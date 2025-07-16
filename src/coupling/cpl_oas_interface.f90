#ifdef COUP_OASIS3MCT

MODULE cpl_oas_interface

!     AUTHOR.
!     -------
!     2020-11: Ha Ho-Hagemann (Hereon): OASIS interface for ICON-CLM

USE cpl_oas_vardef

USE mod_oasis_namcouple        ! OASIS3MCT namcouple variables: e.g. coupling time step
USE mod_oasis,              ONLY: oasis_set_couplcomm

USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message, message_text, finish
USE mo_model_domain,        ONLY: t_patch
USE mo_mpi,                 ONLY: p_pe_work, process_mpi_all_comm, &
                                  my_process_is_stdio, num_work_procs, &
                                  p_io, p_comm_work, p_bcast, abort_mpi, p_barrier, &
                                  get_my_mpi_all_id, mpi_comm_null, &
                                  p_comm_work_only
USE mo_sync,                ONLY: sync_c, sync_patch_array
USE mo_ext_data_types,      ONLY: t_external_data
USE mo_nonhydro_types,      ONLY: t_nh_diag
USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
USE mo_nwp_lnd_types,       ONLY: t_lnd_diag, t_wtr_prog
USE mo_lnd_nwp_config,      ONLY: lseaice, isub_water, isub_seaice, isub_lake, &
                                  ntiles_total, ntiles_lnd, ntiles_water, hice_max, frsi_min
USE mo_impl_constants,      ONLY: min_rlcell_int, grf_bdywidth_c, inwp
USE mo_physical_constants,  ONLY: lh_v=>alv, lh_s=>als
USE mo_parallel_config,     ONLY: idx_1d, nproma
USE mo_loopindices,         ONLY: get_indices_c
USE mo_grid_config,         ONLY: n_dom
USE mo_io_units,            ONLY: find_next_free_unit
USE mo_run_config,          ONLY: msg_level, iforcing, nsteps, ltimer, dtime
USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling, &
       &                          timer_coupling_put, timer_coupling_get,  &
       &                          timer_coupling_init
USE mo_idx_list,            ONLY: t_idx_list_blocked

!==============================================================================

IMPLICIT NONE

!==============================================================================
PUBLIC cpl_oas_send,                &
       cpl_oas_receive,             &
       construct_atmo_coupler_OAS,  &
       destruct_atmo_coupler_OAS

!==============================================================================

CONTAINS


!==============================================================================
SUBROUTINE destruct_atmo_coupler_OAS

INTEGER  :: istat              ! for local error-code

!- End of header
!==============================================================================

  istat = 0

IF ( lpe_cpl ) THEN
  IF (ALLOCATED (ssnd)) THEN
   DEALLOCATE( ssnd, srcv, STAT=istat )
  ENDIF

  IF (ALLOCATED (oas_part)) THEN   
   DEALLOCATE( oas_part, oas_snd_field, oas_rcv_field, oas_rcv_field_icon, STAT=istat )
  ENDIF
ENDIF

!------------------------------------------------------------------------------
!  End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE destruct_atmo_coupler_OAS

!==============================================================================
SUBROUTINE construct_atmo_coupler_OAS(p_patch)

TYPE(t_patch),        TARGET,INTENT(in) :: p_patch(:)    ! pt_patch
!
! local variables
!
TYPE(t_patch), POINTER :: patch_horz

INTEGER :: &
  ji, jj, jg, jn,          & !
  ierrstat,                & !
  izerrstat

INTEGER ::   nuin, ierror

INTEGER :: i_startblk, i_endblk, jb, jc, ic, i_startidx, i_endidx, &
           rl_start, rl_end, c

CHARACTER (LEN=30)         ::   &
  yinput             ! Namelist INPUT file

CHARACTER (LEN=80)         ::   &
  yerrmsg      ! error message

INTEGER :: current_proc
  
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

  debug_oasis = 1

  patch_horz => p_patch(1)

IF ( p_pe_work < patch_horz%proc0) THEN
 lpe_cpl = .FALSE.
 CALL oasis_set_couplcomm(mpi_comm_null)
 WRITE(*,*) 'lpe_cpl=F, mpi_comm_null = ', mpi_comm_null
ELSE
 lpe_cpl = .TRUE.
 CALL oasis_set_couplcomm(p_comm_work_only)
 WRITE(*,*) 'lpe_cpl=T, p_comm_work_only = ', p_comm_work_only
ENDIF
  call flush(6)


  IF (ltimer) CALL timer_start(timer_coupling_init)

! -----------------------------------------------------------------
! ... Setup the OASIS interface
! ----------------------------------------------------------------
 
IF ( my_process_is_stdio().AND. debug_oasis > 15 ) THEN
  print*, ' *****************************************************'
  print*, ' *     ICON-CLM: Setting up OASIS3(-MCT) interface   *'
  print*, ' *****************************************************'
ENDIF

! With OASIS3, all domains involved; for now, this will also be valid for 
! OASIS3-MCT (coupling only on a subset of subdomains becomes important if
! many processors are used)
!IF ( my_process_is_io() .OR. my_process_is_pref() ) THEN

IF ( lpe_cpl ) THEN
! -----------------------------------------------------------------
! ... Define the partition
! -----------------------------------------------------------------

   IF (n_dom > 1) CALL oasis_abort(oas_comp_id, &
      'oas_icon_partition', 'Number of ICON domain > 1 when coupled to CLM.')

    oas_nlat = patch_horz%n_patch_cells_g ! number of points globally
    oas_nlon = 1

    IF ( debug_oasis > 25 ) PRINT *, '++++ ICON-CLM: oas_nlat, oas_nlon=',oas_nlat, oas_nlon
    IF ( debug_oasis > 25 ) PRINT *, '++++ ICON-CLM: p_patch(1)%n_patch_cells=',patch_horz%n_patch_cells
    IF ( debug_oasis > 25 ) PRINT *, '++++ ICON-CLM: grf_bdywidth_c=',grf_bdywidth_c  ! 4
    IF ( debug_oasis > 25 ) PRINT *, '++++ ICON-CLM: min_rlcell_int=',min_rlcell_int  ! -4

!   oas_part(1:N) where N = 2 + number of points
!   oas_part(1)= 4 (indicates a Points partition) !!! Chosen here !!!
!   oas_part(2)= number of points in the partition
!   oas_part(3)= the first global index
!   oas_part(4)= the second global index ...
!   oas_part(N)= the last global index

    !ALLOCATE( oas_part(2+patch_horz%n_patch_cells) )
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = patch_horz%cells%start_blk(rl_start, 1)
    i_endblk   = patch_horz%cells%end_blk(rl_end, MAX(1,patch_horz%n_childdom))

    ! (1) calculate length of oas_part on current proc
    c = 0
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(patch_horz, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      c = c + i_endidx - i_startidx + 1
    END DO
    ALLOCATE( oas_part(c+2) )
    oas_part(1) = 4 ! 4: POINTS partitioning (1: APPLE; 2: BOX; 3: ORANGE; 0: SERIAL)
    oas_part(2) = c  ! number of points in the partition

    ! (2) fill the partition with the cell indices
    c = 0
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(patch_horz, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        c = c + 1
        ic = idx_1d(jc,jb)
        oas_part(c+2) = patch_horz%cells%decomp_info%glb_index(ic)
      END DO
    END DO
    
    CALL oasis_def_partition(oas_part_id, oas_part, oas_error, oas_nlat*oas_nlon)

    IF ( debug_oasis > 25 ) PRINT *, '++++ ICON-CLM: oas_part=', oas_part
    IF ( debug_oasis > 25 ) PRINT *, '++++ ICON-CLM: oas_part_id=', oas_part_id

    IF (oas_error /= 0) &
      CALL oasis_abort(oas_comp_id, oas_comp_name, 'Failure in oasis_def_partition')

ENDIF
! -----------------------------------------------------------------
! ... Initialize some variables
! ----------------------------------------------------------------

  nfld_snd = 0
  nfld_rcv = 0

! -----------------------------------------------------------------
! ... Define list of SENT variables per coupling
! ----------------------------------------------------------------

  IF (my_process_is_stdio() ) THEN
    yinput   = 'namelist_cpl_atm_oce'
    nuin = find_next_free_unit(10,100)

    OPEN (nuin, FILE=TRIM(yinput), FORM='FORMATTED', STATUS='UNKNOWN',  &
    IOSTAT=izerrstat)
    IF (izerrstat /= 0) THEN
     yerrmsg  = ' ERROR    *** Error while opening file namelist_cpl_atm_oce *** '
     ierror   = 2
     RETURN
    ELSE
     print*, "++++ ICON-CLM: opening file ", TRIM(yinput)
    ENDIF
  ENDIF
    
   CALL namelist_cpl (nuin, izerrstat)

   IF (izerrstat > 0) THEN
      yerrmsg  = ' ERROR *** Wrong values occured in NAMELIST group /cpl_nml/ *** '
      ierror   = 3
      RETURN
   ELSEIF (izerrstat < 0) THEN
      PRINT *, ' ERROR while reading NAMELIST group /cpl_nml/ in namelist_cpl '
      ierror   = 4
      RETURN
   ENDIF

   IF (my_process_is_stdio()) THEN
     CLOSE (nuin, STATUS='KEEP', IOSTAT=izerrstat)
     IF (izerrstat /= 0) THEN
      yerrmsg = ' ERROR *** while closing file namelist_cpl_atm_oce *** '
      ierror  = 5
     ENDIF
   ENDIF

! -----------------------------------------------------------------
! ... write variables names for all active couplings to the structures ssnd and srcv
! ----------------------------------------------------------------

IF ( lpe_cpl ) THEN

! initial values
   jprsnd_u10   = 1
   jprsnd_v10   = 2
   jprsnd_swdo  = 3
   jprsnd_swdi  = 4
   jprsnd_lwd   = 5
   jprsnd_t2m   = 6
   jprsnd_q2m   = 7
   jprsnd_pre   = 8
   jprsnd_snw   = 9
   jprsnd_slp   = 10
   jprsnd_clf   = 11
   jprsnd_lhf   = 12 
   jprsnd_shf   = 13
   jprsnd_umoo  = 14
   jprsnd_vmoo  = 15
   jprsnd_umoi  = 16
   jprsnd_vmoi  = 17
   jprsnd_evp   = 18
   jprsnd_sub   = 19
   jprsnd_qnso  = 20
   jprsnd_qnsi  = 21
   jprsnd_dqn   = 22
   jprsnd_wnd   = 23
   jprsnd_ros   = 24
   jprsnd_rog   = 25
   ! Actual number of max snd fields
   nfld_snd_max = 25

   jprrcv_sst   = 1
   jprrcv_ifr   = 2
   jprrcv_ial   = 3
   jprrcv_lhf   = 4
   jprrcv_shf   = 5
   jprrcv_umo   = 6
   jprrcv_vmo   = 7
   jprrcv_tice  = 8
   jprrcv_hice  = 9
   ! Actual number of max rcv fields
   nfld_rcv_max = 9

   ALLOCATE ( ssnd  (nfld_snd_max), STAT=ierrstat )

   ssnd(:)%clname  = 'none'
   ssnd(:)%laction = .FALSE.
   ssnd(:)%nid     = -1

   nfld_snd = 0
! 1
   if (trim(atm_snd_u10) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_u10)%clname = TRIM(atm_snd_u10)
    ssnd(jprsnd_u10)%laction = .TRUE.
   endif
! 2
   if (trim(atm_snd_v10) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_v10)%clname = TRIM(atm_snd_v10)
    ssnd(jprsnd_v10)%laction = .TRUE.
   endif
! 3
   if (trim(atm_snd_swdo) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_swdo)%clname = TRIM(atm_snd_swdo)
    ssnd(jprsnd_swdo)%laction = .TRUE.
   endif
! 4
   if (trim(atm_snd_swdi) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_swdi)%clname = TRIM(atm_snd_swdi)
    ssnd(jprsnd_swdi)%laction = .TRUE.
   endif
! 5
   if (trim(atm_snd_lwd) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_lwd)%clname = TRIM(atm_snd_lwd)
    ssnd(jprsnd_lwd)%laction = .TRUE.
   endif
! 6
   if (trim(atm_snd_t2m) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_t2m)%clname = TRIM(atm_snd_t2m)
    ssnd(jprsnd_t2m)%laction = .TRUE.
   endif
! 7
   if (trim(atm_snd_q2m) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_q2m)%clname = TRIM(atm_snd_q2m)
    ssnd(jprsnd_q2m)%laction = .TRUE.
   endif
! 8
   if (trim(atm_snd_pre) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_pre)%clname = TRIM(atm_snd_pre)
    ssnd(jprsnd_pre)%laction = .TRUE.
   endif
! 9
   if (trim(atm_snd_snw) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_snw)%clname = TRIM(atm_snd_snw)
    ssnd(jprsnd_snw)%laction = .TRUE.
   endif
! 10
   if (trim(atm_snd_slp) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_slp)%clname = TRIM(atm_snd_slp)
    ssnd(jprsnd_slp)%laction = .TRUE.
   endif
! 11
   if (trim(atm_snd_clf) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_clf)%clname = TRIM(atm_snd_clf)
    ssnd(jprsnd_clf)%laction = .TRUE.
   endif
! 12
   if (trim(atm_snd_lhf) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_lhf)%clname = TRIM(atm_snd_lhf)
    ssnd(jprsnd_lhf)%laction = .TRUE.
   endif
! 13
   if (trim(atm_snd_shf) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_shf)%clname = TRIM(atm_snd_shf)
    ssnd(jprsnd_shf)%laction = .TRUE.
   endif
! 14
   if (trim(atm_snd_umoo) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_umoo)%clname = TRIM(atm_snd_umoo)
    ssnd(jprsnd_umoo)%laction = .TRUE.
   endif
! 15
   if (trim(atm_snd_vmoo) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_vmoo)%clname = TRIM(atm_snd_vmoo)
    ssnd(jprsnd_vmoo)%laction = .TRUE.
   endif
! 16
   if (trim(atm_snd_umoi) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_umoi)%clname = TRIM(atm_snd_umoi)
    ssnd(jprsnd_umoi)%laction = .TRUE.
   endif
! 17
   if (trim(atm_snd_vmoi) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_vmoi)%clname = TRIM(atm_snd_vmoi)
    ssnd(jprsnd_vmoi)%laction = .TRUE.
   endif
! 18
   if (trim(atm_snd_evp) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_evp)%clname = TRIM(atm_snd_evp)
    ssnd(jprsnd_evp)%laction = .TRUE.
   endif
! 19
   if (trim(atm_snd_sub) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_sub)%clname = TRIM(atm_snd_sub)
    ssnd(jprsnd_sub)%laction = .TRUE.
   endif
! 20
   if (trim(atm_snd_qnso) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_qnso)%clname = TRIM(atm_snd_qnso)
    ssnd(jprsnd_qnso)%laction = .TRUE.
   endif
! 21
   if (trim(atm_snd_qnsi) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_qnsi)%clname = TRIM(atm_snd_qnsi)
    ssnd(jprsnd_qnsi)%laction = .TRUE.
   endif
! 22
   if (trim(atm_snd_dqn) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_dqn)%clname = TRIM(atm_snd_dqn)
    ssnd(jprsnd_dqn)%laction = .TRUE.
   endif
! 23
   if (trim(atm_snd_wnd) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_wnd)%clname = TRIM(atm_snd_wnd)
    ssnd(jprsnd_wnd)%laction = .TRUE.
   endif

!---------------------------------------------------------------------------
! 24
   if (trim(atm_snd_ros) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_ros)%clname = TRIM(atm_snd_ros)
    ssnd(jprsnd_ros)%laction = .TRUE.
   endif
! 25
   if (trim(atm_snd_rog) .ne. 'none' ) then
    nfld_snd = nfld_snd +1
    ssnd(jprsnd_rog)%clname = TRIM(atm_snd_rog)
    ssnd(jprsnd_rog)%laction = .TRUE.
   endif

! -----------------------------------------------------------------
! ... Define list of RECEIVED variables per coupling
! ----------------------------------------------------------------

   nfld_rcv = 0

   ALLOCATE ( srcv(nfld_rcv_max), STAT=ierrstat )
   srcv(:)%clname  = 'none'
   srcv(:)%laction = .FALSE.
   srcv(:)%nid     = -1

   if (trim(atm_rcv_sst) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_sst)%clname = TRIM(atm_rcv_sst)
    srcv(jprrcv_sst)%laction = .TRUE.
   endif
! 2
   if (trim(atm_rcv_ifr) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_ifr)%clname = TRIM(atm_rcv_ifr)
    srcv(jprrcv_ifr)%laction = .TRUE.
   endif
! 3
   if (trim(atm_rcv_ial) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_ial)%clname = TRIM(atm_rcv_ial)
    srcv(jprrcv_ial)%laction = .TRUE.
   endif
! 4
   if (trim(atm_rcv_lhf) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_lhf)%clname = TRIM(atm_rcv_lhf)
    srcv(jprrcv_lhf)%laction = .TRUE.
   endif
! 5
   if (trim(atm_rcv_shf) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_shf)%clname = TRIM(atm_rcv_shf)
    srcv(jprrcv_shf)%laction = .TRUE.
   endif
! 6
   if (trim(atm_rcv_umo) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_umo)%clname = TRIM(atm_rcv_umo)
    srcv(jprrcv_umo)%laction = .TRUE.
   endif
! 7
   if (trim(atm_rcv_vmo) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_vmo)%clname = TRIM(atm_rcv_vmo)
    srcv(jprrcv_vmo)%laction = .TRUE.
   endif
! 8
   if (trim(atm_rcv_tice) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_tice)%clname = TRIM(atm_rcv_tice)
    srcv(jprrcv_tice)%laction = .TRUE.
   endif
! 9
   if (trim(atm_rcv_hice) .ne. 'none' ) then
    nfld_rcv = nfld_rcv +1
    srcv(jprrcv_hice)%clname = TRIM(atm_rcv_hice)
    srcv(jprrcv_hice)%laction = .TRUE.
   endif

  current_proc = get_my_mpi_all_id()
  IF ( debug_oasis > 15 .AND. current_proc == patch_horz%proc0 ) THEN
    print*, ' nfld_snd=',nfld_snd, ', nfld_rcv=',nfld_rcv
    call flush(6)
  ENDIF


! -----------------------------------------------------------------
! ... Define the shape of the valid region without the halo and overlaps between CPUs
! -----------------------------------------------------------------

    oas_var_nodims(1) = 1   ! rang of field array
    oas_var_nodims(2) = 1   ! 'bundle' always 1 in OASIS3-MCT
    oas_vshape(1)     = 1
    oas_vshape(2)     = oas_part(2)

! -----------------------------------------------------------------
! ... Announce variables to be sent and to be received
! -----------------------------------------------------------------

  ! Announce variables to be sent:
  DO ji = 1, nfld_snd_max
    IF ( ssnd(ji)%laction ) THEN
      CALL oasis_def_var( ssnd(ji)%nid, ssnd(ji)%clname, oas_part_id, &
        oas_var_nodims, OASIS_Out, oas_vshape, OASIS_REAL, oas_error )
      IF ( oas_error /= OASIS_Success ) CALL oasis_abort( ssnd(ji)%nid, &
        'construct_atmo_coupler_OAS', 'Failure in oasis_def_var for '//TRIM(ssnd(ji)%clname) )
      IF(msg_level > 20 ) &
       & WRITE(*,*) 'oasis_def_var with ssnd(',ji,') = ', ssnd(ji)%nid, ssnd(ji)%clname
    ENDIF
  ENDDO
      
  ! Announce variables to be received:
  DO ji = 1, nfld_rcv_max
    IF ( srcv(ji)%laction ) THEN
      CALL oasis_def_var( srcv(ji)%nid, srcv(ji)%clname, oas_part_id, &
        oas_var_nodims, OASIS_In, oas_vshape, OASIS_REAL, oas_error )
      IF ( oas_error /= OASIS_Success ) CALL oasis_abort( srcv(ji)%nid, &
        'construct_atmo_coupler_OAS', 'Failure in oasis_def_var for '//TRIM(srcv(ji)%clname) )
      IF(msg_level > 20 ) &
       & WRITE(*,*) 'oasis_def_var with srcv(',ji,') = ', srcv(ji)%nid, srcv(ji)%clname
    ENDIF
  ENDDO

  ALLOCATE ( oas_snd_field(oas_vshape(1):oas_vshape(2),nfld_snd), stat=oas_error )
  IF (oas_error > 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
      'Failure in allocating icon send buffers' )

  ALLOCATE ( oas_rcv_field(oas_vshape(1):oas_vshape(2),nfld_rcv), stat=oas_error )
  IF (oas_error > 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
      'Failure in allocating icon receive buffers' )

  ALLOCATE ( oas_rcv_field_icon(nproma, patch_horz%nblks_c,nfld_rcv), stat=oas_error )
  IF (oas_error > 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
      'Failure in allocating icon receive buffers' )

  ! initialize buffers
  oas_snd_field      = -1000.
  oas_rcv_field      = -1000.
  oas_rcv_field_icon = -1000._wp
! ------------------------------------------------------------------
! ... End of definition phase (must be called by all processes including 
!     the PEs not involved in the coupling - see also COUP_OASIS3 in io/ ... !! )
! ------------------------------------------------------------------

  CALL oasis_enddef ( oas_error )
  IF ( oas_error /= OASIS_Success ) CALL oasis_abort ( oas_comp_id, &
    'construct_atmo_coupler_OAS', 'Failure in oasis_enddef')

ELSE
  CALL oasis_enddef(oas_error)
  IF ( oas_error /= OASIS_Success ) CALL oasis_abort ( oas_comp_id, &
    'construct_atmo_coupler_OAS', 'Failure in oasis_enddef (on non-worker process)')

ENDIF ! lpe_cpl

  IF (ltimer) CALL timer_stop(timer_coupling_init)

END SUBROUTINE construct_atmo_coupler_OAS

!==============================================================================
!+ input of NAMELIST namelist_cpl_atm_oce
!------------------------------------------------------------------------------

SUBROUTINE namelist_cpl (nuin, ierrstat)

! Parameter list:
  INTEGER , INTENT (IN)      ::        &
    nuin            ! Unit number for Namelist INPUT file

  INTEGER , INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

  CHARACTER (LEN=80)  ::        &
    yerrmsg            ! error message

  CHARACTER (LEN=100)     , ALLOCATABLE ::   charbuf (:)
    
! Local variables: 
  INTEGER   ::  &
    i, iz_err
  
! Define the namelist group
  NAMELIST /cpl_nml/            &
    atm_snd_u10, atm_snd_v10,   &
    atm_snd_swdo, atm_snd_swdi, &
    atm_snd_lwd,                &
    atm_snd_t2m, atm_snd_q2m,   &
    atm_snd_pre, atm_snd_snw,   &
    atm_snd_slp, atm_snd_clf,   &
    atm_snd_lhf, atm_snd_shf,   &
    atm_snd_umoo, atm_snd_vmoo, &
    atm_snd_umoi, atm_snd_vmoi, &
    atm_snd_evp, atm_snd_sub,   &
    atm_snd_qnso, atm_snd_qnsi, &
    atm_snd_dqn,                &
    atm_snd_wnd,                &
    atm_snd_ros, atm_snd_rog,   &
    atm_rcv_sst,                &
    atm_rcv_ifr, atm_rcv_ial,   &
    atm_rcv_lhf, atm_rcv_shf,   &
    atm_rcv_umo, atm_rcv_vmo,   &
    atm_rcv_tice, atm_rcv_hice
  
!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin SUBROUTINE input_oasisctl
!------------------------------------------------------------------------------

ierrstat = 0
iz_err   = 0

!------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!------------------------------------------------------------------------------

IF (my_process_is_stdio()) THEN
  atm_snd_u10 = 'none'
  atm_snd_v10 = 'none'
  atm_snd_swdo = 'none'
  atm_snd_swdi = 'none'
  atm_snd_lwd = 'none'
  atm_snd_t2m = 'none'
  atm_snd_q2m = 'none'
  atm_snd_pre = 'none'
  atm_snd_snw = 'none'
  atm_snd_slp = 'none'
  atm_snd_clf = 'none'
  atm_snd_lhf = 'none'
  atm_snd_shf = 'none'
  atm_snd_umoo = 'none'
  atm_snd_vmoo = 'none'
  atm_snd_umoi = 'none'
  atm_snd_vmoi = 'none'
  atm_snd_evp = 'none'
  atm_snd_sub = 'none'
  atm_snd_qnso = 'none'
  atm_snd_qnsi = 'none'
  atm_snd_dqn = 'none'
  atm_snd_wnd = 'none'
  atm_snd_ros = 'none'
  atm_snd_rog = 'none'
  atm_rcv_sst = 'none'
  atm_rcv_ifr = 'none'
  atm_rcv_ial = 'none'
  atm_rcv_lhf = 'none'
  atm_rcv_shf = 'none'
  atm_rcv_umo = 'none'
  atm_rcv_vmo = 'none'
  atm_rcv_tice = 'none'
  atm_rcv_hice = 'none'
                     
!------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!------------------------------------------------------------------------------

  READ (nuin, cpl_nml, IOSTAT=iz_err)

ENDIF ! my_process_is_stdio()

IF (num_work_procs > 1) THEN
  ! distribute error status to all processors
  CALL p_bcast(iz_err, p_io, p_comm_work)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

!------------------------------------------------------------------------------
!- Section 4: Distribute variables to all nodes
!------------------------------------------------------------------------------

IF (num_work_procs > 1) THEN
  ALLOCATE (charbuf(34))

  IF (my_process_is_stdio()) THEN
   charbuf( 1) = atm_snd_u10
   charbuf( 2) = atm_snd_v10
   charbuf( 3) = atm_snd_swdo
   charbuf( 4) = atm_snd_lwd
   charbuf( 5) = atm_snd_t2m
   charbuf( 6) = atm_snd_q2m
   charbuf( 7) = atm_snd_pre
   charbuf( 8) = atm_snd_snw
   charbuf( 9) = atm_snd_slp
   charbuf(10) = atm_snd_clf
   charbuf(11) = atm_snd_lhf
   charbuf(12) = atm_snd_shf
   charbuf(13) = atm_snd_umoo
   charbuf(14) = atm_snd_vmoo
   charbuf(15) = atm_snd_evp
   charbuf(16) = atm_snd_sub
   charbuf(17) = atm_snd_qnso
   charbuf(18) = atm_snd_dqn
   charbuf(19) = atm_snd_wnd
   charbuf(20) = atm_snd_ros
   charbuf(21) = atm_snd_rog
   charbuf(22) = atm_rcv_sst
   charbuf(23) = atm_rcv_ifr
   charbuf(24) = atm_rcv_ial
   charbuf(25) = atm_rcv_lhf
   charbuf(26) = atm_rcv_shf
   charbuf(27) = atm_rcv_umo
   charbuf(28) = atm_rcv_vmo
   charbuf(29) = atm_snd_swdi
   charbuf(30) = atm_snd_qnsi
   charbuf(31) = atm_snd_umoi
   charbuf(32) = atm_snd_vmoi
   charbuf(33) = atm_rcv_tice
   charbuf(34) = atm_rcv_hice


  ENDIF ! my_process_is_stdio()p_bcast
  
  CALL p_bcast(charbuf, p_io, p_comm_work)

  IF (.NOT. my_process_is_stdio()) THEN
!  IF (my_cart_id /= 0) THEN
   atm_snd_u10  = charbuf( 1)
   atm_snd_v10  = charbuf( 2)
   atm_snd_swdo = charbuf( 3)
   atm_snd_lwd  = charbuf( 4)
   atm_snd_t2m  = charbuf( 5)
   atm_snd_q2m  = charbuf( 6)
   atm_snd_pre  = charbuf( 7)
   atm_snd_snw  = charbuf( 8)
   atm_snd_slp  = charbuf( 9)
   atm_snd_clf  = charbuf(10)
   atm_snd_lhf  = charbuf(11)
   atm_snd_shf  = charbuf(12)
   atm_snd_umoo = charbuf(13)
   atm_snd_vmoo = charbuf(14)
   atm_snd_evp  = charbuf(15)
   atm_snd_sub  = charbuf(16)
   atm_snd_qnso = charbuf(17)
   atm_snd_dqn  = charbuf(18)
   atm_snd_wnd  = charbuf(19)
   atm_snd_ros  = charbuf(20)
   atm_snd_rog  = charbuf(21)
   atm_rcv_sst  = charbuf(22)
   atm_rcv_ifr  = charbuf(23)
   atm_rcv_ial  = charbuf(24)
   atm_rcv_lhf  = charbuf(25)
   atm_rcv_shf  = charbuf(26)
   atm_rcv_umo  = charbuf(27)
   atm_rcv_vmo  = charbuf(28)
   atm_snd_swdi = charbuf(29)
   atm_snd_qnsi = charbuf(30)
   atm_snd_umoi = charbuf(31)
   atm_snd_vmoi = charbuf(32)
   atm_rcv_tice = charbuf(33)
   atm_rcv_hice = charbuf(34)

  ENDIF ! my_process_is_stdio()

  DEALLOCATE (charbuf)

  IF ( ierrstat /= 0 ) THEN
    PRINT *, ' ERROR *** in distributing buffers in namelist_cpl*** '
    RETURN
  ENDIF

ENDIF ! num_work_procs

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE namelist_cpl 

!==============================================================================


!==============================================================================

SUBROUTINE cpl_oas_receive   (p_patch       ,  &
                              p_diag_lnd    ,  &
                              p_prog_wtr_now,  &
                              p_prog_wtr_new,  &
                              p_list_sea_nemo, &
                              sim_time,        &
                              lcpl_hice        )
!!---------------------------------------------------------------------
!!              ***  ROUTINE cpl_oas_receive  ***
!!
!! ** Purpose : - Prepare and receive coupling fields to OASIS
!!    
!!----------------------------------------------------------------------



  TYPE(t_patch),            TARGET,INTENT(in)    :: p_patch           ! p_patch(1)
  TYPE(t_idx_list_blocked), TARGET,INTENT(in)    :: p_list_sea_nemo   ! ext_data(jg)%atm%list_sea_nemo
  TYPE(t_lnd_diag),                INTENT(inout) :: p_diag_lnd
  TYPE(t_wtr_prog),                INTENT(inout) :: p_prog_wtr_now !< prog vars for wtr
  TYPE(t_wtr_prog),                INTENT(inout) :: p_prog_wtr_new !< prog vars for wtr

  REAL(wp), INTENT(IN) :: sim_time         !< elapsed simulation time

  LOGICAL, INTENT(out) :: lcpl_hice

!
! local parameters, variables and arrays
!

INTEGER :: isec, jn, ii, izerror

CHARACTER (LEN=25)   :: yzroutine
CHARACTER (LEN=1024) :: oas_message

INTEGER :: rl_start, rl_end, i_startblk, i_endblk, jb, jc, &
           ic, jlac, nn, i_startidx, i_endidx

INTEGER :: kinfo, nrcvinfo(nfld_rcv)

INTEGER :: current_proc

LOGICAL :: lcpl_tice

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

  IF (ltimer) CALL timer_start(timer_coupling)

  oas_message = ''
  yzroutine= 'cpl_oas_receive'

! Coupling only with working PE's
IF ( lpe_cpl ) THEN

  isec = INT(sim_time)

  current_proc = get_my_mpi_all_id()
  IF ( debug_oasis > 15 .AND. current_proc == p_patch%proc0 ) THEN
    print*, 'in cpl_oas_receive: sim_time =', isec
    call flush(6)
  ENDIF

!------------------------------------------------------------------------------
! 2. Receive all coupling fields (independent of the specific coupling) 
!-------------------------------------------------------------------------------

  nrcvinfo (:) = OASIS_idle

  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int
  i_startblk = p_patch%cells%start_blk(rl_start, 1)
  i_endblk   = p_patch%cells%end_blk(rl_end, MAX(1,p_patch%n_childdom))

  jlac = 0
  lcpl_hice = .false.
  lcpl_tice = .false.

  IF (ltimer) CALL timer_start(timer_coupling_get)

  DO jn = 1, nfld_rcv_max

    IF( srcv(jn)%laction) THEN

      jlac = jlac + 1

      CALL oasis_get(srcv(jn)%nid, isec, oas_rcv_field(:,jlac), kinfo)

      IF (kinfo .NE. OASIS_Ok .AND. kinfo .LT. OASIS_Recvd) THEN
        WRITE(oas_message,*) 'Failure in oasis_get of ', srcv(jn)%clname
        CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
      ENDIF

      IF( kinfo == OASIS_Recvd   .OR. kinfo == OASIS_FromRest    .OR. &
          kinfo == OASIS_RecvOut .OR. kinfo == OASIS_FromRestOut .OR. & 
          kinfo == OASIS_Input) nrcvinfo(jlac) = OASIS_Rcv
       
      ! transform thru oasis recieved variable into something icon can use
      IF( nrcvinfo(jlac) == OASIS_Rcv ) THEN
        ic = 0
        DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, &
           &                 i_endidx, rl_start, rl_end)
          DO jc = i_startidx, i_endidx
            ic = ic + 1
            oas_rcv_field_icon(jc,jb,jlac) = oas_rcv_field(ic,jlac)
          END DO
        END DO

        ! check:
        IF(msg_level > 20 ) THEN
          WRITE(*,*) "oasis_get of ", srcv(jn)%clname, ": min, max = ", &
           & MINVAL(oas_rcv_field(:,jlac)), MAXVAL(oas_rcv_field(:,jlac))
          CALL flush(6)
        ENDIF

      ENDIF ! nrcvinfo(jn) == OASIS_Rcv
       
    ENDIF ! %laction = .TRUE.

  END DO ! jn

  IF (ltimer) CALL timer_stop(timer_coupling_get)


!----------------------------------------------------------------------------
! 2.4 nemoI-specific handling of received fields
!----------------------------------------------------------------------------
  jlac = 0

  DO jn = 1, nfld_rcv_max

! 1
  !____ SST
    IF( jn == jprrcv_sst .AND. srcv(jn)%laction ) THEN

      jlac = jlac + 1

      nn = 0
      DO jb = i_startblk, i_endblk
        !CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO ic = 1, p_list_sea_nemo%ncount(jb)
          jc = p_list_sea_nemo%idx(ic,jb)

          !IF (oas_rcv_field(nn+jc,jlac) .gt. 100._wp .AND. oas_rcv_field(nn+jc,jlac) .lt. 400._wp ) &
          ! & p_diag_lnd%t_seasfc  (jc,jb)  = oas_rcv_field(nn+jc,jlac)
          IF (oas_rcv_field_icon(jc,jb,jlac) .gt. 100._wp .AND. oas_rcv_field_icon(jc,jb,jlac) .lt. 400._wp ) &
           & p_diag_lnd%t_seasfc  (jc,jb)  = oas_rcv_field_icon(jc,jb,jlac)

        END DO
        !nn = nn + (i_endidx-i_startidx+1)
      END DO

      ! vm: needed?
      !CALL sync_patch_array(sync_c, p_patch, p_diag_lnd%t_seasfc(:,:) )

! 2
  !____ fr_ice
!  print*, "++++ ICON-CLM receives: jprrcv_ifr=",jprrcv_ifr,srcv(jprrcv_ifr)%laction
    ELSEIF( jn == jprrcv_ifr .AND. srcv(jn)%laction ) THEN

      jlac = jlac + 1

      DO jb = i_startblk, i_endblk
        !CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO ic = 1, p_list_sea_nemo%ncount(jb)
          jc = p_list_sea_nemo%idx(ic,jb)
          ! frsi_min = 1.0E-10_wp
          p_diag_lnd%fr_seaice (jc,jb) = MIN( MAX(oas_rcv_field_icon(jc,jb,jlac), 1.0E-12_wp), 1._wp)

        END DO
        
        !nn = nn + (i_endidx-i_startidx+1)
      END DO
      !CALL sync_patch_array(sync_c, p_patch, p_diag_lnd%fr_seaice(:,:) )

! 3
  !____ alb_ice
    !!! make sure to set lprog_albsi = true !!!
    ! print*, "++++ ICON-CLM receives: jprrcv_ial=",jprrcv_ial,srcv(jprrcv_ial)%laction
    ELSEIF( jn == jprrcv_ial .AND. srcv(jn)%laction ) THEN

      jlac = jlac + 1

      DO jb = i_startblk, i_endblk

        DO ic = 1, p_list_sea_nemo%ncount(jb)
           jc = p_list_sea_nemo%idx(ic,jb)

           IF (oas_rcv_field_icon(jc,jb,jlac) .gt. 0._wp .AND. oas_rcv_field_icon(jc,jb,jlac) .le. 1._wp) THEN
             p_prog_wtr_now%alb_si(jc,jb)  = oas_rcv_field_icon(jc,jb,jlac)
             p_prog_wtr_new%alb_si(jc,jb)  = p_prog_wtr_now%alb_si(jc,jb)
           END IF
        END DO
      END DO
      !CALL sync_patch_array(sync_c, p_patch, p_prog_wtr_now%alb_si(:,:) )
      !CALL sync_patch_array(sync_c, p_patch, p_prog_wtr_new%alb_si(:,:) )

! 8
  !____ t_ice
    ! print*, "++++ ICON-CLM receives: jprrcv_tice=",jprrcv_tice,srcv(jprrcv_tice)%laction
    ELSEIF( jn == jprrcv_tice .AND. srcv(jn)%laction ) THEN

      jlac = jlac + 1
      lcpl_tice = .true.

      DO jb = i_startblk, i_endblk

        DO ic = 1, p_list_sea_nemo%ncount(jb)
          jc = p_list_sea_nemo%idx(ic,jb)

          IF (oas_rcv_field_icon(jc,jb,jlac) .gt. 100._wp .AND. oas_rcv_field_icon(jc,jb,jlac) .lt. 400._wp)  THEN
            p_prog_wtr_now%t_ice (jc,jb)  = oas_rcv_field_icon(jc,jb,jlac)
            p_prog_wtr_new%t_ice (jc,jb)  = oas_rcv_field_icon(jc,jb,jlac)
           ENDIF

        END DO
      END DO
      !CALL sync_patch_array(sync_c, p_patch, p_prog_wtr_now%t_ice(:,:) )
      !CALL sync_patch_array(sync_c, p_patch, p_prog_wtr_new%t_ice(:,:) )

! 9
  !____ h_ice
    ! print*, "++++ ICON-CLM receives: jprrcv_tice=",jprrcv_tice,srcv(jprrcv_tice)%laction
    ELSEIF( jn == jprrcv_hice .AND. srcv(jn)%laction ) THEN

      jlac = jlac + 1
      lcpl_hice = .true.

      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO ic = 1, p_list_sea_nemo%ncount(jb)
          jc = p_list_sea_nemo%idx(ic,jb)

          IF (oas_rcv_field_icon(jc,jb,jlac) .ge. 0._wp) THEN
            p_prog_wtr_new%h_ice (jc,jb)  = MIN (oas_rcv_field_icon(jc,jb,jlac), hice_max-1.0E-5_wp)
            p_prog_wtr_now%h_ice (jc,jb)  = p_prog_wtr_new%h_ice (jc,jb)
          ENDIF

        END DO
      END DO
      !CALL sync_patch_array(sync_c, p_patch, p_prog_wtr_now%h_ice(:,:) )
      !CALL sync_patch_array(sync_c, p_patch, p_prog_wtr_new%h_ice(:,:) )

    ENDIF

    IF( srcv(jn)%laction .AND. msg_level > 20 ) THEN
      WRITE(*,*) jn, 'Receiving  ', srcv(jn)%clname, 'at t=', isec, &
        &        MINVAL(oas_rcv_field_icon(:,:,jlac)), MAXVAL(oas_rcv_field_icon(:,:,jlac))
      flush(6)
    ENDIF

    ! vm: we could also use p_prog_wtr_new%t_snow_si and p_prog_wtr_new%h_snow_si (OSnwTck);
    !     however, currently these are set by the seaice scheme, but without any effect

  END DO ! jn = 1, nfld_rcv_max

    ! if one of lcpl_hice or lcpl_tice is true, the other has to be true as well
    ! abort if not
    IF (lcpl_hice .OR. lcpl_tice) THEN
      IF (.NOT. (lcpl_hice .AND. lcpl_tice)) THEN
        CALL finish('cpl_oas_interface:','hice and tice can only be coupled together!')
      ENDIF
    ENDIF

ENDIF ! lpe_cpl

  IF (ltimer) CALL timer_stop(timer_coupling)

END SUBROUTINE cpl_oas_receive

!==============================================================================

SUBROUTINE cpl_oas_send      (p_patch       ,&
                              prm_diag      ,&
                              p_diag        ,&
                              p_diag_lnd    ,&
                              ext_data      ,&
                              sim_time       )

!!---------------------------------------------------------------------
!!              ***  ROUTINE cpl_oas_send  ***
!!
!! ** Purpose : Prepare and send coupling fields to OASIS
!!      
!!----------------------------------------------------------------------

  TYPE(t_patch),        TARGET,INTENT(in) :: p_patch    ! pt_patch
  TYPE(t_nwp_phy_diag),        INTENT(in) :: prm_diag   ! prm_diag
  TYPE(t_nh_diag),      TARGET,INTENT(in) :: p_diag     ! pt_diag
  TYPE(t_lnd_diag),            INTENT(in) :: p_diag_lnd ! lnd_diag
  TYPE(t_external_data),       INTENT(in) :: ext_data   ! ext_data

REAL(wp), INTENT(IN) :: sim_time         !< elapsed simulation time

!
! local parameters, variables and arrays
!

INTEGER :: &
  isec,    &
  jn, jb, jc,  &
  ii,      &
  izerror, c, jact

CHARACTER (LEN=25)   :: yzroutine
CHARACTER (LEN=1024) :: oas_message

INTEGER :: rl_start, rl_end, i_startblk, i_endblk
INTEGER , ALLOCATABLE :: i_startidx(:), i_endidx(:)
INTEGER :: isubs
REAL(wp):: runoff_s_inst, runoff_g_inst
INTEGER :: current_proc
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

  oas_message = ''
  yzroutine= 'cpl_oas_send'

!----------------------------------------------------------------------------
! 1 coupling-specific handling of sent fields
!----------------------------------------------------------------------------
  IF (ltimer) CALL timer_start(timer_coupling)

! Coupling only on PE with at least one unmasked grid point
IF ( lpe_cpl ) THEN

  isec = INT(sim_time)

  current_proc = get_my_mpi_all_id()
  IF ( debug_oasis > 15 .AND. current_proc == p_patch%proc0 ) THEN
    print*, 'in cpl_oas_send: sim_time = ', isec
  ENDIF
  call flush(6)

!----------------------------------------------------------------------------
! 1.4 nemoI-specific handling of sent fields
!----------------------------------------------------------------------------
! Ha Ho-Hagemann for nemoI {
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start, 1)
    i_endblk   = p_patch%cells%end_blk(rl_end, MAX(1,p_patch%n_childdom))

    ALLOCATE(i_startidx(i_endblk-i_startblk+1))
    ALLOCATE(i_endidx(i_endblk-i_startblk+1))

    c = 1
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx(c), &
        i_endidx(c), rl_start, rl_end)
      c = c + 1
    ENDDO

    jact = 0

    DO jn = 1, nfld_snd_max

! 1
! U10MtNB
        IF( jn == jprsnd_u10 .AND. ssnd(jprsnd_u10)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%u_10m(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 2
! V10MtNB
        IF( jn == jprsnd_v10 .AND. ssnd(jprsnd_v10)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%v_10m(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 3a
! SWDNOtNB: SOBS = prm_diag%swflxsfc -> vm: solar budget over the ocean for flux coupling, not SW down
        IF( jn == jprsnd_swdo .AND. ssnd(jprsnd_swdo)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%swflxsfc_t(jc,jb,isub_water) ! vm: net SW, ocean tile
            ENDDO
          ENDDO
        ENDIF

! 3b (assign nr 20 while reading from the namelist)
! SWDNItNB: SOBS = prm_diag%swflxsfc -> vm: solar budget over ice
        IF( jn == jprsnd_swdi .AND. ssnd(jprsnd_swdi)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%swflxsfc_t(jc,jb,isub_seaice) ! vm: net SW, seaice tile
            ENDDO
          ENDDO
        ENDIF
! 4
! LWDNtNB
        IF( jn == jprsnd_lwd .AND. ssnd(jprsnd_lwd)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%lwflxsfc(jc,jb) &
                                     + prm_diag%lwflx_up_sfc(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 5
! T2MtNB
        IF( jn == jprsnd_t2m .AND. ssnd(jprsnd_t2m)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%t_2m(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 6
! QV2MtNB
        IF( jn == jprsnd_q2m .AND. ssnd(jprsnd_q2m)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%qv_2m(jc,jb)
            ENDDO
          ENDDO
         ENDIF
! 7
! RAINtNB
        IF( jn == jprsnd_pre .AND. ssnd(jprsnd_pre)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%rain_con_rate(jc,jb) &
                                      + prm_diag%rain_gsp_rate(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 8
! SNOWtNB
        IF( jn == jprsnd_snw .AND. ssnd(jprsnd_snw)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%snow_con_rate(jc,jb) &
                                    + prm_diag%snow_gsp_rate(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 9
! PMSL
        IF( jn == jprsnd_slp .AND. ssnd(jprsnd_slp)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              !oas_snd_field(c,jact) = p_diag%pres_sfc(jc,jb)
              oas_snd_field(c,jact) = p_diag%pres_msl(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 10
! CLCT
        IF( jn == jprsnd_clf .AND. ssnd(jprsnd_clf)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%clct(jc,jb)  
            ENDDO
          ENDDO
        ENDIF
! 11
! LHFLtNB
        IF( jn == jprsnd_lhf .AND. ssnd(jprsnd_lhf)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%lhfl_s(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 12
! SHFLtNB
        IF( jn == jprsnd_shf .AND. ssnd(jprsnd_shf)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%shfl_s(jc,jb)
            ENDDO
          ENDDO
        ENDIF
! 13
! UMFLSOtNB - u momentum flux on ocean tile
        IF( jn == jprsnd_umoo .AND. ssnd(jprsnd_umoo)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%umfl_s_t(jc,jb,isub_water)
            ENDDO
          ENDDO
        ENDIF
! 14
! VMFLSOtNB - v momentum flux on ocean tile
        IF( jn == jprsnd_vmoo .AND. ssnd(jprsnd_vmoo)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%vmfl_s_t(jc,jb,isub_water)
            ENDDO
          ENDDO
        ENDIF
! 13b
! UMFLSItNB - u momentum flux on ice tile
        IF( jn == jprsnd_umoi .AND. ssnd(jprsnd_umoi)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%umfl_s_t(jc,jb,isub_seaice)
            ENDDO
          ENDDO
        ENDIF
! 14b
! VMFLSItNB - v momentum flux on ice tile
        IF( jn == jprsnd_vmoi .AND. ssnd(jprsnd_vmoi)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%vmfl_s_t(jc,jb,isub_seaice)
            ENDDO
          ENDDO
        ENDIF
! 15
! EVAPtNB
        IF( jn == jprsnd_evp .AND. ssnd(jprsnd_evp)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = -prm_diag%lhfl_s(jc,jb)/lh_v
            ENDDO
          ENDDO
        ENDIF
! 16
! SUBLtNB
        IF( jn == jprsnd_sub .AND. ssnd(jprsnd_sub)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = -prm_diag%qhfl_s_t (jc,jb,isub_seaice)
            ENDDO
          ENDDO
        ENDIF
! 17
! NSOBSOtNB = thbs + lhfl_s + shfl_s on ocean tile
        IF( jn == jprsnd_qnso .AND. ssnd(jprsnd_qnso)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%lwflxsfc_t(jc,jb,isub_water) &  ! thbs
                                      + prm_diag%lhfl_s_t(jc,jb,isub_water)     &  ! lhfl_s
                                      + prm_diag%shfl_s_t(jc,jb,isub_water)        ! shfl_s
            ENDDO
          ENDDO
        ENDIF
! 17b
! NSOBSItNB = thbs + lhfl_s + shfl_s on ice tile
        IF( jn == jprsnd_qnsi .AND. ssnd(jprsnd_qnsi)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = prm_diag%lwflxsfc_t(jc,jb,isub_seaice) &  ! thbs
                                      + prm_diag%lhfl_s_t(jc,jb,isub_seaice)     &  ! lhfl_s
                                      + prm_diag%shfl_s_t(jc,jb,isub_seaice)        ! shfl_s
            ENDDO
          ENDDO
        ENDIF
! 18
! NSOSENtNB: nonsolar radiation flux sensitivity at the ground d(Qns)/dT ( W/m2/K)
        IF( jn == jprsnd_dqn .AND. ssnd(jprsnd_dqn)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              ! EM Clearly a problem for SI3 
              oas_snd_field(c,jact) = 0._wp
            ENDDO
          ENDDO
        ENDIF
! 19
! WND10tNB
        IF( jn == jprsnd_wnd .AND. ssnd(jprsnd_wnd)%laction ) THEN
          jact = jact + 1
          c = 0
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              oas_snd_field(c,jact) = &
                   SQRT( prm_diag%u_10m(jc,jb)**2 + prm_diag%v_10m(jc,jb)**2 )
            ENDDO
          ENDDO
        ENDIF
 
! 20
! RUNOFF_S: sum over all surface tiles ; convert to mm/s
       IF( jn == jprsnd_ros .AND. ssnd(jprsnd_ros)%laction ) THEN
          jact = jact + 1
          c = 0
          oas_snd_field(:,jact) = 0._wp
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              DO isubs = 1, ntiles_total
                ! Take P-E over the lake as runoff
                IF ( isubs == isub_lake ) THEN
                  oas_snd_field(c,jact) = oas_snd_field(c,jact) +  &
                     &                    ( prm_diag%tot_prec_rate(jc,jb) + prm_diag%qhfl_s_t(jc,jb,isubs) )  &
                     &                    * ext_data%atm%frac_t(jc,jb,isubs)
                ENDIF
                oas_snd_field(c,jact) = oas_snd_field(c,jact) +  &
                     &                  p_diag_lnd%runoff_s_inst_t (jc, jb, isubs) / dtime * ext_data%atm%frac_t(jc,jb,isubs)
              ENDDO
              ! vm: don't exclude negative runoff as P<E is possible
              ! IF (oas_snd_field(c,jact) .LT. 10e-10_wp) oas_snd_field(c,jact) = 0._wp
            ENDDO
          ENDDO
       ENDIF

! 21
! RUNOFF_G: sum over all surface tiles ; convert to mm/s
       IF( jn == jprsnd_rog .AND. ssnd(jprsnd_rog)%laction ) THEN
          jact = jact + 1
          c = 0
          oas_snd_field(:,jact) = 0._wp
          DO jb = i_startblk, i_endblk
            DO jc = i_startidx(jb-i_startblk+1), i_endidx(jb-i_startblk+1)
              c = c + 1
              DO isubs = 1, ntiles_total
                oas_snd_field(c,jact) = oas_snd_field(c,jact) +  &
                 &                      p_diag_lnd%runoff_g_inst_t (jc, jb, isubs) / dtime * ext_data%atm%frac_t(jc,jb,isubs)
              ENDDO
              IF (oas_snd_field(c,jact) .LT. 10e-10_wp) oas_snd_field(c,jact) = 0._wp
            ENDDO
          ENDDO
       ENDIF

    ENDDO ! loop jn
!
    DEALLOCATE(i_startidx,i_endidx)

  IF (ltimer) CALL timer_start(timer_coupling_put)

    jact = 0

    DO jn = 1, nfld_snd_max
      IF (ssnd(jn)%laction ) THEN
        jact = jact + 1

        CALL oasis_put(ssnd(jn)%nid, isec, oas_snd_field(:,jact), oas_error)
        IF( msg_level > 20 ) THEN
          WRITE(oas_message,*) jn, 'Sending  ', ssnd(jn)%clname, 'at t=', isec, &
            &                  MINVAL(oas_snd_field(:,jact)), MAXVAL(oas_snd_field(:,jact))
          !CALL message(yzroutine, oas_message)
          WRITE(*,*) oas_message
          flush(6)
        ENDIF
        IF (oas_error .NE. OASIS_Ok .AND. oas_error .LT. OASIS_Sent) THEN
          WRITE(oas_message,*) 'Failure in oasis_put of ', ssnd(jn)%clname
          CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
        END IF
      END IF
    ENDDO

  IF (ltimer) CALL timer_stop(timer_coupling_put)

ENDIF ! lpe_cpl

  IF (ltimer) CALL timer_stop(timer_coupling)


END SUBROUTINE cpl_oas_send

!==============================================================================


END MODULE cpl_oas_interface

#endif
