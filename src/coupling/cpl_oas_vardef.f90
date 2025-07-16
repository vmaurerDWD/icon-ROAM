#ifdef COUP_OASIS3MCT

MODULE cpl_oas_vardef

!**** *cpl_oas_vardef * - Controls, definitions and variables
!                         for OASIS communications
!
!     AUTHOR.
!     -------
!     2020-11: Ha Ho-Hagemann (Hereon): OASIS interface for ICON-CLM
!****

USE mo_kind,     ONLY: wp

USE mod_oasis                    ! OASIS3-MCT module

IMPLICIT NONE

PUBLIC

! Debug level of OASIS
!     0 : Minimum debugging
!     1 : Debugging
!     2 : Perfs measurement
!     3 : OASIS restart production

! Stefan Poll {
INTEGER               :: oas_comp_id
!CHARACTER(len=7)      :: oas_comp_name="iconclm"
CHARACTER(len=4)      :: oas_comp_name
INTEGER               :: kl_comm              ! local communicator 
INTEGER               :: oas_error            ! return error code
INTEGER               :: oas_nlat
INTEGER               :: oas_nlon
INTEGER               :: oas_var_nodims(2)
INTEGER               :: oas_vshape(2)
INTEGER               :: oas_part_id
INTEGER, ALLOCATABLE  :: oas_part(:)

!TYPE :: t_oas_field
!  CHARACTER(len = 8)  :: clpname
!  INTEGER             :: vid
!END TYPE t_oas_field
!TYPE(t_oas_field), DIMENSION(11)  :: oas_snd_meta
!TYPE(t_oas_field), DIMENSION(8)  :: oas_rcv_meta

REAL, ALLOCATABLE     :: oas_snd_field(:,:), oas_rcv_field(:,:)
REAL(wp), ALLOCATABLE :: oas_rcv_field_icon(:,:,:)

!PUBLIC :: &
  !oas_comp_id, oas_comp_name, kl_comm, oas_error, oas_nlat, &
  !oas_nlon, oas_var_nodims, oas_vshape, oas_part, &
  !oas_snd_field, oas_rcv_field  !, oas_snd_meta, oas_rcv_meta
!! Stefan Poll }

INTEGER :: &
  debug_oasis                ! Debug level of OASIS

INTEGER :: &
  OASIS_Rcv  = 1,          & ! return code if received field
  OASIS_idle = 0,          & ! return code if nothing was done by OASIS
  OASIS_Success = 0          ! return code if no error in OASIS

!! Ha Ho-Hagemann added {
!   ! OASIS Variables not used. defined only for compilation purpose
!   INTEGER                    ::   OASIS_Out         = -1
!   INTEGER                    ::   OASIS_REAL        = -1
!   INTEGER                    ::   OASIS_Ok          = -1
!   INTEGER                    ::   OASIS_In          = -1
!   INTEGER                    ::   OASIS_Sent        = -1
!   INTEGER                    ::   OASIS_SentOut     = -1
!   INTEGER                    ::   OASIS_ToRest      = -1
!   INTEGER                    ::   OASIS_ToRestOut   = -1
!   INTEGER                    ::   OASIS_Recvd       = -1
!   INTEGER                    ::   OASIS_RecvOut     = -1
!   INTEGER                    ::   OASIS_FromRest    = -1
!   INTEGER                    ::   OASIS_FromRestOut = -1
!! Ha Ho-Hagemann added }

TYPE :: FLD_CPL                        ! Type for coupling field information
  LOGICAL                 :: laction   ! To be coupled or not
  CHARACTER(LEN = 16)     :: clname    ! Name of the coupling field
  INTEGER                 :: nid       ! Id of the field
END TYPE FLD_CPL

TYPE(FLD_CPL), ALLOCATABLE :: &
  srcv(:),                 & ! All fields to be received
  ssnd(:)                    ! All fields to be sent

LOGICAL :: &
  lpe_cpl = .FALSE.


INTEGER :: &
  nfld_snd_max,        & ! total number of fields to be sent
  nfld_rcv_max,        & ! total number of fields to be received
  nfld_snd,            & !
  nfld_rcv               ! number of fields to be received per coupling with lsm

! variable names as defined in the namcouple
! for sent fields, copied to srcv(jprsnd_*)%clname
CHARACTER(len = 10) ::        &
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
  atm_snd_dqn, atm_snd_wnd,   &
  atm_snd_ros, atm_snd_rog
! for received fields, copied to srcv(jprrcv_*)%clname
CHARACTER(len = 10) ::        &
  atm_rcv_sst,                &
  atm_rcv_ifr, atm_rcv_ial,   &
  atm_rcv_lhf, atm_rcv_shf,   &
  atm_rcv_umo, atm_rcv_vmo,   &
  atm_rcv_tice, atm_rcv_hice

! integer list for correct assignment of the internally defined variable list
! and the order in which OASIS is sending/receiving the fields
! (fix numbers in construct_atmo_coupler_OAS)
! for sent fields
INTEGER ::                  &
  jprsnd_u10, jprsnd_v10,   &
  jprsnd_swdo, jprsnd_swdi, & 
  jprsnd_lwd,               &
  jprsnd_t2m, jprsnd_q2m,   &
  jprsnd_pre, jprsnd_snw,   &
  jprsnd_slp, jprsnd_clf,   &
  jprsnd_lhf, jprsnd_shf,   &
  jprsnd_umoo, jprsnd_vmoo, &
  jprsnd_umoi, jprsnd_vmoi, &
  jprsnd_evp, jprsnd_sub,   &
  jprsnd_qnso, jprsnd_qnsi, &
  jprsnd_dqn, jprsnd_wnd,   &
  jprsnd_ros, jprsnd_rog
! for received fields
INTEGER ::                  &
  jprrcv_sst,               &
  jprrcv_ifr, jprrcv_ial,   &
  jprrcv_lhf, jprrcv_shf,   &
  jprrcv_umo, jprrcv_vmo,   &
  jprrcv_tice, jprrcv_hice


END MODULE cpl_oas_vardef

#endif
