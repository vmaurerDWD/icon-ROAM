#ifdef COUP_OASIS3MCT

MODULE cpl_oas_mpi

!     AUTHOR.
!     -------
!     2020-11: Ha Ho-Hagemann (Hereon): OASIS interface for ICON-CLM

USE mo_kind,           ONLY: wp
USE mo_impl_constants, ONLY: modelname  ! 'icon'
USE cpl_oas_vardef

IMPLICIT NONE

SAVE
CHARACTER(len=4)      :: cpl_comp_name

!==============================================================================
PUBLIC :: cpl_oas_init,           &
          cpl_oas_finalize

PUBLIC :: cpl_comp_name

!==============================================================================

CONTAINS

!==============================================================================
SUBROUTINE cpl_oas_init

#if defined COUP_OASIS3MCT
  call oasis_atm_init
  cpl_comp_name = oas_comp_name
#endif
END SUBROUTINE cpl_oas_init

!==============================================================================
SUBROUTINE cpl_oas_finalize

#if defined COUP_OASIS3MCT
  call oasis_atm_finalize
#endif
END SUBROUTINE cpl_oas_finalize

!==============================================================================

SUBROUTINE oasis_atm_init

!**** *INIOASIS*  - Initialize coupled mode communication
!
!     Purpose.
!     --------
!     Initialize coupler to get the MPI communicator
!
!**   Interface.
!     ----------
!       *CALL*  *oasis_atm_init*
!
!     Input:
!     -----
!
!     Output:
!     ------
!       MPLUSERCOMM in MPL_MODULE is updated.
!
!     Method:
!     ------
!       OASIS usage is controlled by environment variables
!
!     Externals:
!     ---------
!       GETENV - Get enviroment variables
!       oasis_init, oasis_init_comp, oasis_get_localcomm, oasis_abort : oasis library
!
!     Reference:
!     ---------
!       S. Valcke, R. Redler, 2007: OASIS4 User Guide ,
!       OASIS Support Initiative Report No 4,
!       CERFACS, Toulouse, France, 60 pp.
!
!     Author:
!     -------
!       R. Redler, NEC Laboratories Europe
!       K. Mogensen, ECMWF
!
!     Modifications.
!     --------------
!       E. Maisonnave : adapted to COSMO
!       A. Dobler: adapted to OASIS3
!       Ha Ho-Hagemann: adapted to OASIS interface for ICON-CLM

INTEGER ::  rank, peer_comm, p_error

!------------------------------------------------------------------
! 1) Initialize the OASIS system for the component
!------------------------------------------------------------------
oas_comp_name = modelname

CALL oasis_init_comp( oas_comp_id, oas_comp_name, oas_error )
IF( oas_error /= OASIS_Success ) THEN
  CALL oasis_abort( oas_comp_id, 'oasis_atm_init', 'Failure in oasis_init_comp' )
ENDIF

!------------------------------------------------------------------
! 2) Get an MPI communicator for ICON-CLM local communication
!------------------------------------------------------------------

CALL oasis_get_localcomm( kl_comm, oas_error )
IF( oas_error /= OASIS_Success ) THEN
  CALL oasis_abort( oas_comp_id, 'oasis_atm_init', 'Failure in oasis_get_localcomm' )
ENDIF

END SUBROUTINE oasis_atm_init
!==============================================================================

SUBROUTINE oasis_atm_finalize

!!---------------------------------------------------------------------
!!              ***  ROUTINE oasis_atm_finalize  ***
!!
!! ** Purpose : - Finalizes the coupling. If MPI_init has not been
!!      called explicitly before oasis_atm_init it will also close
!!      MPI communication.
!!----------------------------------------------------------------------

CALL oasis_terminate ( oas_error )         
IF ( debug_oasis > 1 ) print*, 'Call oasis_terminate in ICON-CLM'

END SUBROUTINE oasis_atm_finalize
!==============================================================================

END MODULE cpl_oas_mpi

#endif
