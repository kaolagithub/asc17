!###############################################################################
!-------------------------------------------------------------------------------

  program masnum_wam_mpi

!-------------------------------------------------------------------------------

  use wammpi_mod
  
  implicit none

!-------------------------------------------------------------------------------

  call wammpi_init
  call readwi_mpi
  call ympi_final
  
!-------------------------------------------------------------------------------

  end program masnum_wam_mpi

!-------------------------------------------------------------------------------
!###############################################################################
