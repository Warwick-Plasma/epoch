MODULE ic_module

  USE shared_data
  USE strings
  USE shunt
  USE dist_fn
  USE probes
  USE iocontrol
  USE laser
  USE helper

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ic_early, ic_late, manual_load

CONTAINS

  !----------------------------------------------------------------------------
  ! This function contains the equilibrium
  !----------------------------------------------------------------------------

  SUBROUTINE ic_early

  END SUBROUTINE ic_early



  SUBROUTINE ic_late

  END SUBROUTINE ic_late



  SUBROUTINE manual_load

  END SUBROUTINE manual_load

END MODULE ic_module
