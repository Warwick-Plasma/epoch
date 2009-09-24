MODULE initial_conditions

  USE shared_data
  USE strings
  USE shunt
  USE dist_fn
  USE probes
  USE iocontrol
  USE laser
  USE simple_io
  USE helper

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: IC_Early,IC_Late,ManualLoad


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE IC_Early

    INTEGER :: iSpecies
    REAL(num) :: target_den=5.0e24_num
    REAL(num) :: spotcentre=0.0_num
    REAL(num) :: spotsize=20.0e-6_num
    TYPE(Laser_Block),POINTER :: Laser



  END SUBROUTINE IC_Early

  SUBROUTINE IC_Late


  END SUBROUTINE IC_Late

  SUBROUTINE ManualLoad
    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies

  END SUBROUTINE ManualLoad

END MODULE initial_conditions
