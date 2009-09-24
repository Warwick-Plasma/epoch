MODULE initial_conditions

  USE shared_data
  USE strings
  USE shunt
  USE dist_fn
  USE iocontrol

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: IC_Early,IC_Late,ManualLoad


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE IC_Early


  END SUBROUTINE IC_Early

  SUBROUTINE IC_Late

    INTEGER :: iSpecies


  END SUBROUTINE IC_Late

  SUBROUTINE ManualLoad
    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies


  END SUBROUTINE ManualLoad

END MODULE initial_conditions
