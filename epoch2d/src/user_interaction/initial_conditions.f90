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
  USE shunt

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


  END SUBROUTINE IC_Late

  SUBROUTINE ManualLoad

  END SUBROUTINE ManualLoad

END MODULE initial_conditions
