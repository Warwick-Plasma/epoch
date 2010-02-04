MODULE custom_parser

  USE shared_data
  USE tokenizer_blocks
  IMPLICIT NONE

CONTAINS

  ! This function is used to register custom functions and constants
  ! Use the register_function and register_constant functions
  SUBROUTINE register_objects

  END SUBROUTINE register_objects



  !----------------------------------------------------------------------------
  ! These functions contain the user functions for the deck parser
  !----------------------------------------------------------------------------

  FUNCTION custom_function(opcode, ix, iy, iz, errcode)

    INTEGER, INTENT(IN) :: opcode, ix, iy, iz
    INTEGER, INTENT(INOUT) :: errcode
    REAL(num) :: custom_function

    ! Leave these lines in place. They cause the code to throw an error if
    ! The opcode is unknown
    custom_function = 0.0_num
    errcode = IOR(errcode, c_err_unknown_element)

  END FUNCTION custom_function



  !----------------------------------------------------------------------------
  ! These functions contain the user constants for the deck parser
  !----------------------------------------------------------------------------

  FUNCTION custom_constant(opcode, ix, iy, iz, errcode)

    INTEGER, INTENT(IN) :: opcode, ix, iy, iz
    INTEGER, INTENT(INOUT) :: errcode
    REAL(num) :: custom_constant

    ! Leave these lines in place. They cause the code to throw an error if
    ! The opcode is unknown
    custom_constant = 0.0_num
    errcode = IOR(errcode, c_err_unknown_element)

  END FUNCTION custom_constant

END MODULE custom_parser
