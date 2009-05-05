MODULE Custom_parser
  USE shared_data
  USE tokenizer_blocks
  USE stack
  IMPLICIT NONE

CONTAINS

  !This function is used to register custom functions and constants
  !Use the RegisterFunction and RegisterConstant functions
  SUBROUTINE RegisterObjects

  END SUBROUTINE RegisterObjects


  !-----------------------------------------------------------------------------
  !These functions contain the user functions for the deck parser
  !-----------------------------------------------------------------------------
  FUNCTION CustomFunction(opcode,ix,iy,errcode)
    INTEGER, INTENT(IN) :: opcode,ix,iy
    INTEGER, INTENT(INOUT) :: errcode
    REAL(num) :: CustomFunction
    REAL(num) :: Values(5)

    !Leave these lines in place. They cause the code to throw an error if
    !The opcode is unknown
    CustomFunction=0.0_num
    errcode=IOR(errcode,ERR_UNKNOWN_ELEMENT)

  END FUNCTION CustomFunction

  !-----------------------------------------------------------------------------
  !These functions contain the user constants for the deck parser
  !-----------------------------------------------------------------------------
  FUNCTION CustomConstant(opcode,ix,iy,errcode)
    INTEGER, INTENT(IN) :: opcode,ix,iy
    INTEGER, INTENT(INOUT) :: errcode
    REAL(num) :: CustomConstant

    !Leave these lines in place. They cause the code to throw an error if
    !The opcode is unknown
    CustomConstant=0.0_num
    errcode=IOR(errcode,ERR_UNKNOWN_ELEMENT)

  END FUNCTION CustomConstant

END MODULE Custom_parser
