! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE custom_parser

  USE shared_data
  USE tokenizer_blocks
  USE stack

  IMPLICIT NONE

CONTAINS

  ! This function is used to register custom functions and constants
  ! Use the register_function and register_constant functions
  SUBROUTINE register_objects

  END SUBROUTINE register_objects



  !----------------------------------------------------------------------------
  ! These functions contain the user functions for the deck parser
  !----------------------------------------------------------------------------

  FUNCTION custom_function(opcode, parameters, errcode)

    INTEGER, INTENT(IN) :: opcode
    TYPE(parameter_pack), INTENT(IN) :: parameters
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

  FUNCTION custom_constant(opcode, parameters, errcode)

    INTEGER, INTENT(IN) :: opcode
    TYPE(parameter_pack), INTENT(IN) :: parameters
    INTEGER, INTENT(INOUT) :: errcode
    REAL(num) :: custom_constant

    ! Leave these lines in place. They cause the code to throw an error if
    ! The opcode is unknown
    custom_constant = 0.0_num
    errcode = IOR(errcode, c_err_unknown_element)

  END FUNCTION custom_constant

END MODULE custom_parser
