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

MODULE custom_deck

  USE constants

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! These functions contain the user input deck elements
  !----------------------------------------------------------------------------

  FUNCTION custom_blocks_handle_element(block_name, element, value) &
      RESULT(errcode)

    CHARACTER(LEN=string_length), INTENT(IN) :: block_name, element, value
    INTEGER :: errcode

    ! The following line must always be present
    errcode = c_err_unknown_block

  END FUNCTION custom_blocks_handle_element



  FUNCTION custom_blocks_check() RESULT(errcode)

    INTEGER :: errcode

    ! This subroutine is to allow you to force the code to bomb out if an
    ! essential element of the input deck is missing. If you either don't
    ! want to check, are not extending the input deck, or all elements are
    ! set then set 'errcode = c_err_none'. Otherwise set the
    ! return value to 'errcode = c_err_missing_elements'.

    errcode = c_err_none

  END FUNCTION custom_blocks_check

END MODULE custom_deck
