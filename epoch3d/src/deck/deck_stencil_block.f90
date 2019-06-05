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

MODULE deck_stencil_block

  USE strings_advanced
  USE fields

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: stencil_deck_initialise, stencil_deck_finalise
  PUBLIC :: stencil_block_start, stencil_block_end
  PUBLIC :: stencil_block_handle_element, stencil_block_check

  LOGICAL :: dt_set

CONTAINS

  SUBROUTINE stencil_deck_initialise

    IF (deck_state == c_ds_first) THEN
      dt_set = .FALSE.
    END IF

  END SUBROUTINE stencil_deck_initialise



  SUBROUTINE stencil_deck_finalise

  END SUBROUTINE stencil_deck_finalise



  SUBROUTINE stencil_block_start

  END SUBROUTINE stencil_block_start



  SUBROUTINE stencil_block_end

  END SUBROUTINE stencil_block_end



  FUNCTION stencil_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'betaxy')) THEN
      betaxy = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'betayx')) THEN
      betayx = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'betaxz')) THEN
      betaxz = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'betazx')) THEN
      betazx = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'betayz')) THEN
      betayz = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'betazy')) THEN
      betazy = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'gammax')) THEN
      gammax = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'gammay')) THEN
      gammay = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'gammaz')) THEN
      gammaz = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'deltax')) THEN
      deltax = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'deltay')) THEN
      deltay = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'deltaz')) THEN
      deltaz = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'dt')) THEN
      dt_custom = as_real_print(value, element, errcode)
      dt_set = .TRUE.
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION stencil_block_handle_element



  FUNCTION stencil_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none

    IF (maxwell_solver /= c_maxwell_solver_custom) RETURN

    IF (.NOT. dt_set) THEN
      errcode = c_err_missing_elements
    END IF

  END FUNCTION stencil_block_check

END MODULE deck_stencil_block
