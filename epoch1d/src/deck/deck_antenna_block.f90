! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2018 Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE deck_antenna_block

  USE evaluator
  USE shunt
  USE antennae
  USE strings_advanced

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: antenna_deck_initialise, antenna_deck_finalise
  PUBLIC :: antenna_block_start, antenna_block_end
  PUBLIC :: antenna_block_handle_element, antenna_block_check

  TYPE(antenna), POINTER :: working_antenna
  LOGICAL :: antenna_ok

CONTAINS

  SUBROUTINE antenna_deck_initialise

  END SUBROUTINE antenna_deck_initialise



  SUBROUTINE antenna_deck_finalise

  END SUBROUTINE antenna_deck_finalise



  SUBROUTINE antenna_block_start

    IF (deck_state == c_ds_first) RETURN
    ALLOCATE(working_antenna)
    antenna_ok = .TRUE.
    CALL initialise_antenna(working_antenna)

  END SUBROUTINE antenna_block_start



  SUBROUTINE antenna_block_end

    IF (deck_state == c_ds_first) RETURN
    antenna_ok = .TRUE.
    CALL add_antenna(working_antenna)

  END SUBROUTINE antenna_block_end



  FUNCTION antenna_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(LEN=*), INTENT(IN) :: element, value
    REAL(num), DIMENSION(:), POINTER :: ranges
    INTEGER :: errcode, nels

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'jx')) THEN
      CALL initialise_stack(working_antenna%jx_expression)
      CALL tokenize(value, working_antenna%jx_expression, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'jy')) THEN
      CALL initialise_stack(working_antenna%jy_expression)
      CALL tokenize(value, working_antenna%jy_expression, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'jz')) THEN
      CALL initialise_stack(working_antenna%jz_expression)
      CALL tokenize(value, working_antenna%jz_expression, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'ranges')) THEN
      CALL initialise_stack(working_antenna%ranges)
      CALL tokenize(value, working_antenna%ranges, errcode)
      ranges => NULL()
      CALL evaluate_and_return_all(working_antenna%ranges, nels, ranges, &
          errcode)
      DEALLOCATE(ranges)
      IF (errcode /= c_err_none .OR. nels /= c_ndims * 2) THEN
        errcode = IOR(errcode, c_err_bad_value)
        antenna_ok = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'omega')) THEN
      CALL initialise_stack(working_antenna%omega)
      CALL tokenize(value, working_antenna%omega, errcode)
      IF (.NOT. working_antenna%omega%is_time_varying) THEN
        working_antenna%omega_value = as_real_print(value, element, errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'start_time')) THEN
      working_antenna%start_time = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'stop_time')) THEN
      working_antenna%stop_time = as_time_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION antenna_block_handle_element



  FUNCTION antenna_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none

  END FUNCTION antenna_block_check

END MODULE deck_antenna_block
