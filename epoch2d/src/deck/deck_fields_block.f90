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

MODULE deck_fields_block

  USE strings_advanced
  USE simple_io

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: fields_deck_initialise, fields_deck_finalise
  PUBLIC :: fields_block_start, fields_block_end
  PUBLIC :: fields_block_handle_element, fields_block_check

  INTEGER(KIND=MPI_OFFSET_KIND) :: offset

CONTAINS

  SUBROUTINE fields_deck_initialise

  END SUBROUTINE fields_deck_initialise



  SUBROUTINE fields_deck_finalise

  END SUBROUTINE fields_deck_finalise



  SUBROUTINE fields_block_start

    offset = 0

  END SUBROUTINE fields_block_start



  SUBROUTINE fields_block_end

  END SUBROUTINE fields_block_end



  FUNCTION fields_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    CHARACTER(LEN=string_length) :: filename
    INTEGER :: err
    LOGICAL :: got_file

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'offset')) THEN
      offset = as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    CALL get_filename(value, filename, got_file, err)

    IF (str_cmp(element, 'ex')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ex, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_ex)
        CALL evaluate_string_in_space(value, ex, 1-ng, nx+ng, 1-ng, ny+ng, &
            errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'ey')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ey, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_ey)
        CALL evaluate_string_in_space(value, ey, 1-ng, nx+ng, 1-ng, ny+ng, &
            errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'ez')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ez, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_ez)
        CALL evaluate_string_in_space(value, ez, 1-ng, nx+ng, 1-ng, ny+ng, &
            errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'bx')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bx, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_bx)
        CALL evaluate_string_in_space(value, bx, 1-ng, nx+ng, 1-ng, ny+ng, &
            errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'by')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, by, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_by)
        CALL evaluate_string_in_space(value, by, 1-ng, nx+ng, 1-ng, ny+ng, &
            errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'bz')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bz, offset, errcode)
      ELSE
        CALL set_tokenizer_stagger(c_stagger_bz)
        CALL evaluate_string_in_space(value, bz, 1-ng, nx+ng, 1-ng, ny+ng, &
            errcode)
        CALL set_tokenizer_stagger(c_stagger_centre)
      END IF
      RETURN
    END IF

  END FUNCTION fields_block_handle_element



  FUNCTION fields_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION fields_block_check

END MODULE deck_fields_block
