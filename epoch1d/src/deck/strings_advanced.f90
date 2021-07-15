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

MODULE strings_advanced

  USE shunt
  USE evaluator

  IMPLICIT NONE

  PRIVATE :: get_allocated_array_rnum, get_allocated_array_int

  INTERFACE get_allocated_array
    MODULE PROCEDURE get_allocated_array_rnum, get_allocated_array_int
  END INTERFACE get_allocated_array

CONTAINS

  SUBROUTINE get_filename(str_in, str_out, got_file, err)

    CHARACTER(*), INTENT(IN) :: str_in
    CHARACTER(*), INTENT(OUT) :: str_out
    LOGICAL, INTENT(OUT) :: got_file
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: str_len, chr, pos, delimiter, c

    str_len = LEN(str_in)
    pos = -1
    got_file = .FALSE.

    c = ICHAR(str_in(1:1))
    ! ACHAR(39) = '\'', ACHAR(34) = '"'
    IF (c /= 39 .AND. c /= 34) RETURN

    delimiter = c

    DO chr = 2, str_len
      c = ICHAR(str_in(chr:chr))
      IF (c == delimiter) THEN
        pos = chr
        EXIT
      END IF
    END DO

    IF (pos < 0) THEN
      err = IOR(err, c_err_bad_value)
      RETURN
    END IF

    got_file = .TRUE.
    str_out = str_in(2:pos-1)

  END SUBROUTINE get_filename



  SUBROUTINE split_off_int(str_in, str_out, int_out, err)

    CHARACTER(*), INTENT(IN) :: str_in
    CHARACTER(*), INTENT(OUT) :: str_out
    INTEGER, INTENT(OUT) :: int_out
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: str_len, chr, pos, c

    str_len = LEN(str_in)
    pos = -1

    DO chr = 1, str_len
      c = ICHAR(str_in(chr:chr))
      ! ACHAR(48)-ACHAR(57) = '0'-'9'
      IF (c > 47 .AND. c < 58) THEN
        pos = chr
        EXIT
      END IF
    END DO

    IF (pos < 0) THEN
      err = IOR(err, c_err_bad_value)
      RETURN
    END IF

    str_out = str_in(1:pos-1)
    int_out = as_integer_simple(str_in(pos:str_len), err)

  END SUBROUTINE split_off_int



  SUBROUTINE split_range(str_in, real1, real2, err)

    CHARACTER(*), INTENT(IN) :: str_in
    REAL(num), INTENT(OUT) :: real1, real2
    INTEGER, INTENT(INOUT) :: err
    TYPE(primitive_stack) :: output
    REAL(num), DIMENSION(2) :: array

    CALL initialise_stack(output)
    CALL tokenize(str_in, output, err)
    IF (err == c_err_none) &
        CALL evaluate_to_array(output, 2, array, err)
    real1 = array(1)
    real2 = array(2)
    CALL deallocate_stack(output)

  END SUBROUTINE split_range



  SUBROUTINE get_vector(str_in, array, err)

    CHARACTER(*), INTENT(IN) :: str_in
    REAL(num), DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(INOUT) :: err
    TYPE(primitive_stack) :: output
    INTEGER :: ndim

    ndim = SIZE(array)
    CALL initialise_stack(output)
    CALL tokenize(str_in, output, err)
    IF (err == c_err_none) &
        CALL evaluate_to_array(output, ndim, array, err)
    CALL deallocate_stack(output)

  END SUBROUTINE get_vector



  SUBROUTINE get_allocated_array_rnum(str_in, array, err)

    CHARACTER(*), INTENT(IN) :: str_in
    REAL(num), POINTER :: array(:)
    INTEGER, INTENT(INOUT) :: err
    TYPE(primitive_stack) :: output
    INTEGER :: ndim
    CHARACTER(LEN=string_length) :: str_tmp

    CALL initialise_stack(output)
    str_tmp = '(' // TRIM(ADJUSTL(str_in)) // ')'

    CALL tokenize(str_tmp, output, err)
    IF (err == c_err_none) THEN
      NULLIFY(array)
      CALL evaluate_and_return_all(output, ndim, array, err)
    END IF
    CALL deallocate_stack(output)

  END SUBROUTINE get_allocated_array_rnum



  SUBROUTINE get_allocated_array_int(str_in, array, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, POINTER :: array(:)
    INTEGER, INTENT(INOUT) :: err
    REAL(num), POINTER :: rarray(:)

    CALL get_allocated_array_rnum(str_in, rarray, err)

    IF (ASSOCIATED(array)) DEALLOCATE(array)
    ALLOCATE(array(SIZE(rarray)))
    array = INT(rarray)
    DEALLOCATE(rarray)

  END SUBROUTINE get_allocated_array_int



  FUNCTION as_integer(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_integer
    as_integer = NINT(as_real(str_in, err))

  END FUNCTION as_integer



  FUNCTION as_long_integer(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER(i8) :: as_long_integer
    as_long_integer = NINT(as_real(str_in, err),i8)

  END FUNCTION as_long_integer



  FUNCTION as_real(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: as_real
    TYPE(primitive_stack) :: output

    CALL initialise_stack(output)
    CALL tokenize(str_in, output, err)
    as_real = evaluate(output, err)
    CALL deallocate_stack(output)

  END FUNCTION as_real



  SUBROUTINE as_list(str_in, array, n_elements, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(OUT) :: array(:)
    INTEGER, INTENT(OUT) :: n_elements
    INTEGER, INTENT(INOUT) :: err
    TYPE(primitive_stack) :: output

    CALL initialise_stack(output)
    CALL tokenize(str_in, output, err)
    IF (err == c_err_none) &
        CALL evaluate_as_list(output, array, n_elements)
    CALL deallocate_stack(output)

  END SUBROUTINE as_list



  SUBROUTINE evaluate_string_in_space(str_in, data_out, x1, x2, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER, INTENT(IN) :: x1, x2
    REAL(num), DIMENSION(:), INTENT(OUT) :: data_out
    TYPE(primitive_stack) :: output
    INTEGER :: ix
    TYPE(parameter_pack) :: parameters

    CALL initialise_stack(output)
    CALL tokenize(str_in, output, err)

    DO ix = x1, x2
      parameters%pack_ix = ix
      data_out(ix-x1+1) = evaluate_with_parameters(output, parameters, err)
    END DO
    CALL deallocate_stack(output)

  END SUBROUTINE evaluate_string_in_space



  FUNCTION as_integer_print(str_in, element, err) RESULT(res)

    CHARACTER(*), INTENT(IN) :: str_in, element
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: res

    res = as_integer(str_in, err)

    IF (.NOT.print_deck_constants .OR. rank /= 0) RETURN

    WRITE(du,'(A,I9)') TRIM(element) // ' = ', res

  END FUNCTION as_integer_print



  FUNCTION as_long_integer_print(str_in, element, err) RESULT(res)

    CHARACTER(*), INTENT(IN) :: str_in, element
    INTEGER, INTENT(INOUT) :: err
    INTEGER(i8) :: res

    res = as_long_integer(str_in, err)

    IF (.NOT.print_deck_constants .OR. rank /= 0) RETURN

    WRITE(du,'(A,I9)') TRIM(element) // ' = ', res

  END FUNCTION as_long_integer_print



  FUNCTION as_real_print(str_in, element, err) RESULT(res)

    CHARACTER(*), INTENT(IN) :: str_in, element
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: res

    res = as_real(str_in, err)

    IF (.NOT.print_deck_constants .OR. rank /= 0) RETURN

    WRITE(du,'(A,G18.11)') TRIM(element) // ' = ', res

  END FUNCTION as_real_print



  FUNCTION as_time(value, err)

    CHARACTER(LEN=*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: as_time

    IF (str_cmp(value, 'start')) THEN
      as_time = 0.0_num
      RETURN
    END IF

    IF (str_cmp(value, 'end')) THEN
      as_time = t_end
      RETURN
    END IF

    as_time = as_real(value, err)

  END FUNCTION as_time



  FUNCTION as_time_print(str_in, element, err) RESULT(res)

    CHARACTER(*), INTENT(IN) :: str_in, element
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: res

    res = as_time(str_in, err)

    IF (.NOT.print_deck_constants .OR. rank /= 0) RETURN

    WRITE(du,'(A,G18.11)') TRIM(element) // ' = ', res

  END FUNCTION as_time_print

END MODULE strings_advanced
