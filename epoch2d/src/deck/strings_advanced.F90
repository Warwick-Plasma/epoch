MODULE strings_advanced


  USE shunt
  USE evaluator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE split_off_int(str_in, str_out, int_out, err)

    CHARACTER(*), INTENT(IN) :: str_in
    CHARACTER(*), INTENT(OUT) :: str_out
    INTEGER, INTENT(OUT) :: int_out
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: str_len, char, pos, c

    str_len = LEN(str_in)
    pos = -1

    DO char = 1, str_len
      c = ICHAR(str_in(char:char))
      IF (c .GT. 47 .AND. c .LT. 58) THEN
        pos = char
        EXIT
      ENDIF
    ENDDO

    IF (pos .LT. 0) THEN
      err = IOR(err, c_err_bad_value)
      RETURN
    ENDIF

    str_out = str_in(1:pos-1)
    int_out = as_integer_simple(str_in(pos:str_len), err)

  END SUBROUTINE split_off_int



  SUBROUTINE split_range(str_in, real1, real2, err)

    CHARACTER(*), INTENT(IN) :: str_in
    REAL(num), INTENT(OUT) :: real1, real2
    INTEGER, INTENT(INOUT) :: err
    TYPE(primitive_stack) :: output
    REAL(num), DIMENSION(2) :: array

    output%stack_point = 0
    CALL tokenize(str_in, output, err)
    IF (err .EQ. c_err_none) &
        CALL evaluate_at_point_to_array(output, 0, 0, 2, array, err)
    real1 = array(1)
    real2 = array(2)

  END SUBROUTINE split_range



  FUNCTION as_integer(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_integer
    as_integer = NINT(as_real(str_in, err))

  END FUNCTION as_integer



  FUNCTION as_long_integer(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER(KIND=8) :: as_long_integer
    as_long_integer = NINT(as_real(str_in, err))

  END FUNCTION as_long_integer



  FUNCTION as_real(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: as_real
    TYPE(primitive_stack) :: output

    output%stack_point = 0
    CALL tokenize(str_in, output, err)
    as_real = evaluate(output, err)

  END FUNCTION as_real



  SUBROUTINE evaluate_string_in_space(str_in, data_out, xrange, yrange, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER, DIMENSION(2), INTENT(IN) :: xrange, yrange
    REAL(num), DIMENSION(1:,1:), INTENT(OUT) :: data_out
    TYPE(primitive_stack) :: output
    INTEGER :: ix, iy

    output%stack_point = 0
    CALL tokenize(str_in, output, err)

    DO iy = yrange(1), yrange(2)
      DO ix = xrange(1), xrange(2)
        data_out(ix-xrange(1)+1, iy-yrange(1)+1) = &
            evaluate_at_point(output, ix, iy, err)
      ENDDO
    ENDDO

  END SUBROUTINE evaluate_string_in_space

END MODULE strings_advanced
