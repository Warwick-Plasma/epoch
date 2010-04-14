MODULE strings

  USE shared_data

  IMPLICIT NONE

  PRIVATE :: integer4_as_string, integer8_as_string

  INTERFACE integer_as_string
    MODULE PROCEDURE integer4_as_string, integer8_as_string
  END INTERFACE integer_as_string

CONTAINS

  SUBROUTINE integer4_as_string(int_in, string)

    INTEGER(KIND=4), INTENT(IN) :: int_in
    CHARACTER(LEN=*), INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(LEN=9) :: numfmt

    IF (int_in .EQ. 0) THEN
      n_nums = 1
    ELSE
      n_nums = 1 + INT(LOG10(REAL(ABS(int_in), num)))
    ENDIF
    WRITE(numfmt, '("(I", I6.6, ")")') n_nums
    WRITE(string, numfmt) int_in

  END SUBROUTINE integer4_as_string



  SUBROUTINE integer8_as_string(int_in, string)

    INTEGER(KIND=8), INTENT(IN) :: int_in
    CHARACTER(LEN=*), INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(LEN=12) :: numfmt

    IF (int_in .EQ. 0) THEN
      n_nums = 1
    ELSE
      n_nums = 1 + INT(LOG10(REAL(ABS(int_in), num)))
    ENDIF
    WRITE(numfmt, '("(I", I9.9, ")")') n_nums
    WRITE(string, numfmt) int_in

  END SUBROUTINE integer8_as_string



  FUNCTION str_cmp(str_in, str_test)

    CHARACTER(*), INTENT(IN) ::  str_in, str_test
    CHARACTER(LEN=string_length) :: str_trim
    LOGICAL :: str_cmp

    str_trim = TRIM(ADJUSTL(str_in))

    IF (LEN(str_test) .GT. LEN(str_trim)) THEN
      str_cmp = .FALSE.
      RETURN
    ENDIF

    IF (LEN(str_test) .LT. string_length) THEN
      IF (str_trim(LEN(str_test)+1:LEN(str_test)+1) .NE. " ") THEN
        str_cmp = .FALSE.
        RETURN
      ENDIF
    ENDIF

    str_cmp = (str_trim(1:LEN(str_test)) .EQ. str_test)

  END FUNCTION str_cmp



  FUNCTION as_real_simple(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: f
    REAL(num) :: as_real_simple
    REAL(num) :: value

    READ(unit=str_in, fmt=*, iostat=f) value
    IF (f .NE. 0) err = IOR(err, c_err_bad_value)
    as_real_simple = value

  END FUNCTION as_real_simple



  FUNCTION as_integer_simple(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_integer_simple, value
    INTEGER :: f

    READ(unit=str_in, fmt=*, iostat=f) value
    IF (f .NE. 0) err = IOR(err, c_err_bad_value)
    as_integer_simple = value

  END FUNCTION as_integer_simple



  FUNCTION as_long_integer_simple(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER(KIND=8) :: as_long_integer_simple, value
    INTEGER :: f

    READ(unit=str_in, fmt=*, iostat=f) value
    IF (f .NE. 0) err = IOR(err, c_err_bad_value)
    as_long_integer_simple = value

  END FUNCTION as_long_integer_simple



  FUNCTION as_boundary(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_boundary

    as_boundary = -1

    IF (str_cmp(str_in, "left") &
        .OR. str_cmp(str_in, "x_min")) as_boundary = c_bd_x_min
    IF (str_cmp(str_in, "right") &
        .OR. str_cmp(str_in, "x_max")) as_boundary = c_bd_x_max
    IF (str_cmp(str_in, "down") &
        .OR. str_cmp(str_in, "y_min")) as_boundary = c_bd_y_min
    IF (str_cmp(str_in, "up") &
        .OR. str_cmp(str_in, "y_max")) as_boundary = c_bd_y_max

    IF (as_boundary .EQ. -1) err = IOR(err, c_err_bad_value)

  END FUNCTION as_boundary



  FUNCTION as_logical(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    LOGICAL :: as_logical

    as_logical = .FALSE.
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "T")) THEN
      as_logical = .TRUE.
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "F")) THEN
      as_logical = .FALSE.
      RETURN
    ENDIF

    err = IOR(err, c_err_bad_value)

  END FUNCTION as_logical



  FUNCTION as_bc(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_bc

    as_bc = -1

    IF (str_cmp(TRIM(ADJUSTL(str_in)), "periodic")) THEN
      as_bc = c_bc_periodic
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "simple_laser")) THEN
      as_bc = c_bc_simple_laser
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "simple_outflow")) THEN
      as_bc = c_bc_simple_outflow
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "other")) THEN
      as_bc = c_bc_other
      RETURN
    ENDIF

    err = IOR(err, c_err_bad_value)

  END FUNCTION as_bc



  FUNCTION as_dump_param(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_dump_param

    as_dump_param = -1

    IF (str_cmp(TRIM(ADJUSTL(str_in)), "never")) THEN
      as_dump_param = c_io_never
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "always")) THEN
      as_dump_param = c_io_always
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "full")) THEN
      as_dump_param = c_io_full
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "restart")) THEN
      as_dump_param = c_io_restartable
      RETURN
    ENDIF

    err = IOR(err, c_err_bad_value)

  END FUNCTION as_dump_param



  FUNCTION as_domain(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_domain

    as_domain = -1

    IF (str_cmp(TRIM(ADJUSTL(str_in)), "decomposed")) THEN
      as_domain = c_do_decomposed
      RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)), "full")) THEN
      as_domain = c_do_full
      RETURN
    ENDIF

    err = IOR(err, c_err_bad_value)

  END FUNCTION as_domain

END MODULE strings
