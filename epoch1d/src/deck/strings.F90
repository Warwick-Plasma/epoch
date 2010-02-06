MODULE strings

  USE shared_data

  IMPLICIT NONE

CONTAINS

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
    INTEGER :: as_integer_simple, value, f

    READ(unit=str_in, fmt=*, iostat=f) value
    IF (f .NE. 0) err = IOR(err, c_err_bad_value)
    as_integer_simple = value

  END FUNCTION as_integer_simple



  FUNCTION as_long_integer_simple(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER(KIND=8) :: as_long_integer_simple, value, f
    READ(unit=str_in, fmt=*, iostat=f) value
    IF (f .NE. 0) err = IOR(err, c_err_bad_value)
    as_long_integer_simple = value

  END FUNCTION as_long_integer_simple



  FUNCTION as_direction(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: as_direction

    as_direction = -1

    IF (str_cmp(str_in, "left")) as_direction = c_bd_left
    IF (str_cmp(str_in, "right")) as_direction = c_bd_right
    IF (str_cmp(str_in, "up")) as_direction = c_bd_up
    IF (str_cmp(str_in, "down")) as_direction = c_bd_down

    IF (as_direction .EQ. -1) err = IOR(err, c_err_bad_value)

  END FUNCTION as_direction



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
