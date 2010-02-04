MODULE strings

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE integer_as_string(INTEGER,string)

    INTEGER,INTENT(IN) :: INTEGER
    CHARACTER(len=*),INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(len=9) :: numfmt

    n_nums=1 + LOG10(REAL(INTEGER,num))
    WRITE(numfmt,'("(I",I6.6,")")') n_nums
    WRITE(string,numfmt) INTEGER

  END SUBROUTINE integer_as_string

  SUBROUTINE integer8_as_string(INTEGER,string)

    INTEGER(KIND=8),INTENT(IN) :: INTEGER
    CHARACTER(len=*),INTENT(OUT) :: string

    INTEGER :: n_nums
    CHARACTER(len=12) :: numfmt

    n_nums=1 + LOG10(REAL(INTEGER,num))
    WRITE(numfmt,'("(I",I9.9,")")') n_nums
    WRITE(string,numfmt) INTEGER

  END SUBROUTINE integer8_as_string

  FUNCTION str_cmp(str_in,str_test)

    CHARACTER(*),INTENT(IN) ::  str_in,str_test
    CHARACTER(len=string_length) :: str_trim
    LOGICAL :: str_cmp

    str_trim=TRIM(ADJUSTL(str_in))

    IF (len(str_test) .GT. len(str_trim)) THEN
       str_cmp=.FALSE.
       RETURN
    ENDIF

    IF (len(str_test) .LT. string_length) THEN
       IF (str_trim(len(str_test)+1:len(str_test)+1) .NE. " ") THEN
          str_cmp=.FALSE.
          RETURN
       ENDIF
    ENDIF

    str_cmp=(str_trim(1:len(str_test)) == str_test)

  END FUNCTION str_cmp

  FUNCTION as_real_simple(str_in,err)
    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER,INTENT(INOUT) :: err
    INTEGER::f
    REAL(num) :: as_real_simple
    REAL(num) :: value

    READ(unit=str_in,fmt=*,iostat=f) value
    IF (f .NE. 0) err=IOR(err,ERR_BAD_VALUE)
    as_real_simple=value
  END FUNCTION as_real_simple

  FUNCTION as_integer_simple(str_in,err)
    CHARACTER(*),INTENT(IN) :: str_in
    INTEGER,INTENT(INOUT) :: err
    INTEGER :: as_integer_simple,value,f
    READ(unit=str_in,fmt=*,iostat=f) value
    IF (f .NE. 0) err=IOR(err,ERR_BAD_VALUE)
    as_integer_simple=value

  END FUNCTION as_integer_simple

  FUNCTION as_direction(str_in,err)
    CHARACTER(*),INTENT(IN) :: str_in
    INTEGER,INTENT(INOUT) :: err
    INTEGER :: as_direction
    as_direction=-1

    IF (str_cmp(str_in,"left")) as_direction=BD_LEFT
    IF (str_cmp(str_in,"right")) as_direction=BD_RIGHT
    IF (str_cmp(str_in,"up")) as_direction=BD_UP
    IF (str_cmp(str_in,"down")) as_direction=BD_DOWN

    IF (as_direction == -1) err=IOR(err,ERR_BAD_VALUE)

  END FUNCTION as_direction

  FUNCTION as_logical(str_in,err)
    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER,INTENT(INOUT) :: err
    LOGICAL :: as_logical

    as_logical=.FALSE.
    IF(str_cmp(TRIM(ADJUSTL(str_in)),"T")) THEN
       as_logical=.TRUE.
       RETURN
    ENDIF
    IF(str_cmp(TRIM(ADJUSTL(str_in)),"F")) THEN
       as_logical=.FALSE.
       RETURN
    ENDIF

    err=IOR(err,ERR_BAD_VALUE)

  END FUNCTION as_logical

  FUNCTION as_bc(str_in,err)
    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER,INTENT(INOUT) :: err
    INTEGER :: as_bc

    as_bc=-1

    IF (str_cmp(TRIM(ADJUSTL(str_in)),"periodic")) THEN
       as_bc=BC_PERIODIC
       RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)),"simple_laser")) THEN
       as_bc=BC_SIMPLE_LASER
       RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)),"simple_outflow")) THEN
       as_bc=BC_SIMPLE_OUTFLOW
       RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)),"other")) THEN
       as_bc=BC_OTHER
       RETURN
    ENDIF

    err=IOR(err,ERR_BAD_VALUE)

  END FUNCTION as_bc

  FUNCTION as_dump_param(str_in,err)
    CHARACTER(*),INTENT(IN) :: str_in
    INTEGER,INTENT(INOUT) :: err
    INTEGER :: as_dump_param

    as_dump_param=-1

    IF (str_cmp(TRIM(ADJUSTL(str_in)),"never")) THEN
       as_dump_param=IO_NEVER
       RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)),"always")) THEN
       as_dump_param=IO_ALWAYS
       RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)),"full")) THEN
       as_dump_param=IO_FULL
       RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)),"restart")) THEN
       as_dump_param=IO_RESTARTABLE
       RETURN
    ENDIF

    err=IOR(err,ERR_BAD_VALUE)

  END FUNCTION as_dump_param

  FUNCTION as_domain(str_in,err)
    CHARACTER(*),INTENT(IN) :: str_in
    INTEGER,INTENT(INOUT) :: err
    INTEGER :: as_domain

    as_domain=-1

    IF (str_cmp(TRIM(ADJUSTL(str_in)),"decomposed")) THEN
       as_domain=DO_DECOMPOSED
       RETURN
    ENDIF
    IF (str_cmp(TRIM(ADJUSTL(str_in)),"full")) THEN
       as_domain=DO_FULL
       RETURN
    ENDIF

    err=IOR(err,ERR_BAD_VALUE)
  END FUNCTION as_domain


END MODULE strings
