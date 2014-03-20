MODULE stack

  USE constants

  IMPLICIT NONE

  TYPE float_stack
    REAL(num) :: entries(1024)
    INTEGER :: flags(1024)
    INTEGER :: stack_point, nvalues
  END TYPE float_stack

  TYPE(float_stack), SAVE :: eval_stack
  INTEGER, PRIVATE :: flag

CONTAINS

  SUBROUTINE eval_reset()

    eval_stack%stack_point = 0
    eval_stack%nvalues = 0
    flag = 0

  END SUBROUTINE eval_reset



  SUBROUTINE push_on_eval(value)

    REAL(num), INTENT(IN) :: value

    eval_stack%stack_point = eval_stack%stack_point + 1
    eval_stack%entries(eval_stack%stack_point) = value
    eval_stack%flags(eval_stack%stack_point) = 0

  END SUBROUTINE push_on_eval



  SUBROUTINE push_eval_flag()

    eval_stack%stack_point = eval_stack%stack_point + 1
    eval_stack%flags(eval_stack%stack_point) = 1

  END SUBROUTINE push_eval_flag



  FUNCTION pop_off_eval() RESULT(value)

    REAL(num) :: value
    INTEGER :: sp

    sp = eval_stack%stack_point
    value = eval_stack%entries(sp)
    flag = flag + eval_stack%flags(sp)
    eval_stack%stack_point = sp - 1

  END FUNCTION pop_off_eval



  SUBROUTINE pop_off_eval_always(value, flag_set)

    REAL(num), INTENT(OUT) :: value
    LOGICAL, INTENT(OUT) :: flag_set
    INTEGER :: sp

    sp = eval_stack%stack_point
    value = eval_stack%entries(sp)
    eval_stack%stack_point = sp - 1

    IF (eval_stack%flags(sp) .EQ. 0) THEN
      flag_set = .FALSE.
    ELSE
      flag_set = .TRUE.
    ENDIF

  END SUBROUTINE pop_off_eval_always



  SUBROUTINE get_values(count, values)

    INTEGER, INTENT(IN) :: count
    REAL(num), DIMENSION(1:) :: values

    INTEGER :: i

    DO i = 1, count
      values(count-i+1) = pop_off_eval()
    ENDDO

    IF (flag .NE. 0) THEN
      eval_stack%nvalues = count
      eval_stack%stack_point = eval_stack%stack_point + count
      flag = 0
    ENDIF

  END SUBROUTINE get_values



  SUBROUTINE stack_point_fix()

    IF (eval_stack%nvalues .GT. 0) THEN
      eval_stack%nvalues = eval_stack%nvalues + 1
      eval_stack%stack_point = eval_stack%stack_point + 1
    ENDIF

  END SUBROUTINE stack_point_fix

END MODULE stack
