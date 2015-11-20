MODULE stack

  USE constants

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER, PRIVATE :: stack_size = 1024

  REAL(num) :: eval_stack_entries(stack_size)
  INTEGER :: eval_stack_flags(stack_size)
  INTEGER :: eval_stack_stack_point, eval_stack_nvalues

  INTEGER, PRIVATE :: flag

CONTAINS

  SUBROUTINE eval_reset()

    eval_stack_stack_point = 0
    eval_stack_nvalues = 0
    flag = 0

  END SUBROUTINE eval_reset



  SUBROUTINE push_on_eval(value)

    REAL(num), INTENT(IN) :: value

    eval_stack_stack_point = eval_stack_stack_point + 1
    eval_stack_entries(eval_stack_stack_point) = value
    eval_stack_flags(eval_stack_stack_point) = 0

  END SUBROUTINE push_on_eval



  SUBROUTINE push_eval_flag()

    eval_stack_stack_point = eval_stack_stack_point + 1
    eval_stack_flags(eval_stack_stack_point) = 1

  END SUBROUTINE push_eval_flag



  FUNCTION pop_off_eval() RESULT(value)

    REAL(num) :: value
    INTEGER :: sp

    sp = eval_stack_stack_point
    value = eval_stack_entries(sp)
    flag = flag + eval_stack_flags(sp)
    eval_stack_stack_point = sp - 1

  END FUNCTION pop_off_eval



  SUBROUTINE get_values(count, values)

    INTEGER, INTENT(IN) :: count
    REAL(num), DIMENSION(1:) :: values

    INTEGER :: i

    DO i = 1, count
      values(count-i+1) = pop_off_eval()
    ENDDO

    IF (flag /= 0) THEN
      eval_stack_nvalues = count
      eval_stack_stack_point = eval_stack_stack_point + count
      flag = 0
    ENDIF

  END SUBROUTINE get_values



  SUBROUTINE stack_point_fix()

    IF (eval_stack_nvalues > 0) THEN
      eval_stack_nvalues = eval_stack_nvalues + 1
      eval_stack_stack_point = eval_stack_stack_point + 1
    ENDIF

  END SUBROUTINE stack_point_fix

END MODULE stack
