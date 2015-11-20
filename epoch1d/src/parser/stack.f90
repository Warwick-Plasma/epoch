MODULE stack

  USE constants

  IMPLICIT NONE

  SAVE

  INTEGER, PRIVATE :: stack_size

  REAL(num), ALLOCATABLE :: eval_stack_entries(:)
  INTEGER, ALLOCATABLE :: eval_stack_flags(:)
  INTEGER :: eval_stack_stack_point, eval_stack_nvalues

  INTEGER, PRIVATE :: flag
  PRIVATE :: eval_stack_grow

CONTAINS

  SUBROUTINE eval_stack_init()

    stack_size = 64
    ALLOCATE(eval_stack_entries(stack_size))
    ALLOCATE(eval_stack_flags(stack_size))

  END SUBROUTINE eval_stack_init



  SUBROUTINE deallocate_eval_stack()

    stack_size = 0
    DEALLOCATE(eval_stack_entries)
    DEALLOCATE(eval_stack_flags)

  END SUBROUTINE deallocate_eval_stack



  SUBROUTINE eval_stack_grow()

    REAL(num), ALLOCATABLE :: eval_stack_entries_tmp(:)
    INTEGER, ALLOCATABLE :: eval_stack_flags_tmp(:)
    INTEGER :: old_size

    ALLOCATE(eval_stack_entries_tmp(stack_size))
    ALLOCATE(eval_stack_flags_tmp(stack_size))

    eval_stack_entries_tmp = eval_stack_entries
    eval_stack_flags_tmp = eval_stack_flags

    old_size = stack_size
    stack_size = 2 * stack_size
    DEALLOCATE(eval_stack_entries)
    DEALLOCATE(eval_stack_flags)
    ALLOCATE(eval_stack_entries(stack_size))
    ALLOCATE(eval_stack_flags(stack_size))

    eval_stack_entries(1:old_size) = eval_stack_entries_tmp
    eval_stack_flags(1:old_size) = eval_stack_flags_tmp

    DEALLOCATE(eval_stack_entries_tmp)
    DEALLOCATE(eval_stack_flags_tmp)

  END SUBROUTINE eval_stack_grow



  SUBROUTINE eval_reset()

    eval_stack_stack_point = 0
    eval_stack_nvalues = 0
    flag = 0

  END SUBROUTINE eval_reset



  SUBROUTINE push_on_eval(value)

    REAL(num), INTENT(IN) :: value

    eval_stack_stack_point = eval_stack_stack_point + 1
    IF (eval_stack_stack_point > stack_size) CALL eval_stack_grow
    eval_stack_entries(eval_stack_stack_point) = value
    eval_stack_flags(eval_stack_stack_point) = 0

  END SUBROUTINE push_on_eval



  SUBROUTINE push_eval_flag()

    eval_stack_stack_point = eval_stack_stack_point + 1
    IF (eval_stack_stack_point > stack_size) CALL eval_stack_grow
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
