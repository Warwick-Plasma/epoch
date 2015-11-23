MODULE evaluator

  USE evaluator_blocks
  USE shunt
  USE utilities

  IMPLICIT NONE

  TYPE stack_list
    TYPE(primitive_stack) :: stack
    TYPE(stack_list), POINTER :: prev, next
  END TYPE stack_list

  TYPE(stack_list), POINTER :: sl_head, sl_tail
  INTEGER :: sl_size

CONTAINS

#ifdef SIMPLIFY_DEBUG
  SUBROUTINE basic_evaluate(input_stack, ix, iy, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i, n_elements
    REAL(num), ALLOCATABLE :: array(:)

    IF (input_stack%should_simplify) THEN
      CALL basic_evaluate_standard(input_stack, ix, iy, err)

      n_elements = eval_stack_stack_point
      ALLOCATE(array(1:n_elements))

      ! Pop off the final answers
      DO i = n_elements,1,-1
        array(i) = pop_off_eval()
      ENDDO

      CALL simplify_stack(input_stack, err)

      CALL basic_evaluate_standard(input_stack, ix, iy, err)

      ! Check the final answers
      DO i = 1, n_elements
        IF (eval_stack_entries(i) /= array(i)) THEN
          PRINT*,i,eval_stack_entries(i),array(i),eval_stack_entries(i)-array(i)
        ENDIF
      ENDDO

      DEALLOCATE(array)
    ELSE
      CALL basic_evaluate_standard(input_stack, ix, iy, err)
    ENDIF

  END SUBROUTINE basic_evaluate



  SUBROUTINE basic_evaluate_standard(input_stack, ix, iy, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i
    TYPE(stack_element) :: block

    CALL eval_reset()

    DO i = 1, input_stack%stack_point
      block = input_stack%entries(i)
      IF (block%ptype == c_pt_variable) THEN
        CALL push_on_eval(block%numerical_data)
      ELSE IF (block%ptype == c_pt_species) THEN
        CALL do_species(block%value, err)
      ELSE IF (block%ptype == c_pt_operator) THEN
        CALL do_operator(block%value, err)
      ELSE IF (block%ptype == c_pt_constant &
          .OR. block%ptype == c_pt_default_constant) THEN
        CALL do_constant(block%value, .FALSE., ix, iy, err)
      ELSE IF (block%ptype == c_pt_function) THEN
        CALL do_functions(block%value, .FALSE., ix, iy, err)
      ENDIF

      IF (err /= c_err_none) THEN
        PRINT *, 'BAD block', err, block%ptype, i, block%value
        CALL abort_code(err)
        STOP
      ENDIF
    ENDDO

  END SUBROUTINE basic_evaluate_standard



#else
  SUBROUTINE basic_evaluate(input_stack, ix, iy, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i
    TYPE(stack_element) :: block

    IF (input_stack%should_simplify) CALL simplify_stack(input_stack, err)

    CALL eval_reset()

    DO i = 1, input_stack%stack_point
      block = input_stack%entries(i)
      IF (block%ptype == c_pt_variable) THEN
        CALL push_on_eval(block%numerical_data)
      ELSE IF (block%ptype == c_pt_species) THEN
        CALL do_species(block%value, err)
      ELSE IF (block%ptype == c_pt_operator) THEN
        CALL do_operator(block%value, err)
      ELSE IF (block%ptype == c_pt_constant &
          .OR. block%ptype == c_pt_default_constant) THEN
        CALL do_constant(block%value, .FALSE., ix, iy, err)
      ELSE IF (block%ptype == c_pt_function) THEN
        CALL do_functions(block%value, .FALSE., ix, iy, err)
      ENDIF

      IF (err /= c_err_none) THEN
        PRINT *, 'BAD block', err, block%ptype, i, block%value
        CALL abort_code(err)
        STOP
      ENDIF
    ENDDO

  END SUBROUTINE basic_evaluate



#endif
  SUBROUTINE sl_append

    TYPE(stack_list), POINTER :: sl_tmp

    IF (sl_size == 0) THEN
      ALLOCATE(sl_head)
      sl_tail => sl_head
    ELSE
      ALLOCATE(sl_tmp)
      sl_tmp%prev => sl_tail
      sl_tail%next => sl_tmp
      sl_tail => sl_tmp
    ENDIF
    sl_size = sl_size + 1
    CALL initialise_stack(sl_tail%stack)

  END SUBROUTINE sl_append



  SUBROUTINE simplify_stack(input_stack, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i
    TYPE(stack_element) :: block
    TYPE(primitive_stack) :: output_stack

    ! Evaluating expressions and push the results onto eval_stack.
    ! When we reach a block which requests a time or space variable,
    ! push a special flag to the eval_stack and create a new stack containing
    ! the block (sl_tail%stack).
    ! If do_* encounters an expression whose arguments contain a space or
    ! time varying expression (eg. gauss(x+y,0,1)), it pushes garbage to the
    ! eval_stack and a flag is set. It is then fixed up by
    ! update_stack_for_block()

    ! Eg. the following expression:
    ! 2*3*x + func1(4+5,x*func2(6*7,8*9*x^2),x+y) + func3(1*2,3,func4(4+5))+1+2
    ! simplifies to:
    ! 6*x + func1(9,x*func2(42,72*x^2),x+y) + f3 + 1 + 2
    !   where f3=func3(2,3,func4(9))

    CALL eval_reset()

    sl_size = 0

    DO i = 1, input_stack%stack_point
      block = input_stack%entries(i)
      IF (block%ptype == c_pt_variable) THEN
        CALL push_on_eval(block%numerical_data)
      ELSE IF (block%ptype == c_pt_species) THEN
        CALL do_species(block%value, err)
      ELSE IF (block%ptype == c_pt_operator) THEN
        CALL do_operator(block%value, err)
        CALL update_stack_for_block(block, err)
      ELSE IF (block%ptype == c_pt_constant &
          .OR. block%ptype == c_pt_default_constant) THEN
        CALL do_constant(block%value, .TRUE., 1, 1, err)
        CALL update_stack_for_block(block, err)
      ELSE IF (block%ptype == c_pt_function) THEN
        CALL do_functions(block%value, .TRUE., 1, 1, err)
        CALL update_stack_for_block(block, err)
      ENDIF

      IF (err /= c_err_none) THEN
        PRINT *, 'BAD block', err, block%ptype, i, block%value
        CALL abort_code(err)
        STOP
      ENDIF
    ENDDO

    ! We may now just be left with a list of values on the eval_stack
    ! If so, push them onto sl_tail%stack
    i = eval_stack_stack_point
    IF (i > 0) CALL update_stack(i)

    ! Now populate output_stack with the simplified expression
    CALL initialise_stack(output_stack)
    output_stack%should_simplify = .FALSE.

    IF (sl_size > 0) THEN
      CALL append_stack(output_stack, sl_tail%stack)
      DEALLOCATE(sl_tail)
      sl_size = 0
      eval_stack_stack_point = 0
    ENDIF

    CALL deallocate_stack(input_stack)
    input_stack = output_stack

  END SUBROUTINE simplify_stack



  SUBROUTINE update_stack(nvalues)

    INTEGER, INTENT(IN) :: nvalues
    INTEGER :: sp, i, n
    REAL(num), ALLOCATABLE :: entries(:)
    INTEGER, ALLOCATABLE :: flags(:)
    TYPE(stack_list), POINTER :: sl_tmp, sl_part
    TYPE(stack_element) :: new_block

    new_block%ptype = c_pt_variable
    new_block%value = 0

    sp = eval_stack_stack_point

    ALLOCATE(entries(nvalues))
    ALLOCATE(flags(nvalues))

    n = nvalues
    DO i = 1, nvalues
      entries(n) = eval_stack_entries(sp)
      flags(n) = eval_stack_flags(sp)
      IF (flags(n) /= 0) THEN
        sl_size = sl_size - 1
        sl_part => sl_tail
        sl_tail => sl_tail%prev
      ENDIF
      n = n - 1
      sp = sp - 1
    ENDDO
    eval_stack_stack_point = sp

    CALL sl_append()
    DO i = 1, nvalues
      IF (flags(i) == 0) THEN
        new_block%numerical_data = entries(i)
        CALL push_to_stack(sl_tail%stack, new_block)
      ELSE
        CALL append_stack(sl_tail%stack, sl_part%stack)
        sl_tmp => sl_part%next
        DEALLOCATE(sl_part)
        sl_part => sl_tmp
      ENDIF
    ENDDO

    DEALLOCATE(entries)
    DEALLOCATE(flags)

  END SUBROUTINE update_stack



  SUBROUTINE update_stack_for_block(block, err)

    TYPE(stack_element), INTENT(INOUT) :: block
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: nvalues

    IF (err == c_err_other) THEN
      err = c_err_none
      ! Operator just pushed a bogus value to stack, so we'll ignore it
      eval_stack_stack_point = eval_stack_stack_point - 1
      CALL push_eval_flag()
      CALL sl_append()
      CALL push_to_stack(sl_tail%stack, block)
      IF (block%value == c_const_time) sl_tail%stack%is_time_varying = .TRUE.
      RETURN
    ENDIF

    ! Number of eval_stack entries consumed by operator
    nvalues = eval_stack_nvalues
    IF (nvalues == 0) RETURN

    eval_stack_nvalues = 0
    ! Operator just pushed a bogus value to stack, so we'll ignore it
    eval_stack_stack_point = eval_stack_stack_point - 1

    CALL update_stack(nvalues)

    CALL push_eval_flag()
    CALL push_to_stack(sl_tail%stack, block)

  END SUBROUTINE update_stack_for_block



  SUBROUTINE evaluate_at_point_to_array(input_stack, ix, iy, n_elements, &
      array, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, INTENT(IN) :: n_elements
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i

    CALL basic_evaluate(input_stack, ix, iy, err)

    IF (eval_stack_stack_point /= n_elements) err = IOR(err, c_err_bad_value)

    ! Pop off the final answers
    DO i = MIN(eval_stack_stack_point,n_elements),1,-1
      array(i) = pop_off_eval()
    ENDDO

  END SUBROUTINE evaluate_at_point_to_array



  SUBROUTINE evaluate_and_return_all(input_stack, ix, iy, n_elements, &
      array, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, INTENT(OUT) :: n_elements
    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i

    IF (ASSOCIATED(array)) DEALLOCATE(array)

    CALL basic_evaluate(input_stack, ix, iy, err)

    n_elements = eval_stack_stack_point
    ALLOCATE(array(1:n_elements))

    ! Pop off the final answers
    DO i = n_elements,1,-1
      array(i) = pop_off_eval()
    ENDDO

  END SUBROUTINE evaluate_and_return_all



  SUBROUTINE evaluate_as_list(input_stack, array, n_elements, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(OUT) :: n_elements
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i
    TYPE(stack_element) :: block

    array(1) = 0
    n_elements = 1

    DO i = 1, input_stack%stack_point
      block = input_stack%entries(i)

      IF (block%ptype == c_pt_subset) THEN
        n_elements = n_elements + 1
        array(n_elements) = block%value
      ELSE IF (block%ptype == c_pt_constant &
          .OR. block%ptype == c_pt_default_constant) THEN
        CALL do_constant(block%value, .FALSE., 1, 1, err)
        array(1) = array(1) + INT(pop_off_eval())
      ELSE IF (block%ptype /= c_pt_operator) THEN
        err = c_err_bad_value
      ENDIF

      IF (err /= c_err_none) THEN
        PRINT *, 'BAD block', err, block%ptype, i, block%value
        CALL abort_code(err)
        STOP
      ENDIF
    ENDDO

  END SUBROUTINE evaluate_as_list



  FUNCTION evaluate_at_point(input_stack, ix, iy, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, INTENT(INOUT) :: err
    REAL(num), DIMENSION(1) :: array
    REAL(num) :: evaluate_at_point

    CALL evaluate_at_point_to_array(input_stack, ix, iy, 1, array, err)
    evaluate_at_point = array(1)

  END FUNCTION evaluate_at_point



  FUNCTION evaluate(input_stack, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate

    evaluate = evaluate_at_point(input_stack, 1, 1, err)

  END FUNCTION evaluate

END MODULE evaluator
