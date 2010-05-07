MODULE evaluator

  USE stack
  USE evaluator_blocks
  !USE shunt

  IMPLICIT NONE

CONTAINS

  SUBROUTINE evaluate_at_point_to_array(input_stack, ix, n_elements, array, err)

    TYPE(primitive_stack), INTENT(IN) :: input_stack
    INTEGER, INTENT(IN) :: ix
    INTEGER, INTENT(IN) :: n_elements
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i
    TYPE(stack_element) :: block

    eval_stack%stack_point = 0

    DO i = 1, input_stack%stack_point
      block = input_stack%data(i)
      IF (block%ptype .EQ. c_pt_variable) THEN
        CALL push_on_eval(block%numerical_data)
      ENDIF

      IF (block%ptype .EQ. c_pt_operator) &
          CALL do_operator(block%data, ix, err)
      IF (block%ptype .EQ. c_pt_constant) &
          CALL do_constant(block%data, ix, err)
      IF (block%ptype .EQ. c_pt_function) &
          CALL do_functions(block%data, ix, err)
      IF (err .NE. c_err_none) THEN
        PRINT *, "BAD block", err, block%ptype, i, block%data
        EXIT
      ENDIF
    ENDDO
    IF (eval_stack%stack_point .NE. n_elements) err = IAND(err, c_err_bad_value)
    ! Pop off the final answers
    DO i = n_elements,1,-1
      array(i) = pop_off_eval()
    ENDDO

  END SUBROUTINE evaluate_at_point_to_array

  FUNCTION evaluate_at_point(input_stack, ix, err)

    TYPE(primitive_stack), INTENT(IN) :: input_stack
    INTEGER, INTENT(IN) :: ix
    INTEGER, INTENT(INOUT) :: err
    REAL(num), DIMENSION(1) :: array
    REAL(num) :: evaluate_at_point

    CALL evaluate_at_point_to_array(input_stack, ix, 1, array, err)
    evaluate_at_point = array(1)

  END FUNCTION evaluate_at_point



  FUNCTION evaluate(input_stack, err)

    TYPE(primitive_stack), INTENT(IN) :: input_stack
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate

    evaluate = evaluate_at_point(input_stack, 0, err)

  END FUNCTION evaluate

END MODULE evaluator
