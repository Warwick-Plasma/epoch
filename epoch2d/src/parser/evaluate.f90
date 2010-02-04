MODULE evaluator
  USE shared_parser_data
  USE stack
  USE custom_parser
  USE evaluator_blocks

  IMPLICIT NONE

CONTAINS

  FUNCTION evaluate_at_point(inputstack,ix,iy,err)

    TYPE(primitive_stack), INTENT(IN) :: inputstack
    INTEGER,INTENT(IN) :: ix,iy
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate_at_point
    INTEGER :: i
    TYPE(stack_element) :: block

    eval_stack%stack_point=0

    DO i=1,inputstack%stack_point
       block=inputstack%data(i)
       IF (block%ptype .EQ. PT_VARIABLE) THEN
          CALL push_on_eval(block%numerical_data)
       ENDIF
       IF (block%ptype .EQ. PT_OPERATOR) CALL do_operator(block%data,ix,iy,err)
       IF (block%ptype .EQ. PT_CONSTANT) CALL do_constant(block%data,ix,iy,err)
       IF (block%ptype .EQ. PT_FUNCTION) CALL do_functions(block%data,ix,iy,err)
       IF (err .NE. ERR_NONE) THEN
          PRINT *,"BAD block",err,block%ptype,i,block%data
          EXIT
       ENDIF
    ENDDO
    IF (eval_stack%stack_point .NE. 1) err=IAND(err,ERR_BAD_VALUE)
    !Just pop off final answer
    evaluate_at_point=pop_off_eval()
  END FUNCTION evaluate_at_point

  FUNCTION evaluate(inputstack,err)
    TYPE(primitive_stack), INTENT(IN) :: inputstack
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate

    evaluate=evaluate_at_point(inputstack,0,0,err)

  END FUNCTION evaluate

  END MODULE evaluator
