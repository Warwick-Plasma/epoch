MODULE evaluator
  USE shared_parser_data
  USE stack
  USE custom_parser
  USE evaluator_blocks
  USE shunt

  IMPLICIT NONE

CONTAINS

  FUNCTION evaluate_at_point(inputstack,ix,iy,iz,err)

    TYPE(primitive_stack), INTENT(IN) :: inputstack
    INTEGER,INTENT(IN) :: ix,iy,iz
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate_at_point
    INTEGER :: i
    TYPE(stack_element) :: block

    eval_stack%stack_point=0
!!$    IF (err .NE. ERR_NONE) THEN
!!$       PRINT *,"STUPID",err
!!$       STOP
!!$    ENDIF

	IF (debug_mode .AND. ix== 0 .AND. iy==0 .AND. iz==0) CALL display_tokens(inputstack)

    DO i=1,inputstack%stack_point
       block=inputstack%data(i)
       IF (block%ptype .EQ. PT_VARIABLE) THEN
          CALL push_on_eval(block%numerical_data)
       ENDIF
       IF (block%ptype .EQ. PT_OPERATOR) CALL do_operator(block%data,ix,iy,iz,err)
       IF (block%ptype .EQ. PT_CONSTANT) CALL do_constant(block%data,ix,iy,iz,err)
       IF (block%ptype .EQ. PT_FUNCTION) CALL do_functions(block%data,ix,iy,iz,err)
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

    evaluate=evaluate_at_point(inputstack,0,0,0,err)

  END FUNCTION evaluate

  END MODULE evaluator
