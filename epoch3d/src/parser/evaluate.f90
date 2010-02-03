MODULE evaluator
  USE shared_parser_data
  USE stack
  USE custom_parser
  USE evaluator_blocks
  USE shunt

  IMPLICIT NONE

CONTAINS

  FUNCTION evaluate_at_point(input_stack,ix,iy,iz,err)

    TYPE(primitive_stack), INTENT(IN) :: input_stack
    INTEGER,INTENT(IN) :: ix,iy,iz
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate_at_point
    INTEGER :: i
    TYPE(stack_element) :: block

    eval_stack%stack_point=0
!!$    IF (err .NE. c_err_none) THEN
!!$       PRINT *,"STUPID",err
!!$       STOP
!!$    ENDIF

  IF (debug_mode .AND. ix== 0 .AND. iy==0 .AND. iz==0) CALL display_tokens(input_stack)

    DO i=1,input_stack%stack_point
       block=input_stack%data(i)
       IF (block%ptype .EQ. c_pt_variable) THEN
          CALL push_on_eval(block%numerical_data)
       ENDIF
       IF (block%ptype .EQ. c_pt_operator) CALL do_operator(block%data,ix,iy,iz,err)
       IF (block%ptype .EQ. c_pt_constant) CALL do_constant(block%data,ix,iy,iz,err)
       IF (block%ptype .EQ. c_pt_function) CALL do_functions(block%data,ix,iy,iz,err)
       IF (err .NE. c_err_none) THEN
          PRINT *,"BAD block",err,block%ptype,i,block%data
          EXIT
       ENDIF
    ENDDO
    IF (eval_stack%stack_point .NE. 1) err=IAND(err,c_err_bad_value)
    !Just pop off final answer
    evaluate_at_point=pop_off_eval()
  END FUNCTION evaluate_at_point

  FUNCTION evaluate(input_stack,err)
    TYPE(primitive_stack), INTENT(IN) :: input_stack
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate

    evaluate=evaluate_at_point(input_stack,0,0,0,err)

  END FUNCTION evaluate

  END MODULE evaluator
