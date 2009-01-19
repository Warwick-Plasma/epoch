MODULE evaluator
  USE shared_parser_data
  USE stack
  USE custom_parser
  USE evaluator_blocks

  IMPLICIT NONE

CONTAINS

  FUNCTION EvaluateAtPoint(inputstack,ix,iy,ERR)

    TYPE(primitivestack), INTENT(IN) :: inputstack
    INTEGER,INTENT(IN) :: ix,iy
    INTEGER, INTENT(INOUT) :: ERR
    REAL(num) :: EvaluateAtPoint
    INTEGER :: i
    TYPE(stackelement) :: block

    EvalStack%StackPoint=0
!!$    IF (ERR .NE. ERR_NONE) THEN
!!$       PRINT *,"STUPID",ERR
!!$       STOP
!!$    ENDIF

    DO i=1,inputstack%stackpoint
       block=inputstack%Data(i)
       IF (block%Type .EQ. PT_VARIABLE) THEN
          CALL PushOnEval(block%NumericalData)
       ENDIF
       IF (block%Type .EQ. PT_OPERATOR) CALL DoOperator(block%Data,ix,iy,err)
       IF (block%Type .EQ. PT_CONSTANT) CALL DoConstant(block%Data,ix,iy,err)
       IF (block%Type .EQ. PT_FUNCTION) CALL DoFunctions(block%Data,ix,iy,err)
       IF (err .NE. ERR_NONE) THEN
          PRINT *,"BAD block",err,Block%Type,i,Block%Data
          EXIT
       ENDIF
    ENDDO
    IF (EvalStack%StackPoint .NE. 1) ERR=IAND(ERR,ERR_BAD_VALUE)
    !Just pop off final answer
    EvaluateAtPoint=PopOffEval()
  END FUNCTION EvaluateAtPoint

  FUNCTION Evaluate(inputstack,ERR)
    TYPE(primitivestack), INTENT(IN) :: inputstack
    INTEGER, INTENT(INOUT) :: ERR
    REAL(num) :: Evaluate

    Evaluate=EvaluateAtPoint(inputstack,0,0,ERR)

  END FUNCTION Evaluate

  END MODULE evaluator
