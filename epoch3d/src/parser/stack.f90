MODULE stack
  USE shared_data

  IMPLICIT NONE

  TYPE :: FloatStack
     REAL(num),DIMENSION(1000) :: Stack
     INTEGER :: StackPoint
  END TYPE FloatStack

  TYPE(FloatStack) :: EvalStack

CONTAINS

  SUBROUTINE PushOnEval(value)

    REAL(num),INTENT(IN) :: value
    EvalStack%StackPoint=EvalStack%StackPoint+1
    EvalStack%Stack(EvalStack%StackPoint)=value

  END SUBROUTINE PushOnEval

  FUNCTION PopOffEval()

    REAL(num) :: PopOffEval

    PopOffEval=EvalStack%Stack(EvalStack%StackPoint)
    EvalStack%StackPoint=EvalStack%StackPoint-1

  END FUNCTION PopOffEval

  SUBROUTINE GetValues(count,values)

    INTEGER,INTENT(IN) :: count
    REAL(num),DIMENSION(1:) :: values

    INTEGER :: i

    DO i=1,count
       values(count-i+1)=PopOffEval()
    ENDDO
  END SUBROUTINE GetValues


END MODULE stack
