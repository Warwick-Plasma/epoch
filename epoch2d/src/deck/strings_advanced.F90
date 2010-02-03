MODULE strings_advanced

  USE strings

  USE shared_parser_data
  USE shunt
  USE evaluator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SplitOffInt(StrIn,StrOut,IntOut,ERR)

    CHARACTER(*),INTENT(IN) :: StrIn
    CHARACTER(*),INTENT(OUT) :: StrOut
    INTEGER,INTENT(OUT) :: IntOut
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: StrLen,Char,pos,C

    StrLen=Len(StrIn)
    pos=-1

    DO Char=1,StrLen
       C=ICHAR(StrIn(Char:Char))
       IF (C .GT. 47 .AND. C .LT. 58) THEN
          pos=Char
          EXIT
       ENDIF
    ENDDO

    IF (pos < 0) THEN 
       ERR=IOR(ERR,ERR_BAD_VALUE)
       RETURN
    ENDIF
    StrOut=StrIn(1:pos-1)
    IntOut=AsIntegerSimple(StrIn(pos:StrLen),ERR)

  END SUBROUTINE SplitOffInt

  SUBROUTINE SplitRange(StrIn,Real1,Real2,ERR)

    CHARACTER(*),INTENT(IN) :: StrIn
    REAL(num),INTENT(OUT) :: Real1,Real2
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: StrLen,Char,pos,C

    StrLen=LEN(TRIM(StrIn))
    pos=-1

    DO Char=1,StrLen
       C=ICHAR(StrIn(Char:Char))
       !Separate on a >
       IF (C .EQ. 62) THEN
          pos=Char
          EXIT
       ENDIF
    ENDDO

    IF (pos < 0) THEN 
       ERR=IOR(ERR,ERR_BAD_VALUE)
       RETURN
    ENDIF
    Real1=AsRealSimple(TRIM(StrIn(1:pos-1)),ERR)
    Real2=AsRealSimple(TRIM(StrIn(pos+1:StrLen)),ERR)

  END SUBROUTINE SplitRange


  FUNCTION AsInteger(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsInteger
    AsInteger=NINT(AsReal(StrIn,ERR))
  END FUNCTION AsInteger

  FUNCTION AsLongInteger(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER(KIND=8) :: AsLongInteger
    AsLongInteger=NINT(AsReal(StrIn,ERR))
  END FUNCTION AsLongInteger


  FUNCTION AsReal(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    REAL(num) :: AsReal
    TYPE(primitivestack) :: output

    output%stackpoint=0
    CALL Tokenize(StrIn,output,ERR)
    AsReal=Evaluate(output,ERR)
  END FUNCTION AsReal

  SUBROUTINE EvaluateStringInSpace(StrIn,DataOut,xrange,yrange,ERR)

    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER,DIMENSION(2),INTENT(IN) :: xrange, yrange
    REAL(num),DIMENSION(1:,1:),INTENT(OUT) :: DataOut 
    TYPE(primitivestack) :: output
    INTEGER :: ix,iy

    output%stackpoint=0
    CALL Tokenize(StrIn,output,ERR)

    DO iy=yrange(1),yrange(2)
       DO ix=xrange(1),xrange(2)
          DataOut(ix-xrange(1)+1,iy-yrange(1)+1) = EvaluateAtPoint(output,ix,iy,ERR)
       ENDDO
    ENDDO

  END SUBROUTINE EvaluateStringInSpace


END MODULE strings_advanced
