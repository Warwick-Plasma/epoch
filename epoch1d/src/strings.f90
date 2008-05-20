MODULE strings

  USE shared_data

  IMPLICIT NONE

CONTAINS
  FUNCTION StrCmp(StrIn,StrTest)

    CHARACTER(*),INTENT(IN) ::  StrIn,StrTest
    CHARACTER(30) :: StrTrim
    LOGICAL :: StrCmp

    StrTrim=TRIM(ADJUSTL(StrIn))

    IF (LEN(StrTest) .GT. LEN(StrIn)) THEN
       StrCmp=.FALSE.
       return
    ENDIF

    IF (StrTrim(LEN(StrTest)+1:LEN(StrTest)+1) .NE. " ") THEN
       StrCmp=.FALSE.
       RETURN
    ENDIF

    StrCmp=StrTrim(1:Len(StrTest)) == StrTest

  END FUNCTION StrCmp

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
    IntOut=AsInteger(StrIn(pos:StrLen),ERR)
    
  END SUBROUTINE SplitOffInt

  FUNCTION AsInteger(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsInteger,Value,f
    READ(unit=StrIn,fmt=*,iostat=f) value
    IF (f .NE. 0) ERR=IOR(ERR,ERR_BAD_VALUE)
    AsInteger=value

  END FUNCTION AsInteger

  FUNCTION AsReal(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER::f
    REAL(num) :: AsReal
    REAL(num) :: Value

    READ(unit=StrIn,fmt=*,iostat=f) Value
    IF (f .NE. 0) ERR=IOR(ERR,ERR_BAD_VALUE)
    AsReal=Value
  END FUNCTION AsReal

  FUNCTION AsLogical(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    LOGICAL :: AsLogical

    AsLogical=.FALSE.
    IF(StrCmp(TRIM(ADJUSTL(StrIn)),"T")) THEN
       AsLogical=.TRUE.
       RETURN
    ENDIF
    IF(StrCmp(TRIM(ADJUSTL(StrIn)),"F")) THEN
       AsLogical=.FALSE.
       RETURN
    ENDIF

    ERR=IOR(ERR,ERR_BAD_VALUE)

  END FUNCTION AsLogical

  FUNCTION AsBC(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsBC

    AsBC=-1

    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"periodic")) THEN
       AsBC=BC_PERIODIC
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"simple_laser")) THEN
       AsBC=BC_SIMPLE_LASER
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"other")) THEN
       AsBC=BC_OTHER
       RETURN
    ENDIF

    ERR=IOR(ERR,ERR_BAD_VALUE)

  END FUNCTION AsBC

  FUNCTION AsDumpParam(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsDumpParam

    AsDumpParam=-1

    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"never")) THEN
       AsDumpParam=IO_NEVER
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"always")) THEN
       AsDumpParam=IO_ALWAYS
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"full")) THEN
       AsDumpParam=IO_FULL
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"restart")) THEN
       AsDumpParam=IO_RESTARTABLE
       RETURN
    ENDIF

    ERR=IOR(ERR,ERR_BAD_VALUE)

  END FUNCTION AsDumpParam

  FUNCTION AsDomain(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsDomain

    AsDomain=-1

    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"decomposed")) THEN
       AsDomain=DO_DECOMPOSED
       RETURN
    ENDIF
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"full")) THEN
       AsDomain=DO_FULL
       RETURN
    ENDIF

    ERR=IOR(ERR,ERR_BAD_VALUE)
  END FUNCTION AsDomain

END MODULE strings
