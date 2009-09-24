MODULE strings

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE IntegerAsString(Integer,String)

    INTEGER,INTENT(IN) :: Integer
    CHARACTER(LEN=*),INTENT(OUT) :: String

    INTEGER :: n_nums
    CHARACTER(LEN=9) :: numfmt

    n_nums=1 + LOG10(REAL(Integer,num))
    WRITE(numfmt,'("(I",I6.6,")")') n_nums
    WRITE(string,numfmt) Integer

  END SUBROUTINE IntegerAsString

  SUBROUTINE Integer8AsString(Integer,String)

    INTEGER(KIND=8),INTENT(IN) :: Integer
    CHARACTER(LEN=*),INTENT(OUT) :: String

    INTEGER :: n_nums
    CHARACTER(LEN=12) :: numfmt

    n_nums=1 + LOG10(REAL(Integer,num))
    WRITE(numfmt,'("(I",I9.9,")")') n_nums
    WRITE(string,numfmt) Integer

  END SUBROUTINE Integer8AsString

  FUNCTION StrCmp(StrIn,StrTest)

    CHARACTER(*),INTENT(IN) ::  StrIn,StrTest
    CHARACTER(LEN=EntryLength) :: StrTrim
    LOGICAL :: StrCmp

    StrTrim=TRIM(ADJUSTL(StrIn))

    IF (LEN(StrTest) .GT. LEN(StrTrim)) THEN
       StrCmp=.FALSE.
       return
    ENDIF

    IF (LEN(StrTest) .LT. EntryLength) THEN
       IF (StrTrim(LEN(StrTest)+1:LEN(StrTest)+1) .NE. " ") THEN
          StrCmp=.FALSE.
          RETURN
       ENDIF
    ENDIF

    StrCmp=(StrTrim(1:Len(StrTest)) == StrTest)

  END FUNCTION StrCmp

  FUNCTION AsRealSimple(StrIn,ERR)
    CHARACTER(*), INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER::f
    REAL(num) :: AsRealSimple
    REAL(num) :: Value

    READ(unit=StrIn,fmt=*,iostat=f) Value
    IF (f .NE. 0) ERR=IOR(ERR,ERR_BAD_VALUE)
    AsRealSimple=Value
  END FUNCTION AsRealSimple

  FUNCTION AsIntegerSimple(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsIntegerSimple,Value,f
    READ(unit=StrIn,fmt=*,iostat=f) value
    IF (f .NE. 0) ERR=IOR(ERR,ERR_BAD_VALUE)
    AsIntegerSimple=value

  END FUNCTION AsIntegerSimple

  FUNCTION AsDirection(StrIn,ERR)
    CHARACTER(*),INTENT(IN) :: StrIn
    INTEGER,INTENT(INOUT) :: ERR
    INTEGER :: AsDirection
    AsDirection=-1

    IF (StrCmp(StrIn,"left")) AsDirection=BD_LEFT
    IF (StrCmp(StrIn,"right")) AsDirection=BD_RIGHT
    IF (StrCmp(StrIn,"up")) AsDirection=BD_UP
    IF (StrCmp(StrIn,"down")) AsDirection=BD_DOWN

    IF (AsDirection == -1) ERR=IOR(ERR,ERR_BAD_VALUE)

  END FUNCTION AsDirection

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
    IF (StrCmp(TRIM(ADJUSTL(StrIn)),"simple_outflow")) THEN
       AsBC=BC_SIMPLE_OUTFLOW
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
