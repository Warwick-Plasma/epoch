MODULE tokenizer_blocks

  USE shared_parser_data
  USE shared_data
  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER,PARAMETER :: MAX_NEW_FUNC=256
  INTEGER :: n_new_func = 0
  TYPE(entry),DIMENSION(MAX_NEW_FUNC) :: new_func_name
  INTEGER, DIMENSION(MAX_NEW_FUNC) :: new_func_code

  INTEGER,PARAMETER :: MAX_NEW_CONST=256
  INTEGER :: n_new_constant = 0
  TYPE(entry),DIMENSION(MAX_NEW_CONST) :: new_constant_name
  INTEGER, DIMENSION(MAX_NEW_CONST) :: new_constant_code
  INTEGER :: Last_Block_Type

CONTAINS

  !Functions to register new functions and constants
  FUNCTION RegisterFunction(name)

    CHARACTER(len=*),INTENT(IN) :: name
    INTEGER :: RegisterFunction

    IF (n_new_func .EQ. MAX_NEW_FUNC) THEN
       RegisterFunction=-1
       RETURN
    ENDIF
    n_new_func=n_new_func+1
    new_func_name(n_new_func)%Value=name
    new_func_code(n_new_func)=FUNC_CUSTOM_LOWBOUND+n_new_func
    RegisterFunction=FUNC_CUSTOM_LOWBOUND+n_new_func

  END FUNCTION RegisterFunction

  FUNCTION RegisterConstant(name)

    CHARACTER(len=*),INTENT(IN) :: name
    INTEGER :: RegisterConstant

    IF (n_new_constant .EQ. MAX_NEW_CONST) THEN
       RegisterConstant=-1
       RETURN
    ENDIF
    n_new_constant=n_new_constant+1
    new_constant_name(n_new_constant)%Value=name
    new_constant_code(n_new_constant)=CONSTANT_CUSTOM_LOWBOUND+n_new_constant
    RegisterConstant=CONSTANT_CUSTOM_LOWBOUND+n_new_constant

  END FUNCTION RegisterConstant


  FUNCTION AsConstant(name)

    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: AsConstant
    INTEGER :: i

    AsConstant=PRC_NOT_THIS_TYPE

    IF (StrCmp(name,"pi")) AsConstant=CONST_PI
    IF (StrCmp(name,"kb")) AsConstant=CONST_KB
    IF (StrCmp(name,"me")) AsConstant=CONST_ME
    IF (StrCmp(name,"qe")) AsConstant=CONST_QE
    IF (StrCmp(name,"c")) AsConstant=CONST_C
    IF (StrCmp(name,"epsilonnought")) AsConstant=CONST_EPS0
    IF (StrCmp(name,"munought"))      AsConstant=CONST_MU0
    IF (StrCmp(name,"ev")) AsConstant=CONST_EV
    IF (StrCmp(name,"kev")) AsConstant=CONST_KEV
    IF (StrCmp(name,"mev")) AsConstant=CONST_MEV
    IF (StrCmp(name,"x"))  AsConstant=CONST_X
    IF (StrCmp(name,"lengthx")) AsConstant=CONST_LX
    IF (StrCmp(name,"dx")) AsConstant=CONST_DX
    IF (StrCmp(name,"x_start")) AsConstant=CONST_START_X
    IF (StrCmp(name,"y_start")) AsConstant=CONST_START_Y
    IF (StrCmp(name,"x_end")) AsConstant=CONST_END_X
    IF (StrCmp(name,"ix")) AsConstant=CONST_IX
    IF (StrCmp(name,"time")) AsConstant=CONST_TIME
    IF (StrCmp(name,"internal_early")) AsConstant=CONST_AUTOEARLY
    IF (StrCmp(name,"internal_late")) AsConstant=CONST_AUTOLATE
    IF (StrCmp(name,"external")) AsConstant=CONST_EXTERNAL
    IF (StrCmp(name,"manual"))   AsConstant=CONST_MANUAL
    IF (StrCmp(name,"restart"))  AsConstant=CONST_RESTART
    IF (StrCmp(name,"never")) AsConstant=CONST_IO_NEVER
    IF (StrCmp(name,"always")) AsConstant=CONST_IO_ALWAYS
    IF (StrCmp(name,"full")) AsConstant=CONST_IO_FULL
    IF (StrCmp(name,"restartable")) AsConstant=CONST_IO_RESTARTABLE
    IF (StrCmp(name,"species")) AsConstant=CONST_IO_SPECIES
    IF (StrCmp(name,"dir_x")) AsConstant=CONST_DIR_X
    IF (StrCmp(name,"dir_px")) AsConstant=CONST_DIR_PX
    IF (StrCmp(name,"dir_py")) AsConstant=CONST_DIR_PY
    IF (StrCmp(name,"dir_pz")) AsConstant=CONST_DIR_PZ


    !User submitted constant using "Register"
    DO i=1,n_new_constant
       IF (StrCmp(TRIM(name),TRIM(new_constant_name(i)%Value))) AsConstant=new_constant_code(i)
    ENDDO

    !Constants set up using the input deck
    DO i=1,n_Deck_Constants
       IF (StrCmp(TRIM(name),TRIM(Deck_Constant_List(i)%Name))) THEN 
          AsConstant=CONSTANT_DECK_LOWBOUND + i
          RETURN
       ENDIF
    ENDDO

  END FUNCTION AsConstant

  FUNCTION AsDeferredExecutionObject(name)
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: AsDeferredExecutionObject
    INTEGER :: i

    AsDeferredExecutionObject=0

    DO i=1,n_Deferred_Execution_Objects
       IF (StrCmp(TRIM(name),TRIM(Deferred_Objects(i)%Name))) THEN 
          AsDeferredExecutionObject=i
          RETURN
       ENDIF
    ENDDO
  END FUNCTION AsDeferredExecutionObject

  FUNCTION AsFunction(name)

    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: AsFunction
    INTEGER :: i
    AsFunction=PRC_NOT_THIS_TYPE

    IF (StrCmp(name,"abs")) AsFunction=FUNC_ABS
    IF (StrCmp(name,"sqrt")) AsFunction=FUNC_SQRT
    IF (StrCmp(name,"sin"))   AsFunction=FUNC_SINE
    IF (StrCmp(name,"cos"))   AsFunction=FUNC_COSINE
    IF (StrCmp(name,"tan"))   AsFunction=FUNC_TAN
    IF (StrCmp(name,"exp"))   AsFunction=FUNC_EXP
    IF (StrCmp(name,"asin"))  AsFunction=FUNC_ARCSINE
    IF (StrCmp(name,"acos"))  AsFunction=FUNC_ARCCOSINE
    IF (StrCmp(name,"atan"))  AsFunction=FUNC_ARCTAN
    IF (StrCmp(name,"-")) AsFunction=FUNC_NEG
    IF (StrCmp(name,"if"))    AsFunction=FUNC_IF
    IF (StrCmp(name,"floor")) AsFunction=FUNC_FLOOR
    IF (StrCmp(name,"ceil"))  AsFunction=FUNC_CEIL
    IF (StrCmp(name,"nint"))  AsFunction=FUNC_NINT
    IF (StrCmp(name,"rho"))   AsFunction=FUNC_RHO
    IF (StrCmp(name,"temp_x"))  AsFunction=FUNC_TEMPX
    IF (StrCmp(name,"temp_y"))  AsFunction=FUNC_TEMPY
    IF (StrCmp(name,"temp_z"))  AsFunction=FUNC_TEMPZ
    IF (StrCmp(name,"interpolate")) AsFunction=FUNC_INTERPOLATE
    IF (StrCmp(name,"tanh")) AsFunction=FUNC_TANH
    IF (StrCmp(name,"sinh")) AsFunction=FUNC_SINH
    IF (StrCmp(name,"cosh")) AsFunction=FUNC_COSH
    IF (StrCmp(name,"ex")) AsFunction=FUNC_EX
    IF (StrCmp(name,"ey")) AsFunction=FUNC_EY
    IF (StrCmp(name,"ez")) AsFunction=FUNC_EZ
    IF (StrCmp(name,"bx")) AsFunction=FUNC_BX
    IF (StrCmp(name,"by")) AsFunction=FUNC_BY
    IF (StrCmp(name,"bz")) AsFunction=FUNC_BZ
    IF (StrCmp(name,"gauss")) AsFunction=FUNC_GAUSS
	 IF (StrCmp(name,"semigauss")) AsFunction=FUNC_SEMIGAUSS
	 IF (StrCmp(name,"critical")) AsFunction=FUNC_CRIT

    DO i=1,n_new_func
       IF (StrCmp(TRIM(name),TRIM(new_func_name(i)%Value))) THEN
          AsFunction=new_func_code(i)
       ENDIF
    ENDDO

  END FUNCTION AsFunction

  FUNCTION AsOperator(name)

    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: AsOperator

    AsOperator=PRC_NOT_THIS_TYPE

    IF (StrCmp(name,"+")) THEN
       AsOperator=OPCODE_PLUS
    ENDIF
    IF (StrCmp(name,"-"))  THEN
       IF((last_block_type .EQ. PT_VARIABLE .OR. last_block_type .EQ. PT_CONSTANT)) THEN
          AsOperator=OPCODE_MINUS
       ELSE
          AsOperator=OPCODE_UNARY_MINUS
       ENDIF
    ENDIF
    IF (StrCmp(name,"*")) THEN
       AsOperator=OPCODE_TIMES
    ENDIF
    IF (StrCmp(name,"/")) THEN 
       AsOperator=OPCODE_DIVIDE
    ENDIF
    IF (StrCmp(name,"^")) THEN
       AsOperator=OPCODE_POWER
    ENDIF
    IF (StrCmp(name,"e")) THEN
       AsOperator=OPCODE_EXPO
    ENDIF
    IF (StrCmp(name,"lt")) AsOperator=OPCODE_LT
    IF (StrCmp(name,"gt")) AsOperator=OPCODE_GT
    IF (StrCmp(name,"eq")) AsOperator=OPCODE_EQ
    IF (StrCmp(name,"and")) AsOperator=OPCODE_AND
    IF (StrCmp(name,"or"))  AsOperator=OPCODE_OR

  END FUNCTION AsOperator


END MODULE tokenizer_blocks
