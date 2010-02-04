MODULE tokenizer_blocks

  USE shared_parser_data
  USE shared_data
  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER,PARAMETER :: MAX_NEW_FUNC=256
  INTEGER :: n_new_func = 0
  TYPE(string_type),DIMENSION(MAX_NEW_FUNC) :: new_func_name
  INTEGER, DIMENSION(MAX_NEW_FUNC) :: new_func_code

  INTEGER,PARAMETER :: MAX_NEW_CONST=256
  INTEGER :: n_new_constant = 0
  TYPE(string_type),DIMENSION(MAX_NEW_CONST) :: new_constant_name
  INTEGER, DIMENSION(MAX_NEW_CONST) :: new_constant_code
  INTEGER :: last_block_type

CONTAINS

  !Functions to register new functions and constants
  FUNCTION register_function(name)

    CHARACTER(len=*),INTENT(IN) :: name
    INTEGER :: register_function

    IF (n_new_func .EQ. MAX_NEW_FUNC) THEN
       register_function=-1
       RETURN
    ENDIF
    n_new_func=n_new_func+1
    new_func_name(n_new_func)%value=name
    new_func_code(n_new_func)=FUNC_CUSTOM_LOWBOUND+n_new_func
    register_function=FUNC_CUSTOM_LOWBOUND+n_new_func

  END FUNCTION register_function

  FUNCTION register_constant(name)

    CHARACTER(len=*),INTENT(IN) :: name
    INTEGER :: register_constant

    IF (n_new_constant .EQ. MAX_NEW_CONST) THEN
       register_constant=-1
       RETURN
    ENDIF
    n_new_constant=n_new_constant+1
    new_constant_name(n_new_constant)%value=name
    new_constant_code(n_new_constant)=CONSTANT_CUSTOM_LOWBOUND+n_new_constant
    register_constant=CONSTANT_CUSTOM_LOWBOUND+n_new_constant

  END FUNCTION register_constant


  FUNCTION as_constant(name)

    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: as_constant
    INTEGER :: i

    as_constant=PRC_NOT_THIS_TYPE

    IF (str_cmp(name,"pi")) as_constant=CONST_PI
    IF (str_cmp(name,"kb")) as_constant=CONST_KB
    IF (str_cmp(name,"me")) as_constant=CONST_ME
    IF (str_cmp(name,"qe")) as_constant=CONST_QE
    IF (str_cmp(name,"c")) as_constant=CONST_C
    IF (str_cmp(name,"epsilonnought")) as_constant=CONST_EPS0
    IF (str_cmp(name,"munought"))      as_constant=CONST_MU0
    IF (str_cmp(name,"x"))  as_constant=CONST_X
    IF (str_cmp(name,"y"))  as_constant=CONST_Y
    IF (str_cmp(name,"z"))  as_constant=CONST_Z
    IF (str_cmp(name,"lengthx") .OR. str_cmp(name,"length_x")) as_constant=CONST_LX
    IF (str_cmp(name,"lengthy") .OR. str_cmp(name,"length_y")) as_constant=CONST_LY
    IF (str_cmp(name,"lengthz") .OR. str_cmp(name,"length_z")) as_constant=CONST_LZ
    IF (str_cmp(name,"dx")) as_constant=CONST_DX
    IF (str_cmp(name,"dy")) as_constant=CONST_DY
    IF (str_cmp(name,"dy")) as_constant=CONST_DZ
    IF (str_cmp(name,"x_start")) as_constant=CONST_START_X
    IF (str_cmp(name,"y_start")) as_constant=CONST_START_Y
    IF (str_cmp(name,"y_start")) as_constant=CONST_START_Z
    IF (str_cmp(name,"x_end")) as_constant=CONST_END_X
    IF (str_cmp(name,"y_end")) as_constant=CONST_END_Y
    IF (str_cmp(name,"z_end")) as_constant=CONST_END_Z
    IF (str_cmp(name,"ix")) as_constant=CONST_IX
    IF (str_cmp(name,"iy")) as_constant=CONST_IY
    IF (str_cmp(name,"iy")) as_constant=CONST_IZ
    IF (str_cmp(name,"time")) as_constant=CONST_TIME
    IF (str_cmp(name,"internal_early")) as_constant=CONST_AUTOEARLY
    IF (str_cmp(name,"internal_late")) as_constant=CONST_AUTOLATE
    IF (str_cmp(name,"external")) as_constant=CONST_EXTERNAL
    IF (str_cmp(name,"manual"))   as_constant=CONST_MANUAL
    IF (str_cmp(name,"restart"))  as_constant=CONST_RESTART
    IF (str_cmp(name,"never")) as_constant=CONST_IO_NEVER
    IF (str_cmp(name,"always")) as_constant=CONST_IO_ALWAYS
    IF (str_cmp(name,"full")) as_constant=CONST_IO_FULL
    IF (str_cmp(name,"restartable")) as_constant=CONST_IO_RESTARTABLE
    IF (str_cmp(name,"species")) as_constant=CONST_IO_SPECIES
    IF (str_cmp(name,"dir_x")) as_constant=CONST_DIR_X
    IF (str_cmp(name,"dir_y")) as_constant=CONST_DIR_Y
    IF (str_cmp(name,"dir_z")) as_constant=CONST_DIR_Z
    IF (str_cmp(name,"dir_px")) as_constant=CONST_DIR_PX
    IF (str_cmp(name,"dir_py")) as_constant=CONST_DIR_PY
    IF (str_cmp(name,"dir_pz")) as_constant=CONST_DIR_PZ
    IF (str_cmp(name,"r_xy")) as_constant=CONST_R_XY
    IF (str_cmp(name,"r_yz")) as_constant=CONST_R_YZ
    IF (str_cmp(name,"r_xz")) as_constant=CONST_R_XZ


    !User submitted constant using "Register"
    DO i=1,n_new_constant
       IF (str_cmp(TRIM(name),TRIM(new_constant_name(i)%value))) as_constant=new_constant_code(i)
    ENDDO

    !Constants set up using the input deck
    DO i=1,n_deck_constants
       IF (str_cmp(TRIM(name),TRIM(deck_constant_list(i)%name))) THEN 
          as_constant=CONSTANT_DECK_LOWBOUND + i
          RETURN
       ENDIF
    ENDDO

  END FUNCTION as_constant

  FUNCTION as_deferred_execution_object(name)
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: as_deferred_execution_object
    INTEGER :: i

    as_deferred_execution_object=0

    DO i=1,n_deferred_execution_objects
       IF (str_cmp(TRIM(name),TRIM(deferred_objects(i)%name))) THEN 
          as_deferred_execution_object=i
          RETURN
       ENDIF
    ENDDO
  END FUNCTION as_deferred_execution_object

  FUNCTION as_function(name)

    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: as_function
    INTEGER :: i
    as_function=PRC_NOT_THIS_TYPE

    IF (str_cmp(name,"sqrt")) as_function=FUNC_SQRT
    IF (str_cmp(name,"sin"))   as_function=FUNC_SINE
    IF (str_cmp(name,"cos"))   as_function=FUNC_COSINE
    IF (str_cmp(name,"tan"))   as_function=FUNC_TAN
    IF (str_cmp(name,"exp"))   as_function=FUNC_EXP
    IF (str_cmp(name,"asin"))  as_function=FUNC_ARCSINE
    IF (str_cmp(name,"acos"))  as_function=FUNC_ARCCOSINE
    IF (str_cmp(name,"atan"))  as_function=FUNC_ARCTAN
    IF (str_cmp(name,"-")) as_function=FUNC_NEG
    IF (str_cmp(name,"if"))    as_function=FUNC_IF
    IF (str_cmp(name,"floor")) as_function=FUNC_FLOOR
    IF (str_cmp(name,"ceil"))  as_function=FUNC_CEIL
    IF (str_cmp(name,"nint"))  as_function=FUNC_NINT
    IF (str_cmp(name,"rho") .OR. str_cmp(name,"number_density"))   as_function=FUNC_RHO
    IF (str_cmp(name,"temp_x"))  as_function=FUNC_TEMPX
    IF (str_cmp(name,"temp_y"))  as_function=FUNC_TEMPY
    IF (str_cmp(name,"temp_z"))  as_function=FUNC_TEMPZ
    IF (str_cmp(name,"interpolate")) as_function=FUNC_INTERPOLATE
    IF (str_cmp(name,"tanh")) as_function=FUNC_TANH
    IF (str_cmp(name,"sinh")) as_function=FUNC_SINH
    IF (str_cmp(name,"cosh")) as_function=FUNC_COSH
    IF (str_cmp(name,"ex")) as_function=FUNC_EX
    IF (str_cmp(name,"ey")) as_function=FUNC_EY
    IF (str_cmp(name,"ez")) as_function=FUNC_EZ
    IF (str_cmp(name,"bx")) as_function=FUNC_BX
    IF (str_cmp(name,"by")) as_function=FUNC_BY
    IF (str_cmp(name,"bz")) as_function=FUNC_BZ
    IF (str_cmp(name,"gauss")) as_function=FUNC_GAUSS
    IF (str_cmp(name,"abs")) as_function=FUNC_ABS

    DO i=1,n_new_func
       IF (str_cmp(TRIM(name),TRIM(new_func_name(i)%value))) THEN
          as_function=new_func_code(i)
       ENDIF
    ENDDO

  END FUNCTION as_function

  FUNCTION as_operator(name)

    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: as_operator

    as_operator=PRC_NOT_THIS_TYPE

    IF (str_cmp(name,"+")) THEN
       as_operator=OPCODE_PLUS
    ENDIF
    IF (str_cmp(name,"-"))  THEN
       IF((last_block_type .EQ. PT_VARIABLE .OR. last_block_type .EQ. PT_CONSTANT) .OR. &
            last_block_type .EQ. PT_PARENTHESIS) THEN
          as_operator=OPCODE_MINUS
       ELSE
          as_operator=OPCODE_UNARY_MINUS
       ENDIF
    ENDIF
    IF (str_cmp(name,"*")) THEN
       as_operator=OPCODE_TIMES
    ENDIF
    IF (str_cmp(name,"/")) THEN 
       as_operator=OPCODE_DIVIDE
    ENDIF
    IF (str_cmp(name,"^")) THEN
       as_operator=OPCODE_POWER
    ENDIF
    IF (str_cmp(name,"e")) THEN
       as_operator=OPCODE_EXPO
    ENDIF
    IF (str_cmp(name,"lt")) as_operator=OPCODE_LT
    IF (str_cmp(name,"gt")) as_operator=OPCODE_GT
    IF (str_cmp(name,"eq")) as_operator=OPCODE_EQ
    IF (str_cmp(name,"and")) as_operator=OPCODE_AND
    IF (str_cmp(name,"or"))  as_operator=OPCODE_OR

  END FUNCTION as_operator


END MODULE tokenizer_blocks
