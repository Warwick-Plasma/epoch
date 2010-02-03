MODULE evaluator_blocks

  USE shared_parser_data
  USE shared_data
  USE strings
  USE custom_parser
  USE stack
  IMPLICIT NONE

CONTAINS


  SUBROUTINE DoOperator(opcode,ix,iy,err)

    INTEGER,INTENT(IN) :: opcode,ix,iy
    INTEGER,INTENT(INOUT) :: err
    REAL(num),DIMENSION(2) :: Values
    REAL(num) :: Result
    LOGICAL :: Comp

    IF (opcode .EQ. OPCODE_PLUS) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Values(1)+Values(2))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_MINUS) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Values(1)-Values(2))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_UNARY_MINUS) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(-Values(1))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_TIMES) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Values(1)*Values(2))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_DIVIDE) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Values(1)/Values(2))
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_POWER) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Values(1)**Values(2))
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_EXPO) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Values(1) * 10.0_num ** Values(2))
    ENDIF
    IF (opcode .EQ. OPCODE_LT) THEN
       CALL GetValues(2,Values)
       Comp=Values(1) .LT. Values(2)
       Result=0.0_num
       IF(Comp) Result=1.0_num
       CALL PushOnEval(Result)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_GT) THEN
       CALL GetValues(2,Values)
       Comp=Values(1) .GT. Values(2)
       Result=0.0_num
       IF(Comp) Result=1.0_num
       CALL PushOnEval(Result)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_EQ) THEN
       CALL GetValues(2,Values)
       Comp=Values(1) .EQ. Values(2)
       Result=0.0_num
       IF(Comp) Result=1.0_num
       CALL PushOnEval(Result)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_OR) THEN
       CALL GetValues(2,Values)
       Comp=(Values(1) .NE. 0 .OR. Values(2) .NE. 0)
       Result=0.0_num
       IF (Comp)Result=1.0_num
       CALL PushOnEval(Result)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_AND) THEN
       CALL GetValues(2,Values)
       Comp=(Values(1) .NE. 0 .AND. Values(2) .NE. 0)
       Result=0.0_num
       IF (Comp)Result=1.0_num
       CALL PushOnEval(Result)
       RETURN
    ENDIF

  END SUBROUTINE DoOperator

  SUBROUTINE DoConstant(opcode,ix,iy,err)

    INTEGER,INTENT(IN) :: opcode,ix,iy
    INTEGER,INTENT(INOUT) :: err
    REAL(num) :: result

    IF (opcode .GE. CONSTANT_CUSTOM_LOWBOUND) THEN
       !Check for custom constants
       Result = CustomConstant(opcode,ix,iy,err)    
       IF(IAND(err,ERR_UNKNOWN_ELEMENT) == 0) CALL PushOnEval(Result)
       RETURN
    ENDIF

    IF (opcode .GE. CONSTANT_DECK_LOWBOUND .AND. opcode .LT. CONSTANT_CUSTOM_LOWBOUND) THEN
       Result = Deck_Constant_List(opcode-CONSTANT_DECK_LOWBOUND)%Value
       CALL PushOnEval(Result)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_PI) THEN
       CALL PushOnEval(pi)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_KB) THEN
       CALL PushOnEval(kb)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_QE) THEN
       CALL PushOnEval(q0)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_C) THEN
       CALL PushOnEval(C)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_ME) THEN
       CALL PushOnEval(m0)
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_EPS0) THEN
       CALL PushOnEval(epsilon0)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_MU0) THEN
       CALL PushOnEval(mu0)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_EV) THEN
       CALL PushOnEval(ev)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_KEV) THEN
       CALL PushOnEval(ev*1000.0_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_MEV) THEN
       CALL PushOnEval(ev*1.0e6_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_LX) THEN
       CALL PushOnEval(length_x)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_LY) THEN
       CALL PushOnEval(length_y)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_X) THEN
       CALL PushOnEval(x(ix)+dx/2.0_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_Y) THEN
       CALL PushOnEval(y(iy)+dy/2.0_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DX) THEN
       CALL PushOnEval(dx)
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_DY) THEN
       CALL PushOnEval(dy)
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_START_X) THEN
       CALL PushOnEval(x_start)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_START_Y) THEN
       CALL PushOnEval(y_start)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_END_X) THEN
       CALL PushOnEval(x_end)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_END_Y) THEN
       CALL PushOnEval(y_end)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IX) THEN
       CALL PushOnEval(REAL(ix,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IY) THEN
       CALL PushOnEval(REAL(iy,num))
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_TIME) THEN
       CALL PushOnEval(time)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_AUTOEARLY) THEN
       CALL PushOnEval(REAL(IC_EARLY_INTERNAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_AUTOLATE) THEN
       CALL PushOnEval(REAL(IC_LATE_INTERNAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_EXTERNAL) THEN
       CALL PushOnEval(REAL(IC_EXTERNAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_MANUAL) THEN
       CALL PushOnEval(REAL(IC_MANUAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_RESTART) THEN
       CALL PushOnEval(REAL(IC_RESTART,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_NEVER) THEN
       CALL PushOnEval(REAL(IO_NEVER,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_ALWAYS) THEN
       CALL PushOnEval(REAL(IO_ALWAYS,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_FULL) THEN
       CALL PushOnEval(REAL(IO_FULL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_RESTARTABLE) THEN
       CALL PushOnEval(REAL(IO_RESTARTABLE,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_SPECIES) THEN
       CALL PushOnEval(REAL(IO_SPECIES,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_X) THEN
       CALL PushOnEval(REAL(DIR_X,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_Y) THEN
       CALL PushOnEval(REAL(DIR_Y,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_PX) THEN
       CALL PushOnEval(REAL(DIR_PX,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_PY) THEN
       CALL PushOnEval(REAL(DIR_PY,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_PZ) THEN
       CALL PushOnEval(REAL(DIR_PZ,num))
       RETURN
    ENDIF

	IF (opcode .EQ. CONST_VAR) THEN
   	CALL PushOnEval(context_variable)
   	RETURN
	ENDIF

  END SUBROUTINE DoConstant

  SUBROUTINE DoFunctions(opcode,ix,iy,err)

    INTEGER,INTENT(IN) :: opcode,ix,iy
    INTEGER,INTENT(INOUT) :: err
    REAL(num),DIMENSION(4) :: Values
    REAL(num) :: Result
    INTEGER :: Count, iPoint
    LOGICAL :: Done
    REAL(num),DIMENSION(:),ALLOCATABLE :: Var_Length_Values
    REAL(num) :: point, T0


    IF (opcode .EQ. FUNC_SQRT) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(SQRT(Values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_SINE) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(SIN(Values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_COSINE) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(COS(Values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TAN) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(TAN(Values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_EXP) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(EXP(Values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_ARCSINE) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(ASIN(Values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_ARCCOSINE) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(ACOS(Values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_ARCTAN) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(ATAN(Values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_NEG) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(-(Values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_IF) THEN
       CALL GetValues(3,Values)
       IF (Values(1) .NE. 0) THEN
          Result=Values(2)
       ELSE
          Result=Values(3)
       ENDIF
       CALL PushOnEval(Result)
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_RHO) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(InitialConditions(NINT(Values(1)))%rho(ix,iy))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TEMPX) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(InitialConditions(NINT(Values(1)))%temp(ix,iy,1))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TEMPY) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(InitialConditions(NINT(Values(1)))%temp(ix,iy,2))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TEMPZ) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(InitialConditions(NINT(Values(1)))%temp(ix,iy,3))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_INTERPOLATE) THEN
       CALL GetValues(1,Values)
       Count=NINT(Values(1))
       ALLOCATE(Var_Length_Values(0:Count*2))
       CALL GetValues(Count*2+1,Var_Length_Values)
       !This could be replaced by a bisection algorithm, change at some point
       !For now, not too bad for small count

       !Var_Length_Values(0) = position in domain
       Done=.FALSE.
       point=Var_Length_Values(0)

       DO ipoint=1,Count-1
          IF ( point .GE. Var_Length_Values(iPoint*2-1) .AND. point .LE. Var_Length_Values(iPoint*2+1)) THEN
             Result=(point-Var_Length_Values(iPoint*2-1))/(Var_Length_Values(iPoint*2+1)-Var_Length_Values(iPoint*2-1)) *&
                  (Var_Length_Values(iPoint*2+2) - Var_Length_Values(iPoint*2)) + Var_Length_Values(iPoint*2)
             Done=.TRUE.
             CALL PushOnEval(Result)
             EXIT
          ENDIF
       ENDDO
       DEALLOCATE(Var_Length_Values)
       IF (.NOT. Done) Err=ERR_BAD_VALUE
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_TANH) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(TANH(Values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_SINH) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(SINH(Values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_COSH) THEN
       CALL GetValues(1,Values)
       CALL PushOnEval(COSH(Values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_EX) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Ex(NINT(Values(1)),NINT(Values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_EY) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Ey(NINT(Values(1)),NINT(Values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_EZ) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Ez(NINT(Values(1)),NINT(Values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_BX) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(BX(NINT(Values(1)),NINT(Values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_BY) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(By(NINT(Values(1)),NINT(Values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_BZ) THEN
       CALL GetValues(2,Values)
       CALL PushOnEval(Bz(NINT(Values(1)),NINT(Values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_GAUSS) THEN
       CALL GetValues(3,Values)
       CALL PushOnEval(EXP(-((Values(1)-Values(2))/Values(3))**2))
       RETURN
    ENDIF

IF (opcode .EQ. FUNC_SEMIGAUSS) THEN
	 CALL GetValues(4,Values)
    !Values are : time, maximum amplitude, amplitude at t=0, characteristic time width
	 T0 = Values(4) * SQRT(-LOG(Values(3)/Values(2)))
	 IF (Values(1) .LE. T0) THEN
	 	CALL PushOnEval(Values(2) * EXP(-((Values(1)-T0)/Values(4))**2))
	 ELSE
		CALL PushOnEval(Values(2))
	 ENDIF
	 RETURN
ENDIF

IF (opcode .EQ. FUNC_CRIT) THEN
	CALL GetValues(1,Values)
	CALL PushOnEval(Values(1)**2 * m0 * epsilon0 / q0**2)
	RETURN
ENDIF

IF (opcode .EQ. FUNC_ABS) THEN
	CALL GetValues(1,Values)
	CALL PushOnEval(ABS(Values(1)))
	RETURN
ENDIF

    !Check for custom functions
    Result = CustomFunction(opcode,ix,iy,err)
    IF(IAND(err,ERR_UNKNOWN_ELEMENT) == 0) THEN
       CALL PushOnEval(Result)
    ENDIF

  END SUBROUTINE DoFunctions


END MODULE evaluator_blocks
