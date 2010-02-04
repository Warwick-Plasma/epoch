MODULE evaluator_blocks

  USE shared_parser_data
  USE shared_data
  USE strings
  USE custom_parser
  USE stack
  IMPLICIT NONE

CONTAINS


  SUBROUTINE do_operator(opcode,ix,iy,err)

    INTEGER,INTENT(IN) :: opcode,ix,iy
    INTEGER,INTENT(INOUT) :: err
    REAL(num),DIMENSION(2) :: values
    REAL(num) :: val
    LOGICAL :: comp

    IF (opcode .EQ. OPCODE_PLUS) THEN
       CALL get_values(2,values)
       CALL push_on_eval(values(1)+values(2))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_MINUS) THEN
       CALL get_values(2,values)
       CALL push_on_eval(values(1)-values(2))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_UNARY_MINUS) THEN
       CALL get_values(1,values)
       CALL push_on_eval(-values(1))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_TIMES) THEN
       CALL get_values(2,values)
       CALL push_on_eval(values(1)*values(2))
       RETURN
    ENDIF

    IF (opcode .EQ. OPCODE_DIVIDE) THEN
       CALL get_values(2,values)
       CALL push_on_eval(values(1)/values(2))
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_POWER) THEN
       CALL get_values(2,values)
       CALL push_on_eval(values(1)**values(2))
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_EXPO) THEN
       CALL get_values(2,values)
       CALL push_on_eval(values(1) * 10.0_num ** values(2))
    ENDIF
    IF (opcode .EQ. OPCODE_LT) THEN
       CALL get_values(2,values)
       comp=values(1) .LT. values(2)
       val=0.0_num
       IF(comp) val=1.0_num
       CALL push_on_eval(val)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_GT) THEN
       CALL get_values(2,values)
       comp=values(1) .GT. values(2)
       val=0.0_num
       IF(comp) val=1.0_num
       CALL push_on_eval(val)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_EQ) THEN
       CALL get_values(2,values)
       comp=values(1) .EQ. values(2)
       val=0.0_num
       IF(comp) val=1.0_num
       CALL push_on_eval(val)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_OR) THEN
       CALL get_values(2,values)
       comp=(values(1) .NE. 0 .OR. values(2) .NE. 0)
       val=0.0_num
       IF (comp)val=1.0_num
       CALL push_on_eval(val)
       RETURN
    ENDIF
    IF (opcode .EQ. OPCODE_AND) THEN
       CALL get_values(2,values)
       comp=(values(1) .NE. 0 .AND. values(2) .NE. 0)
       val=0.0_num
       IF (comp)val=1.0_num
       CALL push_on_eval(val)
       RETURN
    ENDIF

  END SUBROUTINE do_operator

  SUBROUTINE do_constant(opcode,ix,iy,err)

    INTEGER,INTENT(IN) :: opcode,ix,iy
    INTEGER,INTENT(INOUT) :: err
    REAL(num) :: val

    IF (opcode .GE. CONSTANT_CUSTOM_LOWBOUND) THEN
       !Check for custom constants
       val = custom_constant(opcode,ix,iy,err)    
       IF(IAND(err,ERR_UNKNOWN_ELEMENT) == 0) CALL push_on_eval(val)
       RETURN
    ENDIF

    IF (opcode .GE. CONSTANT_DECK_LOWBOUND .AND. opcode .LT. CONSTANT_CUSTOM_LOWBOUND) THEN
       val = deck_constant_list(opcode-CONSTANT_DECK_LOWBOUND)%value
       CALL push_on_eval(val)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_PI) THEN
       CALL push_on_eval(pi)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_KB) THEN
       CALL push_on_eval(kb)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_QE) THEN
       CALL push_on_eval(q0)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_C) THEN
       CALL push_on_eval(c)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_ME) THEN
       CALL push_on_eval(m0)
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_EPS0) THEN
       CALL push_on_eval(epsilon0)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_MU0) THEN
       CALL push_on_eval(mu0)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_EV) THEN
       CALL push_on_eval(ev)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_KEV) THEN
       CALL push_on_eval(ev*1000.0_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_MEV) THEN
       CALL push_on_eval(ev*1.0e6_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_LX) THEN
       CALL push_on_eval(length_x)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_LY) THEN
       CALL push_on_eval(length_y)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_X) THEN
       CALL push_on_eval(x(ix)+dx/2.0_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_Y) THEN
       CALL push_on_eval(y(iy)+dy/2.0_num)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DX) THEN
       CALL push_on_eval(dx)
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_DY) THEN
       CALL push_on_eval(dy)
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_START_X) THEN
       CALL push_on_eval(x_start)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_START_Y) THEN
       CALL push_on_eval(y_start)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_END_X) THEN
       CALL push_on_eval(x_end)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_END_Y) THEN
       CALL push_on_eval(y_end)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IX) THEN
       CALL push_on_eval(REAL(ix,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IY) THEN
       CALL push_on_eval(REAL(iy,num))
       RETURN 
    ENDIF

    IF (opcode .EQ. CONST_TIME) THEN
       CALL push_on_eval(time)
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_AUTOEARLY) THEN
       CALL push_on_eval(REAL(IC_EARLY_INTERNAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_AUTOLATE) THEN
       CALL push_on_eval(REAL(IC_LATE_INTERNAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_EXTERNAL) THEN
       CALL push_on_eval(REAL(IC_EXTERNAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_MANUAL) THEN
       CALL push_on_eval(REAL(IC_MANUAL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_RESTART) THEN
       CALL push_on_eval(REAL(IC_RESTART,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_NEVER) THEN
       CALL push_on_eval(REAL(IO_NEVER,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_ALWAYS) THEN
       CALL push_on_eval(REAL(IO_ALWAYS,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_FULL) THEN
       CALL push_on_eval(REAL(IO_FULL,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_RESTARTABLE) THEN
       CALL push_on_eval(REAL(IO_RESTARTABLE,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_IO_SPECIES) THEN
       CALL push_on_eval(REAL(IO_SPECIES,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_X) THEN
       CALL push_on_eval(REAL(DIR_X,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_Y) THEN
       CALL push_on_eval(REAL(DIR_Y,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_PX) THEN
       CALL push_on_eval(REAL(DIR_PX,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_PY) THEN
       CALL push_on_eval(REAL(DIR_PY,num))
       RETURN
    ENDIF

    IF (opcode .EQ. CONST_DIR_PZ) THEN
       CALL push_on_eval(REAL(DIR_PZ,num))
       RETURN
    ENDIF

	IF (opcode .EQ. CONST_VAR) THEN
   	CALL push_on_eval(context_variable)
   	RETURN
	ENDIF

  END SUBROUTINE do_constant

  SUBROUTINE do_functions(opcode,ix,iy,err)

    INTEGER,INTENT(IN) :: opcode,ix,iy
    INTEGER,INTENT(INOUT) :: err
    REAL(num),DIMENSION(4) :: values
    REAL(num) :: val
    INTEGER :: count, ipoint
    LOGICAL :: done
    REAL(num),DIMENSION(:),ALLOCATABLE :: var_length_values
    REAL(num) :: point, t0


    IF (opcode .EQ. FUNC_SQRT) THEN
       CALL get_values(1,values)
       CALL push_on_eval(SQRT(values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_SINE) THEN
       CALL get_values(1,values)
       CALL push_on_eval(SIN(values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_COSINE) THEN
       CALL get_values(1,values)
       CALL push_on_eval(COS(values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TAN) THEN
       CALL get_values(1,values)
       CALL push_on_eval(TAN(values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_EXP) THEN
       CALL get_values(1,values)
       CALL push_on_eval(EXP(values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_ARCSINE) THEN
       CALL get_values(1,values)
       CALL push_on_eval(ASIN(values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_ARCCOSINE) THEN
       CALL get_values(1,values)
       CALL push_on_eval(ACOS(values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_ARCTAN) THEN
       CALL get_values(1,values)
       CALL push_on_eval(ATAN(values(1)))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_NEG) THEN
       CALL get_values(1,values)
       CALL push_on_eval(-(values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_IF) THEN
       CALL get_values(3,values)
       IF (values(1) .NE. 0) THEN
          val=values(2)
       ELSE
          val=values(3)
       ENDIF
       CALL push_on_eval(val)
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_RHO) THEN
       CALL get_values(1,values)
       CALL push_on_eval(initial_conditions(NINT(values(1)))%rho(ix,iy))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TEMPX) THEN
       CALL get_values(1,values)
       CALL push_on_eval(initial_conditions(NINT(values(1)))%temp(ix,iy,1))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TEMPY) THEN
       CALL get_values(1,values)
       CALL push_on_eval(initial_conditions(NINT(values(1)))%temp(ix,iy,2))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_TEMPZ) THEN
       CALL get_values(1,values)
       CALL push_on_eval(initial_conditions(NINT(values(1)))%temp(ix,iy,3))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_INTERPOLATE) THEN
       CALL get_values(1,values)
       count=NINT(values(1))
       ALLOCATE(var_length_values(0:count*2))
       CALL get_values(count*2+1,var_length_values)
       !This could be replaced by a bisection algorithm, change at some point
       !For now, not too bad for small count

       !var_length_values(0) = position in domain
       done=.FALSE.
       point=var_length_values(0)

       DO ipoint=1,count-1
          IF ( point .GE. var_length_values(ipoint*2-1) .AND. point .LE. var_length_values(ipoint*2+1)) THEN
             val=(point-var_length_values(ipoint*2-1))/(var_length_values(ipoint*2+1)-var_length_values(ipoint*2-1)) *&
                  (var_length_values(ipoint*2+2) - var_length_values(ipoint*2)) + var_length_values(ipoint*2)
             done=.TRUE.
             CALL push_on_eval(val)
             EXIT
          ENDIF
       ENDDO
       DEALLOCATE(var_length_values)
       IF (.NOT. done) err=ERR_BAD_VALUE
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_TANH) THEN
       CALL get_values(1,values)
       CALL push_on_eval(TANH(values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_SINH) THEN
       CALL get_values(1,values)
       CALL push_on_eval(SINH(values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_COSH) THEN
       CALL get_values(1,values)
       CALL push_on_eval(COSH(values(1)))
       RETURN
    ENDIF

    IF (opcode .EQ. FUNC_EX) THEN
       CALL get_values(2,values)
       CALL push_on_eval(ex(NINT(values(1)),NINT(values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_EY) THEN
       CALL get_values(2,values)
       CALL push_on_eval(ey(NINT(values(1)),NINT(values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_EZ) THEN
       CALL get_values(2,values)
       CALL push_on_eval(ez(NINT(values(1)),NINT(values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_BX) THEN
       CALL get_values(2,values)
       CALL push_on_eval(bx(NINT(values(1)),NINT(values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_BY) THEN
       CALL get_values(2,values)
       CALL push_on_eval(by(NINT(values(1)),NINT(values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_BZ) THEN
       CALL get_values(2,values)
       CALL push_on_eval(bz(NINT(values(1)),NINT(values(2))))
       RETURN
    ENDIF
    IF (opcode .EQ. FUNC_GAUSS) THEN
       CALL get_values(3,values)
       CALL push_on_eval(EXP(-((values(1)-values(2))/values(3))**2))
       RETURN
    ENDIF

IF (opcode .EQ. FUNC_SEMIGAUSS) THEN
	 CALL get_values(4,values)
    !values are : time, maximum amplitude, amplitude at t=0, characteristic time width
	 t0 = values(4) * SQRT(-LOG(values(3)/values(2)))
	 IF (values(1) .LE. t0) THEN
	 	CALL push_on_eval(values(2) * EXP(-((values(1)-t0)/values(4))**2))
	 ELSE
		CALL push_on_eval(values(2))
	 ENDIF
	 RETURN
ENDIF

IF (opcode .EQ. FUNC_CRIT) THEN
	CALL get_values(1,values)
	CALL push_on_eval(values(1)**2 * m0 * epsilon0 / q0**2)
	RETURN
ENDIF

    IF (opcode .EQ. FUNC_ABS) THEN
      CALL get_values(1,values)
      CALL push_on_eval(ABS(values(1)))
      RETURN
    ENDIF

    !Check for custom functions
    val = custom_function(opcode,ix,iy,err)
    IF(IAND(err,ERR_UNKNOWN_ELEMENT) == 0) THEN
       CALL push_on_eval(val)
    ENDIF

  END SUBROUTINE do_functions


END MODULE evaluator_blocks
