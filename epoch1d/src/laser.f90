MODULE laser

  ! This module contributed by
  ! Dr. N. J. Sircombe

  USE evaluator
  USE custom_laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_laser(direction, laser)

    INTEGER, INTENT(IN) :: direction
    TYPE(laser_block), INTENT(INOUT) :: laser

    laser%direction = direction
    NULLIFY(laser%next)

  END SUBROUTINE init_laser



  ! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    INTEGER :: direction
    TYPE(laser_block), POINTER :: laser

    direction = laser%direction

    IF (direction .EQ. c_bd_left) THEN
      CALL attach_laser_to_list(laser_left, laser, direction)
    ENDIF
    IF (direction .EQ. c_bd_right) THEN
      CALL attach_laser_to_list(laser_right, laser, direction)
    ENDIF

  END SUBROUTINE attach_laser



  FUNCTION laser_time_profile(laser)

    TYPE(laser_block), POINTER :: laser
    REAL(num) :: laser_time_profile
    INTEGER :: err

    err = 0
    IF (laser%use_time_function) THEN
      laser_time_profile = evaluate(laser%time_function, err)
      RETURN
    ENDIF

    laser_time_profile = custom_laser_time_profile(laser)

  END FUNCTION laser_time_profile



  ! Actually does the attaching of the laser to the correct list
  SUBROUTINE attach_laser_to_list(list, laser, direction)

    TYPE(laser_block), POINTER :: list
    TYPE(laser_block), POINTER :: laser
    INTEGER, INTENT(IN) :: direction
    TYPE(laser_block), POINTER :: current

    IF (ASSOCIATED(list)) THEN
      current=>list
      DO WHILE(ASSOCIATED(current%next))
        current=>current%next
      ENDDO
      current%next=>laser
    ELSE
      list=>laser
    ENDIF

  END SUBROUTINE attach_laser_to_list



  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = 1000000.0_num
    current=>laser_left
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num*pi/current%freq
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO
    current=>laser_right
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num*pi/current%freq
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO

    ! Need at least two iterations per laser period
    ! (Nyquist)
    dt_laser = dt_laser/2.0_num

  END SUBROUTINE set_laser_dt



  ! laser boundary for the left boundary
  SUBROUTINE laser_bcs_left

    REAL(num) :: t_env
    REAL(num) :: lx
    REAL(num) :: fplus

    TYPE(laser_block), POINTER :: current

    lx = dt/dx
    fplus = 0.0_num
    bx(0) =  0.0_num

    current=>laser_left
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus = fplus + t_env * current%amp &
            * SIN(current%freq*time + current%phase) &
            * SIN(current%pol) * COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    by(-1:0) = (1.0_num / (c + lx*c**2)) &
        * (-4.0_num * fplus + 2.0_num * ez(1) - (c - lx*c**2)*by(1) &
        - (dt / epsilon0) * jz(1))

    fplus = 0.0_num
    current=>laser_left
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus = fplus + t_env * current%amp &
            * SIN(current%freq*time + current%phase) * COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(-1:0) = (1.0_num / (c + lx*c**2)) &
        * (4.0_num * fplus - 2.0_num * ey(1) - (c - lx*c**2)*bz(1) &
        + (dt / epsilon0) * jy(1))

  END SUBROUTINE laser_bcs_left



  ! laser boundary for the right boundary
  SUBROUTINE laser_bcs_right

    REAL(num) :: t_env
    REAL(num) :: lx
    REAL(num) :: f_minus

    TYPE(laser_block), POINTER :: current

    lx = dt/dx

    f_minus = 0.0_num
    bx(nx+1) =  0.0_num

    current=>laser_right
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        f_minus = f_minus + t_env * current%amp &
            * SIN(current%freq*time + current%phase) &
            * SIN(current%pol) * COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    by(nx:nx+1) = (1.0_num / (c + lx*c**2)) &
        * (4.0_num * f_minus - 2.0_num * ez(nx) - (c - lx*c**2)*by(nx) &
        + (dt / epsilon0) * jz(nx))

    f_minus = 0.0_num
    current=>laser_right
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        f_minus = f_minus + t_env * current%amp * current%profile &
            * SIN(current%freq*time + current%phase) * COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(nx:nx+1) = (1.0_num / (c + lx*c**2)) &
        * (-4.0_num * f_minus + 2.0_num * ey(nx) - (c - lx*c**2)*bz(nx) &
        - (dt / epsilon0) * jy(nx))

  END SUBROUTINE laser_bcs_right



  SUBROUTINE outflow_bcs_left

    REAL(num) :: lx

    lx = dt/dx
    bx(-1:0) =  0.0_num
    ! Set the y magnetic field
    by(-1:0) = (1.0_num / (c + lx*c**2)) &
        * (2.0_num * ez(1) - (c - lx*c**2)*by(1) - (dt / epsilon0) * jz(1))
    bz(-1:0) = (1.0_num / (c + lx*c**2)) &
        * (-2.0_num * ey(1) - (c - lx*c**2)*bz(1) + (dt / epsilon0) * jy(1))

  END SUBROUTINE outflow_bcs_left



  SUBROUTINE outflow_bcs_right

    REAL(num) :: lx

    lx = dt/dx
    bx(nx) =  0.0_num
    ! Set the y magnetic field
    by(nx:nx+1) = (1.0_num / (c + lx*c**2)) &
        * (- 2.0_num * ez(nx) - (c - lx*c**2)*by(nx) + (dt / epsilon0) * jz(nx))
    bz(nx:nx+1) = (1.0_num / (c + lx*c**2)) &
        * (2.0_num * ey(nx) - (c - lx*c**2)*bz(nx) - (dt / epsilon0) * jy(nx))

  END SUBROUTINE outflow_bcs_right

END MODULE laser
