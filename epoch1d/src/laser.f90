MODULE laser

  ! This module contributed by
  ! Dr. N. J. Sircombe

  USE evaluator
  USE custom_laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_laser(boundary, laser)

    INTEGER, INTENT(IN) :: boundary
    TYPE(laser_block), INTENT(INOUT) :: laser

    laser%boundary = boundary
    laser%id = -1
    laser%use_time_function = .FALSE.
    laser%amp = 0.0_num
    laser%omega = 1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = 0.0_num
    NULLIFY(laser%next)

    laser%profile = 0.0_num
    laser%phase = 0.0_num

  END SUBROUTINE init_laser



  ! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    INTEGER :: boundary
    TYPE(laser_block), POINTER :: laser

    boundary = laser%boundary

    IF (boundary .EQ. c_bd_x_min) THEN
      CALL attach_laser_to_list(laser_x_min, laser)
    ELSE IF (boundary .EQ. c_bd_x_max) THEN
      CALL attach_laser_to_list(laser_x_max, laser)
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
  SUBROUTINE attach_laser_to_list(list, laser)

    TYPE(laser_block), POINTER :: list
    TYPE(laser_block), POINTER :: laser
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

    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO

    ! Need at least two iterations per laser period
    ! (Nyquist)
    dt_laser = dt_laser / 2.0_num

  END SUBROUTINE set_laser_dt



  SUBROUTINE laser_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: lx, sum, diff, dt_eps
    REAL(num) :: fplus
    TYPE(laser_block), POINTER :: current

    lx = c**2 * dt / dx
    sum = 1.0_num / (c + lx)
    diff = c - lx
    dt_eps = dt / epsilon0

    bx(0) =  0.0_num

    fplus = 0.0_num
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus = fplus + t_env * current%amp * current%profile &
            * SIN(current%omega * time + current%phase) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    by(0) = sum * (-4.0_num * fplus &
        + 2.0_num * ez(1) &
        - dt_eps * jz(1) &
        - diff * by(1))

    fplus = 0.0_num
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus = fplus + t_env * current%amp * current%profile &
            * SIN(current%omega * time + current%phase) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(0) = sum * ( 4.0_num * fplus &
        - 2.0_num * ey(1) &
        + dt_eps * jy(1) &
        - diff * bz(1))

  END SUBROUTINE laser_bcs_x_min



  SUBROUTINE laser_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: lx, sum, diff, dt_eps
    REAL(num) :: fneg
    TYPE(laser_block), POINTER :: current

    lx = c**2 * dt / dx
    sum = 1.0_num / (c + lx)
    diff = c - lx
    dt_eps = dt / epsilon0

    bx(nx+1) =  0.0_num

    fneg = 0.0_num
    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg = fneg + t_env * current%amp * current%profile &
            * SIN(current%omega * time + current%phase) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    by(nx) = sum * ( 4.0_num * fneg &
        - 2.0_num * ez(nx) &
        + dt_eps * jz(nx) &
        - diff * by(nx-1))

    fneg = 0.0_num
    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg = fneg + t_env * current%amp * current%profile &
            * SIN(current%omega * time + current%phase) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(nx) = sum * (-4.0_num * fneg &
        + 2.0_num * ey(nx) &
        - dt_eps * jy(nx) &
        - diff * bz(nx-1))

  END SUBROUTINE laser_bcs_x_max



  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: lx, sum, diff, dt_eps

    lx = c**2 * dt / dx
    sum = 1.0_num / (c + lx)
    diff = c - lx
    dt_eps = dt / epsilon0

    bx(0) = 0.0_num
    by(0) = sum * ( 2.0_num * ez(1) &
        - dt_eps * jz(1) - diff * by(1))
    bz(0) = sum * (-2.0_num * ey(1) &
        + dt_eps * jy(1) - diff * bz(1))

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: lx, sum, diff, dt_eps

    lx = c**2 * dt / dx
    sum = 1.0_num / (c + lx)
    diff = c - lx
    dt_eps = dt / epsilon0

    bx(nx+1) = 0.0_num
    by(nx) = sum * (-2.0_num * ez(nx) &
        + dt_eps * jz(nx) - diff * by(nx-1))
    bz(nx) = sum * ( 2.0_num * ey(nx) &
        - dt_eps * jy(nx) - diff * bz(nx-1))

  END SUBROUTINE outflow_bcs_x_max

END MODULE laser
