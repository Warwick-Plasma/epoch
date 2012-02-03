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
    laser%amp = -1.0_num
    laser%omega = -1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = t_end
    NULLIFY(laser%next)

    laser%profile = 1.0_num
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



  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, sum, diff, dt_eps
    REAL(num) :: fplus
    INTEGER :: i
    TYPE(laser_block), POINTER :: current

    i = c_bd_x_min

    dtc2 = dt * c**2
    lx = dtc2 / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    bx(0) = bx_x_min

    fplus = 0.0_num
    IF (add_laser(i)) THEN
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
    ENDIF

    bz(0) = sum * ( 4.0_num * fplus &
        + 2.0_num * (ey_x_min + c * bz_x_min) &
        - 2.0_num * ey(1) &
        + dt_eps * jy(1) &
        + diff * bz(1))

    IF (add_laser(i)) THEN
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
    ENDIF

    by(0) = sum * (-4.0_num * fplus &
        - 2.0_num * (ez_x_min - c * by_x_min) &
        + 2.0_num * ez(1) &
        - dt_eps * jz(1) &
        + diff * by(1))

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, sum, diff, dt_eps
    REAL(num) :: fneg
    INTEGER :: i
    TYPE(laser_block), POINTER :: current

    i = c_bd_x_max

    dtc2 = dt * c**2
    lx = dtc2 / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    bx(nx+1) = bx_x_max

    fneg = 0.0_num
    IF (add_laser(i)) THEN
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
    ENDIF

    bz(nx) = sum * (-4.0_num * fneg &
        - 2.0_num * (ey_x_max - c * bz_x_max) &
        + 2.0_num * ey(nx) &
        - dt_eps * jy(nx) &
        + diff * bz(nx-1))

    IF (add_laser(i)) THEN
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
    ENDIF

    by(nx) = sum * ( 4.0_num * fneg &
        + 2.0_num * (ez_x_max + c * by_x_max) &
        - 2.0_num * ez(nx) &
        + dt_eps * jz(nx) &
        + diff * by(nx-1))

  END SUBROUTINE outflow_bcs_x_max

END MODULE laser
