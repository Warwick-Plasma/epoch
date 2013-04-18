MODULE laser

  USE custom_laser
  USE evaluator

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
      current => list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      ENDDO
      current%next => laser
    ELSE
      list => laser
    ENDIF

  END SUBROUTINE attach_laser_to_list



  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = HUGE(1.0_num)

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    ENDDO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    ENDDO

    ! Need at least two iterations per laser period
    ! (Nyquist)
    dt_laser = dt_laser / 2.0_num

  END SUBROUTINE set_laser_dt



  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, sum, diff, dt_eps, base
    REAL(num) :: source1, source2
    INTEGER :: laserpos, n
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_min

    laserpos = 1
    IF (bc_field(n) .EQ. c_bc_cpml_laser) THEN
      laserpos = cpml_x_min_laser_idx
    ENDIF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos-1) = bx_x_min

    IF (add_laser(n)) THEN
      current => laser_x_min
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
          t_env = laser_time_profile(current) * current%amp
          base = t_env * current%profile &
            * SIN(current%omega * time + current%phase)
          source1 = source1 + base * COS(current%pol_angle)
          source2 = source2 + base * SIN(current%pol_angle)
        ENDIF
        current => current%next
      ENDDO
    ENDIF

    bz(laserpos-1) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ey_x_min + c * bz_x_min) &
        - 2.0_num * ey(laserpos) &
        + dt_eps * jy(laserpos) &
        + diff * bz(laserpos))

    by(laserpos-1) = sum * (-4.0_num * source2 &
        - 2.0_num * (ez_x_min - c * by_x_min) &
        + 2.0_num * ez(laserpos) &
        - dt_eps * jz(laserpos) &
        + diff * by(laserpos))

    IF (dumpmask(c_dump_absorption) .GT. 0) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_min, lasers = laser_x_min)
      ELSE
        CALL calc_absorption(c_bd_x_min)
      ENDIF
    ENDIF

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, sum, diff, dt_eps, base
    REAL(num) :: source1, source2
    INTEGER :: laserpos, n
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_max

    laserpos = nx
    IF (bc_field(n) .EQ. c_bc_cpml_laser) THEN
      laserpos = cpml_x_max_laser_idx
    ENDIF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos+1) = bx_x_max

    IF (add_laser(n)) THEN
      current => laser_x_max
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
          t_env = laser_time_profile(current) * current%amp
          base = t_env * current%profile &
            * SIN(current%omega * time + current%phase)
          source1 = source1 + base * COS(current%pol_angle)
          source2 = source2 + base * SIN(current%pol_angle)
        ENDIF
        current => current%next
      ENDDO
    ENDIF

    bz(laserpos) = sum * (-4.0_num * source1 &
        - 2.0_num * (ey_x_max - c * bz_x_max) &
        + 2.0_num * ey(laserpos) &
        - dt_eps * jy(laserpos) &
        + diff * bz(laserpos-1))

    by(laserpos) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ez_x_max + c * by_x_max) &
        - 2.0_num * ez(laserpos) &
        + dt_eps * jz(laserpos) &
        + diff * by(laserpos-1))

    IF (dumpmask(c_dump_absorption) .GT. 0) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_max, lasers = laser_x_max)
      ELSE
        CALL calc_absorption(c_bd_x_max)
      ENDIF
    ENDIF

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE calc_absorption(bd, lasers)

    TYPE(laser_block), POINTER, OPTIONAL, INTENT(IN) :: lasers
    INTEGER, INTENT(IN) :: bd
    TYPE(laser_block), POINTER :: current
    REAL(num) :: t_env, dir
    REAL(num) :: e1, e2, b1, b2
    INTEGER :: ibc

    ! Note: ideally e1, e2, b1, b2 should be face-centred. However, this is not
    ! possible with 'open' boundaries since E-fields are not defined in the
    ! ghost cell, so we use the cell-centred quantities in the first cell.

    dir = 1.0_num

    ibc = 1
    IF (bd .EQ. c_bd_x_max) THEN
      dir = -1.0_num
      ibc = nx
    ENDIF

    e1 = ey(ibc)
    e2 = ez(ibc)
    b1 = 0.5_num * (bz(ibc-1) + bz(ibc))
    b2 = 0.5_num * (by(ibc-1) + by(ibc))

    laser_absorb_local = laser_absorb_local &
        + dt * dir * (e1 * b1 - e2 * b2) / mu0

    IF (PRESENT(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        t_env = laser_time_profile(current)
        laser_inject_local = laser_inject_local &
            + dir * dt * 0.5_num * epsilon0 * c &
            * (t_env * current%amp * current%profile)**2
        current => current%next
      ENDDO
    ENDIF

  END SUBROUTINE calc_absorption

END MODULE laser
