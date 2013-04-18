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
    laser%use_phase_function = .FALSE.
    laser%amp = -1.0_num
    laser%omega = -1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = t_end
    NULLIFY(laser%profile)
    NULLIFY(laser%phase)
    NULLIFY(laser%next)

    CALL allocate_with_boundary(laser%profile, boundary)
    CALL allocate_with_boundary(laser%phase, boundary)
    laser%profile = 1.0_num
    laser%phase = 0.0_num

  END SUBROUTINE init_laser



  SUBROUTINE delete_laser(laser)

    TYPE(laser_block), INTENT(INOUT) :: laser

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)

  END SUBROUTINE delete_laser



  ! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    INTEGER :: boundary
    TYPE(laser_block), POINTER :: laser

    boundary = laser%boundary

    IF (boundary .EQ. c_bd_x_min) THEN
      CALL attach_laser_to_list(laser_x_min, laser)
    ELSE IF (boundary .EQ. c_bd_x_max) THEN
      CALL attach_laser_to_list(laser_x_max, laser)
    ELSE IF (boundary .EQ. c_bd_y_min) THEN
      CALL attach_laser_to_list(laser_y_min, laser)
    ELSE IF (boundary .EQ. c_bd_y_max) THEN
      CALL attach_laser_to_list(laser_y_max, laser)
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



  SUBROUTINE laser_update_phase(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, err

    err = 0
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO i = 1,ny
          laser%phase(i) = &
              evaluate_at_point(laser%phase_function, 0, i, err)
        ENDDO
      CASE(c_bd_y_min, c_bd_y_max)
        DO i = 1,nx
          laser%phase(i) = &
              evaluate_at_point(laser%phase_function, i, 0, err)
        ENDDO
    END SELECT

  END SUBROUTINE laser_update_phase



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



  SUBROUTINE allocate_with_boundary(array, boundary)

    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: boundary

    IF (boundary .EQ. c_bd_x_min .OR. boundary .EQ. c_bd_x_max) THEN
      ALLOCATE(array(-2:ny+3))
    ELSE IF (boundary .EQ. c_bd_y_min .OR. boundary .EQ. c_bd_y_max) THEN
      ALLOCATE(array(-2:nx+3))
    ENDIF

  END SUBROUTINE allocate_with_boundary



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

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    ENDDO

    current => laser_y_max
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
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_min

    laserpos = 1
    IF (bc_field(n) .EQ. c_bc_cpml_laser) THEN
      laserpos = cpml_x_min_laser_idx
    ENDIF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(ny))
    ALLOCATE(source2(ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos-1, 1:ny) = bx_x_min(1:ny)

    IF (add_laser(n)) THEN
      current => laser_x_min
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,ny
            base = t_env * current%profile(i) &
              * SIN(current%omega * time + current%phase(i))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          ENDDO
        ENDIF
        current => current%next
      ENDDO
    ENDIF

    bz(laserpos-1, 1:ny) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ey_x_min(1:ny) + c * bz_x_min(1:ny)) &
        - 2.0_num * ey(laserpos, 1:ny) &
        + dt_eps * jy(laserpos, 1:ny) &
        + diff * bz(laserpos, 1:ny))

    by(laserpos-1, 1:ny) = sum * (-4.0_num * source2 &
        - 2.0_num * (ez_x_min(1:ny) - c * by_x_min(1:ny)) &
        + 2.0_num * ez(laserpos, 1:ny) &
        - ly * (bx(laserpos, 1:ny) - bx(laserpos, 0:ny-1)) &
        - dt_eps * jz(laserpos, 1:ny) &
        + diff * by(laserpos, 1:ny))

    DEALLOCATE(source1, source2)

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
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_max

    laserpos = nx
    IF (bc_field(n) .EQ. c_bc_cpml_laser) THEN
      laserpos = cpml_x_max_laser_idx
    ENDIF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(ny))
    ALLOCATE(source2(ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos+1, 1:ny) = bx_x_max(1:ny)

    IF (add_laser(n)) THEN
      current => laser_x_max
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,ny
            base = t_env * current%profile(i) &
              * SIN(current%omega * time + current%phase(i))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          ENDDO
        ENDIF
        current => current%next
      ENDDO
    ENDIF

    bz(laserpos, 1:ny) = sum * (-4.0_num * source1 &
        - 2.0_num * (ey_x_max(1:ny) - c * bz_x_max(1:ny)) &
        + 2.0_num * ey(laserpos, 1:ny) &
        - dt_eps * jy(laserpos, 1:ny) &
        + diff * bz(laserpos-1, 1:ny))

    by(laserpos, 1:ny) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ez_x_max(1:ny) + c * by_x_max(1:ny)) &
        - 2.0_num * ez(laserpos, 1:ny) &
        + ly * (bx(laserpos, 1:ny) - bx(laserpos, 0:ny-1)) &
        + dt_eps * jz(laserpos, 1:ny) &
        + diff * by(laserpos-1, 1:ny))

    DEALLOCATE(source1, source2)

    IF (dumpmask(c_dump_absorption) .GT. 0) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_max, lasers = laser_x_max)
      ELSE
        CALL calc_absorption(c_bd_x_max)
      ENDIF
    ENDIF

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE outflow_bcs_y_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_min

    laserpos = 1
    IF (bc_field(n) .EQ. c_bc_cpml_laser) THEN
      laserpos = cpml_y_min_laser_idx
    ENDIF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(nx))
    ALLOCATE(source2(nx))
    source1 = 0.0_num
    source2 = 0.0_num

    by(1:nx, laserpos-1) = by_y_min(1:nx)

    IF (add_laser(n)) THEN
      current => laser_y_min
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,nx
            base = t_env * current%profile(i) &
              * SIN(current%omega * time + current%phase(i))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          ENDDO
        ENDIF
        current => current%next
      ENDDO
    ENDIF

    bx(1:nx, laserpos-1) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ez_y_min(1:nx) + c * bx_y_min(1:nx)) &
        - 2.0_num * ez(1:nx, laserpos) &
        - lx * (by(1:nx, laserpos) - by(0:nx-1, laserpos)) &
        + dt_eps * jz(1:nx, laserpos) &
        + diff * bx(1:nx, laserpos))

    bz(1:nx, laserpos-1) = sum * (-4.0_num * source2 &
        - 2.0_num * (ex_y_min(1:nx) - c * bz_y_min(1:nx)) &
        + 2.0_num * ex(1:nx, laserpos) &
        - dt_eps * jx(1:nx, laserpos) &
        + diff * bz(1:nx, laserpos))

    DEALLOCATE(source1, source2)

    IF (dumpmask(c_dump_absorption) .GT. 0) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_min, lasers = laser_y_min)
      ELSE
        CALL calc_absorption(c_bd_y_min)
      ENDIF
    ENDIF

  END SUBROUTINE outflow_bcs_y_min



  SUBROUTINE outflow_bcs_y_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_max

    laserpos = ny
    IF (bc_field(n) .EQ. c_bc_cpml_laser) THEN
      laserpos = cpml_y_max_laser_idx
    ENDIF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(nx))
    ALLOCATE(source2(nx))
    source1 = 0.0_num
    source2 = 0.0_num

    by(1:nx, laserpos+1) = by_y_max(1:nx)

    IF (add_laser(n)) THEN
      current => laser_y_max
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,nx
            base = t_env * current%profile(i) &
              * SIN(current%omega * time + current%phase(i))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          ENDDO
        ENDIF
        current => current%next
      ENDDO
    ENDIF

    bx(1:nx, laserpos) = sum * (-4.0_num * source1 &
        - 2.0_num * (ez_y_max(1:nx) - c * bx_y_max(1:nx)) &
        + 2.0_num * ez(1:nx, laserpos) &
        + lx * (by(1:nx, laserpos) - by(0:nx-1, laserpos)) &
        - dt_eps * jz(1:nx, laserpos) &
        + diff * bx(1:nx, laserpos-1))

    bz(1:nx, laserpos) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ex_y_max(1:nx) + c * bz_y_max(1:nx)) &
        - 2.0_num * ex(1:nx, laserpos) &
        + dt_eps * jx(1:nx, laserpos) &
        + diff * bz(1:nx, laserpos-1))

    DEALLOCATE(source1, source2)

    IF (dumpmask(c_dump_absorption) .GT. 0) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_max, lasers = laser_y_max)
      ELSE
        CALL calc_absorption(c_bd_y_max)
      ENDIF
    ENDIF

  END SUBROUTINE outflow_bcs_y_max



  SUBROUTINE calc_absorption(bd, lasers)

    TYPE(laser_block), POINTER, OPTIONAL, INTENT(IN) :: lasers
    INTEGER, INTENT(IN) :: bd
    TYPE(laser_block), POINTER :: current
    REAL(num) :: t_env, dir, dd
    REAL(num), DIMENSION(:), ALLOCATABLE :: e1, e2, b1, b2
    INTEGER :: mm, ibc, icell

    ! Note: ideally e1, e2, b1, b2 should be face-centred. However, this is not
    ! possible with 'open' boundaries since E-fields are not defined in the
    ! ghost cell, so we use the cell-centred quantities in the first cell.

    dir = 1.0_num

    SELECT CASE(bd)
      CASE(c_bd_x_min, c_bd_x_max)
        mm = ny
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        dd = dy
        IF (bd .EQ. c_bd_x_max) THEN
          dir = -1.0_num
          ibc = nx
        ENDIF

        e1 = 0.5_num  * (ey(ibc  , 0:ny-1) + ey(ibc, 1:ny  ))
        e2 = ez(ibc, 1:ny)
        b1 = 0.25_num * (bz(ibc-1, 0:ny-1) + bz(ibc, 0:ny-1) &
                       + bz(ibc-1, 1:ny  ) + bz(ibc, 1:ny  ))
        b2 = 0.5_num  * (by(ibc-1, 1:ny  ) + by(ibc, 1:ny  ))

      CASE(c_bd_y_min, c_bd_y_max)
        mm = nx
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        dd = dx
        IF (bd .EQ. c_bd_y_max) THEN
          dir = -1.0_num
          ibc = ny
        ENDIF

        e1 = ez(1:nx, ibc)
        e2 = 0.5_num  * (ez(0:nx-1, ibc  ) + ex(1:nx  , ibc))
        b1 = 0.5_num  * (bx(1:nx  , ibc-1) + bx(1:nx  , ibc))
        b2 = 0.25_num * (bz(0:nx-1, ibc-1) + bz(0:nx-1, ibc) &
                       + bz(1:nx  , ibc-1) + bz(1:nx  , ibc))

      CASE DEFAULT
        dd = 0.0_num
        e1 = 0.0_num
        e2 = 0.0_num
        b1 = 0.0_num
        b2 = 0.0_num
    END SELECT

    laser_absorb_local = laser_absorb_local &
        + dt * dir * SUM(e1 * b1 - e2 * b2) / mu0

    IF (PRESENT(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        t_env = laser_time_profile(current)
        DO icell = 1, mm
          laser_inject_local = laser_inject_local &
              + dir * dt * 0.5_num * epsilon0 * c &
              * (t_env * current%amp * current%profile(icell))**2
        ENDDO
        current => current%next
      ENDDO
    ENDIF

    DEALLOCATE(e1, e2, b1, b2)

  END SUBROUTINE calc_absorption

END MODULE laser
