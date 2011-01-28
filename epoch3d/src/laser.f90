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
    ELSE IF (boundary .EQ. c_bd_z_min) THEN
      CALL attach_laser_to_list(laser_z_min, laser)
    ELSE IF (boundary .EQ. c_bd_z_max) THEN
      CALL attach_laser_to_list(laser_z_max, laser)
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



  SUBROUTINE allocate_with_boundary(array, boundary)

    REAL(num), DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(IN) :: boundary

    IF (boundary .EQ. c_bd_x_min .OR. boundary .EQ. c_bd_x_max) THEN
      ALLOCATE(array(-2:ny+3, -2:nz+3))
    ELSE IF (boundary .EQ. c_bd_y_min .OR. boundary .EQ. c_bd_y_max) THEN
      ALLOCATE(array(-2:nx+3, -2:nz+3))
    ELSE IF (boundary .EQ. c_bd_z_min .OR. boundary .EQ. c_bd_z_max) THEN
      ALLOCATE(array(-2:nx+3, -2:ny+3))
    ENDIF

  END SUBROUTINE allocate_with_boundary



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

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO

    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO

    current=>laser_z_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO

    current=>laser_z_max
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
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    TYPE(laser_block), POINTER :: current

    lx = c**2 * dt / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(fplus(1:ny, 1:nz))
    bx(0, 1:ny, 1:nz) = 0.0_num

    fplus = 0.0_num
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:ny, 1:nz) = fplus(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%omega * time + current%phase(1:ny, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(0, 1:ny, 1:nz) = sum * ( 4.0_num * fplus &
        - 2.0_num * ey(1, 1:ny, 1:nz) &
        + dt_eps * jy(1, 1:ny, 1:nz) &
        + diff * bz(1, 1:ny, 1:nz))

    fplus = 0.0_num
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:ny, 1:nz) = fplus(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%omega * time + current%phase(1:ny, 1:nz)) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    by(0, 1:ny, 1:nz) = sum * (-4.0_num * fplus &
        + 2.0_num * ez(1, 1:ny, 1:nz) &
        - dt_eps * jz(1, 1:ny, 1:nz) &
        + diff * by(1, 1:ny, 1:nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_x_min



  SUBROUTINE laser_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: lx, sum, diff, dt_eps
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fneg
    TYPE(laser_block), POINTER :: current

    lx = c**2 * dt / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(fneg(1:ny, 1:nz))
    bx(nx+1, 1:ny, 1:nz) = 0.0_num

    fneg = 0.0_num
    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg(1:ny, 1:nz) = fneg(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%omega * time + current%phase(1:ny, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(nx, 1:ny, 1:nz) = sum * (-4.0_num * fneg &
        + 2.0_num * ey(nx, 1:ny, 1:nz) &
        - dt_eps * jy(nx, 1:ny, 1:nz) &
        + diff * bz(nx-1, 1:ny, 1:nz))

    fneg = 0.0_num
    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg(1:ny, 1:nz) = fneg(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%omega * time + current%phase(1:ny, 1:nz)) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    by(nx, 1:ny, 1:nz) = sum * ( 4.0_num * fneg &
        - 2.0_num * ez(nx, 1:ny, 1:nz) &
        + dt_eps * jz(nx, 1:ny, 1:nz) &
        + diff * by(nx-1, 1:ny, 1:nz))

    DEALLOCATE(fneg)

  END SUBROUTINE laser_bcs_x_max



  SUBROUTINE laser_bcs_y_min

    REAL(num) :: t_env
    REAL(num) :: ly, sum, diff, dt_eps
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    TYPE(laser_block), POINTER :: current

    ly = c**2 * dt / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(fplus(1:nx, 1:nz))
    by(1:nx, 0, 1:nz) = 0.0_num

    fplus = 0.0_num
    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%omega * time + current%phase(1:nx, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bx(1:nx, 0, 1:nz) = sum * ( 4.0_num * fplus &
        - 2.0_num * ez(1:nx, 1, 1:nz) &
        + dt_eps * jz(1:nx, 1, 1:nz) &
        + diff * bx(1:nx, 1, 1:nz))

    fplus = 0.0_num
    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%omega * time + current%phase(1:nx, 1:nz)) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, 0, 1:nz) = sum * (-4.0_num * fplus &
        + 2.0_num * ex(1:nx, 1, 1:nz) &
        - dt_eps * jx(1:nx, 1, 1:nz) &
        + diff * bz(1:nx, 1, 1:nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_y_min



  SUBROUTINE laser_bcs_y_max

    REAL(num) :: t_env
    REAL(num) :: ly, sum, diff, dt_eps
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fneg
    TYPE(laser_block), POINTER :: current

    ly = c**2 * dt / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(fneg(1:nx, 1:nz))
    by(1:nx, ny+1, 1:nz) = 0.0_num

    fneg = 0.0_num
    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg(1:nx, 1:nz) = fneg(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%omega * time + current%phase(1:nx, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bx(1:nx, ny, 1:nz) = sum * (-4.0_num * fneg &
        + 2.0_num * ez(1:nx, ny, 1:nz) &
        - dt_eps * jz(1:nx, ny, 1:nz) &
        + diff * bx(1:nx, ny-1, 1:nz))

    fneg = 0.0_num
    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg(1:nx, 1:nz) = fneg(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%omega * time + current%phase(1:nx, 1:nz)) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, ny, 1:nz) = sum * ( 4.0_num * fneg &
        - 2.0_num * ex(1:nx, ny, 1:nz) &
        + dt_eps * jx(1:nx, ny, 1:nz) &
        + diff * bz(1:nx, ny-1, 1:nz))

    DEALLOCATE(fneg)

  END SUBROUTINE laser_bcs_y_max



  SUBROUTINE laser_bcs_z_min

    REAL(num) :: t_env
    REAL(num) :: lz, sum, diff, dt_eps
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    TYPE(laser_block), POINTER :: current

    lz = c**2 * dt / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    ALLOCATE(fplus(1:nx, 1:ny))
    bz(1:nx, 1:ny, 0) = 0.0_num

    fplus = 0.0_num
    current=>laser_z_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%omega * time + current%phase(1:nx, 1:ny)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    by(1:nx, 1:ny, 0) = sum * ( 4.0_num * fplus &
        - 2.0_num * ex(1:nx, 1:ny, 1) &
        + dt_eps * jx(1:nx, 1:ny, 1) &
        + diff * by(1:nx, 1:ny, 1))

    fplus = 0.0_num
    current=>laser_z_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%omega * time + current%phase(1:nx, 1:ny)) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bx(1:nx, 1:ny, 0) = sum * (-4.0_num * fplus &
        + 2.0_num * ey(1:nx, 1:ny, 1) &
        - dt_eps * jy(1:nx, 1:ny, 1) &
        + diff * bx(1:nx, 1:ny, 1))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_z_min



  SUBROUTINE laser_bcs_z_max

    REAL(num) :: t_env
    REAL(num) :: lz, sum, diff, dt_eps
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fneg
    TYPE(laser_block), POINTER :: current

    lz = c**2 * dt / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    ALLOCATE(fneg(1:nx, 1:ny))
    bz(1:nx, 1:ny, nz+1) = 0.0_num

    fneg = 0.0_num
    current=>laser_z_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg(1:nx, 1:ny) = fneg(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%omega * time + current%phase(1:nx, 1:ny)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    by(1:nx, 1:ny, nz) = sum * (-4.0_num * fneg &
        + 2.0_num * ex(1:nx, 1:ny, nz) &
        - dt_eps * jx(1:nx, 1:ny, nz) &
        + diff * by(1:nx, 1:ny, nz-1))

    fneg = 0.0_num
    current=>laser_z_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fneg(1:nx, 1:ny) = fneg(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%omega * time + current%phase(1:nx, 1:ny)) &
            * SIN(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bx(1:nx, 1:ny, nz) = sum * ( 4.0_num * fneg &
        - 2.0_num * ey(1:nx, 1:ny, nz) &
        + dt_eps * jy(1:nx, 1:ny, nz) &
        + diff * bx(1:nx, 1:ny, nz-1))

    DEALLOCATE(fneg)

  END SUBROUTINE laser_bcs_z_max



  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: lx, sum, diff, dt_eps

    lx = c**2 * dt / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    bx(0, 1:ny, 1:nz) = 0.0_num
    bz(0, 1:ny, 1:nz) = sum * (-2.0_num * ey(1, 1:ny, 1:nz) &
        + dt_eps * jy(1, 1:ny, 1:nz) + diff * bz(1, 1:ny, 1:nz))
    by(0, 1:ny, 1:nz) = sum * ( 2.0_num * ez(1, 1:ny, 1:nz) &
        - dt_eps * jz(1, 1:ny, 1:nz) + diff * by(1, 1:ny, 1:nz))

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: lx, sum, diff, dt_eps

    lx = c**2 * dt / dx
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    bx(nx+1, 1:ny, 1:nz) = 0.0_num
    bz(nx, 1:ny, 1:nz) = sum * ( 2.0_num * ey(nx, 1:ny, 1:nz) &
        - dt_eps * jy(nx, 1:ny, 1:nz) + diff * bz(nx-1, 1:ny, 1:nz))
    by(nx, 1:ny, 1:nz) = sum * (-2.0_num * ez(nx, 1:ny, 1:nz) &
        + dt_eps * jz(nx, 1:ny, 1:nz) + diff * by(nx-1, 1:ny, 1:nz))

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE outflow_bcs_y_min

    REAL(num) :: ly, sum, diff, dt_eps

    ly = c**2 * dt / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    by(1:nx, 0, 1:nz) = 0.0_num
    bx(1:nx, 0, 1:nz) = sum * (-2.0_num * ez(1:nx, 1, 1:nz) &
        + dt_eps * jz(1:nx, 1, 1:nz) + diff * bx(1:nx, 1, 1:nz))
    bz(1:nx, 0, 1:nz) = sum * ( 2.0_num * ex(1:nx, 1, 1:nz) &
        - dt_eps * jx(1:nx, 1, 1:nz) + diff * bz(1:nx, 1, 1:nz))

  END SUBROUTINE outflow_bcs_y_min



  SUBROUTINE outflow_bcs_y_max

    REAL(num) :: ly, sum, diff, dt_eps

    ly = c**2 * dt / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    by(1:nx, ny+1, 1:nz) = 0.0_num
    bx(1:nx, ny, 1:nz) = sum * ( 2.0_num * ez(1:nx, ny, 1:nz) &
        - dt_eps * jz(1:nx, ny, 1:nz) + diff * bx(1:nx, ny-1, 1:nz))
    bz(1:nx, ny, 1:nz) = sum * (-2.0_num * ex(1:nx, ny, 1:nz) &
        + dt_eps * jx(1:nx, ny, 1:nz) + diff * bz(1:nx, ny-1, 1:nz))

  END SUBROUTINE outflow_bcs_y_max



  SUBROUTINE outflow_bcs_z_min

    REAL(num) :: lz, sum, diff, dt_eps

    lz = c**2 * dt / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    bz(1:nx, 1:ny, 0) = 0.0_num
    by(1:nx, 1:ny, 0) = sum * (-2.0_num * ex(1:nx, 1:ny, 1) &
        + dt_eps * jx(1:nx, 1:ny, 1) + diff * by(1:nx, 1:ny, 1))
    bx(1:nx, 1:ny, 0) = sum * ( 2.0_num * ey(1:nx, 1:ny, 1) &
        - dt_eps * jy(1:nx, 1:ny, 1) + diff * bx(1:nx, 1:ny, 1))

  END SUBROUTINE outflow_bcs_z_min



  SUBROUTINE outflow_bcs_z_max

    REAL(num) :: lz, sum, diff, dt_eps

    lz = c**2 * dt / dz
    sum = 1.0_num / (lz + c)
    diff = lz - c
    dt_eps = dt / epsilon0

    bz(1:nx, 1:ny, nz+1) = 0.0_num
    by(1:nx, 1:ny, nz) = sum * ( 2.0_num * ex(1:nx, 1:ny, nz) &
        - dt_eps * jx(1:nx, 1:ny, nz) + diff * by(1:nx, 1:ny, nz-1))
    bx(1:nx, 1:ny, nz) = sum * (-2.0_num * ey(1:nx, 1:ny, nz) &
        + dt_eps * jy(1:nx, 1:ny, nz) + diff * bx(1:nx, 1:ny, nz-1))

  END SUBROUTINE outflow_bcs_z_max

END MODULE laser
