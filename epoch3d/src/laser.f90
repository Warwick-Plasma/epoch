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

    CALL allocate_with_boundary(laser%profile, boundary)
    CALL allocate_with_boundary(laser%phase, boundary)

    laser%boundary = boundary
    NULLIFY(laser%next)

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

    IF (laser%k .EQ. 0) laser%k = laser%freq

!!$    IF (boundary .EQ. c_bd_x_min .OR. boundary .EQ. c_bd_x_max) THEN
!!$      laser%phase(1:ny, 1:nz) = &
!!$          laser%phase(1:ny) - laser%k * (y(1:ny) * TAN(laser%angle))
!!$    ELSE IF (boundary .EQ. c_bd_y_max .OR. boundary .EQ. c_bd_y_min) THEN
!!$      laser%phase(1:nx, 1:nz) = &
!!$          laser%phase(1:nx) - laser%k * (x(1:nx) * TAN(laser%angle))
!!$    ENDIF

    IF (boundary .EQ. c_bd_x_min) THEN
      CALL attach_laser_to_list(laser_x_min, laser, boundary)
    ENDIF
    IF (boundary .EQ. c_bd_x_max) THEN
      CALL attach_laser_to_list(laser_x_max, laser, boundary)
    ENDIF
    IF (boundary .EQ. c_bd_y_max) THEN
      CALL attach_laser_to_list(laser_y_max, laser, boundary)
    ENDIF
    IF (boundary .EQ. c_bd_y_min) THEN
      CALL attach_laser_to_list(laser_y_min, laser, boundary)
    ENDIF
    IF (boundary .EQ. c_bd_z_max) THEN
      CALL attach_laser_to_list(laser_z_max, laser, boundary)
    ENDIF
    IF (boundary .EQ. c_bd_z_min) THEN
      CALL attach_laser_to_list(laser_z_min, laser, boundary)
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
  SUBROUTINE attach_laser_to_list(list, laser, boundary)

    TYPE(laser_block), POINTER :: list
    TYPE(laser_block), POINTER :: laser
    INTEGER, INTENT(IN) :: boundary
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
    ELSE IF (boundary .EQ. c_bd_y_max .OR. boundary .EQ. c_bd_y_min) THEN
      ALLOCATE(array(-2:nx+3, -2:nz+3))
    ELSE IF (boundary .EQ. c_bd_z_max .OR. boundary .EQ. c_bd_z_min) THEN
      ALLOCATE(array(-2:nx+3, -2:ny+2))
    ENDIF

  END SUBROUTINE allocate_with_boundary



  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = 1000000.0_num
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num*pi/current%freq
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO
    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num*pi/current%freq
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO
    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num*pi/current%freq
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO
    current=>laser_y_min
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
  SUBROUTINE laser_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: lx
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus

    TYPE(laser_block), POINTER :: current

    lx = dt/dx
    ALLOCATE(fplus(1:ny, 1:nz))
    fplus = 0.0_num
    bx(0, 1:ny, 1:nz) =  0.0_num

    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:ny, 1:nz) = fplus(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%freq*time + current%phase(1:ny, 1:nz)) &
            * SIN(current%pol_angle) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    by(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (-4.0_num * fplus + 2.0_num * ez(0, 1:ny, 1:nz) &
        - (c - lx*c**2)*by(0, 1:ny, 1:nz) - (dt / epsilon0) * jz(0, 1:ny, 1:nz))

    fplus = 0.0_num
    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:ny, 1:nz) = fplus(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%freq*time + current%phase(1:ny, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (4.0_num * fplus - 2.0_num * ey(0, 1:ny, 1:nz) &
        - (c - lx*c**2)*bz(0, 1:ny, 1:nz) + (dt / epsilon0) * jy(0, 1:ny, 1:nz))
    DEALLOCATE(fplus)
    ! CALL bfield_bcs

  END SUBROUTINE laser_bcs_x_min



  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: lx

    lx = dt/dx
    bx(0, 1:ny, 1:nz) =  0.0_num
    ! Set the y magnetic field
    by(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (2.0_num * ez(0, 1:ny, 1:nz) - (c - lx*c**2)*by(0, 1:ny, 1:nz) &
        - (dt / epsilon0) * jz(0, 1:ny, 1:nz))
    bz(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (-2.0_num * ey(0, 1:ny, 1:nz) - (c - lx*c**2)*bz(0, 1:ny, 1:nz) &
        + (dt / epsilon0) * jy(0, 1:ny, 1:nz))
    ! CALL bfield_bcs

  END SUBROUTINE outflow_bcs_x_min



  ! laser boundary for the right boundary
  SUBROUTINE laser_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: lx
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: f_minus

    TYPE(laser_block), POINTER :: current

    lx = dt/dx
    ALLOCATE(f_minus(1:ny, 1:nz))
    f_minus = 0.0_num
    bx(nx, 1:ny, 1:nz) =  0.0_num

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        f_minus(1:ny, 1:nz) = f_minus(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%freq*time + current%phase(1:ny, 1:nz)) &
            * SIN(current%pol_angle) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    by(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (4.0_num * f_minus - 2.0_num * ez(nx, 1:ny, 1:nz) &
        - (c - lx*c**2)*by(nx, 1:ny, 1:nz) + (dt / epsilon0) &
        * jz(nx, 1:ny, 1:nz))

    f_minus = 0.0_num
    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        f_minus(1:ny, 1:nz) = f_minus(1:ny, 1:nz) &
            + t_env * current%amp * current%profile(1:ny, 1:nz) &
            * SIN(current%freq*time + current%phase(1:ny, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (-4.0_num * f_minus + 2.0_num * ey(nx, 1:ny, 1:nz) &
        - (c - lx*c**2)*bz(nx, 1:ny, 1:nz) - (dt / epsilon0) &
        * jy(nx, 1:ny, 1:nz))

    DEALLOCATE(f_minus)

  END SUBROUTINE laser_bcs_x_max



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: lx

    lx = dt/dx
    bx(nx, 1:ny, 1:nz) =  0.0_num
    ! Set the y magnetic field
    by(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (- 2.0_num * ez(nx, 1:ny, 1:nz) - (c - lx*c**2)*by(nx, 1:ny, 1:nz) &
        + (dt / epsilon0) * jz(nx, 1:ny, 1:nz))
    bz(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) &
        * (2.0_num * ey(nx, 1:ny, 1:nz) - (c - lx*c**2)*bz(nx, 1:ny, 1:nz) &
        - (dt / epsilon0) * jy(nx, 1:ny, 1:nz))

  END SUBROUTINE outflow_bcs_x_max



  ! laser boundary for the bottom boundary
  SUBROUTINE laser_bcs_y_min

    REAL(num) :: t_env
    REAL(num) :: ly
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus

    TYPE(laser_block), POINTER :: current

    ly = dt/dy
    ALLOCATE(fplus(1:nx, 1:nz))
    fplus = 0.0_num
    by(1:nx, -1, 1:nz) =  0.0_num

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%freq*time + current%phase(1:nx, 1:nz)) &
            * SIN(current%pol_angle) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (-4.0_num * fplus - 2.0_num * ez(1:nx, 0, 1:nz) &
        - (c - ly*c**2)*bx(1:nx, 0, 1:nz) - (dt / epsilon0) * jz(1:nx, 0, 1:nz))

    fplus = 0.0_num
    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%freq*time + current%phase(1:nx, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (4.0_num * fplus + 2.0_num * ex(1:nx, 0, 1:nz) &
        - (c - ly*c**2)*bz(1:nx, 0, 1:nz) + (dt / epsilon0) * jx(1:nx, 0, 1:nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_y_min



  SUBROUTINE outflow_bcs_y_min

    REAL(num) :: ly

    ly = dt/dy
    ! Set the x magnetic field
    bx(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (-2.0_num * ez(1:nx, 0, 1:nz) - (c - ly*c**2)*bx(1:nx, 0, 1:nz) &
        - (dt / epsilon0) * jz(1:nx, 0, 1:nz))
    by(1:nx, 0, 1:nz) =  0.0_num
    bz(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (2.0_num * ex(1:nx, 0, 1:nz) - (c - ly*c**2)*bz(1:nx, 0, 1:nz) &
        + (dt / epsilon0) * jx(1:nx, 0, 1:nz))

  END SUBROUTINE outflow_bcs_y_min



  ! laser boundary for the top boundary
  SUBROUTINE laser_bcs_y_max

    REAL(num) :: t_env
    REAL(num) :: ly
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus

    TYPE(laser_block), POINTER :: current

    ly = dt/dy
    ALLOCATE(fplus(1:nx, 1:nz))
    fplus = 0.0_num
    by(1:nx, ny, 1:nz) =  0.0_num

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%freq*time + current%phase(1:nx, 1:nz)) &
            * SIN(current%pol_angle) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (-4.0_num * fplus + 2.0_num * ez(1:nx, ny, 1:nz) &
        - (c - ly*c**2)*bx(1:nx, ny, 1:nz) - (dt / epsilon0) &
        * jz(1:nx, ny, 1:nz))

    fplus = 0.0_num
    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) &
            + t_env * current%amp * current%profile(1:nx, 1:nz) &
            * SIN(current%freq*time + current%phase(1:nx, 1:nz)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (4.0_num * fplus - 2.0_num * ex(1:nx, ny, 1:nz) &
        - (c - ly*c**2)*bz(1:nx, ny, 1:nz) + (dt / epsilon0) &
        * jx(1:nx, ny, 1:nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_y_max



  SUBROUTINE outflow_bcs_y_max

    REAL(num) :: ly

    ly = dt/dy
    ! Set the x magnetic field
    bx(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (2.0_num * ez(1:nx, ny, 1:nz) - (c - ly*c**2)*bx(1:nx, ny, 1:nz) &
        - (dt / epsilon0) * jz(1:nx, ny, 1:nz))
    by(1:nx, ny, 1:nz) =  0.0_num
    bz(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) &
        * (-2.0_num * ex(1:nx, ny, 1:nz) - (c - ly*c**2)*bz(1:nx, ny, 1:nz) &
        + (dt / epsilon0) * jx(1:nx, ny, 1:nz))

  END SUBROUTINE outflow_bcs_y_max



  ! laser boundary for the back boundary
  SUBROUTINE laser_bcs_z_min

    REAL(num) :: t_env
    REAL(num) :: lz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus

    TYPE(laser_block), POINTER :: current

    lz = dt/dz
    ALLOCATE(fplus(1:nx, 1:ny))
    fplus = 0.0_num
    by(1:nx, 1:ny, 0) =  0.0_num

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%freq*time + current%phase(1:nx, 1:ny)) &
            * SIN(current%pol_angle) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) &
        * (-4.0_num * fplus + 2.0_num * ez(1:nx, 1:ny, 0) &
        - (c - lz*c**2)*bx(1:nx, 1:ny, 0) - (dt / epsilon0) * jz(1:nx, 1:ny, 0))

    fplus = 0.0_num
    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%freq*time + current%phase(1:nx, 1:ny)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) &
        * (4.0_num * fplus - 2.0_num * ex(1:nx, 1:ny, 0) &
        - (c - lz*c**2)*bz(1:nx, 1:ny, 0) + (dt / epsilon0) * jx(1:nx, 1:ny, 0))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_z_min



  SUBROUTINE outflow_bcs_z_min

    REAL(num) :: lz

    lz = dt/dz
    ! Set the x magnetic field
    bx(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) &
        * (2.0_num * ez(1:nx, 1:ny, 0) - (c - lz*c**2)*bx(1:nx, 1:ny, 0) &
        - (dt / epsilon0) * jz(1:nx, 1:ny, 0))
    by(1:nx, 1:ny, 0) =  0.0_num
    bz(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) &
        * (-2.0_num * ex(1:nx, 1:ny, 0) - (c - lz*c**2)*bz(1:nx, 1:ny, 0) &
        + (dt / epsilon0) * jx(1:nx, 1:ny, 0))

  END SUBROUTINE outflow_bcs_z_min



  ! laser boundary for the back boundary
  SUBROUTINE laser_bcs_z_max

    REAL(num) :: t_env
    REAL(num) :: lz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus

    TYPE(laser_block), POINTER :: current

    lz = dt/dz
    ALLOCATE(fplus(1:nx, 1:ny))
    fplus = 0.0_num
    by(1:nx, 1:ny, nz) =  0.0_num

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%freq*time + current%phase(1:nx, 1:ny)) &
            * SIN(current%pol_angle) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) &
        * (-4.0_num * fplus + 2.0_num * ez(1:nx, 1:ny, nz) &
        - (c - lz*c**2)*bx(1:nx, 1:ny, nz) - (dt / epsilon0) &
        * jz(1:nx, 1:ny, nz))

    fplus = 0.0_num
    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) &
            + t_env * current%amp * current%profile(1:nx, 1:ny) &
            * SIN(current%freq*time + current%phase(1:nx, 1:ny)) &
            * COS(current%pol_angle)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) &
        * (4.0_num * fplus - 2.0_num * ex(1:nx, 1:ny, nz) &
        - (c - lz*c**2)*bz(1:nx, 1:ny, nz) + (dt / epsilon0) &
        * jx(1:nx, 1:ny, nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_z_max



  SUBROUTINE outflow_bcs_z_max

    REAL(num) :: lz

    lz = dt/dz
    ! Set the x magnetic field
    bx(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) &
        * (2.0_num * ez(1:nx, 1:ny, nz) - (c - lz*c**2)*bx(1:nx, 1:ny, nz) &
        - (dt / epsilon0) * jz(1:nx, 1:ny, nz))
    by(1:nx, 1:ny, nz) =  0.0_num
    bz(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) &
        * (-2.0_num * ex(1:nx, 1:ny, nz) - (c - lz*c**2)*bz(1:nx, 1:ny, nz) &
        + (dt / epsilon0) * jx(1:nx, 1:ny, nz))

  END SUBROUTINE outflow_bcs_z_max

END MODULE laser
