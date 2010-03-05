MODULE laser

  USE evaluator
  USE custom_laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_laser(direction, laser)

    INTEGER, INTENT(IN) :: direction
    TYPE(laser_block), INTENT(INOUT) :: laser

    CALL allocate_with_direction(laser%profile, direction)
    CALL allocate_with_direction(laser%phase, direction)

    laser%direction = direction
    NULLIFY(laser%next)

  END SUBROUTINE init_laser



  SUBROUTINE delete_laser(laser)

    TYPE(laser_block), INTENT(INOUT) :: laser

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)

  END SUBROUTINE delete_laser



  ! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    INTEGER :: direction
    TYPE(laser_block), POINTER :: laser

    direction = laser%direction

    IF (laser%k .EQ. 0) laser%k = laser%freq

!!$    IF (direction .EQ. c_bd_left .OR. direction .EQ. c_bd_right) THEN
!!$      laser%phase(1:ny, 1:nz) = &
!!$          laser%phase(1:ny) - laser%k * (y(1:ny) * TAN(laser%angle))
!!$    ELSE IF (direction .EQ. c_bd_up .OR. direction .EQ. c_bd_down) THEN
!!$      laser%phase(1:nx, 1:nz) = &
!!$          laser%phase(1:nx) - laser%k * (x(1:nx) * TAN(laser%angle))
!!$    ENDIF

    IF (direction .EQ. c_bd_left) THEN
      CALL attach_laser_to_list(laser_left, laser, direction)
    ENDIF
    IF (direction .EQ. c_bd_right) THEN
      CALL attach_laser_to_list(laser_right, laser, direction)
    ENDIF
    IF (direction .EQ. c_bd_up) THEN
      CALL attach_laser_to_list(laser_up, laser, direction)
    ENDIF
    IF (direction .EQ. c_bd_down) THEN
      CALL attach_laser_to_list(laser_down, laser, direction)
    ENDIF
    IF (direction .EQ. c_bd_front) THEN
      CALL attach_laser_to_list(laser_front, laser, direction)
    ENDIF
    IF (direction .EQ. c_bd_back) THEN
      CALL attach_laser_to_list(laser_back, laser, direction)
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



  SUBROUTINE allocate_with_direction(array, direction)

    REAL(num), DIMENSION(:,:), POINTER :: array
    INTEGER, INTENT(IN) :: direction

    IF (direction .EQ. c_bd_left .OR. direction .EQ. c_bd_right) THEN
      ALLOCATE(array(-2:ny+3, -2:nz+3))
    ELSE IF (direction .EQ. c_bd_up .OR. direction .EQ. c_bd_down) THEN
      ALLOCATE(array(-2:nx+3, -2:nz+3))
    ELSE IF (direction .EQ. c_bd_front .OR. direction .EQ. c_bd_back) THEN
      ALLOCATE(array(-2:nx+3, -2:ny+2))
    ENDIF

  END SUBROUTINE allocate_with_direction



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
    current=>laser_up
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num*pi/current%freq
      dt_laser = MIN(dt_laser, dt_local)
      current=>current%next
    ENDDO
    current=>laser_down
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
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    INTEGER :: err

    TYPE(laser_block), POINTER :: current

    lx = dt/dx
    ALLOCATE(fplus(1:ny, 1:nz))
    fplus = 0.0_num
    err = 0
    bx(0, 1:ny, 1:nz) =  0.0_num

    current=>laser_left
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:ny, 1:nz) = fplus(1:ny, 1:nz) + &
            t_env * current%amp * current%profile(1:ny, 1:nz) * &
            SIN(current%freq*time + current%phase(1:ny, 1:nz)) * &
            SIN(current%pol) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    by(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (-4.0_num * fplus + 2.0_num * ez(0, 1:ny, 1:nz) - &
        (c - lx*c**2)*by(0, 1:ny, 1:nz) - (dt / epsilon0) * jz(0, 1:ny, 1:nz))

    fplus = 0.0_num
    current=>laser_left
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:ny, 1:nz) = fplus(1:ny, 1:nz) + &
            t_env * current%amp * current%profile(1:ny, 1:nz) * &
            SIN(current%freq*time + current%phase(1:ny, 1:nz)) * &
            COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (4.0_num * fplus - 2.0_num * ey(0, 1:ny, 1:nz) - &
        (c - lx*c**2)*bz(0, 1:ny, 1:nz) + (dt / epsilon0) * jy(0, 1:ny, 1:nz))
    DEALLOCATE(fplus)
    ! CALL bfield_bcs

  END SUBROUTINE laser_bcs_left



  SUBROUTINE outflow_bcs_left

    REAL(num) :: lx

    lx = dt/dx
    bx(0, 1:ny, 1:nz) =  0.0_num
    ! Set the y magnetic field
    by(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (2.0_num * ez(0, 1:ny, 1:nz) - (c - lx*c**2)*by(0, 1:ny, 1:nz) - &
        (dt / epsilon0) * jz(0, 1:ny, 1:nz))
    bz(0, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (-2.0_num * ey(0, 1:ny, 1:nz) - (c - lx*c**2)*bz(0, 1:ny, 1:nz) + &
        (dt / epsilon0) * jy(0, 1:ny, 1:nz))
    ! CALL bfield_bcs

  END SUBROUTINE outflow_bcs_left



  ! laser boundary for the right boundary
  SUBROUTINE laser_bcs_right

    REAL(num) :: t_env
    REAL(num) :: lx
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: f_minus
    INTEGER :: err

    TYPE(laser_block), POINTER :: current

    lx = dt/dx
    ALLOCATE(f_minus(1:ny, 1:nz))
    f_minus = 0.0_num
    err = 0
    bx(nx, 1:ny, 1:nz) =  0.0_num

    current=>laser_right
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        f_minus(1:ny, 1:nz) = f_minus(1:ny, 1:nz) + &
            t_env * current%amp * current%profile(1:ny, 1:nz) * &
            SIN(current%freq*time + current%phase(1:ny, 1:nz)) * &
            SIN(current%pol) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    by(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (4.0_num * f_minus - 2.0_num * ez(nx, 1:ny, 1:nz) - &
        (c - lx*c**2)*by(nx, 1:ny, 1:nz) + (dt / epsilon0) * jz(nx, 1:ny, 1:nz))

    f_minus = 0.0_num
    current=>laser_right
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        f_minus(1:ny, 1:nz) = f_minus(1:ny, 1:nz) + &
            t_env * current%amp * current%profile(1:ny, 1:nz) * &
            SIN(current%freq*time + current%phase(1:ny, 1:nz)) * &
            COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (-4.0_num * f_minus + 2.0_num * ey(nx, 1:ny, 1:nz) - &
        (c - lx*c**2)*bz(nx, 1:ny, 1:nz) - (dt / epsilon0) * jy(nx, 1:ny, 1:nz))

    DEALLOCATE(f_minus)

  END SUBROUTINE laser_bcs_right



  SUBROUTINE outflow_bcs_right

    REAL(num) :: lx

    lx = dt/dx
    bx(nx, 1:ny, 1:nz) =  0.0_num
    ! Set the y magnetic field
    by(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (- 2.0_num * ez(nx, 1:ny, 1:nz) - (c - lx*c**2)*by(nx, 1:ny, 1:nz) + &
        (dt / epsilon0) * jz(nx, 1:ny, 1:nz))
    bz(nx, 1:ny, 1:nz) = (1.0_num / (c + lx*c**2)) * &
        (2.0_num * ey(nx, 1:ny, 1:nz) - (c - lx*c**2)*bz(nx, 1:ny, 1:nz) - &
        (dt / epsilon0) * jy(nx, 1:ny, 1:nz))

  END SUBROUTINE outflow_bcs_right



  ! laser boundary for the bottom boundary
  SUBROUTINE laser_bcs_down

    REAL(num) :: t_env
    REAL(num) :: ly
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    INTEGER :: err

    TYPE(laser_block), POINTER :: current

    ly = dt/dy
    ALLOCATE(fplus(1:nx, 1:nz))
    fplus = 0.0_num
    err = 0
    by(1:nx, -1, 1:nz) =  0.0_num

    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) + &
            t_env * current%amp * current%profile(1:nx, 1:nz) * &
            SIN(current%freq*time + current%phase(1:nx, 1:nz)) * &
            SIN(current%pol) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (-4.0_num * fplus - 2.0_num * ez(1:nx, 0, 1:nz) - &
        (c - ly*c**2)*bx(1:nx, 0, 1:nz) - (dt / epsilon0) * jz(1:nx, 0, 1:nz))

    fplus = 0.0_num
    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) + &
            t_env * current%amp * current%profile(1:nx, 1:nz) * &
            SIN(current%freq*time + current%phase(1:nx, 1:nz)) * &
            COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (4.0_num * fplus + 2.0_num * ex(1:nx, 0, 1:nz) - &
        (c - ly*c**2)*bz(1:nx, 0, 1:nz) + (dt / epsilon0) * jx(1:nx, 0, 1:nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_down



  SUBROUTINE outflow_bcs_down

    REAL(num) :: ly

    ly = dt/dy
    ! Set the x magnetic field
    bx(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (-2.0_num * ez(1:nx, 0, 1:nz) - (c - ly*c**2)*bx(1:nx, 0, 1:nz) - &
        (dt / epsilon0) * jz(1:nx, 0, 1:nz))
    by(1:nx, 0, 1:nz) =  0.0_num
    bz(1:nx, 0, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (2.0_num * ex(1:nx, 0, 1:nz) - (c - ly*c**2)*bz(1:nx, 0, 1:nz) + &
        (dt / epsilon0) * jx(1:nx, 0, 1:nz))

  END SUBROUTINE outflow_bcs_down



  ! laser boundary for the top boundary
  SUBROUTINE laser_bcs_up

    REAL(num) :: t_env
    REAL(num) :: ly
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    INTEGER :: err

    TYPE(laser_block), POINTER :: current

    ly = dt/dy
    ALLOCATE(fplus(1:nx, 1:nz))
    fplus = 0.0_num
    err = 0
    by(1:nx, ny, 1:nz) =  0.0_num

    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) + &
            t_env * current%amp * current%profile(1:nx, 1:nz) * &
            SIN(current%freq*time + current%phase(1:nx, 1:nz)) * &
            SIN(current%pol) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (-4.0_num * fplus + 2.0_num * ez(1:nx, ny, 1:nz) - &
        (c - ly*c**2)*bx(1:nx, ny, 1:nz) - (dt / epsilon0) * jz(1:nx, ny, 1:nz))

    fplus = 0.0_num
    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:nz) = fplus(1:nx, 1:nz) + &
            t_env * current%amp * current%profile(1:nx, 1:nz) * &
            SIN(current%freq*time + current%phase(1:nx, 1:nz)) * &
            COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (4.0_num * fplus - 2.0_num * ex(1:nx, ny, 1:nz) - &
        (c - ly*c**2)*bz(1:nx, ny, 1:nz) + (dt / epsilon0) * jx(1:nx, ny, 1:nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_up



  SUBROUTINE outflow_bcs_up

    REAL(num) :: ly

    ly = dt/dy
    ! Set the x magnetic field
    bx(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (2.0_num * ez(1:nx, ny, 1:nz) - (c - ly*c**2)*bx(1:nx, ny, 1:nz) - &
        (dt / epsilon0) * jz(1:nx, ny, 1:nz))
    by(1:nx, ny, 1:nz) =  0.0_num
    bz(1:nx, ny, 1:nz) = (1.0_num / (c + ly*c**2)) * &
        (-2.0_num * ex(1:nx, ny, 1:nz) - (c - ly*c**2)*bz(1:nx, ny, 1:nz) + &
        (dt / epsilon0) * jx(1:nx, ny, 1:nz))

  END SUBROUTINE outflow_bcs_up



  ! laser boundary for the back boundary
  SUBROUTINE laser_bcs_back

    REAL(num) :: t_env
    REAL(num) :: lz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    INTEGER :: err

    TYPE(laser_block), POINTER :: current

    lz = dt/dz
    ALLOCATE(fplus(1:nx, 1:ny))
    fplus = 0.0_num
    err = 0
    by(1:nx, 1:ny, 0) =  0.0_num

    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) + &
            t_env * current%amp * current%profile(1:nx, 1:ny) * &
            SIN(current%freq*time + current%phase(1:nx, 1:ny)) * &
            SIN(current%pol) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) * &
        (-4.0_num * fplus + 2.0_num * ez(1:nx, 1:ny, 0) - &
        (c - lz*c**2)*bx(1:nx, 1:ny, 0) - (dt / epsilon0) * jz(1:nx, 1:ny, 0))

    fplus = 0.0_num
    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) + &
            t_env * current%amp * current%profile(1:nx, 1:ny) * &
            SIN(current%freq*time + current%phase(1:nx, 1:ny)) * &
            COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) * &
        (4.0_num * fplus - 2.0_num * ex(1:nx, 1:ny, 0) - &
        (c - lz*c**2)*bz(1:nx, 1:ny, 0) + (dt / epsilon0) * jx(1:nx, 1:ny, 0))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_back



  SUBROUTINE outflow_bcs_back

    REAL(num) :: lz

    lz = dt/dz
    ! Set the x magnetic field
    bx(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) * &
        (2.0_num * ez(1:nx, 1:ny, 0) - (c - lz*c**2)*bx(1:nx, 1:ny, 0) - &
        (dt / epsilon0) * jz(1:nx, 1:ny, 0))
    by(1:nx, 1:ny, 0) =  0.0_num
    bz(1:nx, 1:ny, 0) = (1.0_num / (c + lz*c**2)) * &
        (-2.0_num * ex(1:nx, 1:ny, 0) - (c - lz*c**2)*bz(1:nx, 1:ny, 0) + &
        (dt / epsilon0) * jx(1:nx, 1:ny, 0))

  END SUBROUTINE outflow_bcs_back



  ! laser boundary for the back boundary
  SUBROUTINE laser_bcs_front

    REAL(num) :: t_env
    REAL(num) :: lz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fplus
    INTEGER :: err

    TYPE(laser_block), POINTER :: current

    lz = dt/dz
    ALLOCATE(fplus(1:nx, 1:ny))
    fplus = 0.0_num
    err = 0
    by(1:nx, 1:ny, nz) =  0.0_num

    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) + &
            t_env * current%amp * current%profile(1:nx, 1:ny) * &
            SIN(current%freq*time + current%phase(1:nx, 1:ny)) * &
            SIN(current%pol) *COS(current%angle)
      ENDIF
      current=>current%next
    ENDDO

    ! Set the y magnetic field
    bx(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) * &
        (-4.0_num * fplus + 2.0_num * ez(1:nx, 1:ny, nz) - &
        (c - lz*c**2)*bx(1:nx, 1:ny, nz) - (dt / epsilon0) * jz(1:nx, 1:ny, nz))

    fplus = 0.0_num
    current=>laser_down
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the laser
      IF (time .GE. current%t_start .AND. time .LE. current%t_end) THEN
        t_env = laser_time_profile(current)
        fplus(1:nx, 1:ny) = fplus(1:nx, 1:ny) + &
            t_env * current%amp * current%profile(1:nx, 1:ny) * &
            SIN(current%freq*time + current%phase(1:nx, 1:ny)) * &
            COS(current%pol)
      ENDIF
      current=>current%next
    ENDDO

    bz(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) * &
        (4.0_num * fplus - 2.0_num * ex(1:nx, 1:ny, nz) - &
        (c - lz*c**2)*bz(1:nx, 1:ny, nz) + (dt / epsilon0) * jx(1:nx, 1:ny, nz))

    DEALLOCATE(fplus)

  END SUBROUTINE laser_bcs_front



  SUBROUTINE outflow_bcs_front

    REAL(num) :: lz

    lz = dt/dz
    ! Set the x magnetic field
    bx(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) * &
        (2.0_num * ez(1:nx, 1:ny, nz) - (c - lz*c**2)*bx(1:nx, 1:ny, nz) - &
        (dt / epsilon0) * jz(1:nx, 1:ny, nz))
    by(1:nx, 1:ny, nz) =  0.0_num
    bz(1:nx, 1:ny, nz) = (1.0_num / (c + lz*c**2)) * &
        (-2.0_num * ex(1:nx, 1:ny, nz) - (c - lz*c**2)*bz(1:nx, 1:ny, nz) + &
        (dt / epsilon0) * jx(1:nx, 1:ny, nz))

  END SUBROUTINE outflow_bcs_front

END MODULE laser
