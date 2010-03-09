MODULE fields

  USE boundary
  USE laser

  IMPLICIT NONE

  REAL(num) :: alpha = 0.0_num

CONTAINS

  SUBROUTINE update_eb_fields_half

    REAL(num) :: lx, cnx, ly, cny, lz, cnz

    lx = dt/dx
    cnx = 0.5_num*lx

    ly = dt/dy
    cny = 0.5_num*ly

    lz = dt/dz
    cnz = 0.5_num*lz

    ! Update ex to t = t0+dt/2

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ex(ix, iy, iz) = ex(ix, iy, iz) &
              + cny*(bz(ix, iy, iz)-bz(ix, iy-1, iz))*c**2 &
              - cnz*(by(ix, iy, iz)-by(ix, iy, iz-1))*c**2 &
              - 0.5_num*dt*jx(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! Update ey to t = t0+dt/2
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ey(ix, iy, iz) = ey(ix, iy, iz) &
              + cnz*(bx(ix, iy, iz)-bx(ix, iy, iz-1))*c**2 &
              - cnx*(bz(ix, iy, iz)-bz(ix-1, iy, iz))*c**2 &
              - 0.5_num*dt*jy(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! Update ez to t = t0+dt/2
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ez(ix, iy, iz) = ez(ix, iy, iz) &
              + cnx*(by(ix, iy, iz)-by(ix-1, iy, iz))*c**2 &
              - cny*(bx(ix, iy, iz)-bx(ix, iy-1, iz))*c**2 &
              - 0.5_num*dt*jz(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! Now have E(t+dt/2), do boundary conditions on E

    CALL efield_bcs

    ! Update B field to t+dt/2 using E(t+dt/2)

    ! bx
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          bx(ix, iy, iz) = bx(ix, iy, iz) &
              + cnz*(ey(ix, iy, iz+1)-ey(ix, iy, iz)) &
              - cny*(ez(ix, iy+1, iz)-ez(ix, iy, iz))
        ENDDO
      ENDDO
    ENDDO

    ! by
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          by(ix, iy, iz) = by(ix, iy, iz) &
              + cnx*(ez(ix+1, iy, iz)-ez(ix, iy, iz)) &
              - cnz*(ex(ix, iy, iz+1)-ex(ix, iy, iz))
        ENDDO
      ENDDO
    ENDDO

    ! bz
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          bz(ix, iy, iz) = bz(ix, iy, iz) &
              - cnx*(ey(ix+1, iy, iz)-ey(ix, iy, iz)) &
              + cny*(ex(ix, iy+1, iz)-ex(ix, iy, iz))
        ENDDO
      ENDDO
    ENDDO

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.FALSE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    REAL(num) :: lx, cnx, ly, cny, lz, cnz

    lx = dt/dx
    cnx = 0.5_num*lx

    ly = dt/dy
    cny = 0.5_num*ly

    lz = dt/dz
    cnz = 0.5_num*lz

    ! bx
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          bx(ix, iy, iz) = bx(ix, iy, iz) &
              + cnz*(ey(ix, iy, iz+1)-ey(ix, iy, iz)) &
              - cny*(ez(ix, iy+1, iz)-ez(ix, iy, iz))
        ENDDO
      ENDDO
    ENDDO

    ! by
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          by(ix, iy, iz) = by(ix, iy, iz) &
              + cnx*(ez(ix+1, iy, iz)-ez(ix, iy, iz)) &
              - cnz*(ex(ix, iy, iz+1)-ex(ix, iy, iz))
        ENDDO
      ENDDO
    ENDDO

    ! bz
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          bz(ix, iy, iz) = bz(ix, iy, iz) &
              - cnx*(ey(ix+1, iy, iz)-ey(ix, iy, iz)) &
              + cny*(ex(ix, iy+1, iz)-ex(ix, iy, iz))
        ENDDO
      ENDDO
    ENDDO

    CALL bfield_bcs(.FALSE.)
    IF (xbc_left .EQ. c_bc_simple_laser .AND. left .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_left
    IF (xbc_left .EQ. c_bc_simple_outflow .AND. left .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_left

    IF (xbc_right .EQ. c_bc_simple_laser .AND. right .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_right
    IF (xbc_right .EQ. c_bc_simple_outflow .AND. right .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_right

    IF (ybc_up .EQ. c_bc_simple_laser .AND. up .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_up
    IF (ybc_up .EQ. c_bc_simple_outflow .AND. up .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_up

    IF (ybc_down .EQ. c_bc_simple_laser .AND. down .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_down
    IF (ybc_down .EQ. c_bc_simple_outflow .AND. down .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_down

    IF (zbc_front .EQ. c_bc_simple_laser .AND. up .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_front
    IF (zbc_back .EQ. c_bc_simple_outflow .AND. up .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_front

    IF (zbc_front .EQ. c_bc_simple_laser .AND. down .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_back
    IF (zbc_back .EQ. c_bc_simple_outflow .AND. down .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_back
    CALL bfield_bcs(.TRUE.)

    ! ex
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ex(ix, iy, iz) = ex(ix, iy, iz) &
              + cny*(bz(ix, iy, iz)-bz(ix, iy-1, iz)) * c**2 &
              - cnz*(by(ix, iy, iz)-by(ix, iy, iz-1)) * c**2 &
              - 0.5*dt*jx(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! ey
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ey(ix, iy, iz) = ey(ix, iy, iz) &
              + cnz*(bx(ix, iy, iz)-bx(ix, iy, iz-1)) * c**2 &
              - cnx*(bz(ix, iy, iz)-bz(ix-1, iy, iz)) * c**2 &
              - 0.5*dt*jy(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! ez
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ez(ix, iy, iz) = ez(ix, iy, iz) &
              + cnx*(by(ix, iy, iz)-by(ix-1, iy, iz))*c**2 &
              - cny*(bx(ix, iy, iz)-bx(ix, iy-1, iz))*c**2 &
              - 0.5*dt*jz(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
