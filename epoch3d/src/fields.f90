MODULE fields

  USE boundary
  USE laser

  IMPLICIT NONE

  REAL(num), DIMENSION(6) :: diff_consts
  INTEGER :: large, small, field_order

CONTAINS

  SUBROUTINE set_field_order(order)

    INTEGER, INTENT(IN) :: order

    field_order = order

    IF (field_order .EQ. 2) THEN
      diff_consts(1:2) = (/ -1.0_num, 1.0_num /)
      large = 1
      small = 0
    ELSE IF (field_order .EQ. 4) THEN
      diff_consts(1:4) = (/ 1.0_num/24.0_num, -9.0_num/8.0_num, &
          9.0_num/8.0_num, -1.0_num/24.0_num /)
      large = 2
      small = 1
    ELSE
      diff_consts(1:6) = (/ 3.0_num/640.0_num, -25.0_num/384.0_num, &
          75.0_num/64.0_num, -75.0_num/64.0_num, 25.0_num/384.0_num, &
          -3.0_num/640.0_num /)
      large = 3
      small = 2
    ENDIF

  END SUBROUTINE set_field_order



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
              + cny*c**2 * SUM(diff_consts * bz(ix, iy-large:iy+small, iz)) &
              - cnz*c**2 * SUM(diff_consts * by(ix, iy, iz-large:iz+small)) &
              - 0.5_num*dt*jx(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! Update ey to t = t0+dt/2
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ey(ix, iy, iz) = ey(ix, iy, iz) &
              + cnz*c**2 * SUM(diff_consts * bx(ix, iy, iz-large:iz+small)) &
              - cnx*c**2 * SUM(diff_consts * bz(ix-large:ix+small, iy, iz)) &
              - 0.5_num*dt*jy(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! Update ez to t = t0+dt/2
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ez(ix, iy, iz) = ez(ix, iy, iz) &
              + cnx*c**2 * SUM(diff_consts * by(ix-large:ix+small, iy, iz)) &
              - cny*c**2 * SUM(diff_consts * bx(ix, iy-large:iy+small, iz)) &
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
              + cnz * SUM(diff_consts * ey(ix, iy, iz-small:iz+large)) &
              - cny * SUM(diff_consts * ez(ix, iy-small:iy+large), iz)
        ENDDO
      ENDDO
    ENDDO

    ! by
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          by(ix, iy, iz) = by(ix, iy, iz) &
              + cnx * SUM(diff_consts * ez(ix-small:ix+large), iy, iz) &
              - cnz * SUM(diff_consts * ex(ix, iy, iz-small:iz+large))
        ENDDO
      ENDDO
    ENDDO

    ! bz
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          bz(ix, iy, iz) = bz(ix, iy, iz) &
              - cnx * SUM(diff_consts * ey(ix-small:ix+large), iy, iz) &
              + cny * SUM(diff_consts * ex(ix, iy-small:iy+large), iz)
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
              + cnz * SUM(diff_consts * ey(ix, iy, iz-small:iz+large)) &
              - cny * SUM(diff_consts * ez(ix, iy-small:iy+large, iz))
        ENDDO
      ENDDO
    ENDDO

    ! by
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          by(ix, iy, iz) = by(ix, iy, iz) &
              + cnx * SUM(diff_consts * ez(ix-small:ix+large, iy, iz)) &
              - cnz * SUM(diff_consts * ex(ix, iy, iz-small:iz+large))
        ENDDO
      ENDDO
    ENDDO

    ! bz
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          bz(ix, iy, iz) = bz(ix, iy, iz) &
              - cnx * SUM(diff_consts * ey(ix-small:ix+large, iy, iz)) &
              + cny * SUM(diff_consts * ex(ix, iy-small:iy+large, iz))
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
              + cny*c**2 * SUM(diff_consts * bz(ix, iy-large:iy+small, iz)) &
              - cnz*c**2 * SUM(diff_consts * by(ix, iy, iz-large:iz+small)) &
              - 0.5*dt*jx(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! ey
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ey(ix, iy, iz) = ey(ix, iy, iz) &
              + cnz*c**2 * SUM(diff_consts * bx(ix, iy, iz-large:iz+small)) &
              - cnx*c**2 * SUM(diff_consts * bz(ix-large:ix+small, iy, iz)) &
              - 0.5*dt*jy(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    ! ez
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          ez(ix, iy, iz) = ez(ix, iy, iz) &
              + cnx*c**2 * SUM(diff_consts * by(ix-large:ix+small, iy, iz)) &
              - cny*c**2 * SUM(diff_consts * bx(ix, iy-large:iy+small, iz)) &
              - 0.5*dt*jz(ix, iy, iz)/epsilon0
        ENDDO
      ENDDO
    ENDDO

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
