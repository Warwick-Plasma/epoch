MODULE fields

  USE boundary

  IMPLICIT NONE

  REAL(num), DIMENSION(6) :: const
  INTEGER :: large, small, order
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx, hdty, hdtz
  REAL(num) :: cnx, cny, cnz

CONTAINS

  SUBROUTINE set_field_order(field_order)

    INTEGER, INTENT(IN) :: field_order

    order = field_order
    large = order / 2
    small = large - 1
    ng = large

    IF (field_order .EQ. 2) THEN
      const(1:2) = (/ -1.0_num, 1.0_num /)
      cfl = 1.0_num
    ELSE IF (field_order .EQ. 4) THEN
      const(1:4) = (/ 1.0_num/24.0_num, -9.0_num/8.0_num, &
          9.0_num/8.0_num, -1.0_num/24.0_num /)
      cfl = 6.0_num / 7.0_num
    ELSE
      const(1:6) = (/ -3.0_num/640.0_num, 25.0_num/384.0_num, &
          -75.0_num/64.0_num, 75.0_num/64.0_num, -25.0_num/384.0_num, &
          3.0_num/640.0_num /)
      cfl = 120.0_num / 149.0_num
    ENDIF

  END SUBROUTINE set_field_order



  SUBROUTINE update_e_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z, j_extra = 0

    IF (cpml_boundaries) THEN
      cpml_x = cnx
      cpml_y = cny
      cpml_z = cnz

      DO iz = 1, nz
        cpml_z = cnz / cpml_kappa_e_dz(iz)
        DO iy = 1, ny
          cpml_y = cny / cpml_kappa_e_dy(iy)
          DO ix = 1, nx
            j_extra = (cpml_e_psixz(ix,iy,iz) - cpml_e_psixy(ix,iy,iz)) / mu0
            ex(ix, iy, iz) = ex(ix, iy, iz) &
                + cpml_y * SUM(const(1:order) * bz(ix, iy-large:iy+small, iz)) &
                - cpml_z * SUM(const(1:order) * by(ix, iy, iz-large:iz+small)) &
                - fac * (jx(ix, iy, iz) + j_extra)
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        cpml_z = cnz / cpml_kappa_e_dz(iz)
        DO iy = 1, ny
          DO ix = 1, nx
            cpml_x = cnx / cpml_kappa_e_dx(ix)
            j_extra = (cpml_e_psiyx(ix,iy,iz) - cpml_e_psiyz(ix,iy,iz)) / mu0
            ey(ix, iy, iz) = ey(ix, iy, iz) &
                + cpml_z * SUM(const(1:order) * bx(ix, iy, iz-large:iz+small)) &
                - cpml_x * SUM(const(1:order) * bz(ix-large:ix+small, iy, iz)) &
                - fac * (jy(ix, iy, iz) + j_extra)
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        DO iy = 1, ny
          cpml_y = cny / cpml_kappa_e_dy(iy)
          DO ix = 1, nx
            cpml_x = cnx / cpml_kappa_e_dx(ix)
            j_extra = (cpml_e_psizy(ix,iy,iz) - cpml_e_psizx(ix,iy,iz)) / mu0
            ez(ix, iy, iz) = ez(ix, iy, iz) &
                + cpml_x * SUM(const(1:order) * by(ix-large:ix+small, iy, iz)) &
                - cpml_y * SUM(const(1:order) * bx(ix, iy-large:iy+small, iz)) &
                - fac * (jz(ix, iy, iz) + j_extra)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ex(ix, iy, iz) = ex(ix, iy, iz) &
                + cny * SUM(const(1:order) * bz(ix, iy-large:iy+small, iz)) &
                - cnz * SUM(const(1:order) * by(ix, iy, iz-large:iz+small)) &
                - fac * jx(ix, iy, iz)
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ey(ix, iy, iz) = ey(ix, iy, iz) &
                + cnz * SUM(const(1:order) * bx(ix, iy, iz-large:iz+small)) &
                - cnx * SUM(const(1:order) * bz(ix-large:ix+small, iy, iz)) &
                - fac * jy(ix, iy, iz)
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ez(ix, iy, iz) = ez(ix, iy, iz) &
                + cnx * SUM(const(1:order) * by(ix-large:ix+small, iy, iz)) &
                - cny * SUM(const(1:order) * bx(ix, iy-large:iy+small, iz)) &
                - fac * jz(ix, iy, iz)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z, j_extra = 0

    IF (cpml_boundaries) THEN
      cpml_x = hdtx
      cpml_y = hdty
      cpml_z = hdtz

      DO iz = 1, nz
        cpml_z = hdtz / cpml_kappa_b_dz(iz)
        DO iy = 1, ny
          cpml_y = hdty / cpml_kappa_b_dy(iy)
          DO ix = 1, nx
            j_extra = cpml_b_psixy(ix,iy,iz) - cpml_b_psixz(ix,iy,iz)
            bx(ix, iy, iz) = bx(ix, iy, iz) &
                - cpml_y * SUM(const(1:order) * ez(ix, iy-small:iy+large, iz)) &
                + cpml_z * SUM(const(1:order) * ey(ix, iy, iz-small:iz+large)) &
                - hdt * j_extra
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        cpml_z = hdtz / cpml_kappa_b_dz(iz)
        DO iy = 1, ny
          DO ix = 1, nx
            cpml_x = hdtx / cpml_kappa_b_dx(ix)
            j_extra = cpml_b_psiyz(ix,iy,iz) - cpml_b_psiyx(ix,iy,iz)
            by(ix, iy, iz) = by(ix, iy, iz) &
                - cpml_z * SUM(const(1:order) * ex(ix, iy, iz-small:iz+large)) &
                + cpml_x * SUM(const(1:order) * ez(ix-small:ix+large, iy, iz)) &
                - hdt * j_extra
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        DO iy = 1, ny
          cpml_y = hdty / cpml_kappa_b_dy(iy)
          DO ix = 1, nx
            cpml_x = hdtx / cpml_kappa_b_dx(ix)
            j_extra = cpml_b_psizx(ix,iy,iz) - cpml_b_psizy(ix,iy,iz)
            bz(ix, iy, iz) = bz(ix, iy, iz) &
                - cpml_x * SUM(const(1:order) * ey(ix-small:ix+large, iy, iz)) &
                + cpml_y * SUM(const(1:order) * ex(ix, iy-small:iy+large, iz)) &
                - hdt * j_extra
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            bx(ix, iy, iz) = bx(ix, iy, iz) &
                - hdty * SUM(const(1:order) * ez(ix, iy-small:iy+large, iz)) &
                + hdtz * SUM(const(1:order) * ey(ix, iy, iz-small:iz+large))
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            by(ix, iy, iz) = by(ix, iy, iz) &
                - hdtz * SUM(const(1:order) * ex(ix, iy, iz-small:iz+large)) &
                + hdtx * SUM(const(1:order) * ez(ix-small:ix+large, iy, iz))
          ENDDO
        ENDDO
      ENDDO

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            bz(ix, iy, iz) = bz(ix, iy, iz) &
                - hdtx * SUM(const(1:order) * ey(ix-small:ix+large, iy, iz)) &
                + hdty * SUM(const(1:order) * ex(ix, iy-small:iy+large, iz))
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy
    hdtz = hdt / dz

    cnx = hdtx * c**2
    cny = hdty * c**2
    cnz = hdtz * c**2

    fac = hdt / epsilon0

    ! Update E field to t+dt/2
    CALL update_e_field

    ! Now have E(t+dt/2), do boundary conditions on E
    CALL efield_bcs

    IF (cpml_boundaries) CALL cpml_advance_b_currents(dt)

    ! Update B field to t+dt/2 using E(t+dt/2)
    CALL update_b_field

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.TRUE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    CALL update_b_field

    CALL bfield_final_bcs

    IF (cpml_boundaries) CALL cpml_advance_e_currents(dt)

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
