MODULE fields

  USE boundary

  IMPLICIT NONE

  REAL(num), DIMENSION(6) :: const
  INTEGER :: large, small, order
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx, hdty
  REAL(num) :: cnx, cny

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

    INTEGER :: ix, iy
    REAL(num) :: cpml_x, cpml_y, j_extra = 0

    IF (cpml_boundaries) THEN
      CALL cpml_advance_e_currents(0.5_num * dt)

      cpml_x = cnx
      cpml_y = cny

      DO iy = 1, ny
        cpml_y = cny / cpml_kappa_ey(iy)
        DO ix = 1, nx
          j_extra = -cpml_psi_exy(ix,iy) / mu0
          ex(ix, iy) = ex(ix, iy) &
              + cpml_y * SUM(const(1:order) * bz(ix, iy-large:iy+small)) &
              - fac * (jx(ix, iy) + j_extra)
        ENDDO
      ENDDO

      DO iy = 1, ny
        DO ix = 1, nx
          cpml_x = cnx / cpml_kappa_ex(ix)
          j_extra = cpml_psi_eyx(ix,iy) / mu0
          ey(ix, iy) = ey(ix, iy) &
              - cpml_x * SUM(const(1:order) * bz(ix-large:ix+small, iy)) &
              - fac * (jy(ix, iy) + j_extra)
        ENDDO
      ENDDO

      DO iy = 1, ny
        cpml_y = cny / cpml_kappa_ey(iy)
        DO ix = 1, nx
          cpml_x = cnx / cpml_kappa_ex(ix)
          j_extra = (cpml_psi_ezy(ix,iy) - cpml_psi_ezx(ix,iy)) / mu0
          ez(ix, iy) = ez(ix, iy) &
              + cpml_x * SUM(const(1:order) * by(ix-large:ix+small, iy)) &
              - cpml_y * SUM(const(1:order) * bx(ix, iy-large:iy+small)) &
              - fac * (jz(ix, iy) + j_extra)
        ENDDO
      ENDDO
    ELSE
      DO iy = 1, ny
        DO ix = 1, nx
          ex(ix, iy) = ex(ix, iy) &
              + cny * SUM(const(1:order) * bz(ix, iy-large:iy+small)) &
              - fac * jx(ix, iy)
        ENDDO
      ENDDO

      DO iy = 1, ny
        DO ix = 1, nx
          ey(ix, iy) = ey(ix, iy) &
              - cnx * SUM(const(1:order) * bz(ix-large:ix+small, iy)) &
              - fac * jy(ix, iy)
        ENDDO
      ENDDO

      DO iy = 1, ny
        DO ix = 1, nx
          ez(ix, iy) = ez(ix, iy) &
              + cnx * SUM(const(1:order) * by(ix-large:ix+small, iy)) &
              - cny * SUM(const(1:order) * bx(ix, iy-large:iy+small)) &
              - fac * jz(ix, iy)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix, iy
    REAL(num) :: cpml_x, cpml_y, j_extra = 0

    IF (cpml_boundaries) THEN
      CALL cpml_advance_b_currents(0.5_num * dt)

      cpml_x = hdtx
      cpml_y = hdty

      DO iy = 1, ny
        cpml_y = hdty / cpml_kappa_by(iy)
        DO ix = 1, nx
          j_extra = cpml_psi_bxy(ix,iy)
          bx(ix, iy) = bx(ix, iy) &
              - cpml_y * SUM(const(1:order) * ez(ix, iy-small:iy+large)) &
              - hdt * j_extra
        ENDDO
      ENDDO

      DO iy = 1, ny
        DO ix = 1, nx
          cpml_x = hdtx / cpml_kappa_bx(ix)
          j_extra = -cpml_psi_byx(ix,iy)
          by(ix, iy) = by(ix, iy) &
              + cpml_x * SUM(const(1:order) * ez(ix-small:ix+large, iy)) &
              - hdt * j_extra
        ENDDO
      ENDDO

      DO iy = 1, ny
        cpml_y = hdty / cpml_kappa_by(iy)
        DO ix = 1, nx
          cpml_x = hdtx / cpml_kappa_bx(ix)
          j_extra = cpml_psi_bzx(ix,iy) - cpml_psi_bzy(ix,iy)
          bz(ix, iy) = bz(ix, iy) &
              - cpml_x * SUM(const(1:order) * ey(ix-small:ix+large, iy)) &
              + cpml_y * SUM(const(1:order) * ex(ix, iy-small:iy+large)) &
              - hdt * j_extra
        ENDDO
      ENDDO
    ELSE
      DO iy = 1, ny
        DO ix = 1, nx
          bx(ix, iy) = bx(ix, iy) &
              - hdty * SUM(const(1:order) * ez(ix, iy-small:iy+large))
        ENDDO
      ENDDO

      DO iy = 1, ny
        DO ix = 1, nx
          by(ix, iy) = by(ix, iy) &
              + hdtx * SUM(const(1:order) * ez(ix-small:ix+large, iy))
        ENDDO
      ENDDO

      DO iy = 1, ny
        DO ix = 1, nx
          bz(ix, iy) = bz(ix, iy) &
              - hdtx * SUM(const(1:order) * ey(ix-small:ix+large, iy)) &
              + hdty * SUM(const(1:order) * ex(ix, iy-small:iy+large))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy

    cnx = hdtx * c**2
    cny = hdty * c**2

    fac = hdt / epsilon0

    ! Update E field to t+dt/2
    CALL update_e_field

    ! Now have E(t+dt/2), do boundary conditions on E
    CALL efield_bcs

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

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
