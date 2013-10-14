MODULE fields

  USE boundary

  IMPLICIT NONE

  INTEGER :: field_order
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx, hdty, hdtz
  REAL(num) :: cnx, cny, cnz

CONTAINS

  SUBROUTINE set_field_order(order)

    INTEGER, INTENT(IN) :: order

    field_order = order
    fng = field_order / 2

    IF (field_order .EQ. 2) THEN
      cfl = 1.0_num
    ELSE IF (field_order .EQ. 4) THEN
      cfl = 6.0_num / 7.0_num
    ELSE
      cfl = 120.0_num / 149.0_num
    ENDIF

  END SUBROUTINE set_field_order



  SUBROUTINE update_e_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3
    REAL(num) :: cz1, cz2, cz3

    IF (cpml_boundaries) THEN
      cpml_x = cnx
      cpml_y = cny
      cpml_z = cnz

      IF (field_order .EQ. 2) THEN
        DO iz = 1, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          DO iy = 1, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            DO ix = 1, nx
              cpml_x = cnx / cpml_kappa_ex(ix)

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cpml_y * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  - cpml_z * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cpml_z * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  - cpml_x * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cpml_x * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  - cpml_y * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order .EQ. 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          DO iy = 1, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            DO ix = 1, nx
              cpml_x = cnx / cpml_kappa_ex(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cpml_z = cnz / cpml_kappa_ez(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          cz3 = c3 * cpml_z
          DO iy = 1, ny
            cpml_y = cny / cpml_kappa_ey(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            cy3 = c3 * cpml_y
            DO ix = 1, nx
              cpml_x = cnx / cpml_kappa_ex(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x
              cx3 = c3 * cpml_x

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  + cy3 * (bz(ix  , iy+2, iz  ) - bz(ix  , iy-3, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - cz3 * (by(ix  , iy  , iz+2) - by(ix  , iy  , iz-3)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  + cz3 * (bx(ix  , iy  , iz+2) - bx(ix  , iy  , iz-3)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - cx3 * (bz(ix+2, iy  , iz  ) - bz(ix-3, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  + cx3 * (by(ix+2, iy  , iz  ) - by(ix-3, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - cy3 * (bx(ix  , iy+2, iz  ) - bx(ix  , iy-3, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL cpml_advance_e_currents(hdt)
    ELSE
      IF (field_order .EQ. 2) THEN
        DO iz = 1, nz
          DO iy = 1, ny
            DO ix = 1, nx

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cny * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  - cnz * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cnz * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  - cnx * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cnx * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  - cny * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order .EQ. 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cz1 = c1 * cnz
          cz2 = c2 * cnz
          DO iy = 1, ny
            cy1 = c1 * cny
            cy2 = c2 * cny
            DO ix = 1, nx
              cx1 = c1 * cnx
              cx2 = c2 * cnx

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cz1 = c1 * cnz
          cz2 = c2 * cnz
          cz3 = c3 * cnz
          DO iy = 1, ny
            cy1 = c1 * cny
            cy2 = c2 * cny
            cy3 = c3 * cny
            DO ix = 1, nx
              cx1 = c1 * cnx
              cx2 = c2 * cnx
              cx3 = c3 * cnx

              ex(ix, iy, iz) = ex(ix, iy, iz) &
                  + cy1 * (bz(ix  , iy  , iz  ) - bz(ix  , iy-1, iz  )) &
                  + cy2 * (bz(ix  , iy+1, iz  ) - bz(ix  , iy-2, iz  )) &
                  + cy3 * (bz(ix  , iy+2, iz  ) - bz(ix  , iy-3, iz  )) &
                  - cz1 * (by(ix  , iy  , iz  ) - by(ix  , iy  , iz-1)) &
                  - cz2 * (by(ix  , iy  , iz+1) - by(ix  , iy  , iz-2)) &
                  - cz3 * (by(ix  , iy  , iz+2) - by(ix  , iy  , iz-3)) &
                  - fac * jx(ix, iy, iz)

              ey(ix, iy, iz) = ey(ix, iy, iz) &
                  + cz1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy  , iz-1)) &
                  + cz2 * (bx(ix  , iy  , iz+1) - bx(ix  , iy  , iz-2)) &
                  + cz3 * (bx(ix  , iy  , iz+2) - bx(ix  , iy  , iz-3)) &
                  - cx1 * (bz(ix  , iy  , iz  ) - bz(ix-1, iy  , iz  )) &
                  - cx2 * (bz(ix+1, iy  , iz  ) - bz(ix-2, iy  , iz  )) &
                  - cx3 * (bz(ix+2, iy  , iz  ) - bz(ix-3, iy  , iz  )) &
                  - fac * jy(ix, iy, iz)

              ez(ix, iy, iz) = ez(ix, iy, iz) &
                  + cx1 * (by(ix  , iy  , iz  ) - by(ix-1, iy  , iz  )) &
                  + cx2 * (by(ix+1, iy  , iz  ) - by(ix-2, iy  , iz  )) &
                  + cx3 * (by(ix+2, iy  , iz  ) - by(ix-3, iy  , iz  )) &
                  - cy1 * (bx(ix  , iy  , iz  ) - bx(ix  , iy-1, iz  )) &
                  - cy2 * (bx(ix  , iy+1, iz  ) - bx(ix  , iy-2, iz  )) &
                  - cy3 * (bx(ix  , iy+2, iz  ) - bx(ix  , iy-3, iz  )) &
                  - fac * jz(ix, iy, iz)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix, iy, iz
    REAL(num) :: cpml_x, cpml_y, cpml_z
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3
    REAL(num) :: cy1, cy2, cy3
    REAL(num) :: cz1, cz2, cz3

    IF (cpml_boundaries) THEN
      cpml_x = hdtx
      cpml_y = hdty
      cpml_z = hdtz

      IF (field_order .EQ. 2) THEN
        DO iz = 1, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          DO iy = 1, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            DO ix = 1, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cpml_y * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  + cpml_z * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cpml_z * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  + cpml_x * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cpml_x * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  + cpml_y * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order .EQ. 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          DO iy = 1, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            DO ix = 1, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cpml_z = hdtz / cpml_kappa_bz(iz)
          cz1 = c1 * cpml_z
          cz2 = c2 * cpml_z
          cz3 = c3 * cpml_z
          DO iy = 1, ny
            cpml_y = hdty / cpml_kappa_by(iy)
            cy1 = c1 * cpml_y
            cy2 = c2 * cpml_y
            cy3 = c3 * cpml_y
            DO ix = 1, nx
              cpml_x = hdtx / cpml_kappa_bx(ix)
              cx1 = c1 * cpml_x
              cx2 = c2 * cpml_x
              cx3 = c3 * cpml_x

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  - cy3 * (ez(ix  , iy+3, iz  ) - ez(ix  , iy-2, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)) &
                  + cz3 * (ey(ix  , iy  , iz+3) - ey(ix  , iy  , iz-2))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  - cz3 * (ex(ix  , iy  , iz+3) - ex(ix  , iy  , iz-2)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )) &
                  + cx3 * (ez(ix+3, iy  , iz  ) - ez(ix-2, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  - cx3 * (ey(ix+3, iy  , iz  ) - ey(ix-2, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )) &
                  + cy3 * (ex(ix  , iy+3, iz  ) - ex(ix  , iy-2, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      CALL cpml_advance_b_currents(hdt)
    ELSE
      IF (field_order .EQ. 2) THEN
        DO iz = 1, nz
          DO iy = 1, ny
            DO ix = 1, nx
              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - hdty * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  + hdtz * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  ))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - hdtz * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  + hdtx * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - hdtx * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  + hdty * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (field_order .EQ. 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO iz = 1, nz
          cz1 = c1 * hdtz
          cz2 = c2 * hdtz
          DO iy = 1, ny
            cy1 = c1 * hdty
            cy2 = c2 * hdty
            DO ix = 1, nx
              cx1 = c1 * hdtx
              cx2 = c2 * hdtx

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO iz = 1, nz
          cz1 = c1 * hdtz
          cz2 = c2 * hdtz
          cz3 = c3 * hdtz
          DO iy = 1, ny
            cy1 = c1 * hdty
            cy2 = c2 * hdty
            cy3 = c3 * hdty
            DO ix = 1, nx
              cx1 = c1 * hdtx
              cx2 = c2 * hdtx
              cx3 = c3 * hdtx

              bx(ix, iy, iz) = bx(ix, iy, iz) &
                  - cy1 * (ez(ix  , iy+1, iz  ) - ez(ix  , iy  , iz  )) &
                  - cy2 * (ez(ix  , iy+2, iz  ) - ez(ix  , iy-1, iz  )) &
                  - cy3 * (ez(ix  , iy+3, iz  ) - ez(ix  , iy-2, iz  )) &
                  + cz1 * (ey(ix  , iy  , iz+1) - ey(ix  , iy  , iz  )) &
                  + cz2 * (ey(ix  , iy  , iz+2) - ey(ix  , iy  , iz-1)) &
                  + cz3 * (ey(ix  , iy  , iz+3) - ey(ix  , iy  , iz-2))

              by(ix, iy, iz) = by(ix, iy, iz) &
                  - cz1 * (ex(ix  , iy  , iz+1) - ex(ix  , iy  , iz  )) &
                  - cz2 * (ex(ix  , iy  , iz+2) - ex(ix  , iy  , iz-1)) &
                  - cz3 * (ex(ix  , iy  , iz+3) - ex(ix  , iy  , iz-2)) &
                  + cx1 * (ez(ix+1, iy  , iz  ) - ez(ix  , iy  , iz  )) &
                  + cx2 * (ez(ix+2, iy  , iz  ) - ez(ix-1, iy  , iz  )) &
                  + cx3 * (ez(ix+3, iy  , iz  ) - ez(ix-2, iy  , iz  ))

              bz(ix, iy, iz) = bz(ix, iy, iz) &
                  - cx1 * (ey(ix+1, iy  , iz  ) - ey(ix  , iy  , iz  )) &
                  - cx2 * (ey(ix+2, iy  , iz  ) - ey(ix-1, iy  , iz  )) &
                  - cx3 * (ey(ix+3, iy  , iz  ) - ey(ix-2, iy  , iz  )) &
                  + cy1 * (ex(ix  , iy+1, iz  ) - ex(ix  , iy  , iz  )) &
                  + cy2 * (ex(ix  , iy+2, iz  ) - ex(ix  , iy-1, iz  )) &
                  + cy3 * (ex(ix  , iy+3, iz  ) - ex(ix  , iy-2, iz  ))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
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

    ! Update B field to t+dt/2 using E(t+dt/2)
    CALL update_b_field

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.TRUE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy
    hdtz = hdt / dz

    cnx = hdtx * c**2
    cny = hdty * c**2
    cnz = hdtz * c**2

    fac = hdt / epsilon0

    CALL update_b_field

    CALL bfield_final_bcs

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
