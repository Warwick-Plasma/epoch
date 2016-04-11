! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE fields

  USE boundary

  IMPLICIT NONE

  INTEGER :: field_order
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx
  REAL(num) :: cnx

CONTAINS

  SUBROUTINE set_field_order(order)

    INTEGER, INTENT(IN) :: order

    field_order = order
    fng = field_order / 2

    IF (field_order == 2) THEN
      cfl = 1.0_num
    ELSE IF (field_order == 4) THEN
      cfl = 6.0_num / 7.0_num
    ELSE
      cfl = 120.0_num / 149.0_num
    ENDIF

  END SUBROUTINE set_field_order



  SUBROUTINE update_e_field

    INTEGER :: ix
    REAL(num) :: cpml_x
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3

    IF (cpml_boundaries) THEN
      cpml_x = cnx

      IF (field_order == 2) THEN
        DO ix = 1, nx
          cpml_x = cnx / cpml_kappa_ex(ix)

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cpml_x * (bz(ix  ) - bz(ix-1)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cpml_x * (by(ix  ) - by(ix-1)) &
              - fac * jz(ix)
        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO ix = 1, nx
          cpml_x = cnx / cpml_kappa_ex(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              - fac * jz(ix)
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO ix = 1, nx
          cpml_x = cnx / cpml_kappa_ex(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x
          cx3 = c3 * cpml_x

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - cx3 * (bz(ix+2) - bz(ix-3)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              + cx3 * (by(ix+2) - by(ix-3)) &
              - fac * jz(ix)
        ENDDO
      ENDIF

      CALL cpml_advance_e_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        DO ix = 1, nx
          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cnx * (bz(ix  ) - bz(ix-1)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cnx * (by(ix  ) - by(ix-1)) &
              - fac * jz(ix)
        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO ix = 1, nx
          cx1 = c1 * cnx
          cx2 = c2 * cnx

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              - fac * jz(ix)
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO ix = 1, nx
          cx1 = c1 * cnx
          cx2 = c2 * cnx
          cx3 = c3 * cnx

          ex(ix) = ex(ix) &
              - fac * jx(ix)

          ey(ix) = ey(ix) &
              - cx1 * (bz(ix  ) - bz(ix-1)) &
              - cx2 * (bz(ix+1) - bz(ix-2)) &
              - cx3 * (bz(ix+2) - bz(ix-3)) &
              - fac * jy(ix)

          ez(ix) = ez(ix) &
              + cx1 * (by(ix  ) - by(ix-1)) &
              + cx2 * (by(ix+1) - by(ix-2)) &
              + cx3 * (by(ix+2) - by(ix-3)) &
              - fac * jz(ix)
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    INTEGER :: ix
    REAL(num) :: cpml_x
    REAL(num) :: c1, c2, c3
    REAL(num) :: cx1, cx2, cx3

    REAL(num) :: alphax, deltax

    INTEGER :: maxwell_solver

    maxwell_solver = 0 ! 0=Yee 1=Lehe

    IF (maxwell_solver == 1) THEN    

    deltax = (1.0_num / 4.0_num)*(1 - ((dx/(c*dt))**2)*(sin(pi*c*dt/(2*dx)))**2)
    alphax = (1 - 3*deltax)

    ENDIF

    IF (cpml_boundaries) THEN
      cpml_x = hdtx

      IF (field_order == 2) THEN
        DO ix = 1, nx
          cpml_x = hdtx / cpml_kappa_bx(ix)

          IF (maxwell_solver == 0) THEN
            by(ix) = by(ix) &
                + cpml_x * (ez(ix+1) - ez(ix  ))

            bz(ix) = bz(ix) &
                - cpml_x * (ey(ix+1) - ey(ix  ))
          ELSE  IF (maxwell_solver == 1) THEN
            by(ix) = by(ix) &
               + cpml_x*(alphax*(ez(ix+1) - ez(ix  )) &
                        +deltax*(ez(ix+2) - ez(ix-1)) &
                        )

            bz(ix) = bz(ix) &
               - cpml_x*(alphax*(ey(ix+1) - ey(ix  )) &
                        +deltax*(ey(ix+2) - ey(ix-1)) &
                        )
          ENDIF

        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO ix = 1, nx
          cpml_x = hdtx / cpml_kappa_bx(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x

          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1))
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO ix = 1, nx
          cpml_x = hdtx / cpml_kappa_bx(ix)
          cx1 = c1 * cpml_x
          cx2 = c2 * cpml_x
          cx3 = c3 * cpml_x

          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1)) &
              + cx3 * (ez(ix+3) - ez(ix-2))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1)) &
              - cx3 * (ey(ix+3) - ey(ix-2))
        ENDDO
      ENDIF

      CALL cpml_advance_b_currents(hdt)
    ELSE
      IF (field_order == 2) THEN
        DO ix = 1, nx
          IF (maxwell_solver == 0) THEN
            by(ix) = by(ix) &
                + hdtx * (ez(ix+1) - ez(ix  ))

            bz(ix) = bz(ix) &
                - hdtx * (ey(ix+1) - ey(ix  ))
          ELSE  IF (maxwell_solver == 1) THEN
            by(ix) = by(ix) &
               + hdtx*(alphax*(ez(ix+1) - ez(ix  )) &
                        +deltax*(ez(ix+2) - ez(ix-1)) &
                        )

            bz(ix) = bz(ix) &
               - hdtx*(alphax*(ey(ix+1) - ey(ix  )) &
                        +deltax*(ey(ix+2) - ey(ix-1)) &
                        )
          ENDIF

        ENDDO
      ELSE IF (field_order == 4) THEN
        c1 = 9.0_num / 8.0_num
        c2 = -1.0_num / 24.0_num

        DO ix = 1, nx
          cx1 = c1 * hdtx
          cx2 = c2 * hdtx

          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1))
        ENDDO
      ELSE
        c1 = 75.0_num / 64.0_num
        c2 = -25.0_num / 384.0_num
        c3 = 3.0_num / 640.0_num

        DO ix = 1, nx
          cx1 = c1 * hdtx
          cx2 = c2 * hdtx
          cx3 = c3 * hdtx

          by(ix) = by(ix) &
              + cx1 * (ez(ix+1) - ez(ix  )) &
              + cx2 * (ez(ix+2) - ez(ix-1)) &
              + cx3 * (ez(ix+3) - ez(ix-2))

          bz(ix) = bz(ix) &
              - cx1 * (ey(ix+1) - ey(ix  )) &
              - cx2 * (ey(ix+2) - ey(ix-1)) &
              - cx3 * (ey(ix+3) - ey(ix-2))
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

    hdt  = 0.5_num * dt
    hdtx = hdt / dx

    cnx = hdtx * c**2

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

    cnx = hdtx * c**2

    fac = hdt / epsilon0

    CALL update_b_field

    CALL bfield_final_bcs

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
