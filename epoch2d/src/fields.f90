MODULE fields

  USE boundary
  USE laser

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

    IF (field_order .EQ. 2) THEN
      const(1:2) = (/ -1.0_num, 1.0_num /)
      large = 1
      small = 0
    ELSE IF (field_order .EQ. 4) THEN
      const(1:4) = (/ 1.0_num/24.0_num, -9.0_num/8.0_num, &
          9.0_num/8.0_num, -1.0_num/24.0_num /)
      large = 2
      small = 1
    ELSE
      const(1:6) = (/ -3.0_num/640.0_num, 25.0_num/384.0_num, &
          -75.0_num/64.0_num, 75.0_num/64.0_num, -25.0_num/384.0_num, &
          3.0_num/640.0_num /)
      large = 3
      small = 2
    ENDIF

  END SUBROUTINE set_field_order



  SUBROUTINE update_e_field

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

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

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
    CALL bfield_bcs(.FALSE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    CALL update_b_field

    CALL bfield_bcs(.FALSE.)

    IF (proc_x_min .EQ. MPI_PROC_NULL) THEN
      IF (bc_x_min_field .EQ. c_bc_simple_laser) THEN
        CALL laser_bcs_x_min
      ELSE IF (bc_x_min_field .EQ. c_bc_simple_outflow) THEN
        CALL outflow_bcs_x_min
      ENDIF
    ENDIF

    IF (proc_x_max .EQ. MPI_PROC_NULL) THEN
      IF (bc_x_max_field .EQ. c_bc_simple_laser) THEN
        CALL laser_bcs_x_max
      ELSE IF (bc_x_max_field .EQ. c_bc_simple_outflow) THEN
        CALL outflow_bcs_x_max
      ENDIF
    ENDIF

    IF (proc_y_min .EQ. MPI_PROC_NULL) THEN
      IF (bc_y_min_field .EQ. c_bc_simple_laser) THEN
        CALL laser_bcs_y_min
      ELSE IF (bc_y_min_field .EQ. c_bc_simple_outflow) THEN
        CALL outflow_bcs_y_min
      ENDIF
    ENDIF

    IF (proc_y_max .EQ. MPI_PROC_NULL) THEN
      IF (bc_y_max_field .EQ. c_bc_simple_laser) THEN
        CALL laser_bcs_y_max
      ELSE IF (bc_y_max_field .EQ. c_bc_simple_outflow) THEN
        CALL outflow_bcs_y_max
      ENDIF
    ENDIF

    CALL bfield_bcs(.TRUE.)

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
