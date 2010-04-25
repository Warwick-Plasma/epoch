MODULE fields

  USE boundary
  USE laser

  IMPLICIT NONE

  REAL(num), DIMENSION(6) :: const
  INTEGER :: large, small, order
  REAL(num) :: hdt, fac
  REAL(num) :: hdtx
  REAL(num) :: cnx

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

    DO ix = 1, nx
      ex(ix) = ex(ix) &
          - fac * jx(ix)
    ENDDO

    DO ix = 1, nx
      ey(ix) = ey(ix) &
          - cnx * SUM(const(1:order) * bz(ix-large:ix+small)) &
          - fac * jy(ix)
    ENDDO

    DO ix = 1, nx
      ez(ix) = ez(ix) &
          + cnx * SUM(const(1:order) * by(ix-large:ix+small)) &
          - fac * jz(ix)
    ENDDO

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    DO ix = 1, nx
      by(ix) = by(ix) &
          + hdtx * SUM(const(1:order) * ez(ix-small:ix+large))
    ENDDO

    DO ix = 1, nx
      bz(ix) = bz(ix) &
          - hdtx * SUM(const(1:order) * ey(ix-small:ix+large))
    ENDDO

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
    CALL bfield_bcs(.FALSE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    CALL update_b_field

    CALL bfield_bcs(.FALSE.)

    IF (bc_x_min .EQ. c_bc_simple_laser .AND. proc_x_min .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_x_min
    IF (bc_x_min .EQ. c_bc_simple_outflow .AND. proc_x_min .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_x_min

    IF (bc_x_max .EQ. c_bc_simple_laser .AND. proc_x_max .EQ. MPI_PROC_NULL) &
        CALL laser_bcs_x_max
    IF (bc_x_max .EQ. c_bc_simple_outflow .AND. proc_x_max .EQ. MPI_PROC_NULL) &
        CALL outflow_bcs_x_max

    CALL bfield_bcs(.TRUE.)

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
