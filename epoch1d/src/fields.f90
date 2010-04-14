MODULE fields

  USE boundary
  USE laser

  IMPLICIT NONE

  REAL(num), DIMENSION(6) :: const
  INTEGER :: large, small, order

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
      const(1:6) = (/ 3.0_num/640.0_num, -25.0_num/384.0_num, &
          75.0_num/64.0_num, -75.0_num/64.0_num, 25.0_num/384.0_num, &
          -3.0_num/640.0_num /)
      large = 3
      small = 2
    ENDIF

  END SUBROUTINE set_field_order



  SUBROUTINE update_eb_fields_half

    REAL(num) :: lx, cnx

    lx = dt/dx
    cnx = 0.5_num*lx

    ! Update ex to t = t0+dt/2
    DO ix = 1, nx
      ex(ix) = ex(ix) &
          - 0.5_num*dt*jx(ix)/epsilon0
    ENDDO

    ! Update ey to t = t0+dt/2
    DO ix = 1, nx
      ey(ix) = ey(ix) &
          - cnx*c**2 * SUM(const(1:order) * bz(ix-large:ix+small)) &
          - 0.5_num*dt*jy(ix)/epsilon0
    ENDDO

    ! Update ez to t = t0+dt/2
    DO ix = 1, nx
      ez(ix) = ez(ix) &
          + cnx*c**2 * SUM(const(1:order) * by(ix-large:ix+small)) &
          - 0.5_num*dt*jz(ix)/epsilon0
    ENDDO

    ! Now have E(t+dt/2), do boundary conditions on E

    CALL efield_bcs

    ! Update B field to t+dt/2 using E(t+dt/2)

    ! bx

    ! by
    DO ix = 1, nx
      by(ix) = by(ix) &
          + cnx * SUM(const(1:order) * ez(ix-small:ix+large))
    ENDDO

    ! bz
    DO ix = 1, nx
      bz(ix) = bz(ix) &
          - cnx * SUM(const(1:order) * ey(ix-small:ix+large))
    ENDDO

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.FALSE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    REAL(num) :: lx, cnx

    lx = dt/dx
    cnx = 0.5_num*lx

    ! bx

    ! by
    DO ix = 1, nx
      by(ix) = by(ix) &
          + cnx * SUM(const(1:order) * ez(ix-small:ix+large))
    ENDDO

    ! bz
    DO ix = 1, nx
      bz(ix) = bz(ix) &
          - cnx * SUM(const(1:order) * ey(ix-small:ix+large))
    ENDDO

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

    ! ex
    DO ix = 1, nx
      ex(ix) = ex(ix) &
          - 0.5*dt*jx(ix)/epsilon0
    ENDDO

    ! ey
    DO ix = 1, nx
      ey(ix) = ey(ix) &
          - cnx*c**2 * SUM(const(1:order) * bz(ix-large:ix+small)) &
          - 0.5*dt*jy(ix)/epsilon0
    ENDDO

    ! ez
    DO ix = 1, nx
      ez(ix) = ez(ix) &
          + cnx*c**2 * SUM(const(1:order) * by(ix-large:ix+small)) &
          - 0.5*dt*jz(ix)/epsilon0
    ENDDO

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields
