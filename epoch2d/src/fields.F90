MODULE field

  USE boundary
  USE laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_eb_fields_half

    REAL(num) :: lx, cnx, ly, cny
    ! Note minus sign
#ifndef ORDER_SIX
    REAL(num), DIMENSION(4) :: diff_consts = (/1.0_num/24.0_num, &
        -9.0_num/8.0_num, 9.0_num/8.0_num, -1.0_num/24.0_num/)
    INTEGER, PARAMETER :: large = 2, small = 1
#else
    REAL(num), DIMENSION(6) :: diff_consts = (/3.0_num/640.0_num, &
        -25.0_num/384.0_num, 75.0_num/64.0_num, -75.0_num/64.0_num, &
        25.0_num/384.0_num, -3.0_num/640.0_num/)
    INTEGER, PARAMETER :: large = 3, small = 2
#endif

    lx = dt/dx
    cnx = 0.5_num*lx

    ly = dt/dy
    cny = 0.5_num*ly

    ! Update ex to t = t0+dt/2
#ifndef HIGH_ORDER_FIELDS
    DO iy = 1, ny
      DO ix = 1, nx
        ex(ix, iy) = ex(ix, iy) + &
            cny*(bz(ix, iy)-bz(ix, iy-1))*c**2 -0.5_num*dt*jx(ix, iy)/epsilon0
      ENDDO
    ENDDO

    ! Update ey to t = t0+dt/2
    DO iy = 1, ny
      DO ix = 1, nx
        ey(ix, iy) = ey(ix, iy) - &
            cnx*(bz(ix, iy)-bz(ix-1, iy))*c**2 -0.5_num*dt*jy(ix, iy)/epsilon0
      ENDDO
    ENDDO

    ! Update ez to t = t0+dt/2
    DO iy = 1, ny
      DO ix = 1, nx
        ez(ix, iy) = ez(ix, iy) + &
            cnx*(by(ix, iy)-by(ix-1, iy))*c**2 - &
            cny*(bx(ix, iy)-bx(ix, iy-1))*c**2 -0.5_num*dt*jz(ix, iy)/epsilon0
      ENDDO
    ENDDO
#else
    DO iy = 1, ny
      DO ix = 1, nx
        ex(ix, iy) = ex(ix, iy) + &
            cny*c**2 * SUM(diff_consts * bz(ix, iy-large:iy+small)) - &
            0.5_num*dt*jx(ix, iy)/epsilon0
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        ey(ix, iy) = ey(ix, iy) - &
            cnx*c**2 * SUM(diff_consts * bz(ix-large:ix+small, iy)) - &
            0.5_num*dt*jy(ix, iy)/epsilon0
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        ez(ix, iy) = ez(ix, iy) + &
            cnx*c**2 * SUM(diff_consts * by(ix-large:ix+small, iy)) - &
            cny*c**2 * SUM(diff_consts * bx(ix, iy-large:iy+small)) - &
            0.5_num*dt*jz(ix, iy)/epsilon0
      ENDDO
    ENDDO
#endif

    ! Now have E(t+dt/2), do boundary conditions on E

    CALL efield_bcs
#ifndef ELECTROSTATIC
    ! Update B field to t+dt/2 using E(t+dt/2)
#ifndef HIGH_ORDER_FIELDS
    ! bx
    DO iy = 1, ny
      DO ix = 1, nx
        bx(ix, iy) = bx(ix, iy) -cny*(ez(ix, iy+1)-ez(ix, iy))
      ENDDO
    ENDDO

    ! by
    DO iy = 1, ny
      DO ix = 1, nx
        by(ix, iy) = by(ix, iy) +cnx*(ez(ix+1, iy)-ez(ix, iy))
      ENDDO
    ENDDO

    ! bz
    DO iy = 1, ny
      DO ix = 1, nx
        bz(ix, iy) = bz(ix, iy) - &
            cnx*(ey(ix+1, iy)-ey(ix, iy)) +cny*(ex(ix, iy+1)-ex(ix, iy))
      ENDDO
    ENDDO
#else
    DO iy = 1, ny
      DO ix = 1, nx
        bx(ix, iy) = bx(ix, iy) - &
            cny * SUM(diff_consts * ez(ix, iy-small:iy+large))
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        by(ix, iy) = by(ix, iy) + &
            cnx * SUM(diff_consts * ez(ix-small:ix+large, iy))
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        bz(ix, iy) = bz(ix, iy) - &
            cnx * SUM(diff_consts * ey(ix-small:ix+large, iy)) + &
            cny * SUM(diff_consts * ex(ix, iy-small:iy+large))
      ENDDO
    ENDDO
#endif

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.FALSE.)
#endif
    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    REAL(num) :: lx, cnx, ly, cny

#ifndef ORDER_SIX
    REAL(num), DIMENSION(4) :: diff_consts = (/1.0_num/24.0_num, &
        -9.0_num/8.0_num, 9.0_num/8.0_num, -1.0_num/24.0_num/)
    INTEGER, PARAMETER :: large = 2, small = 1
#else
    REAL(num), DIMENSION(6) :: diff_consts = (/3.0_num/640.0_num, &
        -25.0_num/384.0_num, 75.0_num/64.0_num, -75.0_num/64.0_num, &
        25.0_num/384.0_num, -3.0_num/640.0_num/)
    INTEGER, PARAMETER :: large = 3, small = 2
#endif

    lx = dt/dx
    cnx = 0.5_num*lx

    ly = dt/dy
    cny = 0.5_num*ly
#ifndef ELECTROSTATIC
#ifndef HIGH_ORDER_FIELDS
    ! bx
    DO iy = 1, ny
      DO ix = 1, nx
        bx(ix, iy) = bx(ix, iy) -cny*(ez(ix, iy+1)-ez(ix, iy))
      ENDDO
    ENDDO

    ! by
    DO iy = 1, ny
      DO ix = 1, nx
        by(ix, iy) = by(ix, iy) +cnx*(ez(ix+1, iy)-ez(ix, iy))
      ENDDO
    ENDDO

    ! bz
    DO iy = 1, ny
      DO ix = 1, nx
        bz(ix, iy) = bz(ix, iy) - &
           cnx*(ey(ix+1, iy)-ey(ix, iy)) +cny*(ex(ix, iy+1)-ex(ix, iy))
      ENDDO
    ENDDO
#else
    DO iy = 1, ny
      DO ix = 1, nx
        bx(ix, iy) = bx(ix, iy) - &
           cny * SUM(diff_consts * ez(ix, iy-small:iy+large))
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        by(ix, iy) = by(ix, iy) + &
           cnx * SUM(diff_consts * ez(ix-small:ix+large, iy))
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        bz(ix, iy) = bz(ix, iy) - &
           cnx * SUM(diff_consts * ey(ix-small:ix+large, iy)) + &
           cny * SUM(diff_consts * ex(ix, iy-small:iy+large))
      ENDDO
    ENDDO
#endif

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
    CALL bfield_bcs(.TRUE.)
#endif
#ifndef HIGH_ORDER_FIELDS
    ! ex
    DO iy = 1, ny
      DO ix = 1, nx
        ex(ix, iy) = ex(ix, iy) + &
           cny*(bz(ix, iy)-bz(ix, iy-1)) * c**2 -0.5*dt*jx(ix, iy)/epsilon0
      ENDDO
    ENDDO

    ! ey
    DO iy = 1, ny
      DO ix = 1, nx
        ey(ix, iy) = ey(ix, iy) - &
           cnx*(bz(ix, iy)-bz(ix-1, iy)) * c **2 -0.5*dt*jy(ix, iy)/epsilon0
      ENDDO
    ENDDO

    ! ez
    DO iy = 1, ny
      DO ix = 1, nx
        ez(ix, iy) = ez(ix, iy) + &
           cnx*(by(ix, iy)-by(ix-1, iy))*c**2 - &
           cny*(bx(ix, iy)-bx(ix, iy-1))*c**2 -0.5*dt*jz(ix, iy)/epsilon0
      ENDDO
    ENDDO
#else
    DO iy = 1, ny
      DO ix = 1, nx
        ex(ix, iy) = ex(ix, iy) + &
           cny*c**2 * SUM(diff_consts * bz(ix, iy-large:iy+small)) - &
           0.5*dt*jx(ix, iy)/epsilon0
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        ey(ix, iy) = ey(ix, iy) - &
           cnx*c**2 * SUM(diff_consts * bz(ix-large:ix+small, iy)) - &
           0.5*dt*jy(ix, iy)/epsilon0
      ENDDO
    ENDDO

    DO iy = 1, ny
      DO ix = 1, nx
        ez(ix, iy) = ez(ix, iy) + &
           cnx*c**2 * SUM(diff_consts * by(ix-large:ix+small, iy)) - &
           cny*c**2 * SUM(diff_consts * bx(ix, iy-large:iy+small)) - &
           0.5_num*dt*jz(ix, iy)/epsilon0
      ENDDO
    ENDDO
#endif

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE field
