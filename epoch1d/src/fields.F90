MODULE field

  USE shared_data
  USE boundary
  USE laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_eb_fields_half

#ifndef ORDER_SIX
    REAL(num), DIMENSION(4) :: diff_consts = (/1.0_num/24.0_num, -9.0_num/8.0_num, 9.0_num/8.0_num, -1.0_num/24.0_num/)
    INTEGER, PARAMETER :: large = 2, small = 1
#else
    REAL(num), DIMENSION(6) :: diff_consts = (/3.0_num/640.0_num, -25.0_num/384.0_num, 75.0_num/64.0_num, -75.0_num/64.0_num, 25.0_num/384.0_num, -3.0_num/640.0_num/)
    INTEGER, PARAMETER :: large = 3, small = 2
#endif

    REAL(num) :: lx, cnx

    lx = dt/dx
    cnx = 0.5_num*lx

    ! Update ex to t = t0+dt/2

#ifndef HIGH_ORDER_FIELDS

    DO ix = 1, nx
      ex(ix) = ex(ix) -0.5_num*dt*jx(ix)/epsilon0
    ENDDO

    ! Update ey to t = t0+dt/2
    DO ix = 1, nx
      ey(ix) = ey(ix) -cnx*(bz(ix)-bz(ix-1))*c**2 -0.5_num*dt*jy(ix)/epsilon0
    ENDDO

    ! Update ez to t = t0+dt/2
    DO ix = 1, nx
      ez(ix) = ez(ix) +cnx*(by(ix)-by(ix-1))*c**2 -0.5_num*dt*jz(ix)/epsilon0
    ENDDO

#else
    DO ix = 1, nx
      ex(ix) = ex(ix) -0.5_num*dt*jx(ix)/epsilon0
    ENDDO

    DO ix = 1, nx
      ey(ix) = ey(ix) -cnx*c**2 * SUM(diff_consts * bz(ix-large:ix+small)) -0.5_num*dt*jy(ix)/epsilon0
    ENDDO

    DO ix = 1, nx
      ez(ix) = ez(ix) +cnx*c**2 * SUM(diff_consts * by(ix-large:ix+small)) -0.5_num*dt*jz(ix)/epsilon0
    ENDDO
#endif

    ! Now have E(t+dt/2), do boundary conditions on E

    CALL efield_bcs

    ! Update B field to t+dt/2 using E(t+dt/2)

#ifndef HIGH_ORDER_FIELDS
    ! bx unchanged in 1D

    ! by
    DO ix = 1, nx
      by(ix) = by(ix) +cnx*(ez(ix+1)-ez(ix))
    ENDDO

    ! bz
    DO ix = 1, nx
      bz(ix) = bz(ix) -cnx*(ey(ix+1)-ey(ix))
    ENDDO
#else
    DO ix = 1, nx
      by(ix) = by(ix) +cnx * SUM(diff_consts * ez(ix-small:ix+large))
    ENDDO

    DO ix = 1, nx
      bz(ix) = bz(ix) -cnx * SUM(diff_consts * ey(ix-small:ix+large))
    ENDDO
#endif

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.FALSE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

#ifndef ORDER_SIX
    REAL(num), DIMENSION(4) :: diff_consts = (/1.0_num/24.0_num, -9.0_num/8.0_num, 9.0_num/8.0_num, -1.0_num/24.0_num/)
    INTEGER, PARAMETER :: large = 2, small = 1
#else
    REAL(num), DIMENSION(6) :: diff_consts = (/3.0_num/640.0_num, -25.0_num/384.0_num, 75.0_num/64.0_num, -75.0_num/64.0_num, 25.0_num/384.0_num, -3.0_num/640.0_num/)
    INTEGER, PARAMETER :: large = 3, small = 2
#endif

    REAL(num) :: lx, cnx

    lx = dt/dx
    cnx = 0.5_num*lx
#ifndef HIGH_ORDER_FIELDS
    ! bx unchanged in 1D

    ! by
    DO ix = 1, nx
      by(ix) = by(ix) +cnx*(ez(ix+1)-ez(ix))
    ENDDO

    ! bz
    DO ix = 1, nx
      bz(ix) = bz(ix) -cnx*(ey(ix+1)-ey(ix))
    ENDDO
#else
    DO ix = 1, nx
      by(ix) = by(ix) +cnx * SUM(diff_consts * ez(ix-small:ix+large))
    ENDDO

    DO ix = 1, nx
      bz(ix) = bz(ix) -cnx * SUM(diff_consts * ey(ix-small:ix+large))
    ENDDO
#endif

    CALL bfield_bcs(.FALSE.)
    IF(xbc_left == c_bc_simple_laser .AND. left == MPI_PROC_NULL) CALL laser_bcs_left
    IF(xbc_left == c_bc_simple_outflow .AND. left == MPI_PROC_NULL) CALL outflow_bcs_left

    IF(xbc_right == c_bc_simple_laser .AND. right == MPI_PROC_NULL) CALL laser_bcs_right
    IF(xbc_right == c_bc_simple_outflow .AND. right == MPI_PROC_NULL) CALL outflow_bcs_right
    CALL bfield_bcs(.TRUE.)

#ifndef HIGH_ORDER_FIELDS
    ! ex
    DO ix = 1, nx
      ex(ix) = ex(ix) -0.5*dt*jx(ix)/epsilon0
    ENDDO

    ! ey
    DO ix = 1, nx
      ey(ix) = ey(ix) -cnx*(bz(ix)-bz(ix-1)) * c**2 -0.5*dt*jy(ix)/epsilon0
    ENDDO

    ! ez
    DO ix = 1, nx
      ez(ix) = ez(ix) +cnx*(by(ix)-by(ix-1))*c**2 -0.5*dt*jz(ix)/epsilon0
    ENDDO
#else
    DO ix = 1, nx
      ex(ix) = ex(ix) -0.5_num*dt*jx(ix)/epsilon0
    ENDDO

    DO ix = 1, nx
      ey(ix) = ey(ix) -cnx*c**2 * SUM(diff_consts * bz(ix-large:ix+small)) -0.5_num*dt*jy(ix)/epsilon0
    ENDDO

    DO ix = 1, nx
      ez(ix) = ez(ix) +cnx*c**2 * SUM(diff_consts * by(ix-large:ix+small)) -0.5_num*dt*jz(ix)/epsilon0
    ENDDO
#endif

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE field
