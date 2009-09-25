MODULE field
  USE shared_data
  USE boundary
  USE laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_eb_fields_half

    REAL(num) :: lx,cnx,ly,cny
    !Note minus sign
#ifndef ORDER_SIX
    REAL(num),DIMENSION(4) :: Diff_Consts&
         =(/1.0_num/24.0_num,-9.0_num/8.0_num,9.0_num/8.0_num,-1.0_num/24.0_num/)
    INTEGER,PARAMETER :: Large=2, Small=1
#else
    REAL(num),DIMENSION(6) :: Diff_Consts&
         =(/3.0_num/640.0_num, -25.0_num/384.0_num, 75.0_num/64.0_num,&
         -75.0_num/64.0_num,25.0_num/384.0_num,-3.0_num/640.0_num/)
    INTEGER,PARAMETER :: Large=3, Small=2
#endif


    lx=dt/dx
    cnx=0.5_num*lx

    ly=dt/dy
    cny=0.5_num*ly

    !Update Ex to t=t0+dt/2
#ifndef HIGH_ORDER_FIELDS
    DO iy=1,ny
       DO ix=1,nx
          Ex(ix,iy)=Ex(ix,iy)&
               +cny*(bz(ix,iy)-bz(ix,iy-1))*c**2&
               -0.5_num*dt*Jx(ix,iy)/epsilon0
       ENDDO
    ENDDO

    !Update Ey to t=t0+dt/2
    DO iy=1,ny
       DO ix=1,nx
          Ey(ix,iy)=Ey(ix,iy)&
               -cnx*(bz(ix,iy)-bz(ix-1,iy))*c**2&
               -0.5_num*dt*Jy(ix,iy)/epsilon0
       ENDDO
    ENDDO

    !Update Ez to t=t0+dt/2
    DO iy=1,ny
       DO ix=1,nx
          Ez(ix,iy)=Ez(ix,iy)&
               +cnx*(by(ix,iy)-by(ix-1,iy))*c**2&
               -cny*(bx(ix,iy)-bx(ix,iy-1))*c**2&
               -0.5_num*dt*Jz(ix,iy)/epsilon0
       ENDDO
    ENDDO
#else
    DO iy=1,ny
       DO ix=1,nx
          Ex(ix,iy)=Ex(ix,iy)&
               +cny*c**2 * SUM(Diff_Consts * Bz(ix,iy-Large:iy+Small))&
               -0.5_num*dt*Jx(ix,iy)/epsilon0
       ENDDO
    ENDDO
    DO iy=1,ny
       DO ix=1,nx
          Ey(ix,iy)=Ey(ix,iy)&
               -cnx*c**2 * SUM(Diff_Consts * Bz(ix-Large:ix+Small,iy))&
               -0.5_num*dt*Jy(ix,iy)/epsilon0
       ENDDO
    ENDDO

    DO iy=1,ny
       DO ix=1,nx
          Ez(ix,iy)=Ez(ix,iy)&
               +cnx*c**2 * SUM(Diff_Consts * By(ix-Large:ix+Small,iy))&
               -cny*c**2 * SUM(Diff_Consts * Bx(ix,iy-Large:iy+Small))&
               -0.5_num*dt*Jz(ix,iy)/epsilon0
       ENDDO
    ENDDO
#endif

    !Now have E(t+dt/2), do boundary conditions on E

    CALL Efield_bcs
#ifndef ELECTROSTATIC
    !Update B field to t+dt/2 using E(t+dt/2)
#ifndef HIGH_ORDER_FIELDS
    !Bx
    DO iy=1,ny
       DO ix=1,nx
          Bx(ix,iy) = Bx(ix,iy)&
               -cny*(Ez(ix,iy+1)-Ez(ix,iy))
       ENDDO
    ENDDO


    !By
    DO iy=1,ny
       DO ix=1,nx
          By(ix,iy)=By(ix,iy)&
               +cnx*(Ez(ix+1,iy)-Ez(ix,iy))
       ENDDO
    ENDDO

    !Bz
    DO iy=1,ny
       DO ix=1,nx
          Bz(ix,iy)=Bz(ix,iy)&
               -cnx*(Ey(ix+1,iy)-Ey(ix,iy))&
               +cny*(Ex(ix,iy+1)-Ex(ix,iy))
       ENDDO
    ENDDO
#else
    DO iy=1,ny
       DO ix=1,nx
          Bx(ix,iy)=Bx(ix,iy)&
               -cny * SUM(Diff_Consts * Ez(ix,iy-Small:iy+Large))
       ENDDO
    ENDDO
    DO iy=1,ny
       DO ix=1,nx
          By(ix,iy)=By(ix,iy)&
               +cnx * SUM(Diff_Consts * Ez(ix-Small:ix+Large,iy))
       ENDDO
    ENDDO
    DO iy=1,ny
       DO ix=1,nx
          Bz(ix,iy)=Bz(ix,iy)&
               -cnx * SUM(Diff_Consts * Ey(ix-Small:ix+Large,iy))&
               +cny * SUM(Diff_Consts * Ex(ix,iy-Small:iy+Large))
       ENDDO
    ENDDO
#endif

    !Now have B field at t+dt/2. Do boundary conditions on B
    CALL Bfield_bcs(.FALSE.)
#endif
    !Now have E&B fields at t=t+dt/2
    !Move to particle pusher

  END SUBROUTINE update_eb_fields_half

  SUBROUTINE update_eb_fields_final

    REAL(num) :: lx,cnx,ly,cny

#ifndef ORDER_SIX
    REAL(num),DIMENSION(4) :: Diff_Consts&
         =(/1.0_num/24.0_num,-9.0_num/8.0_num,9.0_num/8.0_num,-1.0_num/24.0_num/)
    INTEGER,PARAMETER :: Large=2, Small=1
#else
    REAL(num),DIMENSION(6) :: Diff_Consts&
         =(/3.0_num/640.0_num, -25.0_num/384.0_num, 75.0_num/64.0_num,&
         -75.0_num/64.0_num,25.0_num/384.0_num,-3.0_num/640.0_num/)
    INTEGER,PARAMETER :: Large=3, Small=2
#endif

    lx=dt/dx
    cnx=0.5_num*lx

    ly=dt/dy
    cny=0.5_num*ly
#ifndef ELECTROSTATIC
#ifndef HIGH_ORDER_FIELDS
    !Bx 
    DO iy=1,ny
       DO ix=1,nx
          Bx(ix,iy) = Bx(ix,iy)&
               -cny*(Ez(ix,iy+1)-Ez(ix,iy))
       ENDDO
    ENDDO

    !By
    DO iy=1,ny
       DO ix=1,nx
          By(ix,iy)=By(ix,iy)&
               +cnx*(Ez(ix+1,iy)-Ez(ix,iy))
       ENDDO
    ENDDO

    !Bz
    DO iy=1,ny
       DO ix=1,nx
          Bz(ix,iy)=Bz(ix,iy)&
               -cnx*(Ey(ix+1,iy)-Ey(ix,iy))&
               +cny*(Ex(ix,iy+1)-Ex(ix,iy))
       ENDDO
    ENDDO
#else
    DO iy=1,ny
       DO ix=1,nx
          Bx(ix,iy)=Bx(ix,iy)&
               -cny * SUM(Diff_Consts * Ez(ix,iy-Small:iy+Large))
       ENDDO
    ENDDO
    DO iy=1,ny
       DO ix=1,nx
          By(ix,iy)=By(ix,iy)&
               +cnx * SUM(Diff_Consts * Ez(ix-Small:ix+Large,iy))
       ENDDO
    ENDDO
    DO iy=1,ny
       DO ix=1,nx
          Bz(ix,iy)=Bz(ix,iy)&
               -cnx * SUM(Diff_Consts * Ey(ix-Small:ix+Large,iy))&
               +cny * SUM(Diff_Consts * Ex(ix,iy-Small:iy+Large))
       ENDDO
    ENDDO
#endif

    CALL BField_BCS(.FALSE.)
    IF(xbc_left == BC_SIMPLE_LASER .AND. left == MPI_PROC_NULL) CALL laser_bcs_left
    IF(xbc_left == BC_SIMPLE_OUTFLOW .AND. left == MPI_PROC_NULL) CALL outflow_bcs_left

    IF(xbc_right == BC_SIMPLE_LASER .AND. right == MPI_PROC_NULL) CALL laser_bcs_right
    IF(xbc_right == BC_SIMPLE_OUTFLOW .AND. right == MPI_PROC_NULL) CALL outflow_bcs_right

    IF(ybc_up == BC_SIMPLE_LASER .AND. up == MPI_PROC_NULL) CALL laser_bcs_up
    IF(ybc_up == BC_SIMPLE_OUTFLOW .AND. up == MPI_PROC_NULL) CALL outflow_bcs_up

    IF(ybc_down == BC_SIMPLE_LASER .AND. down == MPI_PROC_NULL) CALL laser_bcs_down
    IF(ybc_down == BC_SIMPLE_OUTFLOW .AND. down == MPI_PROC_NULL) CALL outflow_bcs_down
    CALL BField_BCS(.TRUE.)
#endif
#ifndef HIGH_ORDER_FIELDS
    !Ex
    DO iy=1,ny
       DO ix=1,nx
          Ex(ix,iy)=Ex(ix,iy)&
               +cny*(Bz(ix,iy)-Bz(ix,iy-1)) * c**2&
               -0.5*dt*Jx(ix,iy)/epsilon0
       ENDDO
    ENDDO

    !Ey
    DO iy=1,ny
       DO ix=1,nx
          Ey(ix,iy)=Ey(ix,iy)&
               -cnx*(Bz(ix,iy)-Bz(ix-1,iy)) * c **2&
               -0.5*dt*Jy(ix,iy)/epsilon0
       ENDDO
    ENDDO



    !Ez
    DO iy=1,ny
       DO ix=1,nx
          Ez(ix,iy)=Ez(ix,iy)&
               +cnx*(By(ix,iy)-By(ix-1,iy))*c**2&
               -cny*(Bx(ix,iy)-Bx(ix,iy-1))*c**2&
               -0.5*dt*Jz(ix,iy)/epsilon0
       ENDDO
    ENDDO
#else
    DO iy=1,ny
       DO ix=1,nx
          Ex(ix,iy)=Ex(ix,iy)&
               +cny*c**2 * SUM(Diff_Consts * Bz(ix,iy-Large:iy+Small))&
               -0.5*dt*Jx(ix,iy)/epsilon0
       ENDDO
    ENDDO
    DO iy=1,ny
       DO ix=1,nx
          Ey(ix,iy)=Ey(ix,iy)&
               -cnx*c**2 * SUM(Diff_Consts * Bz(ix-Large:ix+Small,iy))&
               -0.5*dt*Jy(ix,iy)/epsilon0
       ENDDO
    ENDDO
    DO iy=1,ny
       DO ix=1,nx
          Ez(ix,iy)=Ez(ix,iy)&
               +cnx*c**2 * SUM(Diff_Consts * By(ix-Large:ix+Small,iy))&
               -cny*c**2 * SUM(Diff_Consts * Bx(ix,iy-Large:iy+Small))&
               -0.5_num*dt*Jz(ix,iy)/epsilon0
       ENDDO
    ENDDO
#endif

    CALL Efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE field














