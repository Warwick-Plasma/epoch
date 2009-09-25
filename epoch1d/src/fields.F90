MODULE field
  USE shared_data
  USE boundary
  USE laser

  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_eb_fields_half

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


    REAL(num) :: lx,cnx

    lx=dt/dx
    cnx=0.5_num*lx

    !Update Ex to t=t0+dt/2

#ifndef HIGH_ORDER_FIELDS

    DO ix=1,nx
       Ex(ix)=Ex(ix)&
            -0.5_num*dt*Jx(ix)/epsilon0
    ENDDO

    !Update Ey to t=t0+dt/2
    DO ix=1,nx
       Ey(ix)=Ey(ix)&
            -cnx*(bz(ix)-bz(ix-1))*c**2&
            -0.5_num*dt*Jy(ix)/epsilon0
    ENDDO

    !Update Ez to t=t0+dt/2
    DO ix=1,nx
       Ez(ix)=Ez(ix)&
            +cnx*(by(ix)-by(ix-1))*c**2&
            -0.5_num*dt*Jz(ix)/epsilon0
    ENDDO

#else
   DO ix=1,nx
     	Ex(ix)=Ex(ix)&
		-0.5_num*dt*Jx(ix)/epsilon0
 	ENDDO

   DO ix=1,nx
      Ey(ix)=Ey(ix)&
           -cnx*c**2 * SUM(Diff_Consts * Bz(ix-Large:ix+Small))&
			  -0.5_num*dt*Jy(ix)/epsilon0
   ENDDO

   DO ix=1,nx
      Ez(ix)=Ez(ix)&
           +cnx*c**2 * SUM(Diff_Consts * By(ix-Large:ix+Small))&
           -0.5_num*dt*Jz(ix)/epsilon0
   ENDDO
#endif


    !Now have E(t+dt/2), do boundary conditions on E

    CALL Efield_bcs

    !Update B field to t+dt/2 using E(t+dt/2)

#ifndef HIGH_ORDER_FIELDS
    !Bx unchanged in 1D

    !By
    DO ix=1,nx
       By(ix)=By(ix)&
            +cnx*(Ez(ix+1)-Ez(ix))
    ENDDO

    !Bz
    DO ix=1,nx
       Bz(ix)=Bz(ix)&
            -cnx*(Ey(ix+1)-Ey(ix))
    ENDDO
#else
	DO ix=1,nx
   	By(ix)=By(ix)&
        	+cnx * SUM(Diff_Consts * Ez(ix-Small:ix+Large))
	ENDDO
	
	DO ix=1,nx
   	Bz(ix)=Bz(ix)&
        	-cnx * SUM(Diff_Consts * Ey(ix-Small:ix+Large))
	ENDDO
#endif

    !Now have B field at t+dt/2. Do boundary conditions on B
    CALL Bfield_bcs(.FALSE.)

    !Now have E&B fields at t=t+dt/2
    !Move to particle pusher

  END SUBROUTINE update_eb_fields_half

  SUBROUTINE update_eb_fields_final

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

    REAL(num) :: lx,cnx

    lx=dt/dx
    cnx=0.5_num*lx
#ifndef HIGH_ORDER_FIELDS
    !Bx unchanged in 1D

    !By
    DO ix=1,nx
       By(ix)=By(ix)&
            +cnx*(Ez(ix+1)-Ez(ix))
    ENDDO

    !Bz
    DO ix=1,nx
       Bz(ix)=Bz(ix)&
            -cnx*(Ey(ix+1)-Ey(ix))
    ENDDO
#else
	DO ix=1,nx
		By(ix)=By(ix)&
     		+cnx * SUM(Diff_Consts * Ez(ix-Small:ix+Large))
	ENDDO

	DO ix=1,nx
		Bz(ix)=Bz(ix)&
     		-cnx * SUM(Diff_Consts * Ey(ix-Small:ix+Large))
	ENDDO
#endif


    CALL BField_BCS(.FALSE.)
    IF(xbc_left == BC_SIMPLE_LASER .AND. left == MPI_PROC_NULL) CALL laser_bcs_left
    IF(xbc_left == BC_SIMPLE_OUTFLOW .AND. left == MPI_PROC_NULL) CALL outflow_bcs_left

    IF(xbc_right == BC_SIMPLE_LASER .AND. right == MPI_PROC_NULL) CALL laser_bcs_right
    IF(xbc_right == BC_SIMPLE_OUTFLOW .AND. right == MPI_PROC_NULL) CALL outflow_bcs_right
    CALL BField_BCS(.TRUE.)

#ifndef HIGH_ORDER_FIELDS
    !Ex
    DO ix=1,nx
       Ex(ix)=Ex(ix)&
            -0.5*dt*Jx(ix)/epsilon0
    ENDDO

    !Ey
    DO ix=1,nx
       Ey(ix)=Ey(ix)&
            -cnx*(Bz(ix)-Bz(ix-1)) * c**2&
            -0.5*dt*Jy(ix)/epsilon0
    ENDDO

    !Ez
    DO ix=1,nx
       Ez(ix)=Ez(ix)&
            +cnx*(By(ix)-By(ix-1))*c**2&
            -0.5*dt*Jz(ix)/epsilon0
    ENDDO
#else
	DO ix=1,nx
  		Ex(ix)=Ex(ix)&
		-0.5_num*dt*Jx(ix)/epsilon0
	ENDDO

	DO ix=1,nx
   	Ey(ix)=Ey(ix)&
        	-cnx*c**2 * SUM(Diff_Consts * Bz(ix-Large:ix+Small))&
		  	-0.5_num*dt*Jy(ix)/epsilon0
	ENDDO

	DO ix=1,nx
   	Ez(ix)=Ez(ix)&
        	+cnx*c**2 * SUM(Diff_Consts * By(ix-Large:ix+Small))&
        	-0.5_num*dt*Jz(ix)/epsilon0
	ENDDO
#endif

    CALL Efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE field














