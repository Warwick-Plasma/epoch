MODULE field
  USE shared_data
  USE boundary

  IMPLICIT NONE

SAVE

CONTAINS

  SUBROUTINE update_eb_fields_half

    REAL(num) :: lx,cnx

    lx=dt/dx
    cnx=0.5_num*lx

    !Update Ex to t=t0+dt/2

    DO ix=0,nx+1
       Ex(ix)=Ex(ix)&
            -0.5_num*dt*Jx(ix)/epsilon0
    ENDDO

    !Update Ey to t=t0+dt/2
    DO ix=0,nx+1
       Ey(ix)=Ey(ix)&
            -cnx*(Bz(ix)-Bz(ix-1))*c**2&
            -0.5_num*dt*Jy(ix)/epsilon0
    ENDDO

    !Update Ez to t=t0+dt/2
    DO ix=0,nx+1
       Ez(ix)=Ez(ix)&
            +cnx*(By(ix)-By(ix-1))*c**2&
            -0.5_num*dt*Jz(ix)/epsilon0
    ENDDO

    !Now have E(t+dt/2), do boundary conditions on E

    CALL Efield_bcs

    !Update B field to t+dt/2 using E(t+dt/2)
    !Bx is unchanged in 1D


    !By
    DO ix=0,nx+1
       By(ix)=By(ix)&
            +cnx*(Ez(ix+1)-Ez(ix))
    ENDDO

    !Bz
    DO ix=0,nx+1
       Bz(ix)=Bz(ix)&
            -cnx*(Ey(ix+1)-Ey(ix))
    ENDDO

    !Now have B field at t+dt/2. Do boundary conditions on B
    CALL Bfield_bcs

    !Now have E&B fields at t=t+dt/2
    !Move to particle pusher

  END SUBROUTINE update_eb_fields_half

  SUBROUTINE update_eb_fields_final

    REAL(num) :: lx,cnx

    lx=dt/dx
    cnx=0.5_num*lx

    !Bx unchanged

    !By
    DO ix=0,nx+1
       By(ix)=By(ix)&
            +cnx*(Ez(ix+1)-Ez(ix))
    ENDDO

    !Bz
    DO ix=0,nx+1
       Bz(ix)=Bz(ix)&
            -cnx*(Ey(ix+1)-Ey(ix))
    ENDDO


    CALL Bfield_bcs

    !Ex
    DO ix=0,nx+1
       Ex(ix)=Ex(ix)&
            -0.5*dt*Jx(ix)/epsilon0
    ENDDO

    !Ey
    DO ix=0,nx+1
       Ey(ix)=Ey(ix)&
            -cnx*(Bz(ix)-Bz(ix-1))*c**2&
            -0.5*dt*Jy(ix)/epsilon0
    ENDDO

    !Ez
    DO ix=0,nx+1
       Ez(ix)=Ez(ix)&
            +cnx*(By(ix)-By(ix-1))*c**2&
            -0.5*dt*Jz(ix)/epsilon0
    ENDDO

    CALL Efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE field














