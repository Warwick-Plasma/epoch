MODULE field
  USE shared_data
  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_eb_fields_half

    REAL(num) :: lx,cnx,ly,cny

    lx=dt/dx
    cnx=0.5_num*lx

    ly=dt/dy
    cny=0.5_num*ly

    !Update Ex to t=t0+dt/2

    DO iy=0,ny+1
       DO ix=0,nx+1
          Ex(ix,iy)=Ex(ix,iy)&
               +cny*(bz(ix,iy)-bz(ix,iy-1))*c**2&
               -0.5_num*dt*Jx(ix,iy)/epsilon0
       ENDDO
    ENDDO

    !Update Ey to t=t0+dt/2
    DO iy=0,ny+1
       DO ix=0,nx+1
          Ey(ix,iy)=Ey(ix,iy)&
               -cnx*(bz(ix,iy)-bz(ix-1,iy))*c**2&
               -0.5_num*dt*Jy(ix,iy)/epsilon0
       ENDDO
    ENDDO

    !Update Ez to t=t0+dt/2
    DO iy=0,ny+1
       DO ix=0,nx+1
          Ez(ix,iy)=Ez(ix,iy)&
               +cnx*(by(ix,iy)-by(ix-1,iy))*c**2&
               -cny*(bx(ix,iy)-bx(ix,iy-1))*c**2&
               -0.5_num*dt*Jz(ix,iy)/epsilon0
       ENDDO
    ENDDO

    !Now have E(t+dt/2), do boundary conditions on E

    CALL Efield_bcs

    !Update B field to t+dt/2 using E(t+dt/2)

    !Bx
    DO iy=0,ny+1
       DO ix=0,nx+1
          Bx(ix,iy) = Bx(ix,iy)&
               -cny*(Ez(ix,iy+1)-Ez(ix,iy))
       ENDDO
    ENDDO


    !By
    DO iy=0,ny+1
       DO ix=0,nx+1
          By(ix,iy)=By(ix,iy)&
               +cnx*(Ez(ix+1,iy)-Ez(ix,iy))
       ENDDO
    ENDDO

    !Bz
    DO iy=0,ny+1
       DO ix=0,nx+1
          Bz(ix,iy)=Bz(ix,iy)&
               -cnx*(Ey(ix+1,iy)-Ey(ix,iy))&
               +cny*(Ex(ix,iy+1)-Ex(ix,iy))
       ENDDO
    ENDDO

!    Bz=0.0_num

    !Now have B field at t+dt/2. Do boundary conditions on B
    CALL Bfield_bcs

    !Now have E&B fields at t=t+dt/2
    !Move to particle pusher

  END SUBROUTINE update_eb_fields_half

  SUBROUTINE update_eb_fields_final

    REAL(num) :: lx,cnx,ly,cny

    lx=dt/dx
    cnx=0.5_num*lx

    ly=dt/dy
    cny=0.5_num*ly

    !Bx 
    DO iy=0,ny+1
       DO ix=0,nx+1
          Bx(ix,iy) = Bx(ix,iy)&
               -cny*(Ez(ix,iy+1)-Ez(ix,iy))
       ENDDO
    ENDDO

    !By
    DO iy=0,ny+1
       DO ix=0,nx+1
          By(ix,iy)=By(ix,iy)&
               +cnx*(Ez(ix+1,iy)-Ez(ix,iy))
       ENDDO
    ENDDO

    !Bz
    DO iy=0,ny+1
       DO ix=0,nx+1
          Bz(ix,iy)=Bz(ix,iy)&
               -cnx*(Ey(ix+1,iy)-Ey(ix,iy))&
               +cny*(Ex(ix,iy+1)-Ex(ix,iy))
       ENDDO
    ENDDO


!    Bz=0.0_num
    CALL Bfield_bcs

    !Ex
    DO iy=0,ny+1
       DO ix=0,nx+1
          Ex(ix,iy)=Ex(ix,iy)&
               +cny*(Bz(ix,iy)-Bz(ix,iy-1)) * c**2&
               -0.5*dt*Jx(ix,iy)/epsilon0
       ENDDO
    ENDDO

    !Ey
    DO iy=0,ny+1
       DO ix=0,nx+1
          Ey(ix,iy)=Ey(ix,iy)&
               -cnx*(Bz(ix,iy)-Bz(ix-1,iy)) * c **2&
               -0.5*dt*Jy(ix,iy)/epsilon0
       ENDDO
    ENDDO



    !Ez
    DO iy=0,ny+1
       DO ix=0,nx+1
          Ez(ix,iy)=Ez(ix,iy)&
               +cnx*(By(ix,iy)-By(ix-1,iy))*c**2&
               -cny*(Bx(ix,iy)-Bx(ix,iy-1))*c**2&
               -0.5*dt*Jz(ix,iy)/epsilon0
       ENDDO
    ENDDO

    CALL Efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE field














