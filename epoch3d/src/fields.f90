MODULE field
  USE shared_data
  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_eb_fields_half

    REAL(num) :: lx, cnx, ly, cny, lz, cnz

    lx=dt/dx
    cnx=0.5_num*lx

    ly=dt/dy
    cny=0.5_num*ly

    lz=dt/dz
    cnz=0.5_num*lz

    !Update Ex to t=t0+dt/2

    DO iz=0,nz
       DO iy=0,ny
          DO ix=0,nx
             Ex(ix,iy,iz)=Ex(ix,iy,iz)&
                  +cny*(bz(ix,iy,iz)-bz(ix,iy-1,iz))*c**2&
                  -cnz*(by(ix,iy,iz)-by(ix,iy,iz-1))*c**2&
                  -0.5_num*dt*Jx(ix,iy,iz)/epsilon0
          ENDDO
       ENDDO
    ENDDO

    !Update Ey to t=t0+dt/2
    DO iz=0,nz
       DO iy=0,ny
          DO ix=0,nx
             Ey(ix,iy,iz)=Ey(ix,iy,iz)&
                  +cnz*(bx(ix,iy,iz)-bx(ix,iy,iz-1))*c**2&
                  -cnx*(bz(ix,iy,iz)-bz(ix-1,iy,iz))*c**2&
                  -0.5_num*dt*Jy(ix,iy,iz)/epsilon0
          ENDDO
       ENDDO
    ENDDO

    !Update Ez to t=t0+dt/2
    DO iz=0,nz
       DO iy=0,ny
          DO ix=0,nx
             Ez(ix,iy,iz)=Ez(ix,iy,iz)&
                  +cnx*(by(ix,iy,iz)-by(ix-1,iy,iz))*c**2&
                  -cny*(bx(ix,iy,iz)-bx(ix,iy-1,iz))*c**2&
                  -0.5_num*dt*Jz(ix,iy,iz)/epsilon0
          ENDDO
       ENDDO
    ENDDO

    !Now have E(t+dt/2), do boundary conditions on E

    CALL Efield_bcs

    !Update B field to t+dt/2 using E(t+dt/2)

    !Bx
    DO iz=0,nz
       DO iy=0,ny
          DO ix=0,nx
             Bx(ix,iy,iz) = Bx(ix,iy,iz)&
                  +cnz*(Ey(ix,iy,iz+1)-Ey(ix,iy,iz))&
                  -cny*(Ez(ix,iy+1,iz)-Ez(ix,iy,iz))
          ENDDO
       ENDDO
    ENDDO

    !By
    DO iz=0,nz
       DO iy=0,ny
          DO ix=0,nx
             By(ix,iy,iz)=By(ix,iy,iz)&
                  +cnx*(Ez(ix+1,iy,iz)-Ez(ix,iy,iz))&
                  -cnz*(Ex(ix,iy,iz+1)-Ex(ix,iy,iz))
          ENDDO
       ENDDO
    ENDDO

    !Bz
    DO iz=0,nz
       DO iy=0,ny
          DO ix=0,nx
             Bz(ix,iy,iz)=Bz(ix,iy,iz)&
                  -cnx*(Ey(ix+1,iy,iz)-Ey(ix,iy,iz))&
                  +cny*(Ex(ix,iy+1,iz)-Ex(ix,iy,iz))
          ENDDO
       ENDDO
    ENDDO

    !Now have B field at t+dt/2. Do boundary conditions on B
    CALL Bfield_bcs

    !Now have E&B fields at t=t+dt/2
    !Move to particle pusher

  END SUBROUTINE update_eb_fields_half

  SUBROUTINE update_eb_fields_final

    REAL(num) :: lx,cnx,ly,cny, lz, cnz

    lx=dt/dx
    cnx=0.5_num*lx

    ly=dt/dy
    cny=0.5_num*ly

    lz=dt/dz
    cnz=0.5_num*lz

    !Bx 
    DO iz=0,nz+1
       DO iy=0,ny+1
          DO ix=0,nx+1
             Bx(ix,iy,iz) = Bx(ix,iy,iz)&
                  +cnz*(Ey(ix,iy,iz+1)-Ey(ix,iy,iz))&
                  -cny*(Ez(ix,iy+1,iz)-Ez(ix,iy,iz))
          ENDDO
       ENDDO
    ENDDO

    !By
    DO iz=0,nz+1
       DO iy=0,ny+1
          DO ix=0,nx+1
             By(ix,iy,iz)=By(ix,iy,iz)&
                  +cnx*(Ez(ix+1,iy,iz)-Ez(ix,iy,iz))&
                  -cnz*(Ex(ix,iy,iz+1)-Ex(ix,iy,iz))
          ENDDO
       ENDDO
    ENDDO

    !Bz
    DO iz=0,nz+1
       DO iy=0,ny+1
          DO ix=0,nx+1
             Bz(ix,iy,iz)=Bz(ix,iy,iz)&
                  -cnx*(Ey(ix+1,iy,iz)-Ey(ix,iy,iz))&
                  +cny*(Ex(ix,iy+1,iz)-Ex(ix,iy,iz))
          ENDDO
       ENDDO
    ENDDO

    CALL Bfield_bcs

    !Ex
    DO iz=0,nz+1
       DO iy=0,ny+1
          DO ix=0,nx+1
             Ex(ix,iy,iz)=Ex(ix,iy,iz)&
                  +cny*(Bz(ix,iy,iz)-Bz(ix,iy-1,iz)) * c**2&
                  -cnz*(By(ix,iy,iz)-By(ix,iy,iz-1)) * c**2&
                  -0.5*dt*Jx(ix,iy,iz)/epsilon0
          ENDDO
       ENDDO
    ENDDO

    !Ey
    DO iz=0,nz+1
       DO iy=0,ny+1
          DO ix=0,nx+1
             Ey(ix,iy,iz)=Ey(ix,iy,iz)&
                  +cnz*(Bx(ix,iy,iz)-Bx(ix,iy,iz-1)) * c **2&
                  -cnx*(Bz(ix,iy,iz)-Bz(ix-1,iy,iz)) * c **2&
                  -0.5*dt*Jy(ix,iy,iz)/epsilon0
          ENDDO
       ENDDO
    ENDDO



    !Ez
    DO iz=0,nz+1
       DO iy=0,ny+1
          DO ix=0,nx+1
             Ez(ix,iy,iz)=Ez(ix,iy,iz)&
                  +cnx*(By(ix,iy,iz)-By(ix-1,iy,iz))*c**2&
                  -cny*(Bx(ix,iy,iz)-Bx(ix,iy-1,iz))*c**2&
                  -0.5*dt*Jz(ix,iy,iz)/epsilon0
          ENDDO
       ENDDO
    ENDDO

    CALL Efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE field














