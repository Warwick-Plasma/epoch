MODULE Solve_Gauss

  USE shared_data
  USE multigrid
  USE boundary

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE Do_Gauss

      !Contains the integer cell position of the particle in x,y,z
      INTEGER :: Cell_x,Cell_y,Cell_z

      !Properties of the current particle. Copy out of particle arrays for speed
      REAL(num) :: part_x,part_y,part_z,part_px,part_py,part_pz,part_q,part_m

      !Contains the floating point version of the cell number (never actually used)
      REAL(num) :: cell_x_r,cell_y_r,cell_z_r

      !The fraction of a cell between the particle position and the cell boundary
      REAL(num) :: cell_frac_x,cell_frac_y,cell_frac_z

      !The weight of a particle
      REAL(num) :: l_weight

      !Weighting factors as Eqn 4.77 page 25 of manual
      !Eqn 4.77 would be written as
      !F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
      !Defined at the particle position

      REAL(num),DIMENSION(-1:1) :: gx,gy
      !The data to be weighted onto the grid
      REAL(num) :: Data
      REAL(num),DIMENSION(:,:),ALLOCATABLE :: rho,phi
      TYPE(Particle),POINTER :: Current
      INTEGER :: iSpecies, nGrids_local,nGrids
      LOGICAL :: force_converged=.TRUE.

      ALLOCATE(rho(-2:nx+3,-2:ny+3),phi(-2:nx+3,-2:ny+3))
      rho=0.0_num
      phi=0.0_num

      l_weight=weight
      DO iSpecies=1,nspecies
         Current=>ParticleSpecies(iSpecies)%AttachedList%Head
         DO WHILE (ASSOCIATED(Current))

            !Copy the particle properties out for speed
            part_x  = Current%Part_pos(1) - x_start_local
            part_y  = Current%Part_pos(2) - y_start_local
            part_px = Current%Part_P(1)
            part_py = Current%Part_P(2)
            part_pz = Current%Part_P(3)
#ifdef PER_PARTICLE_CHARGEMASS
            part_q  = Current%Charge
            part_m  = Current%Mass
#else
            part_q  = ParticleSpecies(iSpecies)%Charge
            part_m  = ParticleSpecies(iSpecies)%Mass
#endif

#ifdef PER_PARTICLE_WEIGHT
            l_weight=Current%Weight
#endif

            cell_x_r = part_x / dx
            cell_x  = NINT(cell_x_r)
            cell_frac_x = REAL(cell_x,num) - cell_x_r
            cell_x=cell_x+1

            cell_y_r = part_y / dy
            cell_y  = NINT(cell_y_r)
            cell_frac_y = REAL(cell_y,num) - cell_y_r
            cell_y=cell_y+1



            gx(-1) = 0.5_num * (1.5_num - ABS(cell_frac_x - 1.0_num))**2
            gx( 0) = 0.75_num - ABS(cell_frac_x)**2
            gx( 1) = 0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2

            gy(-1) = 0.5_num * (1.5_num - ABS(cell_frac_y - 1.0_num))**2
            gy( 0) = 0.75_num - ABS(cell_frac_y)**2
            gy( 1) = 0.5_num * (1.5_num - ABS(cell_frac_y + 1.0_num))**2

            DO iy=-1,1
               DO ix=-1,1
                  Data=part_q * l_weight / (dx*dy)
                  rho(cell_x+ix,cell_y+iy) = rho(cell_x+ix,cell_y+iy) + &
                       gx(ix) * gy(iy) * Data
               ENDDO
            ENDDO
            Current=>Current%Next
         ENDDO
      ENDDO

      CALL Processor_Summation_BCS(rho)

      !Now know charge density (rho), now use the multigrid scheme

      !Calculate maximum number of allowed grids
      nGrids_local=FLOOR(MIN(LOG(REAL(nx,num)),LOG(REAL(ny,num)))/LOG(2.0_num))-2
      CALL MPI_ALLREDUCE(nGrids_local,nGrids,1,MPI_INTEGER,MPI_MIN,comm,errcode)

      CALL Setup_Multigrid(nx,ny,nGrids,dx,dy)
      CALL Setup_MPI(left,right,up,down,comm)
      CALL Run_Multigrid(phi(1:nx,1:ny),rho(1:nx,1:ny),30,10,1000,force_converged)
      CALL Field_BC(phi)

      DO iy=1,ny
         DO ix=1,nx
            Ex(ix,iy)=(phi(ix+1,iy)-phi(ix,iy))/(epsilon0 * dx)
            Ey(ix,iy)=(phi(ix,iy+1)-phi(ix,iy))/(epsilon0 * dy)
         ENDDO
      ENDDO

    END SUBROUTINE Do_Gauss

END MODULE Solve_Gauss
