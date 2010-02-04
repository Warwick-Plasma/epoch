MODULE solve_gauss

  USE shared_data
  USE multigrid
  USE boundary

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE do_gauss

      !Contains the integer cell position of the particle in x,y,z
      INTEGER :: cell_x,cell_y

      !Properties of the current particle. Copy out of particle arrays for speed
      REAL(num) :: part_x,part_y,part_px,part_py,part_pz,part_q,part_m

      !Contains the floating point version of the cell number (never actually used)
      REAL(num) :: cell_x_r,cell_y_r

      !The fraction of a cell between the particle position and the cell boundary
      REAL(num) :: cell_frac_x,cell_frac_y

      !The weight of a particle
      REAL(num) :: l_weight

      !Weighting factors as Eqn 4.77 page 25 of manual
      !Eqn 4.77 would be written as
      !F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
      !Defined at the particle position

      REAL(num),DIMENSION(-1:1) :: gx,gy
      !The data to be weighted onto the grid
      REAL(num) :: data
      REAL(num),DIMENSION(:,:),ALLOCATABLE :: rho,phi
      TYPE(particle),POINTER :: current
      INTEGER :: ispecies, ngrids_local,ngrids
      LOGICAL :: force_converged=.TRUE.

      ALLOCATE(rho(-2:nx+3,-2:ny+3),phi(-2:nx+3,-2:ny+3))
      rho=0.0_num
      phi=0.0_num

      l_weight=weight
      DO ispecies=1,n_species
         current=>particle_species(ispecies)%attached_list%head
         DO WHILE (ASSOCIATED(current))

            !Copy the particle properties out for speed
            part_x  = current%part_pos(1) - x_start_local
            part_y  = current%part_pos(2) - y_start_local
            part_px = current%part_p(1)
            part_py = current%part_p(2)
            part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
            part_q  = current%charge
            part_m  = current%mass
#else
            part_q  = particle_species(ispecies)%charge
            part_m  = particle_species(ispecies)%mass
#endif

#ifdef PER_PARTICLE_WEIGHT
            l_weight=current%weight
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
                  data=part_q * l_weight / (dx*dy)
                  rho(cell_x+ix,cell_y+iy) = rho(cell_x+ix,cell_y+iy) + &
                       gx(ix) * gy(iy) * data
               ENDDO
            ENDDO
            current=>current%next
         ENDDO
      ENDDO

      CALL processor_summation_bcs(rho)

      !Now know charge density (rho), now use the multigrid scheme

      !Calculate maximum number of allowed grids
      ngrids_local=FLOOR(MIN(LOG(REAL(nx,num)),LOG(REAL(ny,num)))/LOG(2.0_num))-2
      CALL MPI_ALLREDUCE(ngrids_local,ngrids,1,MPI_INTEGER,MPI_MIN,comm,errcode)

      CALL setup_multigrid(nx,ny,ngrids,dx,dy)
      CALL setup_mpi(left,right,up,down,comm)
      CALL run_multigrid(phi(1:nx,1:ny),rho(1:nx,1:ny),30,10,1000,force_converged)
      CALL field_bc(phi)

      DO iy=1,ny
         DO ix=1,nx
            ex(ix,iy)=(phi(ix+1,iy)-phi(ix,iy))/(epsilon0 * dx)
            ey(ix,iy)=(phi(ix,iy+1)-phi(ix,iy))/(epsilon0 * dy)
         ENDDO
      ENDDO

    END SUBROUTINE do_gauss

END MODULE solve_gauss
