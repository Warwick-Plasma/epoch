MODULE calc_df
  USE shared_data
  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_mass_density(DataArray)

    !Contains the integer cell position of the particle in x,y,z
    INTEGER :: Cell_x,Cell_y,Cell_z

    !Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x,part_y,part_z,part_px,part_py,part_pz,part_q,part_m

    !Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r,cell_y_r,cell_z_r

    !The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x,cell_frac_y,cell_frac_z

    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    !Defined at the particle position

    REAL(num),DIMENSION(-1:1) :: gx,gy
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-1:,-1:),INTENT(INOUT) :: DataArray
    REAL(num),DIMENSION(:,:), ALLOCATABLE :: Temp

    DataArray=0.0_num
    ALLOCATE(Temp(-1:nx+1,-1:ny+1))


    DO ipart=1,npart

       IF (part_species(ipart) == 0) CYCLE

       !Copy the particle properties out for speed
       part_x  = Part_pos(ipart,1) - x_start
       part_y  = Part_pos(ipart,2) - y_start
       part_px = Part_P(ipart,1)
       part_py = Part_P(ipart,2)
       part_pz = Part_P(ipart,3)
       !Use a lookup table for charge and mass to save memory
       !No reason not to do this (I think), check properly later
       part_q  = species(Part_species(ipart),1)
       part_m  = species(Part_species(ipart),2)



       cell_x_r = part_x / dx
       cell_x  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       cell_y_r = part_y / dy
       cell_y  = NINT(cell_y_r)
       cell_frac_y = REAL(cell_y,num) - cell_y_r
       cell_y=cell_y+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

       gy(-1)=0.5_num * (0.5_num + cell_frac_y)**2
       gy(0)=0.75_num - cell_frac_y**2
       gy(1)=0.5_num * (0.5_num - cell_frac_y)**2

       DO iy=-1,1
          DO ix=-1,1
             Data=part_m
             DataArray(cell_x+ix,cell_y+iy) = DataArray(cell_x+ix,cell_y+iy) + &
                  gx(ix) * gy(iy) * Data
          ENDDO
       ENDDO

    ENDDO

    Temp=0.0_num
    CALL MPI_ALLREDUCE(DataArray(1:nx,1:ny),Temp(1:nx,1:ny),nx*ny,mpireal,MPI_SUM,comm,errcode)

    DataArray(-1:nx+1,-1:ny+1)=Temp
    DEALLOCATE(Temp)

    CALL Periodic_Summation_bcs(DataArray)

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_charge_density(DataArray)

    !Contains the integer cell position of the particle in x,y,z
    INTEGER :: Cell_x,Cell_y,Cell_z

    !Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x,part_y,part_z,part_px,part_py,part_pz,part_q,part_m

    !Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r,cell_y_r,cell_z_r

    !The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x,cell_frac_y,cell_frac_z

    !Weighting factors as Eqn 4.77 page 25 of manual
    !Eqn 4.77 would be written as
    !F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    !Defined at the particle position

    REAL(num),DIMENSION(-1:1) :: gx,gy
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: DataArray
    REAL(num),DIMENSION(:,:), ALLOCATABLE :: Temp

    DataArray=0.0_num
    ALLOCATE(Temp(-1:nx+1,-1:ny+1))


    DO ipart=1,npart

       IF (part_species(ipart) == 0) CYCLE

       !Copy the particle properties out for speed
       part_x  = Part_pos(ipart,1) - x_start
       part_y  = Part_pos(ipart,2) - y_start
       part_px = Part_P(ipart,1)
       part_py = Part_P(ipart,2)
       part_pz = Part_P(ipart,3)
       !Use a lookup table for charge and mass to save memory
       !No reason not to do this (I think), check properly later
       part_q  = species(Part_species(ipart),1)
       part_m  = species(Part_species(ipart),2)


       cell_x_r = part_x / dx
       cell_x  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       cell_y_r = part_y / dy
       cell_y  = NINT(cell_y_r)
       cell_frac_y = REAL(cell_y,num) - cell_y_r
       cell_y=cell_y+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

       gy(-1)=0.5_num * (0.5_num + cell_frac_y)**2
       gy(0)=0.75_num - cell_frac_y**2
       gy(1)=0.5_num * (0.5_num - cell_frac_y)**2

       DO iy=-1,1
          DO ix=-1,1
             Data=part_q
             DataArray(cell_x+ix,cell_y+iy) = DataArray(cell_x+ix,cell_y+iy) + &
                  gx(ix) * gy(iy) * Data
          ENDDO
       ENDDO

    ENDDO

    Temp=0.0_num
    CALL MPI_ALLREDUCE(DataArray(1:nx,1:ny),Temp(1:nx,1:ny),nx*ny,mpireal,MPI_SUM,comm,errcode)

    DataArray(-1:nx+1,-1:ny+1)=Temp
    DEALLOCATE(Temp)

    CALL Periodic_Summation_bcs(DataArray)

  END SUBROUTINE calc_charge_density

END MODULE calc_df















