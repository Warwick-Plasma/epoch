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

    REAL(num),DIMENSION(-1:1) :: gx
    !The data to be weighted onto the grid
    REAL(num) :: Data,Spt,l_weight

    REAL(num),DIMENSION(-1:),INTENT(INOUT) :: DataArray
    REAL(num),DIMENSION(:), ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Cur

    DataArray=0.0_num

    Cur=>MainRoot%Head
#ifndef PER_PARTICLE_WEIGHT
    l_weight=weight
#endif
    DO WHILE (ASSOCIATED(Cur))

       !Copy the particle properties out for speed
       part_x  = Cur%Part_pos - x_start_local
       part_px = Cur%Part_P(1)
       part_py = Cur%Part_P(2)
       part_pz = Cur%Part_P(3)
#ifdef PER_PARTICLE_CHARGEMASS
       part_q  = Cur%Charge
       part_m  = Cur%Mass
#else
       part_q  = species(Cur%Part_species,1)
       part_m  = species(Cur%Part_species,2)
#endif

#ifdef PER_PARTICLE_WEIGHT
       l_weight=Cur%Weight
#endif

       cell_x_r = part_x / dx
       cell_x  = NINT(cell_x_r) 
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

       DO ix=-1,1
          Data=part_m / dx*l_weight !Times weight to get real particle density, divide by dx to give real density, not density per cell
          DataArray(cell_x+ix) = DataArray(cell_x+ix) + &
               gx(ix) * Data
       ENDDO

       Cur=>Cur%Next
    ENDDO

    IF (domain == DO_FULL) THEN
       CALL Field_Reduction(DataArray)
       CALL Periodic_Summation_BCS(DataArray)
    ELSE
       CALL Processor_Summation_BCS(DataArray)
    ENDIF

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

    REAL(num),DIMENSION(-1:1) :: gx
    !The data to be weighted onto the grid
    REAL(num) :: Data,Spt,l_weight

    REAL(num),DIMENSION(-1:),INTENT(INOUT) :: DataArray
    REAL(num),DIMENSION(:), ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Cur

    DataArray=0.0_num


    Cur=>MainRoot%Head
#ifndef PER_PARTICLE_WEIGHT
    l_weight=weight
#endif
    DO WHILE (ASSOCIATED(Cur))

       !Copy the particle properties out for speed
       part_x  = Cur%Part_pos - x_start_local
       part_px = Cur%Part_P(1)
       part_py = Cur%Part_P(2)
       part_pz = Cur%Part_P(3)

#ifdef PER_PARTICLE_CHARGEMASS
       part_q  = Cur%Charge
       part_m  = Cur%Mass
#else
       part_q  = species(Cur%Part_species,1)
       part_m  = species(Cur%Part_species,2)
#endif

#ifdef PER_PARTICLE_WEIGHT
       l_weight=Cur%Weight
#endif


       cell_x_r = part_x / dx
       cell_x  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2


       DO ix=-1,1
          Data=part_q / dx * l_weight
          DataArray(cell_x+ix) = DataArray(cell_x+ix) + &
               gx(ix) * Data
       ENDDO
       Cur=>Cur%Next
    ENDDO

    IF (domain == DO_FULL) THEN
       CALL Field_Reduction(DataArray)
       CALL Periodic_Summation_BCS(DataArray)
    ELSE
       CALL Processor_Summation_BCS(DataArray)
    ENDIF

  END SUBROUTINE calc_charge_density


  SUBROUTINE calc_number_density(DataArray,iSpecies)

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

    REAL(num),DIMENSION(-1:1) :: gx
    !The data to be weighted onto the grid
    REAL(num) :: Data,Spt,l_weight

    REAL(num),DIMENSION(-1:),INTENT(INOUT) :: DataArray
    INTEGER,INTENT(IN) :: iSpecies
    REAL(num),DIMENSION(:), ALLOCATABLE :: Temp
    TYPE(Particle),POINTER :: Cur

    DataArray=0.0_num

    Cur=>MainRoot%Head
#ifndef PER_PARTICLE_WEIGHT
    l_weight=weight
#endif
    DO WHILE (ASSOCIATED(Cur))

       !Copy the particle properties out for speed
       part_x  = Cur%Part_pos - x_start_local
       part_px = Cur%Part_P(1)
       part_py = Cur%Part_P(2)
       part_pz = Cur%Part_P(3)
#ifdef PER_PARTICLE_CHARGEMASS
       part_q  = Cur%Charge
       part_m  = Cur%Mass
#else
       part_q  = species(Cur%Part_species,1)
       part_m  = species(Cur%Part_species,2)
#endif

#ifdef PER_PARTICLE_WEIGHT
       l_weight=Cur%Weight
#endif

       cell_x_r = part_x / dx
       cell_x  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

       DO ix=-1,1
          IF (Cur%Part_Species .EQ. iSpecies) THEN
             Data=l_weight/dx !Times weight to get real particle density, divide by dx to give real density, not density per cell
          ELSE
             Data=0.0_num
          ENDIF
          DataArray(cell_x+ix) = DataArray(cell_x+ix) + &
               gx(ix) * Data
       ENDDO

       Cur=>Cur%Next
    ENDDO

    IF (domain == DO_FULL) THEN
       CALL Field_Reduction(DataArray)
       CALL Periodic_Summation_BCS(DataArray)
    ELSE
       CALL Processor_Summation_BCS(DataArray)
    ENDIF

  END SUBROUTINE calc_number_density

END MODULE calc_df















