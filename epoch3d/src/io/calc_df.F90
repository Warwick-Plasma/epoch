MODULE calc_df
  USE shared_data
  USE boundary
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_mass_density(DataArray,CurSpecies)

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

	!Particle Weight factors as described in the manual (FIXREF)
    REAL(num),DIMENSION(-2:2) :: gx,gy,gz
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: DataArray
    INTEGER,INTENT(IN) :: CurSpecies

    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies, spec_start,spec_end

    DataArray=0.0_num

    l_weight=weight
    spec_start=CurSpecies
    spec_end=CurSpecies

    IF (CurSpecies .LE. 0) THEN
       spec_start=1
       spec_end=nSpecies
    ENDIF

    DO iSpecies=spec_start,spec_end
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE (ASSOCIATED(Current))

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_z  = Current%Part_pos(3) - z_start_local
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

          cell_z_r = part_z / dz
          cell_z  = NINT(cell_z_r)
          cell_frac_z = REAL(cell_z,num) - cell_z_r
          cell_z=cell_z+1

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)
			 CALL ParticleToGrid(cell_frac_z,gz)			

          DO iz=-sf_order,sf_order
             DO iy=-sf_order,sf_order
                DO ix=-sf_order,sf_order
                   Data=part_m * l_weight / (dx*dy*dz)
                   DataArray(cell_x+ix,cell_y+iy,cell_z+iz) = DataArray(cell_x+ix,cell_y+iy,cell_z+iz) + &
                        gx(ix) * gy(iy) * gz(iz) * Data
                ENDDO
             ENDDO
          ENDDO

          Current=>Current%Next
       ENDDO
    ENDDO


    CALL Processor_Summation_BCS(DataArray)
  END SUBROUTINE calc_mass_density


 SUBROUTINE calc_ekbar(DataArray,CurSpecies)

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

	!Particle Weight factors as described in the manual (FIXREF)
    REAL(num),DIMENSION(-2:2) :: gx,gy,gz
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: DataArray
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: ct
    INTEGER,INTENT(IN) :: CurSpecies

    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies, spec_start,spec_end

	 ALLOCATE(ct(-2:nx+3,-2:ny+3,-2:nz+3))
    DataArray=0.0_num
	 ct=0.0_num
	

    l_weight=weight

    spec_start=CurSpecies
    spec_end=CurSpecies

    IF (CurSpecies .LE. 0) THEN
       spec_start=1
       spec_end=nSpecies
    ENDIF

    DO iSpecies=spec_start,spec_end
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE (ASSOCIATED(Current))

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_y  = Current%Part_pos(3) - z_start_local
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

          cell_z_r = part_z / dz
          cell_z  = NINT(cell_z_r)
          cell_frac_z = REAL(cell_z,num) - cell_z_r
          cell_z=cell_z+1

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)
			 CALL ParticleToGrid(cell_frac_z,gz)

			DO iz=-sf_order,sf_order
          DO iy=-sf_order,sf_order
             DO ix=-sf_order,sf_order
                Data=SQRT(((part_px*l_weight)**2+(part_py*l_weight)**2+(part_pz*l_weight)**2)*c**2 &
						+ (part_m*l_weight)**2*c**4) - (part_m*l_weight)*c**2
                DataArray(cell_x+ix,cell_y+iy,cell_z+iz) = DataArray(cell_x+ix,cell_y+iy,cell_z+iz) + &
                     gx(ix) * gy(iy) * gz(iz) * Data
	             ct(cell_x+ix,cell_y+iy,cell_z+iz) = ct(cell_x+ix,cell_y+iy,cell_z+iz) + &
	                  gx(ix) * gy(iy) * gz(iz) * l_weight
             ENDDO
          ENDDO
			ENDDO

          Current=>Current%Next
       ENDDO
    ENDDO


    CALL Processor_Summation_BCS(DataArray)
	 CALL Processor_Summation_BCS(ct)
	
	 DataArray = DataArray / MAX(ct,none_zero)
    CALL Field_Zero_Gradient(DataArray,.TRUE.)
	 DEALLOCATE(ct)

  END SUBROUTINE calc_ekbar

  SUBROUTINE calc_charge_density(DataArray,CurSpecies)
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

    REAL(num),DIMENSION(-2:2) :: gx,gy,gz
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: DataArray
    INTEGER,INTENT(IN) :: CurSpecies

    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies, spec_start, spec_end

    DataArray=0.0_num

    l_weight=weight

    spec_start=CurSpecies
    spec_end=CurSpecies

    IF (CurSpecies .LE. 0) THEN
       spec_start=1
       spec_end=nSpecies
    ENDIF

    DO iSpecies=spec_start,spec_end
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE (ASSOCIATED(Current))

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_z  = Current%Part_pos(3) - z_start_local
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

          cell_z_r = part_z / dz
          cell_z  = NINT(cell_z_r)
          cell_frac_z = REAL(cell_z,num) - cell_z_r
          cell_z=cell_z+1

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)
			 CALL ParticleToGrid(cell_frac_z,gz)			

          DO iz=-sf_order,sf_order
             DO iy=-sf_order,sf_order
                DO ix=-sf_order,sf_order
                   Data=part_q * l_weight / (dx*dy*dz)
                   DataArray(cell_x+ix,cell_y+iy,cell_z+iz) = DataArray(cell_x+ix,cell_y+iy,cell_z+iz) + &
                        gx(ix) * gy(iy) * gz(iz) * Data
                ENDDO
             ENDDO
          ENDDO

          Current=>Current%Next
       ENDDO
    ENDDO


    CALL Processor_Summation_BCS(DataArray)


  END SUBROUTINE calc_charge_density

  SUBROUTINE calc_number_density(DataArray,CurSpecies)
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

    REAL(num),DIMENSION(-2:2) :: gx,gy,gz
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: DataArray
    INTEGER,INTENT(IN) :: CurSpecies

    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies, spec_start, spec_end

    DataArray=0.0_num

    l_weight=weight

    spec_start=CurSpecies
    spec_end=CurSpecies

    IF (CurSpecies .LE. 0) THEN
       spec_start=1
       spec_end=nSpecies
    ENDIF

    DO iSpecies=spec_start,spec_end
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE (ASSOCIATED(Current))

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_z  = Current%Part_pos(3) - z_start_local
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

          cell_z_r = part_z / dz
          cell_z  = NINT(cell_z_r)
          cell_frac_z = REAL(cell_z,num) - cell_z_r
          cell_z=cell_z+1

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)
			 CALL ParticleToGrid(cell_frac_z,gz)			

          DO iz=-sf_order,sf_order
             DO iy=-sf_order,sf_order
                DO ix=-sf_order,sf_order
                   Data=l_weight / (dx*dy*dz)
                   DataArray(cell_x+ix,cell_y+iy,cell_z+iz) = DataArray(cell_x+ix,cell_y+iy,cell_z+iz) + &
                        gx(ix) * gy(iy) * gz(iz) * Data
                ENDDO
             ENDDO
          ENDDO

          Current=>Current%Next
       ENDDO
    ENDDO

    CALL Processor_Summation_BCS(DataArray)

  END SUBROUTINE calc_number_density

  SUBROUTINE calc_on_grid_with_evaluator(DataArray,CurSpecies,evaluator)
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

    REAL(num),DIMENSION(-2:2) :: gx,gy,gz
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:,-2:),INTENT(INOUT) :: DataArray
    INTEGER,INTENT(IN) :: CurSpecies

    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies, spec_start, spec_end

    INTERFACE
       FUNCTION evaluator(aParticle,species_eval)
         USE shared_data
         TYPE(particle), POINTER :: aParticle
         INTEGER,INTENT(IN) :: species_eval
         REAL(num) :: evaluator
       END FUNCTION evaluator
    END INTERFACE

    DataArray=0.0_num

    l_weight=weight

    spec_start=CurSpecies
    spec_end=CurSpecies

    IF (CurSpecies .LE. 0) THEN
       spec_start=1
       spec_end=nSpecies
    ENDIF

    DO iSpecies=spec_start,spec_end
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE (ASSOCIATED(Current))

          !Copy the particle properties out for speed
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_z  = Current%Part_pos(3) - z_start_local
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

          cell_z_r = part_z / dz
          cell_z  = NINT(cell_z_r)
          cell_frac_z = REAL(cell_z,num) - cell_z_r
          cell_z=cell_z+1

 			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)
			 CALL ParticleToGrid(cell_frac_z,gz)
			
          DO iz=-sf_order,sf_order
             DO iy=-sf_order,sf_order
                DO ix=-sf_order,sf_order
                   Data=evaluator(Current,iSpecies)
                   DataArray(cell_x+ix,cell_y+iy,cell_z+iz) = DataArray(cell_x+ix,cell_y+iy,cell_z+iz) + &
                        gx(ix) * gy(iy) * gz(iz) * Data
                ENDDO
             ENDDO
          ENDDO

          Current=>Current%Next
       ENDDO
    ENDDO


    CALL Processor_Summation_BCS(DataArray)


  END SUBROUTINE calc_on_grid_with_evaluator

END MODULE calc_df















