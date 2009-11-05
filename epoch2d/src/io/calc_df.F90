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
    REAL(num),DIMENSION(-2:2) :: gx,gy
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: DataArray
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

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)

          DO iy=-sf_order,sf_order
             DO ix=-sf_order,sf_order
                Data=part_m * l_weight / (dx*dy)
                DataArray(cell_x+ix,cell_y+iy) = DataArray(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * Data
             ENDDO
          ENDDO

          Current=>Current%Next
       ENDDO
    ENDDO


    CALL Processor_Summation_BCS(DataArray)
    CALL Field_Zero_Gradient(DataArray,.TRUE.)

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
    REAL(num),DIMENSION(-2:2) :: gx,gy
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: DataArray
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: ct
    INTEGER,INTENT(IN) :: CurSpecies

    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies, spec_start,spec_end

	 ALLOCATE(ct(-2:nx+3,-2:ny+3))
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

          IF (cell_y .LT. -2) PRINT *,cell_y,rank,part_y

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)

          DO iy=-sf_order,sf_order
             DO ix=-sf_order,sf_order
                Data=SQRT(((part_px*l_weight)**2+(part_py*l_weight)**2+(part_pz*l_weight)**2)*c**2 &
						+ (part_m*l_weight)**2*c**4) - (part_m*l_weight)*c**2
                DataArray(cell_x+ix,cell_y+iy) = DataArray(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * Data
	             ct(cell_x+ix,cell_y+iy) = ct(cell_x+ix,cell_y+iy) + &
	                  gx(ix) * gy(iy) * l_weight
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

	!Particle Weight factors as described in the manual (FIXREF)
    REAL(num),DIMENSION(-2:2) :: gx,gy
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: DataArray
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

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)

          DO iy=-sf_order,sf_order
             DO ix=-sf_order,sf_order
                Data=part_q * l_weight / (dx*dy)
                DataArray(cell_x+ix,cell_y+iy) = DataArray(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * Data
             ENDDO
          ENDDO

          Current=>Current%Next
       ENDDO
    ENDDO

    CALL Processor_Summation_BCS(DataArray)
    CALL Field_Zero_Gradient(DataArray,.TRUE.)

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

	!Particle Weight factors as described in the manual (FIXREF)
    REAL(num),DIMENSION(-2:2) :: gx,gy
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: DataArray
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


			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)

          DO iy=-sf_order,sf_order
             DO ix=-sf_order,sf_order
                Data=1.0_num * l_weight / (dx*dy)
                DataArray(cell_x+ix,cell_y+iy) = DataArray(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * Data
             ENDDO
          ENDDO
          Current=>Current%Next
       ENDDO
    ENDDO

    CALL Processor_Summation_BCS(DataArray)
    CALL Field_Zero_Gradient(DataArray,.TRUE.)

  END SUBROUTINE calc_number_density

  SUBROUTINE calc_temperature(DataArray,CurSpecies)

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
    REAL(num),DIMENSION(-2:2) :: gx,gy
    !The data to be weighted onto the grid
    REAL(num) :: Data

    REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: DataArray
    REAL(num),DIMENSION(:,:),ALLOCATABLE ::  part_count, p_max, p_min, mass, sigma, mean
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

!!$    ALLOCATE(p_max(-2:nx+3,-2:ny+3),p_min(-2:nx+3,-2:ny+3))
    ALLOCATE(mean(-2:nx+3,-2:ny+3),part_count(-2:nx+3,-2:ny+3))
    ALLOCATE(mass(-2:nx+3,-2:ny+3), sigma(-2:nx+3,-2:ny+3))
    DataArray=0.0_num
    mean=0.0_num
    part_count=0.0_num
    mass=0.0_num
    sigma=0.0_num

!!$    p_min=1.0e13_num
!!$    p_max=-1.0e13_num

    DO iSpecies=spec_start,spec_end
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
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

          cell_x_r = part_x / dx !
          cell_x  = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x,num) - cell_x_r
          cell_x=cell_x+1

          cell_y_r = part_y / dy !
          cell_y  = NINT(cell_y_r)
          cell_frac_y = REAL(cell_y,num) - cell_y_r
          cell_y=cell_y+1

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)

          DO iy=-sf_order,sf_order
             DO ix=-sf_order,sf_order
                Data=SQRT(part_px**2+part_py**2+part_pz**2) * l_weight
!!$                IF (Data .GE. p_min(cell_x+ix,cell_y+iy) .AND. Data .LE. p_max(cell_x+ix,cell_y+iy)) THEN
                mean(cell_x+ix,cell_y+iy) = mean(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * Data
                Data=l_weight
                part_count(cell_x+ix,cell_y+iy) = part_count(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * Data
                Data=part_m * l_weight
                mass(cell_x+ix,cell_y+iy) = mass(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * Data
!!$                ENDIF
             ENDDO
          ENDDO
          Current=>Current%Next
       ENDDO
    ENDDO

    CALL Processor_Summation_BCS(mean)
    CALL Field_Zero_Gradient(mean,.TRUE.)
    CALL Processor_Summation_BCS(Part_Count)
    CALL Field_Zero_Gradient(Part_Count,.TRUE.)
    CALL Processor_Summation_BCS(mass)
    CALL Field_Zero_Gradient(mass,.TRUE.)

    mean=mean/part_count
    mass=mass/part_count

    WHERE (part_count .LT. 1.0_num) mean=0.0_num
    WHERE (part_count .LT. 1.0_num) mass=0.0_num

    part_count=0.0_num
    DO iSpecies=spec_start,spec_end
       Current=>ParticleSpecies(iSpecies)%AttachedList%Head
       DO WHILE(ASSOCIATED(Current))
          part_x  = Current%Part_pos(1) - x_start_local
          part_y  = Current%Part_pos(2) - y_start_local
          part_px = Current%Part_P(1)
          part_py = Current%Part_P(2)
          part_pz = Current%Part_P(3)

          cell_x_r = part_x / dx !
          cell_x  = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x,num) - cell_x_r
          cell_x=cell_x+1

          cell_y_r = part_y / dy !
          cell_y  = NINT(cell_y_r)
          cell_frac_y = REAL(cell_y,num) - cell_y_r
          cell_y=cell_y+1

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)

          DO iy=-sf_order,sf_order
             DO ix=-sf_order,sf_order
                Data=SQRT(part_px**2+part_py**2+part_pz**2)
!!$                IF (Data .GE. p_min(cell_x+ix,cell_y+iy) .AND. Data .LE. p_max(cell_x+ix,cell_y+iy)) THEN
                sigma(cell_x+ix,cell_y+iy) = sigma(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy) * (Data-mean(cell_x+ix,cell_y+iy))**2
                part_count(cell_x+ix,cell_y+iy) = part_count(cell_x+ix,cell_y+iy) + &
                     gx(ix) * gy(iy)
!!$                ENDIF
             ENDDO
          ENDDO
          Current=>Current%Next
       ENDDO
    ENDDO
    CALL Processor_Summation_BCS(sigma)
    CALL Field_Zero_Gradient(sigma,.TRUE.)
    CALL Processor_Summation_BCS(Part_Count)
    CALL Field_Zero_Gradient(Part_Count,.TRUE.)
    sigma=sigma/MAX(part_count,0.1_num)
    WHERE (part_count .LT. 1.0) sigma=0.0_num

    DataArray=sigma/(kb * MAX(mass,M0))

    DEALLOCATE(part_count,mean,mass,sigma)
!!$    DEALLOCATE(p_min,p_max)

    CALL Field_BC(DataArray)
    CALL Field_Zero_Gradient(DataArray,.TRUE.)

  END SUBROUTINE calc_temperature


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

	!Particle Weight factors as described in the manual (FIXREF)
       REAL(num),DIMENSION(-2:2) :: gx,gy
       !The data to be weighted onto the grid
       REAL(num) :: Data

       REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: DataArray
       INTEGER,INTENT(IN) :: CurSpecies

       TYPE(Particle),POINTER :: Current
       INTEGER :: iSpecies, spec_start,spec_end

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

             part_x  = Current%Part_pos(1) - x_start_local
             part_y  = Current%Part_pos(2) - y_start_local

             cell_x_r = part_x / dx !
             cell_x  = NINT(cell_x_r)
             cell_frac_x = REAL(cell_x,num) - cell_x_r
             cell_x=cell_x+1

             cell_y_r = part_y / dy !
             cell_y  = NINT(cell_y_r)
             cell_frac_y = REAL(cell_y,num) - cell_y_r
             cell_y=cell_y+1

			 CALL ParticleToGrid(cell_frac_x,gx)
			 CALL ParticleToGrid(cell_frac_y,gy)

             DO iy=-sf_order,sf_order
                DO ix=-sf_order,sf_order
                   Data=evaluator(Current,iSpecies)
                   DataArray(cell_x+ix,cell_y+iy) = DataArray(cell_x+ix,cell_y+iy) + &
                        gx(ix) * gy(iy) * Data
                ENDDO
             ENDDO

             Current=>Current%Next
          ENDDO
       ENDDO

       CALL Processor_Summation_BCS(DataArray)
       CALL Field_Zero_Gradient(DataArray,.TRUE.)

     END SUBROUTINE calc_on_grid_with_evaluator

END MODULE calc_df















