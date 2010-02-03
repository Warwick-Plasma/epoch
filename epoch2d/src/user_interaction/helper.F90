MODULE helper

  USE shared_data
  USE partlist
  USE boundary
  USE shape_functions
  IMPLICIT NONE

  SAVE
  REAL(num) :: max_rand

CONTAINS

  SUBROUTINE AutoLoad
    INTEGER :: iSpecies
    INTEGER :: clock,idum
    TYPE(ParticleFamily),POINTER :: PartFam

    CALL SYSTEM_CLOCK(clock)
    idum=-(clock+rank+1)
    DO iSpecies=1,nspecies
       PartFam=>ParticleSpecies(iSpecies)
       IF (move_window) THEN
          ParticleSpecies(iSpecies)%Density=InitialConditions(iSpecies)%Rho(nx,:)
          ParticleSpecies(iSpecies)%Temperature=InitialConditions(iSpecies)%Temp(nx,:,:)
       ENDIF
#ifdef PER_PARTICLE_WEIGHT
		       CALL SetupParticleDensity(InitialConditions(iSpecies)%Rho,Partfam,InitialConditions(iSpecies)%MinRho,InitialConditions(iSpecies)%MaxRho,idum)
#else
		       CALL NonUniformLoadParticles(InitialConditions(iSpecies)%Rho,Partfam,InitialConditions(iSpecies)%MinRho,InitialConditions(iSpecies)%MaxRho,idum)
#endif
       CALL SetupParticleTemperature(InitialConditions(iSpecies)%Temp(:,:,1)&
			,DIR_X,PartFam,InitialConditions(iSpecies)%drift,idum)
       CALL SetupParticleTemperature(InitialConditions(iSpecies)%Temp(:,:,2)&
			,DIR_Y,PartFam,InitialConditions(iSpecies)%drift,idum)
       CALL SetupParticleTemperature(InitialConditions(iSpecies)%Temp(:,:,3)&
			,DIR_Z,PartFam,InitialConditions(iSpecies)%drift,idum)
    ENDDO
  END SUBROUTINE AutoLoad

  SUBROUTINE AllocateIC
    INTEGER :: iSpecies

    ALLOCATE(InitialConditions(1:nspecies))
    DO iSpecies=1,nSpecies
       ALLOCATE(InitialConditions(iSpecies)%Rho(-2:nx+3,-2:ny+3))
       ALLOCATE(InitialConditions(iSpecies)%Temp(-2:nx+3,-2:nx+3,1:3))

       InitialConditions(iSpecies)%Rho=1.0_num
       InitialConditions(iSpecies)%Temp=0.0_num
       InitialConditions(iSpecies)%MinRho=0.0_num
       InitialConditions(iSpecies)%MaxRho=0.0_num
		 InitialConditions(iSpecies)%drift=0.0_num
    ENDDO
    Ex=0.0_num
    Ey=0.0_num
    Ez=0.0_num

    Bx=0.0_num
    By=0.0_num
    Bz=0.0_num
  END SUBROUTINE AllocateIC

  SUBROUTINE DeallocateIC
    INTEGER :: iSpecies, ix, iy
	 REAL(num) :: min_dt, omega, k_max

	 min_dt=1000000.0_num
	 k_max=2.0_num * pi / MIN(dx,dy)
	 !Identify the plasma frequency (Bohm-Gross dispersion relation)
	 !Note that this doesn't get strongly relativistic plasmas right
	DO iSpecies=1,nSpecies
   	DO iy=1,ny
      	DO ix=1,nx
         	omega=SQRT((InitialConditions(iSpecies)%Rho(ix,iy) * Q0**2)/(ParticleSpecies(iSpecies)%Mass * epsilon0) + &
              	3.0_num * k_max**2 * kb * MAXVAL(InitialConditions(iSpecies)%Temp(ix,iy,:))/(ParticleSpecies(iSpecies)%Mass))
         	IF (2.0_num * pi/omega .LT. min_dt) min_dt=2.0_num * pi /omega
      	ENDDO
   	ENDDO
	ENDDO
	CALL MPI_ALLREDUCE(min_dt,dt_plasma_frequency,1,mpireal,MPI_MIN,comm,errcode)
	!Must resolve plasma frequency
	dt_plasma_frequency=dt_plasma_frequency/2.0_num

    DO iSpecies=1,nspecies
       DEALLOCATE(InitialConditions(iSpecies)%Rho)
       DEALLOCATE(InitialConditions(iSpecies)%Temp)
    ENDDO
    DEALLOCATE(InitialConditions)
  END SUBROUTINE DeallocateIC

  SUBROUTINE NonUniformLoadParticles(Density,SpeciesList,minrho,maxrho,idum)

    REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: Density
    TYPE(ParticleFamily),POINTER :: SpeciesList
    REAL(num),INTENT(INOUT) :: minrho,maxrho
    INTEGER,INTENT(INOUT) :: idum

    INTEGER(KIND=8) :: num_valid_cells, num_valid_cells_global
    INTEGER(KIND=8) :: npart_per_cell_average
    INTEGER(KIND=8) :: npart_per_cell
    REAL(num) :: Density_Total, Density_Total_global, Density_Average
    INTEGER(KIND=8) :: npart_this_proc_new, ipart, npart_this_species
    REAL(num) :: rpos

    TYPE(ParticleList),POINTER :: PartList
    TYPE(Particle),POINTER :: Current, Next

    PartList=>SpeciesList%AttachedList

    num_valid_cells=0
    Density_Total=0.0_num
    DO iy=1,ny
       DO ix=1,nx
          IF (Density(ix,iy) .GE. minrho) THEN
             num_valid_cells=num_valid_cells+1
             Density_Total=Density_Total + Density(ix,iy)
          ENDIF
          IF (Density(ix,iy) .GT. maxrho .AND. maxrho .GT. 0.0_num) Density(ix,iy) = maxrho
       ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(num_valid_cells,num_valid_cells_global,1,MPI_INTEGER8,MPI_MAX,comm,errcode)
    npart_per_cell_average=SpeciesList%Count/num_valid_cells_global
    IF (npart_per_cell_average == 0) npart_per_cell_average=1

    CALL MPI_ALLREDUCE(Density_Total,Density_Total_global,1,mpireal,MPI_SUM,comm,errcode)
    Density_Average=Density_Total_global/REAL(num_valid_cells_global,num)

    !Assume that a cell with the average density has the average number of particles per cell
    !Now calculate the new minimum density 
    minrho = Density_Average/REAL(npart_per_cell_average,num)
    !Set the particle weight
    weight = minrho * dx * dy

    !Recalculate the number of valid cells and the summed density
    num_valid_cells=0
    Density_Total=0.0_num
    DO iy=1,ny
       DO ix=1,nx
          IF (Density(ix,iy) .GE. minrho) THEN
             num_valid_cells=num_valid_cells+1
             Density_Total=Density_Total + Density(ix,iy)
          ENDIF
       ENDDO
    ENDDO

    npart_this_proc_new=Density_Total/Density_Average * REAL(npart_per_cell_average,num)

    CALL Destroy_PartList(PartList)
    CALL Create_Allocated_PartList(PartList,npart_this_proc_new)
    Current=>PartList%Head
    DO ix=1,nx
       DO iy=1,ny
          ipart=0
          npart_per_cell = Density(ix,iy)/Density_Average * REAL(npart_per_cell_average,num)
          DO WHILE(ASSOCIATED(Current) .AND. ipart .LT. npart_per_cell)
#ifdef PER_PARTICLE_CHARGEMASS
             !Even if particles have per particle charge and mass, assume that 
             !initially they all have the same charge and mass (user can easily override)
             Current%Charge=SpeciesList%Charge
             Current%Mass=SpeciesList%Mass
#endif
             rpos=random(idum)-0.5_num
             rpos=(rpos*dx)+x(ix)
             Current%Part_Pos(1)=rpos
             rpos=random(idum)-0.5_num
             rpos=(rpos*dy)+y(iy)
             Current%Part_Pos(2)=rpos
             ipart=ipart+1
             Current=>Current%Next
          ENDDO
       ENDDO
    ENDDO

    DO WHILE(ASSOCIATED(Current))
       Next=>Current%Next
       CALL Remove_Particle_From_PartList(PartList,Current)
       DEALLOCATE(Current)
       Current=>Next
    ENDDO
    CALL MPI_REDUCE(PartList%Count,npart_this_species,1,MPI_INTEGER8,MPI_SUM,0,comm,errcode)
    SpeciesList%Count=npart_this_species
    IF (rank .EQ. 0) THEN
       PRINT *,"Loaded",npart_this_species,"particles of species ",TRIM(SpeciesList%Name)
       WRITE(20,*) "Loaded",npart_this_species,"particles of species ",TRIM(SpeciesList%Name)
    ENDIF
    CALL Particle_BCS


  END SUBROUTINE NonUniformLoadParticles

  !This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE LoadParticles(SpeciesList,loadlist,idum)

    INTEGER, INTENT(INOUT) :: idum
    TYPE(ParticleFamily), POINTER,INTENT(INOUT) :: SpeciesList
    LOGICAL,DIMENSION(-2:,-2:),INTENT(IN) :: LoadList
    TYPE(ParticleList),POINTER :: PartList
    TYPE(particle),POINTER :: Current,Next
    INTEGER(KIND=8) :: ipart,npart_per_cell,ipart_total,ncell_per_part
    INTEGER(KIND=8) :: num_valid_cells,num_valid_cells_local,npart_this_species
    INTEGER(KIND=8) :: num_new_particles,npart_left
    REAL(num) :: valid_cell_frac
    REAL(dbl) :: rpos
    INTEGER :: upper_x, upper_y, lower_x, lower_y, cell_x,cell_y
    REAL(num) :: cell_x_r,cell_frac_x
    REAL(num) :: cell_y_r,cell_frac_y
    INTEGER(KIND=8) :: i
    INTEGER :: j


    upper_x=nx
    upper_y=ny
    lower_x=1
    lower_y=1

    IF (Coordinates(2) .EQ. nprocx-1) upper_x=nx+1
    IF (Coordinates(2) .EQ. 0) lower_x=0
    IF (Coordinates(1) .EQ. nprocy-1) upper_y=ny+1
    IF (Coordinates(1) .EQ. 0) lower_y=0

    PartList=>SpeciesList%AttachedList

    npart_this_species=SpeciesList%Count
    IF (npart_this_species .LT. 0) THEN
       IF (rank .EQ. 0) PRINT *,"Unable to continue, species ",TRIM(SpeciesList%Name)," has not had a number of particles set"
       CALL MPI_ABORT(comm,errcode)
    ENDIF
    IF (npart_this_species .EQ. 0) RETURN
    num_valid_cells_local=0
    DO iy=lower_y,upper_y
       DO ix=lower_x,upper_x
          IF (loadlist(ix,iy)) num_valid_cells_local=num_valid_cells_local+1
       ENDDO
    ENDDO

    !Calculate global number of particles per cell
    CALL MPI_ALLREDUCE(num_valid_cells_local,num_valid_cells,1,MPI_INTEGER8,MPI_SUM,comm,errcode)

    IF (num_valid_cells .EQ. 0) THEN
       IF (rank .EQ. 0) THEN
          PRINT *,"Intial condition settings mean that there are no cells where particles may validly be placed for at least one species. Code terminates."
          CALL MPI_ABORT(comm,errcode)
       ENDIF
    ENDIF


    valid_cell_frac=REAL(num_valid_cells_local,num)/REAL(num_valid_cells,num)
    num_new_particles=INT(npart_this_species*valid_cell_frac,KIND=8)
    CALL Destroy_PartList(PartList)
    CALL Create_Allocated_PartList(PartList,num_new_particles)


    npart_per_cell=npart_this_species/num_valid_cells
    ncell_per_part=num_valid_cells/npart_this_species
    SpeciesList%npart_per_cell=npart_per_cell
    IF (SpeciesList%npart_per_cell .EQ. 0) SpeciesList%npart_per_cell=1

    ipart=0
    ipart_total=0
    ix=1
    iy=1

    npart_left=npart_this_species
    Current=>PartList%Head
    IF (npart_per_cell .GT. 0) THEN

       DO ix=lower_x,upper_x
          DO iy=lower_y,upper_y
             ipart=0
             IF (loadlist(ix,iy)) THEN
                DO WHILE(ASSOCIATED(Current) .AND. ipart .LT. npart_per_cell)
#ifdef PER_PARTICLE_CHARGEMASS
                   !Even if particles have per particle charge and mass, assume that 
                   !initially they all have the same charge and mass (user can easily override)
                   Current%Charge=SpeciesList%Charge
                   Current%Mass=SpeciesList%Mass
#endif
                   rpos=random(idum)-0.5_num
                   rpos=(rpos*dx)+x(ix)
                   Current%Part_Pos(1)=rpos
                   rpos=random(idum)-0.5_num
                   rpos=(rpos*dy)+y(iy)
                   Current%Part_Pos(2)=rpos
                   ipart=ipart+1
                   Current=>Current%Next
                   !One particle sucessfully placed
                   npart_left=npart_left-1
                ENDDO
             ENDIF
          ENDDO
       ENDDO

    ENDIF


    DO i=1,npart_left*4
       IF (.NOT. ASSOCIATED(Current)) EXIT
       DO j=1,200
          rpos=random(idum)*(x(nx)-x(1) + dx) - dx/2.0_num
          Current%Part_Pos(1)=rpos+x(1)
          rpos=random(idum)*(y(ny)-y(1) + dy) - dy/2.0_num
          Current%Part_Pos(2)=rpos+y(1)
          cell_x_r = (Current%Part_Pos(1)-x_start_local)/dx - 0.5_num
          cell_x = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x,num) - cell_x_r
          cell_x=cell_x+1

          cell_y_r = (Current%Part_Pos(2)-y_start_local)/dy -0.5_num
          cell_y = NINT(cell_y_r)
          cell_frac_y = REAL(cell_y,num) - cell_y_r
          cell_y=cell_y+1
          IF (loadlist(cell_x,cell_y)) THEN
             EXIT
          ENDIF
       ENDDO
       Current=>Current%Next
    ENDDO

    DO WHILE(ASSOCIATED(Current))
       Next=>Current%Next
       CALL Remove_Particle_From_PartList(PartList,Current)
       DEALLOCATE(Current)
       Current=>Next
    ENDDO
    CALL MPI_REDUCE(PartList%Count,npart_this_species,1,MPI_INTEGER8,MPI_SUM,0,comm,errcode)
    SpeciesList%Count=npart_this_species
    IF (rank .EQ. 0) THEN
       PRINT *,"Loaded",npart_this_species,"particles of species ",TRIM(SpeciesList%Name)
       WRITE(20,*) "Loaded",npart_this_species,"particles of species ",TRIM(SpeciesList%Name)
    ENDIF
    CALL Particle_BCS


  END SUBROUTINE LoadParticles

!!$  !Subroutine to initialise a thermal particle distribution
  SUBROUTINE SetupParticleTemperature(Temperature,Direction,PartFamily,drift,idum)

    REAL(num),DIMENSION(-2:,-2:), INTENT(IN) :: Temperature
    INTEGER, INTENT(IN) :: Direction
    TYPE(ParticleFamily),POINTER,INTENT(INOUT) :: PartFamily
	 REAL(num), DIMENSION(3), INTENT(IN) :: drift
    INTEGER, INTENT(INOUT) :: idum
    TYPE(ParticleList),POINTER :: PartList
    REAL(num) :: mass,temp_local
    REAL(num) :: cell_x_r,cell_frac_x,cell_y_r,cell_frac_y
    REAL(num),DIMENSION(-2:2) :: gx, gy
    TYPE(particle),POINTER :: Current
    INTEGER :: cell_x,cell_y
    INTEGER(KIND=8) :: ipart

    PartList=>PartFamily%AttachedList
    Current=>PartList%Head
    ipart=0
    DO WHILE(ipart < PartList%Count)
#ifdef PER_PARTICLE_CHARGEMASS
       mass=Current%Mass
#else
       mass=PartFamily%Mass
#endif

       !Assume that temperature is cell centred
       cell_x_r = (Current%Part_Pos(1)-x_start_local)/dx - 0.5_num
       cell_x = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       cell_y_r = (Current%Part_Pos(2)-y_start_local)/dy -0.5_num
       cell_y = NINT(cell_y_r)
       cell_frac_y = REAL(cell_y,num) - cell_y_r
       cell_y=cell_y+1

		 CALL GridToParticle(cell_frac_x,gx)
		 CALL GridToParticle(cell_frac_y,gy)

		 temp_local=0.0_num
		 DO ix=-sf_order,sf_order
			DO iy=-sf_order,sf_order
				temp_local=temp_local+gx(ix)*gy(iy)*Temperature(cell_x+ix,cell_y+iy)
			ENDDO
 		 ENDDO

       IF (IAND(Direction,DIR_X) .NE. 0) Current%Part_P(1)&
			=MomentumFromTemperature(mass,temp_local,idum) + drift(1)

       IF (IAND(Direction,DIR_Y) .NE. 0) Current%Part_P(2)&
			=MomentumFromTemperature(mass,temp_local,idum) + drift(2)

       IF (IAND(Direction,DIR_Z) .NE. 0) Current%Part_P(3)&
			=MomentumFromTemperature(mass,temp_local,idum) + drift(3)

       Current=>Current%Next
       ipart=ipart+1
    ENDDO

  END SUBROUTINE SetupParticleTemperature


  SUBROUTINE SetupParticleDensity(DensityIn,PartFamily,min_density,max_density,idum)

    REAL(num),DIMENSION(-2:,-2:), INTENT(IN) :: DensityIn
    TYPE(ParticleFamily),POINTER :: PartFamily
    REAL(num),INTENT(IN) :: min_density,max_density
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: weight_local
    REAL(num) :: cell_x_r,cell_frac_x
    REAL(num) :: cell_y_r,cell_frac_y
    TYPE(particle),POINTER :: Current
    INTEGER :: cell_x,cell_y
    INTEGER(KIND=8) :: ipart
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Weight_Fn,Temp
    REAL(num),DIMENSION(-2:2) :: gx
    REAL(num),DIMENSION(-2:2) :: gy
    REAL(num) :: Data
    TYPE(ParticleList),POINTER :: PartList
    INTEGER :: iSubx,iSuby
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Density
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: DensityMap

    ALLOCATE(Density(-2:nx+3,-2:ny+3),DensityMap(-2:nx+3,-2:ny+3))
    Density=0.0_num
    Density=DensityIn
    CALL Field_BC(Density)

    DensityMap=.FALSE.
    DO iy=-2,ny+3
       DO ix=-2,nx+3
          IF (Density(ix,iy) .GT. min_density) THEN
             DensityMap(ix,iy)=.TRUE.
          ENDIF
          IF (Density(ix,iy) .GT. max_density  .AND. max_density .GT. 0.0_num) THEN
             Density(ix,iy)=max_density
          ENDIF
       ENDDO
    ENDDO

#ifdef PER_PARTICLE_WEIGHT
    !Uniformly load particles in space
    CALL LoadParticles(PartFamily,DensityMap,idum)
    DEALLOCATE(DensityMap)

    ALLOCATE(Weight_Fn(-2:nx+3,-2:ny+3),Temp(-2:nx+3,-2:ny+3))
    CALL MPI_BARRIER(comm,errcode)
    Weight_Fn=0.0_num
    Temp=0.0_num

    PartList=>PartFamily%AttachedList
    !If using per particle weighing then use the weight function to match the uniform pseudoparticle density to the 
    !Real particle density
    Current=>PartList%Head
    ipart=0
    !First loop converts number density into weight function
    DO WHILE(ipart < PartList%Count)
       IF (.NOT. ASSOCIATED(Current)) PRINT *,"Bad Particle"
       cell_x_r = (Current%Part_Pos(1)-x_start_local) / dx - 0.5_num
       cell_x=NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       cell_y_r = (Current%Part_Pos(2)-y_start_local) / dy - 0.5_num
       cell_y=NINT(cell_y_r)
       cell_frac_y = REAL(cell_y,num) - cell_y_r
       cell_y=cell_y+1

		 CALL ParticleToGrid(cell_frac_x,gx)
		 CALL ParticleToGrid(cell_frac_y,gy)

       Data=1.0_num/(dx*dy) !Simply want to count particles per metre^2
       DO iSuby=-sf_order,sf_order
          DO iSubx=-sf_order,sf_order
             Weight_Fn(cell_x+iSubx,cell_y+iSuby) = Weight_Fn(cell_x+iSubx,cell_y+iSuby) + & 
                  gx(iSubx) * gy(iSuby) * Data
          ENDDO
       ENDDO

       Current=>Current%Next
       ipart=ipart+1
    ENDDO
    CALL Processor_Summation_BCS(Weight_Fn)
    IF (left  .EQ. MPI_PROC_NULL) Weight_Fn(0 ,:)=Weight_Fn(1,:)
    IF (right .EQ. MPI_PROC_NULL) Weight_Fn(nx,:)=Weight_Fn(nx-1,:)
    IF (down  .EQ. MPI_PROC_NULL) Weight_Fn(: ,0)=Weight_Fn(:,1)
    IF (up    .EQ. MPI_PROC_NULL) Weight_Fn(:,ny)=Weight_Fn(:,ny-1)
    CALL Field_Zero_Gradient(Weight_Fn,.TRUE.)
    DO iy=-2,ny+2
       DO ix=-2,nx+2
          IF (Weight_Fn(ix,iy) .GT. 0.0_num) THEN
             Weight_Fn(ix,iy)=Density(ix,iy)/Weight_Fn(ix,iy)
          ELSE
             Weight_Fn(ix,iy)=0.0_num
          ENDIF
       ENDDO
    ENDDO
    IF (left  .EQ. MPI_PROC_NULL) Weight_Fn(0 ,:)=Weight_Fn(1,:)
    IF (right .EQ. MPI_PROC_NULL) Weight_Fn(nx,:)=Weight_Fn(nx-1,:)
    IF (down  .EQ. MPI_PROC_NULL) Weight_Fn(: ,0)=Weight_Fn(:,1)
    IF (up    .EQ. MPI_PROC_NULL) Weight_Fn(:,ny)=Weight_Fn(:,ny-1)
    CALL Field_Zero_Gradient(Weight_Fn,.TRUE.)


    PartList=>PartFamily%AttachedList
    !Second loop actually assigns weights to particles
    !Again assumes linear interpolation
    Current=>PartList%Head
    ipart=0
    DO WHILE(ipart < PartList%Count)
       cell_x_r = (Current%Part_Pos(1)-x_start_local) / dx !-0.5_num
       cell_x=NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       cell_y_r = (Current%Part_Pos(2)-y_start_local) / dy !-0.5_num
       cell_y=NINT(cell_y_r)
       cell_frac_y = REAL(cell_y,num) - cell_y_r
       cell_y=cell_y+1

		 CALL GridToParticle(cell_frac_x,gx)
		 CALL GridToParticle(cell_frac_y,gy)

       weight_local=0.0_num
       DO iSuby=-sf_order,sf_order
          DO iSubx=-sf_order,+sf_order
             weight_local=weight_local+gx(iSubx)*gy(iSuby)*Weight_Fn(cell_x+iSubx,cell_y+iSuby)
          ENDDO
       ENDDO
       Current%Weight=weight_local
       Current=>Current%Next
       ipart=ipart+1
    ENDDO
#endif
    DEALLOCATE(Weight_Fn)
    DEALLOCATE(Density)
  END SUBROUTINE SetupParticleDensity


  FUNCTION MomentumFromTemperature(mass,temperature,idum)

    REAL(num), INTENT(IN) :: mass,temperature
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: MomentumFromTemperature

    REAL(num) :: stdev
    REAL(num) :: rand1,rand2,w

    !This is a basic polar Box-Muller transform
    !It generates gaussian distributed random numbers
    !The standard deviation (stdev) is related to temperature

    stdev=SQRT(2.0_num * temperature*kb*mass)

    DO
       rand1=Random(idum)
       rand2=Random(idum)

       rand1=2.0_num*rand1 - 1.0_num
       rand2=2.0_num*rand2 - 1.0_num

       w=rand1**2 + rand2**2

       IF (w .LT. 1.0_num) EXIT
    ENDDO

    w = SQRT((-2.0_num * LOG(w) )/w)

    MomentumFromTemperature = rand1 * w * stdev

  END FUNCTION MomentumFromTemperature

  FUNCTION Sample_Dist_Function(axis,Dist_Fn,idum)
    REAL(num),DIMENSION(:),INTENT(IN) :: axis,Dist_Fn
    INTEGER, INTENT(INOUT) :: idum
    REAL(num),DIMENSION(:),ALLOCATABLE :: CDF
    REAL(num) :: Position,d_cdf
    INTEGER :: n_points, iPoint, start, endpoint, current
    REAL(num) :: Sample_Dist_Function

    n_points=SIZE(Dist_Fn)
    ALLOCATE(CDF(1:n_points))
    DO iPoint=1,n_points
       CDF(iPoint)=SUM(Dist_Fn(1:iPoint))
    ENDDO
    CDF=CDF/SUM(Dist_Fn)

    Position=Random(iDum)
	sample_dist_function=0.0_num

    start=1
    endpoint=n_points
    current=(start+endpoint)/2

    DO current=1,n_points-1
       IF (CDF(Current) .LE. Position .AND. CDF(Current+1) .GE. Position) THEN
          d_cdf=CDF(Current+1)-CDF(Current)
          Sample_Dist_Function=(Axis(Current)*(Position-CDF(Current))/d_cdf + &
               Axis(Current+1)*(CDF(Current+1)-Position)/d_cdf)
          EXIT
       ENDIF
    ENDDO
    DEALLOCATE(CDF)


  END FUNCTION Sample_Dist_Function

  FUNCTION Random(idum)

    INTEGER,INTENT(INOUT) :: idum
    REAL(num) :: Random
    INTEGER,PARAMETER :: IA=16807, IM=2147483647, IQ=127773
    INTEGER,PARAMETER :: IR=2836, MASK=123459876
    REAL(dbl),PARAMETER :: AM=1.0_dbl/2147483647.0_dbl

    INTEGER :: k

    idum=XOR(idum,mask)
    k=idum/IQ

    idum=IA*(idum-k*IQ)-IR*k
    IF (idum .LT. 0) THEN 
       idum=idum+IM
    ENDIF

    Random=AM*idum
    idum=XOR(idum,mask)

    IF (random .GT. max_rand) max_rand=random

  END FUNCTION Random

END MODULE helper
