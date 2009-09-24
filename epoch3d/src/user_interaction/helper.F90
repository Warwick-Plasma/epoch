MODULE helper

  USE shared_data
  USE partlist
  USE boundary
  USE strings
  IMPLICIT NONE

  SAVE
  REAL(num) :: max_rand

CONTAINS

  SUBROUTINE AutoLoad
    INTEGER :: iSpecies
    INTEGER :: clock,idum
    TYPE(ParticleFamily),POINTER :: PartFam
	REAL(num) :: max_rho_local,max_rho
	
	dt_plasma=1000000.0_num

    CALL SYSTEM_CLOCK(clock)
    idum=-(clock+rank+1)
    DO iSpecies=1,nspecies
       PartFam=>ParticleSpecies(iSpecies)
		 !Calculate the inverse plasma frequency for stability
		 max_rho_local=MAXVAL(InitialConditions(iSpecies)%Rho)
		 CALL MPI_ALLREDUCE(max_rho_local,max_rho,1,mpireal,MPI_MAX,comm,errcode)
		 dt_plasma=MIN(dt_plasma,1.0_num/SQRT(max_rho*Q0**2/(M0*EPSILON0)))
       !ParticleSpecies(iSpecies)%Window_Density=InitialConditions(iSpecies)%Rho(nx,:,:)
       !ParticleSpecies(iSpecies)%Window_Temperature=InitialConditions(iSpecies)%Temp(nx,:,:,:)

       CALL SetupParticleDensity(InitialConditions(iSpecies)%Rho,Partfam,InitialConditions(iSpecies)%MinRho,InitialConditions(iSpecies)%MaxRho,idum)
       CALL SetupParticleTemperature(InitialConditions(iSpecies)%Temp(:,:,:,1),DIR_X,PartFam,idum)
       CALL SetupParticleTemperature(InitialConditions(iSpecies)%Temp(:,:,:,2),DIR_Y,PartFam,idum)
       CALL SetupParticleTemperature(InitialConditions(iSpecies)%Temp(:,:,:,3),DIR_Z,PartFam,idum)
    ENDDO

	IF (rank .EQ. 0) PRINT *,"Plasma characteristic timescale ",dt_plasma
  END SUBROUTINE AutoLoad

  SUBROUTINE AllocateIC
    INTEGER :: iSpecies

    ALLOCATE(InitialConditions(1:nspecies))
    DO iSpecies=1,nSpecies
       ALLOCATE(InitialConditions(iSpecies)%Rho(-2:nx+3,-2:ny+3,-2:nz+3))
       ALLOCATE(InitialConditions(iSpecies)%Temp(-2:nx+3,-2:nx+3,-2:nz+3,1:3))

       InitialConditions(iSpecies)%Rho=0.0_num
       InitialConditions(iSpecies)%Temp=0.0_num
       InitialConditions(iSpecies)%MinRho=0.0_num
       InitialConditions(iSpecies)%MaxRho=0.0_num
    ENDDO
    Ex=0.0_num
    Ey=0.0_num
    Ez=0.0_num

    Bx=0.0_num
    By=0.0_num
    Bz=0.0_num
  END SUBROUTINE AllocateIC

  SUBROUTINE DeallocateIC
    INTEGER :: iSpecies

    DO iSpecies=1,nspecies
       DEALLOCATE(InitialConditions(iSpecies)%Rho)
       DEALLOCATE(InitialConditions(iSpecies)%Temp)
    ENDDO
    DEALLOCATE(InitialConditions)
  END SUBROUTINE DeallocateIC

  !This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE LoadParticles(SpeciesList,loadlist,idum)

    INTEGER, INTENT(INOUT) :: idum
    TYPE(ParticleFamily), POINTER,INTENT(INOUT) :: SpeciesList
    LOGICAL,DIMENSION(-2:,-2:,-2:),INTENT(IN) :: LoadList
    TYPE(ParticleList),POINTER :: PartList
    TYPE(particle),POINTER :: Current,PrevHead,Next
    INTEGER :: curspecies
    INTEGER(KIND=8) :: ipart,npart_per_cell,npart,ipart_total,ncell_per_part,iproc
    INTEGER(KIND=8) :: num_valid_cells,temp,num_valid_cells_local,npart_this_species
    INTEGER(KIND=8) :: num_new_particles,npart_left
    REAL(num) :: valid_cell_frac
    REAL(dbl) :: frac, frac_to_pos,rpos
    LOGICAL :: fullrand
    INTEGER :: cell_x,cell_y,cell_z
    REAL(num) :: cell_x_r,dcell_x,cell_frac_x
    REAL(num) :: cell_y_r,dcell_y,cell_frac_y
    REAL(num) :: cell_z_r,dcell_z,cell_frac_z
    INTEGER(KIND=8) :: i
    INTEGER :: j
    CHARACTER(LEN=15) :: string

    PartList=>SpeciesList%AttachedList

    npart_this_species=SpeciesList%Count
    IF (npart_this_species .LT. 0) THEN
       IF (rank .EQ. 0) PRINT *,"Unable to continue, species ",TRIM(SpeciesList%Name)," has not had a number of particles set"
       CALL MPI_ABORT(comm,errcode)
    ENDIF
    IF (npart_this_species .EQ. 0) RETURN
    num_valid_cells_local=0
    DO iz=1,nz
       DO iy=1,ny
          DO ix=1,nx
             IF (loadlist(ix,iy,iz)) num_valid_cells_local=num_valid_cells_local+1
          ENDDO
       ENDDO
    ENDDO

    !Calculate global number of particles per cell
    CALL MPI_ALLREDUCE(num_valid_cells_local,num_valid_cells,1,MPI_INTEGER8,MPI_SUM,comm,errcode)

    valid_cell_frac=REAL(num_valid_cells_local,num)/REAL(num_valid_cells,num)
    num_new_particles=INT(npart_this_species*valid_cell_frac,KIND=8)
    CALL Destroy_PartList(PartList)
    CALL Create_Allocated_PartList(PartList,num_new_particles)


    npart_per_cell=npart_this_species/num_valid_cells
    ncell_per_part=num_valid_cells/npart_this_species
    SpeciesList%window_npart_per_cell=npart_per_cell
    IF (SpeciesList%window_npart_per_cell .EQ. 0) SpeciesList%window_npart_per_cell=1

    ipart=0
    ipart_total=0
    ix=1
    iy=1

    npart_left=npart_this_species
    Current=>PartList%Head
    IF (npart_per_cell .GT. 0) THEN

       DO ix=1,nx
          DO iy=1,ny
             DO iz=1,nz
                ipart=0
                IF (loadlist(ix,iy,iz)) THEN
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
                      rpos=random(idum)-0.5_num
                      rpos=(rpos*dz)+z(iz)
                      Current%Part_Pos(3)=rpos
                      ipart=ipart+1
                      Current=>Current%Next
                      !One particle sucessfully placed
                      npart_left=npart_left-1
                   ENDDO
                ENDIF
             ENDDO
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
          rpos=random(idum)*(z(nz)-z(1) + dz) - dz/2.0_num
          Current%Part_Pos(3)=rpos+z(1)
          cell_x_r = (Current%Part_Pos(1)-x_start_local)/dx - 0.5_num
          cell_x = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x,num) - cell_x_r
          cell_x=cell_x+1

          cell_y_r = (Current%Part_Pos(2)-y_start_local)/dy -0.5_num
          cell_y = NINT(cell_y_r)
          cell_frac_y = REAL(cell_y,num) - cell_y_r
          cell_y=cell_y+1

          cell_z_r = (Current%Part_Pos(3)-z_start_local)/dz -0.5_num
          cell_z = NINT(cell_z_r)
          cell_frac_z = REAL(cell_z,num) - cell_z_r
          cell_z=cell_z+1

          IF (loadlist(cell_x,cell_y,cell_z)) THEN
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
       CALL Integer8AsString(npart_this_species,string)
       PRINT *,"Loaded ",TRIM(ADJUSTL(string))," particles of species ",TRIM(SpeciesList%Name)
       WRITE(20,*) "Loaded ",TRIM(ADJUSTL(string))," particles of species ",TRIM(SpeciesList%Name)
    ENDIF


  END SUBROUTINE LoadParticles

!!$
!!$  !Subroutine to initialise a thermal particle distribution
!!$  !Assumes linear interpolation of temperature between cells
  SUBROUTINE SetupParticleTemperature(Temperature,Direction,PartFamily,idum)

    REAL(num),DIMENSION(-2:,-2:,-2:), INTENT(IN) :: Temperature
    INTEGER, INTENT(IN) :: Direction
    TYPE(ParticleFamily),POINTER,INTENT(INOUT) :: PartFamily
    INTEGER, INTENT(INOUT) :: idum
    TYPE(ParticleList),POINTER :: PartList
    REAL(num) :: mass,temp_local
    REAL(num) :: cell_x_r,cell_frac_x
    REAL(num) :: cell_y_r,cell_frac_y
    REAL(num) :: cell_z_r,cell_frac_z
    REAL(num) :: part_weight
    REAL(num) :: g0x,gpx,gmx
    REAL(num) :: g0y,gpy,gmy
    REAL(num) :: g0z,gpz,gmz
    TYPE(particle),POINTER :: Current
    INTEGER :: cell_x,cell_y,cell_z
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

       cell_z_r = (Current%Part_Pos(3)-z_start_local)/dz -0.5_num
       cell_z = NINT(cell_z_r)
       cell_frac_z = REAL(cell_z,num) - cell_z_r
       cell_z=cell_z+1

       gmx=0.5_num * (0.5_num + cell_frac_x)**2
       g0x=0.75_num - cell_frac_x**2
       gpx=0.5_num * (0.5_num - cell_frac_x)**2

       gmy=0.5_num * (0.5_num + cell_frac_y)**2
       g0y=0.75_num - cell_frac_y**2
       gpy=0.5_num * (0.5_num - cell_frac_y)**2

       gmz=0.5_num * (0.5_num + cell_frac_z)**2
       g0z=0.75_num - cell_frac_z**2
       gpz=0.5_num * (0.5_num - cell_frac_z)**2


       temp_local = gmz * (&
            gmy * (gmx * Temperature(cell_x-1,cell_y-1,cell_z-1) + g0x * &
            Temperature(cell_x,cell_y-1,cell_z-1) + gpx * Temperature(cell_x+1,cell_y-1,cell_z-1)) + &
            g0y * (gmx * Temperature(cell_x-1,cell_y  ,cell_z-1) + g0x * &
            Temperature(cell_x,cell_y  ,cell_z-1) + gpx * Temperature(cell_x+1,cell_y  ,cell_z-1)) + &
            gpy * (gmx * Temperature(cell_x-1,cell_y+1,cell_z-1) + g0x * &
            Temperature(cell_x,cell_y+1,cell_z-1) + gpx * Temperature(cell_x+1,cell_y+1,cell_z-1))) + &
            g0z *(&
            gmy * (gmx * Temperature(cell_x-1,cell_y-1,cell_z  ) + g0x * &
            Temperature(cell_x,cell_y-1,cell_z  ) + gpx * Temperature(cell_x+1,cell_y-1,cell_z  )) + &
            g0y * (gmx * Temperature(cell_x-1,cell_y  ,cell_z  ) + g0x * &
            Temperature(cell_x,cell_y  ,cell_z  ) + gpx * Temperature(cell_x+1,cell_y  ,cell_z  )) + &
            gpy * (gmx * Temperature(cell_x-1,cell_y+1,cell_z  ) + g0x * &
            Temperature(cell_x,cell_y+1,cell_z  ) + gpx * Temperature(cell_x+1,cell_y+1,cell_z  ))) + &
            gpz *(&
            gmy * (gmx * Temperature(cell_x-1,cell_y-1,cell_z+1) + g0x * &
            Temperature(cell_x,cell_y-1,cell_z+1) + gpx * Temperature(cell_x+1,cell_y-1,cell_z+1)) + &
            g0y * (gmx * Temperature(cell_x-1,cell_y  ,cell_z+1) + g0x * &
            Temperature(cell_x,cell_y  ,cell_z+1) + gpx * Temperature(cell_x+1,cell_y  ,cell_z+1)) + &
            gpy * (gmx * Temperature(cell_x-1,cell_y+1,cell_z+1) + g0x * &
            Temperature(cell_x,cell_y+1,cell_z+1) + gpx * Temperature(cell_x+1,cell_y+1,cell_z+1)))

       IF (IAND(Direction,DIR_X) .NE. 0) Current%Part_P(1)=MomentumFromTemperature(mass,temp_local,idum)

       IF (IAND(Direction,DIR_Y) .NE. 0) Current%Part_P(2)=MomentumFromTemperature(mass,temp_local,idum)

       IF (IAND(Direction,DIR_Z) .NE. 0) Current%Part_P(3)=MomentumFromTemperature(mass,temp_local,idum)
       Current=>Current%Next
       ipart=ipart+1
    ENDDO

  END SUBROUTINE SetupParticleTemperature


  SUBROUTINE SetupParticleDensity(DensityIn,PartFamily,min_density,max_density,idum)

    REAL(num),DIMENSION(-2:,-2:,-2:), INTENT(IN) :: DensityIn
    TYPE(ParticleFamily),POINTER :: PartFamily
    REAL(num),INTENT(IN) :: min_density,max_density
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: weight_local
    REAL(num) :: cell_x_r,dcell_x,cell_frac_x
    REAL(num) :: cell_y_r,dcell_y,cell_frac_y
    REAL(num) :: cell_z_r,dcell_z,cell_frac_z
    TYPE(particle),POINTER :: Current
    INTEGER :: cell_x,cell_y,cell_z
    INTEGER(KIND=8) :: ipart
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: Weight_Fn,Temp
    REAL(num),DIMENSION(-1:1) :: gx,gy,gz
    REAL(num) :: Data,rpos
    TYPE(ParticleList),POINTER :: PartList
    INTEGER :: iSubx,iSuby,iSubz
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: Density
    LOGICAL,DIMENSION(:,:,:),ALLOCATABLE :: DensityMap

    ALLOCATE(Density(-2:nx+3,-2:ny+3,-2:nz+3),DensityMap(-2:nx+3,-2:ny+3,-2:nz+3))
    Density=DensityIn

    CALL Field_BC(Density)

    DensityMap=.FALSE.
    DO iz=-2,nz+3
       DO iy=-2,ny+3
          DO ix=-2,nx+3
             IF (Density(ix,iy,iz) .GT. min_density) THEN
                DensityMap(ix,iy,iz)=.TRUE.
             ENDIF
             IF (Density(ix,iy,iz) .GT. max_density  .AND. max_density .GT. 0.0_num) THEN
                Density(ix,iy,iz)=max_density
             ENDIF
          ENDDO
       ENDDO
    ENDDO

#ifdef PER_PARTICLE_WEIGHT
    !Uniformly load particles in space
    CALL LoadParticles(PartFamily,DensityMap,idum)
    DEALLOCATE(DensityMap)

    ALLOCATE(Weight_Fn(-2:nx+3,-2:ny+3,-2:nz+3),Temp(-2:nx+3,-2:ny+3,-2:nz+3))
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
       cell_x_r = (Current%Part_Pos(1)-x_start_local) / dx !- 0.5_num
       cell_x=NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       cell_y_r = (Current%Part_Pos(2)-y_start_local) / dy !- 0.5_num
       cell_y=NINT(cell_y_r)
       cell_frac_y = REAL(cell_y,num) - cell_y_r
       cell_y=cell_y+1

       cell_z_r = (Current%Part_Pos(3)-z_start_local) / dz !- 0.5_num
       cell_z=NINT(cell_z_r)
       cell_frac_z = REAL(cell_z,num) - cell_z_r
       cell_z=cell_z+1
!       PRINT *,rank,cell_z

       gx(-1) = 0.5_num * (1.5_num - ABS(cell_frac_x - 1.0_num))**2
       gx( 0) = 0.75_num - ABS(cell_frac_x)**2
       gx( 1) = 0.5_num * (1.5_num - ABS(cell_frac_x + 1.0_num))**2

       gy(-1) = 0.5_num * (1.5_num - ABS(cell_frac_y - 1.0_num))**2
       gy( 0) = 0.75_num - ABS(cell_frac_y)**2
       gy( 1) = 0.5_num * (1.5_num - ABS(cell_frac_y + 1.0_num))**2

       gz(-1) = 0.5_num * (1.5_num - ABS(cell_frac_z - 1.0_num))**2
       gz( 0) = 0.75_num - ABS(cell_frac_z)**2
       gz( 1) = 0.5_num * (1.5_num - ABS(cell_frac_z + 1.0_num))**2
       Data=1.0_num/(dx*dy*dz) !Simply want to count particles per metre^2
       DO iSubz=-1,1
          DO iSuby=-1,1
             DO iSubx=-1,1
                Weight_Fn(cell_x+iSubx,cell_y+iSuby,cell_z+iSubz) = Weight_Fn(cell_x+iSubx,cell_y+iSuby,cell_z+iSubZ) + & 
                     gx(iSubx) * gy(iSuby) * gz(iSubz) * Data
             ENDDO
          ENDDO
       ENDDO
       Current=>Current%Next
       ipart=ipart+1
    ENDDO
    CALL Processor_Summation_BCS(Weight_Fn)
!!$    IF (left  .EQ. MPI_PROC_NULL) Weight_Fn(0 ,:)=Weight_Fn(1,:)
!!$    IF (right .EQ. MPI_PROC_NULL) Weight_Fn(nx,:)=Weight_Fn(nx-1,:)
!!$    IF (down  .EQ. MPI_PROC_NULL) Weight_Fn(: ,0)=Weight_Fn(:,1)
!!$    IF (up    .EQ. MPI_PROC_NULL) Weight_Fn(:,ny)=Weight_Fn(:,ny-1)
    CALL Field_Zero_Gradient(Weight_Fn,.TRUE.)
    DO iz=-2,nz+2
       DO iy=-2,ny+2
          DO ix=-2,nx+2
             IF (Weight_Fn(ix,iy,iz) .GT. 0.0_num) THEN
                Weight_Fn(ix,iy,iz)=Density(ix,iy,iz)/Weight_Fn(ix,iy,iz)
             ELSE
                Weight_Fn(ix,iy,iz)=0.0_num
             ENDIF
          ENDDO
       ENDDO
    ENDDO
!!$    IF (left  .EQ. MPI_PROC_NULL) Weight_Fn(0 ,:)=Weight_Fn(1,:)
!!$    IF (right .EQ. MPI_PROC_NULL) Weight_Fn(nx,:)=Weight_Fn(nx-1,:)
!!$    IF (down  .EQ. MPI_PROC_NULL) Weight_Fn(: ,0)=Weight_Fn(:,1)
!!$    IF (up    .EQ. MPI_PROC_NULL) Weight_Fn(:,ny)=Weight_Fn(:,ny-1)
    CALL Field_Zero_Gradient(Weight_Fn,.TRUE.)


    PartList=>PartFamily%AttachedList
    !Second loop actually assigns weights to particles
    !Again assumes linear interpolation
    Current=>PartList%Head
    ipart=0
    DO WHILE(ipart < PartList%Count)
       cell_x_r = (Current%Part_Pos(1)-x_start_local) / dx -0.5_num
       cell_x=NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       cell_y_r = (Current%Part_Pos(2)-y_start_local) / dy -0.5_num
       cell_y=NINT(cell_y_r)
       cell_frac_y = REAL(cell_y,num) - cell_y_r
       cell_y=cell_y+1

       cell_z_r = (Current%Part_Pos(3)-z_start_local) / dz -0.5_num
       cell_z=NINT(cell_z_r)
       cell_frac_z = REAL(cell_z,num) - cell_z_r
       cell_z=cell_z+1


       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

       gy(-1) = 0.5_num * (0.5_num + cell_frac_y)**2
       gy( 0) = 0.75_num - cell_frac_y**2
       gy( 1) = 0.5_num * (0.5_num - cell_frac_y)**2

       gz(-1) = 0.5_num * (0.5_num + cell_frac_z)**2
       gz( 0) = 0.75_num - cell_frac_z**2
       gz( 1) = 0.5_num * (0.5_num - cell_frac_z)**2

       weight_local=0.0_num
       DO iSubz=-1,+1
          DO iSuby=-1,+1
             DO iSubx=-1,+1
                weight_local=weight_local+gx(iSubx)*gy(iSuby)*gz(iSubz)*Weight_Fn(cell_x+iSubx,cell_y+iSuby,cell_z+iSubz)
             ENDDO
          ENDDO
       ENDDO
       Current%Weight=weight_local
       Current=>Current%Next
       ipart=ipart+1
    ENDDO
#else
    IF (rank .EQ. 0) PRINT *,"Autoloader only available when using per particle weighting"
    CALL MPI_ABORT(comm,errcode)
#endif
    DEALLOCATE(Weight_Fn)
    DEALLOCATE(Density)
  END SUBROUTINE SetupParticleDensity


  FUNCTION MomentumFromTemperature(mass,temperature,idum)

    REAL(num), INTENT(IN) :: mass,temperature
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: MomentumFromTemperature

    REAL(num) :: mean=0.0_num
    REAL(num) :: stdev
    REAL(num) :: rand1,rand2,w

    !This is a basic polar Box-Muller transform
    !It generates gaussian distributed random numbers
    !The standard deviation (stdev) is related to temperature

    stdev=SQRT(temperature*kb*mass)

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
