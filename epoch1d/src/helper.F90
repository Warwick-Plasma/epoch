MODULE helper

  USE shared_data
  USE partlist
  IMPLICIT NONE

  SAVE
  !Direction parameters
  INTEGER, PARAMETER :: DIR_X=1, DIR_Y=2, DIR_Z=4, DIR_UX=8, DIR_UY=16, DIR_UZ=32
  REAL(num) :: max_rand

CONTAINS

  !This subroutine automatically loads a uniform density of pseudoparticles
  !Uses the "fraction" property of the species information to assign particle type
  SUBROUTINE LoadParticles(idum,SpeciesLists)

    INTEGER, INTENT(INOUT) :: idum
    TYPE(ParticlePointer), DIMENSION(1:nspecies), INTENT(INOUT) :: SpeciesLists
    TYPE(particle),POINTER :: Current,PrevHead
    INTEGER(KIND=8) :: ipart
    INTEGER:: curspecies
    REAL(dbl) :: frac, frac_to_pos,rpos

    DO ipart=1,nspecies
       CALL Create_Empty_PartList(SpeciesLists(ipart))
    ENDDO

    Current=>MainRoot%Head
    ipart=0
    frac_to_pos=0.0_num
    curspecies=1

    PrevHead=>Current
    DO WHILE(ASSOCIATED(Current))

       !Move to next particle species when this one exhausted
       frac=REAL(ipart,dbl)/REAL(MainRoot%Count,dbl)
       IF (frac .GE. frac_to_pos + species(curspecies,3)) THEN
          !Current will then be the FIRST instance of curspecies+1
          CALL Create_Unsafe_PartList_By_Tail(SpeciesLists(curspecies),PrevHead,Current%Prev)
          frac_to_pos=frac_to_pos+species(curspecies,3)
          curspecies=curspecies+1
          PrevHead=>Current%Next
       ENDIF

       !Random positions
       rpos=x_start-dx
       DO WHILE(rpos .LT. x_start .OR. rpos .GT. x_end)
          rpos=random(idum)
          rpos=(rpos*(x_end-x_start))+x_start
       ENDDO
       Current%Part_pos=rpos

       Current%Part_Species=curspecies
#ifdef PER_PARTICLE_CHARGEMASS
       Current%Charge=Species(Current%Part_Species,1)
       Current%Mass=Species(Current%Part_Species,2)
#endif

       Current%Part_P=0.0_num

       ipart=ipart+1
       Current=>Current%Next
    ENDDO
    CALL Create_Unsafe_PartList_By_Tail(SpeciesLists(nspecies),PrevHead,MainRoot%Tail)

  END SUBROUTINE LoadParticles

  !Wrapper routine that applies the same temperature to all particles
  SUBROUTINE SetupParticlesByTemperature(Temperature,Direction,idum)

    REAL(num),DIMENSION(0:), INTENT(IN) :: Temperature
    INTEGER, INTENT(IN) :: Direction
    INTEGER, INTENT(INOUT) :: idum

    CALL SetupParticleListByTemperature(Temperature,Direction,MainRoot,idum)

  END SUBROUTINE SetupParticlesByTemperature

  !Subroutine to initialise a thermal particle distribution
  !Assumes linear interpolation of temperature between cells
  SUBROUTINE SetupParticleListByTemperature(Temperature,Direction,PartList,idum)

    REAL(num),DIMENSION(-1:), INTENT(IN) :: Temperature
    INTEGER, INTENT(IN) :: Direction
    TYPE(ParticlePointer),INTENT(IN) :: PartList
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: mass,temp_local,cell_x_r,cell_frac_x,part_weight
    REAL(num) :: g0x,gpx,gmx
    TYPE(particle),POINTER :: Current
    INTEGER :: cell_x,cell_xp
    INTEGER(KIND=8) :: ipart

    Current=>PartList%Head
    ipart=0
    DO WHILE(ipart < PartList%Count)
#ifdef PER_PARTICLE_CHARGEMASS
       mass=Current%Mass
#else
       mass=species(Current%part_species,2)
#endif

       !Particle weighting function
#ifdef PER_PARTICLE_WEIGHT
       part_weight=Current%Weight
#else
       part_weight=weight
#endif


       !Assume that temperature is cell centred
       cell_x_r = (Current%Part_Pos-x_start)/dx
       cell_x=NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       gmx=0.5_num * (0.5_num + cell_frac_x)**2
       g0x=0.75_num - cell_frac_x**2
       gpx=0.5_num * (0.5_num - cell_frac_x)**2

       temp_local=gmx*Temperature(cell_x-1) + g0x*Temperature(cell_x) + gpx*Temperature(cell_x+1)

       IF (IAND(Direction,DIR_X) .NE. 0) Current%Part_P(1)=MomentumFromTemperature(mass,temp_local,idum)

       IF (IAND(Direction,DIR_Y) .NE. 0) Current%Part_P(2)=MomentumFromTemperature(mass,temp_local,idum)

       IF (IAND(Direction,DIR_Z) .NE. 0) Current%Part_P(3)=MomentumFromTemperature(mass,temp_local,idum)
       Current=>Current%Next
       ipart=ipart+1
    ENDDO

  END SUBROUTINE SetupParticleListByTemperature

  SUBROUTINE SetupParticlesByDensity(Density,idum)

    REAL(num),DIMENSION(-1:),INTENT(IN) :: Density
    INTEGER,INTENT(INOUT) :: idum

    CALL SetupParticleListByDensity(Density,MainRoot,idum)

  END SUBROUTINE SetupParticlesByDensity

  SUBROUTINE SetupParticleListByDensity(Density,PartList,idum)

    REAL(num),DIMENSION(-1:), INTENT(IN) :: Density
    TYPE(ParticlePointer),INTENT(IN) :: PartList
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: weight_local,cell_x_r,dcell_x,cell_frac_x
    TYPE(particle),POINTER :: Current
    INTEGER :: cell_x,cell_xp
    INTEGER(KIND=8) :: ipart
    REAL(num),DIMENSION(:),ALLOCATABLE :: Weight_Fn,Temp
    REAL(num),DIMENSION(-1:1) :: gx
    REAL(num) :: Data,rpos
    INTEGER :: Low, High, Point, OldPoint

#ifdef PER_PARTICLE_WEIGHT
    !If using per particle weighing then use the weight function to match the uniform pseudoparticle density to the 
    !Real particle density
    Current=>PartList%Head
    ipart=0
    !First loop converts number density into weight function
    ALLOCATE(Weight_Fn(-1:nx_global+2))
    Weight_Fn=0.0_num
    DO WHILE(ipart < PartList%Count)
       cell_x_r = (Current%Part_Pos-x_start) / dx
       cell_x  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

       DO ix=-1,1
          Data=1.0_num/dx !Simply want to count particles per metre
          Weight_Fn(cell_x+ix) = Weight_Fn(cell_x+ix) + &
               gx(ix) * Data
       ENDDO
       Current=>Current%Next
       ipart=ipart+1
    ENDDO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,Weight_Fn,nx_global+4,mpireal,MPI_SUM,comm,errcode)

    DO ix=-1,nx_global+1
       IF (Weight_Fn(ix) .GT. 0.0_num) THEN
          Weight_Fn(ix)=Density(ix)/Weight_Fn(ix)
       ENDIF
    ENDDO


    !Second loop actually assigns weights to particles
    !Again assumes linear interpolation
    Current=>PartList%Head
    ipart=0
    DO WHILE(ipart < PartList%Count)
!!$       cell_x_r=(Current%Part_Pos-x_start)/dx
!!$       cell_x = INT(cell_x_r)
!!$       cell_xp = cell_x + 1
!!$
!!$       cell_frac_x=cell_x_r-REAL(cell_x,num)
!!$       weight_local=Weight_Fn(cell_x) * (1.0_num - cell_frac_x) + Weight_Fn(cell_xp) * cell_frac_x

       cell_x_r = (Current%Part_Pos-x_start) / dx
       cell_x  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

       weight_local=0.0_num
       DO ix=-1,+1
          weight_local=weight_local+gx(ix)*Weight_Fn(cell_x+ix)
       ENDDO

       Current%Weight=weight_local

       Current=>Current%Next
       ipart=ipart+1
    ENDDO
    PRINT *,weight_local
#else
    !If not then use a CDF type scheme to redistibute the pseudoparticle density
    !This isn't perfect because particles are still split evenly across processors
    ALLOCATE(Weight_Fn(1:nx_global),temp(-1:nx_global+2))
    DO ix=1,nx_global
       Weight_Fn(ix)=SUM(Density(1:ix)/SUM(Density))
    ENDDO
    Current=>PartList%Head
    ipart=0
    DO WHILE(ipart < PartList%Count)
       !Random number between 0 and 1
       rpos=random(idum)
       rpos=(rpos*(MAXVAL(Weight_Fn)-MINVAL(Weight_Fn)))+MINVAL(Weight_Fn)

       Low=1
       High=nx_global
       Point=(Low+High)/2
       OldPoint=Point
       DO
          IF (rpos .GE. Weight_Fn(Point) .AND. rpos .LE. Weight_Fn(Point+1)) EXIT
          IF (rpos .GT. Weight_Fn(Point+1)) Low=Point+1
          IF (rpos .LT. Weight_Fn(Point)) High=Point
          Point=(Low+High)/2
          IF (Point .EQ. OldPoint) THEN
             PRINT *,"Particle Locked, Point=",Point
          ELSE
             OldPoint=Point
          ENDIF
       ENDDO

       !Uniform random position within the cell
       rpos=random(idum)
       rpos=(rpos*(x_global(Point+1)-x_global(Point)))+x_global(Point)

       Current%Part_Pos=rpos

       ipart=ipart+1
       Current=>Current%Next
    ENDDO
    DEALLOCATE(Weight_Fn)
    !Pseudoparticles are now positioned as best they can
    !Now set weight constant to match to real particle number density
    ALLOCATE(Weight_Fn(-1:nx_global+2))
    Weight_Fn=0.0_num
    Current=>PartList%Head
    ipart=0
    DO WHILE(ipart < PartList%Count)
       cell_x_r = (Current%Part_Pos-x_start) / dx
       cell_x  = NINT(cell_x_r)
       cell_frac_x = REAL(cell_x,num) - cell_x_r
       cell_x=cell_x+1

       gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
       gx( 0) = 0.75_num - cell_frac_x**2
       gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2


       DO ix=-1,1
          Data=1.0_num/dx !Simply want to count particles per metre
          Weight_Fn(cell_x+ix) = Weight_Fn(cell_x+ix) + &
               gx(ix) * Data
       ENDDO
       Current=>Current%Next
       ipart=ipart+1
    ENDDO
    CALL MPI_ALLREDUCE(Weight_Fn,Temp,nx_global+4,mpireal,MPI_SUM,comm,errcode)
    Weight_Fn=Temp
    DO ix=1,nx_global
       IF (Weight_Fn(ix) .GT. 0.0_num) THEN
          Weight_Fn(ix)=Density(ix)/Weight_Fn(ix)
       ENDIF
    ENDDO
    weight=SUM(Weight_Fn(1:nx_global))/REAL(nx_global,num)
#endif
    DEALLOCATE(Weight_Fn)
  END SUBROUTINE SetupParticleListByDensity


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
    REAL(num),PARAMETER :: AM=1.0/REAL(IM,8)

    INTEGER :: k

    idum=XOR(idum,mask)
    k=idum/IQ

    idum=IA*(idum-k*IQ)-IR*k
    IF (idum .LT. 0) idum=idum+IM

    Random=AM*idum
    idum=XOR(idum,mask)

    IF (random .GT. max_rand) max_rand=random

  END FUNCTION Random

END MODULE helper
