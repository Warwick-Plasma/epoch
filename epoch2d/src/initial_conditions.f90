MODULE initial_conditions

  USE shared_data
  USE strings

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: HandleCustomBlock,Equilibrium, CheckCustomBlocks


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE Equilibrium

    REAL(num) :: rpos,rvel,rtype,n0,n_real=1.0e0_num

    REAL(num) :: TMax,TMin
    REAL(num) :: x1,x2,t
    REAL(num) :: alpha=0.2_num
    x1=256.0_num
    x2=324.0_num
    TMax=1.0_num
    TMin=0.01_num


    Ex=0.0_num
    Ey=0.0_num
    Ez=0.0_num
    Bx=0.0_num
    By=0.0_num
    Bz=0.0_num
    Jx=0.0_num
    Jy=0.0_num
    Jz=0.0_num

    !n0 is the pseudoparticle density
    n0=REAL(npart_global,num)/(length_x * length_y)
    !set n_real to the real particle density
    n_real=1.0e15_num
    !weight is then the effective number of real particles represented by a single pseudoparticle
    weight=n_real/n0

    !charge=species(part_species(ipart),1)
    !mass  =species(part_species(ipart),2)
    !Note that this code uses arrays rather than particle types
    DO ipart=1,npart
       Part_p(ipart,:)=0.0_num
       CALL RANDOM_NUMBER(rpos)
       rpos=(rpos*length_x)+x_start
       Part_pos(ipart,1)=rpos

       CALL RANDOM_NUMBER(rpos)
       rpos=(rpos*length_y)+y_start
       Part_pos(ipart,2)=rpos
       !Electron
       IF (ipart .LE. npart/2) THEN
          part_species(ipart) = 2
       ELSE
          part_species(ipart) = 1
       ENDIF
       t=1.0e6_num * EXP((-(Part_pos(ipart,1)-0.017877_num/2.0_num)**2 - (Part_pos(ipart,2)-0.017877_num/2.0_num)**2)*100000.0_num)
       part_p(ipart,1)=MomentumFromTemperature(species(part_species(ipart),2),t)
       part_p(ipart,2)=MomentumFromTemperature(species(part_species(ipart),2),t)
    ENDDO

  END SUBROUTINE Equilibrium

  FUNCTION Gaussian(stdev)

    REAL(num),INTENT(IN) :: stdev
    REAL(num) :: Gaussian

    REAL(num) :: mean=0.0_num
    REAL(num) :: rand1,rand2,w

    !This is a basic polar Box-Muller transform
    !It generates gaussian distributed random numbers
    !The standard deviation (stdev) is related to temperature

    DO
       CALL RANDOM_NUMBER(rand1)
       CALL RANDOM_NUMBER(rand2)

       rand1=2.0_num*rand1 - 1.0_num
       rand2=2.0_num*rand2 - 1.0_num

       w=rand1**2 + rand2**2

       IF (w .LT. 1.0_num) EXIT
    ENDDO

    w = SQRT((-2.0_num * LOG(w) )/w)


    Gaussian = rand1 * w * stdev

  END FUNCTION Gaussian

  FUNCTION MomentumFromTemperature(mass,temperature)

    REAL(num),INTENT(IN) :: mass,temperature
    REAL(num) :: MomentumFromTemperature

    REAL(num) :: mean=0.0_num
    REAL(num) :: stdev
    REAL(num) :: rand1,rand2,w

    !This is a basic polar Box-Muller transform
    !It generates gaussian distributed random numbers
    !The standard deviation (stdev) is related to temperature

    stdev=SQRT(kb * temperature * mass)

    DO
       CALL RANDOM_NUMBER(rand1)
       CALL RANDOM_NUMBER(rand2)

       rand1=2.0_num*rand1 - 1.0_num
       rand2=2.0_num*rand2 - 1.0_num

       w=rand1**2 + rand2**2

       IF (w .LT. 1.0_num) EXIT
    ENDDO

    w = SQRT((-2.0_num * LOG(w) )/w)


    MomentumFromTemperature = rand1 * w * stdev

  END FUNCTION MomentumFromTemperature
  !-----------------------------------------------------------------------------
  !These functions contain the user input deck elements
  !-----------------------------------------------------------------------------

  FUNCTION HandleCustomBlock(blockname,Element,Value)

    CHARACTER(len=30),INTENT(IN)::blockname,Element,Value
    INTEGER :: HandleCustomBlock

    !The following line must always be present
    HandleCustomBlock=ERR_UNKNOWN_BLOCK

  END FUNCTION HandleCustomBlock

  FUNCTION CheckCustomBlocks

    INTEGER :: CheckCustomBlocks

    !This subroutine is to allow you to force the code to bomb out if an essential element
    !Of the input deck is missing. If you either don't want to check ,are not extending the
    !Input deck, or all elements are set then set "CheckCustomBlocks = ERR_NONE". Otherwise
    !Set the return value to "CheckCustomBlocks = ERR_MISSING_ELEMENTS". If you want an error
    !Message then use the function "Report_Deck_Error(ErrorString)

    CheckCustomBlocks=ERR_NONE

  END FUNCTION CheckCustomBlocks


END MODULE initial_conditions
