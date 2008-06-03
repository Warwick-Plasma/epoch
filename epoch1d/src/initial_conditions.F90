MODULE initial_conditions

  USE shared_data
  USE strings
  USE helper
  USE partlist


  IMPLICIT NONE

  PRIVATE

  PUBLIC:: HandleCustomBlock,Equilibrium,CheckCustomBlocks


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE Equilibrium


    REAL(num) :: rpos,rvel,rtype,n0,n_real=1.0e0_num

    REAL(num) :: TMax,TMin,lambda_D
    REAL(num) :: x1,x2
    INTEGER :: i
    TYPE(Particle), POINTER :: Cur
    INTEGER :: n, clock,rate=1000,max
    INTEGER :: idum
    REAL(num),DIMENSION(-1:nx_global+2) :: temp
    TYPE(ParticlePointer),DIMENSION(1:nspecies) :: SpeciesList


    max_rand=0.0_num
    max=100000000
    CALL SYSTEM_CLOCK(clock,rate,max)
    idum=-(clock+rank)

    lambda_d=dx/0.4_num

    x1=256.0_num * dx
    x2=324.0_num * dx
    TMax=1.0e8_num
    TMin=1.0e6_num


    Ex=0.0_num
    Ey=0.0_num
    Ez=0.0_num
    Bx=0.0_num
    By=0.0_num
    Bz=0.0_num
    Jx=0.0_num
    Jy=0.0_num
    Jz=0.0_num

    Cur=>MainRoot%Head

!!$    weight=1.0_num
!!$    DO WHILE(ASSOCIATED(Cur))
!!$       Cur%Part_pos=(random(idum)*(x_end-x_start))+x_start
!!$       Cur%Part_species=1
!!$       Cur%Part_p=0.0_num
!!$       Cur%Part_p(1) = 1.0e8_num * m0
!!$       Cur=>Cur%Next
!!$    ENDDO



    CALL LoadParticles(idum,speciesList)
    temp=(TMax * kb * epsilon0)/(lambda_d**2 * Q0**2)
    CALL SetupParticleListByDensity(temp,SpeciesList(2),idum)
    CALL SetupParticleListByDensity(temp,SpeciesList(1),idum)

    IF (rank .EQ. 0) PRINT *,"lambda_d=",lambda_d

    DO ix=-1,nx_global+2
       IF (x_global(ix) .LE. x1) THEN
          temp(ix)=Tmax
       ELSE IF (x_global(ix) .GT. x1 .AND. x_global(ix) .LT. x2) THEN
          temp(ix)=(Tmax-Tmin)/(x1-x2) * (x_global(ix)-x2) + Tmin
       ELSE
          temp(ix)=TMin
       ENDIF
    ENDDO

    CALL SetupParticleListByTemperature(temp,DIR_X,SpeciesList(1),idum)


    temp=TMin
    CALL SetupParticleListByTemperature(temp,DIR_X,SpeciesList(2),idum)
    temp=0.0_num
    CALL SetupParticlesByTemperature(temp,DIR_Y+DIR_Z,idum)



  END SUBROUTINE Equilibrium

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
 !Set the return value to "CheckCustomBlocks = ERR_MISSING_ELEMENTS".

 CheckCustomBlocks=ERR_NONE

END FUNCTION CheckCustomBlocks

END MODULE initial_conditions
