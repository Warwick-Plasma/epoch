MODULE initial_conditions

  USE shared_data
  USE strings
  USE helper

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
    REAL(num) :: x1,x2
    REAL(num) :: alpha=0.2_num
    TYPE(PARTICLE), POINTER :: Cur
    INTEGER :: clock,idum


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

    CALL SYSTEM_CLOCK(COUNT=clock)
    idum=clock+rank

    !n0 is the pseudoparticle density
    n0=REAL(npart_global,num)/(length_x * length_y * length_z)
    !set n_real to the real particle density
    n_real=1.0e15_num
    !weight is then the effective number of real particles represented by a single pseudoparticle
    weight=n_real/n0

    !charge=species(Cur%part_species,1)
    !mass  =species(Cur%part_species,2)


    Cur=>Head
    ipart=0
    DO WHILE(ASSOCIATED(Cur))
       Cur%Part_p=0.0_num
       rpos=random(idum)

       rpos=(rpos*length_x_local)+x_start_local
       Cur%Part_pos(1)=rpos

       CALL RANDOM_NUMBER(rpos)
       rpos=(rpos*length_y_local)+y_start_local
       Cur%Part_pos(2)=rpos

       CALL RANDOM_NUMBER(rpos)
       rpos=(rpos*length_z_local)+z_start_local
       Cur%Part_pos(3)=rpos

       IF (Cur%Part_pos(1) .GT. x_end_local) Cur%part_pos(1) =  x_end_local - (Cur%part_pos(1)-x_end_local)
       IF (Cur%Part_pos(2) .GT. y_end_local) Cur%part_pos(2) =  y_end_local - (Cur%part_pos(2)-y_end_local)
       IF (Cur%Part_pos(3) .GT. z_end_local) Cur%part_pos(3) =  z_end_local - (Cur%part_pos(3)-z_end_local)

       !Electron
       IF (ipart .LE. npart/2) THEN
          Cur%part_species = 2
          Cur%Part_P(1)=0.0_num!Cur%Part_pos(1)!MomentumFromTemperature(species(2,2),EXP(-(Cur%Part_pos(1)*100.0_num)**2 + &
          !(Cur%Part_pos(2)*100.0_num)**2 + (Cur%Part_pos(3)*100.0_num)**2))
       ELSE
          Cur%part_species = 1
       ENDIF
       !       Cur%Label=rank
       Cur=>Cur%Next
       ipart=ipart+1
    ENDDO

  END SUBROUTINE Equilibrium

  !-----------------------------------------------------------------------------
  !These functions contain the user input deck elements
  !-----------------------------------------------------------------------------

  FUNCTION HandleCustomBlock(blockname,Element,Value)

    CHARACTER(len=30),INTENT(IN)::blockname,Element,Value
    INTEGER :: HandleCustomBlock

    !The following line must always be present
    HandleCustomBlock=ERR_UNKNOWN_BLOCK

    IF (StrCmp(blockname,"laser")) THEN

       HandleCustomBlock=ERR_NONE

       IF (StrCmp(Element,"amplitude")) THEN
          las_amp=AsReal(Value,HandleCustomBlock)
          RETURN
       ENDIF
       IF (StrCmp(Element,"period")) THEN
          las_omega=2.0_num * pi /AsReal(Value,HandleCustomBlock)
          RETURN
       ENDIF

       IF (StrCmp(Element,"t_drop")) THEN
          las_drop=AsReal(Value,HandleCustomBlock)
          RETURN
       ENDIF

       HandleCustomBlock=ERR_UNKNOWN_ELEMENT

    ENDIF


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
