MODULE deck_io_block

  USE shared_data
  USE strings
  USE strings_advanced

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: N_Var_Special = 7
  INTEGER,PARAMETER :: IOBlockElements=N_Var_Special+num_vars_to_dump
  LOGICAL, DIMENSION(IOBlockElements)  :: IOBlockDone
  CHARACTER(len=EntryLength),DIMENSION(IOBlockElements) :: IOBlockName=(/"dt_snapshot","full_dump_every","restart_dump_every","force_final_to_be_restartable","use_offset_grid","use_extended_io",&
       "extended_io_file","particles","grid","px","py","pz","vx","vy","vz",&
       "ex","ey","ez","bx","by","bz","jx","jy","jz","charge","mass",&
       "ekbar","mass_density","charge_density","number_density","particle_weight","species_id","distribution_functions","particle_probes","temperature","ejected_particles"/)

CONTAINS



  FUNCTION HandleIODeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleIODeck
    INTEGER :: loop,elementselected

    HandleIODeck=ERR_UNKNOWN_ELEMENT

    elementselected=0

    DO loop=1,IOBlockElements
       IF(StrCmp(Element,TRIM(ADJUSTL(IOBlockName(loop))))) THEN
          elementselected=loop
          EXIT
       ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (IOBlockDone(elementselected)) THEN
       HandleIODeck=ERR_PRESET_ELEMENT
       RETURN
    ENDIF
    IOBlockDone(elementselected)=.TRUE.
    HandleIODeck=ERR_NONE

    SELECT CASE (elementselected)
    CASE(1)
       dt_snapshots=AsReal(Value,HandleIODeck)
    CASE(2)
       full_dump_every=AsInteger(Value,HandleIODeck)
    CASE(3)
       restart_dump_every=AsInteger(Value,HandleIODeck)
    CASE(4)
       force_final_to_be_restartable=AsLogical(Value,HandleIODeck)
    CASE(5)
       use_offset_grid=AsLogical(Value,HandleIODeck)
    CASE(6)
       use_extended_io=AsLogical(Value,HandleIODeck)
    CASE(7)
       extended_io_file=TRIM(Value)
    END SELECT

    IF (elementselected .LE. N_Var_Special) RETURN
    IF (elementselected .GT. N_Var_Special) DumpMask(elementselected-N_Var_Special)=AsReal(Value,HandleIODeck)

    !If setting dumpmask for particle probes then report if the code wasn't compiled for particle probes
#ifndef PARTICLE_PROBES   
    IF (elementselected-N_Var_Special .EQ. 27) THEN
       HandleIODeck=ERR_PP_OPTIONS_WRONG
       Extended_Error_String="-DPARTICLE_PROBES"
    ENDIF
#endif

  END FUNCTION HandleIODeck


  FUNCTION CheckIOBlock()

    INTEGER :: CheckIOBlock,Index

    !Just assume that anything not included except for the compulsory elements is not wanted
    CheckIOBlock=ERR_NONE

    !If not requesting extended io then don't check for extended_io_file
    IF (.NOT. IOBlockDone(6) .OR. .NOT. use_extended_io) IOBlockDone(6:7)=.TRUE.

    !Particle Positions
    DumpMask(1:5) = IOR(DumpMask(1:5),IO_RESTARTABLE)
    !Fields
    DumpMask(9:14) = IOR(DumpMask(9:14),IO_RESTARTABLE) 
    !Weight and species info
    DumpMask(24:25) = IOR(DumpMask(24:25),IO_RESTARTABLE)

    DO index=1,N_Var_Special
       IF (.NOT. IOBlockDone(index)) THEN
          IF (rank .EQ. 0) THEN
             PRINT *,"***ERROR***"
             PRINT *,"Required output block element ",TRIM(ADJUSTL(IOBlockName(index))), " absent. Please create this entry in the input deck"
             WRITE(40,*) ""
             WRITE(40,*) "***ERROR***"
             WRITE(40,*) "Required output block element ",TRIM(ADJUSTL(IOBlockName(index))), " absent. Please create this entry in the input deck"   
          ENDIF
          CheckIOBlock=ERR_MISSING_ELEMENTS
       ENDIF
    ENDDO

!!$    IF (.NOT. IOBlockDone(26) .AND. rank .EQ. 0) PRINT *,"***WARNING*** You have selected not to output particle species information. This means that it is impossible to restart the code from these data dumps"
!!$    IF (.NOT. IOBlockDone(25) .AND. rank .EQ. 0) PRINT *,"***WARNING*** You have selected not to output particle weight information. This means that it is impossible to restart the code from these data dumps"

  END FUNCTION CheckIOBlock
END MODULE deck_io_block
