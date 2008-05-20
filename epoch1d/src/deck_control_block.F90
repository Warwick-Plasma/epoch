MODULE deck_control_block
  USE shared_data
  USE strings
IMPLICIT NONE
SAVE 
  INTEGER,PARAMETER :: ControlBlockElements =13
  LOGICAL, DIMENSION(ControlBlockElements) :: ControlBlockDone =.FALSE.
  CHARACTER(len=30),DIMENSION(ControlBlockElements) :: ControlBlockName = (/"nx","npart",&
       "nsteps","t_end","x_start","x_end",&
         "dt_multiplier",&
         "data_dir","restart","restart_snapshot","domain","dlb","dlb_threshold"/)

CONTAINS

  FUNCTION HandleControlDeck(Element,Value)

    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleControlDeck
    INTEGER :: loop,elementselected
    HandleControlDeck=ERR_UNKNOWN_ELEMENT

    elementselected=0

    DO loop=1,ControlBlockElements
       IF(StrCmp(Element,TRIM(ADJUSTL(ControlBlockName(loop))))) THEN
          elementselected=loop
          EXIT
       ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (ControlBlockDone(elementselected) .AND. elementselected .NE. 8) THEN
       HandleControlDeck=ERR_PRESET_ELEMENT
       RETURN
    ENDIF
    ControlBlockDone(elementselected)=.TRUE.
    HandleControlDeck=ERR_NONE

#ifdef INPUTFILE
    !If using INPUTFILE input specifier then data_dir is already known
    ControlBlockDone(8)=.TRUE.
#endif

    SELECT CASE (elementselected)
    CASE(1)
       nx_global=AsInteger(Value,HandleControlDeck)
    CASE(2)
       npart_global=AsInteger(Value,HandleControlDeck)
    CASE(3)
       nsteps=AsInteger(Value,HandleControlDeck)
    CASE(4)
       t_end=AsReal(Value,HandleControlDeck)
    CASE(5)
       x_start=AsReal(Value,HandleControlDeck)
    CASE(6)
       x_end=AsReal(Value,HandleControlDeck)
    CASE(7)
       dt_multiplier=AsReal(Value,HandleControlDeck)
    CASE(8)
#ifdef INPUTFILE
       IF (rank .EQ. 0) PRINT *,"***WARNING*** Input deck attempting to reset output directory, ignoring"
#else
       data_dir=Value(1:MIN(LEN(Value),Data_Dir_Max_Length))
#endif
    CASE(9)
       restart=AsLogical(Value,HandleControlDeck)
    CASE(10)
       restart_snapshot=AsInteger(Value,HandleControlDeck)
    CASE(11)
       Domain=AsDomain(Value,HandleControlDeck)
    CASE(12)
       DLB=AsLogical(Value,HandleControlDeck)
    CASE(13)
       DLB_Threshold=AsReal(Value,HandleControlDeck)
    END SELECT

  END FUNCTION HandleControlDeck

  FUNCTION CheckControlBlock
    INTEGER :: CheckControlBlock,Index

    CheckControlBlock=ERR_NONE

    DO index=1,ControlBlockElements
       IF (.NOT. ControlBlockDone(index)) THEN
          IF (rank .EQ. 0) THEN
             PRINT *,"***ERROR***"
             PRINT *,"Required control block element ",TRIM(ADJUSTL(ControlBlockName(index))), " absent. Please create this entry in the input deck"
             WRITE(40,*) ""
             WRITE(40,*) "***ERROR***"
             WRITE(40,*) "Required control block element ",TRIM(ADJUSTL(ControlBlockName(index))), " absent. Please create this entry in the input deck"             
          ENDIF
          CheckControlBlock=ERR_MISSING_ELEMENTS
       ENDIF
    ENDDO
  END FUNCTION CheckControlBlock

END MODULE
