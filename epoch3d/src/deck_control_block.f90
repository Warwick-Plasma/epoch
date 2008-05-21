MODULE deck_control_block
  USE shared_data
  USE strings
IMPLICIT NONE
SAVE 
  INTEGER,PARAMETER :: ControlBlockElements =19
  LOGICAL, DIMENSION(ControlBlockElements) :: ControlBlockDone =.FALSE.
  CHARACTER(len=30),DIMENSION(ControlBlockElements) :: ControlBlockName = (/"nx","ny","nz","npart",&
       "nsteps","t_end","x_start","x_end","y_start","y_end","z_start","z_end",&
         "dt_multiplier",&
         "data_dir","restart","restart_snapshot","nprocx","nprocy","nprocz"/)

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
    IF (ControlBlockDone(elementselected)) THEN
       HandleControlDeck=ERR_PRESET_ELEMENT
       RETURN
    ENDIF
    ControlBlockDone(elementselected)=.TRUE.
    HandleControlDeck=ERR_NONE

    SELECT CASE (elementselected)
    CASE(1)
       nx_global=AsInteger(Value,HandleControlDeck)
    CASE(2)
       ny_global=AsInteger(Value,HandleControlDeck)
    CASE(3)
       nz_global=AsInteger(Value,HandleControlDeck)
    CASE(4)
       npart_global=AsInteger(Value,HandleControlDeck)
    CASE(5)
       nsteps=AsInteger(Value,HandleControlDeck)
    CASE(6)
       t_end=AsReal(Value,HandleControlDeck)
    CASE(7)
       x_start=AsReal(Value,HandleControlDeck)
    CASE(8)
       x_end=AsReal(Value,HandleControlDeck)
    CASE(9)
       y_start=AsReal(Value,HandleControlDeck)
    CASE(10)
       y_end=AsReal(Value,HandleControlDeck)
    CASE(11)
       z_start=AsReal(Value,HandleControlDeck)
    CASE(12)
       z_end=AsReal(Value,HandleControlDeck)
    CASE(13)
       dt_multiplier=AsReal(Value,HandleControlDeck)
    CASE(14)
       data_dir=Value(1:MIN(LEN(Value),Data_Dir_Max_Length))
    CASE(15)
       restart=AsLogical(Value,HandleControlDeck)
    CASE(16)
       restart_snapshot=AsInteger(Value,HandleControlDeck)
    CASE(17)
       nprocx=AsInteger(Value,HandleControlDeck)
    CASE(18)
       nprocy=AsInteger(Value,HandleControlDeck)
    CASE(19)
       nprocz=AsInteger(Value,HandleControlDeck)
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
