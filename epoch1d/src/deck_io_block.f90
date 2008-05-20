MODULE deck_io_block

  USE shared_data
  USE strings

  IMPLICIT NONE

  SAVE
  INTEGER, PARAMETER :: N_Var_Special = 3
  INTEGER,PARAMETER :: IOBlockElements=num_vars_to_dump + N_Var_Special
  LOGICAL, DIMENSION(IOBlockElements)  :: IOBlockDone
  CHARACTER(len=30),DIMENSION(IOBlockElements) :: IOBlockName=(/"dt_snapshot","full_dump_every","force_final_to_be_restartable","particles","grid","px","py","pz","vx","vy","vz",&
       "ex","ey","ez","bx","by","bz","jx","jy","jz","charge","mass","weight",&
       "ekbar","mass_density","charge_density","species","number_density"/)

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

    IF (elementselected .EQ. 0) THEN
       HandleIODeck=ERR_UNKNOWN_ELEMENT
       RETURN
    ENDIF
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
       IF (AsLogical(Value,HandleIODeck)) THEN
          DumpMask(1:5) = IOR(DumpMask(1:5),IO_RESTARTABLE)
          DumpMask(9:14) = IOR(DumpMask(9:14),IO_RESTARTABLE) 
          DumpMask(20) = IOR(DumpMask(20),IO_RESTARTABLE)
          DumpMask(24) = IOR(DumpMask(24),IO_RESTARTABLE)
       ENDIF
    END SELECT

    IF (elementselected .LE. N_Var_Special) RETURN
    IF (elementselected .GT. N_Var_Special) DumpMask(elementselected-N_Var_Special)=IOR(DumpMask(elementselected-N_Var_Special),AsDumpParam(Value,HandleIODeck))

  END FUNCTION HandleIODeck


  FUNCTION CheckIOBlock

    INTEGER :: CheckIOBlock,Index

    !Just assume that anything not included except for the first two elements is not wanted

    CheckIOBlock=ERR_NONE

    DO index=1,3
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


  END FUNCTION CheckIOBlock
END MODULE deck_io_block
