MODULE deck_io_block

  USE shared_data
  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER,PARAMETER :: IOBlockElements=26
  LOGICAL, DIMENSION(IOBlockElements)  :: IOBlockDone
  CHARACTER(len=30),DIMENSION(IOBlockElements) :: IOBlockName=(/"dt_snapshot","full_dump_every","particles","grid","px","py","pz","vx","vy","vz",&
       "ex","ey","ez","bx","by","bz","jx","jy","jz","charge","mass",&
       "temperatures","mass_density","charge_density","particle_weight","species_id"/)

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
    END SELECT

    IF (elementselected .LE. 2) RETURN
    IF (elementselected .GT. 2) DumpMask(elementselected-2)=AsDumpParam(Value,HandleIODeck)

  END FUNCTION HandleIODeck


     FUNCTION CheckIOBlock

       INTEGER :: CheckIOBlock,Index

       !Just assume that anything not included except for the first two elements is not wanted

       CheckIOBlock=ERR_NONE

       DO index=1,2
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

       IF (.NOT. IOBlockDone(26) .AND. rank .EQ. 0) PRINT *,"***WARNING*** You have selected not to output particle species information. This means that it is impossible to restart the code from these data dumps"
       IF (.NOT. IOBlockDone(25) .AND. rank .EQ. 0) PRINT *,"***WARNING*** You have selected not to output particle weight information. This means that it is impossible to restart the code from these data dumps"

     END FUNCTION CheckIOBlock
   END MODULE deck_io_block
