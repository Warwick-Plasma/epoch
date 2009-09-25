MODULE deck_boundaries_block

  USE shared_data
  USE strings
  USE strings_advanced

  IMPLICIT NONE

  SAVE

  INTEGER,PARAMETER :: BoundaryBlockElements=2
  LOGICAL, DIMENSION(BoundaryBlockElements)  :: BoundaryBlockDone
  CHARACTER(len=EntryLength),DIMENSION(BoundaryBlockElements) :: BoundaryBlockName=(/"xbc_left","xbc_right"/)

CONTAINS



  FUNCTION HandleBoundaryDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleBoundaryDeck
    INTEGER :: loop,elementselected

    HandleBoundaryDeck=ERR_UNKNOWN_ELEMENT

    elementselected=0

    DO loop=1,BoundaryBlockElements
       IF(StrCmp(Element,TRIM(ADJUSTL(BoundaryBlockName(loop))))) THEN
          elementselected=loop
          EXIT
       ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (BoundaryBlockDone(elementselected)) THEN
       HandleBoundaryDeck=ERR_PRESET_ELEMENT
       RETURN
    ENDIF
    BoundaryBlockDone(elementselected)=.TRUE.
    HandleBoundaryDeck=ERR_NONE

    SELECT CASE (elementselected)
    CASE(1)
       xbc_left=AsBC(Value,HandleBoundaryDeck)
    CASE(2)
       xbc_right=AsBC(Value,HandleBoundaryDeck)
    END SELECT

  END FUNCTION HandleBoundaryDeck

  FUNCTION CheckBoundaryBlock()

    INTEGER :: CheckBoundaryBlock
    INTEGER :: index

    CheckBoundaryBlock=ERR_NONE

    DO index=1,BoundaryBlockElements
       IF (.NOT. BoundaryBlockDone(index)) THEN
          IF (rank .EQ. 0) THEN
             PRINT *,"***ERROR***"
             PRINT *,"Required boundary block element ",TRIM(ADJUSTL(BoundaryBlockName(index))), " absent. Please create this entry in the input deck"
             WRITE(40,*) ""
             WRITE(40,*) "***ERROR***"
             WRITE(40,*) "Required boundary block element ",TRIM(ADJUSTL(BoundaryBlockName(index))), " absent. Please create this entry in the input deck"   
          ENDIF
          CheckBoundaryBlock = ERR_MISSING_ELEMENTS
       ENDIF
    ENDDO

  END FUNCTION CheckBoundaryBlock

END MODULE deck_boundaries_block
