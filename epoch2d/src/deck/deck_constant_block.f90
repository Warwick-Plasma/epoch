MODULE deck_constant_block

  USE shared_data
  USE shared_parser_data
  USE strings_advanced
  USE strings

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION HandleConstantDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleConstantDeck
    INTEGER :: ix,iy,ERR
    TYPE(Deck_Constant),DIMENSION(:),ALLOCATABLE :: Buffer

    HandleConstantDeck=ERR_NONE

    IF (Value .EQ. Blank) RETURN

    !First check whether constant already exists
    DO ix=1,n_Deck_Constants
       IF (StrCmp(TRIM(Element),TRIM(Deck_Constant_List(ix)%Name))) THEN
          Deck_Constant_List(ix)%Value=AsReal(Value,HandleConstantDeck)
          RETURN
       ENDIF
    ENDDO

    !If we're here then then named constant doesn't yet exist, so create it

    !Take a copy of the old list
    IF (n_Deck_Constants .GT. 0) THEN
       ALLOCATE(Buffer(1:n_Deck_Constants))
       Buffer=Deck_Constant_List
       DEALLOCATE(Deck_Constant_List)
    ENDIF
    !Allocate the new list
    n_Deck_Constants=n_Deck_Constants+1
    ALLOCATE(Deck_Constant_List(1:n_Deck_Constants))
    !If old list not empty then
    IF (n_Deck_Constants .GT. 1) THEN
       Deck_Constant_List(1:n_Deck_Constants-1)=Buffer
       DEALLOCATE(Buffer)
    ENDIF
    !Add the new value
    Deck_Constant_List(n_Deck_Constants)%Value=AsReal(Value,HandleConstantDeck)
    Deck_Constant_List(n_Deck_Constants)%Name=Element
  END FUNCTION HandleConstantDeck

END MODULE deck_constant_block
