MODULE deck_deo_block

  USE shared_data
  USE shared_parser_data
  USE strings_advanced
  USE strings
  USE shunt

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION HandleDEODeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleDEODeck
    INTEGER :: ix,iy,ERR
    TYPE(Deferred_Execution_Object),DIMENSION(:),ALLOCATABLE :: Buffer
    TYPE(PrimitiveStack) :: Temp

    HandleDEODeck=ERR_NONE
    Temp%Stackpoint=0

    IF (Value .EQ. Blank) RETURN
    CALL Tokenize(Value,Temp,ERR)
    IF (ERR .NE. ERR_NONE) THEN
       HandleDEODeck=ERR
       RETURN
    ENDIF

    !First check whether constant already exists
    DO ix=1,n_Deferred_Execution_Objects
       IF (StrCmp(TRIM(Element),TRIM(Deferred_Objects(ix)%Name))) THEN
          Deferred_Objects(ix)%Execution_Stream=Temp
          RETURN
       ENDIF
    ENDDO

    !If we're here then then named constant doesn't yet exist, so create it

    !Take a copy of the old list
    IF (n_Deferred_Execution_Objects .GT. 0) THEN
       ALLOCATE(Buffer(1:n_Deferred_Execution_Objects))
       Buffer=Deferred_Objects
       DEALLOCATE(Deferred_Objects)
    ENDIF
    !Allocate the new list
    n_Deferred_Execution_Objects=n_Deferred_Execution_Objects+1
    ALLOCATE(Deferred_Objects(1:n_Deferred_Execution_Objects))
    !If old list not empty then
    IF (n_Deferred_Execution_Objects .GT. 1) THEN
       Deferred_Objects(1:n_Deferred_Execution_Objects-1)=Buffer
       DEALLOCATE(Buffer)
    ENDIF
    !Add the new value
    Deferred_Objects(n_Deferred_Execution_Objects)%Execution_Stream=Temp
    Deferred_Objects(n_Deferred_Execution_Objects)%Name=Element
  END FUNCTION HandleDEODeck

END MODULE deck_deo_block
