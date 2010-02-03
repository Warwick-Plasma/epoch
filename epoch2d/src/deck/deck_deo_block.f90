MODULE deck_deo_block

  USE shared_data
  USE shared_parser_data
  USE strings_advanced
  USE strings
  USE shunt

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION handle_deo_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_deo_deck
    INTEGER :: ix,err
    TYPE(deferred_execution_object),DIMENSION(:),ALLOCATABLE :: buffer
    TYPE(primitive_stack) :: temp

    handle_deo_deck=c_err_none

    IF (value .EQ. blank) RETURN
   temp%stack_point=0
    CALL tokenize(value,temp,err)
    IF (err .NE. c_err_none) THEN
       handle_deo_deck=err
       RETURN
    ENDIF

    !First check whether constant already exists
    DO ix=1,n_deferred_execution_objects
       IF (str_cmp(TRIM(element),TRIM(deferred_objects(ix)%name))) THEN
          deferred_objects(ix)%execution_stream=temp
          RETURN
       ENDIF
    ENDDO

    !If we're here then then named constant doesn't yet exist, so create it

    !Take a copy of the old list
    IF (n_deferred_execution_objects .GT. 0) THEN
       ALLOCATE(buffer(1:n_deferred_execution_objects))
       buffer=deferred_objects
       DEALLOCATE(deferred_objects)
    ENDIF
    !Allocate the new list
    n_deferred_execution_objects=n_deferred_execution_objects+1
    ALLOCATE(deferred_objects(1:n_deferred_execution_objects))
    !If old list not empty then
    IF (n_deferred_execution_objects .GT. 1) THEN
       deferred_objects(1:n_deferred_execution_objects-1)=buffer
       DEALLOCATE(buffer)
    ENDIF
    !Add the new value
    deferred_objects(n_deferred_execution_objects)%execution_stream=temp
    deferred_objects(n_deferred_execution_objects)%name=element
  END FUNCTION handle_deo_deck

END MODULE deck_deo_block
