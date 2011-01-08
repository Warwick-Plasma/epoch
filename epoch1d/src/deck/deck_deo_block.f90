MODULE deck_deo_block

  USE strings
  USE shunt

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION handle_deo_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_deo_deck
    INTEGER :: ix, err
    TYPE(deferred_execution_object), DIMENSION(:), ALLOCATABLE :: buffer
    TYPE(primitive_stack) :: temp

    handle_deo_deck = c_err_none

    IF (value .EQ. blank) RETURN

    CALL initialise_stack(temp)
    CALL tokenize(value, temp, err)
    IF (err .NE. c_err_none) THEN
      handle_deo_deck = err
      CALL deallocate_stack(temp)
      RETURN
    ENDIF

    ! First check whether constant already exists
    DO ix = 1, n_deferred_execution_objects
      IF (str_cmp(TRIM(element), TRIM(deferred_objects(ix)%name))) THEN
        deferred_objects(ix)%execution_stream = temp
        RETURN
      ENDIF
    ENDDO

    ! If we're here then then named constant doesn't yet exist, so create it

    ! Take a copy of the old list
    IF (n_deferred_execution_objects .GT. 0) THEN
      ALLOCATE(buffer(1:n_deferred_execution_objects))
      buffer = deferred_objects

      ! Allocate the new list
      DEALLOCATE(deferred_objects)
      ALLOCATE(deferred_objects(1:n_deferred_execution_objects+1))

      deferred_objects(1:n_deferred_execution_objects) = buffer
      DEALLOCATE(buffer)
    ELSE
      ! Allocate the new list
      ALLOCATE(deferred_objects(1:n_deferred_execution_objects+1))
    ENDIF

    ! Add the new value
    deferred_objects(n_deferred_execution_objects)%execution_stream = temp
    deferred_objects(n_deferred_execution_objects)%name = element

    n_deferred_execution_objects = n_deferred_execution_objects + 1

  END FUNCTION handle_deo_deck

END MODULE deck_deo_block
