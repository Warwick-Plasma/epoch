MODULE deck_constant_block

  USE strings
  USE shunt

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION handle_constant_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_constant_deck
    INTEGER :: ix, err
    TYPE(deck_constant), DIMENSION(:), ALLOCATABLE :: buffer
    TYPE(primitive_stack) :: temp

    handle_constant_deck = c_err_none

    IF (value .EQ. blank) RETURN

    CALL initialise_stack(temp)
    CALL tokenize(value, temp, err)
    IF (err .NE. c_err_none) THEN
      handle_constant_deck = err
      CALL deallocate_stack(temp)
      RETURN
    ENDIF

    ! First check whether constant already exists
    DO ix = 1, n_deck_constants
      IF (str_cmp(TRIM(element), TRIM(deck_constant_list(ix)%name))) THEN
        deck_constant_list(ix)%execution_stream = temp
        RETURN
      ENDIF
    ENDDO

    ! If we're here then then named constant doesn't yet exist, so create it

    ! Take a copy of the old list
    IF (n_deck_constants .GT. 0) THEN
      ALLOCATE(buffer(1:n_deck_constants))
      buffer = deck_constant_list

      ! Allocate the new list
      DEALLOCATE(deck_constant_list)
      ALLOCATE(deck_constant_list(1:n_deck_constants+1))

      deck_constant_list(1:n_deck_constants) = buffer
      DEALLOCATE(buffer)
    ELSE
      ! Allocate the new list
      ALLOCATE(deck_constant_list(1:n_deck_constants+1))
    ENDIF

    ! Add the new value
    deck_constant_list(n_deck_constants+1)%execution_stream = temp
    deck_constant_list(n_deck_constants+1)%name = element

    n_deck_constants = n_deck_constants + 1

  END FUNCTION handle_constant_deck

END MODULE deck_constant_block
