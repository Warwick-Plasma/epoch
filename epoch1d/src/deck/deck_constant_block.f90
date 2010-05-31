MODULE deck_constant_block

  USE strings_advanced

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION handle_constant_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_constant_deck
    INTEGER :: ix
    TYPE(deck_constant), DIMENSION(:), ALLOCATABLE :: buffer

    handle_constant_deck = c_err_none

    IF (value .EQ. blank) RETURN

    ! First check whether constant already exists
    DO ix = 1, n_deck_constants
      IF (str_cmp(TRIM(element), TRIM(deck_constant_list(ix)%name))) THEN
        deck_constant_list(ix)%value = as_real(value, handle_constant_deck)
        RETURN
      ENDIF
    ENDDO

    ! If we're here then then named constant doesn't yet exist, so create it

    ! Take a copy of the old list
    IF (n_deck_constants .GT. 0) THEN
      ALLOCATE(buffer(1:n_deck_constants))
      buffer = deck_constant_list
      DEALLOCATE(deck_constant_list)
    ENDIF

    ! Allocate the new list
    ALLOCATE(deck_constant_list(1:n_deck_constants+1))

    ! If old list not empty then
    IF (n_deck_constants .GE. 1) THEN
      deck_constant_list(1:n_deck_constants) = buffer
      DEALLOCATE(buffer)
    ENDIF

    ! Add the new value
    deck_constant_list(n_deck_constants+1)%value = &
        as_real(value, handle_constant_deck)
    deck_constant_list(n_deck_constants+1)%name = element

    n_deck_constants = n_deck_constants+1

  END FUNCTION handle_constant_deck

END MODULE deck_constant_block
