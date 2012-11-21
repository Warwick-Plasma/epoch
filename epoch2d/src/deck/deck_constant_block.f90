MODULE deck_constant_block

  USE strings
  USE shunt

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: constant_deck_initialise, constant_deck_finalise
  PUBLIC :: constant_block_start, constant_block_end
  PUBLIC :: constant_block_handle_element, constant_block_check

CONTAINS

  SUBROUTINE constant_deck_initialise

  END SUBROUTINE constant_deck_initialise



  SUBROUTINE constant_deck_finalise

  END SUBROUTINE constant_deck_finalise



  SUBROUTINE constant_block_start

  END SUBROUTINE constant_block_start



  SUBROUTINE constant_block_end

  END SUBROUTINE constant_block_end



  FUNCTION constant_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: ix, io
    TYPE(deck_constant), DIMENSION(:), ALLOCATABLE :: buffer
    TYPE(primitive_stack) :: temp
    TYPE(stack_element) :: block

    errcode = c_err_none

    IF (value .EQ. blank) RETURN

    ! First check whether constant already exists
    DO ix = 1, n_deck_constants
      IF (str_cmp(TRIM(element), TRIM(deck_constant_list(ix)%name))) RETURN
    ENDDO

    ! If we're here then then named constant doesn't yet exist, so create it

    ! First issue a warning message if the name overrides a built-in one
    CALL load_block(element, block)
    IF (block%ptype .NE. c_pt_bad .AND. block%ptype .NE. c_pt_null) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'The constant block variable "' // TRIM(element) &
              // '" conflicts with a'
          WRITE(io,*) 'built-in constant name. It will be ignored.'
          WRITE(io,*)
        ENDDO
      ENDIF
    ENDIF

    CALL initialise_stack(temp)
    CALL tokenize(value, temp, errcode)
    IF (errcode .NE. c_err_none) THEN
      CALL deallocate_stack(temp)
      RETURN
    ENDIF

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

  END FUNCTION constant_block_handle_element



  FUNCTION constant_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION constant_block_check

END MODULE deck_constant_block
