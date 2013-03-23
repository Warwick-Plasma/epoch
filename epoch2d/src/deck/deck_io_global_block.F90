MODULE deck_io_global_block

  USE strings

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: io_global_deck_initialise, io_global_deck_finalise
  PUBLIC :: io_global_block_start, io_global_block_end
  PUBLIC :: io_global_block_handle_element, io_global_block_check

CONTAINS

  SUBROUTINE io_global_deck_initialise

  END SUBROUTINE io_global_deck_initialise



  SUBROUTINE io_global_deck_finalise

  END SUBROUTINE io_global_deck_finalise



  SUBROUTINE io_global_block_start

    INTEGER :: io, ierr
    CHARACTER(LEN=c_max_string_length) :: list_filename

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (.NOT. new_style_io_block) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'You cannot use an "output_global" block in ', &
              'conjunction with ', 'unnamed "output" blocks.'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

  END SUBROUTINE io_global_block_start



  SUBROUTINE io_global_block_end

  END SUBROUTINE io_global_block_end



  FUNCTION io_global_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, 'force_first_to_be_restartable')) THEN
      force_first_to_be_restartable = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'force_final_to_be_restartable')) THEN
      force_final_to_be_restartable = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'use_offset_grid')) THEN
      use_offset_grid = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'dump_source_code')) THEN
      dump_source_code = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'dump_input_decks')) THEN
      dump_input_decks = as_logical(value, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION io_global_block_handle_element



  FUNCTION io_global_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none

  END FUNCTION io_global_block_check

END MODULE deck_io_global_block
