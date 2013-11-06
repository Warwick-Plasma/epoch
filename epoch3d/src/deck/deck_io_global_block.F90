MODULE deck_io_global_block

  USE strings_advanced

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: io_global_deck_initialise, io_global_deck_finalise
  PUBLIC :: io_global_block_start, io_global_block_end
  PUBLIC :: io_global_block_handle_element, io_global_block_check
  LOGICAL, SAVE :: dump_first, dump_last, got_dump_first, got_dump_last

CONTAINS

  SUBROUTINE io_global_deck_initialise

    time_start  = -1.0_num
    time_stop   = HUGE(1.0_num)
    nstep_start = -1
    nstep_stop  = HUGE(1)
    got_dump_first = .FALSE.
    got_dump_last = .FALSE.
    ! Set point data buffer size to 64MB by default.
    sdf_buffer_size = 64 * 1024 * 1024
    filesystem = ''

  END SUBROUTINE io_global_deck_initialise



  SUBROUTINE io_global_deck_finalise

    INTEGER :: i

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (got_dump_first) THEN
      DO i = 1, n_io_blocks
        io_block_list(i)%dump_first = dump_first
      ENDDO
    ENDIF

    IF (got_dump_last) THEN
      DO i = 1, n_io_blocks
        io_block_list(i)%dump_last = dump_last
      ENDDO
    ENDIF

  END SUBROUTINE io_global_deck_finalise



  SUBROUTINE io_global_block_start

    INTEGER :: io, iu, ierr

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (.NOT. new_style_io_block) THEN
      IF (rank .EQ. 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
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
      force_first_to_be_restartable = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'force_final_to_be_restartable') &
        .OR. str_cmp(element, 'force_last_to_be_restartable')) THEN
      force_final_to_be_restartable = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'use_offset_grid')) THEN
      use_offset_grid = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'time_start')) THEN
      time_start = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'time_stop')) THEN
      time_stop = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nstep_start')) THEN
      nstep_start = as_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'nstep_stop')) THEN
      nstep_stop = as_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'dump_first')) THEN
      got_dump_first = .TRUE.
      dump_first = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'dump_last') &
        .OR. str_cmp(element, 'dump_final')) THEN
      got_dump_last = .TRUE.
      dump_last = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'sdf_buffer_size')) THEN
      sdf_buffer_size = as_long_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'filesystem')) THEN
      filesystem = TRIM(value) // ':'
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION io_global_block_handle_element



  FUNCTION io_global_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none

  END FUNCTION io_global_block_check

END MODULE deck_io_global_block
