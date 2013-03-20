MODULE deck_fields_block

  USE strings_advanced
  USE simple_io

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: fields_deck_initialise, fields_deck_finalise
  PUBLIC :: fields_block_start, fields_block_end
  PUBLIC :: fields_block_handle_element, fields_block_check

  INTEGER(KIND=MPI_OFFSET_KIND) :: offset

CONTAINS

  SUBROUTINE fields_deck_initialise

  END SUBROUTINE fields_deck_initialise



  SUBROUTINE fields_deck_finalise

  END SUBROUTINE fields_deck_finalise



  SUBROUTINE fields_block_start

    offset = 0

  END SUBROUTINE fields_block_start



  SUBROUTINE fields_block_end

  END SUBROUTINE fields_block_end



  FUNCTION fields_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    CHARACTER(LEN=string_length) :: filename
    INTEGER :: err
    LOGICAL :: got_file

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, 'offset')) THEN
      offset = as_long_integer_simple(value, errcode)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, err)

    IF (str_cmp(element, 'ex')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ex, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, ex, &
            -2, nx+3, -2, ny+3, -2, nz+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'ey')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ey, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, ey, &
            -2, nx+3, -2, ny+3, -2, nz+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'ez')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, ez, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, ez, &
            -2, nx+3, -2, ny+3, -2, nz+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'bx')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bx, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, bx, &
            -2, nx+3, -2, ny+3, -2, nz+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'by')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, by, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, by, &
            -2, nx+3, -2, ny+3, -2, nz+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'bz')) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, bz, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, bz, &
            -2, nx+3, -2, ny+3, -2, nz+3, errcode)
      ENDIF
      RETURN
    ENDIF

  END FUNCTION fields_block_handle_element



  FUNCTION fields_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION fields_block_check

END MODULE deck_fields_block
