MODULE deck_window_block

  USE strings_advanced

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: window_deck_initialise, window_deck_finalise
  PUBLIC :: window_block_start, window_block_end
  PUBLIC :: window_block_handle_element, window_block_check

CONTAINS

  SUBROUTINE window_deck_initialise

  END SUBROUTINE window_deck_initialise



  SUBROUTINE window_deck_finalise

  END SUBROUTINE window_deck_finalise



  SUBROUTINE window_block_start

    IF (deck_state .NE. c_ds_first) RETURN

    bc_x_min_after_move = bc_field(c_bd_x_min)
    bc_x_max_after_move = bc_field(c_bd_x_max)

  END SUBROUTINE window_block_start



  SUBROUTINE window_block_end

  END SUBROUTINE window_block_end



  FUNCTION window_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state .NE. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, 'move_window')) THEN
      move_window = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'window_v_x')) THEN
      window_v_x = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'window_start_time')) THEN
      window_start_time = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'bc_x_min_after_move') &
        .OR. str_cmp(element, 'xbc_left_after_move')) THEN
      bc_x_min_after_move = as_bc(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'bc_x_max_after_move') &
        .OR. str_cmp(element, 'xbc_right_after_move')) THEN
      bc_x_max_after_move = as_bc(value, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION window_block_handle_element



  FUNCTION window_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION window_block_check

END MODULE deck_window_block
