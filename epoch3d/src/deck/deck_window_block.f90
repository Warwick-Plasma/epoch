MODULE deck_window_block

  USE shared_data
  USE strings_advanced

  IMPLICIT NONE

CONTAINS

  FUNCTION handle_window_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_window_deck

    handle_window_deck=c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element,"move_window")) THEN
       move_window=as_logical(value,handle_window_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"window_v_x")) THEN
       window_v_x=as_real(value,handle_window_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"window_start_time")) THEN
       window_start_time=as_real(value,handle_window_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"xbc_left_after_move")) THEN
       xbc_left_after_move=as_bc(value,handle_window_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"xbc_right_after_move")) THEN
       xbc_right_after_move=as_bc(value,handle_window_deck)
       RETURN
    ENDIF

    handle_window_deck=c_err_unknown_element

  END FUNCTION handle_window_deck

  FUNCTION check_window_block()

    INTEGER :: check_window_block

    !Should do error checking but can't be bothered at the moment
    check_window_block=c_err_none

  END FUNCTION check_window_block

  SUBROUTINE window_start
    xbc_left_after_move=xbc_left
    xbc_right_after_move=xbc_right
  END SUBROUTINE window_start

END MODULE deck_window_block
