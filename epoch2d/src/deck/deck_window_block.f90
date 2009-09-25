MODULE deck_window_block

  USE shared_data
  USE strings_advanced

  IMPLICIT NONE

CONTAINS

  FUNCTION HandleWindowDeck(Element,Value)
    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleWindowDeck

    HandleWindowDeck=ERR_NONE
    IF (Element .EQ. blank .OR. Value .EQ. blank) RETURN

    IF (StrCmp(Element,"move_window")) THEN
       move_window=AsLogical(Value,HandleWindowDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"window_v_x")) THEN
       window_v_x=AsReal(Value,HandleWindowDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"window_start_time")) THEN
       window_start_time=AsReal(Value,HandleWindowDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"xbc_left_after_move")) THEN
       xbc_left_after_move=AsBC(Value,HandleWindowDeck)
       RETURN
    ENDIF

    IF (StrCmp(Element,"xbc_right_after_move")) THEN
       xbc_right_after_move=AsBC(Value,HandleWindowDeck)
       RETURN
    ENDIF

    HandleWindowDeck=ERR_UNKNOWN_ELEMENT

  END FUNCTION HandleWindowDeck

  FUNCTION CheckWindowBlock()

    INTEGER :: CheckWindowBlock

    !Should do error checking but can't be bothered at the moment
    CheckWindowBlock=ERR_NONE

  END FUNCTION CheckWindowBlock

  SUBROUTINE Window_Start
    xbc_left_after_move=xbc_left
    xbc_right_after_move=xbc_right
  END SUBROUTINE Window_Start

END MODULE deck_window_block
