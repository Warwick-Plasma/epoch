MODULE custom_deck

  USE shared_data

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! These functions contain the user input deck elements
  !----------------------------------------------------------------------------

  FUNCTION custom_blocks_handle_element(block_name, element, value) &
      RESULT(errcode)

    CHARACTER(LEN=string_length), INTENT(IN) :: block_name, element, value
    INTEGER :: errcode

    ! The following line must always be present
    errcode = c_err_unknown_block

  END FUNCTION custom_blocks_handle_element



  FUNCTION custom_blocks_check() RESULT(errcode)

    INTEGER :: errcode

    ! This subroutine is to allow you to force the code to bomb out if an
    ! essential element of the input deck is missing. If you either don't
    ! want to check, are not extending the input deck, or all elements are
    ! set then set 'errcode = c_err_none'. Otherwise set the
    ! return value to 'errcode = c_err_missing_elements'.

    errcode = c_err_none

  END FUNCTION custom_blocks_check

END MODULE custom_deck
