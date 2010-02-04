MODULE custom_deck

  USE shared_data
  USE strings
  USE strings_advanced
  USE shunt
  USE evaluator
  USE shared_parser_data
  IMPLICIT NONE

CONTAINS

  ! -----------------------------------------------------------------------------
  ! These functions contain the user input deck elements
  ! -----------------------------------------------------------------------------

  FUNCTION handle_custom_block(block_name, element, value)

    CHARACTER(LEN=string_length), INTENT(IN) :: block_name, element, value
    INTEGER :: handle_custom_block

    ! The following line must always be present
    handle_custom_block = c_err_unknown_block

  END FUNCTION handle_custom_block



  FUNCTION check_custom_blocks()

    INTEGER :: check_custom_blocks

    ! This subroutine is to allow you to force the code to bomb out if an essential element
    ! Of the input deck is missing. If you either don't want to check, are not extending the
    ! Input deck, or all elements are set then set "check_custom_blocks = c_err_none". Otherwise
    ! Set the return value to "check_custom_blocks = c_err_missing_elements".

    check_custom_blocks = c_err_none

  END FUNCTION check_custom_blocks

END MODULE custom_deck
