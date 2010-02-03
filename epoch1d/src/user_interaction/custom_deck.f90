MODULE Custom_deck
  USE shared_data
  USE strings
  USE strings_advanced
  USE shunt
  USE evaluator
  USE shared_parser_data
  IMPLICIT NONE

CONTAINS


  !-----------------------------------------------------------------------------
  !These functions contain the user input deck elements
  !-----------------------------------------------------------------------------

  FUNCTION HandleCustomBlock(blockname,Element,Value)

    CHARACTER(len=entrylength),INTENT(IN)::blockname,Element,Value
    INTEGER :: HandleCustomBlock

    !The following line must always be present
    HandleCustomBlock=ERR_UNKNOWN_BLOCK

  END FUNCTION HandleCustomBlock


  FUNCTION CheckCustomBlocks()

    INTEGER :: CheckCustomBlocks

    !This subroutine is to allow you to force the code to bomb out if an essential element
    !Of the input deck is missing. If you either don't want to check ,are not extending the
    !Input deck, or all elements are set then set "CheckCustomBlocks = ERR_NONE". Otherwise
    !Set the return value to "CheckCustomBlocks = ERR_MISSING_ELEMENTS".

    CheckCustomBlocks=ERR_NONE

  END FUNCTION CheckCustomBlocks

END MODULE Custom_deck
