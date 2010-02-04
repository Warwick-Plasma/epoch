MODULE input_functions

  USE shared_data
  USE iocommon

  IMPLICIT NONE

  SAVE

CONTAINS

  SUBROUTINE cfd_skip_block

    ! Minimal subroutine to skip past current block
    ! Assumes that the file is at a point where the block header has been read
    ! If it's at the start of a block call cfd_skip_block_header first
    current_displacement = block_header_end + block_length

  END SUBROUTINE cfd_skip_block



  SUBROUTINE cfd_skip_block_header

    ! Minimal subroutine used to skip block header
    current_displacement = block_header_end

  END SUBROUTINE cfd_skip_block_header



  SUBROUTINE cfd_skip_block_metadata

    ! Minimal subroutine to skip straight to the start of the real data in
    ! the current block
    current_displacement = block_header_end+block_md_length

  END SUBROUTINE cfd_skip_block_metadata

END MODULE input_functions
