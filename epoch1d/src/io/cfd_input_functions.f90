MODULE cfd_input_functions

  USE cfd_common

  IMPLICIT NONE

  SAVE

CONTAINS

  SUBROUTINE cfd_skip_block(h)

    TYPE(cfd_file_handle) :: h

    ! Minimal subroutine to skip past current block
    ! Assumes that the file is at a point where the block header has been read
    ! If it's at the start of a block call cfd_skip_block_header first
    h%current_displacement = h%block_header_end + h%block_length

  END SUBROUTINE cfd_skip_block



  SUBROUTINE cfd_skip_block_header(h)

    TYPE(cfd_file_handle) :: h

    ! Minimal subroutine used to skip block header
    h%current_displacement = h%block_header_end

  END SUBROUTINE cfd_skip_block_header



  SUBROUTINE cfd_skip_block_metadata(h)

    TYPE(cfd_file_handle) :: h

    ! Minimal subroutine to skip straight to the start of the real data in
    ! the current block
    h%current_displacement = h%block_header_end + h%block_md_length

  END SUBROUTINE cfd_skip_block_metadata

END MODULE cfd_input_functions
