MODULE cfd_input_functions

  USE cfd_common

  IMPLICIT NONE

  SAVE

CONTAINS

  SUBROUTINE cfd_skip_block(h)

    TYPE(cfd_file_handle) :: h

    ! Minimal subroutine to skip to the start of the next block
    h%current_displacement = &
        h%block_start + h%block_header_size + h%block_length

  END SUBROUTINE cfd_skip_block



  SUBROUTINE cfd_skip_block_header(h)

    TYPE(cfd_file_handle) :: h

    ! Minimal subroutine to skip past the block header in
    ! the current block
    h%current_displacement = h%block_start + h%block_header_size

  END SUBROUTINE cfd_skip_block_header



  SUBROUTINE cfd_skip_block_metadata(h)

    TYPE(cfd_file_handle) :: h

    ! Minimal subroutine to skip to the start of the real data in
    ! the current block
    h%current_displacement = &
        h%block_start + h%block_header_size + h%block_md_length

  END SUBROUTINE cfd_skip_block_metadata

END MODULE cfd_input_functions
