MODULE sdf_input_functions

  USE sdf_common

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_skip_block(h)

    TYPE(sdf_file_handle) :: h

    ! Minimal subroutine to skip to the start of the next block
    h%current_displacement = &
        h%block_start + h%block_header_size + h%block_length

  END SUBROUTINE sdf_skip_block



  SUBROUTINE sdf_skip_block_header(h)

    TYPE(sdf_file_handle) :: h

    ! Minimal subroutine to skip past the block header in
    ! the current block
    h%current_displacement = h%block_start + h%block_header_size

  END SUBROUTINE sdf_skip_block_header



  SUBROUTINE sdf_skip_block_info(h)

    TYPE(sdf_file_handle) :: h

    ! Minimal subroutine to skip to the start of the real data in
    ! the current block
    h%current_displacement = &
        h%block_start + h%block_header_size + h%block_md_length

  END SUBROUTINE sdf_skip_block_info

END MODULE sdf_input_functions
