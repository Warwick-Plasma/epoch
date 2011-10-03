MODULE sdf_input

  USE mpi
  USE sdf_common
  USE sdf_input_ru
  USE sdf_input_r4
  USE sdf_input_r8

  IMPLICIT NONE

  INTERFACE sdf_read_header
    MODULE PROCEDURE &
        read_header_ru, &
        read_header_r4, &
        read_header_r8
  END INTERFACE sdf_read_header

  INTERFACE sdf_read_srl
    MODULE PROCEDURE &
        read_constant_real_r4, &
        read_constant_real_r8, &
        read_constant_integer, &
        read_constant_logical, &
        read_1d_array_real_r4, &
        read_2d_array_real_r4, &
        read_1d_array_real_r8, &
        read_2d_array_real_r8, &
        read_1d_array_integer, &
        read_2d_array_integer, &
        read_1d_array_logical, &
        read_2d_array_character
  END INTERFACE sdf_read_srl

END MODULE sdf_input
