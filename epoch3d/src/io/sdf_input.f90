MODULE sdf_input

  USE mpi
  USE sdf_common
  USE sdf_input_ru
  USE sdf_input_r8

  IMPLICIT NONE

  INTERFACE sdf_read_srl
    MODULE PROCEDURE &
        sdf_read_constant_real, &
        sdf_read_constant_integer, &
        sdf_read_constant_logical, &
        sdf_read_1d_array_real, &
        sdf_read_2d_array_real, &
        sdf_read_1d_array_integer, &
        sdf_read_2d_array_integer, &
        sdf_read_1d_array_logical, &
        sdf_read_2d_array_character
  END INTERFACE sdf_read_srl

END MODULE sdf_input
