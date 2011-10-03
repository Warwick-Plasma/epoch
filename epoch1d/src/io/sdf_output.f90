MODULE sdf_output

  USE mpi
  USE sdf_common
  USE sdf_output_ru
  USE sdf_output_r8

  IMPLICIT NONE

  INTERFACE sdf_write_srl
    MODULE PROCEDURE &
        sdf_write_constant_real, &
        sdf_write_constant_integer, &
        sdf_write_constant_logical, &
        sdf_write_1d_array_real, &
        sdf_write_2d_array_real, &
        sdf_write_2d_array_real_spec, &
        sdf_write_1d_array_integer, &
        sdf_write_2d_array_integer, &
        sdf_write_2d_array_integer_spec, &
        sdf_write_1d_array_logical, &
        sdf_write_2d_array_character
  END INTERFACE sdf_write_srl

END MODULE sdf_output
