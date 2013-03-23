MODULE sdf_input

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
        read_3d_array_real_r4, &
        read_1d_array_real_r8, &
        read_2d_array_real_r8, &
        read_3d_array_real_r8, &
        read_1d_array_integer, &
        read_2d_array_integer, &
        read_1d_array_logical, &
        read_2d_array_character
  END INTERFACE sdf_read_srl

  INTERFACE sdf_read_srl_cpu_split
    MODULE PROCEDURE &
        read_srl_cpu_split, &
        read_srl_cpu_split_part
  END INTERFACE sdf_read_srl_cpu_split

  INTERFACE sdf_read_run_info
    MODULE PROCEDURE &
        read_run_info, &
        read_run_info_old, &
        read_run_info_minor
  END INTERFACE sdf_read_run_info

END MODULE sdf_input
