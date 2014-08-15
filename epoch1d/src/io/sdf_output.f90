MODULE sdf_output

  USE sdf_output_ru
  USE sdf_output_r4
  USE sdf_output_r8

  IMPLICIT NONE

  INTERFACE sdf_write_header
    MODULE PROCEDURE &
        write_header_r4, &
        write_header_r8
  END INTERFACE sdf_write_header

  INTERFACE sdf_write_srl
    MODULE PROCEDURE &
        write_constant_real_r4, &
        write_constant_real_r8, &
        write_constant_integer, &
        write_constant_logical, &
        write_1d_array_real_r4, &
        write_1d_array_real_spec_r4, &
        write_2d_array_real_r4, &
        write_2d_array_real_spec_r4, &
        write_1d_array_real_r8, &
        write_1d_array_real_spec_r8, &
        write_2d_array_real_r8, &
        write_2d_array_real_spec_r8, &
        write_1d_array_integer, &
        write_1d_array_integer_spec, &
        write_2d_array_integer, &
        write_2d_array_integer_spec, &
        write_1d_array_logical, &
        write_1d_array_logical_spec, &
        write_2d_array_character
  END INTERFACE sdf_write_srl

  INTERFACE sdf_write_cpu_split
    MODULE PROCEDURE &
        write_cpu_split_1d, &
        write_cpu_split_1d_spec, &
        write_cpu_split_mix, &
        write_cpu_split_3d, &
        write_cpu_split_part
  END INTERFACE sdf_write_cpu_split

  INTERFACE sdf_write_run_info
    MODULE PROCEDURE &
        write_run_info_old, &
        write_run_info_minor
  END INTERFACE sdf_write_run_info

  INTERFACE sdf_write_datablock
    MODULE PROCEDURE &
        write_datablock_bytearray, &
        write_datablock_chararray
  END INTERFACE sdf_write_datablock

  INTERFACE sdf_write_namevalue
    MODULE PROCEDURE &
        write_namevalue_i4, &
        write_namevalue_i8, &
        write_namevalue_r4, &
        write_namevalue_r8, &
        write_namevalue_logical, &
        write_namevalue_string
  END INTERFACE sdf_write_namevalue

END MODULE sdf_output
