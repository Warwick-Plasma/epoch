MODULE sdf_input_point

  USE sdf_input_point_r4
  USE sdf_input_point_r8

  IMPLICIT NONE

  INTERFACE sdf_read_point_mesh_info
    MODULE PROCEDURE &
        read_point_mesh_info_ru, &
        read_point_mesh_info_i4_ru, &
        read_point_mesh_info_r4, &
        read_point_mesh_info_r8
  END INTERFACE sdf_read_point_mesh_info

  INTERFACE sdf_read_point_variable_info
    MODULE PROCEDURE &
        read_point_variable_info_ru, &
        read_point_variable_info_r4, &
        read_point_variable_info_r8
  END INTERFACE sdf_read_point_variable_info

  INTERFACE sdf_read_point_mesh
    MODULE PROCEDURE &
        read_point_mesh_r4, &
        read_point_mesh_r8
  END INTERFACE sdf_read_point_mesh

  INTERFACE sdf_read_point_variable
    MODULE PROCEDURE &
        read_point_variable_r4, &
        read_point_variable_r8, &
        read_point_variable_i4, &
        read_point_variable_i8
  END INTERFACE sdf_read_point_variable

  INTERFACE sdf_read_srl_point_mesh
    MODULE PROCEDURE &
        read_srl_1d_pt_mesh_array_r4, &
        read_srl_2d_pt_mesh_array_r4, &
        read_srl_3d_pt_mesh_array_r4, &
        read_srl_1d_pt_mesh_array_r8, &
        read_srl_2d_pt_mesh_array_r8, &
        read_srl_3d_pt_mesh_array_r8
  END INTERFACE sdf_read_srl_point_mesh

  INTERFACE sdf_read_srl_point_variable
    MODULE PROCEDURE &
        read_srl_pt_var_flt_array_r4, &
        read_srl_pt_var_flt_array_r8, &
        read_srl_pt_var_int_array, &
        read_srl_pt_var_logical_array
  END INTERFACE sdf_read_srl_point_variable

END MODULE sdf_input_point
