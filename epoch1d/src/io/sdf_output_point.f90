MODULE sdf_output_point

  USE sdf_output_point_r4
  USE sdf_output_point_r8

  IMPLICIT NONE

  INTERFACE sdf_write_point_mesh
    MODULE PROCEDURE &
        write_point_mesh_r4, &
        write_point_mesh_r8, &
        write_point_mesh_gen_r4, &
        write_point_mesh_gen_r8, &
        write_nospec_point_mesh_r4, &
        write_nospec_point_mesh_r8, &
        write_nospec_point_mesh_gen_r4, &
        write_nospec_point_mesh_gen_r8
  END INTERFACE sdf_write_point_mesh

  INTERFACE sdf_write_point_variable
    MODULE PROCEDURE &
        write_point_variable_r4, &
        write_point_variable_r8, &
        write_point_variable_i4, &
        write_point_variable_i8, &
        write_point_variable_gen_r4, &
        write_point_variable_gen_r8, &
        write_point_variable_gen_i4, &
        write_point_variable_gen_i8, &
        write_nospec_point_variable_r4, &
        write_nospec_point_variable_r8, &
        write_nospec_point_variable_i4, &
        write_nospec_point_variable_i8, &
        write_nospec_point_variable_gen_r4, &
        write_nospec_point_variable_gen_r8, &
        write_nospec_point_variable_gen_i4, &
        write_nospec_point_variable_gen_i8
  END INTERFACE sdf_write_point_variable

  INTERFACE sdf_write_srl_point_mesh
    MODULE PROCEDURE &
        write_srl_1d_pt_mesh_i4_r4, &
        write_srl_2d_pt_mesh_i4_r4, &
        write_srl_3d_pt_mesh_i4_r4, &
        write_srl_1d_pt_mesh_i4_r8, &
        write_srl_2d_pt_mesh_i4_r8, &
        write_srl_3d_pt_mesh_i4_r8, &
        write_srl_1d_pt_mesh_i8_r4, &
        write_srl_2d_pt_mesh_i8_r4, &
        write_srl_3d_pt_mesh_i8_r4, &
        write_srl_1d_pt_mesh_i8_r8, &
        write_srl_2d_pt_mesh_i8_r8, &
        write_srl_3d_pt_mesh_i8_r8, &
        write_nospec_srl_1d_pt_mesh_i4_r4, &
        write_nospec_srl_2d_pt_mesh_i4_r4, &
        write_nospec_srl_3d_pt_mesh_i4_r4, &
        write_nospec_srl_1d_pt_mesh_i4_r8, &
        write_nospec_srl_2d_pt_mesh_i4_r8, &
        write_nospec_srl_3d_pt_mesh_i4_r8, &
        write_nospec_srl_1d_pt_mesh_i8_r4, &
        write_nospec_srl_2d_pt_mesh_i8_r4, &
        write_nospec_srl_3d_pt_mesh_i8_r4, &
        write_nospec_srl_1d_pt_mesh_i8_r8, &
        write_nospec_srl_2d_pt_mesh_i8_r8, &
        write_nospec_srl_3d_pt_mesh_i8_r8
  END INTERFACE sdf_write_srl_point_mesh

  INTERFACE sdf_write_srl_point_variable
    MODULE PROCEDURE &
        write_srl_pt_var_flt_i4_r4, &
        write_srl_pt_var_flt_i4_r8, &
        write_srl_pt_var_flt_i8_r4, &
        write_srl_pt_var_flt_i8_r8, &
        write_srl_pt_var_int_i4_r4, &
        write_srl_pt_var_int_i4_r8, &
        write_srl_pt_var_int_i8_r4, &
        write_srl_pt_var_int_i8_r8, &
        write_srl_pt_var_logical_i4_r4, &
        write_srl_pt_var_logical_i4_r8, &
        write_srl_pt_var_logical_i8_r4, &
        write_srl_pt_var_logical_i8_r8, &
        write_nospec_srl_pt_var_flt_i4_r4, &
        write_nospec_srl_pt_var_flt_i4_r8, &
        write_nospec_srl_pt_var_flt_i8_r4, &
        write_nospec_srl_pt_var_flt_i8_r8, &
        write_nospec_srl_pt_var_int_i4_r4, &
        write_nospec_srl_pt_var_int_i4_r8, &
        write_nospec_srl_pt_var_int_i8_r4, &
        write_nospec_srl_pt_var_int_i8_r8, &
        write_nospec_srl_pt_var_logical_i4_r4, &
        write_nospec_srl_pt_var_logical_i4_r8, &
        write_nospec_srl_pt_var_logical_i8_r4, &
        write_nospec_srl_pt_var_logical_i8_r8
  END INTERFACE sdf_write_srl_point_variable

END MODULE sdf_output_point
