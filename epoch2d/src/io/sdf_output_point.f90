MODULE sdf_output_point

  USE mpi
  USE sdf_common
  USE sdf_output
  USE sdf_output_point_ru
  USE sdf_output_point_r8

  IMPLICIT NONE

  INTERFACE sdf_write_srl_point_mesh
    MODULE PROCEDURE &
        sdf_write_srl_1d_pt_mesh_array, &
        sdf_write_srl_2d_pt_mesh_array, &
        sdf_write_srl_3d_pt_mesh_array, &
        sdf_write_srl_1d_pt_mesh_array4, &
        sdf_write_srl_2d_pt_mesh_array4, &
        sdf_write_srl_3d_pt_mesh_array4
  END INTERFACE sdf_write_srl_point_mesh

  INTERFACE sdf_write_srl_point_variable
    MODULE PROCEDURE &
        sdf_write_srl_pt_var_int_array, &
        sdf_write_srl_pt_var_int_array4, &
        sdf_write_srl_pt_var_flt_array, &
        sdf_write_srl_pt_var_flt_array4
  END INTERFACE sdf_write_srl_point_variable

END MODULE sdf_output_point
