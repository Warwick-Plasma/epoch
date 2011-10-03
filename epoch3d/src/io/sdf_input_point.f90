MODULE sdf_input_point

  USE mpi
  USE sdf_common
  USE sdf_input
  USE sdf_input_point_ru
  USE sdf_input_point_r8

  IMPLICIT NONE

  INTERFACE sdf_read_srl_point_mesh
    MODULE PROCEDURE &
        sdf_read_srl_1d_pt_mesh_array, &
        sdf_read_srl_2d_pt_mesh_array, &
        sdf_read_srl_3d_pt_mesh_array
  END INTERFACE sdf_read_srl_point_mesh

  INTERFACE sdf_read_srl_point_variable
    MODULE PROCEDURE &
        sdf_read_srl_pt_var_int_array, &
        sdf_read_srl_pt_var_flt_array
  END INTERFACE sdf_read_srl_point_variable

END MODULE sdf_input_point
