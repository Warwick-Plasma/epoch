MODULE sdf_input_cartesian

  USE mpi
  USE sdf_common
  USE sdf_input
  USE sdf_input_cartesian_ru
  USE sdf_input_cartesian_r4
  USE sdf_input_cartesian_r8

  IMPLICIT NONE

  INTERFACE sdf_read_plain_mesh_info
    MODULE PROCEDURE &
        read_plain_mesh_info_ru, &
        read_plain_mesh_info_r4, &
        read_plain_mesh_info_r8
  END INTERFACE sdf_read_plain_mesh_info

  INTERFACE sdf_read_plain_variable_info
    MODULE PROCEDURE &
        read_plain_variable_info_ru, &
        read_plain_variable_info_r4, &
        read_plain_variable_info_r8
  END INTERFACE sdf_read_plain_variable_info

  INTERFACE sdf_read_srl_plain_mesh
    MODULE PROCEDURE &
        read_srl_1d_mesh_r4, &
        read_srl_2d_mesh_r4, &
        read_srl_3d_mesh_r4, &
        read_srl_1d_mesh_r8, &
        read_srl_2d_mesh_r8, &
        read_srl_3d_mesh_r8
  END INTERFACE sdf_read_srl_plain_mesh

  INTERFACE sdf_read_plain_mesh
    MODULE PROCEDURE &
        read_1d_mesh_r4, &
        read_2d_mesh_r4, &
        read_3d_mesh_r4, &
        read_1d_mesh_r8, &
        read_2d_mesh_r8, &
        read_3d_mesh_r8
  END INTERFACE sdf_read_plain_mesh

  INTERFACE sdf_read_plain_variable
    MODULE PROCEDURE &
        read_1d_float_r4, &
        read_1d_float_r8, &
        ! 2d_float and 3d_float share the same arguments as 1d_material
        ! and 2d_material so we need to disambiguate
        read_2d_variable_r4, &
        read_3d_variable_r4, &
        read_3d_material_r4, &
        read_2d_variable_r8, &
        read_3d_variable_r8, &
        read_3d_material_r8, &
        read_1d_integer, &
        read_2d_integer, &
        read_3d_integer, &
        read_1d_character, &
        read_2d_character, &
        read_3d_character
  END INTERFACE sdf_read_plain_variable

END MODULE sdf_input_cartesian
