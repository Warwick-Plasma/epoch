MODULE sdf_input_cartesian

  USE mpi
  USE sdf_common
  USE sdf_input
  USE sdf_input_cartesian_ru
  USE sdf_input_cartesian_r8

  IMPLICIT NONE

  INTERFACE sdf_read_srl_plain_mesh
    MODULE PROCEDURE &
        sdf_read_srl_1d_mesh, &
        sdf_read_srl_2d_mesh, &
        sdf_read_srl_3d_mesh
  END INTERFACE sdf_read_srl_plain_mesh

  INTERFACE sdf_read_plain_variable
    MODULE PROCEDURE &
        sdf_read_1d_float, &
        ! 2d_float and 3d_float share the same arguments as 1d_material
        ! and 2d_material so we need to disambiguate
        sdf_read_2d_variable, &
        sdf_read_3d_variable, &
        sdf_read_3d_material, &
        sdf_read_1d_integer, &
        sdf_read_2d_integer, &
        sdf_read_3d_integer, &
        sdf_read_1d_character, &
        sdf_read_2d_character, &
        sdf_read_3d_character
  END INTERFACE sdf_read_plain_variable

END MODULE sdf_input_cartesian
