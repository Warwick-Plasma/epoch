MODULE sdf_output_cartesian

  USE mpi
  USE sdf_common
  USE sdf_output
  USE sdf_output_cartesian_ru
  USE sdf_output_cartesian_r8

  IMPLICIT NONE

  INTERFACE sdf_write_srl_plain_mesh
    MODULE PROCEDURE &
        sdf_write_srl_1d_mesh, &
        sdf_write_srl_2d_mesh, &
        sdf_write_srl_3d_mesh
  END INTERFACE sdf_write_srl_plain_mesh

  INTERFACE sdf_write_plain_mesh
    MODULE PROCEDURE &
        sdf_write_1d_mesh, &
        sdf_write_2d_mesh, &
        sdf_write_3d_mesh
  END INTERFACE sdf_write_plain_mesh

  INTERFACE sdf_write_plain_variable
    MODULE PROCEDURE &
        sdf_write_1d_float, &
        sdf_write_2d_float, &
        sdf_write_3d_float, &
        sdf_write_1d_character, &
        sdf_write_2d_character, &
        sdf_write_3d_character, &
        sdf_write_1d_integer, &
        sdf_write_2d_integer, &
        sdf_write_3d_integer
  END INTERFACE sdf_write_plain_variable

  INTERFACE sdf_write_material
    MODULE PROCEDURE &
        sdf_write_1d_material, &
        sdf_write_2d_material, &
        sdf_write_3d_material
  END INTERFACE sdf_write_material

  INTERFACE sdf_write_matvar
    MODULE PROCEDURE &
        sdf_write_1d_matvar, &
        sdf_write_2d_matvar, &
        sdf_write_3d_matvar
  END INTERFACE sdf_write_matvar

  INTERFACE sdf_write_species
    MODULE PROCEDURE &
        sdf_write_1d_species, &
        sdf_write_2d_species, &
        sdf_write_3d_species
  END INTERFACE sdf_write_species

END MODULE sdf_output_cartesian
