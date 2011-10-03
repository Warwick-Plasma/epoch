MODULE sdf_output_cartesian

  USE mpi
  USE sdf_common
  USE sdf_output
  USE sdf_output_cartesian_ru
  USE sdf_output_cartesian_r8

  IMPLICIT NONE

  INTERFACE sdf_write_srl_plain_mesh
    MODULE PROCEDURE &
        write_srl_1d_mesh_r8, &
        write_srl_2d_mesh_r8, &
        write_srl_3d_mesh_r8
  END INTERFACE sdf_write_srl_plain_mesh

  INTERFACE sdf_write_plain_mesh
    MODULE PROCEDURE &
        write_1d_mesh_r8, &
        write_2d_mesh_r8, &
        write_3d_mesh_r8
  END INTERFACE sdf_write_plain_mesh

  INTERFACE sdf_write_plain_variable
    MODULE PROCEDURE &
        write_1d_float_r8, &
        write_2d_float_r8, &
        write_3d_float_r8, &
        write_1d_character_r8, &
        write_2d_character_r8, &
        write_3d_character_r8, &
        write_1d_integer_r8, &
        write_2d_integer_r8, &
        write_3d_integer_r8
  END INTERFACE sdf_write_plain_variable

  INTERFACE sdf_write_material
    MODULE PROCEDURE &
        write_1d_material_r8, &
        write_2d_material_r8, &
        write_3d_material_r8
  END INTERFACE sdf_write_material

  INTERFACE sdf_write_matvar
    MODULE PROCEDURE &
        write_1d_matvar_r8, &
        write_2d_matvar_r8, &
        write_3d_matvar_r8
  END INTERFACE sdf_write_matvar

  INTERFACE sdf_write_species
    MODULE PROCEDURE &
        write_1d_species_r8, &
        write_2d_species_r8, &
        write_3d_species_r8
  END INTERFACE sdf_write_species

END MODULE sdf_output_cartesian
