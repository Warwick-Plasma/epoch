MODULE sdf

  USE mpi
  USE sdf_common
  USE sdf_control
  USE sdf_input
  USE sdf_input_cartesian
  USE sdf_input_point
  USE sdf_input_util
  USE sdf_output
  USE sdf_output_cartesian
  USE sdf_output_point
  USE sdf_output_util

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sdf_file_handle

  PUBLIC :: c_sdf_read
  PUBLIC :: c_sdf_write
  PUBLIC :: c_id_length

  PUBLIC :: c_blocktype_scrubbed
  PUBLIC :: c_blocktype_null
  PUBLIC :: c_blocktype_plain_mesh
  PUBLIC :: c_blocktype_point_mesh
  PUBLIC :: c_blocktype_plain_variable
  PUBLIC :: c_blocktype_point_variable
  PUBLIC :: c_blocktype_constant
  PUBLIC :: c_blocktype_array
  PUBLIC :: c_blocktype_run_info
  PUBLIC :: c_blocktype_source
  PUBLIC :: c_blocktype_stitched_tensor
  PUBLIC :: c_blocktype_stitched_material
  PUBLIC :: c_blocktype_stitched_matvar
  PUBLIC :: c_blocktype_stitched_species
  PUBLIC :: c_blocktype_species
  PUBLIC :: c_blocktype_plain_derived
  PUBLIC :: c_blocktype_point_derived
  PUBLIC :: c_blocktype_multi_tensor
  PUBLIC :: c_blocktype_multi_material
  PUBLIC :: c_blocktype_multi_matvar
  PUBLIC :: c_blocktype_multi_species
  PUBLIC :: c_blocktype_cpu_split

  PUBLIC :: c_datatype_null
  PUBLIC :: c_datatype_integer4
  PUBLIC :: c_datatype_integer8
  PUBLIC :: c_datatype_real4
  PUBLIC :: c_datatype_real8
  PUBLIC :: c_datatype_real16
  PUBLIC :: c_datatype_character
  PUBLIC :: c_datatype_logical
  PUBLIC :: c_datatype_other

  PUBLIC :: c_geometry_null
  PUBLIC :: c_geometry_cartesian
  PUBLIC :: c_geometry_cylindrical
  PUBLIC :: c_geometry_spherical

  PUBLIC :: c_dimension_irrelevant
  PUBLIC :: c_dimension_1d
  PUBLIC :: c_dimension_2d
  PUBLIC :: c_dimension_3d

  PUBLIC :: c_stagger_cell_centre
  PUBLIC :: c_stagger_face_x
  PUBLIC :: c_stagger_face_y
  PUBLIC :: c_stagger_face_z
  PUBLIC :: c_stagger_edge_x
  PUBLIC :: c_stagger_edge_y
  PUBLIC :: c_stagger_edge_z
  PUBLIC :: c_stagger_vertex

  PUBLIC :: sdf_open
  PUBLIC :: sdf_close
  PUBLIC :: sdf_set_string_length
  PUBLIC :: sdf_set_default_rank
  PUBLIC :: sdf_seek_start
  PUBLIC :: sdf_read_header
  PUBLIC :: sdf_read_blocklist
  PUBLIC :: sdf_read_nblocks
  PUBLIC :: sdf_read_jobid
  PUBLIC :: sdf_read_next_block_header
  PUBLIC :: sdf_read_srl
  PUBLIC :: sdf_read_array_info
  PUBLIC :: sdf_read_material_info
  PUBLIC :: sdf_read_plain_mesh_info
  PUBLIC :: sdf_read_plain_variable
  PUBLIC :: sdf_read_plain_variable_info
  PUBLIC :: sdf_read_point_mesh
  PUBLIC :: sdf_read_point_mesh_info
  PUBLIC :: sdf_read_point_variable
  PUBLIC :: sdf_read_point_variable_info
  PUBLIC :: sdf_read_srl_plain_mesh
  PUBLIC :: sdf_read_srl_point_mesh
  PUBLIC :: sdf_read_srl_point_variable
  PUBLIC :: sdf_write_header
  PUBLIC :: sdf_write_run_info
  PUBLIC :: sdf_write_source_code
  PUBLIC :: sdf_write_material
  PUBLIC :: sdf_write_matvar
  PUBLIC :: sdf_write_species
  PUBLIC :: sdf_write_srl
  PUBLIC :: sdf_write_plain_mesh
  PUBLIC :: sdf_write_plain_variable
  PUBLIC :: sdf_write_point_mesh
  PUBLIC :: sdf_write_point_variable
  PUBLIC :: sdf_write_srl_plain_mesh
  PUBLIC :: sdf_write_srl_point_mesh
  PUBLIC :: sdf_write_srl_point_variable
  PUBLIC :: sdf_write_stitched_tensor
  PUBLIC :: sdf_write_stitched_tensor_mat
  PUBLIC :: sdf_write_stitched_material
  PUBLIC :: sdf_write_stitched_matvar
  PUBLIC :: sdf_write_stitched_species
  PUBLIC :: sdf_write_multi_tensor
  PUBLIC :: sdf_write_multi_tensor_mat
  PUBLIC :: sdf_write_multi_material
  PUBLIC :: sdf_write_multi_matvar
  PUBLIC :: sdf_write_multi_species
  PUBLIC :: sdf_write_cpu_split

END MODULE sdf
