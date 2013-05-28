MODULE sdf

  USE sdf_control
  USE sdf_input
  USE sdf_input_cartesian
  USE sdf_input_point
  USE sdf_input_station
  USE sdf_input_util
  USE sdf_output
  USE sdf_output_cartesian
  USE sdf_output_point
  USE sdf_output_station

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sdf_file_handle

  PUBLIC :: c_sdf_read
  PUBLIC :: c_sdf_write
  PUBLIC :: c_sdf_append
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
  PUBLIC :: c_blocktype_contiguous_tensor
  PUBLIC :: c_blocktype_contiguous_material
  PUBLIC :: c_blocktype_contiguous_matvar
  PUBLIC :: c_blocktype_contiguous_species
  PUBLIC :: c_blocktype_cpu_split
  PUBLIC :: c_blocktype_stitched_obstacle_group
  PUBLIC :: c_blocktype_unstructured_mesh
  PUBLIC :: c_blocktype_stitched
  PUBLIC :: c_blocktype_contiguous
  PUBLIC :: c_blocktype_lagrangian_mesh
  PUBLIC :: c_blocktype_station
  PUBLIC :: c_blocktype_max

  PUBLIC :: c_datatype_null
  PUBLIC :: c_datatype_integer4
  PUBLIC :: c_datatype_integer8
  PUBLIC :: c_datatype_real4
  PUBLIC :: c_datatype_real8
  PUBLIC :: c_datatype_real16
  PUBLIC :: c_datatype_character
  PUBLIC :: c_datatype_logical
  PUBLIC :: c_datatype_other
  PUBLIC :: c_datatype_max

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

  PUBLIC :: c_blocktypes_char
  PUBLIC :: c_datatypes_char
  PUBLIC :: c_type_sizes

  PUBLIC :: sdf_open
  PUBLIC :: sdf_close
  PUBLIC :: sdf_set_string_length
  PUBLIC :: sdf_set_default_rank
  PUBLIC :: sdf_set_point_array_size
  PUBLIC :: sdf_get_point_array_size
  PUBLIC :: sdf_seek_start
  PUBLIC :: sdf_seek_block
  PUBLIC :: sdf_find_block_by_id
  PUBLIC :: sdf_get_block_id
  PUBLIC :: sdf_get_data_location
  PUBLIC :: sdf_set_data_location
  PUBLIC :: sdf_get_all_stations
  PUBLIC :: sdf_station_seek_time
  PUBLIC :: sdf_flush
  PUBLIC :: sdf_read_header
  PUBLIC :: sdf_read_blocklist
  PUBLIC :: sdf_read_nblocks
  PUBLIC :: sdf_read_jobid
  PUBLIC :: sdf_read_block_header
  PUBLIC :: sdf_read_next_block_header
  PUBLIC :: sdf_read_srl
  PUBLIC :: sdf_read_array_info
  PUBLIC :: sdf_read_material_info
  PUBLIC :: sdf_read_stitched_info
  PUBLIC :: sdf_read_plain_mesh
  PUBLIC :: sdf_read_plain_mesh_info
  PUBLIC :: sdf_read_obstacle_group_info
  PUBLIC :: sdf_read_plain_variable
  PUBLIC :: sdf_read_plain_variable_info
  PUBLIC :: sdf_read_point_mesh
  PUBLIC :: sdf_read_point_mesh_info
  PUBLIC :: sdf_read_point_variable
  PUBLIC :: sdf_read_point_variable_info
  PUBLIC :: sdf_read_lagrangian_mesh
  PUBLIC :: sdf_read_srl_plain_mesh
  PUBLIC :: sdf_read_srl_point_mesh
  PUBLIC :: sdf_read_srl_point_variable
  PUBLIC :: sdf_read_cpu_split_info
  PUBLIC :: sdf_read_srl_cpu_split
  PUBLIC :: sdf_read_run_info
  PUBLIC :: sdf_read_station_info
  PUBLIC :: sdf_read_station_info_arrays
  PUBLIC :: sdf_read_station_info_arrays_all
  PUBLIC :: sdf_read_station_array
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
  PUBLIC :: sdf_write_lagrangian_mesh
  PUBLIC :: sdf_write_srl_plain_mesh
  PUBLIC :: sdf_write_srl_point_mesh
  PUBLIC :: sdf_write_srl_point_variable
  PUBLIC :: sdf_write_station_header
  PUBLIC :: sdf_write_station_array
  PUBLIC :: sdf_write_station_material
  PUBLIC :: sdf_write_station_matvar
  PUBLIC :: sdf_write_station_species
  PUBLIC :: sdf_write_stitched
  PUBLIC :: sdf_write_stitched_var
  PUBLIC :: sdf_write_stitched_tensor
  PUBLIC :: sdf_write_stitched_tensor_mat
  PUBLIC :: sdf_write_stitched_material
  PUBLIC :: sdf_write_stitched_matvar
  PUBLIC :: sdf_write_stitched_species
  PUBLIC :: sdf_write_stitched_obstacle_group
  PUBLIC :: sdf_write_cpu_split
  PUBLIC :: sdf_errorcode

END MODULE sdf
