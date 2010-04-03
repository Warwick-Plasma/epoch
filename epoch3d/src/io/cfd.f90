MODULE cfd

  USE cfd_common
  USE cfd_control
  USE cfd_input
  USE cfd_input_cartesian
  USE cfd_input_particle
  USE cfd_input_functions
  USE cfd_output
  USE cfd_output_cartesian
  USE cfd_output_particle
  USE mpi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cfd_file_handle

  PUBLIC :: c_cfd_read
  PUBLIC :: c_cfd_write

  PUBLIC :: c_type_scribble
  PUBLIC :: c_type_additional
  PUBLIC :: c_type_mesh
  PUBLIC :: c_type_mesh_variable
  PUBLIC :: c_type_snapshot
  PUBLIC :: c_type_stitched_vector
  PUBLIC :: c_type_stitched_magnitude
  PUBLIC :: c_type_constant
  PUBLIC :: c_type_arb_db
  PUBLIC :: c_type_integerarray
  PUBLIC :: c_type_info

  PUBLIC :: c_mesh_cartesian
  PUBLIC :: c_mesh_particle
  PUBLIC :: c_particle_cartesian
  PUBLIC :: c_particle_polar
  PUBLIC :: c_particle_cylindrical
  PUBLIC :: c_var_cartesian
  PUBLIC :: c_var_particle

  PUBLIC :: c_dimension_irrelevant
  PUBLIC :: c_dimension_1d
  PUBLIC :: c_dimension_2d
  PUBLIC :: c_dimension_3d

  PUBLIC :: cfd_open
  PUBLIC :: cfd_close
  PUBLIC :: cfd_set_max_string_length
  PUBLIC :: cfd_set_default_rank
  PUBLIC :: cfd_get_nblocks
  PUBLIC :: cfd_get_jobid
  PUBLIC :: cfd_get_next_block_info_all
  PUBLIC :: cfd_get_common_meshtype_metadata_all
  PUBLIC :: cfd_get_snapshot
  PUBLIC :: cfd_get_real_constant
  PUBLIC :: cfd_get_1d_cartesian_grid_all
  PUBLIC :: cfd_get_2d_cartesian_grid_all
  PUBLIC :: cfd_get_3d_cartesian_grid_all
  PUBLIC :: cfd_get_nd_cartesian_grid_metadata_all
  PUBLIC :: cfd_get_nd_cartesian_variable_metadata_all
  PUBLIC :: cfd_get_1d_cartesian_variable_parallel
  PUBLIC :: cfd_get_2d_cartesian_variable_parallel
  PUBLIC :: cfd_get_3d_cartesian_variable_parallel
  PUBLIC :: cfd_skip_block
  PUBLIC :: cfd_get_nd_particle_grid_metadata_all
  PUBLIC :: cfd_get_nd_particle_grid_parallel_with_iterator
  PUBLIC :: cfd_get_nd_particle_grid_parallel
  PUBLIC :: cfd_get_nd_particle_variable_metadata_all
  PUBLIC :: cfd_get_nd_particle_variable_parallel
  PUBLIC :: cfd_get_nd_particle_variable_parallel_with_iterator
  PUBLIC :: cfd_write_job_info
  PUBLIC :: cfd_write_source_code
  PUBLIC :: cfd_write_1d_cartesian_grid
  PUBLIC :: cfd_write_2d_cartesian_grid
  PUBLIC :: cfd_write_3d_cartesian_grid
  PUBLIC :: cfd_write_1d_cartesian_variable_parallel
  PUBLIC :: cfd_write_2d_cartesian_variable_parallel
  PUBLIC :: cfd_write_3d_cartesian_variable_parallel
  PUBLIC :: cfd_write_nd_particle_grid_with_iterator_all
  PUBLIC :: cfd_write_nd_particle_variable_with_iterator_all
  PUBLIC :: cfd_write_stitched_vector
  PUBLIC :: cfd_write_stitched_magnitude
  PUBLIC :: cfd_write_real_constant
  !PUBLIC :: cfd_open_read
  !PUBLIC :: cfd_skip_block_header
  !PUBLIC :: cfd_skip_block_metadata
  !PUBLIC :: cfd_open_clobber
  !PUBLIC :: cfd_write_block_header
  !PUBLIC :: cfd_write_meshtype_header
  !PUBLIC :: cfd_safe_write_string

END MODULE cfd
