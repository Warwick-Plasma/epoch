MODULE iocommon

  USE shared_data
  USE job_info

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: c_cfd_read = 0
  INTEGER, PARAMETER :: c_cfd_write = 1

  INTEGER(4), PARAMETER :: c_type_scribble = -1
  INTEGER(4), PARAMETER :: c_type_additional = 0
  INTEGER(4), PARAMETER :: c_type_mesh = 1
  INTEGER(4), PARAMETER :: c_type_mesh_variable = 2
  INTEGER(4), PARAMETER :: c_type_snapshot = 3
  INTEGER(4), PARAMETER :: c_type_stitched_vector = 4
  INTEGER(4), PARAMETER :: c_type_stitched_magnitude = 5
  INTEGER(4), PARAMETER :: c_type_constant = 6
  INTEGER(4), PARAMETER :: c_type_arb_db = 7
  INTEGER(4), PARAMETER :: c_type_integerarray = 8
  INTEGER(4), PARAMETER :: c_type_info = 9

  INTEGER(4), PARAMETER :: c_mesh_cartesian = 0
  INTEGER(4), PARAMETER :: c_mesh_particle = 1
  INTEGER(4), PARAMETER :: c_particle_cartesian = 0
  INTEGER(4), PARAMETER :: c_particle_polar = 1
  INTEGER(4), PARAMETER :: c_particle_cylindrical = 2
  INTEGER(4), PARAMETER :: c_var_cartesian = 0
  INTEGER(4), PARAMETER :: c_var_particle = 1

  ! c_dimension_irrelevant is used where the dimensionality isn't needed, as
  ! with particle variables still keep dimensionality as a common quantity
  ! because other than this, they really are very alike
  INTEGER(4), PARAMETER :: c_dimension_irrelevant = 0
  INTEGER(4), PARAMETER :: c_dimension_1d = 1
  INTEGER(4), PARAMETER :: c_dimension_2d = 2
  INTEGER(4), PARAMETER :: c_dimension_3d = 3

  INTEGER(KIND=MPI_OFFSET_KIND) :: current_displacement
  INTEGER :: cfd_filehandle = -1, cfd_rank, cfd_comm
  INTEGER(4) :: cfd_nblocks
  INTEGER(4), PARAMETER :: cfd_version = 1, cfd_revision = 1

  INTEGER(4) :: max_string_len = 60
  INTEGER :: default_rank = 0

  INTEGER, PARAMETER :: header_offset_this_version = 10 * 4 + 8 + 3
  INTEGER, PARAMETER :: nblocks_offset_this_version = 5 * 4 + 3

  ! This cannot be changed without a major revision
  ! If you want to add more to every meshtype, tough luck
  ! You'll either have to tag it to every class or
  ! submit it for inclusion in the next major revision
  ! (This shouldn't ever happen, meshtype covers too many things,
  ! The only thing in common is that they include spatial information)
  INTEGER, PARAMETER :: meshtype_header_offset = 3 * 4

  INTEGER(4), PARAMETER :: soi  = 4 ! Size of integer
  INTEGER(4), PARAMETER :: soi8 = 8 ! Size of long (normally 8 byte integer)

  INTEGER(4) :: block_header_size, header_offset

  INTEGER :: cfd_errcode, cfd_status(MPI_STATUS_SIZE), cfd_mode
  LOGICAL :: cfd_writing

  ! Current block info
  INTEGER(8) :: block_length, block_md_length
  INTEGER(KIND=MPI_OFFSET_KIND) :: block_header_start, block_header_end

  TYPE(jobid_type) :: cfd_jobid

END MODULE iocommon
