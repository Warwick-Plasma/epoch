MODULE iocommon

  USE shared_data

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: c_type_scribble=-1, c_type_additional=0, c_type_mesh=1, c_type_mesh_variable=2, c_type_snapshot=3
  INTEGER, PARAMETER :: c_type_stitched_vector=4, c_type_stitched_magnitude=5, c_type_constant=6, c_type_arb_db=7
  INTEGER, PARAMETER :: c_type_integerarray=8

  INTEGER, PARAMETER :: c_mesh_cartesian=0, c_mesh_particle=1
  INTEGER, PARAMETER :: c_particle_cartesian=0, c_particle_polar=1 ,c_particle_cylindrical=2
  INTEGER, PARAMETER :: c_var_cartesian=0, c_var_particle=1

  !c_dimension_irrelevant is used where the dimensionality isn't needed, as with particle variables
  !Still keep dimensionality as a common quantity because other than this, they really are very
  !Alike
  INTEGER, PARAMETER :: c_dimension_irrelevant=0, c_dimension_1d=1, c_dimension_2d=2, c_dimension_3d=3

  INTEGER(KIND=MPI_OFFSET_KIND) :: current_displacement
  INTEGER :: cfd_filehandle=-1,cfd_rank,cfd_comm,nblocks
  INTEGER,PARAMETER :: cfd_version=1, cfd_revision=0

  INTEGER :: max_string_len = 60, default_rank=0

  INTEGER,PARAMETER :: header_offset_this_version = 6*4+3
  !This cannot be changed without a major revision
  !If you want to add more to every meshtype, tough luck
  !You'll either have to tag it to every class or 
  !Submit it for inclusion in the next major revision
  !(This shouldn't ever happen, meshtype covers too many things,
  !The only thing in common is that they include spatial information)
  INTEGER,PARAMETER :: meshtype_header_offset = 3 * 4

  INTEGER,PARAMETER :: soi  = 4 !Size of INTEGER
  INTEGER,PARAMETER :: soi8 = 8 !Size of long (normally 8 byte INTEGER)

  INTEGER  :: block_header_size, header_offset

  INTEGER :: cfd_errcode,cfd_status(MPI_STATUS_SIZE),cfd_mode
  LOGICAL :: cfd_writing, cfd_reading

  !current block info
  INTEGER(KIND=8) :: block_length,block_md_length
  INTEGER(KIND=MPI_OFFSET_KIND) :: block_header_start,block_header_end 


END MODULE iocommon
