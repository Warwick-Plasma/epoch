MODULE iocommon

  USE shared_data

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: TYPE_SCRIBBLE=-1, TYPE_ADDITIONAL=0, TYPE_MESH=1, TYPE_MESH_VARIABLE=2, TYPE_SNAPSHOT=3
  INTEGER, PARAMETER :: TYPE_STITCHED_VECTOR=4, TYPE_STITCHED_MAGNITUDE=5, TYPE_CONSTANT=6, TYPE_ARB_DB=7
  INTEGER, PARAMETER :: TYPE_INTEGERARRAY=8

  INTEGER, PARAMETER :: MESH_CARTESIAN=0, MESH_PARTICLE=1
  INTEGER, PARAMETER :: PARTICLE_CARTESIAN=0, PARTICLE_POLAR=1 ,PARTICLE_CYLINDRICAL=2
  INTEGER, PARAMETER :: VAR_CARTESIAN=0, VAR_PARTICLE=1

  !Dimension_Irrelevant is used where the dimensionality isn't needed, as with particle variables
  !Still keep dimensionality as a common quantity because other than this, they really are very
  !Alike
  INTEGER, PARAMETER :: DIMENSION_IRRELEVANT=0, DIMENSION_1D=1, DIMENSION_2D=2, DIMENSION_3D=3

  INTEGER(KIND=MPI_OFFSET_KIND) :: current_displacement
  INTEGER :: cfd_filehandle=-1,cfd_rank,cfd_comm,nblocks
  INTEGER,PARAMETER :: CFD_Version=1, CFD_Revision=0

  INTEGER :: MaxStringLen = 40, default_rank=0

  INTEGER,PARAMETER :: header_offset_this_version = 6*4+3
  !This cannot be changed without a major revision
  !If you want to add more to every meshtype, tough luck
  !You'll either have to tag it to every class or 
  !Submit it for inclusion in the next major revision
  !(This shouldn't ever happen, meshtype covers too many things,
  !The only thing in common is that they include spatial information)
  INTEGER,PARAMETER :: MeshType_Header_offset = 3 * 4

  INTEGER,PARAMETER :: SoI  = 4 !Size of integer
  INTEGER,PARAMETER :: SoI8 = 8 !Size of long (normally 8 byte integer)

  INTEGER  :: block_header_size, header_offset

  INTEGER :: cfd_errcode,cfd_status(MPI_STATUS_SIZE),cfd_mode
  LOGICAL :: cfd_writing, cfd_reading

  !Current block info
  INTEGER :: c_block_Type
  INTEGER(kind=8) :: block_length,block_md_length
  INTEGER(kind=MPI_OFFSET_KIND) :: block_header_start,block_header_end 


END MODULE iocommon
