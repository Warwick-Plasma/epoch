MODULE input_cartesian

  USE shared_data
  USE iocommon
  USE input_functions

  IMPLICIT NONE

  SAVE

CONTAINS

  !Grid loading functions
  SUBROUTINE cfd_get_nd_cartesian_grid_metadata_all(ndims,dims,extents)

    INTEGER,DIMENSION(:),INTENT(OUT) :: dims
    REAL(num),DIMENSION(:),INTENT(OUT) :: extents
    INTEGER,INTENT(IN) :: ndims
    !this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, dims, ndims, MPI_INTEGER, cfd_status, cfd_errcode)
    current_displacement = current_displacement + ndims * soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, ndims*2, mpireal, cfd_status, cfd_errcode)

    !After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    !start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_cartesian_grid_metadata_all

  SUBROUTINE cfd_get_1d_cartesian_grid_all(x)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: x
    INTEGER :: nx

    CALL cfd_skip_block_metadata()
    nx=SIZE(x)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal,cfd_status,cfd_errcode)
    !That should be it, so now skip to end of block
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_1d_cartesian_grid_all

  SUBROUTINE cfd_get_2d_cartesian_grid_all(x,y)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: x,y
    INTEGER :: nx, ny

    CALL cfd_skip_block_metadata()
    nx=SIZE(x)
    ny=SIZE(y)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, y, ny, mpireal,cfd_status,cfd_errcode)
    !That should be it, so now skip to end of block
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_2d_cartesian_grid_all

  SUBROUTINE cfd_get_3d_cartesian_grid_all(x,y,z)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: x,y,z
    INTEGER :: nx,ny,nz

    nx=SIZE(x)
    ny=SIZE(y)
    nz=SIZE(z)

    CALL cfd_skip_block_metadata()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, y, ny, mpireal,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, z, nz, mpireal,cfd_status,cfd_errcode)
    !That should be it, so now skip to end of block
    CALL cfd_skip_block

  END SUBROUTINE cfd_get_3d_cartesian_grid_all


  !variable loading functions
  SUBROUTINE cfd_get_nd_cartesian_variable_metadata_all(ndims,dims,extents,stagger,mesh_name,mesh_class)

    INTEGER,DIMENSION(:),INTENT(OUT) :: dims
    REAL(num),DIMENSION(:),INTENT(OUT) :: extents
    REAL(num),DIMENSION(:),INTENT(OUT) :: stagger
    INTEGER,INTENT(IN) :: ndims
    CHARACTER(len=*),INTENT(INOUT) :: mesh_name,mesh_class
    CHARACTER(len=max_string_len) :: mesh_name_file,mesh_class_file
    INTEGER :: len_name,len_class

    len_name=len(mesh_name)
    len_class=len(mesh_class)

    !this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, dims, ndims, MPI_INTEGER, cfd_status, cfd_errcode)
    current_displacement = current_displacement + ndims*soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    !Read grid stagger
    CALL MPI_FILE_READ_ALL(cfd_filehandle, stagger, ndims, mpireal, cfd_status, cfd_errcode)
    !Read data range
    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, 2, mpireal, cfd_status, cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh_name_file, max_string_len, MPI_CHARACTER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, mesh_class_file, max_string_len, MPI_CHARACTER, cfd_status, cfd_errcode)

    mesh_name=mesh_name_file(1:MIN(max_string_len,len_name))
    mesh_class=mesh_class_file(1:MIN(max_string_len,len_class))

    !After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    !start of data

    CALL cfd_skip_block_metadata()

  END SUBROUTINE cfd_get_nd_cartesian_variable_metadata_all

  SUBROUTINE cfd_get_1d_cartesian_variable_parallel(variable,subtype)

    REAL(num),INTENT(IN),DIMENSION(:) :: variable
    INTEGER,INTENT(IN) :: subtype

    CALL cfd_skip_block_metadata()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_1d_cartesian_variable_parallel

  SUBROUTINE cfd_get_1d_cartesian_variable_all(variable)

    REAL(num),INTENT(IN),DIMENSION(:) :: variable

    CALL cfd_skip_block_metadata()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_1d_cartesian_variable_all

  SUBROUTINE cfd_get_2d_cartesian_variable_parallel(variable,subtype)

    REAL(num),INTENT(IN),DIMENSION(:,:) :: variable
    INTEGER,INTENT(IN) :: subtype

    CALL cfd_skip_block_metadata()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_2d_cartesian_variable_parallel

  SUBROUTINE cfd_get_2d_cartesian_variable_all(variable)

    REAL(num),INTENT(IN),DIMENSION(:,:) :: variable

    CALL cfd_skip_block_metadata()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_2d_cartesian_variable_all

  SUBROUTINE cfd_get_3d_cartesian_variable_parallel(variable,subtype)

    REAL(num),INTENT(IN),DIMENSION(:,:,:) :: variable
    INTEGER,INTENT(IN) :: subtype

    CALL cfd_skip_block_metadata()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_3d_cartesian_variable_parallel

  SUBROUTINE cfd_get_3d_cartesian_variable_all(variable)

    REAL(num),INTENT(IN),DIMENSION(:,:,:) :: variable

    CALL cfd_skip_block_metadata()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_3d_cartesian_variable_all

END MODULE input_cartesian
