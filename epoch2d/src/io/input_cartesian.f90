MODULE input_cartesian

  USE shared_data
  USE iocommon
  USE inputfunctions

  IMPLICIT NONE

  SAVE

CONTAINS

  !Grid loading functions
  SUBROUTINE cfd_Get_nD_Cartesian_Grid_MetaData_All(ndims,dims,extents)

    INTEGER,DIMENSION(:),INTENT(OUT) :: dims
    REAL(num),DIMENSION(:),INTENT(OUT) :: extents
    INTEGER,INTENT(IN) :: ndims
    !this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, dims, ndims, MPI_INTEGER, cfd_status, cfd_errcode)
    current_displacement = current_displacement + ndims * SoI 

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, ndims*2, mpireal, cfd_status, cfd_errcode)

    !After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    !Start of Data

    CALL cfd_Skip_Block_MetaData()

  END SUBROUTINE cfd_Get_nD_Cartesian_Grid_MetaData_All

  SUBROUTINE cfd_Get_1D_Cartesian_Grid_All(x)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: x
    INTEGER :: nx

    CALL cfd_Skip_Block_MetaData()
    nx=SIZE(x)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal,cfd_status,cfd_errcode)
    !That should be it, so now skip to end of block
    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_1D_Cartesian_Grid_All

  SUBROUTINE cfd_Get_2D_Cartesian_Grid_All(x,y)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: x,y
    INTEGER :: nx, ny

    CALL cfd_Skip_Block_MetaData()
    nx=SIZE(x)
    ny=SIZE(y)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, y, ny, mpireal,cfd_status,cfd_errcode)
    !That should be it, so now skip to end of block
    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_2D_Cartesian_Grid_All

  SUBROUTINE cfd_Get_3D_Cartesian_Grid_All(x,y,z)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: x,y,z
    INTEGER :: nx,ny,nz

    nx=SIZE(x)
    ny=SIZE(y)
    nz=SIZE(z)

    CALL cfd_Skip_Block_MetaData()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, x, nx, mpireal,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, y, ny, mpireal,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, z, nz, mpireal,cfd_status,cfd_errcode)
    !That should be it, so now skip to end of block
    CALL cfd_Skip_Block

  END SUBROUTINE cfd_Get_3D_Cartesian_Grid_All


  !Variable loading functions
  SUBROUTINE cfd_Get_nD_Cartesian_Variable_MetaData_All(ndims,dims,extents,stagger,meshname,meshclass)

    INTEGER,DIMENSION(:),INTENT(OUT) :: dims
    REAL(num),DIMENSION(:),INTENT(OUT) :: extents
    REAL(num),DIMENSION(:),INTENT(OUT) :: stagger
    INTEGER,INTENT(IN) :: ndims
    CHARACTER(len=*),INTENT(INOUT) :: meshname,meshclass
    CHARACTER(len=MaxStringLen) :: Meshname_file,meshclass_file
    INTEGER :: len_name,len_class

    len_name=LEN(meshname)
    len_class=LEN(meshclass)

    !this subroutine MUST be called after the call to Get_Common_Mesh_MetaData_All or it will break everything
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, dims, ndims, MPI_INTEGER, cfd_status, cfd_errcode)
    current_displacement = current_displacement + ndims*SoI

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    !Read grid stagger
    CALL MPI_FILE_READ_ALL(cfd_filehandle, stagger, ndims, mpireal, cfd_status, cfd_errcode)
    !Read data range
    CALL MPI_FILE_READ_ALL(cfd_filehandle, extents, 2, mpireal, cfd_status, cfd_errcode)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, Meshname_file, MaxStringLen, MPI_CHARACTER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, Meshclass_file, MaxStringLen, MPI_CHARACTER, cfd_status, cfd_errcode)

    meshname=meshname_file(1:MIN(MaxStringLen,len_name))
    meshclass=meshclass_file(1:MIN(MaxStringLen,len_class))

    !After this subroutine, all the metadata should be read in, so to make sure, just jump to known
    !start of Data

    CALL cfd_Skip_Block_MetaData()

  END SUBROUTINE cfd_Get_nD_Cartesian_Variable_MetaData_All

  SUBROUTINE cfd_Get_1D_Cartesian_Variable_Parallel(variable,subtype)

    REAL(num),INTENT(IN),DIMENSION(:) :: Variable
    INTEGER,INTENT(IN) :: subtype

    CALL cfd_Skip_Block_MetaData()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_1D_Cartesian_Variable_Parallel

  SUBROUTINE cfd_Get_1D_Cartesian_Variable_All(variable)

    REAL(num),INTENT(IN),DIMENSION(:) :: Variable

    CALL cfd_Skip_Block_MetaData()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_1D_Cartesian_Variable_All

  SUBROUTINE cfd_Get_2D_Cartesian_Variable_Parallel(variable,subtype)

    REAL(num),INTENT(IN),DIMENSION(:,:) :: Variable
    INTEGER,INTENT(IN) :: subtype

    CALL cfd_Skip_Block_MetaData()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_2D_Cartesian_Variable_Parallel

  SUBROUTINE cfd_Get_2D_Cartesian_Variable_All(variable)

    REAL(num),INTENT(IN),DIMENSION(:,:) :: Variable

    CALL cfd_Skip_Block_MetaData()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_2D_Cartesian_Variable_All

  SUBROUTINE cfd_Get_3D_Cartesian_Variable_Parallel(variable,subtype)

    REAL(num),INTENT(IN),DIMENSION(:,:,:) :: Variable
    INTEGER,INTENT(IN) :: subtype

    CALL cfd_Skip_Block_MetaData()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, subtype,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_3D_Cartesian_Variable_Parallel

  SUBROUTINE cfd_Get_3D_Cartesian_Variable_All(variable)

    REAL(num),INTENT(IN),DIMENSION(:,:,:) :: Variable

    CALL cfd_Skip_Block_MetaData()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, variable, SIZE(variable), mpireal, cfd_status, cfd_errcode)

    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_3D_Cartesian_Variable_All

END MODULE input_cartesian
