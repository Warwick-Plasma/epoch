
MODULE input

  USE shared_data
  USE iocommon
  USE inputfunctions
  USE input_cartesian

  IMPLICIT NONE

  SAVE

CONTAINS


  SUBROUTINE cfd_Open_Read(filename)

    CHARACTER(len=*),INTENT(IN) :: filename
    CHARACTER(len=3) :: CFD

    INTEGER:: FileVersion,FileRevision

    CALL MPI_BARRIER(cfd_comm,cfd_errcode)
    CALL MPI_FILE_OPEN(cfd_comm, TRIM(filename), cfd_mode, &
         MPI_INFO_NULL, cfd_filehandle, cfd_errcode)

    current_displacement=0
       CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
          "native", MPI_INFO_NULL, cfd_errcode)
       !Read the header
       CALL MPI_FILE_READ_ALL(cfd_filehandle, CFD, 3, MPI_CHARACTER, cfd_status, cfd_errcode)
       !If this isn't "CFD" then this isn't a CFD file
       IF (CFD /= "CFD") THEN
          CALL MPI_FILE_CLOSE(cfd_filehandle,cfd_errcode)
          PRINT *,"The specified file is not a valid CFD file"
          CALL MPI_ABORT(cfd_comm,cfd_errcode)
       ENDIF
       current_displacement=3
       CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
           "native", MPI_INFO_NULL, cfd_errcode)
       !Read in the basic file info. Should check version info, but this is version 1, so 
       !Let's not worry about it
       CALL MPI_FILE_READ_ALL(cfd_filehandle, header_offset, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, block_header_size, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, FileVersion, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, FileRevision, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, MaxStringLen, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, nBlocks, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    IF (FileVersion .GT. cfd_Version) THEN
       IF (rank == default_rank) PRINT *,"Version number incompatible"
       CALL MPI_ABORT(cfd_comm,cfd_errcode)
    ENDIF

    IF (FileRevision .GT. cfd_Revision) THEN
       IF (rank == default_rank) PRINT *,"Revision number of file is too high. Writing disabled"
       cfd_writing=.FALSE.
    ENDIF

    current_displacement = header_offset

  END SUBROUTINE cfd_Open_Read

  SUBROUTINE cfd_Get_Next_Block_Info_All(name,class,Type)

    CHARACTER(len=*),INTENT(INOUT) :: name,class
    CHARACTER(len=MaxStringLen) :: name_l,class_l
    INTEGER, INTENT(OUT) :: Type
    INTEGER :: len_name,len_class

    len_name=LEN(name)
    len_class=LEN(name)

    block_header_start=current_displacement

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, name_l, MaxStringLen, MPI_CHARACTER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, class_l, MaxStringLen, MPI_CHARACTER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + 2 * MaxStringLen
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, Type, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    c_Block_Type=Type

    name=name_l(1:MIN(len_name,MaxStringLen))
    class=class_l(1:MIN(len_class,MaxStringLen))

    current_displacement = current_displacement +  4
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, MPI_INTEGER8,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle,block_md_length,1,MPI_INTEGER8,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle,block_length,1,MPI_INTEGER8,cfd_status,cfd_errcode)

    !Skip past the header block
    current_displacement = block_header_start + block_header_size 

    block_header_end = current_displacement

  END SUBROUTINE cfd_Get_Next_Block_Info_All



  SUBROUTINE cfd_Get_Common_MeshType_MetaData_All(Type,nd,sof)

    !Mesh and mesh variables (and other types such as multimat objects start in the same way
    !An integer type and a dimensionality, so just have one routine
    INTEGER,INTENT(INOUT) :: Type,nd,sof

    CALL cfd_Skip_Block_Header()
    !Now at start of metadata
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, Type, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, nd, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, sof, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + 3 * SoI
  END SUBROUTINE cfd_Get_Common_MeshType_MetaData_All

  SUBROUTINE cfd_Get_Snapshot(time,snap)

    REAL(KIND=8),INTENT(OUT) :: time
    INTEGER,INTENT(OUT) :: snap

    CALL cfd_Skip_Block_Header()
    !Now at start of metadata
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, snap, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement=current_displacement + SoI

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, time, 1, mpireal, cfd_status, cfd_errcode)

    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_Snapshot


  SUBROUTINE cfd_Get_Real_Constant(Value)
    REAL(num),INTENT(OUT) :: Value

    CALL cfd_Skip_Block_Header()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL,cfd_errcode)  
    CALL MPI_FILE_READ_ALL(cfd_filehandle, Value, 1, mpireal, cfd_status, cfd_errcode)
    CALL cfd_Skip_Block()

  END SUBROUTINE cfd_Get_Real_Constant

END MODULE input
