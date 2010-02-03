
MODULE output

  USE shared_data
  USE iocommon

  IMPLICIT NONE

  SAVE

  PRIVATE


  PUBLIC :: cfd_Open_Clobber,cfd_Write_Block_Header,cfd_Write_MeshType_Header
  PUBLIC :: cfd_Safe_Write_String,cfd_Write_Snapshot_Data,cfd_Write_Stitched_Vector
  PUBLIC :: cfd_Write_Stitched_Magnitude,cfd_Write_Real_Constant
  PUBLIC :: cfd_Write_VisIT_Expression

CONTAINS

  SUBROUTINE cfd_Open_Clobber(filename)

    CHARACTER(len=*),INTENT(IN) :: filename

    !Set the block header
    block_header_size = MaxStringLen * 2 + 4 + 2 * 8

    !Delete file and wait
    IF (cfd_rank == default_rank) CALL MPI_FILE_DELETE(TRIM(filename), MPI_INFO_NULL, cfd_errcode)
    CALL MPI_BARRIER(cfd_comm,cfd_errcode)
    CALL MPI_FILE_OPEN(cfd_comm, TRIM(filename), cfd_mode, &
         MPI_INFO_NULL, cfd_filehandle, cfd_errcode)
	 CALL MPI_FILE_SET_ATOMICITY(cfd_filehandle,0,cfd_errcode)

    !

    IF (cfd_rank == default_rank) THEN
       !Write the header
       CALL MPI_FILE_WRITE(cfd_filehandle, "CFD", 3, MPI_CHARACTER, cfd_status, cfd_errcode)
       !This goes next so that stuff can be added to the global header without breaking 
       !Everything
       CALL MPI_FILE_WRITE(cfd_filehandle, header_offset, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, block_header_size, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, CFD_Version, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, CFD_Revision, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, MaxStringLen, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       !This is where the nblocks variable will go, put a zero for now
       CALL MPI_FILE_WRITE(cfd_filehandle, 0, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    !Currently no blocks written
    nblocks=0
    !Current displacement is just the header
    current_displacement=header_offset

  END SUBROUTINE cfd_Open_Clobber

  SUBROUTINE cfd_safe_write_string(string)

    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=MaxStringLen) :: output
    CHARACTER,DIMENSION(MaxStringLen) :: mask
    INTEGER :: len_s

    len_s=LEN(string)

    !This subroutine expects that the record marker is in place and that
    !the view is set correctly. Call it only on the node which is doing the writing
    !You still have to advance the file pointer yourself on all nodes

    output(1:MIN(MaxStringLen,len_s))=string(1:MIN(MaxStringLen,len_s))
    !If this isn't the full string length then tag in a ACHAR(0) to help 
    !With C++ string handling
    IF (len_s+1 < MaxStringLen) output(len_s+1:MaxStringLen)=ACHAR(0)
    CALL MPI_FILE_WRITE(cfd_filehandle,output,MaxStringLen,MPI_CHARACTER,cfd_status,cfd_errcode)

  END SUBROUTINE cfd_safe_write_string

  SUBROUTINE cfd_write_block_header(blockname,blockclass,blocktype,blocklength,blockmetadatalength,rank_write)

   CHARACTER(len=*), INTENT(IN) :: blockname, blockclass
    INTEGER, INTENT(IN) :: blocktype,rank_write
    INTEGER(KIND=8),INTENT(IN) :: blocklength,blockmetadatalength
    INTEGER :: len_bn, len_bc
    CHARACTER(len=MaxStringLen) :: output

    len_bn=LEN(blockname)
    len_bc=LEN(blockclass)


    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
       CALL cfd_Safe_Write_String(blockname)
       CALL cfd_Safe_Write_String(blockclass)
    ENDIF
    current_displacement = current_displacement + 2 * MaxStringLen

    !Write the block type
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write)&
         CALL MPI_FILE_WRITE(cfd_filehandle, blocktype, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + 4

    !Write the block skip and metadata skip data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, MPI_INTEGER8,&
         "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank == rank_write) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, blockmetadatalength, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, blocklength, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * 8

    nblocks=nblocks+1

  END SUBROUTINE cfd_write_block_header

  SUBROUTINE cfd_Write_MeshType_Header(type,dim,sof,rank_write)
    !MeshTypes (Meshes, fluid variables, multimat blocks etc)
    !All have a common header, this is what writes that (although the content
    !Of type will depend on what meshtype you're using)

    INTEGER, INTENT(IN) :: type,dim,rank_write,sof

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)

    IF (cfd_rank == rank_write) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, type, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, dim , 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, sof , 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF

    Current_Displacement = Current_Displacement + MeshType_Header_Offset

  END SUBROUTINE cfd_Write_MeshType_Header
 
  SUBROUTINE cfd_Write_Snapshot_Data(time,cycle,rank_write)

    INTEGER, INTENT(IN) :: rank_write,cycle
    INTEGER(8) :: mdlength
    REAL(8),INTENT(IN) :: time

    mdlength = SoI + num

    CALL cfd_Write_Block_Header("Snapshot","Snapshot",TYPE_SNAPSHOT,mdlength,mdlength,rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == rank_write) &
         CALL MPI_FILE_WRITE(cfd_filehandle, cycle, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement=current_displacement + SoI

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == rank_write) &
         CALL MPI_FILE_WRITE(cfd_filehandle, time, 1, MPI_DOUBLE_PRECISION, cfd_status, cfd_errcode)

    current_displacement=current_displacement + 8

  END SUBROUTINE cfd_Write_Snapshot_Data

  SUBROUTINE cfd_Write_Stitched_Vector(vector_name,vector_class,mesh_name,mesh_class,name,class,rank_write)

    CHARACTER(len=*),DIMENSION(:), INTENT(IN) :: name,class
    CHARACTER(len=*),INTENT(IN) :: vector_name,vector_class,mesh_name,mesh_class
    INTEGER,INTENT(IN) :: rank_write
    INTEGER(8) :: n_Dims,mdlength,blocklength
    INTEGER :: iLoop

    n_Dims=SIZE(name)

    mdlength = 2 * MaxStringLen + SoI
    blocklength=mdlength + n_Dims * 2 * MaxStringLen

    CALL cfd_Write_Block_Header(vector_name,vector_class,TYPE_STITCHED_VECTOR,blocklength,mdlength,rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
       CALL cfd_Safe_Write_String(mesh_name)
       CALL cfd_Safe_Write_String(mesh_class)
    ENDIF
    current_displacement = current_displacement + 2 * MaxStringLen

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)  
    IF (cfd_rank == rank_write) &
         CALL MPI_FILE_WRITE(cfd_filehandle, n_Dims, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + SoI

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL,cfd_errcode)  
    IF (cfd_rank == rank_write) THEN
       DO iLoop=1,n_Dims
          CALL cfd_Safe_Write_String(name(iLoop))
          CALL cfd_Safe_Write_String(class(iLoop))
       ENDDO
    ENDIF
    current_displacement=current_displacement + 2 * n_Dims* MaxStringLen

  END SUBROUTINE cfd_Write_Stitched_Vector

  SUBROUTINE cfd_Write_Stitched_Magnitude(Magn_name,magn_class,mesh_name,mesh_class,name,class,rank_write)

    CHARACTER(len=*),DIMENSION(:), INTENT(IN) :: name,class
    CHARACTER(len=*),INTENT(IN) :: magn_name,magn_class,mesh_name,mesh_class
    INTEGER,INTENT(IN) :: rank_write
    INTEGER(8) :: n_Dims,mdlength,blocklength
    INTEGER :: iLoop

    n_Dims=SIZE(name)

    mdlength = 2 * MaxStringLen + SoI
    blocklength=mdlength + n_Dims * 2 * MaxStringLen

    CALL cfd_Write_Block_Header(magn_name,magn_class,TYPE_STITCHED_MAGNITUDE,blocklength,mdlength,rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == rank_write) THEN
       CALL cfd_Safe_Write_String(mesh_name)
       CALL cfd_Safe_Write_String(mesh_class)
    ENDIF
    current_displacement = current_displacement + 2 * MaxStringLen

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)  
    IF (cfd_rank == rank_write) &
         CALL MPI_FILE_WRITE(cfd_filehandle, n_Dims, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + SoI

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL,cfd_errcode)  
    IF (cfd_rank == rank_write) THEN
       DO iLoop=1,n_Dims
          CALL cfd_Safe_Write_String(name(iLoop))
          CALL cfd_Safe_Write_String(class(iLoop))
       ENDDO
    ENDIF
    current_displacement=current_displacement + 2 * n_Dims* MaxStringLen

  END SUBROUTINE cfd_Write_Stitched_Magnitude

  SUBROUTINE cfd_Write_Real_Constant(Name,Class,Value,rank_write)
    CHARACTER(len=*), INTENT(IN) :: name,class
    REAL(num),INTENT(IN) :: Value
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: mdlength

    mdlength=num

    CALL cfd_Write_Block_Header(name,class,TYPE_CONSTANT,mdlength,mdlength,rank_write)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL,cfd_errcode)  
    IF (cfd_rank == rank_write) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, Value, 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    Current_displacement=Current_displacement + num
    

  END SUBROUTINE cfd_Write_Real_Constant

  SUBROUTINE cfd_Write_1D_Integer_Array(Name,Class,Values,rank_write)
    CHARACTER(len=*), INTENT(IN) :: name,class
    INTEGER, DIMENSION(:), INTENT(IN) :: Values
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: mdlength

    mdlength=2 * SoI

    CALL cfd_Write_Block_Header(name,class,TYPE_INTEGERARRAY,mdlength,mdlength,rank_write)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)  
    IF (cfd_rank == rank_write) THEN
       !1D
       CALL MPI_FILE_WRITE(cfd_filehandle, 1, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       !Size of array
       CALL MPI_FILE_WRITE(cfd_filehandle, 1, SIZE(Values), MPI_INTEGER, cfd_status, cfd_errcode)
       !Actual Array
       CALL MPI_FILE_WRITE(cfd_filehandle, Values, SIZE(Values), MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    Current_displacement=Current_displacement + mdlength

  END SUBROUTINE cfd_Write_1D_Integer_Array

  SUBROUTINE cfd_Write_VisIT_Expression(ExpressionName,ExpressionClass,Expression)

    CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: ExpressionName,ExpressionClass,Expression

    PRINT *,LEN(Expression(1)),LEN(Expression(2))

  END SUBROUTINE cfd_Write_VisIT_Expression

END MODULE output







