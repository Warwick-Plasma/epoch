
MODULE input

  USE shared_data
  USE iocommon
  USE input_functions
  USE input_cartesian

  IMPLICIT NONE

  SAVE

CONTAINS


  SUBROUTINE cfd_open_read(filename)

    CHARACTER(len=*),INTENT(IN) :: filename
    CHARACTER(len=3) :: cfd

    INTEGER:: file_version,file_revision

    CALL MPI_BARRIER(cfd_comm,cfd_errcode)
    CALL MPI_FILE_OPEN(cfd_comm, TRIM(filename), cfd_mode, &
         MPI_INFO_NULL, cfd_filehandle, cfd_errcode)

    current_displacement=0
       CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
          "native", MPI_INFO_NULL, cfd_errcode)
       !Read the header
       CALL MPI_FILE_READ_ALL(cfd_filehandle, cfd, 3, MPI_CHARACTER, cfd_status, cfd_errcode)
       !If this isn't "CFD" then this isn't a CFD file
       IF (cfd /= "CFD") THEN
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
       CALL MPI_FILE_READ_ALL(cfd_filehandle, file_version, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, file_revision, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, max_string_len, 1, MPI_INTEGER, cfd_status, cfd_errcode)
       CALL MPI_FILE_READ_ALL(cfd_filehandle, nblocks, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    IF (file_version .GT. cfd_version) THEN
       IF (rank == default_rank) PRINT *,"Version number incompatible"
       CALL MPI_ABORT(cfd_comm,cfd_errcode)
    ENDIF

    IF (file_revision .GT. cfd_revision) THEN
       IF (rank == default_rank) PRINT *,"Revision number of file is too high. Writing disabled"
       cfd_writing=.FALSE.
    ENDIF

    current_displacement = header_offset

  END SUBROUTINE cfd_open_read

  SUBROUTINE cfd_get_next_block_info_all(name,class,block_type)

    CHARACTER(len=*),INTENT(INOUT) :: name,class
    CHARACTER(len=max_string_len) :: name_l,class_l
    INTEGER, INTENT(OUT) :: block_type
    INTEGER :: len_name,len_class

    len_name=len(name)
    len_class=len(name)

    block_header_start=current_displacement

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, name_l, max_string_len, MPI_CHARACTER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, class_l, max_string_len, MPI_CHARACTER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + 2 * max_string_len
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, block_type, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    name=name_l(1:MIN(len_name,max_string_len))
    class=class_l(1:MIN(len_class,max_string_len))

    current_displacement = current_displacement +  4
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, MPI_INTEGER8,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle,block_md_length,1,MPI_INTEGER8,cfd_status,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle,block_length,1,MPI_INTEGER8,cfd_status,cfd_errcode)

    !Skip past the header block
    current_displacement = block_header_start + block_header_size

    block_header_end = current_displacement

  END SUBROUTINE cfd_get_next_block_info_all



  SUBROUTINE cfd_get_common_meshtype_metadata_all(meshtype,nd,sof)

    !Mesh and mesh variables (and other types such as multimat objects start in the same way
    !an integer type and a dimensionality, so just have one routine
    INTEGER,INTENT(INOUT) :: meshtype,nd,sof

    CALL cfd_skip_block_header()
    !Now at start of metadata
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, meshtype, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, nd, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, sof, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + 3 * soi
  END SUBROUTINE cfd_get_common_meshtype_metadata_all

  SUBROUTINE cfd_get_snapshot(time,snap)

    REAL(KIND=8),INTENT(OUT) :: time
    INTEGER,INTENT(OUT) :: snap

    CALL cfd_skip_block_header()
    !Now at start of metadata
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, snap, 1, MPI_INTEGER, cfd_status, cfd_errcode)

    current_displacement=current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, time, 1, mpireal, cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_snapshot


  SUBROUTINE cfd_get_real_constant(value)
    REAL(num),INTENT(OUT) :: value

    CALL cfd_skip_block_header()
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_READ_ALL(cfd_filehandle, value, 1, mpireal, cfd_status, cfd_errcode)
    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_real_constant

END MODULE input
