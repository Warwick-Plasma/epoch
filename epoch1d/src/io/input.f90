MODULE input

  USE input_functions

  IMPLICIT NONE

  SAVE

CONTAINS

  SUBROUTINE cfd_open_read(filename, step, time)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(4), INTENT(OUT) :: step
    DOUBLE PRECISION, INTENT(OUT) :: time
    CHARACTER(LEN=3) :: cfd
    INTEGER(4) :: file_version, file_revision
    INTEGER(4) :: endianness
    INTEGER :: ierr

    step = -1

    CALL MPI_BARRIER(cfd_comm, cfd_errcode)

    CALL MPI_FILE_OPEN(cfd_comm, TRIM(filename), cfd_mode, MPI_INFO_NULL, &
        cfd_filehandle, cfd_errcode)

    current_displacement = 0
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    ! Read the header
    CALL MPI_FILE_READ_ALL(cfd_filehandle, cfd, 3, MPI_CHARACTER, &
        cfd_status, cfd_errcode)

    ! If this isn't "CFD" then this isn't a CFD file
    IF (cfd .NE. "CFD") THEN
      CALL MPI_FILE_CLOSE(cfd_filehandle, cfd_errcode)
      IF (rank .EQ. default_rank) &
          PRINT *, "The specified file is not a valid CFD file"
      CALL MPI_ABORT(cfd_comm, cfd_errcode, ierr)
    ENDIF

    current_displacement = 3

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    ! Read in the basic file info. Should check version info, but this is
    ! version 1, so let's not worry about it
    CALL MPI_FILE_READ_ALL(cfd_filehandle, header_offset, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, block_header_size, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, file_version, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, file_revision, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, max_string_len, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, cfd_nblocks, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, endianness, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, cfd_jobid%start_seconds, 1, &
        MPI_INTEGER4, cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, cfd_jobid%start_milliseconds, 1, &
        MPI_INTEGER4, cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, step, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, time, 1, MPI_DOUBLE_PRECISION, &
        cfd_status, cfd_errcode)

    IF (file_version .GT. cfd_version) THEN
      IF (rank .EQ. default_rank) PRINT *, "Version number incompatible"
      CALL MPI_ABORT(cfd_comm, cfd_errcode, ierr)
    ENDIF

    IF (file_revision .GT. cfd_revision) THEN
      IF (rank .EQ. default_rank) &
          PRINT *, "Revision number of file is too high. Writing disabled"
      cfd_writing = .FALSE.
    ENDIF

    current_displacement = header_offset

  END SUBROUTINE cfd_open_read



  SUBROUTINE cfd_get_next_block_info_all(name, class, block_type)

    CHARACTER(LEN=*), INTENT(INOUT) :: name, class
    CHARACTER(LEN=max_string_len) :: name_l, class_l
    INTEGER, INTENT(OUT) :: block_type
    INTEGER(4) :: block_type4
    INTEGER :: len_name, len_class

    len_name = LEN(name)
    len_class = LEN(name)

    block_header_start = current_displacement

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, name_l, max_string_len, &
        MPI_CHARACTER, cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, class_l, max_string_len, &
        MPI_CHARACTER, cfd_status, cfd_errcode)

    current_displacement = current_displacement + 2 * max_string_len
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, block_type4, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    block_type = block_type4

    name = name_l(1:MIN(len_name, max_string_len))
    class = class_l(1:MIN(len_class, max_string_len))

    current_displacement = current_displacement +  4
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, block_md_length, 1, &
        MPI_INTEGER8, cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, block_length, 1, &
        MPI_INTEGER8, cfd_status, cfd_errcode)

    ! Skip past the header block
    current_displacement = block_header_start + block_header_size

    block_header_end = current_displacement

  END SUBROUTINE cfd_get_next_block_info_all



  SUBROUTINE cfd_get_common_meshtype_metadata_all(meshtype, nd, sof)

    ! Mesh and mesh variables (and other types such as multimat objects start
    ! in the same way). An integer type and a dimensionality, so just have one
    ! routine

    INTEGER, INTENT(INOUT) :: meshtype, nd, sof
    INTEGER(4) :: meshtype4, nd4, sof4

    CALL cfd_skip_block_header()

    ! Now at start of metadata
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, meshtype4, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, nd4, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, sof4, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    meshtype = meshtype4
    nd = nd4
    sof = sof4

    current_displacement = current_displacement + 3 * soi

  END SUBROUTINE cfd_get_common_meshtype_metadata_all



  SUBROUTINE cfd_get_snapshot(time, snap)

    REAL(8), INTENT(OUT) :: time
    INTEGER, INTENT(OUT) :: snap
    INTEGER(4) :: snap4

    CALL cfd_skip_block_header()

    ! Now at start of metadata
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, snap4, 1, MPI_INTEGER4, &
        cfd_status, cfd_errcode)

    snap = snap4

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, "native", &
        MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, time, 1, MPI_DOUBLE_PRECISION, &
        cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_snapshot



  SUBROUTINE cfd_get_real_constant(value)

    REAL(num), INTENT(OUT) :: value

    CALL cfd_skip_block_header()

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    CALL MPI_FILE_READ_ALL(cfd_filehandle, value, 1, mpireal, &
        cfd_status, cfd_errcode)

    CALL cfd_skip_block()

  END SUBROUTINE cfd_get_real_constant

END MODULE input
