MODULE cfd_input

  USE cfd_common
  USE cfd_input_functions
  USE mpi

  IMPLICIT NONE

  SAVE

CONTAINS

  SUBROUTINE cfd_open_read(h, filename, step, time)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(4), INTENT(OUT) :: step
    DOUBLE PRECISION, INTENT(OUT) :: time
    CHARACTER(LEN=3) :: cfd
    INTEGER(4) :: file_version, file_revision
    INTEGER(4) :: endianness
    INTEGER :: errcode, ierr

    CALL MPI_FILE_OPEN(h%comm, TRIM(filename), h%mode, MPI_INFO_NULL, &
        h%filehandle, errcode)

    h%current_displacement = 0
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    ! Read the header
    CALL MPI_FILE_READ_ALL(h%filehandle, cfd, 3, MPI_CHARACTER, &
        MPI_STATUS_IGNORE, errcode)

    ! If this isn't "CFD" then this isn't a CFD file
    IF (cfd .NE. "CFD") THEN
      CALL MPI_FILE_CLOSE(h%filehandle, errcode)
      IF (h%rank .EQ. h%default_rank) &
          PRINT *, "The specified file is not a valid CFD file"
      CALL MPI_ABORT(h%comm, errcode, ierr)
    ENDIF

    h%current_displacement = 3

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    ! Read in the basic file info. Should check version info, but this is
    ! version 1, so let's not worry about it
    CALL MPI_FILE_READ_ALL(h%filehandle, h%header_offset, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, h%block_header_size, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, file_version, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, file_revision, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, h%max_string_len, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, h%nblocks, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, endianness, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, cfd_jobid%start_seconds, 1, &
        MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, cfd_jobid%start_milliseconds, 1, &
        MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, step, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, time, 1, MPI_DOUBLE_PRECISION, &
        MPI_STATUS_IGNORE, errcode)

    IF (file_version .GT. cfd_version) THEN
      IF (h%rank .EQ. h%default_rank) PRINT *, "Version number incompatible"
      CALL MPI_ABORT(h%comm, errcode, ierr)
    ENDIF

    IF (file_revision .GT. cfd_revision) THEN
      IF (h%rank .EQ. h%default_rank) &
          PRINT *, "Revision number of file is too high. Writing disabled"
      h%writing = .FALSE.
    ENDIF

    h%current_displacement = h%header_offset

  END SUBROUTINE cfd_open_read



  SUBROUTINE cfd_get_next_block_info_all(h, name, class, block_type)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(INOUT) :: name, class
    INTEGER, INTENT(OUT) :: block_type
    CHARACTER(LEN=h%max_string_len) :: name_l, class_l
    INTEGER(4) :: block_type4
    INTEGER :: len_name, len_class, errcode
    INTEGER(KIND=MPI_OFFSET_KIND) :: block_header_start

    len_name = LEN(name)
    len_class = LEN(name)

    block_header_start = h%current_displacement

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, name_l, h%max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, class_l, h%max_string_len, &
        MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)

    h%current_displacement = h%current_displacement + 2 * h%max_string_len
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, block_type4, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    block_type = block_type4

    name = name_l(1:MIN(len_name, h%max_string_len))
    class = class_l(1:MIN(len_class, h%max_string_len))

    h%current_displacement = h%current_displacement +  4
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, h%block_md_length, 1, &
        MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, h%block_length, 1, &
        MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

    ! Skip past the header block
    h%current_displacement = block_header_start + h%block_header_size

    h%block_header_end = h%current_displacement

  END SUBROUTINE cfd_get_next_block_info_all



  SUBROUTINE cfd_get_common_meshtype_metadata_all(h, meshtype, nd, sof)

    ! Mesh and mesh variables (and other types such as multimat objects start
    ! in the same way). An integer type and a dimensionality, so just have one
    ! routine

    TYPE(cfd_file_handle) :: h
    INTEGER, INTENT(INOUT) :: meshtype, nd, sof
    INTEGER(4) :: meshtype4, nd4, sof4
    INTEGER :: errcode

    CALL cfd_skip_block_header(h)

    ! Now at start of metadata
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, meshtype4, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, nd4, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, sof4, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    meshtype = meshtype4
    nd = nd4
    sof = sof4

    h%current_displacement = h%current_displacement + 3 * soi

  END SUBROUTINE cfd_get_common_meshtype_metadata_all



  SUBROUTINE cfd_get_snapshot(h, time, snap)

    TYPE(cfd_file_handle) :: h
    REAL(8), INTENT(OUT) :: time
    INTEGER, INTENT(OUT) :: snap
    INTEGER(4) :: snap4
    INTEGER :: errcode

    CALL cfd_skip_block_header(h)

    ! Now at start of metadata
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, snap4, 1, MPI_INTEGER4, &
        MPI_STATUS_IGNORE, errcode)

    snap = snap4

    h%current_displacement = h%current_displacement + soi

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, "native", &
        MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, time, 1, MPI_DOUBLE_PRECISION, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_snapshot



  SUBROUTINE cfd_get_real_constant(h, value)

    TYPE(cfd_file_handle) :: h
    REAL(num), INTENT(OUT) :: value
    INTEGER :: errcode

    CALL cfd_skip_block_header(h)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, value, 1, h%mpireal, &
        MPI_STATUS_IGNORE, errcode)

    CALL cfd_skip_block(h)

  END SUBROUTINE cfd_get_real_constant

END MODULE cfd_input
