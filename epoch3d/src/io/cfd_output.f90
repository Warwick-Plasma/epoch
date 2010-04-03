MODULE cfd_output

  USE cfd_common
  USE cfd_job_info
  USE shared_data
  USE version_data
  USE mpi

  IMPLICIT NONE

  SAVE

CONTAINS

  SUBROUTINE cfd_open_clobber(h, filename, step, time)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: step
    DOUBLE PRECISION, INTENT(IN) :: time
    INTEGER(4) :: step4, endianness
    INTEGER(8) :: block_header_size
    INTEGER :: errcode

    ! Set the block header
    block_header_size = h%max_string_len * 2_4 + soi + 2_4 * soi8

    ! Delete file and wait
    IF (h%rank .EQ. h%default_rank) &
        CALL MPI_FILE_DELETE(TRIM(filename), MPI_INFO_NULL, errcode)

    CALL MPI_FILE_OPEN(h%comm, TRIM(filename), h%mode, MPI_INFO_NULL, &
        h%filehandle, errcode)

    endianness = 16911887

    ! Currently no blocks written
    h%nblocks = 0

    IF (h%rank .EQ. h%default_rank) THEN
      ! Write the header
      CALL MPI_FILE_WRITE(h%filehandle, "CFD", 3, MPI_CHARACTER, &
          MPI_STATUS_IGNORE, errcode)

      ! This goes next so that stuff can be added to the global header without
      ! breaking everything
      CALL MPI_FILE_WRITE(h%filehandle, c_header_offset, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, block_header_size, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, cfd_version, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, cfd_revision, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, h%max_string_len, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, h%nblocks, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, endianness, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, cfd_jobid%start_seconds, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, cfd_jobid%start_milliseconds, 1, &
          MPI_INTEGER4, MPI_STATUS_IGNORE, errcode)
      step4 = INT(step, 4)
      CALL MPI_FILE_WRITE(h%filehandle, step4, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, time, 1, MPI_DOUBLE_PRECISION, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    ! Current displacement is just the header
    h%current_displacement = c_header_offset

  END SUBROUTINE cfd_open_clobber



  SUBROUTINE cfd_safe_write_string(h, string)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=h%max_string_len) :: output
    INTEGER :: len_s, errcode

    len_s = LEN(string)

    IF (h%max_string_len .LT. len_s .AND. h%rank .EQ. h%default_rank) THEN
      PRINT*, '***WARNING***'
      PRINT*, 'Output string "' // string // '" has been truncated'
    ENDIF

    ! This subroutine expects that the record marker is in place and that
    ! the view is set correctly. Call it only on the node which is doing the
    ! writing. You still have to advance the file pointer yourself on all nodes

    output(1:MIN(h%max_string_len, len_s)) = &
        string(1:MIN(h%max_string_len, len_s))

    ! If this isn't the full string length then tag in a ACHAR(0) to help
    ! With C++ string handling
    IF (len_s + 1 .LT. h%max_string_len) &
        output(len_s+1:h%max_string_len) = ACHAR(0)

    CALL MPI_FILE_WRITE(h%filehandle, output, h%max_string_len, MPI_CHARACTER, &
        MPI_STATUS_IGNORE, errcode)

  END SUBROUTINE cfd_safe_write_string



  SUBROUTINE cfd_write_block_header(h, block_name, block_class, block_type, &
      block_length, block_md_length, rank_write)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: block_name, block_class
    INTEGER(4), INTENT(IN) :: block_type
    INTEGER(8), INTENT(IN) :: block_length, block_md_length
    INTEGER, INTENT(IN) :: rank_write
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(h, block_name)
      CALL cfd_safe_write_string(h, block_class)
    ENDIF
    h%current_displacement = h%current_displacement + 2 * h%max_string_len

    ! Write the block type
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER4, MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(h%filehandle, block_type, 1, MPI_INTEGER4, &
            MPI_STATUS_IGNORE, errcode)

    h%current_displacement = h%current_displacement + 4

    ! Write the block skip and metadata skip data
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(h%filehandle, block_md_length, 1, MPI_INTEGER8, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, block_length, 1, MPI_INTEGER8, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + 2 * 8

    h%nblocks = h%nblocks + 1_4

  END SUBROUTINE cfd_write_block_header



  SUBROUTINE cfd_write_meshtype_header(h, meshtype, dim, sof, rank_write)

    ! MeshTypes (Meshes, fluid variables, multimat blocks etc)
    ! All have a common header, this is what writes that (although the content
    ! Of type will depend on what meshtype you're using)

    TYPE(cfd_file_handle) :: h
    INTEGER(4), INTENT(IN) :: meshtype, dim
    INTEGER(8), INTENT(IN) :: sof
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(4) :: sof4
    INTEGER :: errcode

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(h%filehandle, meshtype, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, dim, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      sof4 = INT(sof,4)
      CALL MPI_FILE_WRITE(h%filehandle, sof4, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + c_meshtype_header_offset

  END SUBROUTINE cfd_write_meshtype_header



  SUBROUTINE cfd_write_job_info(h, code_io_version, version, revision, &
      defines, compile_date, run_date, restart_flag, code_name, commit_id, &
      sha1sum, compile_machine, compile_flags, rank_write)

    TYPE(cfd_file_handle) :: h
    INTEGER(4), INTENT(IN) :: code_io_version, version, revision, defines
    INTEGER(4), INTENT(IN) :: compile_date, run_date
    INTEGER, INTENT(IN) :: restart_flag
    CHARACTER(LEN=*), INTENT(IN) :: code_name, commit_id, sha1sum
    CHARACTER(LEN=*), INTENT(IN) :: compile_machine, compile_flags
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length
    INTEGER(4) :: io_date, restart_flag4
    INTEGER :: errcode

    io_date = get_unix_time()

    md_length = 8 * soi + 4 * h%max_string_len

    CALL cfd_write_block_header(h, code_name, "Job_info", c_type_info, &
        md_length, md_length, rank_write)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(h%filehandle, code_io_version, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, version, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, revision, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + 3 * soi

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(h, commit_id)
      CALL cfd_safe_write_string(h, sha1sum)
      CALL cfd_safe_write_string(h, compile_machine)
      CALL cfd_safe_write_string(h, compile_flags)
    ENDIF

    h%current_displacement = h%current_displacement + 4 * h%max_string_len

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(h%filehandle, defines, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, compile_date, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, run_date, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, io_date, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      restart_flag4 = INT(restart_flag, 4)
      CALL MPI_FILE_WRITE(h%filehandle, restart_flag4, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + 5 * soi

  END SUBROUTINE cfd_write_job_info



  SUBROUTINE cfd_write_stitched_vector(h, vector_name, vector_class, &
      mesh_name, mesh_class, name, class, rank_write)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: vector_name, vector_class
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: name, class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length, block_length
    INTEGER(4) :: ndims
    INTEGER :: iloop, errcode

    ndims = INT(SIZE(name),4)

    md_length = 2 * h%max_string_len + soi
    block_length = md_length + ndims * 2 * h%max_string_len

    CALL cfd_write_block_header(h, vector_name, vector_class, &
        c_type_stitched_vector, block_length, md_length, rank_write)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(h, mesh_name)
      CALL cfd_safe_write_string(h, mesh_class)
    ENDIF

    h%current_displacement = h%current_displacement + 2 * h%max_string_len

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(h%filehandle, ndims, 1, MPI_INTEGER4, &
            MPI_STATUS_IGNORE, errcode)

    h%current_displacement = h%current_displacement + soi

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      DO iloop = 1, ndims
        CALL cfd_safe_write_string(h, name(iloop))
        CALL cfd_safe_write_string(h, class(iloop))
      ENDDO
    ENDIF

    h%current_displacement = &
        h%current_displacement + 2 * ndims * h%max_string_len

  END SUBROUTINE cfd_write_stitched_vector



  SUBROUTINE cfd_write_stitched_magnitude(h, magn_name, magn_class, mesh_name, &
      mesh_class, name, class, rank_write)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: magn_name, magn_class
    CHARACTER(LEN=*), INTENT(IN) :: mesh_name, mesh_class
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: name, class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length, block_length
    INTEGER(4) :: ndims
    INTEGER :: iloop, errcode

    ndims = INT(SIZE(name),4)

    md_length = 2 * h%max_string_len + soi
    block_length = md_length + ndims * 2 * h%max_string_len

    CALL cfd_write_block_header(h, magn_name, magn_class, &
        c_type_stitched_magnitude, block_length, md_length, rank_write)

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(h, mesh_name)
      CALL cfd_safe_write_string(h, mesh_class)
    ENDIF

    h%current_displacement = h%current_displacement + 2 * h%max_string_len

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(h%filehandle, ndims, 1, MPI_INTEGER4, &
            MPI_STATUS_IGNORE, errcode)

    h%current_displacement = h%current_displacement + soi

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      DO iloop = 1, ndims
        CALL cfd_safe_write_string(h, name(iloop))
        CALL cfd_safe_write_string(h, class(iloop))
      ENDDO
    ENDIF

    h%current_displacement = &
        h%current_displacement + 2 * ndims * h%max_string_len

  END SUBROUTINE cfd_write_stitched_magnitude



  SUBROUTINE cfd_write_real_constant(h, name, class, value, rank_write)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length
    INTEGER :: errcode

    md_length = sof

    CALL cfd_write_block_header(h, name, class, c_type_constant, md_length, &
        md_length, rank_write)
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, h%mpireal, &
        h%mpireal, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(h%filehandle, value, 1, h%mpireal, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + sof

  END SUBROUTINE cfd_write_real_constant



  SUBROUTINE cfd_write_source_code(h, name, class, array, last, rank_write)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: array
    CHARACTER(LEN=*), INTENT(IN) :: last
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length, sz, len1, len2
    INTEGER :: i, errcode

    IF (h%rank .EQ. rank_write) THEN
      sz   = SIZE(array)
      len1 = LEN(array)
      len2 = LEN(last)
      md_length = sz*len1 + len2
    ENDIF

    CALL MPI_BCAST(md_length, 1, MPI_INTEGER8, 0, h%comm, errcode)

    CALL cfd_write_block_header(h, name, class, c_type_constant, md_length, &
        md_length, rank_write)
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, &
        MPI_CHARACTER,  MPI_CHARACTER, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      DO i = 1, sz
        CALL MPI_FILE_WRITE(h%filehandle, array(i), len1, &
            MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)
      ENDDO
      CALL MPI_FILE_WRITE(h%filehandle, last, len2, &
          MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + md_length

  END SUBROUTINE cfd_write_source_code



  SUBROUTINE cfd_write_1d_integer_array(h, name, class, values, rank_write)

    TYPE(cfd_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, DIMENSION(:), INTENT(IN) :: values
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length
    INTEGER(4) :: sz
    INTEGER :: errcode

    md_length = 3 * soi

    CALL cfd_write_block_header(h, name, class, c_type_integerarray, &
        md_length, md_length, rank_write)
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_displacement, MPI_INTEGER4, &
        MPI_INTEGER4, "native", MPI_INFO_NULL, errcode)

    IF (h%rank .EQ. rank_write) THEN
      ! 1D
      CALL MPI_FILE_WRITE(h%filehandle, 1, c_dimension_1d, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      ! INTEGER kind
      sz = KIND(values)
      CALL MPI_FILE_WRITE(h%filehandle, 1, sz, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      ! Size of array
      CALL MPI_FILE_WRITE(h%filehandle, 1, SIZE(values), MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
      ! Actual array
      CALL MPI_FILE_WRITE(h%filehandle, values, SIZE(values), MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    ENDIF

    h%current_displacement = h%current_displacement + md_length

  END SUBROUTINE cfd_write_1d_integer_array

END MODULE cfd_output
