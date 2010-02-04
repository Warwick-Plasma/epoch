MODULE output

  USE iocommon

  IMPLICIT NONE

  SAVE

  PRIVATE

  PUBLIC :: cfd_open_clobber, cfd_write_block_header, cfd_write_meshtype_header
  PUBLIC :: cfd_safe_write_string, cfd_write_snapshot_data
  PUBLIC :: cfd_write_stitched_vector
  PUBLIC :: cfd_write_stitched_magnitude, cfd_write_real_constant
  PUBLIC :: cfd_write_visit_expression

CONTAINS

  SUBROUTINE cfd_open_clobber(filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename

    ! Set the block header
    block_header_size = max_string_len * 2 + 4 + 2 * 8

    ! Delete file and wait
    IF (cfd_rank .EQ. default_rank) &
        CALL MPI_FILE_DELETE(TRIM(filename), MPI_INFO_NULL, cfd_errcode)

    CALL MPI_BARRIER(cfd_comm, cfd_errcode)
    CALL MPI_FILE_OPEN(cfd_comm, TRIM(filename), cfd_mode, MPI_INFO_NULL, &
        cfd_filehandle, cfd_errcode)

    IF (cfd_rank .EQ. default_rank) THEN
      ! Write the header
      CALL MPI_FILE_WRITE(cfd_filehandle, "CFD", 3, MPI_CHARACTER, &
          cfd_status, cfd_errcode)

      ! This goes next so that stuff can be added to the global header without
      ! breaking everything
      CALL MPI_FILE_WRITE(cfd_filehandle, header_offset, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, block_header_size, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, cfd_version, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, cfd_revision, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, max_string_len, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)

      ! This is where the nblocks variable will go, put a zero for now
      CALL MPI_FILE_WRITE(cfd_filehandle, 0, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
    ENDIF

    ! Currently no blocks written
    nblocks = 0
    ! Current displacement is just the header
    current_displacement = header_offset

  END SUBROUTINE cfd_open_clobber



  SUBROUTINE cfd_safe_write_string(string)

    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=max_string_len) :: output
    INTEGER :: len_s

    len_s = LEN(string)

    ! This subroutine expects that the record marker is in place and that
    ! the view is set correctly. Call it only on the node which is doing the
    ! writing. You still have to advance the file pointer yourself on all nodes

    output(1:MIN(max_string_len, len_s)) = string(1:MIN(max_string_len, len_s))

    ! If this isn't the full string length then tag in a ACHAR(0) to help
    ! With C++ string handling
    IF (len_s+1 .LT. max_string_len) output(len_s+1:max_string_len) = ACHAR(0)

    CALL MPI_FILE_WRITE(cfd_filehandle, output, max_string_len, MPI_CHARACTER, &
        cfd_status, cfd_errcode)

  END SUBROUTINE cfd_safe_write_string



  SUBROUTINE cfd_write_block_header(block_name, block_class, block_type, &
      block_length, block_md_length, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: block_name, block_class
    INTEGER, INTENT(IN) :: block_type, rank_write
    INTEGER(KIND=8), INTENT(IN) :: block_length, block_md_length
    INTEGER :: len_bn, len_bc

    len_bn = LEN(block_name)
    len_bc = LEN(block_class)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(block_name)
      CALL cfd_safe_write_string(block_class)
    ENDIF
    current_displacement = current_displacement + 2 * max_string_len

    ! Write the block type
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, block_type, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + 4

    ! Write the block skip and metadata skip data
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_INTEGER8, MPI_INTEGER8, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, block_md_length, 1, MPI_INTEGER8, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, block_length, 1, MPI_INTEGER8, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + 2 * 8

    nblocks = nblocks + 1

  END SUBROUTINE cfd_write_block_header



  SUBROUTINE cfd_write_meshtype_header(meshtype, dim, sof, rank_write)

    ! MeshTypes (Meshes, fluid variables, multimat blocks etc)
    ! All have a common header, this is what writes that (although the content
    ! Of type will depend on what meshtype you're using)

    INTEGER, INTENT(IN) :: meshtype, dim, rank_write, sof

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, meshtype, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, dim, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      CALL MPI_FILE_WRITE(cfd_filehandle, sof, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + meshtype_header_offset

  END SUBROUTINE cfd_write_meshtype_header



  SUBROUTINE cfd_write_snapshot_data(time, step, rank_write)

    INTEGER, INTENT(IN) :: rank_write, step
    INTEGER(8) :: md_length
    REAL(8), INTENT(IN) :: time

    md_length = soi + num

    CALL cfd_write_block_header("Snapshot", "Snapshot", c_type_snapshot, &
        md_length, md_length, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, step, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, "native", MPI_INFO_NULL, &
        cfd_errcode)

    IF (cfd_rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, time, 1, MPI_DOUBLE_PRECISION, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + 8

  END SUBROUTINE cfd_write_snapshot_data



  SUBROUTINE cfd_write_stitched_vector(vector_name, vector_class, mesh_name, &
      mesh_class, name, class, rank_write)

    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: name, class
    CHARACTER(LEN=*), INTENT(IN) :: vector_name, vector_class, mesh_name
    CHARACTER(LEN=*), INTENT(IN) :: mesh_class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: ndims, md_length, block_length
    INTEGER :: iloop

    ndims = SIZE(name)

    md_length = 2 * max_string_len + soi
    block_length = md_length + ndims * 2 * max_string_len

    CALL cfd_write_block_header(vector_name, vector_class, &
        c_type_stitched_vector, block_length, md_length, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, ndims, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO iloop = 1, ndims
        CALL cfd_safe_write_string(name(iloop))
        CALL cfd_safe_write_string(class(iloop))
      ENDDO
    ENDIF

    current_displacement = current_displacement + 2 * ndims * max_string_len

  END SUBROUTINE cfd_write_stitched_vector



  SUBROUTINE cfd_write_stitched_magnitude(magn_name, magn_class, mesh_name, &
      mesh_class, name, class, rank_write)

    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: name, class
    CHARACTER(LEN=*), INTENT(IN) :: magn_name, magn_class, mesh_name, mesh_class
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: ndims, md_length, block_length
    INTEGER :: iloop

    ndims = SIZE(name)

    md_length = 2 * max_string_len + soi
    block_length = md_length + ndims * 2 * max_string_len

    CALL cfd_write_block_header(magn_name, magn_class, &
        c_type_stitched_magnitude, block_length, md_length, rank_write)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL cfd_safe_write_string(mesh_name)
      CALL cfd_safe_write_string(mesh_class)
    ENDIF

    current_displacement = current_displacement + 2 * max_string_len

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) &
        CALL MPI_FILE_WRITE(cfd_filehandle, ndims, 1, MPI_INTEGER, &
            cfd_status, cfd_errcode)

    current_displacement = current_displacement + soi

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, &
        MPI_CHARACTER, MPI_CHARACTER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      DO iloop = 1, ndims
        CALL cfd_safe_write_string(name(iloop))
        CALL cfd_safe_write_string(class(iloop))
      ENDDO
    ENDIF

    current_displacement = current_displacement + 2 * ndims * max_string_len

  END SUBROUTINE cfd_write_stitched_magnitude



  SUBROUTINE cfd_write_real_constant(name, class, value, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    REAL(num), INTENT(IN) :: value
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length

    md_length = num

    CALL cfd_write_block_header(name, class, c_type_constant, md_length, &
        md_length, rank_write)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, &
        mpireal, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      CALL MPI_FILE_WRITE(cfd_filehandle, value, 1, mpireal, cfd_status, &
          cfd_errcode)
    ENDIF

    current_displacement = current_displacement + num

  END SUBROUTINE cfd_write_real_constant



  SUBROUTINE cfd_write_1d_integer_array(name, class, values, rank_write)

    CHARACTER(LEN=*), INTENT(IN) :: name, class
    INTEGER, DIMENSION(:), INTENT(IN) :: values
    INTEGER, INTENT(IN) :: rank_write
    INTEGER(8) :: md_length

    md_length = 2 * soi

    CALL cfd_write_block_header(name, class, c_type_integerarray, md_length, &
        md_length, rank_write)
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, &
        MPI_INTEGER, "native", MPI_INFO_NULL, cfd_errcode)

    IF (cfd_rank .EQ. rank_write) THEN
      ! 1D
      CALL MPI_FILE_WRITE(cfd_filehandle, 1, 1, MPI_INTEGER, &
          cfd_status, cfd_errcode)
      ! Size of array
      CALL MPI_FILE_WRITE(cfd_filehandle, 1, SIZE(values), MPI_INTEGER, &
          cfd_status, cfd_errcode)
      ! Actual array
      CALL MPI_FILE_WRITE(cfd_filehandle, values, SIZE(values), MPI_INTEGER, &
          cfd_status, cfd_errcode)
    ENDIF

    current_displacement = current_displacement + md_length

  END SUBROUTINE cfd_write_1d_integer_array



  SUBROUTINE cfd_write_visit_expression(expression_name, expression_class, &
      expression)

    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: expression_name
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: expression_class, expression

    PRINT *, LEN(expression(1)), LEN(expression(2))

  END SUBROUTINE cfd_write_visit_expression

END MODULE output
