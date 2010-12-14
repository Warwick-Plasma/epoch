MODULE sdf_input_util

  USE sdf_input
  USE sdf_input_cartesian
  USE sdf_input_point

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_read_blocklist(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: i, step, code_io_version, blocktype, ndims, datatype
    REAL(num) :: time
    CHARACTER(LEN=c_id_length) :: code_name, id
    CHARACTER(LEN=512) :: name
    LOGICAL :: restart, other_domains
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: errcode, buflen

    NULLIFY(h%current_block)
    IF (.NOT. h%done_header) CALL sdf_read_header(h)
    IF (ASSOCIATED(h%blocklist)) RETURN

    buflen = h%summary_size
    ALLOCATE(h%buffer(buflen))

    IF (h%rank .EQ. h%rank_master) THEN
      h%current_location = h%summary_location
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      CALL MPI_FILE_READ(h%filehandle, h%buffer, buflen, &
          MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL MPI_BCAST(h%buffer, buflen, MPI_CHARACTER, h%rank_master, &
        h%comm, errcode)

    h%start_location = h%summary_location
    DO i = 1,h%nblocks
      CALL sdf_read_block_info(h)
    ENDDO

    DEALLOCATE(h%buffer)
    NULLIFY(h%current_block)

  END SUBROUTINE sdf_read_blocklist



  SUBROUTINE sdf_read_block_info(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_read_next_block_header(h)
    b => h%current_block

    IF (b%blocktype .EQ. c_blocktype_plain_mesh) THEN
      CALL sdf_read_plain_mesh_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_point_mesh) THEN
      CALL sdf_read_point_mesh_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_plain_variable) THEN
      CALL sdf_read_plain_variable_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_point_variable) THEN
      CALL sdf_read_point_variable_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_constant) THEN
      CALL sdf_read_constant(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_array) THEN
      CALL sdf_read_array_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_run_info) THEN
      CALL sdf_read_run_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_tensor) THEN
      CALL sdf_read_stitched_tensor(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_material) THEN
      CALL sdf_read_stitched_material(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_matvar) THEN
      CALL sdf_read_stitched_matvar(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_species) THEN
      CALL sdf_read_stitched_species(h)
    ENDIF

  END SUBROUTINE sdf_read_block_info

END MODULE sdf_input_util
