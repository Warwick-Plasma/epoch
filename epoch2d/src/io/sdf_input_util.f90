MODULE sdf_input_util

  USE sdf_input
  USE sdf_input_cartesian
  USE sdf_input_point
  USE sdf_input_station
  USE sdf_output_station_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_read_blocklist(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: errcode, buflen, i

    NULLIFY(h%current_block)
    IF (.NOT. h%done_header) CALL sdf_read_header(h)
    IF (ASSOCIATED(h%blocklist)) RETURN

    h%current_location = h%summary_location

    buflen = h%summary_size
    IF (buflen .GT. 0) THEN
      ALLOCATE(h%buffer(buflen))

      IF (h%rank .EQ. h%rank_master) THEN
        CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, h%buffer, &
            buflen, MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)
      ENDIF

      CALL MPI_BCAST(h%buffer, buflen, MPI_CHARACTER, h%rank_master, &
          h%comm, errcode)
    ENDIF

    h%start_location = h%summary_location
    DO i = 1,h%nblocks
      CALL sdf_read_next_block_info(h)
    ENDDO

    IF (ASSOCIATED(h%buffer)) DEALLOCATE(h%buffer)
    NULLIFY(h%current_block)

  END SUBROUTINE sdf_read_blocklist



  SUBROUTINE sdf_read_next_block_info(h)

    TYPE(sdf_file_handle) :: h

    CALL sdf_read_next_block_header(h)
    CALL sdf_read_block_info(h)

  END SUBROUTINE sdf_read_next_block_info



  SUBROUTINE sdf_read_block_info(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b

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
    ELSE IF (b%blocktype .EQ. c_blocktype_cpu_split) THEN
      CALL sdf_read_cpu_split_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_run_info) THEN
      CALL sdf_read_run_info(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched &
        .OR. b%blocktype .EQ. c_blocktype_contiguous &
        .OR. b%blocktype .EQ. c_blocktype_stitched_tensor &
        .OR. b%blocktype .EQ. c_blocktype_contiguous_tensor) THEN
      CALL sdf_read_stitched(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_material &
        .OR. b%blocktype .EQ. c_blocktype_contiguous_material) THEN
      CALL sdf_read_stitched_material(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_matvar &
        .OR. b%blocktype .EQ. c_blocktype_contiguous_matvar) THEN
      CALL sdf_read_stitched_matvar(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_species &
        .OR. b%blocktype .EQ. c_blocktype_contiguous_species) THEN
      CALL sdf_read_stitched_species(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_stitched_obstacle_group) THEN
      CALL sdf_read_stitched_obstacle_group(h)
    ELSE IF (b%blocktype .EQ. c_blocktype_station) THEN
      CALL sdf_read_station_info(h)
    ENDIF

  END SUBROUTINE sdf_read_block_info



  SUBROUTINE sdf_read_stitched_info(h, variable_ids)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: variable_ids
    INTEGER :: iloop
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    CALL sdf_read_stitched(h)
    b => h%current_block

    DO iloop = 1, b%ndims
      CALL safe_copy_string(b%variable_ids(iloop), variable_ids(iloop))
    ENDDO

  END SUBROUTINE sdf_read_stitched_info



  FUNCTION sdf_station_seek_time(h, time) RESULT(found)

    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(IN) :: time
    LOGICAL :: found
    INTEGER :: i, ns, nstep, errcode, mpireal
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    REAL(r4) :: time4
    TYPE(sdf_block_type), POINTER :: b, station_block

    IF (.NOT.ASSOCIATED(h%blocklist)) CALL sdf_read_blocklist(h)

    found = .FALSE.
    b => h%blocklist
    DO i = 1,h%nblocks
      IF (b%blocktype .EQ. c_blocktype_station) THEN
        station_block => b
        found = .TRUE.
      ENDIF
      b => b%next_block
    ENDDO

    IF (.NOT.found) RETURN

    found = .TRUE.

    b => station_block
    ns = INT((time - b%time) * b%nelements / (h%time - b%time))
    IF (ns .GE. b%nelements) THEN
      ns = INT(b%nelements-1)
    ELSE IF (ns .LT. 0) THEN
      ns = 0
    ENDIF
    nstep = ns

    offset = b%data_location + ns * b%type_size
    mpireal = MPI_REAL4

    CALL MPI_FILE_READ_AT(h%filehandle, offset, time4, 1, mpireal, &
        MPI_STATUS_IGNORE, errcode)
    IF (time4 .GT. time) THEN
      DO i = ns-1, 0, -1
        offset = offset - b%type_size
        CALL MPI_FILE_READ_AT(h%filehandle, offset, time4, 1, mpireal, &
            MPI_STATUS_IGNORE, errcode)
        nstep = i
        IF (time4 .LT. time) EXIT
      ENDDO
    ELSE
      DO i = ns+1, INT(b%nelements)
        offset = offset + b%type_size
        CALL MPI_FILE_READ_AT(h%filehandle, offset, time4, 1, mpireal, &
            MPI_STATUS_IGNORE, errcode)
        IF (time4 .GT. time) EXIT
        nstep = i
      ENDDO
    ENDIF

    h%current_block => b
    b%nelements = nstep
    b%data_length = b%nelements * b%type_size
    h%step = b%step + nstep
    h%time = time4

    CALL write_station_update(h)

  END FUNCTION sdf_station_seek_time

END MODULE sdf_input_util
