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



  FUNCTION sdf_station_seek_time(h, time, overwrite_in) RESULT(found)

    ! overwrite_in .EQ. .TRUE. if we intend to write new data after
    ! after the current time step. Used for restarting a run.
    !
    ! overwrite_in .EQ. .FALSE. otherwise. Used for quickly finding
    ! a timestep to read.

    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(IN) :: time
    LOGICAL, INTENT(IN), OPTIONAL :: overwrite_in
    LOGICAL :: found
    INTEGER :: i, ns, nstep, errcode, mpireal
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    REAL(r8) :: time8
    REAL(r4) :: real4
    TYPE(sdf_block_type), POINTER :: b, station_block
    LOGICAL :: overwrite

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
    mpireal = MPI_REAL8
    IF (b%variable_types(1) .EQ. c_datatype_real4) mpireal = MPI_REAL4

    CALL MPI_FILE_READ_AT(h%filehandle, offset, time8, 1, mpireal, &
        MPI_STATUS_IGNORE, errcode)

    IF (mpireal .EQ. MPI_REAL4) time8 = REAL(TRANSFER(time8, real4), r8)

    IF (time8 .GT. time) THEN
      DO i = ns-1, 0, -1
        offset = offset - b%type_size
        CALL MPI_FILE_READ_AT(h%filehandle, offset, time8, 1, mpireal, &
            MPI_STATUS_IGNORE, errcode)
        IF (mpireal .EQ. MPI_REAL4) time8 = REAL(TRANSFER(time8, real4), r8)
        nstep = i
        IF (time8 .LE. time) EXIT
      ENDDO
    ELSE IF (time8 .LT. time) THEN
      DO i = ns+1, INT(b%nelements)
        offset = offset + b%type_size
        CALL MPI_FILE_READ_AT(h%filehandle, offset, time8, 1, mpireal, &
            MPI_STATUS_IGNORE, errcode)
        IF (mpireal .EQ. MPI_REAL4) time8 = REAL(TRANSFER(time8, real4), r8)
        IF (time8 .GT. time) EXIT
        nstep = i
      ENDDO
    ENDIF

    h%current_block => b
    overwrite = .TRUE.
    IF (PRESENT(overwrite_in)) overwrite = overwrite_in

    IF (overwrite) THEN
       b%nelements = nstep + 1
       b%data_length = b%nelements * b%type_size
       h%step = b%step + nstep
       h%time = time8
       CALL write_station_update(h)
    ELSE
       CALL MPI_FILE_SEEK(h%filehandle, offset, MPI_SEEK_SET, errcode)
       h%current_location = offset
    ENDIF


  END FUNCTION sdf_station_seek_time



  SUBROUTINE sdf_get_all_stations(h, nstations, station_ids)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT), OPTIONAL :: nstations
    CHARACTER(LEN=c_id_length), POINTER, OPTIONAL :: station_ids(:)
    INTEGER :: i, m, n, nstat_max, nstat_list_max = 0, nextra
    TYPE(sdf_block_type), POINTER :: b, next
    CHARACTER(LEN=c_id_length), POINTER :: ctmp(:)
    LOGICAL :: found
    LOGICAL, ALLOCATABLE :: found_id(:)
    INTEGER, ALLOCATABLE :: extra_id(:)

    IF (.NOT.ASSOCIATED(h%blocklist)) CALL sdf_read_blocklist(h)

    IF (.NOT.ASSOCIATED(h%station_ids)) THEN
      nstat_max = 2
      ALLOCATE(extra_id(nstat_max))
      ALLOCATE(found_id(1))

      nstat_list_max = 0
      next => h%blocklist
      DO i = 1,h%nblocks
        h%current_block => next
        next => h%current_block%next_block
        IF (h%current_block%blocktype .NE. c_blocktype_station) CYCLE

        b => h%current_block
        CALL sdf_read_block_info(h)

        IF (b%nstations .GT. nstat_max) THEN
          nstat_max = b%nstations * 11 / 10 + 2
          DEALLOCATE(extra_id)
          ALLOCATE(extra_id(nstat_max))
        ENDIF

        ALLOCATE(b%station_index(b%nstations))
        found_id = .FALSE.
        nextra = 0
        DO n = 1, b%nstations
          found = .FALSE.
          DO m = 1, nstat_list_max
            IF (found_id(m)) CYCLE
            IF (b%station_ids(n) .NE. h%station_ids(m)) CYCLE
            b%station_index(n) = m
            found = .TRUE.
            found_id(m) = .TRUE.
            EXIT
          ENDDO
          IF (.NOT.found) THEN
            nextra = nextra + 1
            extra_id(nextra) = n
            b%station_index(n) = nstat_list_max + n
          ENDIF
        ENDDO

        IF (nextra .EQ. 0) CYCLE

        IF (nstat_list_max .EQ. 0) THEN
          ALLOCATE(h%station_ids(nstat_list_max+nextra))
        ELSE
          ALLOCATE(ctmp(nstat_list_max))
          ctmp(1:nstat_list_max) = h%station_ids(1:nstat_list_max)
          DEALLOCATE(h%station_ids)
          ALLOCATE(h%station_ids(nstat_list_max+nextra))
          h%station_ids(1:nstat_list_max) = ctmp(1:nstat_list_max)
          DEALLOCATE(ctmp)
        ENDIF

        DEALLOCATE(found_id)
        ALLOCATE(found_id(nstat_list_max+nextra))

        m = nstat_list_max
        DO n = 1,nextra
          m = m + 1
          h%station_ids(m) = b%station_ids(extra_id(n))
        ENDDO

        nstat_list_max = nstat_list_max + nextra
      ENDDO

      h%nstations = nstat_list_max

      DEALLOCATE(extra_id, found_id)
    ENDIF

    IF (PRESENT(nstations)) nstations = nstat_list_max

    IF (PRESENT(station_ids) .AND. h%nstations .GT. 0) THEN
      ALLOCATE(station_ids(h%nstations))
      DO n = 1,h%nstations
        station_ids(n) = h%station_ids(n)
      ENDDO
    ENDIF

  END SUBROUTINE sdf_get_all_stations

END MODULE sdf_input_util
