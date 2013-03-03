MODULE sdf_output_station_ru

  USE sdf_output_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_station_header_gen(ndims, h, id, name, step, &
      step_increment, time, time_increment, nstations, station_ids, &
      station_names, station_nvars, station_move, station_x, station_y, &
      station_z, variable_ids, variable_names, variable_types, variable_units, &
      variable_mults )

    INTEGER, INTENT(IN) :: ndims
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    INTEGER, INTENT(IN), OPTIONAL :: step, step_increment
    REAL(r8), INTENT(IN), OPTIONAL :: time, time_increment
    INTEGER, INTENT(IN), OPTIONAL :: nstations
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: station_ids(:), station_names(:)
    INTEGER, INTENT(IN), OPTIONAL :: station_nvars(:), station_move(:)
    REAL(r8), INTENT(IN), OPTIONAL :: station_x(:), station_y(:), station_z(:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: variable_ids(:), variable_names(:)
    INTEGER, INTENT(IN), OPTIONAL :: variable_types(:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: variable_units(:)
    REAL(r8), INTENT(IN), OPTIONAL :: variable_mults(:)
    INTEGER :: i, use_mult, errcode
    CHARACTER(LEN=1) :: use_mult_char
    CHARACTER(LEN=4) :: padding
    TYPE(sdf_block_type), POINTER :: b

    IF (nstations .LE. 0) RETURN

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%geometry = c_geometry_cartesian
    b%ndims = ndims
    b%blocktype = c_blocktype_station
    b%nelements = 0
    b%datatype = c_datatype_other
    h%station_file = .TRUE.

    IF (PRESENT(id)) THEN
      b%nstations = nstations

      b%step = step
      b%time = time
      b%step_increment = step_increment
      b%time_increment = time_increment

      ALLOCATE(b%station_ids(b%nstations))
      ALLOCATE(b%station_names(b%nstations))
      ALLOCATE(b%station_nvars(b%nstations))
      ALLOCATE(b%station_move(b%nstations))
      ALLOCATE(b%station_grid(b%nstations,b%ndims))

      b%nvariables = 0
      IF (b%step_increment .NE. 1) b%nvariables = b%nvariables + 1
      IF (b%time_increment .NE. 1.d0) b%nvariables = b%nvariables + 1
      DO i = 1,b%nstations
        CALL safe_copy_id(h, station_ids(i), b%station_ids(i))
        CALL safe_copy_string(station_names(i), b%station_names(i))
        b%station_nvars(i) = station_nvars(i)
        b%station_move(i) = station_move(i)
        b%station_grid(i,1) = station_x(i)
        IF (b%ndims .GT. 1) b%station_grid(i,2) = station_y(i)
        IF (b%ndims .GT. 2) b%station_grid(i,3) = station_z(i)
        b%nvariables = b%nvariables + b%station_nvars(i)
      ENDDO

      b%use_mult = .FALSE.
      IF (PRESENT(variable_mults)) b%use_mult = .TRUE.

      ALLOCATE(b%variable_ids(b%nvariables))
      ALLOCATE(b%dim_units(b%nvariables))
      ALLOCATE(b%material_names(b%nvariables))
      ALLOCATE(b%variable_types(b%nvariables))
      IF (b%use_mult) ALLOCATE(b%dim_mults(b%nvariables))

      b%type_size = 0
      DO i = 1,b%nvariables
        CALL safe_copy_id(h, variable_ids(i), b%variable_ids(i))
        CALL safe_copy_id(h, variable_units(i), b%dim_units(i))
        CALL safe_copy_string(variable_names(i), b%material_names(i))
        b%variable_types(i) = variable_types(i)
        IF (b%use_mult) b%dim_mults(i) = variable_mults(i)
        b%type_size = b%type_size + c_type_sizes(b%variable_types(i))
      ENDDO
    ENDIF

    use_mult = 0
    IF (b%use_mult) use_mult = 1

    ! Metadata is
    ! - nelements INTEGER(i8)
    ! - entry_len INTEGER(i4)
    ! - nstations INTEGER(i4)
    ! - nvars     INTEGER(i4)
    ! - step0     INTEGER(i4)
    ! - step_inc  INTEGER(i4)
    ! - time0     REAL(r8)
    ! - time_inc  REAL(r8)
    ! - use_mult  CHARACTER(1)
    ! - padding   CHARACTER(3)
    ! - statids   CHARACTER(id_length), DIMENSION(nstations)
    ! - statnames CHARACTER(string_length), DIMENSION(nstations)
    ! - statnvars INTEGER(i4), DIMENSION(nstations)
    ! - statmove  INTEGER(i4), DIMENSION(nstations)
    ! - statx0    REAL(r8), DIMENSION(nstations*ndims)
    ! - varids    CHARACTER(id_length), DIMENSION(nvars)
    ! - varnames  CHARACTER(string_length), DIMENSION(nvars)
    ! - vartypes  INTEGER(i4), DIMENSION(nvars)
    ! - varunits  CHARACTER(id_length), DIMENSION(nvars)
    ! - varmults  REAL(r8), DIMENSION(use_mult*nvars)

    b%info_length = h%block_header_length + soi8 &
        + (5 + 2*b%nstations + b%nvariables) * soi4 &
        + (2 + b%ndims*b%nstations + use_mult*b%nvariables) * sof8 &
        + (b%nstations + 2*b%nvariables) * c_id_length &
        + (b%nstations + b%nvariables) * h%string_length + 4
    b%data_length = b%nelements * b%type_size

    IF (h%rank .EQ. h%rank_master) THEN
      IF (PRESENT(id)) THEN
        CALL sdf_write_block_header(h, id, name)
      ELSE
        CALL write_block_header(h)
      ENDIF

      CALL MPI_FILE_WRITE(h%filehandle, b%nelements, 1, MPI_INTEGER8, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%type_size, 1, MPI_INTEGER, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%nstations, 1, MPI_INTEGER, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%nvariables, 1, MPI_INTEGER, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%step, 1, MPI_INTEGER, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%step_increment, 1, MPI_INTEGER, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%time, 1, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%time_increment, 1, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)
      use_mult_char = ACHAR(use_mult)
      CALL MPI_FILE_WRITE(h%filehandle, use_mult_char, 1, MPI_CHARACTER, &
          MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, padding, 3, MPI_CHARACTER, &
          MPI_STATUS_IGNORE, errcode)
      DO i = 1,b%nstations
        CALL sdf_safe_write_id(h, b%station_ids(i))
      ENDDO
      DO i = 1,b%nstations
        CALL sdf_safe_write_string(h, b%station_names(i))
      ENDDO
      CALL MPI_FILE_WRITE(h%filehandle, b%station_nvars, nstations, &
          MPI_INTEGER, MPI_STATUS_IGNORE, errcode)
      CALL MPI_FILE_WRITE(h%filehandle, b%station_move, nstations, &
          MPI_INTEGER, MPI_STATUS_IGNORE, errcode)
      DO i = 1,b%ndims
        CALL MPI_FILE_WRITE(h%filehandle, b%station_grid(1,i), nstations, &
            MPI_REAL8, MPI_STATUS_IGNORE, errcode)
      ENDDO
      DO i = 1,b%nvariables
        CALL sdf_safe_write_id(h, b%variable_ids(i))
      ENDDO
      DO i = 1,b%nvariables
        CALL sdf_safe_write_string(h, b%material_names(i))
      ENDDO
      CALL MPI_FILE_WRITE(h%filehandle, b%variable_types, b%nvariables, &
          MPI_INTEGER, MPI_STATUS_IGNORE, errcode)
      DO i = 1,b%nvariables
        CALL sdf_safe_write_id(h, b%dim_units(i))
      ENDDO
      IF (b%use_mult) THEN
        CALL MPI_FILE_WRITE(h%filehandle, b%dim_mults, b%nvariables, &
            MPI_REAL8, MPI_STATUS_IGNORE, errcode)
      ENDIF
    ENDIF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_station_header_gen



  SUBROUTINE write_station_header_1d_r4(h, id, name, step, &
      step_increment, time, time_increment, nstations, station_ids, &
      station_names, station_nvars, station_move, station_x, &
      variable_ids, variable_names, variable_types, variable_units, &
      variable_mults)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: step, step_increment
    REAL(r4), INTENT(IN) :: time, time_increment
    INTEGER, INTENT(IN) :: nstations
    CHARACTER(LEN=*), INTENT(IN) :: station_ids(:), station_names(:)
    INTEGER, INTENT(IN) :: station_nvars(:), station_move(:)
    REAL(r4), INTENT(IN) :: station_x(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_ids(:), variable_names(:)
    INTEGER, INTENT(IN) :: variable_types(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_units(:)
    REAL(r4), INTENT(IN), OPTIONAL :: variable_mults(:)
    REAL(r8), ALLOCATABLE :: mults(:)
    REAL(r8) :: dummy(1)

    IF (PRESENT(variable_mults)) THEN
      ALLOCATE(mults(SIZE(variable_types)))
      mults = REAL(variable_mults,r8)
      CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
          REAL(time,r8), REAL(time_increment,r8), nstations, station_ids, &
          station_names, station_nvars, station_move, REAL(station_x,r8), &
          dummy, dummy, variable_ids, &
          variable_names, variable_types, variable_units, mults)
      DEALLOCATE(mults)
    ELSE
      CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
          REAL(time,r8), REAL(time_increment,r8), nstations, station_ids, &
          station_names, station_nvars, station_move, REAL(station_x,r8), &
          dummy, dummy, variable_ids, &
          variable_names, variable_types, variable_units)
    ENDIF

  END SUBROUTINE write_station_header_1d_r4



  SUBROUTINE write_station_header_2d_r4(h, id, name, step, &
      step_increment, time, time_increment, nstations, station_ids, &
      station_names, station_nvars, station_move, station_x, station_y, &
      variable_ids, variable_names, variable_types, variable_units, &
      variable_mults)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: step, step_increment
    REAL(r4), INTENT(IN) :: time, time_increment
    INTEGER, INTENT(IN) :: nstations
    CHARACTER(LEN=*), INTENT(IN) :: station_ids(:), station_names(:)
    INTEGER, INTENT(IN) :: station_nvars(:), station_move(:)
    REAL(r4), INTENT(IN) :: station_x(:), station_y(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_ids(:), variable_names(:)
    INTEGER, INTENT(IN) :: variable_types(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_units(:)
    REAL(r4), INTENT(IN), OPTIONAL :: variable_mults(:)
    REAL(r8), ALLOCATABLE :: mults(:)
    REAL(r8) :: dummy(1)

    IF (PRESENT(variable_mults)) THEN
      ALLOCATE(mults(SIZE(variable_types)))
      mults = REAL(variable_mults,r8)
      CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
          REAL(time,r8), REAL(time_increment,r8), nstations, station_ids, &
          station_names, station_nvars, station_move, REAL(station_x,r8), &
          REAL(station_y,r8), dummy, variable_ids, &
          variable_names, variable_types, variable_units, mults)
      DEALLOCATE(mults)
    ELSE
      CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
          REAL(time,r8), REAL(time_increment,r8), nstations, station_ids, &
          station_names, station_nvars, station_move, REAL(station_x,r8), &
          REAL(station_y,r8), dummy, variable_ids, &
          variable_names, variable_types, variable_units)
    ENDIF

  END SUBROUTINE write_station_header_2d_r4



  SUBROUTINE write_station_header_3d_r4(h, id, name, step, &
      step_increment, time, time_increment, nstations, station_ids, &
      station_names, station_nvars, station_move, station_x, station_y, &
      station_z, variable_ids, variable_names, variable_types, variable_units, &
      variable_mults)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: step, step_increment
    REAL(r4), INTENT(IN) :: time, time_increment
    INTEGER, INTENT(IN) :: nstations
    CHARACTER(LEN=*), INTENT(IN) :: station_ids(:), station_names(:)
    INTEGER, INTENT(IN) :: station_nvars(:), station_move(:)
    REAL(r4), INTENT(IN) :: station_x(:), station_y(:), station_z(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_ids(:), variable_names(:)
    INTEGER, INTENT(IN) :: variable_types(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_units(:)
    REAL(r4), INTENT(IN), OPTIONAL :: variable_mults(:)
    REAL(r8), ALLOCATABLE :: mults(:)

    IF (PRESENT(variable_mults)) THEN
      ALLOCATE(mults(SIZE(variable_types)))
      mults = REAL(variable_mults,r8)
      CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
          REAL(time,r8), REAL(time_increment,r8), nstations, station_ids, &
          station_names, station_nvars, station_move, REAL(station_x,r8), &
          REAL(station_y,r8), REAL(station_z,r8), variable_ids, &
          variable_names, variable_types, variable_units, mults)
      DEALLOCATE(mults)
    ELSE
      CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
          REAL(time,r8), REAL(time_increment,r8), nstations, station_ids, &
          station_names, station_nvars, station_move, REAL(station_x,r8), &
          REAL(station_y,r8), REAL(station_z,r8), variable_ids, &
          variable_names, variable_types, variable_units)
    ENDIF

  END SUBROUTINE write_station_header_3d_r4



  SUBROUTINE write_station_header_1d_r8(h, id, name, step, &
      step_increment, time, time_increment, nstations, station_ids, &
      station_names, station_nvars, station_move, station_x, &
      variable_ids, variable_names, variable_types, variable_units, &
      variable_mults)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: step, step_increment
    REAL(r8), INTENT(IN) :: time, time_increment
    INTEGER, INTENT(IN) :: nstations
    CHARACTER(LEN=*), INTENT(IN) :: station_ids(:), station_names(:)
    INTEGER, INTENT(IN) :: station_nvars(:), station_move(:)
    REAL(r8), INTENT(IN) :: station_x(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_ids(:), variable_names(:)
    INTEGER, INTENT(IN) :: variable_types(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_units(:)
    REAL(r8), INTENT(IN), OPTIONAL :: variable_mults(:)

    CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
        time, time_increment, nstations, station_ids, &
        station_names, station_nvars, station_move, station_x, &
        station_x, station_x, variable_ids, &
        variable_names, variable_types, variable_units, variable_mults)

  END SUBROUTINE write_station_header_1d_r8



  SUBROUTINE write_station_header_2d_r8(h, id, name, step, &
      step_increment, time, time_increment, nstations, station_ids, &
      station_names, station_nvars, station_move, station_x, station_y, &
      variable_ids, variable_names, variable_types, variable_units, &
      variable_mults)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: step, step_increment
    REAL(r8), INTENT(IN) :: time, time_increment
    INTEGER, INTENT(IN) :: nstations
    CHARACTER(LEN=*), INTENT(IN) :: station_ids(:), station_names(:)
    INTEGER, INTENT(IN) :: station_nvars(:), station_move(:)
    REAL(r8), INTENT(IN) :: station_x(:), station_y(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_ids(:), variable_names(:)
    INTEGER, INTENT(IN) :: variable_types(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_units(:)
    REAL(r8), INTENT(IN), OPTIONAL :: variable_mults(:)

    CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
        time, time_increment, nstations, station_ids, &
        station_names, station_nvars, station_move, station_x, &
        station_y, station_x, variable_ids, &
        variable_names, variable_types, variable_units, variable_mults)

  END SUBROUTINE write_station_header_2d_r8



  SUBROUTINE write_station_header_3d_r8(h, id, name, step, &
      step_increment, time, time_increment, nstations, station_ids, &
      station_names, station_nvars, station_move, station_x, station_y, &
      station_z, variable_ids, variable_names, variable_types, variable_units, &
      variable_mults)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    INTEGER, INTENT(IN) :: step, step_increment
    REAL(r8), INTENT(IN) :: time, time_increment
    INTEGER, INTENT(IN) :: nstations
    CHARACTER(LEN=*), INTENT(IN) :: station_ids(:), station_names(:)
    INTEGER, INTENT(IN) :: station_nvars(:), station_move(:)
    REAL(r8), INTENT(IN) :: station_x(:), station_y(:), station_z(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_ids(:), variable_names(:)
    INTEGER, INTENT(IN) :: variable_types(:)
    CHARACTER(LEN=*), INTENT(IN) :: variable_units(:)
    REAL(r8), INTENT(IN), OPTIONAL :: variable_mults(:)

    CALL write_station_header_gen(ndims, h, id, name, step, step_increment, &
        time, time_increment, nstations, station_ids, &
        station_names, station_nvars, station_move, station_x, &
        station_y, station_z, variable_ids, &
        variable_names, variable_types, variable_units, variable_mults)

  END SUBROUTINE write_station_header_3d_r8



  SUBROUTINE write_station_update(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: errcode
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset

    IF (h%rank .EQ. h%rank_master) THEN
      b => h%current_block
      b%next_block_location = b%data_location + b%data_length

      offset = b%block_start
      CALL MPI_FILE_WRITE_AT(h%filehandle, offset, b%next_block_location, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

      offset = b%block_start + 2*soi8 + c_id_length
      CALL MPI_FILE_WRITE_AT(h%filehandle, offset, b%data_length, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)

      offset = b%block_start + h%block_header_length
      CALL MPI_FILE_WRITE_AT(h%filehandle, offset, b%nelements, 1, &
          MPI_INTEGER8, MPI_STATUS_IGNORE, errcode)
    ENDIF

    CALL sdf_update(h)

  END SUBROUTINE write_station_update



  SUBROUTINE sdf_write_station_material(h, id, name, units, dims, nmat, &
      stagger, mesh_id, material_names)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    INTEGER :: i
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_write_stitched_material(h, id, name, mesh_id, stagger, &
        material_names, variable_ids, nmat)

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_station_material



  SUBROUTINE sdf_write_station_matvar(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_names)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    INTEGER :: i
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_write_stitched_matvar(h, id, name, mesh_id, stagger, &
        material_id, variable_ids, nmat)

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_station_matvar



  SUBROUTINE sdf_write_station_species(h, id, name, units, dims, nmat, &
      stagger, mesh_id, material_id, material_name, specnames)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_name
    CHARACTER(LEN=*), INTENT(IN) :: specnames(:)
    INTEGER :: i
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    ALLOCATE(variable_ids(nmat))

    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) .EQ. 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(specnames(i)), variable_ids(i))
      ENDIF
    ENDDO

    CALL sdf_write_stitched_species(h, id, name, mesh_id, stagger, &
        material_id, material_name, specnames, variable_ids, nmat)

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE sdf_write_station_species

END MODULE sdf_output_station_ru
