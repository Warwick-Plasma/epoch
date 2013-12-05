MODULE sdf_input_station_ru

  USE sdf_input_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_read_station_info(h, nelements, type_size, nstations, &
      nvariables, step, step_increment, time, time_increment, use_mult, &
      variable_ids, variable_names, station_ids)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT), OPTIONAL :: nelements
    INTEGER, INTENT(OUT), OPTIONAL :: type_size, nstations, nvariables
    INTEGER, INTENT(OUT), OPTIONAL :: step, step_increment
    REAL(r8), INTENT(OUT), OPTIONAL :: time, time_increment
    LOGICAL, INTENT(OUT), OPTIONAL :: use_mult
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: &
        variable_ids(:), variable_names(:), station_ids(:)
    INTEGER :: i
    CHARACTER(LEN=4) :: padding
    TYPE(sdf_block_type), POINTER :: b

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

    IF (sdf_info_init(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) THEN
      CALL read_entry_int8(h, b%nelements)
      CALL read_entry_int4(h, b%type_size)
      CALL read_entry_int4(h, b%nstations)
      CALL read_entry_int4(h, b%nvariables)
      CALL read_entry_int4(h, b%step)
      CALL read_entry_int4(h, b%step_increment)
      CALL read_entry_real8(h, b%time)
      CALL read_entry_real8(h, b%time_increment)
      CALL read_entry_logical(h, b%use_mult)
      CALL read_entry_stringlen(h, padding, 3)

      ALLOCATE(b%station_ids(b%nstations))
      ALLOCATE(b%station_names(b%nstations))
      ALLOCATE(b%station_nvars(b%nstations))
      ALLOCATE(b%station_move(b%nstations))
      ALLOCATE(b%station_grid(b%nstations,b%ndims))
      ALLOCATE(b%variable_ids(b%nvariables))
      ALLOCATE(b%dim_units(b%nvariables))
      ALLOCATE(b%material_names(b%nvariables))
      ALLOCATE(b%variable_types(b%nvariables))

      DO i = 1,b%nstations
        CALL read_entry_id(h, b%station_ids(i))
      ENDDO
      DO i = 1,b%nstations
        CALL read_entry_string(h, b%station_names(i))
      ENDDO
      CALL read_entry_array_int4(h, b%station_nvars, b%nstations)
      CALL read_entry_array_int4(h, b%station_move, b%nstations)
      DO i = 1,b%ndims
        CALL read_entry_array_real8(h, b%station_grid(:,i), b%nstations)
      ENDDO
      DO i = 1,b%nvariables
        CALL read_entry_id(h, b%variable_ids(i))
      ENDDO
      DO i = 1,b%nvariables
        CALL read_entry_string(h, b%material_names(i))
      ENDDO
      CALL read_entry_array_int4(h, b%variable_types, b%nvariables)
      DO i = 1,b%nvariables
        CALL read_entry_id(h, b%dim_units(i))
      ENDDO

      IF (b%use_mult) THEN
        ALLOCATE(b%dim_mults(b%nvariables))
        CALL read_entry_array_real8(h, b%dim_mults, b%nvariables)
      ENDIF
    ENDIF

    IF (PRESENT(nelements     )) nelements      = b%nelements
    IF (PRESENT(type_size     )) type_size      = b%type_size
    IF (PRESENT(nstations     )) nstations      = b%nstations
    IF (PRESENT(nvariables    )) nvariables     = b%nvariables
    IF (PRESENT(step          )) step           = b%step
    IF (PRESENT(step_increment)) step_increment = b%step_increment
    IF (PRESENT(time          )) time           = b%time
    IF (PRESENT(time_increment)) time_increment = b%time_increment
    IF (PRESENT(use_mult      )) use_mult       = b%use_mult
    IF (PRESENT(variable_ids)) THEN
       i = MIN(SIZE(variable_ids),b%nvariables)
       variable_ids(1:i) = b%variable_ids(1:i)
    ENDIF
    IF (PRESENT(variable_names)) THEN
       i = MIN(SIZE(variable_names),b%nvariables)
       variable_names(1:i) = b%material_names(1:i)
    ENDIF
    IF (PRESENT(station_ids)) THEN
       i = MIN(SIZE(station_ids),b%nvariables)
       station_ids(1:i) = b%station_ids(1:i)
    ENDIF

    h%current_location = b%data_location
    b%done_info = .TRUE.

  END SUBROUTINE sdf_read_station_info



  SUBROUTINE sdf_read_station_info_arrays(h, station_move, station_nvars)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT), OPTIONAL :: station_move(:)
    INTEGER, INTENT(OUT), OPTIONAL :: station_nvars(:)
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block
    IF (PRESENT(station_move)) &
        station_move(1:b%nstations) = b%station_move(1:b%nstations)
    IF (PRESENT(station_nvars)) &
        station_nvars(1:b%nstations) = b%station_nvars(1:b%nstations)

  END SUBROUTINE sdf_read_station_info_arrays



  SUBROUTINE sdf_read_station_info_arrays_all(h, station_move, station_nvars)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT), OPTIONAL :: station_move(:)
    INTEGER, INTENT(OUT), OPTIONAL :: station_nvars(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block
    IF (PRESENT(station_move)) THEN
      station_move = -1
      DO i = 1,b%nstations
        station_move(b%station_index(i)) = b%station_move(i)
      ENDDO
    ENDIF
    IF (PRESENT(station_nvars)) THEN
      station_nvars = 0
      DO i = 1,b%nstations
        station_nvars(b%station_index(i)) = b%station_nvars(i)
      ENDDO
    ENDIF

  END SUBROUTINE sdf_read_station_info_arrays_all

END MODULE sdf_input_station_ru
