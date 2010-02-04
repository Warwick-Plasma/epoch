MODULE deck_ic_external_block

  USE shared_data
  USE strings_advanced
  USE mpi_subtype_control
  USE boundary

  IMPLICIT NONE

  SAVE

  INTEGER :: species_loaded

CONTAINS

  FUNCTION handle_ic_external_species_deck(species_id, element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER, INTENT(IN) :: species_id
    INTEGER :: handle_ic_external_species_deck

    handle_ic_external_species_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (species_id .LT. 0 .OR. species_id .GT. n_species) THEN
      IF (rank .EQ. 0) &
          PRINT *, "Attempting to set non-existent species initial &
              &conditions. Ignoring."
      RETURN
    ENDIF

    IF (str_cmp(element, "minrho")) THEN
      initial_conditions(species_id)%minrho = &
          as_real(value, handle_ic_external_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "maxrho")) THEN
      initial_conditions(species_id)%maxrho = &
          as_real(value, handle_ic_external_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "rho")) THEN
      CALL load_data_file(value, initial_conditions(species_id)%rho, &
          handle_ic_external_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp")) THEN
      CALL load_data_file(value, initial_conditions(species_id)%temp(:,:,:,1), &
          handle_ic_external_species_deck)
      initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, -1:nz+2, 2) = &
          initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, -1:nz+2, 1)
      initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, -1:nz+2, 3) = &
          initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, -1:nz+2, 1)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_x")) THEN
      CALL load_data_file(value, initial_conditions(species_id)%temp(:,:,:,1), &
          handle_ic_external_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_y")) THEN
      CALL load_data_file(value, initial_conditions(species_id)%temp(:,:,:,2), &
          handle_ic_external_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_z")) THEN
      CALL load_data_file(value, initial_conditions(species_id)%temp(:,:,:,3), &
          handle_ic_external_species_deck)
      RETURN
    ENDIF

  END FUNCTION handle_ic_external_species_deck



  FUNCTION handle_ic_external_fields_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_ic_external_fields_deck

    handle_ic_external_fields_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, "ex")) THEN
      CALL load_data_file(value, ex, handle_ic_external_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "ey")) THEN
      CALL load_data_file(value, ey, handle_ic_external_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "ez")) THEN
      CALL load_data_file(value, ez, handle_ic_external_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "bx")) THEN
      CALL load_data_file(value, bx, handle_ic_external_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "by")) THEN
      CALL load_data_file(value, by, handle_ic_external_fields_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "bz")) THEN
      CALL load_data_file(value, bz, handle_ic_external_fields_deck)
      RETURN
    ENDIF

  END FUNCTION handle_ic_external_fields_deck



  SUBROUTINE load_data_file(filename, array, err)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: subtype, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, &
        fh, errcode)
    IF (errcode .NE. 0) THEN
      IF (rank .EQ. 0) PRINT *, "file ", TRIM(filename), " does not exist."
      err = IOR(err, c_err_bad_value)
      RETURN
    ENDIF
    subtype = create_current_field_subtype()
    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype, "native", &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_READ_ALL(fh, array(1:nx, 1:ny, 1:nz), nx*ny*nz, mpireal, &
        status, errcode)
    CALL MPI_FILE_CLOSE(fh, errcode)
    CALL MPI_TYPE_FREE(subtype, errcode)

    CALL field_bc(array)
    CALL field_zero_gradient(array, .TRUE.)

  END SUBROUTINE load_data_file

END MODULE deck_ic_external_block
