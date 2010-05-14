MODULE deck_ic_species_block

  USE strings_advanced
  USE simple_io

  IMPLICIT NONE

  SAVE

  INTEGER :: species_loaded
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0

CONTAINS

  SUBROUTINE species_start

    offset = 0

  END SUBROUTINE species_start



  FUNCTION handle_ic_species_deck(species_id, element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER, INTENT(IN) :: species_id
    INTEGER :: handle_ic_species_deck
    CHARACTER(LEN=string_length) :: filename
    INTEGER :: err
    LOGICAL :: got_file

    handle_ic_species_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (species_id .LT. 0 .OR. species_id .GT. n_species) THEN
      IF (rank .EQ. 0) &
          PRINT *, "Attempting to set non-existent species initial &
              &conditions. Ignoring."
      RETURN
    ENDIF

    IF (str_cmp(element, "offset")) THEN
      offset = as_long_integer_simple(value, handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "minrho")) THEN
      initial_conditions(species_id)%minrho = &
          as_real(value, handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "maxrho")) THEN
      initial_conditions(species_id)%maxrho = &
          as_real(value, handle_ic_species_deck)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, err)

    IF (str_cmp(element, "rho") .OR. str_cmp(element, "number_density")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%rho, offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%rho(-2:nx+3, -2:ny+3), &
            (/-2, nx+3/), (/-2, ny+3/), handle_ic_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "mass_density")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%rho, offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%rho(-2:nx+3, -2:ny+3), &
            (/-2, nx+3/), (/-2, ny+3/), handle_ic_species_deck)
      ENDIF
      initial_conditions(species_id)%rho = &
          initial_conditions(species_id)%rho / particle_species(species_id)%mass
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_x")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%drift(:,:,1), offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(-1:nx+2, -1:ny+2, 1), &
            (/-1, nx+2/), (/-1, ny+2/), handle_ic_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_y")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%drift(:,:,2), offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(-1:nx+2, -1:ny+2, 2), &
            (/-1, nx+2/), (/-1, ny+2/), handle_ic_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_z")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%drift(:,:,3), offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(-1:nx+2, -1:ny+2, 3), &
            (/-1, nx+2/), (/-1, ny+2/), handle_ic_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%temp(:,:,1), offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 1), &
            (/-1, nx+2/), (/-1, ny+2/), handle_ic_species_deck)
      ENDIF
      debug_mode = .FALSE.
      initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 2) = &
          initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 1)
      initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 3) = &
          initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 1)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_x")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%temp(:,:,1), offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 1), &
            (/-1, nx+2/), (/-1, ny+2/), handle_ic_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_y")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%temp(:,:,2), offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 2), &
            (/-1, nx+2/), (/-1, ny+2/), handle_ic_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_z")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_data_file(filename, &
            initial_conditions(species_id)%temp(:,:,3), offset, &
            handle_ic_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(-1:nx+2, -1:ny+2, 3), &
            (/-1, nx+2/), (/-1, ny+2/), handle_ic_species_deck)
      ENDIF
      RETURN
    ENDIF

  END FUNCTION handle_ic_species_deck



  FUNCTION check_ic_species_block()

    INTEGER :: check_ic_species_block

    ! Should do error checking but can't be bothered at the moment
    check_ic_species_block = c_err_none

  END FUNCTION check_ic_species_block

END MODULE deck_ic_species_block
