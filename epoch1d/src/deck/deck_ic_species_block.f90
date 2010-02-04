MODULE deck_ic_species_block

  USE shared_data
  USE strings_advanced
  USE shared_parser_data
  USE strings
  USE shunt
  USE evaluator

  IMPLICIT NONE

  SAVE

  INTEGER :: species_loaded

CONTAINS

  FUNCTION handle_ic_species_deck(species_id, element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER, INTENT(IN) :: species_id
    INTEGER :: handle_ic_species_deck

    handle_ic_species_deck = c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (species_id .LT. 0 .OR. species_id .GT. n_species) THEN
      IF (rank .EQ. 0) PRINT *, "Attempting to set non-existent species initial conditions. Ignoring."
      RETURN
    ENDIF

    IF (str_cmp(element, "minrho")) THEN
      initial_conditions(species_id)%minrho = as_real(value, handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "maxrho")) THEN
      initial_conditions(species_id)%maxrho = as_real(value, handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "rho") .OR. str_cmp(element, "number_density")) THEN
      CALL evaluate_string_in_space(value, initial_conditions(species_id)%rho(-2:nx+3), (/-2, nx+3/), handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "mass_density")) THEN
      initial_conditions(species_id)%rho = &
          initial_conditions(species_id)%rho / particle_species(species_id)%mass
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_x")) THEN
      initial_conditions(species_id)%drift(1) = &
          as_real(value, handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_y")) THEN
      initial_conditions(species_id)%drift(2) = &
          as_real(value, handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_z")) THEN
      initial_conditions(species_id)%drift(3) = &
          as_real(value, handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp")) THEN
      CALL evaluate_string_in_space(value, initial_conditions(species_id)%temp(-1:nx+2, 1), (/-1, nx+2/), handle_ic_species_deck)
      initial_conditions(species_id)%temp(-1:nx+2, 2) = initial_conditions(species_id)%temp(-1:nx+2, 1)
      initial_conditions(species_id)%temp(-1:nx+2, 3) = initial_conditions(species_id)%temp(-1:nx+2, 1)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_x")) THEN
      CALL evaluate_string_in_space(value, initial_conditions(species_id)%temp(-1:nx+2, 1), (/-1, nx+2/), handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_y")) THEN
      CALL evaluate_string_in_space(value, initial_conditions(species_id)%temp(-1:nx+2, 2), (/-1, nx+2/), handle_ic_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_z")) THEN
      CALL evaluate_string_in_space(value, initial_conditions(species_id)%temp(-1:nx+2, 3), (/-1, nx+2/), handle_ic_species_deck)
      RETURN
    ENDIF

  END FUNCTION handle_ic_species_deck



  FUNCTION check_ic_species_block()

    INTEGER :: check_ic_species_block

    ! Should do error checking but can't be bothered at the moment
    check_ic_species_block = c_err_none

  END FUNCTION check_ic_species_block

END MODULE deck_ic_species_block
