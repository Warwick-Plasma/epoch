MODULE deck_species_block

  USE strings_advanced
  USE setup

  IMPLICIT NONE

  SAVE

  INTEGER :: species_loaded, species_id

CONTAINS

  FUNCTION handle_species_deck(element_in, value)

    CHARACTER(*), INTENT(IN) :: element_in, value
    CHARACTER(LEN=string_length) :: element
    INTEGER :: handle_species_deck
    LOGICAL :: handled

    handle_species_deck = c_err_none
    IF (value .EQ. blank .OR. element .EQ. blank) RETURN

    handle_species_deck = c_err_unknown_element

    handled = .FALSE.
    IF (str_cmp(element_in, "n_species")) THEN
      n_species = as_integer(value, handle_species_deck)
      IF (n_species .GT. 0) THEN
        IF (rank .EQ. 0) &
            PRINT '("Code running with ", i2, " species")', n_species
        CALL setup_species
      ENDIF
      handle_species_deck = c_err_none
      handled = .TRUE.
    ENDIF

    IF (n_species .LE. 0) THEN
      IF (rank .EQ. 0) THEN
        PRINT *, "Invalid number of species specified"
      ENDIF
      handle_species_deck = c_err_bad_value
    ENDIF

    IF (.NOT. ASSOCIATED(particle_species)) THEN
      extended_error_string = "n_species"
      IF (rank .EQ. 0) THEN
        PRINT *, "Attempting to set species data before setting n_species"
      ENDIF
      handle_species_deck = c_err_required_element_not_set
      RETURN
    ENDIF

    IF (handled) RETURN

    CALL split_off_int(element_in, element, species_id, handle_species_deck)

    IF (species_id .LT. 1 .OR. species_id .GT. n_species) THEN
      IF (rank .EQ. 0) &
          PRINT '("Ignoring attempt to set property ", a, &
              &" for non existent species ", i2)', TRIM(ADJUSTL(element)), &
              species_id
      handle_species_deck = c_err_none
      RETURN
    ENDIF

    IF (str_cmp(element, "name")) THEN
      particle_species(species_id)%name = TRIM(value)
      IF (rank .EQ. 0) &
          PRINT '("Name of species ", i3, " is ", a)', species_id, TRIM(value)
      handle_species_deck = c_err_none
      RETURN
    ENDIF

    IF (str_cmp(element, "mass")) THEN
      handle_species_deck = c_err_none
      particle_species(species_id)%mass = &
          as_real(value, handle_species_deck) * m0
      RETURN
    ENDIF

    IF (str_cmp(element, "charge")) THEN
      handle_species_deck = c_err_none
      particle_species(species_id)%charge = &
          as_real(value, handle_species_deck) * q0
      RETURN
    ENDIF

    IF (str_cmp(element, "frac") .OR. str_cmp(element, "fraction")) THEN
      handle_species_deck = c_err_none
      IF (npart_global .GE. 0) THEN
        particle_species(species_id)%count = &
            INT(as_real(value, handle_species_deck) * npart_global)
      ELSE
        extended_error_string = "npart"
        handle_species_deck = c_err_required_element_not_set
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "npart")) THEN
      handle_species_deck = c_err_none
      particle_species(species_id)%count = &
          as_long_integer(value, handle_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "dump")) THEN
      handle_species_deck = c_err_none
      particle_species(species_id)%dump = &
          as_logical(value, handle_species_deck)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for tracer particles
    ! *************************************************************
    IF (str_cmp(element, "tracer")) THEN
      handle_species_deck = c_err_none
#ifdef TRACER_PARTICLES
      particle_species(species_id)%tracer = &
          as_logical(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DTRACER_PARTICLES"
#endif
    ENDIF

    ! *************************************************************
    ! This section sets properties for particle splitting
    ! *************************************************************
    IF (str_cmp(element, "split")) THEN
      handle_species_deck = c_err_none
#ifdef SPLIT_PARTICLES_AFTER_PUSH
      particle_species(species_id)%split = &
          as_logical(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DSPLIT_PARTICLES_AFTER_PUSH"
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, "npart_max")) THEN
      handle_species_deck = c_err_none
#ifdef SPLIT_PARTICLES_AFTER_PUSH
      particle_species(species_id)%npart_max = &
          as_long_integer(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DSPLIT_PARTICLES_AFTER_PUSH"
#endif
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for ionisation
    ! *************************************************************
    IF (str_cmp(part1, "ionise")) THEN
      handle_species_deck = c_err_none
#ifdef PARTICLE_IONISE
      particle_species(part2)%ionise = as_logical(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

    IF (str_cmp(part1, "ionise_to_species")) THEN
      handle_species_deck = c_err_none
#ifdef PARTICLE_IONISE
      particle_species(part2)%ionise_to_species = &
          as_integer(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

    IF (str_cmp(part1, "release_species_on_ionise")) THEN
      handle_species_deck = c_err_none
#ifdef PARTICLE_IONISE
      particle_species(part2)%release_species = &
          as_integer(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

    IF (str_cmp(part1, "ionisation_energy")) THEN
      handle_species_deck = c_err_none
#ifdef PARTICLE_IONISE
      particle_species(part2)%ionisation_energy = &
          as_real(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

  END FUNCTION handle_species_deck



  FUNCTION check_species_block()

    INTEGER :: check_species_block

    ! Should do error checking but isn't really necessary at the moment
    check_species_block = c_err_none

  END FUNCTION check_species_block

END MODULE deck_species_block
