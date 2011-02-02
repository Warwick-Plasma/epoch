MODULE deck_species_block

  USE strings_advanced
  USE setup
  USE simple_io
  USE utilities

  IMPLICIT NONE

  SAVE

  INTEGER :: species_id, current_block
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: species_names
  INTEGER, DIMENSION(:), POINTER :: species_blocks
  LOGICAL, DIMENSION(:), ALLOCATABLE :: species_charge_set
  LOGICAL :: got_name
  INTEGER :: check_block = c_err_none

CONTAINS

  SUBROUTINE species_initialise

    current_block = 0
    IF (deck_state .EQ. c_ds_deck) THEN
      n_species = 0
      ALLOCATE(species_names(4))
      ALLOCATE(species_blocks(4))
    ENDIF

  END SUBROUTINE species_initialise



  SUBROUTINE species_finalise

    INTEGER :: i
    CHARACTER(LEN=8) :: string

    IF (deck_state .EQ. c_ds_deck) THEN
      CALL setup_species
      ALLOCATE(species_charge_set(n_species))
      species_charge_set = .FALSE.
      DO i = 1, n_species
        particle_species(i)%name = species_names(i)
        IF (rank .EQ. 0) THEN
          CALL integer_as_string(i, string)
          PRINT*,'Name of species ', TRIM(ADJUSTL(string)), ' is ', &
              TRIM(species_names(i))
        ENDIF
      ENDDO
      IF (rank .EQ. 0) PRINT*
      DEALLOCATE(species_names)
    ENDIF
    ! Ugly hack.
    ! Deallocate species_blocks after the last call to read_deck.
    ! The last call depends on whether it is a restart dump or not.
    IF ((ic_from_restart .AND. deck_state .EQ. c_ds_eio) &
        .OR. deck_state .EQ. c_ds_ic) THEN
      DEALLOCATE(species_blocks)
      DEALLOCATE(species_charge_set)
    ENDIF

  END SUBROUTINE species_finalise



  SUBROUTINE species_start

    current_block = current_block + 1
    got_name = .FALSE.
    IF (deck_state .EQ. c_ds_deck) RETURN
    species_id = species_blocks(current_block)
    offset = 0

  END SUBROUTINE species_start



  SUBROUTINE species_end

    CHARACTER(LEN=8) :: id_string

    IF (.NOT.got_name) THEN
      IF (rank .EQ. 0) THEN
        CALL integer_as_string(current_block, id_string)
        PRINT*, '*** ERROR ***'
        PRINT*, 'Species block number ' // TRIM(id_string) &
            // ' has no "name" element.'
        WRITE(40,*) '*** ERROR ***'
        WRITE(40,*) 'Species block number ' // TRIM(id_string) &
            // ' has no "name" element.'
      ENDIF

      check_block = c_err_missing_elements
    ENDIF

  END SUBROUTINE species_end



  FUNCTION species_number_from_name(name)

    CHARACTER(*), INTENT(IN) :: name
    INTEGER :: species_number_from_name
    INTEGER :: i

    DO i = 1, n_species
      IF (str_cmp(name, species_names(i))) THEN
        species_number_from_name = i
        RETURN
      ENDIF
    ENDDO
    n_species = n_species + 1
    CALL grow_array(species_names, n_species)
    species_names(n_species) = TRIM(name)
    species_number_from_name = n_species

  END FUNCTION species_number_from_name



  FUNCTION handle_species_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_species_deck
    CHARACTER(LEN=string_length) :: filename
    LOGICAL :: got_file

    handle_species_deck = c_err_none
    IF (value .EQ. blank .OR. element .EQ. blank) RETURN

    IF (str_cmp(element, "name")) THEN
      IF (deck_state .NE. c_ds_deck) RETURN
      IF (got_name) THEN
        handle_species_deck = c_err_preset_element
        RETURN
      ENDIF
      CALL grow_array(species_blocks, current_block)
      species_blocks(current_block) = species_number_from_name(value)
      got_name = .TRUE.
      RETURN
    ENDIF

    IF (deck_state .EQ. c_ds_deck) RETURN

    IF (str_cmp(element, "mass")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
      particle_species(species_id)%mass = &
          as_real(value, handle_species_deck) * m0
      RETURN
    ENDIF

    IF (str_cmp(element, "charge")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
      particle_species(species_id)%charge = &
          as_real(value, handle_species_deck) * q0
      species_charge_set(species_id) = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, "frac") .OR. str_cmp(element, "fraction")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
      IF (npart_global .GE. 0) THEN
        particle_species(species_id)%count = &
            INT(as_real(value, handle_species_deck) * npart_global)
      ELSE
        particle_species(species_id)%count = 0
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "npart")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
      particle_species(species_id)%count = &
          as_long_integer(value, handle_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "dump")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
      particle_species(species_id)%dump = &
          as_logical(value, handle_species_deck)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for tracer particles
    ! *************************************************************
    IF (str_cmp(element, "tracer")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
#ifdef TRACER_PARTICLES
      particle_species(species_id)%tracer = &
          as_logical(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DTRACER_PARTICLES"
#endif
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for particle splitting
    ! *************************************************************
    IF (str_cmp(element, "split")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
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
      IF (deck_state .NE. c_ds_eio) RETURN
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
    IF (str_cmp(element, "ionise")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
#ifdef PARTICLE_IONISE
      particle_species(species_id)%ionise = &
          as_logical(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, "ionise_to_species")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
#ifdef PARTICLE_IONISE
      particle_species(species_id)%ionise_to_species = &
          as_integer(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, "release_species_on_ionise")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
#ifdef PARTICLE_IONISE
      particle_species(species_id)%release_species = &
          as_integer(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, "ionisation_energy")) THEN
      IF (deck_state .NE. c_ds_eio) RETURN
#ifdef PARTICLE_IONISE
      particle_species(species_id)%ionisation_energy = &
          as_real(value, handle_species_deck)
#else
      handle_species_deck = c_err_pp_options_wrong
      extended_error_string = "-DPARTICLE_IONISE"
#endif
      RETURN
    ENDIF

    ! Initial conditions

    IF (str_cmp(element, "offset")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      offset = as_long_integer_simple(value, handle_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "minrho")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      initial_conditions(species_id)%minrho = &
          as_real(value, handle_species_deck)
      RETURN
    ENDIF

    IF (str_cmp(element, "maxrho")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      initial_conditions(species_id)%maxrho = &
          as_real(value, handle_species_deck)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, handle_species_deck)

    IF (str_cmp(element, "rho") .OR. str_cmp(element, "number_density")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%rho(:), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%rho(:), &
            -2, nx+3, handle_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "mass_density")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%rho(:), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%rho(:), &
            -2, nx+3, handle_species_deck)
      ENDIF
      initial_conditions(species_id)%rho = &
          initial_conditions(species_id)%rho / particle_species(species_id)%mass
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_x")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%drift(:,1), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(:,1), &
            -2, nx+3, handle_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_y")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%drift(:,2), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(:,2), &
            -2, nx+3, handle_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_z")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%drift(:,3), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(:,3), &
            -2, nx+3, handle_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,1), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,1), &
            -2, nx+3, handle_species_deck)
      ENDIF
      debug_mode = .FALSE.
      initial_conditions(species_id)%temp(-2:nx+3, 2) = 0.0_num
      initial_conditions(species_id)%temp(-2:nx+3, 3) = 0.0_num
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_x")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,1), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,1), &
            -2, nx+3, handle_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_y")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,2), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,2), &
            -2, nx+3, handle_species_deck)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_z")) THEN
      IF (deck_state .NE. c_ds_ic) RETURN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,3), offset, &
            handle_species_deck)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,3), &
            -2, nx+3, handle_species_deck)
      ENDIF
      RETURN
    ENDIF

    handle_species_deck = c_err_unknown_element

  END FUNCTION handle_species_deck



  FUNCTION check_species_block()

    INTEGER :: check_species_block
    INTEGER :: i

    check_species_block = check_block

    IF (deck_state .NE. c_ds_ic) RETURN

    DO i = 1, n_species
      IF (particle_species(i)%mass .LT. 0) THEN
        IF (rank .EQ. 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'No mass specified for particle species "', &
              TRIM(particle_species(i)%name),'"'
          WRITE(40,*) '*** ERROR ***'
          WRITE(40,*) 'No mass specified for particle species "', &
              TRIM(particle_species(i)%name),'"'
        ENDIF
        check_species_block = c_err_missing_elements
      ENDIF
      IF (.NOT. species_charge_set(i)) THEN
        IF (rank .EQ. 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'No charge specified for particle species "', &
              TRIM(particle_species(i)%name),'"'
          WRITE(40,*) '*** ERROR ***'
          WRITE(40,*) 'No charge specified for particle species "', &
              TRIM(particle_species(i)%name),'"'
        ENDIF
        check_species_block = c_err_missing_elements
      ENDIF
    ENDDO

  END FUNCTION check_species_block

END MODULE deck_species_block
