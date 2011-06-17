MODULE deck_species_block

  USE strings_advanced
  USE setup
  USE simple_io
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: species_deck_initialise, species_deck_finalise
  PUBLIC :: species_block_start, species_block_end
  PUBLIC :: species_block_handle_element, species_block_check

  INTEGER :: species_id, current_block
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: species_names
  INTEGER, DIMENSION(:), POINTER :: species_blocks
  LOGICAL, DIMENSION(:), ALLOCATABLE :: species_charge_set
  LOGICAL :: got_name
  INTEGER :: check_block = c_err_none

CONTAINS

  SUBROUTINE species_deck_initialise

    current_block = 0
    IF (deck_state .EQ. c_ds_first) THEN
      n_species = 0
      ALLOCATE(species_names(4))
      ALLOCATE(species_blocks(4))
    ENDIF

  END SUBROUTINE species_deck_initialise



  SUBROUTINE species_deck_finalise

    INTEGER :: i
    CHARACTER(LEN=8) :: string

    IF (deck_state .EQ. c_ds_first) THEN
      CALL setup_species
      ALLOCATE(species_charge_set(n_species))
      species_charge_set = .FALSE.

      DO i = 1, n_species
        species_list(i)%name = species_names(i)
        IF (rank .EQ. 0) THEN
          CALL integer_as_string(i, string)
          PRINT*,'Name of species ', TRIM(ADJUSTL(string)), ' is ', &
              TRIM(species_names(i))
        ENDIF
      ENDDO
      IF (rank .EQ. 0) PRINT*
      DEALLOCATE(species_names)
    ELSE
      DEALLOCATE(species_blocks)
      DEALLOCATE(species_charge_set)

      IF (dumpmask(c_dump_ejected_particles) .NE. c_io_never) THEN
        ALLOCATE(ejected_list(n_species))
        ejected_list = species_list
        DO i = 1, n_species
          ejected_list(i)%count = 0
          ejected_list(i)%name = &
              'ejected_' // TRIM(ADJUSTL(species_list(i)%name))
          CALL create_empty_partlist(ejected_list(i)%attached_list)
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE species_deck_finalise



  SUBROUTINE species_block_start

    current_block = current_block + 1
    got_name = .FALSE.
    IF (deck_state .EQ. c_ds_first) RETURN
    species_id = species_blocks(current_block)
    offset = 0

  END SUBROUTINE species_block_start



  SUBROUTINE species_block_end

    CHARACTER(LEN=8) :: id_string
    INTEGER :: io

    IF (.NOT.got_name) THEN
      IF (rank .EQ. 0) THEN
        CALL integer_as_string(current_block, id_string)
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Species block number ', TRIM(id_string), &
              ' has no "name" element.'
        ENDDO
      ENDIF

      check_block = c_err_missing_elements
    ENDIF

  END SUBROUTINE species_block_end



  FUNCTION species_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dmin
    CHARACTER(LEN=string_length) :: filename
    LOGICAL :: got_file, dump
    INTEGER :: io

    errcode = c_err_none
    IF (value .EQ. blank .OR. element .EQ. blank) RETURN

    IF (str_cmp(element, "name")) THEN
      IF (got_name) THEN
        errcode = c_err_preset_element
        RETURN
      ENDIF
      got_name = .TRUE.
      IF (deck_state .NE. c_ds_first) RETURN
      CALL grow_array(species_blocks, current_block)
      species_blocks(current_block) = species_number_from_name(value)
      RETURN
    ENDIF

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (str_cmp(element, "mass")) THEN
      species_list(species_id)%mass = as_real(value, errcode) * m0
      IF (species_list(species_id)%mass .LT. 0) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle species cannot have negative mass.'
          ENDDO
        ENDIF
        errcode = c_err_bad_value
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "charge")) THEN
      species_list(species_id)%charge = as_real(value, errcode) * q0
      species_charge_set(species_id) = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, "frac") .OR. str_cmp(element, "fraction")) THEN
      IF (npart_global .GE. 0) THEN
        species_list(species_id)%count = &
            INT(as_real(value, errcode) * npart_global)
      ELSE
        species_list(species_id)%count = 0
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "npart")) THEN
      species_list(species_id)%count = as_long_integer(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "dump")) THEN
      dump = as_logical(value, errcode)
      IF (dump) THEN
        species_list(species_id)%dumpmask = c_io_always
      ELSE
        species_list(species_id)%dumpmask = c_io_never
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "dumpmask")) THEN
      species_list(species_id)%dumpmask = as_integer(value, errcode)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for tracer particles
    ! *************************************************************
    IF (str_cmp(element, "tracer")) THEN
#ifdef TRACER_PARTICLES
      species_list(species_id)%tracer = as_logical(value, errcode)
#else
      IF (as_logical(value, errcode)) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = "-DTRACER_PARTICLES"
      ENDIF
#endif
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for particle splitting
    ! *************************************************************
    IF (str_cmp(element, "split")) THEN
      species_list(species_id)%split = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "npart_max")) THEN
      species_list(species_id)%npart_max = as_long_integer(value, errcode)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for ionisation
    ! *************************************************************
    IF (str_cmp(element, "ionise")) THEN
#ifdef PARTICLE_IONISE
      species_list(species_id)%ionise = as_logical(value, errcode)
#else
      IF (as_logical(value, errcode)) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = "-DPARTICLE_IONISE"
      ENDIF
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, "ionise_to_species")) THEN
#ifdef PARTICLE_IONISE
      species_list(species_id)%ionise_to_species = as_integer(value, errcode)
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, "release_species_on_ionise")) THEN
#ifdef PARTICLE_IONISE
      species_list(species_id)%release_species = as_integer(value, errcode)
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, "ionisation_energy")) THEN
#ifdef PARTICLE_IONISE
      species_list(species_id)%ionisation_energy = as_real(value, errcode)
#endif
      RETURN
    ENDIF

    IF (ic_from_restart) RETURN

    ! Initial conditions

    IF (str_cmp(element, "offset")) THEN
      offset = as_long_integer_simple(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, "density_min") .OR. str_cmp(element, "minrho")) THEN
      dmin = as_real(value, errcode)
      IF (dmin .LE. 0.0_num) dmin = EPSILON(1.0_num)
      initial_conditions(species_id)%density_min = dmin
      RETURN
    ENDIF

    IF (str_cmp(element, "density_max") .OR. str_cmp(element, "maxrho")) THEN
      initial_conditions(species_id)%density_max = as_real(value, errcode)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, errcode)

    IF (str_cmp(element, "density") .OR. str_cmp(element, "rho")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%density(:,:), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%density(:,:), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "mass_density")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%density(:,:), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%density(:,:), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      initial_conditions(species_id)%density = &
          initial_conditions(species_id)%density &
              / species_list(species_id)%mass
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_x")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%drift(:,:,1), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(:,:,1), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_y")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%drift(:,:,2), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(:,:,2), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "drift_z")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%drift(:,:,3), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%drift(:,:,3), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp") .OR. str_cmp(element, "temp_k") &
        .OR. str_cmp(element, "temp_ev")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,:,1), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,:,1), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      IF (str_cmp(element, "temp_ev")) THEN
        initial_conditions(species_id)%temp(:,:,1) = ev / kb * &
            initial_conditions(species_id)%temp(:,:,1)
      ENDIF
      debug_mode = .FALSE.
      initial_conditions(species_id)%temp(:,:,2) = &
          initial_conditions(species_id)%temp(:,:,1)
      initial_conditions(species_id)%temp(:,:,3) = 0.0_num
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_x") .OR. str_cmp(element, "temp_x_k") &
        .OR. str_cmp(element, "temp_x_ev")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,:,1), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,:,1), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      IF (str_cmp(element, "temp_x_ev")) THEN
        initial_conditions(species_id)%temp(:,:,1) = ev / kb * &
            initial_conditions(species_id)%temp(:,:,1)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_y") .OR. str_cmp(element, "temp_y_k") &
        .OR. str_cmp(element, "temp_y_ev")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,:,2), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,:,2), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      IF (str_cmp(element, "temp_y_ev")) THEN
        initial_conditions(species_id)%temp(:,:,2) = ev / kb * &
            initial_conditions(species_id)%temp(:,:,2)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, "temp_z") .OR. str_cmp(element, "temp_z_k") &
        .OR. str_cmp(element, "temp_z_ev")) THEN
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, &
            initial_conditions(species_id)%temp(:,:,3), offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, &
            initial_conditions(species_id)%temp(:,:,3), &
            -2, nx+3, -2, ny+3, errcode)
      ENDIF
      IF (str_cmp(element, "temp_z_ev")) THEN
        initial_conditions(species_id)%temp(:,:,3) = ev / kb * &
            initial_conditions(species_id)%temp(:,:,3)
      ENDIF
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION species_block_handle_element



  FUNCTION species_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: i, io

    errcode = check_block

    IF (deck_state .EQ. c_ds_first) RETURN

    DO i = 1, n_species
      IF (species_list(i)%mass .LT. 0) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No mass specified for particle species "', &
                TRIM(species_list(i)%name),'"'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
      IF (.NOT. species_charge_set(i)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No charge specified for particle species "', &
                TRIM(species_list(i)%name),'"'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
    ENDDO

  END FUNCTION species_block_check



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

END MODULE deck_species_block
