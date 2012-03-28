MODULE deck_species_block

  USE mpi
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
  INTEGER :: n_secondary_species_in_block
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
          PRINT*, 'Name of species ', TRIM(ADJUSTL(string)), ' is ', &
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

    n_secondary_species_in_block = 0
    current_block = current_block + 1
    got_name = .FALSE.
    IF (deck_state .EQ. c_ds_first) RETURN
    species_id = species_blocks(current_block)
    offset = 0

  END SUBROUTINE species_block_start



  SUBROUTINE species_block_end

    CHARACTER(LEN=8) :: id_string
    CHARACTER(LEN=string_length) :: name
    INTEGER :: io, i, itmp, block_species_id

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

    IF (deck_state .EQ. c_ds_first) THEN
      block_species_id = n_species
      DO i = 1, n_secondary_species_in_block
        CALL integer_as_string(i, id_string)
        name = TRIM(TRIM(species_names(block_species_id))//id_string)
        itmp = species_number_from_name(name)
      ENDDO
    ENDIF

  END SUBROUTINE species_block_end



  FUNCTION species_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
#ifdef PARTICLE_IONISE
    TYPE(primitive_stack) :: stack
    TYPE(particle_species), POINTER :: base_species, species
    REAL(num), DIMENSION(:), POINTER :: dat
    INTEGER :: i
#endif
    REAL(num), DIMENSION(:,:), POINTER :: array
    REAL(num) :: dmin, mult
    CHARACTER(LEN=string_length) :: filename, mult_string
    LOGICAL :: got_file, dump
    INTEGER :: io, n

    errcode = c_err_none
    IF (value .EQ. blank .OR. element .EQ. blank) RETURN

    IF (str_cmp(element, 'name')) THEN
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

    IF (str_cmp(element, 'ionisation_energies')) THEN
#ifdef PARTICLE_IONISE
      IF (deck_state .EQ. c_ds_first) THEN
        NULLIFY(dat)
        stack%stack_point = 0
        CALL tokenize(value, stack, errcode)
        CALL evaluate_and_return_all(stack, 0, 0, &
            n_secondary_species_in_block, dat, errcode)
        DEALLOCATE(dat)
        RETURN
      ENDIF
#endif
    ENDIF

    IF (deck_state .EQ. c_ds_first) RETURN

    ! *************************************************************
    ! This section identifies a species. Generic
    ! but currently only used in photon production
    ! *************************************************************
    IF (str_cmp(element,'identify')) THEN
      CALL identify_species(value,errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mass')) THEN
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

    IF (str_cmp(element, 'charge')) THEN
      species_list(species_id)%charge = as_real(value, errcode) * q0
      species_charge_set(species_id) = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'frac') .OR. str_cmp(element, 'fraction')) THEN
      IF (npart_global .GE. 0) THEN
        species_list(species_id)%count = &
            INT(as_real(value, errcode) * npart_global)
      ELSE
        species_list(species_id)%count = 0
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'npart')) THEN
      species_list(species_id)%count = as_long_integer(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'npart_per_cell')) THEN
      species_list(species_id)%npart_per_cell = as_real(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'dump')) THEN
      dump = as_logical(value, errcode)
      IF (dump) THEN
        species_list(species_id)%dumpmask = c_io_always
      ELSE
        species_list(species_id)%dumpmask = c_io_never
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'dumpmask')) THEN
      species_list(species_id)%dumpmask = as_integer(value, errcode)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for tracer particles
    ! *************************************************************
    IF (str_cmp(element, 'tracer')) THEN
#ifdef TRACER_PARTICLES
      species_list(species_id)%tracer = as_logical(value, errcode)
#else
      IF (as_logical(value, errcode)) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-DTRACER_PARTICLES'
      ENDIF
#endif
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for particle splitting
    ! *************************************************************
    IF (str_cmp(element, 'split')) THEN
      species_list(species_id)%split = as_logical(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'npart_max')) THEN
      species_list(species_id)%npart_max = as_long_integer(value, errcode)
      RETURN
    ENDIF

    ! *************************************************************
    ! This section sets properties for ionisation
    ! *************************************************************
    IF (str_cmp(element, 'ionisation_electron_species') &
        .OR. str_cmp(element, 'electron')) THEN
#ifdef PARTICLE_IONISE
      species_list(species_id)%release_species = as_integer(value, errcode)
#endif
      RETURN
    ENDIF

    IF (str_cmp(element, 'ionisation_energies')) THEN
#ifdef PARTICLE_IONISE
      IF (species_list(species_id)%release_species .LT. 0) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Attempting to set ionisation energies without', &
                ' specifying electron species'
          ENDDO
        ENDIF
        errcode = c_err_required_element_not_set
        extended_error_string = 'ionisation_electron_species'
        RETURN
      ENDIF

      IF (species_list(species_id)%ionise) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Attempting to set ionisation energies twice for'
            WRITE(io,*) 'species "' &
                // TRIM(species_list(species_id)%name) // '"'
          ENDDO
        ENDIF
        errcode = c_err_preset_element
        RETURN
      ENDIF

      species_list(species_id)%ionise = .TRUE.
      NULLIFY(dat)
      stack%stack_point = 0
      CALL tokenize(value, stack, errcode)
      CALL evaluate_and_return_all(stack, 0, 0, &
          n_secondary_species_in_block, dat, errcode)

      base_species=>species_list(species_id)
      species=>species_list(species_id)

      DO i = 1, n_secondary_species_in_block
        CALL create_ion_subspecies(base_species, species, i, dat(i))
      ENDDO
      DEALLOCATE(dat)
#endif
      RETURN
    ENDIF

    IF (ic_from_restart) RETURN

    ! Initial conditions

    IF (str_cmp(element, 'offset')) THEN
      offset = as_long_integer_simple(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'density_min') .OR. str_cmp(element, 'minrho')) THEN
      dmin = as_real(value, errcode)
      IF (dmin .LE. 0.0_num) dmin = EPSILON(1.0_num)
      initial_conditions(species_id)%density_min = dmin
      RETURN
    ENDIF

    IF (str_cmp(element, 'density_max') .OR. str_cmp(element, 'maxrho')) THEN
      initial_conditions(species_id)%density_max = as_real(value, errcode)
      RETURN
    ENDIF

    CALL get_filename(value, filename, got_file, errcode)

    mult = 1.0_num

    IF (str_cmp(element, 'density') .OR. str_cmp(element, 'rho') &
        .OR. str_cmp(element, 'mass_density')) THEN

      IF (str_cmp(element, 'mass_density')) THEN
        mult = 1.0_num / species_list(species_id)%mass
        mult_string = '/ species_list(species_id)%mass'
      ENDIF

      CALL fill_array(species_list(species_id)%density_function, &
          initial_conditions(species_id)%density, &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_x')) THEN
      array => initial_conditions(species_id)%drift(:,:,1)
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, array, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, array, -2, nx+3, -2, ny+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_y')) THEN
      array => initial_conditions(species_id)%drift(:,:,2)
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, array, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, array, -2, nx+3, -2, ny+3, errcode)
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_z')) THEN
      array => initial_conditions(species_id)%drift(:,:,3)
      IF (got_file) THEN
        CALL load_single_array_from_file(filename, array, offset, errcode)
      ELSE
        CALL evaluate_string_in_space(value, array, -2, nx+3, -2, ny+3, errcode)
      ENDIF
      RETURN
    ENDIF

    mult_string = '* ev / kb'

    IF (str_cmp(element, 'temp') .OR. str_cmp(element, 'temp_k') &
        .OR. str_cmp(element, 'temp_ev')) THEN
      IF (str_cmp(element, 'temp_ev')) mult = ev / kb

      n = 1
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,n), &
          mult, mult_string, element, value, filename, got_file)

      species_list(species_id)%temperature_function(2) = &
          species_list(species_id)%temperature_function(n)
      CALL set_stack_zero(species_list(species_id)%temperature_function(3))

      debug_mode = .FALSE.
      initial_conditions(species_id)%temp(:,:,2) = &
          initial_conditions(species_id)%temp(:,:,n)
      initial_conditions(species_id)%temp(:,:,3) = 0.0_num
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_x') .OR. str_cmp(element, 'temp_x_k') &
        .OR. str_cmp(element, 'temp_x_ev')) THEN
      IF (str_cmp(element, 'temp_x_ev')) mult = ev / kb

      n = 1
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_y') .OR. str_cmp(element, 'temp_y_k') &
        .OR. str_cmp(element, 'temp_y_ev')) THEN
      IF (str_cmp(element, 'temp_y_ev')) mult = ev / kb

      n = 2
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,n), &
          mult, mult_string, element, value, filename, got_file)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_z') .OR. str_cmp(element, 'temp_z_k') &
        .OR. str_cmp(element, 'temp_z_ev')) THEN
      IF (str_cmp(element, 'temp_z_ev')) mult = ev / kb

      n = 3
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          initial_conditions(species_id)%temp(:,:,n), &
          mult, mult_string, element, value, filename, got_file)
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
                TRIM(species_list(i)%name), '"'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
      IF (.NOT. species_charge_set(i)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No charge specified for particle species "', &
                TRIM(species_list(i)%name), '"'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
      IF (species_list(i)%npart_per_cell .GE. 0) THEN
        IF (species_list(i)%count .GE. 0 .AND. rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Two forms of npart used for particle species "', &
                TRIM(species_list(i)%name), '"'
            WRITE(io,*) 'Just using "npart_per_cell".'
          ENDDO
        ENDIF
        species_list(i)%count = species_list(i)%npart_per_cell
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



  SUBROUTINE fill_array(output, array, mult, mult_string, element, value, &
      filename, got_file)

    TYPE(primitive_stack), INTENT(INOUT) :: output
    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: array
    REAL(num), INTENT(IN) :: mult
    CHARACTER(LEN=*), INTENT(IN) :: mult_string, element, value, filename
    LOGICAL, INTENT(IN) :: got_file
    TYPE(stack_element) :: block
    TYPE(primitive_stack) :: stack
    INTEGER :: io, ix, iy

    CALL initialise_stack(stack)
    IF (got_file) THEN
      IF (move_window) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Cannot load from file whilst using a moving window.'
          ENDDO
        ENDIF
        errcode = c_err_bad_value
        RETURN
      ENDIF

      CALL load_single_array_from_file(filename, array, offset, errcode)

      CALL load_block(species_list(species_id)%name, block)
      CALL push_to_stack(stack, block)
      CALL load_block(element, block)
      CALL push_to_stack(stack, block)
      IF (mult .NE. 1) array = mult * array
    ELSE
      CALL tokenize(value, stack, errcode)
      IF (mult .NE. 1) CALL tokenize(mult_string, stack, errcode)

      DO iy = -2, ny+3
        DO ix = -2, nx+3
          array(ix,iy) = evaluate_at_point(stack, ix, iy, errcode)
        ENDDO
      ENDDO
    ENDIF

    CALL deallocate_stack(output)
    output = stack

  END SUBROUTINE fill_array



  SUBROUTINE identify_species(value,errcode)
    CHARACTER(*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: errcode

    !Just a plain old electron
    IF (str_cmp(value,'electron')) THEN
      species_list(species_id)%charge=-q0
      species_list(species_id)%mass=m0
      species_charge_set(species_id)=.TRUE.
      species_list(species_id)%species_type=c_species_id_electron
      RETURN
    ENDIF
    !trident process electron
    IF (str_cmp(value,'trident_electron')) THEN
#ifdef PHOTONS
#ifdef TRIDENT_PHOTONS
      species_list(species_id)%charge=-q0
      species_list(species_id)%mass=m0
      species_list(species_id)%species_type=c_species_id_electron
      species_charge_set(species_id)=.TRUE.
      trident_electron_species = species_id
#else
      !built with photon support but not trident support
      errcode = c_err_generic_error
      extended_error_string = "Cannot identify species " // &
        TRIM(species_list(species_id)%name) // " as " //&
        TRIM(value) // " as compiler option -DTRIDENT_PHOTONS "//&
        " has not been set."
#endif
#else
      !built with neither photon support nor trident support
      errcode = c_err_generic_error
      extended_error_string = "Cannot identify species " // &
        TRIM(species_list(species_id)%name) // " as " //&
        TRIM(value) // " as compiler option -DPHOTONS " //&
        "-DTRIDENT_PHOTONS has not been set."
#endif
      RETURN
    ENDIF

    !Breit Wheeler process electron
    IF (str_cmp(value,'bw_electron')) THEN
#ifdef PHOTONS
      species_list(species_id)%charge=-q0
      species_list(species_id)%mass=m0
      species_list(species_id)%species_type=c_species_id_electron
      species_charge_set(species_id)=.TRUE.
      breit_wheeler_electron_species = species_id
#else
      !built without photon support
      errcode = c_err_generic_error
      extended_error_string = "Cannot identify species " // &
        TRIM(species_list(species_id)%name) // " as " //&
        TRIM(value) // " as compiler options -DPHOTONS has not been set."
#endif
      RETURN
    ENDIF

    IF (str_cmp(value,'proton')) THEN
      species_list(species_id)%charge=q0
      species_list(species_id)%mass=m0*1836.2_num
      species_charge_set(species_id)=.TRUE.
      species_list(species_id)%species_type=c_species_id_proton
      RETURN
    ENDIF

    IF (str_cmp(value,'positron')) THEN
      species_list(species_id)%charge=q0
      species_list(species_id)%mass=m0
      species_charge_set(species_id)=.TRUE.
      species_list(species_id)%species_type=c_species_id_positron
      RETURN
    ENDIF

    IF (str_cmp(value,'bw_positron')) THEN
      species_list(species_id)%charge=q0
      species_list(species_id)%mass=m0
      species_charge_set(species_id)=.TRUE.
      species_list(species_id)%species_type=c_species_id_positron
#ifdef PHOTONS
      breit_wheeler_positron_species=species_id
#endif
      RETURN
    ENDIF

    !trident process positron
    IF (str_cmp(value,'trident_positron')) THEN
#ifdef PHOTONS
#ifdef TRIDENT_PHOTONS
      species_list(species_id)%charge=q0
      species_list(species_id)%mass=m0
      species_list(species_id)%species_type=c_species_id_positron
      species_charge_set(species_id)=.TRUE.
      trident_electron_species = species_id
#else
      !built with photon support but not trident support
      errcode = c_err_generic_error
      extended_error_string = "Cannot identify species " // &
        TRIM(species_list(species_id)%name) // " as " //&
        TRIM(value) // " as compiler option -DTRIDENT_PHOTONS has not been set."
#endif
#else
      !built with neither photon support nor trident support
      errcode = c_err_generic_error
      extended_error_string = "Cannot identify species " // &
        TRIM(species_list(species_id)%name) // " as " //&
        TRIM(value) // " as compiler options -DPHOTONS " //&
        "-DTRIDENT_PHOTONS have not been set."
#endif
      RETURN
    ENDIF

    IF (str_cmp(value,'photon')) THEN
#ifdef PHOTONS
      species_list(species_id)%charge=0.0_num
      species_list(species_id)%mass=0.0_num
      species_list(species_id)%species_type=c_species_id_photon
      species_charge_set(species_id)=.TRUE.
      IF (photon_species .EQ. -1) photon_species=species_id
#else
      errcode = c_err_generic_error
      extended_error_string = "Cannot identify species " // &
        TRIM(species_list(species_id)%name) // " as " //&
        TRIM(value) // " as compiler option -DPHOTONS has not been set."
#endif
      RETURN
    ENDIF

    errcode = IAND(errcode,c_err_bad_value)

  END SUBROUTINE identify_species



#ifdef PARTICLE_IONISE
  SUBROUTINE create_ion_subspecies(base_species, current_species, index, energy)

    TYPE(particle_species), POINTER :: base_species, current_species
    INTEGER, INTENT(IN) :: index
    REAL(num), INTENT(IN) :: energy
    TYPE(particle_species), POINTER :: working_species
    INTEGER :: working_id, err
    CHARACTER(LEN=string_length) :: name, indexstring

    CALL integer_as_string(index, indexstring)
    name = TRIM(TRIM(base_species%name)//indexstring)

    working_id = as_integer(name, err)
    working_species => species_list(working_id)

    IF (rank .EQ. 0) THEN
      PRINT '(''Autocreated ionised form of '', a, '' as species '', i3, &
          & '' named '', a)', TRIM(base_species%name), working_id, TRIM(name)
    ENDIF

    working_species%charge = current_species%charge &
        - species_list(base_species%release_species)%charge
    working_species%mass = current_species%mass &
        - species_list(base_species%release_species)%mass
    working_species%dumpmask = base_species%dumpmask
    working_species%split = base_species%split
    working_species%npart_max = base_species%npart_max
    working_species%count = 0

    current_species%ionise = .TRUE.
    current_species%ionise_to_species = working_species%id
    current_species%release_species = base_species%release_species
    current_species%ionisation_energy = energy

    current_species=>working_species

  END SUBROUTINE create_ion_subspecies
#endif

END MODULE deck_species_block
