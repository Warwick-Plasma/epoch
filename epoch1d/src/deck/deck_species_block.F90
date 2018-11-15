! Copyright (C) 2010-2016 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE deck_species_block

  USE strings_advanced
  USE setup
  USE simple_io
  USE utilities
  USE partlist

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: species_deck_initialise, species_deck_finalise
  PUBLIC :: species_block_start, species_block_end
  PUBLIC :: species_block_handle_element, species_block_check
  PUBLIC :: species_number_from_name

  INTEGER :: species_id, current_block
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
  CHARACTER(LEN=string_length), DIMENSION(:), ALLOCATABLE :: species_names
  INTEGER, DIMENSION(:), ALLOCATABLE :: species_blocks
  LOGICAL :: got_name
  INTEGER :: check_block = c_err_none
  LOGICAL, DIMENSION(:), ALLOCATABLE :: species_charge_set
  INTEGER :: n_secondary_species_in_block
  CHARACTER(LEN=string_length) :: release_species_list
  CHARACTER(LEN=string_length), DIMENSION(:), ALLOCATABLE :: release_species
  REAL(num), DIMENSION(:), ALLOCATABLE :: species_ionisation_energies
  REAL(num), DIMENSION(:), ALLOCATABLE :: ionisation_energies
  REAL(num), DIMENSION(:), ALLOCATABLE :: mass, charge
  INTEGER, DIMENSION(:), ALLOCATABLE :: principle, angular, part_count
  INTEGER, DIMENSION(:), ALLOCATABLE :: ionise_to_species, dumpmask_array
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: bc_particle_array
  REAL(num) :: species_mass, species_charge
  INTEGER :: species_dumpmask
  INTEGER, DIMENSION(2*c_ndims) :: species_bc_particle

CONTAINS

  SUBROUTINE species_deck_initialise

    current_block = 0
    IF (deck_state == c_ds_first) THEN
      n_species = 0
      ALLOCATE(species_names(4))
      ALLOCATE(species_blocks(4))
      ! All the following information is required during c_ds_first so that the
      ! derived ionisation species can be correctly set up
      ALLOCATE(ionise_to_species(4))
      ALLOCATE(release_species(4))
      ALLOCATE(ionisation_energies(4))
      ALLOCATE(mass(4))
      ALLOCATE(charge(4))
      ALLOCATE(principle(4))
      ALLOCATE(angular(4))
      ALLOCATE(part_count(4))
      ALLOCATE(dumpmask_array(4))
      ALLOCATE(bc_particle_array(2*c_ndims,4))
      release_species = ''
    END IF

  END SUBROUTINE species_deck_initialise



  SUBROUTINE species_deck_finalise

    INTEGER :: i, j, idx, io, iu, nlevels, nrelease
    CHARACTER(LEN=8) :: string
    INTEGER :: errcode, bc
    TYPE(primitive_stack) :: stack
    INTEGER, DIMENSION(2*c_ndims) :: bc_species
    LOGICAL :: error

    IF (deck_state == c_ds_first) THEN
      CALL setup_species
      ALLOCATE(species_charge_set(n_species))
      species_charge_set = .FALSE.

      DO i = 1, n_species
        species_list(i)%name = species_names(i)
        IF (rank == 0) THEN
          CALL integer_as_string(i, string)
          PRINT*, 'Name of species ', TRIM(ADJUSTL(string)), ' is ', &
              TRIM(species_names(i))
        END IF
        ! This would usually be set after c_ds_first but all of this is required
        ! during setup of derived ionisation species
        species_list(i)%ionise_to_species = ionise_to_species(i)
        species_list(i)%ionisation_energy = ionisation_energies(i)
        species_list(i)%n = principle(i)
        species_list(i)%l = angular(i)
        species_list(i)%mass = mass(i)
        species_list(i)%charge = charge(i)
        species_list(i)%count = INT(part_count(i),i8)
        species_list(i)%dumpmask = dumpmask_array(i)
        species_list(i)%bc_particle = bc_particle_array(:,i)
        IF (species_list(i)%ionise_to_species > 0) &
            species_list(i)%ionise = .TRUE.
      END DO

      DEALLOCATE(bc_particle_array)
      DEALLOCATE(dumpmask_array)
      DEALLOCATE(part_count)
      DEALLOCATE(principle)
      DEALLOCATE(angular)
      DEALLOCATE(charge)
      DEALLOCATE(mass)
      DEALLOCATE(ionisation_energies)

      DO i = 1, n_species
        IF (TRIM(release_species(i)) == '') CYCLE

        CALL initialise_stack(stack)
        CALL tokenize(release_species(i), stack, errcode)
        nlevels = 0
        j = i
        ! Count number of ionisation levels of species i
        DO WHILE(species_list(j)%ionise)
          nlevels = nlevels + 1
          j = species_list(j)%ionise_to_species
        END DO

        ! Count number of release species listed for species i; we need to do
        ! this because sometimes extra values are returned on the stack
        nrelease = 0
        DO j = 1, SIZE(stack%entries)
          IF (stack%entries(j)%value > 0 &
              .AND. stack%entries(j)%value <= n_species) &
                  nrelease = nrelease + 1
        END DO

        ! If there's only one release species use it for all ionisation levels
        IF (SIZE(stack%entries) == 1) THEN
          j = i
          species_list(stack%entries(1)%value)%electron = .TRUE.
          DO WHILE(species_list(j)%ionise)
            species_list(j)%release_species = stack%entries(1)%value
            j = species_list(j)%ionise_to_species
          END DO
        ! If there's a list of release species use it
        ELSE IF (nlevels == nrelease) THEN
          nlevels = 1
          j = i
          DO WHILE(species_list(j)%ionise)
            species_list(j)%release_species = stack%entries(nlevels)%value
            species_list(stack%entries(nlevels)%value)%electron = .TRUE.
            nlevels = nlevels + 1
            j = species_list(j)%ionise_to_species
          END DO
        ! If there's too many or not enough release species specified use the
        ! first one only and throw an error
        ELSE
          j = i
          species_list(stack%entries(1)%value)%electron = .TRUE.
          DO WHILE(species_list(j)%ionise)
            species_list(j)%release_species = stack%entries(1)%value
            j = species_list(j)%ionise_to_species
          END DO
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Incorrect number of release species specified ', &
                  'for ', TRIM(species_names(i)), '. Using only first ', &
                  'specified.'
            END DO
          END IF
        END IF

        CALL deallocate_stack(stack)
      END DO
      DEALLOCATE(release_species)
      DEALLOCATE(ionise_to_species)
      DEALLOCATE(species_names)

      ! Sanity check on periodic boundaries
      DO i = 1, n_species
        ! First, set the per-species boundary condition to the same value
        ! as bc_particle if it hasn't been set yet
        DO idx = 1, 2*c_ndims
          IF (species_list(i)%bc_particle(idx) == c_bc_null) &
              species_list(i)%bc_particle(idx) = bc_particle(idx)
        END DO

        bc_species = species_list(i)%bc_particle

        error = .FALSE.
        DO idx = 1, c_ndims
          IF (bc_species(2*idx-1) == c_bc_periodic) THEN
            IF (bc_species(2*idx) /= c_bc_periodic) &
                error = .TRUE.
          ELSE
            IF (bc_species(2*idx) == c_bc_periodic) &
                error = .TRUE.
          END IF
        END DO

        IF (error) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Periodic boundaries must be specified on both', &
                  ' sides of the domain.'
            END DO
          END IF
          CALL abort_code(c_err_bad_value)
        END IF

        ! Finally, set bc_allspecies. This will be equal to the per-species
        ! value if they are all compatible or c_bc_mixed otherwise
        DO idx = 1, 2*c_ndims
          bc = bc_species(idx)

          IF (bc == c_bc_reflect) THEN
            bc = c_bc_reflect
          ELSE IF (bc /= c_bc_periodic) THEN
            bc = c_bc_open
          END IF

          IF (i == 1) THEN
            bc_allspecies(idx) = bc
          ELSE
            IF (bc_allspecies(idx) /= bc) bc_allspecies(idx) = c_bc_mixed
          END IF
        END DO
      END DO
    ELSE
      DEALLOCATE(species_charge_set)
      DEALLOCATE(species_blocks)

      ! Sanity check
      DO i = 1, n_species
        IF (species_list(i)%species_type == c_species_id_photon) CYCLE
        IF (species_list(i)%mass > c_tiny) CYCLE
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'The species named "' // TRIM(species_list(i)%name) &
                // '" must have a positive mass.'
          END DO
        END IF
        CALL abort_code(c_err_bad_value)
      END DO

      IF (track_ejected_particles) THEN
        ALLOCATE(ejected_list(n_species))
        ejected_list = species_list
        DO i = 1, n_species
          ejected_list(i)%count = 0
          ejected_list(i)%name = &
              'ejected_' // TRIM(ADJUSTL(species_list(i)%name))
          CALL create_empty_partlist(ejected_list(i)%attached_list)
        END DO
      END IF
    END IF

    IF (use_field_ionisation) need_random_state = .TRUE.

  END SUBROUTINE species_deck_finalise



  SUBROUTINE species_block_start

    n_secondary_species_in_block = 0
    current_block = current_block + 1
    got_name = .FALSE.
    species_dumpmask = c_io_always
    species_bc_particle = c_bc_null
    IF (deck_state == c_ds_first) RETURN
    species_id = species_blocks(current_block)
    offset = 0

  END SUBROUTINE species_block_start



  SUBROUTINE species_block_end

    CHARACTER(LEN=8) :: id_string
    CHARACTER(LEN=string_length) :: name
    INTEGER :: i, io, iu, block_species_id

    IF (.NOT.got_name) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(current_block, id_string)
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Species block number ', TRIM(id_string), &
              ' has no "name" element.'
        END DO
      END IF

      check_block = c_err_missing_elements
    END IF

    IF (deck_state == c_ds_first) THEN
      block_species_id = n_species
      charge(n_species) = species_charge
      mass(n_species) = species_mass
      bc_particle_array(:, n_species) = species_bc_particle
      IF (n_secondary_species_in_block > 0) THEN
        ! Create an empty species for each ionisation energy listed in species
        ! block
        release_species(n_species) = release_species_list
        DO i = 1, n_secondary_species_in_block
          CALL integer_as_string(i, id_string)
          name = TRIM(TRIM(species_names(block_species_id))//id_string)
          CALL create_ionisation_species_from_name(name, &
              species_ionisation_energies(i), &
              n_secondary_species_in_block + 1 - i)
        END DO
        DEALLOCATE(species_ionisation_energies)
      END IF
    END IF

  END SUBROUTINE species_block_end



  FUNCTION species_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    TYPE(primitive_stack) :: stack
    REAL(num) :: dmin, mult
    CHARACTER(LEN=string_length) :: filename, mult_string
    LOGICAL :: got_file, dump
    INTEGER :: i, j, io, iu, n
    TYPE(initial_condition_block), POINTER :: ic

    errcode = c_err_none
    IF (value == blank .OR. element == blank) RETURN

    IF (str_cmp(element, 'name')) THEN
      IF (got_name) THEN
        errcode = c_err_preset_element
        RETURN
      END IF
      got_name = .TRUE.
      IF (deck_state /= c_ds_first) RETURN
      CALL grow_array(species_blocks, current_block)
      species_blocks(current_block) = create_species_number_from_name(value)
      RETURN
    END IF

    ! Collect ionisation energies for the species
    IF (str_cmp(element, 'ionisation_energies') &
        .OR. str_cmp(element, 'ionization_energies')) THEN
      IF (deck_state == c_ds_first) THEN
        CALL initialise_stack(stack)
        CALL tokenize(value, stack, errcode)
        CALL evaluate_and_return_all(stack, &
            n_secondary_species_in_block, species_ionisation_energies, errcode)
        CALL deallocate_stack(stack)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'ionisation_electron_species') &
        .OR. str_cmp(element, 'ionization_electron_species') &
        .OR. str_cmp(element, 'electron_species') &
        .OR. str_cmp(element, 'electron')) THEN
      release_species_list = value
      RETURN
    END IF

    IF (str_cmp(element, 'mass')) THEN
      species_mass = as_real_print(value, element, errcode) * m0
    END IF

    IF (str_cmp(element, 'charge')) THEN
      species_charge = as_real_print(value, element, errcode) * q0
    END IF

    IF (str_cmp(element, 'dump')) THEN
      dump = as_logical_print(value, element, errcode)
      IF (dump) THEN
        species_dumpmask = c_io_always
      ELSE
        species_dumpmask = c_io_never
      END IF
    END IF

    IF (str_cmp(element, 'dumpmask')) THEN
      species_dumpmask = as_integer_print(value, element, errcode)
    END IF

    IF (str_cmp(element, 'bc_x_min')) THEN
      species_bc_particle(c_bd_x_min) = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_x_max')) THEN
      species_bc_particle(c_bd_x_max) = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (deck_state == c_ds_first) RETURN

    ! *************************************************************
    ! This section identifies a species. Generic
    ! but currently only used in photon production
    ! *************************************************************
    IF (str_cmp(element, 'identify')) THEN
      CALL identify_species(value, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mass')) THEN
      species_list(species_id)%mass = species_mass
      ! Find the release species for each ionising species and subtract the
      ! release mass from the ionising species and each child species. Doing it
      ! like this ensures the right number of electron masses is removed for
      ! each ion.
      DO i = 1, n_species
        IF (species_id == species_list(i)%release_species) THEN
          j = species_list(i)%ionise_to_species
          DO WHILE(j > 0)
            species_list(j)%mass = species_list(j)%mass &
                - species_list(species_id)%mass
            j = species_list(j)%ionise_to_species
          END DO
        END IF
      END DO
      IF (species_list(species_id)%mass < 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle species cannot have negative mass.'
          END DO
        END IF
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'charge')) THEN
      species_list(species_id)%charge = species_charge
      species_charge_set(species_id) = .TRUE.
      ! Find the release species for each ionising species and subtract the
      ! release charge from the ionising species and each child species. Doing
      ! it like this ensures the right number of electron charges is removed for
      ! each ion. The species charge is considered set for the derived ionised
      ! species if it is touched in this routine.
      DO i = 1, n_species
        IF (species_id == species_list(i)%release_species) THEN
          j = species_list(i)%ionise_to_species
          DO WHILE(j > 0)
            species_list(j)%charge = species_list(j)%charge &
                - species_list(species_id)%charge
            species_charge_set(j) = .TRUE.
            j = species_list(j)%ionise_to_species
          END DO
        END IF
      END DO
      RETURN
    END IF

    IF (str_cmp(element, 'frac') .OR. str_cmp(element, 'fraction')) THEN
      IF (npart_global >= 0) THEN
        species_list(species_id)%count = &
            INT(as_real_print(value, element, errcode) * npart_global, i8)
      ELSE
        species_list(species_id)%count = 0
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'npart')) THEN
      species_list(species_id)%count = &
          as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'npart_per_cell')) THEN
      species_list(species_id)%npart_per_cell = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'dump') .OR. str_cmp(element, 'dumpmask')) THEN
      species_list(species_id)%dumpmask = species_dumpmask
      RETURN
    END IF

    IF (str_cmp(element, 'immobile')) THEN
      species_list(species_id)%immobile = &
          as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'meet_injectors') &
        .OR. str_cmp(element, 'load_up_to_injectors')) THEN
      species_list(species_id)%fill_ghosts = &
          as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_x_min')) THEN
      species_list(species_id)%bc_particle(c_bd_x_min) = &
          species_bc_particle(c_bd_x_min)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_x_max')) THEN
      species_list(species_id)%bc_particle(c_bd_x_max) = &
          species_bc_particle(c_bd_x_max)
      RETURN
    END IF

    ! *************************************************************
    ! This section sets properties for tracer particles
    ! *************************************************************
    IF (str_cmp(element, 'tracer')) THEN
#ifndef NO_TRACER_PARTICLES
      species_list(species_id)%tracer = &
          as_logical_print(value, element, errcode)
#else
      IF (as_logical_print(value, element, errcode)) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-DNO_TRACER_PARTICLES'
      END IF
#endif
      RETURN
    END IF

    ! *************************************************************
    ! This section sets properties for particle splitting
    ! *************************************************************
    IF (str_cmp(element, 'split')) THEN
      species_list(species_id)%split = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'npart_max')) THEN
      species_list(species_id)%npart_max = &
          as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    ! *************************************************************
    ! This section sets properties for migration
    ! *************************************************************

    IF (str_cmp(element, 'migrate')) THEN
      species_list(species_id)%migrate%this_species = &
          as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'promote_to')) THEN
      species_list(species_id)%migrate%promote_to_species = &
          as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'demote_to')) THEN
      species_list(species_id)%migrate%demote_to_species = &
          as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'promote_multiplier')) THEN
      species_list(species_id)%migrate%promotion_energy_factor = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'demote_multiplier')) THEN
      species_list(species_id)%migrate%demotion_energy_factor = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'promote_density')) THEN
      species_list(species_id)%migrate%promotion_density = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'demote_density')) THEN
      species_list(species_id)%migrate%demotion_density = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    ! Initial conditions

    IF (str_cmp(element, 'offset')) THEN
      offset = as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density_min') .OR. str_cmp(element, 'minrho')) THEN
      dmin = as_real_print(value, element, errcode)
      IF (dmin <= 0.0_num) dmin = EPSILON(1.0_num)
      species_list(species_id)%initial_conditions%density_min = dmin
      RETURN
    END IF

    IF (str_cmp(element, 'density_max') .OR. str_cmp(element, 'maxrho')) THEN
      species_list(species_id)%initial_conditions%density_max = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density_back')) THEN
      species_list(species_id)%initial_conditions%density_back = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    CALL get_filename(value, filename, got_file, errcode)

    mult = 1.0_num

    IF (str_cmp(element, 'density') .OR. str_cmp(element, 'rho') &
        .OR. str_cmp(element, 'mass_density')) THEN

      IF (str_cmp(element, 'mass_density')) THEN
        mult = 1.0_num / species_list(species_id)%mass
        mult_string = '/ species_list(species_id)%mass'
      END IF

      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%density)) THEN
        ALLOCATE(ic%density(1-ng:nx+ng))
      END IF

      CALL fill_array(species_list(species_id)%density_function, &
          ic%density, mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_x')) THEN
      n = 1
      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%drift)) THEN
        ALLOCATE(ic%drift(1-ng:nx+ng,3))
      END IF

      CALL fill_array(species_list(species_id)%drift_function(n), &
          ic%drift(:,n), mult, mult_string, element, value, filename, &
          got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_y')) THEN
      n = 2
      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%drift)) THEN
        ALLOCATE(ic%drift(1-ng:nx+ng,3))
      END IF

      CALL fill_array(species_list(species_id)%drift_function(n), &
          ic%drift(:,n), mult, mult_string, element, value, filename, &
          got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_z')) THEN
      n = 3
      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%drift)) THEN
        ALLOCATE(ic%drift(1-ng:nx+ng,3))
      END IF

      CALL fill_array(species_list(species_id)%drift_function(n), &
          ic%drift(:,n), mult, mult_string, element, value, filename, &
          got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_x_back')) THEN
      species_list(species_id)%initial_conditions%drift_back(1) = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_y_back')) THEN
      species_list(species_id)%initial_conditions%drift_back(2) = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_z_back')) THEN
      species_list(species_id)%initial_conditions%drift_back(3) = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    mult_string = '* ev / kb'

    IF (str_cmp(element, 'temp') .OR. str_cmp(element, 'temp_k') &
        .OR. str_cmp(element, 'temp_ev')) THEN
      IF (str_cmp(element, 'temp_ev')) mult = ev / kb

      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%temp)) THEN
        ALLOCATE(ic%temp(1-ng:nx+ng,3))
      END IF

      n = 1
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          ic%temp(:,n), mult, mult_string, element, value, filename, got_file)
      n = 2
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          ic%temp(:,n), mult, mult_string, element, value, filename, got_file)
      n = 3
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          ic%temp(:,n), mult, mult_string, element, value, filename, got_file)

      debug_mode = .FALSE.
      RETURN
    END IF

    IF (str_cmp(element, 'temp_back') .OR. str_cmp(element, 'temp_back_k') &
         .OR. str_cmp(element, 'temp_back_ev')) THEN
      IF (str_cmp(element, 'temp_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(1) = &
          as_real_print(value, element, errcode) * mult
      species_list(species_id)%initial_conditions%temp_back(2) = &
          as_real_print(value, element, errcode) * mult
      species_list(species_id)%initial_conditions%temp_back(3) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_x_back') &
         .OR. str_cmp(element, 'temp_x_back_ev')) THEN
      IF (str_cmp(element, 'temp_x_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(1) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_y_back') &
         .OR. str_cmp(element, 'temp_y_back_ev')) THEN
      IF (str_cmp(element, 'temp_y_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(2) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_z_back') &
         .OR. str_cmp(element, 'temp_z_back_ev')) THEN
      IF (str_cmp(element, 'temp_z_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(3) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_x') .OR. str_cmp(element, 'temp_x_k') &
        .OR. str_cmp(element, 'temp_x_ev')) THEN
      IF (str_cmp(element, 'temp_x_ev')) mult = ev / kb

      n = 1
      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%temp)) THEN
        ALLOCATE(ic%temp(1-ng:nx+ng,3))
      END IF

      CALL fill_array(species_list(species_id)%temperature_function(n), &
          ic%temp(:,n), mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'temp_y') .OR. str_cmp(element, 'temp_y_k') &
        .OR. str_cmp(element, 'temp_y_ev')) THEN
      IF (str_cmp(element, 'temp_y_ev')) mult = ev / kb

      n = 2
      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%temp)) THEN
        ALLOCATE(ic%temp(1-ng:nx+ng,3))
      END IF

      CALL fill_array(species_list(species_id)%temperature_function(n), &
          ic%temp(:,n), mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'temp_z') .OR. str_cmp(element, 'temp_z_k') &
        .OR. str_cmp(element, 'temp_z_ev')) THEN
      IF (str_cmp(element, 'temp_z_ev')) mult = ev / kb

      n = 3
      ic => species_list(species_id)%initial_conditions
      IF (got_file .AND. .NOT. ALLOCATED(ic%temp)) THEN
        ALLOCATE(ic%temp(1-ng:nx+ng,3))
      END IF

      CALL fill_array(species_list(species_id)%temperature_function(n), &
          ic%temp(:,n), mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION species_block_handle_element



  FUNCTION species_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: i, io, iu

    errcode = check_block

    IF (deck_state == c_ds_first) RETURN

    DO i = 1, n_species
      IF (species_list(i)%mass < 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No mass specified for particle species "', &
                TRIM(species_list(i)%name), '"'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
      IF (.NOT. species_charge_set(i)) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No charge specified for particle species "', &
                TRIM(species_list(i)%name), '"'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
      IF (species_list(i)%npart_per_cell >= 0) THEN
        IF (species_list(i)%count >= 0 .AND. rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Two forms of npart used for particle species "', &
                TRIM(species_list(i)%name), '"'
            WRITE(io,*) 'Just using "npart_per_cell".'
          END DO
        END IF
        species_list(i)%count = INT(species_list(i)%npart_per_cell, i8)
      END IF
    END DO

  END FUNCTION species_block_check



  FUNCTION create_species_number_from_name(name)

    CHARACTER(*), INTENT(IN) :: name
    INTEGER :: create_species_number_from_name
    INTEGER :: i, io, iu
    TYPE(stack_element) :: iblock

    DO i = 1, n_species
      IF (str_cmp(name, species_names(i))) THEN
        create_species_number_from_name = i
        RETURN
      END IF
    END DO

    ! If we're here then then named species doesn't yet exist

    ! First issue a warning message if the name overrides a built-in one
    CALL load_block(name, iblock)
    IF (iblock%ptype /= c_pt_bad .AND. iblock%ptype /= c_pt_null) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The species name "' // TRIM(name) // '" is not valid.'
          WRITE(io,*) 'Please choose a different name and try again.'
        END DO
      END IF
      CALL abort_code(c_err_bad_value)
    END IF

    n_species = n_species + 1
    create_species_number_from_name = n_species

    CALL grow_array(species_names, n_species)
    CALL grow_array(ionise_to_species, n_species)
    CALL grow_array(release_species, n_species)
    CALL grow_array(mass, n_species)
    CALL grow_array(charge, n_species)
    CALL grow_array(ionisation_energies, n_species)
    CALL grow_array(principle, n_species)
    CALL grow_array(angular, n_species)
    CALL grow_array(part_count, n_species)
    CALL grow_array(dumpmask_array, n_species)
    CALL grow_array(bc_particle_array, 2*c_ndims, n_species)

    species_names(n_species) = TRIM(name)
    ionise_to_species(n_species) = -1
    release_species(n_species) = ''
    mass(n_species) = -1.0_num
    charge(n_species) = 0.0_num
    ionisation_energies(n_species) = HUGE(0.0_num)
    principle(n_species) = -1
    angular(n_species) = -1
    part_count(n_species) = -1
    dumpmask_array(n_species) = species_dumpmask
    bc_particle_array(:,n_species) = species_bc_particle

    RETURN

  END FUNCTION create_species_number_from_name



  SUBROUTINE create_ionisation_species_from_name(name, ionisation_energy, &
      n_electrons)

    CHARACTER(*), INTENT(IN) :: name
    REAL(num), INTENT(IN) :: ionisation_energy
    INTEGER, INTENT(IN) :: n_electrons
    INTEGER :: i, n, l

    DO i = 1, n_species
      IF (str_cmp(name, species_names(i))) RETURN
    END DO
    ! This calculates the principle and angular quantum number based on the
    ! assumption that shells are filled as they would be in the ground state
    ! e.g. 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, 7s, etc
    n = 0
    l = 0
    i = 0
    DO WHILE(n_electrons > i)
      n = n + 1
      DO l = (n - 1) / 2, 0, -1
        i = i + 4 * l + 2
        IF (n_electrons <= i) THEN
          n = n - l
          EXIT
        END IF
      END DO
    END DO
    principle(n_species) = n
    angular(n_species) = l
    ionisation_energies(n_species) = ionisation_energy
    ionise_to_species(n_species) = n_species + 1
    n_species = n_species + 1
    CALL grow_array(species_names, n_species)
    species_names(n_species) = TRIM(name)
    CALL grow_array(ionise_to_species, n_species)
    ionise_to_species(n_species) = -1
    CALL grow_array(release_species, n_species)
    release_species(n_species) = ''
    CALL grow_array(mass, n_species)
    mass(n_species) = species_mass
    CALL grow_array(charge, n_species)
    charge(n_species) = species_charge
    CALL grow_array(ionisation_energies, n_species)
    ionisation_energies(n_species) = HUGE(0.0_num)
    CALL grow_array(principle, n_species)
    principle(n_species) = -1
    CALL grow_array(angular, n_species)
    angular(n_species) = -1
    CALL grow_array(part_count, n_species)
    part_count(n_species) = 0
    CALL grow_array(dumpmask_array, n_species)
    dumpmask_array(n_species) = species_dumpmask
    CALL grow_array(bc_particle_array, 2*c_ndims, n_species)
    bc_particle_array(:,n_species) = species_bc_particle
    RETURN

  END SUBROUTINE create_ionisation_species_from_name



  FUNCTION species_number_from_name(name)

    CHARACTER(*), INTENT(IN) :: name
    INTEGER :: species_number_from_name
    INTEGER :: i

    DO i = 1, n_species
      IF (str_cmp(name, species_list(i)%name)) THEN
        species_number_from_name = i
        RETURN
      END IF
    END DO
    species_number_from_name = -1
    RETURN

  END FUNCTION species_number_from_name



  SUBROUTINE fill_array(output, array, mult, mult_string, element, value, &
      filename, got_file)

    TYPE(primitive_stack), INTENT(INOUT) :: output
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: array
    REAL(num), INTENT(IN) :: mult
    CHARACTER(LEN=*), INTENT(IN) :: mult_string, element, value, filename
    LOGICAL, INTENT(IN) :: got_file
    TYPE(stack_element) :: iblock
    TYPE(primitive_stack) :: stack
    INTEGER :: io, iu
    REAL(num) :: tmp

    CALL initialise_stack(stack)
    IF (got_file) THEN
      IF (move_window) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Cannot load from file whilst using a moving window.'
          END DO
        END IF
        errcode = c_err_bad_value
        RETURN
      END IF

      CALL load_single_array_from_file(filename, array, offset, errcode)

      CALL load_block(species_list(species_id)%name, iblock)
      CALL push_to_stack(stack, iblock)
      CALL load_block(element, iblock)
      CALL push_to_stack(stack, iblock)
      IF (ABS(mult - 1.0_num) > c_tiny) array = mult * array
    ELSE
      CALL tokenize(value, stack, errcode)
      IF (ABS(mult - 1.0_num) > c_tiny) &
          CALL tokenize(mult_string, stack, errcode)

      ! Sanity check
      tmp = evaluate(stack, errcode)
      IF (errcode /= c_err_none) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Unable to parse input deck.'
          END DO
        END IF
        CALL abort_code(errcode)
      END IF
    END IF

    CALL deallocate_stack(output)
    output = stack

  END SUBROUTINE fill_array



  SUBROUTINE identify_species(value, errcode)

    CHARACTER(*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: errcode

    ! Just a plain old electron
    IF (str_cmp(value, 'electron')) THEN
      species_list(species_id)%charge = -q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_electron
      species_list(species_id)%electron = .TRUE.
      RETURN
    END IF

    IF (str_cmp(value, 'proton')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0 * 1836.2_num
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_proton
      RETURN
    END IF

    IF (str_cmp(value, 'positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_positron
      RETURN
    END IF

#ifndef PHOTONS
    extended_error_string = 'Cannot identify species "' &
        // TRIM(species_list(species_id)%name) // '" as "' // TRIM(value) &
        // '" because' // CHAR(10) &
        // ' compiler option -DPHOTONS has not been set.'
#else
#ifndef TRIDENT_PHOTONS
    extended_error_string = 'Cannot identify species "' &
        // TRIM(species_list(species_id)%name) // '" as "' // TRIM(value) &
        // '" because' // CHAR(10) &
        // ' compiler option -DTRIDENT_PHOTONS has not been set.'
#endif
#endif

    ! trident process electron
    IF (str_cmp(value, 'trident_electron')) THEN
      species_list(species_id)%charge = -q0
      species_list(species_id)%mass = m0
      species_list(species_id)%species_type = c_species_id_electron
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%electron = .TRUE.
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
      trident_electron_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    END IF

    ! Breit Wheeler process electron
    IF (str_cmp(value, 'bw_electron')) THEN
      species_list(species_id)%charge = -q0
      species_list(species_id)%mass = m0
      species_list(species_id)%species_type = c_species_id_electron
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%electron = .TRUE.
#ifdef PHOTONS
      breit_wheeler_electron_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    END IF

    IF (str_cmp(value, 'bw_positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_positron
#ifdef PHOTONS
      breit_wheeler_positron_species = species_id
#endif
      RETURN
    END IF

    ! trident process positron
    IF (str_cmp(value, 'trident_positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_list(species_id)%species_type = c_species_id_positron
      species_charge_set(species_id) = .TRUE.
#if defined(PHOTONS) && defined(TRIDENT_PHOTONS)
      trident_positron_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    END IF

    IF (str_cmp(value, 'photon')) THEN
      species_list(species_id)%charge = 0.0_num
      species_list(species_id)%mass = 0.0_num
      species_list(species_id)%species_type = c_species_id_photon
      species_charge_set(species_id) = .TRUE.
#ifdef PHOTONS
      IF (photon_species == -1) photon_species = species_id
#else
      errcode = c_err_generic_warning
#endif
      RETURN
    END IF

    errcode = IOR(errcode, c_err_bad_value)

  END SUBROUTINE identify_species

END MODULE deck_species_block
