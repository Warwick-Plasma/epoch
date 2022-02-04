! Copyright (C) 2009-2019 University of Warwick
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
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: species_names
  INTEGER, DIMENSION(:), POINTER :: species_blocks
  LOGICAL :: got_name
  INTEGER :: check_block = c_err_none
  LOGICAL, DIMENSION(:), ALLOCATABLE :: species_charge_set
  INTEGER, DIMENSION(:), ALLOCATABLE :: species_n, species_l
  LOGICAL :: use_ionise
  INTEGER :: n_secondary_species_in_block, n_secondary_limit
  LOGICAL :: unique_electrons
  CHARACTER(LEN=string_length) :: release_species_list
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: release_species
  REAL(num), DIMENSION(:), POINTER :: species_ionisation_energies
  REAL(num), DIMENSION(:), POINTER :: ionisation_energies
  REAL(num), DIMENSION(:), POINTER :: mass, charge
  LOGICAL, DIMENSION(:), POINTER :: auto_electrons
  INTEGER, DIMENSION(:), POINTER :: principle, angular, part_count
  INTEGER, DIMENSION(:), POINTER :: ionise_to_species, dumpmask_array
  INTEGER, DIMENSION(:,:), POINTER :: bc_particle_array
  REAL(num) :: species_mass, species_charge
  INTEGER :: species_dumpmask
  INTEGER :: species_atomic_number
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
      ALLOCATE(auto_electrons(4))
      release_species = ''
      release_species_list = ''
    END IF

  END SUBROUTINE species_deck_initialise



  SUBROUTINE species_deck_finalise

    INTEGER :: i, j, idx, io, iu, nlevels, nrelease, n_species_chain
    CHARACTER(LEN=8) :: string
    INTEGER :: errcode, bc
    TYPE(primitive_stack) :: stack
    INTEGER, DIMENSION(2*c_ndims) :: bc_species
    LOGICAL, ALLOCATABLE :: release_species_set(:)
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

      ALLOCATE(release_species_set(n_species))
      release_species_set = .FALSE.

      ! Scan for ionising species with automatically generated electron
      ! populations
      DO i = 1, n_species
        IF (auto_electrons(i)) THEN
          ! Deduce number of ionisation states attributed to the base state
          j = i
          DO WHILE(species_list(j)%ionise_to_species > 0)
            j = j + 1
          END DO
          ! Number of ionisation states, including the base state itself
          n_species_chain = j - i + 1

          ! Set release species for all ions. If auto-generation is used for a
          ! list of N species, then (N+1) to (2N-1) are the species ID for the
          ! release electrons (final ion in chain has no release)
          DO j = i, i + n_species_chain-2
            species_list(j)%release_species = j + n_species_chain
            release_species_set(j) = .TRUE.
          END DO
        END IF
      END DO

      DEALLOCATE(bc_particle_array)
      DEALLOCATE(dumpmask_array)
      DEALLOCATE(part_count)
      DEALLOCATE(principle)
      DEALLOCATE(angular)
      DEALLOCATE(charge)
      DEALLOCATE(mass)
      DEALLOCATE(ionisation_energies)

      ! Set release species of species_list elements which have been
      ! user-defined
      DO i = 1, n_species
        ! Release species is already present
        IF (release_species_set(i)) CYCLE

        ! No release species needed for a non-ionising species
        IF (.NOT. species_list(i)%ionise) CYCLE

        ! Error if no release species has been provided
        IF (TRIM(release_species(i)) == '') THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) ''
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Missing release species for ', TRIM(species_names(i))
              WRITE(io,*) ''
            END DO
            CALL abort_code(c_err_missing_elements)
          END IF
        END IF

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
            release_species_set(j) = .TRUE.
            j = species_list(j)%ionise_to_species
          END DO
        ! If there's a list of release species use it
        ELSE IF (nlevels == nrelease) THEN
          nlevels = 1
          j = i
          DO WHILE(species_list(j)%ionise)
            species_list(j)%release_species = stack%entries(nlevels)%value
            release_species_set(j) = .TRUE.
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
            release_species_set(j) = .TRUE.
            j = species_list(j)%ionise_to_species
          END DO
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) ''
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'Incorrect number of release species specified ', &
                  'for ', TRIM(species_names(i)), '. Using only first ', &
                  'specified.'
              WRITE(io,*) ''
            END DO
          END IF
        END IF

        CALL deallocate_stack(stack)
      END DO
      DEALLOCATE(auto_electrons)
      DEALLOCATE(release_species)
      DEALLOCATE(ionise_to_species)
      DEALLOCATE(species_names)
      DEALLOCATE(release_species_set)

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
              WRITE(io,*) ''
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
            WRITE(io,*) ''
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'The species named "' // TRIM(species_list(i)%name) &
                // '" must have a positive mass.'
            WRITE(io,*) ''
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

    use_ionise = .FALSE.
    unique_electrons = .FALSE.
    n_secondary_species_in_block = 0
    n_secondary_limit = 200  ! 200 allows all ionisations from any table element
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
    INTEGER :: max_ionisation, species_ionisation_state
    INTEGER :: i, io, iu, block_species_id
    INTEGER :: i_el, i_ion

    IF (.NOT.got_name) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(current_block, id_string)
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) ''
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Species block number ', TRIM(id_string), &
              ' has no "name" element.'
          WRITE(io,*) ''
        END DO
      END IF

      check_block = c_err_missing_elements
    END IF

    IF (deck_state == c_ds_first) THEN

      ! On first pass, read ionisation tables if using ionisation
      IF (use_ionise) THEN

        ! Ensure the user has entered an atomic number for this species
        IF (species_atomic_number < 1 .OR. species_atomic_number > 100) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) ''
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Ionising species must specify an atomic number'
              WRITE(io,*) ''
            END DO
          END IF
          check_block = c_err_missing_elements
          RETURN
        END IF

        ! Number of possible ionisation states
        species_ionisation_state = NINT(species_charge / q0)
        max_ionisation = species_atomic_number - species_ionisation_state

        ! User can ignore species above a certain ionisation-state
        n_secondary_species_in_block = MIN(max_ionisation, n_secondary_limit)

        ! Populate the species_ionisation_energies array
        IF (n_secondary_species_in_block > 0) THEN
          ALLOCATE(species_ionisation_energies(n_secondary_species_in_block))
          ALLOCATE(species_n(n_secondary_species_in_block))
          ALLOCATE(species_l(n_secondary_species_in_block))
          CALL read_ionisation_data(species_atomic_number, &
              species_ionisation_state, n_secondary_species_in_block, &
              species_ionisation_energies, species_l, species_n)
        END IF
      END IF

      block_species_id = n_species
      charge(n_species) = species_charge
      mass(n_species) = species_mass
      bc_particle_array(:, n_species) = species_bc_particle
      IF (n_secondary_species_in_block > 0) THEN
        ! Create an empty species for each ionisation level considered
        release_species(n_species) = release_species_list
        DO i = 1, n_secondary_species_in_block
          CALL integer_as_string(i, id_string)
          name = TRIM(TRIM(species_names(block_species_id))//id_string)
          CALL create_ionisation_species_from_name(name, &
              species_ionisation_energies(i), &
              n_secondary_species_in_block + 1 - i, species_n(i), species_l(i))
        END DO
        DEALLOCATE(species_ionisation_energies)

        ! Auto-generate unique electron release species if requested
        IF (unique_electrons) THEN
          auto_electrons(block_species_id) = .TRUE.
          DO i = 0, n_secondary_species_in_block-1
            ! Name of electron species, e.g., for a Carbon species, these are
            ! electron_from_Carbon
            ! electron_from_Carbon1
            IF (i == 0) THEN
              name = TRIM(TRIM('electron_from_' &
                  //species_names(block_species_id)))
            ELSE
              CALL integer_as_string(i, id_string)
              name = TRIM(TRIM('electron_from_' &
                  // species_names(block_species_id)) // id_string)
            END IF

            ! Create this species
            CALL create_electron_species_from_name(name, block_species_id, i)
          END DO
        END IF

        release_species(block_species_id) = release_species_list
      END IF

      IF (use_ionise) THEN
        DEALLOCATE(species_n, species_l)
      END IF

    ELSE
      ! On second pass, species have been defined - but auto-generated electrons
      ! do not appear in the input deck, so we must set their properties
      ! manually
      IF (unique_electrons) THEN
        i = species_id
        ! Loop over all ionising species in this chain
        DO WHILE(species_list(i)%ionise_to_species > 0)
          ! Electron charge was defined on creation
          i_el = species_list(i)%release_species
          species_charge_set(i_el) = .TRUE.

          ! Set properties of ionised species
          i_ion = species_list(i)%ionise_to_species
          species_list(i_ion)%charge = species_list(i)%charge - &
              species_list(i_el)%charge
          species_charge_set(i_ion) = .TRUE.
          species_list(i_ion)%mass = species_list(i)%mass - &
              species_list(i_el)%mass

          i = i_ion
        END DO
      END IF

    END IF

  END SUBROUTINE species_block_end



  FUNCTION species_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dmin, mult
    REAL(num), TARGET :: dummy(1,1)
    REAL(num), POINTER :: array(:,:)
    CHARACTER(LEN=string_length) :: filename, mult_string
    LOGICAL :: got_file, dump
    LOGICAL, SAVE :: warn_tracer = .TRUE.
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

    ! If set to T, then atomic number and charge state is used to deduce how
    ! many secondary particles there are
    IF (str_cmp(element, 'ionise') &
        .OR. str_cmp(element, 'ionize')) THEN
      use_ionise = as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! If using ionise, this can restrict the number of secondary particles to
    ! consider
    IF (str_cmp(element, 'ionise_limit') &
        .OR. str_cmp(element, 'ionize_limit')) THEN
      n_secondary_limit = as_integer_print(value, element, errcode)
      RETURN
    END IF

    ! If using ionise, this can restrict the number of secondary particles to
    ! consider
    IF (str_cmp(element, 'unique_electron_species')) THEN
      unique_electrons = as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! Support for manual writing of ionisation energies has been dropped. Issue
    ! warning
    IF (str_cmp(element, 'ionisation_energies') &
        .OR. str_cmp(element, 'ionization_energies')) THEN
      IF (deck_state == c_ds_first) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) ''
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Ionisation energies are now known up to Z=100'
            WRITE(io,*) 'EPOCH no longer supports manual entry of energies'
            WRITE(io,*) 'Ionisation of a species with atomic number Z is ', &
                'now activated by adding these'
            WRITE(io,*) 'lines to the species block:'
            WRITE(io,*) 'ionise = T'
            WRITE(io,*) 'atomic_no = Z # Replace Z with atomic number'
            WRITE(io,*) ''
          END DO
        END IF
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

    IF (str_cmp(element, 'atomic_no') &
        .OR. str_cmp(element, 'atomic_number')) THEN
      species_atomic_number = as_integer_print(value, element, errcode)
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

    IF (str_cmp(element, 'bc_y_min')) THEN
      species_bc_particle(c_bd_y_min) = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_y_max')) THEN
      species_bc_particle(c_bd_y_max) = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_z_min')) THEN
      RETURN
    END IF

    IF (str_cmp(element, 'bc_z_max')) THEN
      RETURN
    END IF

    IF (deck_state == c_ds_first) RETURN

    ! This sets up whether or not to use the MJ sampler for a species.
    ! It could go in the first deck pass, but that requires more temporary
    ! variables and seems unnecessary
    IF (str_cmp(element, 'use_maxwell_juettner') &
        .OR. str_cmp(element, 'use_maxwell_juttner')) THEN
      IF (as_logical_print(value, element, errcode)) THEN
        species_list(species_id)%ic_df_type = c_ic_df_relativistic_thermal
      ELSE
        species_list(species_id)%ic_df_type = c_ic_df_thermal
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'fractional_tail_cutoff')) THEN
      species_list(species_id)%fractional_tail_cutoff = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    ! *************************************************************
    ! This section identifies a species. Generic
    ! but currently only used in photon production
    ! *************************************************************
    IF (str_cmp(element, 'identify')) THEN
      CALL identify_species(value, errcode)

      ! If this particle is the release species of an ionising species, then
      ! subtract the charge and mass of this species from the ionising species
      ! to get the charge and mass of the child species.
      DO i = 1, n_species
        IF (species_id == species_list(i)%release_species) THEN
          j = species_list(i)%ionise_to_species
          DO WHILE(j > 0)
            species_list(j)%mass = species_list(j)%mass &
                - species_list(species_id)%mass
            species_list(j)%charge = species_list(j)%charge &
                - species_list(species_id)%charge
            species_charge_set(j) = .TRUE.
            j = species_list(j)%ionise_to_species
          END DO
        END IF
      END DO
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
            WRITE(io,*) ''
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
            WRITE(io,*) 'Particle species cannot have negative mass.'
            WRITE(io,*) ''
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

    IF (str_cmp(element, 'npart') &
        .OR. str_cmp(element, 'nparticles')) THEN
      species_list(species_id)%count = &
          as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'npart_per_cell') &
        .OR. str_cmp(element, 'nparticles_per_cell')) THEN
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

    IF (str_cmp(element, 'bc_y_min')) THEN
      species_list(species_id)%bc_particle(c_bd_y_min) = &
          species_bc_particle(c_bd_y_min)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_y_max')) THEN
      species_list(species_id)%bc_particle(c_bd_y_max) = &
          species_bc_particle(c_bd_y_max)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_z_min')) THEN
      RETURN
    END IF

    IF (str_cmp(element, 'bc_z_max')) THEN
      RETURN
    END IF

    IF (str_cmp(element, 'background_species') &
        .OR. str_cmp(element, 'background')) THEN
      species_list(species_id)%background_species = &
          as_logical_print(value, element, errcode)
      RETURN
    END IF

    ! *************************************************************
    ! This section sets properties for bremsstrahlung emission
    ! *************************************************************
    IF (str_cmp(element, 'atomic_no') &
        .OR. str_cmp(element, 'atomic_number')) THEN

      species_list(species_id)%atomic_no = species_atomic_number
      species_list(species_id)%atomic_no_set = .TRUE.

      ! Identify if the current species ionises to another species
      j = species_list(species_id)%ionise_to_species
      DO WHILE(j > 0)
        species_list(j)%atomic_no = species_list(species_id)%atomic_no
        species_list(j)%atomic_no_set = .TRUE.
        j = species_list(j)%ionise_to_species
      END DO
      RETURN
    END IF

    ! *************************************************************
    ! This section sets properties for zero_current particles
    ! *************************************************************
    IF (str_cmp(element, 'zero_current') .OR. str_cmp(element, 'tracer')) THEN
#ifndef NO_TRACER_PARTICLES
      species_list(species_id)%zero_current = &
          as_logical_print(value, element, errcode)
#else
      IF (as_logical_print(value, element, errcode)) THEN
        errcode = c_err_pp_options_wrong
        extended_error_string = '-DNO_TRACER_PARTICLES'
      END IF
#endif
      IF (warn_tracer .AND. rank == 0 .AND. str_cmp(element, 'tracer')) THEN
        warn_tracer = .FALSE.
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) ''
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'The "tracer" species do not behave in the way that ', &
                      'many users expect them to'
          WRITE(io,*) 'and can lead to unexpected and undesirable results. ', &
                      'Please see the'
          WRITE(io,*) 'documentation for further details.'
          WRITE(io,*) 'For this reason, the "tracer" flag is being renamed ', &
                      'to "zero_current".'
          WRITE(io,*) 'As of version 5.0, the "tracer" flag will be removed ', &
                      'entirely.'
          WRITE(io,*)
        END DO
      END IF
      RETURN
    END IF

    ! *************************************************************
    ! This section sets properties for particle splitting
    ! *************************************************************
    IF (str_cmp(element, 'split')) THEN
      species_list(species_id)%split = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'npart_max') &
        .OR. str_cmp(element, 'nparticles_max')) THEN
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

    IF (str_cmp(element, 'promote_density') &
        .OR. str_cmp(element, 'promote_number_density')) THEN
      species_list(species_id)%migrate%promotion_density = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'demote_density') &
        .OR. str_cmp(element, 'demote_number_density')) THEN
      species_list(species_id)%migrate%demotion_density = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    ! Initial conditions

    IF (str_cmp(element, 'offset')) THEN
      offset = as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density_min') .OR. str_cmp(element, 'minrho') &
        .OR. str_cmp(element, 'number_density_min')) THEN
      dmin = as_real_print(value, element, errcode)
      IF (dmin <= 0.0_num) dmin = EPSILON(1.0_num)
      species_list(species_id)%initial_conditions%density_min = dmin
      RETURN
    END IF

    IF (str_cmp(element, 'density_max') .OR. str_cmp(element, 'maxrho') &
        .OR. str_cmp(element, 'number_density_max')) THEN
      species_list(species_id)%initial_conditions%density_max = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density_back') &
        .OR. str_cmp(element, 'number_density_back') &
        .OR. str_cmp(element, 'density_background') &
        .OR. str_cmp(element, 'number_density_background')) THEN
      species_list(species_id)%initial_conditions%density_back = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    CALL get_filename(value, filename, got_file, errcode)

    mult = 1.0_num

    IF (str_cmp(element, 'density') .OR. str_cmp(element, 'rho') &
        .OR. str_cmp(element, 'number_density') &
        .OR. str_cmp(element, 'mass_density')) THEN

      IF (str_cmp(element, 'mass_density')) THEN
        mult = 1.0_num / species_list(species_id)%mass
        WRITE(mult_string, '(''*'', e23.15e3)') mult
        i = INDEX(mult_string, 'E')
        mult_string(i:i) = 'e'
      END IF

      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%density)) THEN
          ALLOCATE(ic%density(1-ng:nx+ng,1-ng:ny+ng))
        END IF
        array => ic%density
      ELSE
        array => dummy
      END IF

      CALL fill_array(species_list(species_id)%density_function, &
          array, mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_x') .OR. str_cmp(element, 'drift_px')) THEN
      n = 1
      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%drift)) THEN
          ALLOCATE(ic%drift(1-ng:nx+ng,1-ng:ny+ng,3))
          ic%drift = 0.0_num
        END IF
        array => ic%drift(:,:,n)
      ELSE
        array => dummy
      END IF

      CALL fill_array(species_list(species_id)%drift_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_y') .OR. str_cmp(element, 'drift_py')) THEN
      n = 2
      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%drift)) THEN
          ALLOCATE(ic%drift(1-ng:nx+ng,1-ng:ny+ng,3))
          ic%drift = 0.0_num
        END IF
        array => ic%drift(:,:,n)
      ELSE
        array => dummy
      END IF

      CALL fill_array(species_list(species_id)%drift_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_z') .OR. str_cmp(element, 'drift_pz')) THEN
      n = 3
      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%drift)) THEN
          ALLOCATE(ic%drift(1-ng:nx+ng,1-ng:ny+ng,3))
          ic%drift = 0.0_num
        END IF
        array => ic%drift(:,:,n)
      ELSE
        array => dummy
      END IF

      CALL fill_array(species_list(species_id)%drift_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'dist_fn')) THEN
      species_list(species_id)%ic_df_type = c_ic_df_arbitrary
      CALL initialise_stack(species_list(species_id)%dist_fn)
      CALL tokenize(value, species_list(species_id)%dist_fn, errcode)
      species_list(species_id)%dist_fn%should_simplify = .FALSE.
      RETURN
    END IF

    IF (str_cmp(element, 'dist_fn_px_range')) THEN
      CALL deallocate_stack(species_list(species_id)%dist_fn_range(1))
      CALL initialise_stack(species_list(species_id)%dist_fn_range(1))
      CALL tokenize(value, species_list(species_id)%dist_fn_range(1), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'dist_fn_py_range')) THEN
      CALL deallocate_stack(species_list(species_id)%dist_fn_range(2))
      CALL initialise_stack(species_list(species_id)%dist_fn_range(2))
      CALL tokenize(value, species_list(species_id)%dist_fn_range(2), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'dist_fn_pz_range')) THEN
      CALL deallocate_stack(species_list(species_id)%dist_fn_range(3))
      CALL initialise_stack(species_list(species_id)%dist_fn_range(3))
      CALL tokenize(value, species_list(species_id)%dist_fn_range(3), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_x_back') &
        .OR. str_cmp(element, 'drift_px_back') &
        .OR. str_cmp(element, 'drift_x_background') &
        .OR. str_cmp(element, 'drift_px_background')) THEN
      species_list(species_id)%initial_conditions%drift_back(1) = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_y_back') &
        .OR. str_cmp(element, 'drift_py_back') &
        .OR. str_cmp(element, 'drift_y_background') &
        .OR. str_cmp(element, 'drift_py_background')) THEN
      species_list(species_id)%initial_conditions%drift_back(2) = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_z_back') &
        .OR. str_cmp(element, 'drift_pz_back') &
        .OR. str_cmp(element, 'drift_z_background') &
        .OR. str_cmp(element, 'drift_pz_background')) THEN
      species_list(species_id)%initial_conditions%drift_back(3) = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    mult_string = '* ev / kb'

    IF (str_cmp(element, 'temp') &
        .OR. str_cmp(element, 'temp_k') &
        .OR. str_cmp(element, 'temp_ev') &
        .OR. str_cmp(element, 'temperature') &
        .OR. str_cmp(element, 'temperature_k') &
        .OR. str_cmp(element, 'temperature_ev')) THEN
      IF (str_cmp(element, 'temperature_ev') &
          .OR. str_cmp(element, 'temp_ev')) mult = ev / kb

      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%temp)) THEN
          ALLOCATE(ic%temp(1-ng:nx+ng,1-ng:ny+ng,3))
          ic%temp = 0.0_num
        END IF
      ELSE
        array => dummy
      END IF

      n = 1
      IF (got_file) array => ic%temp(:,:,n)
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
      n = 2
      IF (got_file) array => ic%temp(:,:,n)
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
      n = 3
      IF (got_file) array => ic%temp(:,:,n)
      CALL fill_array(species_list(species_id)%temperature_function(n), &
          array, mult, mult_string, element, value, filename, got_file)

      debug_mode = .FALSE.
      RETURN
    END IF

    IF (str_cmp(element, 'temp_back') &
        .OR. str_cmp(element, 'temp_back_k') &
        .OR. str_cmp(element, 'temp_back_ev') &
        .OR. str_cmp(element, 'temperature_background') &
        .OR. str_cmp(element, 'temperature_background_k') &
        .OR. str_cmp(element, 'temperature_background_ev')) THEN
      IF (str_cmp(element, 'temperature_background_ev') &
          .OR. str_cmp(element, 'temp_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(1) = &
          as_real_print(value, element, errcode) * mult
      species_list(species_id)%initial_conditions%temp_back(2) = &
          as_real_print(value, element, errcode) * mult
      species_list(species_id)%initial_conditions%temp_back(3) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_x_back') &
        .OR. str_cmp(element, 'temp_x_back_k') &
        .OR. str_cmp(element, 'temp_x_back_ev') &
        .OR. str_cmp(element, 'temperature_x_background') &
        .OR. str_cmp(element, 'temperature_x_background_k') &
        .OR. str_cmp(element, 'temperature_x_background_ev')) THEN
      IF (str_cmp(element, 'temperature_x_background_ev') &
          .OR. str_cmp(element, 'temp_x_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(1) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_y_back') &
        .OR. str_cmp(element, 'temp_y_back_k') &
        .OR. str_cmp(element, 'temp_y_back_ev') &
        .OR. str_cmp(element, 'temperature_y_background') &
        .OR. str_cmp(element, 'temperature_y_background_k') &
        .OR. str_cmp(element, 'temperature_y_background_ev')) THEN
      IF (str_cmp(element, 'temperature_y_background_ev') &
          .OR. str_cmp(element, 'temp_y_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(2) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_z_back') &
        .OR. str_cmp(element, 'temp_z_back_k') &
        .OR. str_cmp(element, 'temp_z_back_ev') &
        .OR. str_cmp(element, 'temperature_z_background') &
        .OR. str_cmp(element, 'temperature_z_background_k') &
        .OR. str_cmp(element, 'temperature_z_background_ev')) THEN
      IF (str_cmp(element, 'temperature_z_background_ev') &
          .OR. str_cmp(element, 'temp_z_back_ev')) mult = ev / kb

      species_list(species_id)%initial_conditions%temp_back(3) = &
          as_real_print(value, element, errcode) * mult
      RETURN
    END IF

    IF (str_cmp(element, 'temp_x') &
        .OR. str_cmp(element, 'temp_x_k') &
        .OR. str_cmp(element, 'temp_x_ev') &
        .OR. str_cmp(element, 'temperature_x') &
        .OR. str_cmp(element, 'temperature_x_k') &
        .OR. str_cmp(element, 'temperature_x_ev')) THEN
      IF (str_cmp(element, 'temperature_x_ev') &
          .OR. str_cmp(element, 'temp_x_ev')) mult = ev / kb

      n = 1
      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%temp)) THEN
          ALLOCATE(ic%temp(1-ng:nx+ng,1-ng:ny+ng,3))
          ic%temp = 0.0_num
        END IF
        array => ic%temp(:,:,n)
      ELSE
        array => dummy
      END IF

      CALL fill_array(species_list(species_id)%temperature_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'temp_y') &
        .OR. str_cmp(element, 'temp_y_k') &
        .OR. str_cmp(element, 'temp_y_ev') &
        .OR. str_cmp(element, 'temperature_y') &
        .OR. str_cmp(element, 'temperature_y_k') &
        .OR. str_cmp(element, 'temperature_y_ev')) THEN
      IF (str_cmp(element, 'temperature_y_ev') &
          .OR. str_cmp(element, 'temp_y_ev')) mult = ev / kb

      n = 2
      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%temp)) THEN
          ALLOCATE(ic%temp(1-ng:nx+ng,1-ng:ny+ng,3))
          ic%temp = 0.0_num
        END IF
        array => ic%temp(:,:,n)
      ELSE
        array => dummy
      END IF

      CALL fill_array(species_list(species_id)%temperature_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
      RETURN
    END IF

    IF (str_cmp(element, 'temp_z') &
        .OR. str_cmp(element, 'temp_z_k') &
        .OR. str_cmp(element, 'temp_z_ev') &
        .OR. str_cmp(element, 'temperature_z') &
        .OR. str_cmp(element, 'temperature_z_k') &
        .OR. str_cmp(element, 'temperature_z_ev')) THEN
      IF (str_cmp(element, 'temperature_z_ev') &
          .OR. str_cmp(element, 'temp_z_ev')) mult = ev / kb

      n = 3
      ic => species_list(species_id)%initial_conditions
      IF (got_file) THEN
        IF (.NOT. ASSOCIATED(ic%temp)) THEN
          ALLOCATE(ic%temp(1-ng:nx+ng,1-ng:ny+ng,3))
          ic%temp = 0.0_num
        END IF
        array => ic%temp(:,:,n)
      ELSE
        array => dummy
      END IF

      CALL fill_array(species_list(species_id)%temperature_function(n), &
          array, mult, mult_string, element, value, filename, got_file)
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
            WRITE(io,*) ''
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No mass specified for particle species "', &
                TRIM(species_list(i)%name), '"'
            WRITE(io,*) ''
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
      IF (.NOT. species_charge_set(i)) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) ''
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No charge specified for particle species "', &
                TRIM(species_list(i)%name), '"'
            WRITE(io,*) ''
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
      IF (species_list(i)%npart_per_cell >= 0) THEN
        IF (species_list(i)%count >= 0 .AND. rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) ''
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Two forms of nparticles used for particle species "', &
                TRIM(species_list(i)%name), '"'
            WRITE(io,*) 'Just using "nparticles_per_cell".'
            WRITE(io,*) ''
          END DO
        END IF
        species_list(i)%count = INT(species_list(i)%npart_per_cell, i8)
      END IF
    END DO

    ! Atomic numbers are only mandatory if running with bremsstrahlung
    IF (.NOT.use_bremsstrahlung) RETURN

    ! Have all species been assigned an atomic number?
    DO i = 1, n_species
      IF (species_list(i)%atomic_no_set) CYCLE

      ! Does this species ionise to another species? If so, the charge cannot
      ! be the atomic number, and we cannot run the bremsstrahlung module
      IF (species_list(i)%ionise_to_species > 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) ''
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) TRIM(species_list(i)%name), ' missing atomic number'
            WRITE(io,*) ''
          END DO
        END IF
        errcode = c_err_missing_elements
      ELSE
        species_list(i)%atomic_no = NINT(species_list(i)%charge/q0)
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) ''
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'No atomic number has been specified for species: ', &
                TRIM(species_list(i)%name)
            WRITE(io,*) 'Atomic number has been set to species particle charge'
            WRITE(io,*) ''
          END DO
        END IF
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
          WRITE(io,*) ''
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The species name "' // TRIM(name) // '" is not valid.'
          WRITE(io,*) 'Please choose a different name and try again.'
          WRITE(io,*) ''
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
    CALL grow_array(auto_electrons, n_species)
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
    auto_electrons(n_species) = .FALSE.
    bc_particle_array(:,n_species) = species_bc_particle

    RETURN

  END FUNCTION create_species_number_from_name



  SUBROUTINE read_ionisation_data(atomic_no, ion_state, ionise_num, &
        ionise_energy, ion_l, ion_n)

    ! Populates the array ionise_energy with energies taken from the file
    ! "ionisation_energies.table". The table lists ionisation energies [eV],
    ! with each line referring to an element of the corresponding atomic number
    ! (line 1 for H, line 2 for He, etc). Energies are listed in ascending
    ! order, and the ionise_energy array is filled starting from the ionisation
    ! state of the parent species ("ion_state"), and holds the next "ionise_num"
    ! energies.
    !
    ! A set of n and l quantum numbers are also read for the release electron,
    ! based on the ground-state configuration of element ions (taken from NIST).
    ! The release electron is assumed to be the electron missing when comparing
    ! the ground state electron energy configurations of subsequent ions. In
    ! some cases, two electrons will change position between subsequent ion
    ! ground-states - one removed and one changing orbitals. Here, we still use
    ! (n,l) of the vanishing electron. The format of "ion_l.table" and
    ! "ion_n.table" matches "ionisation_energies.table"

    INTEGER, INTENT(IN) :: atomic_no, ion_state, ionise_num
    REAL(num), INTENT(OUT) :: ionise_energy(:)
    INTEGER, INTENT(OUT) :: ion_l(:), ion_n(:)
    REAL(num), ALLOCATABLE :: full_line_energy(:)
    INTEGER, ALLOCATABLE :: full_line_l(:), full_line_n(:)
    INTEGER :: i_file, io, iu
    LOGICAL :: exists

    IF (atomic_no < 1 .OR. atomic_no > 100) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) ''
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'Ionising species must have an atomic number between'
        WRITE(io,*) '1 and 100'
        WRITE(io,*) ''
      END DO
      CALL abort_code(c_err_bad_value)
    END IF

    ! Check if the tables can be seen, issue warning if not
    INQUIRE(FILE='src/physics_packages/TABLES/ionisation_energies.table', &
        EXIST=exists)
    IF (.NOT.exists) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) ''
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'Unable to find the file:'
        WRITE(io,*) 'src/physics_packages/TABLES/ionisation_energies.table'
        WRITE(io,*) ''
      END DO
      CALL abort_code(c_err_io_error)
    END IF

    INQUIRE(FILE='src/physics_packages/TABLES/ion_l.table', &
        EXIST=exists)
    IF (.NOT.exists) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) ''
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'Unable to find the file:'
        WRITE(io,*) 'src/physics_packages/TABLES/ion_l.table'
        WRITE(io,*) ''
      END DO
      CALL abort_code(c_err_io_error)
    END IF

    INQUIRE(FILE='src/physics_packages/TABLES/ion_n.table', &
        EXIST=exists)
    IF (.NOT.exists) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) ''
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'Unable to find the file:'
        WRITE(io,*) 'src/physics_packages/TABLES/ion_n.table'
        WRITE(io,*) ''
      END DO
      CALL abort_code(c_err_io_error)
    END IF

    OPEN(UNIT = lu, &
        FILE = 'src/physics_packages/TABLES/ionisation_energies.table', &
        STATUS = 'OLD')
    OPEN(UNIT = lu + 1, &
        FILE = 'src/physics_packages/TABLES/ion_l.table', &
        STATUS = 'OLD')
    OPEN(UNIT = lu + 2, &
        FILE = 'src/physics_packages/TABLES/ion_n.table', &
        STATUS = 'OLD')

    ! Keep reading each file until the correct line is reached
    DO i_file = 1, atomic_no-1
      READ(lu,*)
      READ(lu+1,*)
      READ(lu+2,*)
    END DO

    ! Read the full line matching the current atomic number
    ALLOCATE(full_line_energy(atomic_no))
    ALLOCATE(full_line_l(atomic_no))
    ALLOCATE(full_line_n(atomic_no))
    READ(lu,*) full_line_energy(1:atomic_no)
    READ(lu+1,*) full_line_l(1:atomic_no)
    READ(lu+2,*) full_line_n(1:atomic_no)
    CLOSE(lu)
    CLOSE(lu+1)
    CLOSE(lu+2)

    ! Only consider data for ions starting at the current ion_state, and up
    ! to ion_state + ionise_num. Note ion_state=0 corresponds to table index 1.
    ! Also convert to [J] for ionise_energy
    ionise_energy = full_line_energy(ion_state+1:ion_state+ionise_num)*q0
    ion_l = full_line_l(ion_state+1:ion_state+ionise_num)
    ion_n = full_line_n(ion_state+1:ion_state+ionise_num)
    DEALLOCATE(full_line_energy, full_line_l, full_line_n)

  END SUBROUTINE read_ionisation_data



  SUBROUTINE create_ionisation_species_from_name(name, ionisation_energy, &
      n_electrons, n_in, l_in)

    ! The subroutine saves the variables of the current species to the temporary
    ! species arrays. Note that "name" refers to the next ionised state, but
    ! "n_species" refers to the current state until the line
    ! n_species = n_species + 1
    !
    ! E.g. if we have species like Carbon, Carbon1, Carbon2, Carbon3 ...
    ! Then if "name" is Carbon2, "species_names(n_species)" would initially be
    ! Carbon1

    CHARACTER(*), INTENT(IN) :: name
    REAL(num), INTENT(IN) :: ionisation_energy
    INTEGER, INTENT(IN) :: n_electrons
    INTEGER, INTENT(IN) :: n_in, l_in
    INTEGER :: i

    DO i = 1, n_species
      IF (str_cmp(name, species_names(i))) RETURN
    END DO

    ! Use quantum numbers of the electron which vanishes between subsequent
    ! ion groundstate electron configurations (from NIST)
    principle(n_species) = n_in
    angular(n_species) = l_in

    ! Set ionisation energy of the current species
    ionisation_energies(n_species) = ionisation_energy

    ! Append a new species to the species list, for the current species to
    ! ionise to
    ionise_to_species(n_species) = n_species + 1
    n_species = n_species + 1

    ! Ensure the temporary arrays for species information are large enough to
    ! contain values for the new species using "grow array". Initialise values
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
    CALL grow_array(auto_electrons, n_species)
    auto_electrons(n_species) = .FALSE.
    CALL grow_array(bc_particle_array, 2*c_ndims, n_species)
    bc_particle_array(:,n_species) = species_bc_particle
    RETURN

  END SUBROUTINE create_ionisation_species_from_name



  SUBROUTINE create_electron_species_from_name(name, block_species_id, i_el)

    ! The subroutine creates a release electron species with the name "name".
    ! The electron species is also set as the release species of the relevant
    ! ion, where the species ID is calculated using block_species_id and i_el

    CHARACTER(*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: block_species_id, i_el
    INTEGER :: i

    ! Check the species doesn't already exist
    DO i = 1, n_species
      IF (str_cmp(name, species_names(i))) RETURN
    END DO

    ! Append to number of species
    n_species = n_species + 1

    ! Ensure the temporary arrays for species information are large enough to
    ! contain values for the new species using "grow array". Initialise values
    CALL grow_array(species_names, n_species)
    species_names(n_species) = TRIM(name)
    CALL grow_array(ionise_to_species, n_species)
    ionise_to_species(n_species) = -1
    CALL grow_array(release_species, n_species)
    release_species(n_species) = ''
    CALL grow_array(mass, n_species)
    mass(n_species) = m0
    CALL grow_array(charge, n_species)
    charge(n_species) = -q0
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
    CALL grow_array(auto_electrons, n_species)
    auto_electrons(n_species) = .FALSE.
    CALL grow_array(bc_particle_array, 2*c_ndims, n_species)
    bc_particle_array(:,n_species) = species_bc_particle
    RETURN

  END SUBROUTINE create_electron_species_from_name



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
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(INOUT) :: array
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
            WRITE(io,*) ''
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Cannot load from file whilst using a moving window.'
            WRITE(io,*) ''
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
      CALL tokenize(value, stack, errcode, species_id)
      IF (ABS(mult - 1.0_num) > c_tiny) &
          CALL tokenize(mult_string, stack, errcode)

      ! Sanity check
      tmp = evaluate(stack, errcode)
      IF (errcode /= c_err_none) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) ''
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Unable to parse input deck.'
            WRITE(io,*) ''
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
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
      RETURN
    END IF

    IF (str_cmp(value, 'proton')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0 * 1836.2_num
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_proton
      species_list(species_id)%atomic_no = 1
      species_list(species_id)%atomic_no_set = .TRUE.
      RETURN
    END IF

    IF (str_cmp(value, 'positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_positron
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
      RETURN
    END IF

#if !defined(PHOTONS) && !defined(BREMSSTRAHLUNG)
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
      IF (use_qed .OR. use_bremsstrahlung) errcode = c_err_generic_warning
#endif
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
      RETURN
    END IF

    ! Breit Wheeler process electron
    IF (str_cmp(value, 'bw_electron') &
        .OR. str_cmp(value, 'breit_wheeler_electron')) THEN
      species_list(species_id)%charge = -q0
      species_list(species_id)%mass = m0
      species_list(species_id)%species_type = c_species_id_electron
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%electron = .TRUE.
#ifdef PHOTONS
      breit_wheeler_electron_species = species_id
#else
      IF (use_qed .OR. use_bremsstrahlung) errcode = c_err_generic_warning
#endif
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
      RETURN
    END IF

    IF (str_cmp(value, 'bw_positron') &
        .OR. str_cmp(value, 'breit_wheeler_positron')) THEN
      species_list(species_id)%charge = q0
      species_list(species_id)%mass = m0
      species_charge_set(species_id) = .TRUE.
      species_list(species_id)%species_type = c_species_id_positron
#ifdef PHOTONS
      breit_wheeler_positron_species = species_id
#endif
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
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
      IF (use_qed .OR. use_bremsstrahlung) errcode = c_err_generic_warning
#endif
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
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
      IF (use_qed .OR. use_bremsstrahlung) errcode = c_err_generic_warning
#endif
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
      RETURN
    END IF

    ! Bremsstrahlung photon
    IF (str_cmp(value, 'brem_photon')) THEN
      species_list(species_id)%charge = 0.0_num
      species_list(species_id)%mass = 0.0_num
      species_list(species_id)%species_type = c_species_id_photon
      species_charge_set(species_id) = .TRUE.
#ifdef BREMSSTRAHLUNG
      IF (bremsstrahlung_photon_species == -1) &
          bremsstrahlung_photon_species = species_id
#else
      IF (use_bremsstrahlung) errcode = c_err_generic_warning
#endif
      species_list(species_id)%atomic_no = 0
      species_list(species_id)%atomic_no_set = .TRUE.
      RETURN
    END IF

    errcode = IOR(errcode, c_err_bad_value)

  END SUBROUTINE identify_species

END MODULE deck_species_block
