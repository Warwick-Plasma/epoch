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

MODULE deck_part_from_file_block

  USE strings_advanced
  USE utilities
  USE deck_species_block

  IMPLICIT NONE
  SAVE

  PRIVATE

  PUBLIC :: part_from_file_deck_initialise, part_from_file_deck_finalise
  PUBLIC :: part_from_file_block_start, part_from_file_block_end
  PUBLIC :: part_from_file_block_handle_element, part_from_file_block_check

  LOGICAL :: got_species
  INTEGER :: check_block
  INTEGER(KIND=8) :: current_offset
  INTEGER :: current_block_num
  INTEGER :: current_loader_id
  INTEGER, DIMENSION(:), POINTER :: loader_block_ids
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: loader_species

CONTAINS

  SUBROUTINE part_from_file_deck_initialise

    check_block = c_err_none
    current_loader_id = -1
    current_block_num = 0
    current_offset = 0
    IF (deck_state == c_ds_first) THEN
      n_custom_loaders = 0
      ALLOCATE(loader_block_ids(1))
      ALLOCATE(loader_species(1))
    END IF

  END SUBROUTINE part_from_file_deck_initialise



  SUBROUTINE part_from_file_deck_finalise

    ! First pass: Validate loader species and abort on error
    ! Allocate loader data structures for second pass
    IF (deck_state == c_ds_first) THEN
      IF (crosscheck_loader_species() /= c_err_none) THEN
        CALL abort_code(c_err_bad_value)
      END IF
      CALL setup_custom_loaders_list
      RETURN
    END IF
    DEALLOCATE(loader_block_ids)
    DEALLOCATE(loader_species)

  END SUBROUTINE part_from_file_deck_finalise



  SUBROUTINE part_from_file_block_start

#ifdef PER_PARTICLE_CHARGE_MASS
    INTEGER :: io, iu

    IF (rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) '"particles_from_file" block not supported when used ', &
            'with "PER_PARTICLE_CHARGE_MASS"'
      END DO
    END IF
    CALL abort_code(c_err_bad_setup)
#endif
    current_block_num = current_block_num + 1
    got_species = .FALSE.
    IF (deck_state == c_ds_last) THEN
      current_loader_id = loader_block_ids(current_block_num)
    END IF

  END SUBROUTINE part_from_file_block_start



  SUBROUTINE part_from_file_block_end

    CHARACTER(LEN=8) :: id_string
    INTEGER :: io, iu

    IF (.NOT.got_species) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(current_block_num, id_string)
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Particle block number ', TRIM(id_string), &
              ' has no "species" element.'
        END DO
      END IF
      check_block = c_err_missing_elements
    END IF

  END SUBROUTINE part_from_file_block_end



  FUNCTION part_from_file_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(LEN=string_length), INTENT(IN) :: element, value
    INTEGER :: errcode, filename_error_ignore
    CHARACTER(LEN=string_length) :: filename
    LOGICAL :: got_filename
    INTEGER :: io, iu

    errcode = c_err_none
    IF (value == blank .OR. element == blank) RETURN

    IF (str_cmp(element, 'species')) THEN
      IF (got_species) THEN
        errcode = c_err_preset_element
        RETURN
      END IF
      got_species = .TRUE.
      IF (deck_state /= c_ds_first) RETURN
      CALL grow_array(loader_block_ids, current_block_num)
      loader_block_ids(current_block_num) = index_by_species(value)
    END IF

    ! First pass just setup for each particles_from_file block
    IF (deck_state == c_ds_first) RETURN

    ! Second pass: Continue to read elements and write to loader structure

    IF (str_cmp(element, 'offset')) THEN
      current_offset = as_long_integer_print(value, element, errcode)
      RETURN
    END IF

    ! below here only expect to be dealing with filenames
    CALL get_filename(value, filename, got_filename, filename_error_ignore)

    IF (str_cmp(element, 'x_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%x_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%x_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'y_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%y_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%y_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'z_data')) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) '"z_data" was ignored'
          WRITE(io,*)
        END DO
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'px_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%px_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%px_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%px_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'py_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%py_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%py_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%py_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'pz_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%pz_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%pz_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%pz_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
    IF (str_cmp(element, 'w_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%w_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%w_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    IF (str_cmp(element, 'id4_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%id_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%id_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%id_data_4byte = .TRUE.
        custom_loaders_list(current_loader_id)%id_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF
    IF (str_cmp(element, 'id8_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%id_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%id_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%id_data_4byte = .FALSE.
        custom_loaders_list(current_loader_id)%id_data_offset = current_offset
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF
#endif

    ! If we got to here, then element is not specified
    errcode = c_err_unknown_element

  END FUNCTION part_from_file_block_handle_element



  FUNCTION part_from_file_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: i, io, iu
    CHARACTER(LEN=string_length) :: id_string
    TYPE(particle_species), POINTER :: species

    errcode = check_block

    IF (deck_state == c_ds_first) RETURN

    DO i = 1, n_custom_loaders
      species => species_list(custom_loaders_list(i)%species_id)

      IF (str_cmp(custom_loaders_list(i)%x_data, '')) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(current_block_num, id_string)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block number ', TRIM(id_string), &
                ' (species = "', TRIM(species%name), '") ', &
                'has no "x_data" element.'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF

      IF (str_cmp(custom_loaders_list(i)%y_data, '')) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(current_block_num, id_string)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block number ', TRIM(id_string), &
                ' (species = "', TRIM(species%name), '") ', &
                'has no "y_data" element.'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF

#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
      IF (str_cmp(custom_loaders_list(i)%w_data, '')) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(current_block_num, id_string)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block number ', TRIM(id_string), &
                ' (species = "', TRIM(species%name), '") ', &
                'has no "w_data" element.'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
#endif
    END DO

  END FUNCTION part_from_file_block_check



  FUNCTION index_by_species(species)

    CHARACTER(LEN=string_length), INTENT(IN) :: species
    INTEGER :: index_by_species
    INTEGER :: i

    DO i = 1, n_custom_loaders
      IF (str_cmp(species, loader_species(i))) THEN
        index_by_species = i
        RETURN
      END IF
    END DO

    n_custom_loaders = n_custom_loaders + 1
    index_by_species = n_custom_loaders

    CALL grow_array(loader_species, n_custom_loaders)
    loader_species(n_custom_loaders) = TRIM(species)

  END FUNCTION index_by_species



  FUNCTION crosscheck_loader_species() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: i, io, iu

    errcode = c_err_bad_value

    DO i = 1, n_custom_loaders
      IF (species_number_from_name(loader_species(i)) <= 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units  ! Print to all registered output buffers
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block with species "' &
                // TRIM(loader_species(i)) // '" ', &
                'must have a matching species defined'
          END DO
        END IF
        RETURN
      END IF
    END DO

    errcode = c_err_none

  END FUNCTION crosscheck_loader_species



  SUBROUTINE setup_custom_loaders_list

    INTEGER :: i

    IF (n_custom_loaders < 1) RETURN

    ALLOCATE(custom_loaders_list(n_custom_loaders))
    DO i = 1, n_custom_loaders
      custom_loaders_list(i)%species_id = &
          species_number_from_name(loader_species(i))
      custom_loaders_list(i)%x_data = ''
      custom_loaders_list(i)%x_data_offset = 0
      custom_loaders_list(i)%y_data = ''
      custom_loaders_list(i)%y_data_offset = 0
      custom_loaders_list(i)%px_data = ''
      custom_loaders_list(i)%px_data_offset = 0
      custom_loaders_list(i)%px_data_given = .FALSE.
      custom_loaders_list(i)%py_data = ''
      custom_loaders_list(i)%py_data_offset = 0
      custom_loaders_list(i)%py_data_given = .FALSE.
      custom_loaders_list(i)%pz_data = ''
      custom_loaders_list(i)%pz_data_offset = 0
      custom_loaders_list(i)%pz_data_given = .FALSE.
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
      custom_loaders_list(i)%w_data = ''
      custom_loaders_list(i)%w_data_offset = 0
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      custom_loaders_list(i)%id_data = ''
      custom_loaders_list(i)%id_data_offset = 0
      custom_loaders_list(i)%id_data_given = .FALSE.
#endif
    END DO

  END SUBROUTINE setup_custom_loaders_list

END MODULE deck_part_from_file_block
