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

MODULE deck_particles_block

  USE strings_advanced
  USE setup
  USE simple_io
  USE utilities
  USE partlist
  USE deck_species_block

  IMPLICIT NONE
  SAVE

  PRIVATE

  PUBLIC :: particles_deck_initialise, particles_deck_finalise
  PUBLIC :: particles_block_start, particles_block_end
  PUBLIC :: particles_block_handle_element, particles_block_check

  LOGICAL :: got_target
  INTEGER :: check_block
  INTEGER(KIND=8) :: current_offset
  INTEGER :: current_block_num
  INTEGER :: current_loader_id 
  INTEGER, DIMENSION(:), POINTER :: loader_block_ids
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: loader_targets
  
CONTAINS

  SUBROUTINE particles_deck_initialise

    check_block = c_err_none
    current_loader_id = -1
    current_block_num = 0
    current_offset = 0
    IF (deck_state == c_ds_first) THEN
      n_custom_loaders = 0
      ALLOCATE(loader_block_ids(1))
      ALLOCATE(loader_targets(1))
    ENDIF

  END SUBROUTINE particles_deck_initialise


  SUBROUTINE particles_deck_finalise
  
    ! First pass: Validate loader targets and abort on error
    ! Allocate loader data structures for second pass
    IF (deck_state == c_ds_first) THEN
      IF (crosscheck_loader_species() /= c_err_none) THEN
        CALL abort_code(c_err_bad_value) 
      ENDIF
      CALL setup_custom_loaders_list
    ENDIF
    IF (deck_state == c_ds_last) THEN
      CALL loader_debug_print
    ENDIF

  END SUBROUTINE particles_deck_finalise


  SUBROUTINE particles_block_start

#ifdef PER_PARTICLE_CHARGE_MASS
    INTEGER :: io iu

    IF (rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) '"particles" block not supported when used with "PER_PARTICLE_CHARGE_MASS'
      ENDDO
    ENDIF
    CALL abort_code(c_err_bad_setup)
#endif
    current_block_num = current_block_num + 1
    got_target = .FALSE.
    IF (deck_state == c_ds_last) THEN
      current_loader_id = loader_block_ids(current_block_num)
    ENDIF

  END SUBROUTINE particles_block_start


  SUBROUTINE particles_block_end

    CHARACTER(LEN=8) :: id_string
    INTEGER :: io, iu

     IF (.NOT.got_target) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(current_block_num, id_string)
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Particle block number ', TRIM(id_string), &
              ' has no "target" element.'
        ENDDO
      ENDIF
      check_block = c_err_missing_elements
    ENDIF

  END SUBROUTINE particles_block_end


  FUNCTION particles_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(LEN=string_length), INTENT(IN) :: element, value
    INTEGER :: errcode, filename_error_ignore
    CHARACTER(LEN=string_length) :: filename
    LOGICAL :: got_filename

    errcode = c_err_none
    IF (value == blank .OR. element == blank) RETURN

    IF (str_cmp(element, 'target')) THEN
      IF (got_target) THEN
        errcode = c_err_preset_element
        RETURN
      ENDIF
      got_target = .TRUE.
      IF (deck_state /= c_ds_first) RETURN
      CALL grow_array(loader_block_ids, current_block_num)
      loader_block_ids(current_block_num) = index_by_target_as_needed(value)
    ENDIF

    ! First pass just setup for each particles block
    IF (deck_state == c_ds_first) RETURN 

    ! Second pass: Continue to read elements and write to loader structure 
    
    IF (str_cmp(element, 'offset')) THEN
        current_offset = as_long_integer_print(value, element, errcode)
      RETURN
    ENDIF
      
    ! below here only expect to be dealing with filenames
    CALL get_filename(value, filename, got_filename, filename_error_ignore)
    PRINT*, '"', filename, '"'

    IF (str_cmp(element, 'x_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%x_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%x_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
 
    IF (str_cmp(element, 'y_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%y_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%y_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
  
    IF (str_cmp(element, 'z_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%z_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%z_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
 
    IF (str_cmp(element, 'px_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%px_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%px_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%px_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
 
    IF (str_cmp(element, 'py_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%py_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%py_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%py_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
 
    IF (str_cmp(element, 'pz_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%pz_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%pz_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%pz_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
 
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
    IF (str_cmp(element, 'w_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%w_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%w_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    IF (str_cmp(element, 'id_data')) THEN
      IF (got_filename) THEN
        custom_loaders_list(current_loader_id)%id_data = TRIM(filename)
        custom_loaders_list(current_loader_id)%id_data_given = .TRUE.
        custom_loaders_list(current_loader_id)%id_data_offset = current_offset
        RETURN
      ENDIF
      errcode = c_err_bad_value
      RETURN
    ENDIF
#endif

    ! If we got to here, then element is not specified
    errcode = c_err_unknown_element

  END FUNCTION particles_block_handle_element


  FUNCTION particles_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: i, io, iu
    CHARACTER(LEN=string_length) :: id_string
  
    errcode = check_block

    IF (deck_state == c_ds_first) RETURN

    DO i = 1, n_custom_loaders
      IF (str_cmp(custom_loaders_list(i)%x_data, '')) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(current_block_num, id_string)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block number ', TRIM(id_string), &
                ' (target = "', TRIM(custom_loaders_list(i)%target_name), '")', &
                ' has no "x_data" element.'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF      
      IF (str_cmp(custom_loaders_list(i)%y_data, '')) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(current_block_num, id_string)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block number ', TRIM(id_string), &
                ' (target = "', TRIM(custom_loaders_list(i)%target_name), '")', &
                ' has no "y_data" element.'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
      IF (str_cmp(custom_loaders_list(i)%z_data, '')) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(current_block_num, id_string)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block number ', TRIM(id_string), &
                ' (target = "', TRIM(custom_loaders_list(i)%target_name), '")', &
                ' has no "z_data" element.'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF      
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS) 
      IF (str_cmp(custom_loaders_list(i)%w_data, '')) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(current_block_num, id_string)
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block number ', TRIM(id_string), &
                ' (target = "', TRIM(custom_loaders_list(i)%target_name), '")', &
                ' has no "w_data" element.'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
#endif
    ENDDO

  END FUNCTION particles_block_check


  FUNCTION index_by_target_as_needed(target)

    CHARACTER(LEN=string_length), INTENT(IN) :: target
    INTEGER :: index_by_target_as_needed
    INTEGER :: i

    DO i = 1, n_custom_loaders
      IF (str_cmp(target, loader_targets(i))) THEN
        index_by_target_as_needed = i
        RETURN
      ENDIF
    ENDDO
  
    n_custom_loaders = n_custom_loaders + 1
    index_by_target_as_needed = n_custom_loaders

    CALL grow_array(loader_targets, n_custom_loaders)
    loader_targets(n_custom_loaders) = TRIM(target)
    
  END FUNCTION index_by_target_as_needed

  FUNCTION crosscheck_loader_species() RESULT(errcode)
    
    INTEGER :: errcode
    INTEGER :: i, io, iu

    errcode = c_err_bad_value

    DO i = 1, n_custom_loaders
      IF (species_number_from_name(loader_targets(i)) <= 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units  ! Print to all registered output buffers
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle block with target "' // TRIM(loader_targets(i)) &
              // '" must have a matching species defined/'
          ENDDO
        ENDIF
        RETURN
      ENDIF
    ENDDO

    errcode = c_err_none

  END FUNCTION crosscheck_loader_species

  SUBROUTINE setup_custom_loaders_list

    INTEGER :: i
    
    ALLOCATE(custom_loaders_list(n_custom_loaders))
    DO i = 1, n_custom_loaders
      custom_loaders_list(i)%target_name = loader_targets(i)
      custom_loaders_list(i)%target_id = species_number_from_name(loader_targets(i)) 
      custom_loaders_list(i)%x_data = ''
      custom_loaders_list(i)%x_data_offset = 0
      custom_loaders_list(i)%y_data = ''
      custom_loaders_list(i)%y_data_offset = 0
      custom_loaders_list(i)%z_data = ''
      custom_loaders_list(i)%z_data_offset = 0
      custom_loaders_list(i)%px_data = ''
      custom_loaders_list(i)%px_data_offset = 0
      custom_loaders_list(i)%px_data_given = .FALSE.
      custom_loaders_list(i)%py_data = ''
      custom_loaders_list(i)%py_data_offset = 0
      custom_loaders_list(i)%py_data_given = .FALSE.
      custom_loaders_list(i)%pz_data = ''
      custom_loaders_list(i)%pz_data_offset = 0
      custom_loaders_list(i)%pz_data_given = .FALSE.
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
      custom_loaders_list(i)%w_data = ''
      custom_loaders_list(i)%w_data_offset = 0
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      custom_loaders_list(i)%id_data = ''
      custom_loaders_list(i)%id_data_offset = 0
      custom_loaders_list(i)%id_data_given = .FALSE.
#endif
    ENDDO

  END SUBROUTINE setup_custom_loaders_list
 
  SUBROUTINE loader_debug_print

    INTEGER :: i
    CHARACTER(LEN=string_length) :: offset_as_str

    IF (rank /= 0) RETURN

    PRINT*, "*** Custom Loader Debug***"
    PRINT*, ""

    DO i = 1, n_custom_loaders
      CALL integer_as_string(i, offset_as_str)
      PRINT*, "** Loader ", TRIM(offset_as_str), " target: ", &
              TRIM(custom_loaders_list(i)%target_name)
      CALL integer_as_string(custom_loaders_list(i)%x_data_offset, offset_as_str)
      PRINT*, "x_data: ", TRIM(custom_loaders_list(i)%x_data), ":", TRIM(offset_as_str)
      CALL integer_as_string(custom_loaders_list(i)%y_data_offset, offset_as_str)
      PRINT*, "y_data: ", TRIM(custom_loaders_list(i)%y_data), ":", TRIM(offset_as_str)
      CALL integer_as_string(custom_loaders_list(i)%z_data_offset, offset_as_str)
      PRINT*, "z_data: ", TRIM(custom_loaders_list(i)%z_data), ":", TRIM(offset_as_str)
#if !defined(PER_SPECIES_WEIGHT) || defined(PHOTONS)
      CALL integer_as_string(custom_loaders_list(i)%w_data_offset, offset_as_str)
      PRINT*, "w_data: ", TRIM(custom_loaders_list(i)%w_data), ":", TRIM(offset_as_str)
#endif
      IF (custom_loaders_list(i)%px_data_given) THEN
        CALL integer_as_string(custom_loaders_list(i)%px_data_offset, offset_as_str)
        PRINT*, "px_data: ", custom_loaders_list(i)%px_data, TRIM(offset_as_str)
      ELSE
        PRINT*, "px_data: NULL"
      ENDIF
      IF (custom_loaders_list(i)%py_data_given) THEN
        CALL integer_as_string(custom_loaders_list(i)%py_data_offset, offset_as_str)
        PRINT*, "py_data: ", custom_loaders_list(i)%px_data, TRIM(offset_as_str)
      ELSE
        PRINT*, "py_data: NULL"
      ENDIF
      IF (custom_loaders_list(i)%pz_data_given) THEN
        CALL integer_as_string(custom_loaders_list(i)%pz_data_offset, offset_as_str)
        PRINT*, "pz_data: ", custom_loaders_list(i)%pz_data, TRIM(offset_as_str)
      ELSE
        PRINT*, "pz_data: NULL"
      ENDIF
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      IF (custom_loaders_list(i)%id_data_given) THEN
        CALL integer_as_string(custom_loaders_list(i)%id_data_offset, offset_as_str)
        PRINT*, "id_data: ", custom_loaders_list(i)%id_data, TRIM(offset_as_str)
      ELSE
        PRINT*, "id_data: NULL"
      ENDIF
#endif
      PRINT*, ""
    ENDDO

    PRINT*, "*** End Debug ***"
    PRINT*, ""

  END SUBROUTINE loader_debug_print

END MODULE deck_particles_block

