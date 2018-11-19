! Copyright (C) 2011-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
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

MODULE deck_subset_block

  USE strings_advanced
  USE utilities
  USE particle_id_hash_mod

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: subset_deck_initialise, subset_deck_finalise
  PUBLIC :: subset_block_start, subset_block_end
  PUBLIC :: subset_block_handle_element, subset_block_check

  INTEGER :: subset_id, current_block
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
  CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: subset_names
  INTEGER, DIMENSION(:), POINTER :: subset_blocks
  LOGICAL :: got_name
  INTEGER :: check_block = c_err_none

CONTAINS

  SUBROUTINE subset_deck_initialise

    INTEGER :: i

    current_block = 0
    IF (deck_state == c_ds_first) THEN
      n_subsets = 0
      ALLOCATE(subset_names(4))
      ALLOCATE(subset_blocks(4))
    ELSE
      DO i = 1, n_subsets
        ALLOCATE(subset_list(i)%use_species(n_species))
        subset_list(i)%use_species = .FALSE.
      END DO
    END IF

  END SUBROUTINE subset_deck_initialise



  SUBROUTINE subset_deck_finalise

    INTEGER :: i
    TYPE(subset), POINTER :: sub

    IF (deck_state == c_ds_first) THEN
      CALL setup_subsets

      DO i = 1, n_subsets
        subset_list(i)%name = subset_names(i)
      END DO
      DEALLOCATE(subset_names)
    ELSE
      DEALLOCATE(subset_blocks)
      DO i = 1, n_subsets
        sub => subset_list(i)

        sub%skip = (SUM(sub%skip_dir - 1) /= 0)

        ! Check for any spatial restrictions in place
        sub%space_restrictions = sub%use_x_min .OR. sub%use_x_max &
            .OR. sub%use_y_min .OR. sub%use_y_max &
            .OR. sub%use_z_min .OR. sub%use_z_max
        IF (sub%skip .AND. sub%space_restrictions) THEN
          IF (rank == 0) THEN
            PRINT*, 'Skip and spatial restrictions specified for ', &
                TRIM(sub%name), &
                ': field variables will not be trimmmed'
          END IF
          sub%space_restrictions = .FALSE.
        END IF
      END DO
    END IF

  END SUBROUTINE subset_deck_finalise



  SUBROUTINE subset_block_start

    current_block = current_block + 1
    got_name = .FALSE.
    IF (deck_state == c_ds_first) RETURN
    subset_id = subset_blocks(current_block)
    offset = 0

  END SUBROUTINE subset_block_start



  SUBROUTINE subset_block_end

    CHARACTER(LEN=8) :: id_string
    INTEGER :: io, iu

    IF (.NOT.got_name) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(current_block, id_string)
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Restriction block number ', TRIM(id_string), &
              ' has no "name" element.'
        END DO
      END IF

      check_block = c_err_missing_elements
    END IF

  END SUBROUTINE subset_block_end



  FUNCTION subset_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: io, iu, ispecies
    TYPE(subset), POINTER :: sub
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    TYPE(particle_id_hash), POINTER :: current_hash
#endif

    errcode = c_err_none
    IF (value == blank .OR. element == blank) RETURN

    IF (str_cmp(element, 'name')) THEN
      IF (got_name) THEN
        errcode = c_err_preset_element
        RETURN
      END IF
      got_name = .TRUE.
      IF (deck_state /= c_ds_first) RETURN
      CALL grow_array(subset_blocks, current_block)
      subset_blocks(current_block) = subset_number_from_name(value)
      RETURN
    END IF

    IF (deck_state == c_ds_first) RETURN

    sub => subset_list(subset_id)

    IF (str_cmp(element, 'random_fraction')) THEN
      sub%random_fraction = as_real_print(value, element, errcode)
      sub%use_random = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'gamma_min')) THEN
      sub%gamma_min = as_real_print(value, element, errcode)
      sub%use_gamma_min = .TRUE.
      sub%use_gamma = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'gamma_max')) THEN
      sub%gamma_max = as_real_print(value, element, errcode)
      sub%use_gamma_max = .TRUE.
      sub%use_gamma = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'x_min')) THEN
      sub%x_min = as_real_print(value, element, errcode)
      sub%use_x_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'x_max')) THEN
      sub%x_max = as_real_print(value, element, errcode)
      sub%use_x_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'y_min')) THEN
      sub%y_min = as_real_print(value, element, errcode)
      sub%use_y_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'y_max')) THEN
      sub%y_max = as_real_print(value, element, errcode)
      sub%use_y_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'z_min')) THEN
      sub%z_min = as_real_print(value, element, errcode)
      sub%use_z_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'z_max')) THEN
      sub%z_max = as_real_print(value, element, errcode)
      sub%use_z_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'px_min')) THEN
      sub%px_min = as_real_print(value, element, errcode)
      sub%use_px_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'px_max')) THEN
      sub%px_max = as_real_print(value, element, errcode)
      sub%use_px_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'py_min')) THEN
      sub%py_min = as_real_print(value, element, errcode)
      sub%use_py_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'py_max')) THEN
      sub%py_max = as_real_print(value, element, errcode)
      sub%use_py_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'pz_min')) THEN
      sub%pz_min = as_real_print(value, element, errcode)
      sub%use_pz_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'pz_max')) THEN
      sub%pz_max = as_real_print(value, element, errcode)
      sub%use_pz_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'weight_min')) THEN
      sub%weight_min = as_real_print(value, element, errcode)
      sub%use_weight_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'weight_max')) THEN
      sub%weight_max = as_real_print(value, element, errcode)
      sub%use_weight_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'charge_min')) THEN
      sub%charge_min = as_real_print(value, element, errcode)
      sub%use_charge_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'charge_max')) THEN
      sub%charge_max = as_real_print(value, element, errcode)
      sub%use_charge_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'mass_min')) THEN
      sub%mass_min = as_real_print(value, element, errcode)
      sub%use_mass_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'mass_max')) THEN
      sub%mass_max = as_real_print(value, element, errcode)
      sub%use_mass_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'id_min')) THEN
      sub%id_min = as_integer_print(value, element, errcode)
      sub%use_id_min = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'id_max')) THEN
      sub%id_max = as_integer_print(value, element, errcode)
      sub%use_id_max = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'dumpmask')) THEN
      sub%mask = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'skip')) THEN
      sub%skip_dir = as_integer_print(value, element, errcode) + 1
      RETURN
    END IF

    IF (str_cmp(element, 'skip_x')) THEN
      sub%skip_dir(1) = as_integer_print(value, element, errcode) + 1
      RETURN
    END IF

    IF (str_cmp(element, 'skip_y')) THEN
      sub%skip_dir(2) = as_integer_print(value, element, errcode) + 1
      RETURN
    END IF

    IF (str_cmp(element, 'skip_z')) THEN
      sub%skip_dir(3) = as_integer_print(value, element, errcode) + 1
      RETURN
    END IF

    IF (str_cmp(element, 'include_species')) THEN
      ispecies = as_integer_print(value, element, errcode)
      IF (errcode == c_err_none) THEN
        IF (ispecies > 0 .AND. ispecies <= n_species) THEN
          sub%use_species(ispecies) = .TRUE.
        ELSE
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Unable to apply subset to non existant ', &
                  'species ', ispecies
            END DO
          END IF
          errcode = c_err_bad_value
        END IF
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'persist_start_time') &
        .OR. str_cmp(element, 'persist_after_time')) THEN
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      sub%persistent = .TRUE.
      sub%persist_start_time = as_real_print(value, element, errcode)
      current_hash => id_registry%get_hash(sub%name)
      IF (ASSOCIATED(current_hash)) THEN
        CALL current_hash%init(1000)
      ELSE
        IF (rank == 0) PRINT*, 'Can only have 64 persistent subsets'
        errcode = c_err_bad_value
        RETURN
      END IF
#else
      errcode = c_err_pp_options_missing
      extended_error_string = '-DPARTICLE_ID'
#endif
      RETURN
    END IF

    IF (str_cmp(element, 'from_file')) THEN
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      sub%persistent = .TRUE.
      sub%filename = TRIM(value)
      sub%from_file = .TRUE.
      current_hash => id_registry%get_hash(sub%name)
      IF (ASSOCIATED(current_hash)) THEN
        CALL current_hash%init(1000)
      ELSE
        IF (rank == 0) PRINT*, 'Can only have 64 persistent subsets'
        errcode = c_err_bad_value
        RETURN
      END IF
#else
      errcode = c_err_pp_options_missing
      extended_error_string = '-DPARTICLE_ID'
#endif
      RETURN
    END IF

    IF (str_cmp(element, 'file_overrides_restart')) THEN
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      sub%add_after_restart = str_cmp(element, 'from_file_on_restart')
#else
      errcode = c_err_pp_options_missing
      extended_error_string = '-DPARTICLE_ID'
#endif
      RETURN
    END IF

    IF (str_cmp(element, 'sorted_file')) THEN
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
      sub%file_sorted = as_logical_print(value, element, errcode)
#else
      errcode = c_err_pp_options_missing
      extended_error_string = '-DPARTICLE_ID'
#endif
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION subset_block_handle_element



  FUNCTION subset_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = check_block

  END FUNCTION subset_block_check



  FUNCTION subset_number_from_name(name)

    CHARACTER(*), INTENT(IN) :: name
    INTEGER :: subset_number_from_name
    INTEGER :: i

    DO i = 1, n_subsets
      IF (str_cmp(name, subset_names(i))) THEN
        subset_number_from_name = i
        RETURN
      END IF
    END DO
    n_subsets = n_subsets + 1
    CALL grow_array(subset_names, n_subsets)
    subset_names(n_subsets) = TRIM(name)
    subset_number_from_name = n_subsets

  END FUNCTION subset_number_from_name



  SUBROUTINE setup_subsets

    INTEGER :: i

    ALLOCATE(subset_list(n_subsets))

    DO i = 1, n_subsets
      subset_list(i)%name = blank
      subset_list(i)%use_random     = .FALSE.
      subset_list(i)%use_gamma      = .FALSE.
      subset_list(i)%use_gamma_min  = .FALSE.
      subset_list(i)%use_gamma_max  = .FALSE.
      subset_list(i)%use_x_min      = .FALSE.
      subset_list(i)%use_x_max      = .FALSE.
      subset_list(i)%use_y_min      = .FALSE.
      subset_list(i)%use_y_max      = .FALSE.
      subset_list(i)%use_z_min      = .FALSE.
      subset_list(i)%use_z_max      = .FALSE.
      subset_list(i)%use_px_min     = .FALSE.
      subset_list(i)%use_px_max     = .FALSE.
      subset_list(i)%use_py_min     = .FALSE.
      subset_list(i)%use_py_max     = .FALSE.
      subset_list(i)%use_pz_min     = .FALSE.
      subset_list(i)%use_pz_max     = .FALSE.
      subset_list(i)%use_weight_min = .FALSE.
      subset_list(i)%use_weight_max = .FALSE.
      subset_list(i)%use_charge_min = .FALSE.
      subset_list(i)%use_charge_max = .FALSE.
      subset_list(i)%use_mass_min   = .FALSE.
      subset_list(i)%use_mass_max   = .FALSE.
      subset_list(i)%use_id_min     = .FALSE.
      subset_list(i)%use_id_max     = .FALSE.
      subset_list(i)%skip           = .FALSE.
      subset_list(i)%dump_field_grid = .FALSE.
      subset_list(i)%space_restrictions = .FALSE.
      subset_list(i)%skip_dir       = 1
      subset_list(i)%random_fraction = 0.0_num
      subset_list(i)%gamma_min  = -HUGE(1.0_num)
      subset_list(i)%gamma_max  =  HUGE(1.0_num)
      subset_list(i)%x_min      = -HUGE(1.0_num)
      subset_list(i)%x_max      =  HUGE(1.0_num)
      subset_list(i)%y_min      = -HUGE(1.0_num)
      subset_list(i)%y_max      =  HUGE(1.0_num)
      subset_list(i)%z_min      = -HUGE(1.0_num)
      subset_list(i)%z_max      =  HUGE(1.0_num)
      subset_list(i)%px_min     = -HUGE(1.0_num)
      subset_list(i)%px_max     =  HUGE(1.0_num)
      subset_list(i)%py_min     = -HUGE(1.0_num)
      subset_list(i)%py_max     =  HUGE(1.0_num)
      subset_list(i)%pz_min     = -HUGE(1.0_num)
      subset_list(i)%pz_max     =  HUGE(1.0_num)
      subset_list(i)%weight_min = -HUGE(1.0_num)
      subset_list(i)%weight_max =  HUGE(1.0_num)
      subset_list(i)%charge_min = -HUGE(1.0_num)
      subset_list(i)%charge_max =  HUGE(1.0_num)
      subset_list(i)%mass_min   = -HUGE(1.0_num)
      subset_list(i)%mass_max   =  HUGE(1.0_num)
      subset_list(i)%id_min     = -HUGE(1)
      subset_list(i)%id_max     =  HUGE(1)
      subset_list(i)%mask = c_io_always
      ALLOCATE(subset_list(i)%dumpmask(n_io_blocks,num_vars_to_dump))
      subset_list(i)%dumpmask = c_io_none
      subset_list(i)%persistent = .FALSE.
      subset_list(i)%persist_start_time = -1.0_num
      subset_list(i)%locked = .FALSE.
      subset_list(i)%from_file = .FALSE.
      subset_list(i)%file_sorted = .FALSE.
    END DO

  END SUBROUTINE setup_subsets

END MODULE deck_subset_block
