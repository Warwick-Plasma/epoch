! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2010 Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE deck_hy_laser_block

  USE strings_advanced
  USE hy_laser
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: hy_laser_deck_initialise, hy_laser_deck_finalise
  PUBLIC :: hy_laser_block_start, hy_laser_block_end
  PUBLIC :: hy_laser_block_handle_element, hy_laser_block_check

  TYPE(hy_laser_block), POINTER :: working_laser
  LOGICAL :: hy_boundary_set = .FALSE.
  INTEGER :: hy_boundary

CONTAINS

  SUBROUTINE hy_laser_deck_initialise

  END SUBROUTINE hy_laser_deck_initialise



  SUBROUTINE hy_laser_deck_finalise

  END SUBROUTINE hy_laser_deck_finalise



  SUBROUTINE hy_laser_block_start

    IF (deck_state == c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_laser)
    working_laser%use_time_function = .FALSE.
    working_laser%use_profile_function = .FALSE.
    working_laser%use_omega_function = .FALSE.

  END SUBROUTINE hy_laser_block_start



  SUBROUTINE hy_laser_block_end

    IF (deck_state == c_ds_first) RETURN

    ! If we have user defined energy and user defined weight, then the laser
    ! parameters intensity, omega and eficiency will be unused.
    IF (working_laser%e_dist == e_dist_mono_weight .AND. &
        working_laser%mean == c_mean_E_val) working_laser%ignore_las = .TRUE.

    CALL attach_hy_laser(working_laser)
    hy_boundary_set = .FALSE.

  END SUBROUTINE hy_laser_block_end



  FUNCTION hy_laser_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dummy
    INTEGER :: io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'boundary')) THEN
      ! If the hy_boundary has already been set, simply ignore further calls
      IF (hy_boundary_set) RETURN
      hy_boundary = as_boundary_print(value, element, errcode)
      hy_boundary_set = .TRUE.
      CALL init_hy_laser(hy_boundary, working_laser)
      RETURN
    END IF

    IF (.NOT. hy_boundary_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Cannot set hy_laser properties before boundary is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'hy_boundary'
      errcode = c_err_required_element_not_set
      RETURN
    END IF

    IF (str_cmp(element, 'ppc')) THEN
      working_laser%ppc = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'profile')) THEN
      CALL initialise_stack(working_laser%profile_function)
      CALL tokenize(value, working_laser%profile_function, errcode)
      working_laser%profile = 0.0_num
      working_laser%use_profile_function = .TRUE.
      CALL hy_laser_update_profile(working_laser)
      IF (.NOT. working_laser%profile_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%profile_function)
        working_laser%use_profile_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'omega')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_omega
      working_laser%use_omega_function = .TRUE.
      CALL hy_laser_update_omega(working_laser)
      IF (.NOT. working_laser%omega_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%omega_function)
        working_laser%use_omega_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'frequency')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_freq
      working_laser%use_omega_function = .TRUE.
      CALL hy_laser_update_omega(working_laser)
      IF (.NOT. working_laser%omega_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%omega_function)
        working_laser%use_omega_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'lambda') .OR. str_cmp(element, 'wavelength')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_lambda
      working_laser%use_omega_function = .TRUE.
      CALL hy_laser_update_omega(working_laser)
      IF (.NOT. working_laser%omega_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%omega_function)
        working_laser%use_omega_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 't_profile')) THEN
      working_laser%use_time_function = .TRUE.
      CALL initialise_stack(working_laser%time_function)
      CALL tokenize(value, working_laser%time_function, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_laser%time_function, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'intensity')) THEN
      ! [W/m²]
      working_laser%intensity = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_start')) THEN
      working_laser%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      working_laser%t_end = as_time_print(value, element, errcode)
      working_laser%has_t_end = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'profile_min')) THEN
      working_laser%profile_min = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'efficiency')) THEN
      working_laser%efficiency = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'species')) THEN
      working_laser%species = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mean_energy')) THEN
      IF (str_cmp(value, 'a0')) THEN
        working_laser%mean = c_mean_a0
      ELSE IF (str_cmp(value, 'Wilks') .OR. str_cmp(value, 'wilks')) THEN
        working_laser%mean = c_mean_wilks
      ELSE IF (str_cmp(value, 'E_val') .OR. str_cmp(value, 'e_val')) THEN
        working_laser%mean = c_mean_E_val
      ELSE
        extended_error_string = 'Unrecognised mean energy model'
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'energy_dist')) THEN
      IF (str_cmp(value, 'exp')) THEN
        working_laser%e_dist = e_dist_exp
      ELSE IF (str_cmp(value, 'mono')) THEN
        working_laser%e_dist = e_dist_mono
      ELSE IF (str_cmp(value, 'tophat') .OR. str_cmp(value, 'top_hat')) THEN
        working_laser%e_dist = e_dist_tophat
      ELSE IF (str_cmp(value, 'exp_weight')) THEN
        working_laser%e_dist = e_dist_exp_weight
      ELSE IF (str_cmp(value, 'mono_weight')) THEN
        working_laser%e_dist = e_dist_mono_weight
      ELSE IF (str_cmp(value, 'mono_las_weight')) THEN
        working_laser%e_dist = e_dist_mono_las_weight
      ELSE
        extended_error_string = 'Unrecognised energy distribution model'
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'angular_dist')) THEN
      IF (str_cmp(value, 'uniform')) THEN
        working_laser%ang_dist = c_ang_uniform
      ELSE IF (str_cmp(value, 'cos')) THEN
        working_laser%ang_dist = c_ang_cos
      ELSE IF (str_cmp(value, 'beam')) THEN
        working_laser%ang_dist = c_ang_beam
      ELSE
        extended_error_string = 'Unrecognised angular distribution model'
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'mono_weight')) THEN
      working_laser%user_weight = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mean_E')) THEN
      working_laser%user_mean_KE = as_real_print(value, element, errcode) - mc2
      RETURN
    END IF

    IF (str_cmp(element, 'mean_KE')) THEN
      working_laser%user_mean_KE = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'las_weight_KE')) THEN
      working_laser%las_weight_KE = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'theta_max')) THEN
      working_laser%user_theta_max = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'cos_n_power')) THEN
      working_laser%cos_n_power = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'top_hat_L')) THEN
      working_laser%top_hat_L = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mean_mult')) THEN
      working_laser%mean_mult = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_sheng_dir')) THEN
      working_laser%use_sheng_dir = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sheng_angle')) THEN
      working_laser%sheng_angle = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_moore_max')) THEN
      working_laser%use_moore_max = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'theta')) THEN
      working_laser%theta_mean = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'phi')) THEN
      working_laser%phi_mean = as_real_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION hy_laser_block_handle_element



  FUNCTION hy_laser_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(hy_laser_block), POINTER :: current
    INTEGER :: error, io, iu

    errcode = c_err_none

    error = 0
    current => hy_laser_x_min
    DO WHILE(ASSOCIATED(current))
      error = single_hy_laser_check(current, error)
      current => current%next
    END DO

    current => hy_laser_x_max
    DO WHILE(ASSOCIATED(current))
      error = single_hy_laser_check(current, error)
      current => current%next
    END DO

    current => hy_laser_y_min
    DO WHILE(ASSOCIATED(current))
      error = single_hy_laser_check(current, error)
      current => current%next
    END DO

    current => hy_laser_y_max
    DO WHILE(ASSOCIATED(current))
      error = single_hy_laser_check(current, error)
      current => current%next
    END DO

    current => hy_laser_z_min
    DO WHILE(ASSOCIATED(current))
      error = single_hy_laser_check(current, error)
      current => current%next
    END DO

    current => hy_laser_z_max
    DO WHILE(ASSOCIATED(current))
      error = single_hy_laser_check(current, error)
      current => current%next
    END DO

    IF (IAND(error, 2**0) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "ppc" for every hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**1) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "species" for every hy_laser to inject to.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**2) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must set a model for mean electron energy with', &
              ' "mean_energy" for every'
          WRITE(io,*) 'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**3) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must set a model for the energy distribution with', &
              ' "energy_dist" for every'
          WRITE(io,*) 'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**4) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must set a model for the angular distribution with', &
              ' "angular_dist" for every'
          WRITE(io,*) 'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**5) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Top-hat fractional-width L must satisfy 0 < L < 1.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 2**6) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Incident (Sheng) angle cannot exceed pi/2 radians.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 2**7) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'mean_mult must be greater than zero.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 2**8) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Electron kinetic energy must be greater than zero.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 2**9) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Electron weight must be greater than zero.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 2**10) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "lambda" or "omega" for every hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**11) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "intensity" for every hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**12) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "efficiency" (laser energy to electron', &
              ' energy) for every'
          WRITE(io,*)'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2**13) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Efficiency cannot exceed 1.0 (100% conversion from', &
              ' laser energy to e- energy).'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 2**14) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Only user-defined weights are compatable with the', &
             ' -DPER_SPECIES_WEIGHT flag'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

  END FUNCTION hy_laser_block_check



  FUNCTION single_hy_laser_check(current, error)

    TYPE(hy_laser_block), POINTER :: current
    INTEGER :: error, single_hy_laser_check

    IF (current%ppc < 0) error = IOR(error, 2**0)
    IF (current%species < 0) error = IOR(error, 2**1)
    IF (current%mean < 0) error = IOR(error, 2**2)
    IF (current%e_dist < 0) error = IOR(error, 2**3)
    IF (current%ang_dist < 0) error = IOR(error, 2**4)
    IF (ABS(current%top_hat_L-0.5_num) >= 0.5_num) error = IOR(error, 2**5)
    IF (ABS(current%sheng_angle) > 0.5_num*pi) error = IOR(error, 2**6)
    IF (current%mean_mult < 0.0_num) error = IOR(error, 2**7)
    IF (current%mean == c_mean_E_val) THEN
      IF (current%user_mean_KE < 0.0_num) error = IOR(error, 2**8)
    END IF
    IF (current%e_dist == e_dist_mono_weight) THEN
      IF (current%user_weight < 0.0_num) error = IOR(error, 2**9)
    END IF
    IF (.NOT. current%ignore_las) THEN
      IF (current%omega < 0.0_num) error = IOR(error, 2**10)
      IF (current%intensity < 0.0_num) error = IOR(error, 2**11)
      IF (current%efficiency < 0.0_num) error = IOR(error, 2**12)
      IF (current%efficiency > 1.0_num) error = IOR(error, 2**13)
    END IF

#ifdef PER_SPECIES_WEIGHT
    IF (.NOT. current%e_dist == e_dist_mono_weight) error = IOR(error, 2**14)
#endif

   ! Output error-code
   single_hy_laser_check = error

  END FUNCTION

END MODULE deck_hy_laser_block
