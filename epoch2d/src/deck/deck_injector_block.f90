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

MODULE deck_injector_block

  USE strings_advanced
  USE shunt
  USE injectors
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: injector_deck_initialise, injector_deck_finalise
  PUBLIC :: injector_block_start, injector_block_end
  PUBLIC :: injector_block_handle_element, injector_block_check

  TYPE(injector_block), POINTER :: working_injector
  LOGICAL :: boundary_set = .FALSE.
  INTEGER :: boundary

CONTAINS

  SUBROUTINE injector_deck_initialise

    IF (deck_state /= c_ds_first) RETURN

    NULLIFY(injector_x_min)
    NULLIFY(injector_x_max)
    NULLIFY(injector_y_min)
    NULLIFY(injector_y_max)

  END SUBROUTINE injector_deck_initialise



  SUBROUTINE injector_deck_finalise

  END SUBROUTINE injector_deck_finalise



  SUBROUTINE injector_block_start

    IF (deck_state == c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_injector)

  END SUBROUTINE injector_block_start



  SUBROUTINE injector_block_end

    IF (deck_state == c_ds_first) RETURN

    CALL attach_injector(working_injector)
    boundary_set = .FALSE.

  END SUBROUTINE injector_block_end



  FUNCTION injector_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, i
    INTEGER :: io, iu
    REAL(num) :: mult
    CHARACTER(LEN=string_length) :: mult_string

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'boundary')) THEN
      ! If the boundary has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary = as_boundary_print(value, element, errcode)
      boundary_set = .TRUE.
      CALL init_injector(boundary, working_injector)
      RETURN
    END IF

    IF (.NOT. boundary_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Cannot set injector properties before boundary is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'boundary'
      errcode = c_err_required_element_not_set
      RETURN
    END IF

    IF (str_cmp(element, 'species')) THEN
      working_injector%species = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_start')) THEN
      working_injector%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      working_injector%t_end = as_time_print(value, element, errcode)
      working_injector%has_t_end = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'use_flux_maxwellian')) THEN
      working_injector%use_flux_injector = as_logical_print(value, element, &
          errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'npart_per_cell') &
        .OR. str_cmp(element, 'nparticles_per_cell')) THEN
      working_injector%npart_per_cell = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density_min') .OR. str_cmp(element, 'minrho') &
        .OR. str_cmp(element, 'number_density_min')) THEN
      working_injector%density_min = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density_max') .OR. str_cmp(element, 'maxrho') &
        .OR. str_cmp(element, 'number_density_max')) THEN
      working_injector%density_max = as_real_print(value, element, errcode)
      IF (working_injector%density_max < 0.0_num) THEN
        working_injector%density_max = HUGE(1.0_num)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'density') &
        .OR. str_cmp(element, 'number_density')) THEN
      CALL initialise_stack(working_injector%density_function)
      CALL tokenize(value, working_injector%density_function, errcode)
      RETURN
    END IF

    mult_string = '* ev / kb'
    mult = ev / kb

    IF (str_cmp(element, 'temp') &
        .OR. str_cmp(element, 'temp_k') &
        .OR. str_cmp(element, 'temp_ev') &
        .OR. str_cmp(element, 'temperature') &
        .OR. str_cmp(element, 'temperature_k') &
        .OR. str_cmp(element, 'temperature_ev')) THEN
      DO i = 1, 3
        CALL initialise_stack(working_injector%temperature_function(i))
        CALL tokenize(value, working_injector%temperature_function(i), errcode)
        IF (str_cmp(element, 'temperature_ev') &
            .OR. str_cmp(element, 'temp_ev')) THEN
          CALL tokenize(mult_string, working_injector%temperature_function(i), &
                        errcode)
        END IF
      END DO
      RETURN
    END IF

    IF (str_cmp(element, 'temp_x') &
        .OR. str_cmp(element, 'temp_x_k') &
        .OR. str_cmp(element, 'temp_x_ev') &
        .OR. str_cmp(element, 'temperature_x') &
        .OR. str_cmp(element, 'temperature_x_k') &
        .OR. str_cmp(element, 'temperature_x_ev')) THEN
      i = 1
      CALL initialise_stack(working_injector%temperature_function(i))
      CALL tokenize(value, working_injector%temperature_function(i), errcode)
      IF (str_cmp(element, 'temperature_x_ev') &
          .OR. str_cmp(element, 'temp_x_ev')) THEN
        CALL tokenize(mult_string, working_injector%temperature_function(i), &
                      errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'temp_y') &
        .OR. str_cmp(element, 'temp_y_k') &
        .OR. str_cmp(element, 'temp_y_ev') &
        .OR. str_cmp(element, 'temperature_y') &
        .OR. str_cmp(element, 'temperature_y_k') &
        .OR. str_cmp(element, 'temperature_y_ev')) THEN
      i = 2
      CALL initialise_stack(working_injector%temperature_function(i))
      CALL tokenize(value, working_injector%temperature_function(i), errcode)
      IF (str_cmp(element, 'temperature_y_ev') &
          .OR. str_cmp(element, 'temp_y_ev')) THEN
        CALL tokenize(mult_string, working_injector%temperature_function(i), &
                      errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'temp_z') &
        .OR. str_cmp(element, 'temp_z_k') &
        .OR. str_cmp(element, 'temp_z_ev') &
        .OR. str_cmp(element, 'temperature_z') &
        .OR. str_cmp(element, 'temperature_z_k') &
        .OR. str_cmp(element, 'temperature_z_ev')) THEN
      i = 3
      CALL initialise_stack(working_injector%temperature_function(i))
      CALL tokenize(value, working_injector%temperature_function(i), errcode)
      IF (str_cmp(element, 'temperature_z_ev') &
          .OR. str_cmp(element, 'temp_z_ev')) THEN
        CALL tokenize(mult_string, working_injector%temperature_function(i), &
                      errcode)
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'drift_x') .OR. str_cmp(element, 'drift_px')) THEN
      i = 1
      CALL initialise_stack(working_injector%drift_function(i))
      CALL tokenize(value, working_injector%drift_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_y') .OR. str_cmp(element, 'drift_py')) THEN
      i = 2
      CALL initialise_stack(working_injector%drift_function(i))
      CALL tokenize(value, working_injector%drift_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_z') .OR. str_cmp(element, 'drift_pz')) THEN
      i = 3
      CALL initialise_stack(working_injector%drift_function(i))
      CALL tokenize(value, working_injector%drift_function(i), errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION injector_block_handle_element



  FUNCTION injector_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(injector_block), POINTER :: current
    INTEGER :: io, iu
    LOGICAL :: error

    use_injectors = .FALSE.
    error = .FALSE.
    errcode = c_err_none

    current => injector_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = .TRUE.
      use_injectors = .TRUE.
      current => current%next
    END DO

    current => injector_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = .TRUE.
      use_injectors = .TRUE.
      current => current%next
    END DO

    current => injector_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = .TRUE.
      use_injectors = .TRUE.
      current => current%next
    END DO

    current => injector_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = .TRUE.
      use_injectors = .TRUE.
      current => current%next
    END DO

    IF (error) THEN
      use_injectors = .FALSE.
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a species for every injector'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

  END FUNCTION injector_block_check

END MODULE deck_injector_block
