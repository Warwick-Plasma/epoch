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

MODULE deck_window_block

  USE strings_advanced

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: window_deck_initialise, window_deck_finalise
  PUBLIC :: window_block_start, window_block_end
  PUBLIC :: window_block_handle_element, window_block_check

CONTAINS

  SUBROUTINE window_deck_initialise

    window_stop_time = c_largest_number

  END SUBROUTINE window_deck_initialise



  SUBROUTINE window_deck_finalise

    INTEGER :: i, io, iu, iwarn, bc(2)
    LOGICAL :: warn

    IF (.NOT.move_window) RETURN

    need_random_state = .TRUE.

    IF (deck_state == c_ds_first) RETURN

    IF (bc_x_min_after_move == c_bc_null) &
        bc_x_min_after_move = bc_field(c_bd_x_min)
    IF (bc_x_max_after_move == c_bc_null) &
        bc_x_max_after_move = bc_field(c_bd_x_max)
    IF (bc_y_min_after_move == c_bc_null) &
        bc_y_min_after_move = bc_field(c_bd_y_min)
    IF (bc_y_max_after_move == c_bc_null) &
        bc_y_max_after_move = bc_field(c_bd_y_max)
    IF (bc_z_min_after_move == c_bc_null) &
        bc_z_min_after_move = bc_field(c_bd_z_min)
    IF (bc_z_max_after_move == c_bc_null) &
        bc_z_max_after_move = bc_field(c_bd_z_max)

    CALL check_injector_boundaries(iwarn)

    IF (rank /= 0) RETURN

    ! Issue warnings about unsupported boundary conditions

    bc(1) = bc_x_min_after_move
    bc(2) = bc_x_max_after_move

    warn = .FALSE.
    DO i = 1, 2
      IF (bc(i) == c_bc_simple_laser &
          .OR. bc(i) == c_bc_cpml_laser &
          .OR. bc(i) == c_bc_cpml_outflow) THEN
        warn = .TRUE.
      END IF
    END DO

    IF (warn) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'You have specified lasers and/or CPML boundary ', &
                    'conditions for an X boundary'
        WRITE(io,*) 'after the moving window begins. These boundary ', &
                    'conditions are not compatible'
        WRITE(io,*) 'with moving windows and are unlikely to give correct ', &
                    'results.'
        WRITE(io,*)
      END DO
    END IF

    warn = .FALSE.
    DO i = 2, 2*c_ndims
      IF (bc_field(i) == c_bc_simple_laser &
          .OR. bc_field(i) == c_bc_cpml_laser &
          .OR. bc_field(i) == c_bc_cpml_outflow) THEN
        warn = .TRUE.
      END IF
    END DO

    IF (warn) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'You have specified lasers and/or CPML boundary ', &
                    'conditions for the Y or'
        WRITE(io,*) 'Z boundaries. These boundary conditions are not fully ', &
                    'compatible with moving'
        WRITE(io,*) 'windows and might give incorrect results.'
        WRITE(io,*)
      END DO
    END IF

    IF (n_custom_loaders > 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'You have specified particle loading from file in ', &
                    'conjunction with moving'
        WRITE(io,*) 'windows. The file contents will be ignored for new ', &
                    'particles entering the'
        WRITE(io,*) 'domain once the window begins to move.'
        WRITE(io,*)
      END DO
    END IF

    IF (iwarn == 1) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'You have specified injectors in conjunction with the ', &
                    'moving window.'
        WRITE(io,*) 'The t_end time of the injectors exceeds the ', &
                    'window_stop_time. The injectors'
        WRITE(io,*) 'will be disabled once the moving window starts.'
        WRITE(io,*)
      END DO
    ELSE IF (iwarn == 2) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'You have specified injectors in conjunction with the ', &
                    'moving window.'
        WRITE(io,*) 'These can only be used if they are explicitly ', &
                    'disabled before the window'
        WRITE(io,*) 'start time.'
        WRITE(io,*)
      END DO
      CALL abort_code(c_err_bad_value)
    END IF

  END SUBROUTINE window_deck_finalise



  SUBROUTINE window_block_start

  END SUBROUTINE window_block_start



  SUBROUTINE window_block_end

  END SUBROUTINE window_block_end



  FUNCTION window_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    LOGICAL, SAVE :: alloc = .FALSE.

    errcode = c_err_none
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'move_window')) THEN
      move_window = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'window_v_x')) THEN
      window_expression = .TRUE.
      IF (window_v_x_stack%init) CALL deallocate_stack(window_v_x_stack)
      CALL initialise_stack(window_v_x_stack)
      alloc = .TRUE.
      CALL tokenize(value, window_v_x_stack, errcode)
      ! evaluate it once to check that it's a valid block
      window_v_x = evaluate(window_v_x_stack, errcode)
      use_window_stack = window_v_x_stack%is_time_varying
      IF (.NOT.use_window_stack) THEN
        CALL deallocate_stack(window_v_x_stack)
        alloc = .FALSE.
      END IF
      window_expression = .FALSE.
      RETURN
    END IF

    IF (str_cmp(element, 'window_start_time')) THEN
      window_start_time = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'window_stop_time')) THEN
      window_stop_time = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_x_min_after_move') &
        .OR. str_cmp(element, 'xbc_left_after_move')) THEN
      bc_x_min_after_move = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_x_max_after_move') &
        .OR. str_cmp(element, 'xbc_right_after_move')) THEN
      bc_x_max_after_move = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_y_min_after_move') &
        .OR. str_cmp(element, 'ybc_down_after_move')) THEN
      bc_y_min_after_move = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_y_max_after_move') &
        .OR. str_cmp(element, 'ybc_up_after_move')) THEN
      bc_y_max_after_move = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_z_min_after_move') &
        .OR. str_cmp(element, 'zbc_back_after_move')) THEN
      bc_z_min_after_move = as_bc_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bc_z_max_after_move') &
        .OR. str_cmp(element, 'zbc_front_after_move')) THEN
      bc_z_max_after_move = as_bc_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION window_block_handle_element



  FUNCTION window_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION window_block_check



  SUBROUTINE check_injector_boundaries(iwarn)

    INTEGER, INTENT(OUT) :: iwarn
    INTEGER :: ierr, warn_buf(2), warn_sum(2)
    TYPE(injector_block), POINTER :: current
    LOGICAL :: got_t_end
    REAL(num) :: t_end

    iwarn = 0

    IF (ASSOCIATED(injector_list)) THEN
      t_end = HUGE(1.0_num)
      got_t_end = .FALSE.
      current => injector_list
      DO WHILE(ASSOCIATED(current))
        IF (current%has_t_end) THEN
          got_t_end = .TRUE.
          IF (current%t_end < t_end) t_end = current%t_end
        END IF
        current => current%next
      END DO

      IF (got_t_end) THEN
        IF (t_end > window_start_time) THEN
          iwarn = 1
        END IF
      ELSE
        iwarn = 2
      END IF
    END IF

    warn_buf(:) = 0
    IF (iwarn > 0) warn_buf(iwarn) = 1
    CALL MPI_REDUCE(warn_buf, warn_sum, 2, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    IF (warn_sum(2) > 0) THEN
      iwarn = 2
    ELSE IF (warn_sum(1) > 0) THEN
      iwarn = 1
    ELSE
      iwarn = 0
    END IF

  END SUBROUTINE check_injector_boundaries

END MODULE deck_window_block
