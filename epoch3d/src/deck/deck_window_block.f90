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

    INTEGER :: i, io, iu, bc(2)
    LOGICAL :: warn, warn_window, warn_no_t_end

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

    CALL check_injector_boundaries(warn_window, warn_no_t_end)

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

    IF (warn_window .OR. warn_no_t_end) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'You have specified injectors in conjunction with the ', &
                    'moving window.'
        IF (warn_window) THEN
          WRITE(io,*) 'The t_end time of one or more injectors exceeds the ', &
                      'window_start_time.'
          WRITE(io,*) 'These injectors will continue to run once the moving ', &
                      'window starts, but '
          WRITE(io,*) 'care should be taken when interpreting these results.'
          WRITE(io,*)
        END IF
        IF (warn_no_t_end) THEN
          WRITE(io,*) 'One or more injectors has no explicit end time.'
          WRITE(io,*) 'These injectors will be disabled once the moving ', &
                      'window starts.'
          WRITE(io,*)
        END IF
      END DO
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



  SUBROUTINE check_injector_boundaries(warn_window, warn_no_t_end)

    LOGICAL, INTENT(OUT) :: warn_window, warn_no_t_end
    LOGICAL, DIMENSION(2) :: warn_loc, warn_glob
    INTEGER :: ierr
    TYPE(injector_block), POINTER :: current
    LOGICAL :: got_t_end
    REAL(num) :: t_end

    warn_loc = .FALSE.

    IF (ASSOCIATED(injector_list)) THEN
      t_end = HUGE(1.0_num)
      got_t_end = .FALSE.
      current => injector_list
      DO WHILE(ASSOCIATED(current))
        IF (current%has_t_end) THEN
          got_t_end = .TRUE.
          IF (current%t_end < t_end) t_end = current%t_end
        ELSE
          warn_loc(2) = .TRUE.
        END IF
        current => current%next
      END DO

      IF (got_t_end) THEN
        IF (t_end > window_start_time) THEN
          warn_loc(1) = .TRUE.
        END IF
      END IF
    END IF

    CALL MPI_REDUCE(warn_loc, warn_glob, 2, MPI_LOGICAL, MPI_LOR, 0, comm, ierr)

    warn_window = warn_glob(1)
    warn_no_t_end = warn_glob(2)

  END SUBROUTINE check_injector_boundaries

END MODULE deck_window_block
