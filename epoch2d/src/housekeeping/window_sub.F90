! Copyright (C) 2009-2019 University of Warwick
! Copyright (C) 2020 Juelich Supercomputing Center
!                    Forschungszentrum Juelich GmbH
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

SUBMODULE (window) window_sub
  ! This SUBMODULE is used to circumvent the circular module dependency between
  ! diagnostics and window (see comments therein)
  USE diagnostics

CONTAINS

  MODULE PROCEDURE moving_window
     USE diagnostics

#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: window_shift_real, window_shift_steps
    INTEGER :: window_shift_cells, errcode = 0
    INTEGER :: i, nchunks, nremainder
    logical :: print_arrays(1:SIZE(file_prefixes))
#endif

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    IF (.NOT. window_started) THEN
      IF (time >= window_start_time .AND. time < window_stop_time) THEN
        bc_field(c_bd_x_min) = bc_x_min_after_move
        bc_field(c_bd_x_max) = bc_x_max_after_move
        bc_field(c_bd_y_min) = bc_y_min_after_move
        bc_field(c_bd_y_max) = bc_y_max_after_move
        CALL setup_boundaries
        IF (.NOT.ic_from_restart) window_shift_fraction = 0.0_num
        window_started = .TRUE.
      END IF
    END IF

    ! If we have a moving window then update the window position
    IF (window_started) THEN
      IF (time >= window_stop_time) RETURN
      IF (use_window_stack) window_v_x = evaluate(window_v_x_stack, errcode)
      IF (window_v_x <= 0.0_num) RETURN
      window_shift_fraction = window_shift_fraction + dt * window_v_x / dx
      window_shift_cells = FLOOR(window_shift_fraction)
      ! CHECK FOR IO TAKING PLACE IN NEXT STEP...
      print_arrays = .false.
      DO i = 1, SIZE(file_prefixes)
        CALL io_test(i, future_step, print_arrays(i), force_dump, prefix_first_call)
      END DO

      ! Allow for posibility of having jumped two cells at once
      IF ( (window_shift_cells > ng - 1) .OR. ANY(print_arrays) ) THEN
        window_shift_real = REAL(window_shift_cells, num)
        window_offset = window_offset + window_shift_real * dx
        nremainder = MOD(window_shift_cells, ng)
        DO i = ng, window_shift_cells, ng  ! CHECK IF THIS LOOP IS CALLED IF window_shift_cells < ng
          CALL shift_window(ng)
        END DO
        IF (ANY(print_arrays)) then
          CALL shift_window(nremainder)
          nremainder = 0
        END IF
        CALL setup_bc_lists
        CALL particle_bcs
        window_shift_fraction = window_shift_fraction - window_shift_real &
                                + REAL(nremainder, num)
      END IF
    END IF
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    END IF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END PROCEDURE moving_window

END SUBMODULE window_sub
