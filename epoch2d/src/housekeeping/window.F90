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

MODULE window

  USE boundary
  USE partlist
  USE evaluator

  IMPLICIT NONE

  REAL(num), ALLOCATABLE :: density(:), temperature(:,:), drift(:,:)
  REAL(num), SAVE :: window_shift_fraction

CONTAINS

  SUBROUTINE initialise_window

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    ALLOCATE(density(0:ny+1))
    ALLOCATE(temperature(0:ny+1, 1:3))
    ALLOCATE(drift(0:ny+1, 1:3))
    window_started = .FALSE.
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    END IF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE initialise_window



  SUBROUTINE deallocate_window

    IF (use_window_stack) CALL deallocate_stack(window_v_x_stack)
    IF (ALLOCATED(density)) DEALLOCATE(density)
    IF (ALLOCATED(temperature)) DEALLOCATE(temperature)
    IF (ALLOCATED(drift)) DEALLOCATE(drift)

  END SUBROUTINE deallocate_window



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE shift_window(window_shift_cells)

    INTEGER, INTENT(IN) :: window_shift_cells
    INTEGER :: iwindow, ix
    REAL(num) :: xb_min

    ! Shift the window round # of window shift cells at a time.
    CALL insert_particles(window_shift_cells)

    x_grid_min = x_global(1) + window_shift_cells * dx
    xb_min = xb_global(1) + window_shift_cells * dx
    x_min = xb_min + window_shift_cells * dx * cpml_thickness

    ! Setup global grid
    DO ix = 1-ng, nx_global + ng
      x_global(ix) = x_grid_min + (ix - 1) * dx
      xb_global(ix) = xb_min + (ix - 1) * dx
    END DO
    x_grid_max = x_global(nx_global)
    x_max = xb_global(nx_global+1) - dx * cpml_thickness

    CALL setup_grid_x

    CALL remove_particles

    ! Shift fields around
    CALL shift_fields(window_shift_cells)

  END SUBROUTINE shift_window


  SUBROUTINE shift_fields(window_shift_cells)

    INTEGER :: j
    INTEGER, INTENT(IN) :: window_shift_cells

    CALL shift_field(ex, ng, window_shift_cells)
    CALL shift_field(ey, ng, window_shift_cells)
    CALL shift_field(ez, ng, window_shift_cells)

    CALL moving_window_field_bc(ex, ey, ez, ng, nx, ny)

    CALL shift_field(bx, ng, window_shift_cells)
    CALL shift_field(by, ng, window_shift_cells)
    CALL shift_field(bz, ng, window_shift_cells)

    CALL moving_window_field_bc(bx, by, bz, ng ,nx, ny)

    CALL shift_field(jx, jng, window_shift_cells)
    CALL shift_field(jy, jng, window_shift_cells)
    CALL shift_field(jz, jng, window_shift_cells)

    CALL moving_window_field_bc(jx, jy, jz, ng, nx, ny)

    IF (cpml_boundaries) THEN
      CALL shift_field(cpml_psi_eyx, ng, &
              window_shift_cells)
      CALL shift_field(cpml_psi_ezx, ng, &
              window_shift_cells)
      CALL shift_field(cpml_psi_byx, ng, &
              window_shift_cells)
      
      CALL moving_window_field_bc(cpml_psi_eyx, &
              cpml_psi_ezx, cpml_psi_byx, &
              ng, nx, ny)
      
      CALL shift_field(cpml_psi_bzx, ng, &
              window_shift_cells) 
      
      CALL field_bc(cpml_psi_bzx, ng)

      CALL shift_field(cpml_psi_exy, ng, &
              window_shift_cells)
      CALL shift_field(cpml_psi_ezy, ng, &
              window_shift_cells)
      CALL shift_field(cpml_psi_bxy, ng, &
              window_shift_cells)
      CALL moving_window_field_bc(cpml_psi_exy, &
              cpml_psi_ezy, cpml_psi_ezy, &
              ng, nx, ny)
      
      CALL shift_field(cpml_psi_bzy, ng)
      CALL field_bc(cpml_psi_bzy, ng)
    END IF

    IF (x_max_boundary) THEN
      DO j = 1-ng, ny+ng
        ! Fix incoming field cell.
        ex(nx,j)   = ex_x_max(j)
        ex(nx+1,j) = ex_x_max(j)
        ey(nx+1,j) = ey_x_max(j)
        ez(nx+1,j) = ez_x_max(j)
        ex(nx-1,j) = 0.5_num * (ex(nx-2,j) + ex(nx,j))
        ey(nx,j)   = 0.5_num * (ey(nx-1,j) + ey(nx+1,j))
        ez(nx,j)   = 0.5_num * (ez(nx-1,j) + ez(nx+1,j))
        bx(nx+1,j) = bx_x_max(j)
        by(nx,j)   = by_x_max(j)
        bz(nx,j)   = bz_x_max(j)
        bx(nx,j)   = 0.5_num * (bx(nx-1,j) + bx(nx+1,j))
        by(nx-1,j) = 0.5_num * (by(nx-2,j) + by(nx,j))
        bz(nx-1,j) = 0.5_num * (bz(nx-2,j) + bz(nx,j))
      END DO

      IF (cpml_boundaries) THEN
        DO j = 1-ng, ny+ng
          cpml_psi_eyx(nx:nx+1,j) = cpml_psi_eyx(nx,j)
          cpml_psi_ezx(nx:nx+1,j) = cpml_psi_ezx(nx,j)
          cpml_psi_byx(nx:nx+1,j) = cpml_psi_byx(nx,j)
          cpml_psi_bzx(nx:nx+1,j) = cpml_psi_bzx(nx,j)

          cpml_psi_exy(nx:nx+1,j) = cpml_psi_exy(nx,j)
          cpml_psi_ezy(nx:nx+1,j) = cpml_psi_ezy(nx,j)
          cpml_psi_bxy(nx:nx+1,j) = cpml_psi_bxy(nx,j)
          cpml_psi_bzy(nx:nx+1,j) = cpml_psi_bzy(nx,j)
        END DO
      END IF
    END IF

  END SUBROUTINE shift_fields


  SUBROUTINE shift_field(field, ng, window_shift_cells)

    INTEGER, INTENT(IN) :: ng, window_shift_cells
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, j

    ! Shift field to the left by window_shift_cells
    DO j = 1-ng, ny+ng
    DO i = 1-ng, nx+ng-window_shift_cells
      field(i,j) = field(i+window_shift_cells, j)
    END DO
    END DO
    !CALL field_bc(field, ng)

  END SUBROUTINE shift_field


  SUBROUTINE moving_window_field_bc(fieldx, fieldy, fieldz, ng, nx_local, &
        ny_local)


    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(INOUT) :: fieldx, fieldy, fieldz
    INTEGER, INTENT(IN) :: nx_local, ny_local
    INTEGER, DIMENSION(c_ndims) :: sizes, subsizes
    INTEGER :: basetype, sz, szmax, i, j, k, n
    REAL(num), ALLOCATABLE :: field_left(:)
    REAL(num), ALLOCATABLE :: temp(:)

    basetype = mpireal

    sizes(1) = nx_local + 2 * ng
    sizes(2) = ny_local + 2 * ng
    starts = 1

    szmax = 3 * sizes(1) * ng
    sz = 3 * sizes(2) * ng
    IF (sz > szmax) szmax = sz

    ALLOCATE(temp(szmax))

    subsizes(1) = ng
    subsizes(2) = sizes(2)

    sz = 3 * subsizes(1) * subsizes(2)

  
    ALLOCATE(field_left(szmax))
    n = 1
    DO k = 1, 3
    DO j = 1-ng, subsizes(2)-ng
    DO i = 1,ng
      select case(k)
      case(1)
      field_left(n) = fieldx(i,j)
      case(2)
      field_left(n) = fieldy(i,j)
      case(3)
      field_left(n) = fieldz(i,j)
      end select
      n = n + 1
    END DO
    END DO
    END DO

    CALL MPI_SENDRECV(field_left, sz, basetype, proc_x_min, &
        tag, temp, sz, basetype, proc_x_max, tag, comm, status, errcode)

    IF (.NOT. x_max_boundary .OR. bc_field(c_bd_x_max)==c_bc_periodic) THEN
      n = 1
      DO k = 1, 3
      DO j = 1-ng, subsizes(2)-ng
      DO i = nx_local+1, subsizes(1)+nx_local
        select case(k)
        case(1)
        fieldx(i,j) = temp(n)
        case(2)
        fieldy(i,j) = temp(n)
        case(3)
        fieldz(i,j) = temp(n)
        end select
        n = n + 1
      END DO
      END DO
      END DO
    END IF

    DEALLOCATE(field_left)
    DEALLOCATE(temp)

  END SUBROUTINE moving_window_field_bc


  SUBROUTINE insert_particles(window_shift_cells)

    TYPE(particle), POINTER :: current
    TYPE(particle_list) :: append_list
    INTEGER :: ispecies, i, iy, isuby, errcode
    INTEGER(i8) :: ipart, npart_per_cell, n_frac
    INTEGER, INTENT(IN) :: window_shift_cells
    REAL(num) :: cell_frac_y, cy2
    REAL(num), DIMENSION(-1:1) :: gy
    REAL(num), DIMENSION(c_ndirs) :: temp_local, drift_local
    REAL(num) :: npart_frac
    REAL(num) :: weight_local, x0, dmin, dmax, wdata
    TYPE(parameter_pack) :: parameters
    TYPE(particle_species), POINTER :: species
    REAL(num), DIMENSION(c_ndirs, 2) :: ranges

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    IF (nproc > 1) THEN
      IF (SIZE(density) /= ny+2) THEN
        DEALLOCATE(density, temperature, drift)
        ALLOCATE(density(0:ny+1))
        ALLOCATE(temperature(0:ny+1, 1:3))
        ALLOCATE(drift(0:ny+1, 1:3))
      END IF
    END IF

    errcode = c_err_none

    DO ispecies = 1, n_species
      species => species_list(ispecies)

      CALL create_empty_partlist(append_list)
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
      npart_frac = species%npart_per_cell - npart_per_cell

      dmin = species%initial_conditions%density_min
      dmax = species%initial_conditions%density_max

      parameters%pack_ix = nx
      DO i = 1, 3
        DO iy = 0, ny+1
          parameters%pack_iy = iy
          temperature(iy,i) = evaluate_with_parameters( &
              species%temperature_function(i), parameters, errcode)
          drift(iy,i) = evaluate_with_parameters( &
              species%drift_function(i), parameters, errcode)
        END DO
      END DO
      DO iy = 0, ny+1
        parameters%pack_iy = iy
        density(iy) = evaluate_with_parameters(species%density_function, &
            parameters, errcode)
        IF (density(iy) > dmax) density(iy) = dmax
        IF (density(iy) < dmin) density(iy) = 0.0_num
      END DO

      x0 = x_grid_max + 0.5_num * dx
      DO iy = 1, ny
        IF (density(iy) < dmin) CYCLE

        ! Place extra particle based on probability
        n_frac = 0
        IF (npart_frac > 0.0_num) THEN
          IF (random() < npart_frac) n_frac = 1
        END IF

        wdata = dx * dy / (npart_per_cell + n_frac)

        DO ipart = 1, window_shift_cells * (npart_per_cell + n_frac)
          CALL create_particle(current)
          cell_frac_y = 0.5_num - random()
          current%part_pos(1) = x0 + random() * window_shift_cells * dx
          current%part_pos(2) = y(iy) - cell_frac_y * dy

          ! Always use the triangle particle weighting for simplicity
          cy2 = cell_frac_y**2
          gy(-1) = 0.5_num * (0.25_num + cy2 + cell_frac_y)
          gy( 0) = 0.75_num - cy2
          gy( 1) = 0.5_num * (0.25_num + cy2 - cell_frac_y)

          temp_local = 0.0_num
          drift_local = 0.0_num
          DO i = 1, c_ndirs
            DO isuby = -1, 1
              temp_local(i) = temp_local(i) + gy(isuby) &
                  * temperature(iy+isuby, i)
              drift_local(i) = drift_local(i) + gy(isuby) &
                  * drift(iy+isuby, i)
            END DO
          END DO

          IF (species%ic_df_type == c_ic_df_thermal) THEN
            DO i = 1, c_ndirs
              current%part_p(i) = momentum_from_temperature(species%mass, &
                  temp_local(i), drift_local(i))
            END DO
          ELSE IF (species%ic_df_type == c_ic_df_relativistic_thermal) THEN
            current%part_p = momentum_from_temperature_relativistic(&
                species%mass, temp_local, species%fractional_tail_cutoff, &
                drift_local)
          ELSE IF (species%ic_df_type == c_ic_df_arbitrary) THEN
            parameters%use_grid_position = .FALSE.
            parameters%pack_pos = current%part_pos
            parameters%pack_iy = iy
            errcode = c_err_none
            CALL evaluate_with_parameters_to_array(species%dist_fn_range(1), &
                parameters, 2, ranges(1,:), errcode)
            CALL evaluate_with_parameters_to_array(species%dist_fn_range(2), &
                parameters, 2, ranges(2,:), errcode)
            CALL evaluate_with_parameters_to_array(species%dist_fn_range(3), &
                parameters, 2, ranges(3,:), errcode)

            CALL sample_from_deck_expression(current, &
                species%dist_fn, parameters, ranges, species%mass, drift)
          END IF

          weight_local = 0.0_num
          DO isuby = -1, 1
            weight_local = weight_local + gy(isuby) * density(iy+isuby)
          END DO

          current%weight = weight_local * wdata
#ifdef PARTICLE_DEBUG
          current%processor = rank
          current%processor_at_t0 = rank
#endif
          CALL add_particle_to_partlist(append_list, current)
        END DO
      END DO

      CALL append_partlist(species%attached_list, append_list)
    END DO

  END SUBROUTINE insert_particles



  SUBROUTINE remove_particles

    TYPE(particle), POINTER :: current, next
    INTEGER :: ispecies

    ! Only processors on the left need do anything
    IF (x_min_boundary) THEN
      DO ispecies = 1, n_species
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          next => current%next
          IF (current%part_pos(1) < x_min) THEN
            CALL remove_particle_from_partlist(&
                species_list(ispecies)%attached_list, current)
            CALL destroy_particle(current)
          END IF
          current => next
        END DO
      END DO
    END IF

  END SUBROUTINE remove_particles
#endif



  SUBROUTINE moving_window

#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: window_shift_real
    INTEGER :: window_shift_cells, errcode = 0
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
      ! Allow for posibility of having jumped two cells at once
      IF (window_shift_cells > 0) THEN
        window_shift_real = REAL(window_shift_cells, num)
        window_offset = window_offset + window_shift_real * dx
        CALL shift_window(window_shift_cells)
        CALL setup_bc_lists
        CALL particle_bcs
        window_shift_fraction = window_shift_fraction - window_shift_real
      END IF
    END IF
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    END IF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE moving_window

END MODULE window
