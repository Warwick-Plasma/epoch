! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

  LOGICAL, SAVE :: window_started
  REAL(num), ALLOCATABLE :: density(:), temperature(:,:), drift(:,:)
  REAL(num), SAVE :: window_shift_fraction
  INTEGER(i8), ALLOCATABLE :: indices(:)
  LOGICAL, ALLOCATABLE :: got_index(:)

CONTAINS

  SUBROUTINE initialise_window

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    ALLOCATE(density(0:ny+1))
    ALLOCATE(temperature(0:ny+1, 1:3))
    ALLOCATE(drift(0:ny+1, 1:3))
    ALLOCATE(indices(ny))
    ALLOCATE(got_index(ny))
    window_started = .FALSE.
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    ENDIF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE initialise_window



  SUBROUTINE deallocate_window

    IF (ALLOCATED(density)) DEALLOCATE(density)
    IF (ALLOCATED(temperature)) DEALLOCATE(temperature)
    IF (ALLOCATED(drift)) DEALLOCATE(drift)
    IF (ALLOCATED(indices)) DEALLOCATE(indices)
    IF (ALLOCATED(got_index)) DEALLOCATE(got_index)

  END SUBROUTINE deallocate_window



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE shift_window(window_shift_cells)

    INTEGER, INTENT(IN) :: window_shift_cells
    INTEGER :: iwindow, ix, iproc
    REAL(num) :: xb_min

    ! Shift the window round one cell at a time.
    ! Inefficient, but it works
    DO iwindow = 1, window_shift_cells
      CALL insert_particles

      ! Shift the box around
      x_grid_min = x_global(1) + dx
      xb_min = xb_global(1) + dx
      x_min = xb_min + dx * cpml_thickness

      ! Setup global grid
      DO ix = 1-ng, nx_global + ng
        x_global(ix) = x_grid_min + (ix - 1) * dx
        xb_global(ix) = xb_min + (ix - 1) * dx
      ENDDO
      x_grid_max = x_global(nx_global)
      x_max = xb_global(nx_global+1) - dx * cpml_thickness

      DO iproc = 0, nprocx-1
        x_grid_mins(iproc) = x_global(cell_x_min(iproc+1))
        x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
      ENDDO

      x_grid_min_local = x_grid_mins(x_coords)
      x_grid_max_local = x_grid_maxs(x_coords)

      x_min_local = x_grid_min_local + (cpml_x_min_offset - 0.5_num) * dx
      x_max_local = x_grid_max_local - (cpml_x_max_offset - 0.5_num) * dx

      x(1-ng:nx+ng) = x_global(nx_global_min-ng:nx_global_max+ng)
      xb(1-ng:nx+ng) = xb_global(nx_global_min-ng:nx_global_max+ng)

      CALL remove_particles

      ! Shift fields around
      CALL shift_fields
    ENDDO

  END SUBROUTINE shift_window



  SUBROUTINE shift_fields

    INTEGER :: j

    CALL shift_field(ex, ng)
    CALL shift_field(ey, ng)
    CALL shift_field(ez, ng)

    CALL shift_field(bx, ng)
    CALL shift_field(by, ng)
    CALL shift_field(bz, ng)

    CALL shift_field(jx, jng)
    CALL shift_field(jy, jng)
    CALL shift_field(jz, jng)

    IF (cpml_boundaries) THEN
      CALL shift_field(cpml_psi_eyx, ng)
      CALL shift_field(cpml_psi_ezx, ng)
      CALL shift_field(cpml_psi_byx, ng)
      CALL shift_field(cpml_psi_bzx, ng)

      CALL shift_field(cpml_psi_exy, ng)
      CALL shift_field(cpml_psi_ezy, ng)
      CALL shift_field(cpml_psi_bxy, ng)
      CALL shift_field(cpml_psi_bzy, ng)
    ENDIF

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
      ENDDO

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
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE shift_fields



  SUBROUTINE shift_field(field, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, j

    ! Shift field to the left by one cell
    DO j = 1-ng, ny+ng
    DO i = 1-ng, nx+ng-1
      field(i,j) = field(i+1,j)
    ENDDO
    ENDDO

    CALL field_bc(field, ng)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    TYPE(particle_list) :: append_list
    INTEGER :: ispecies, i, iy, isuby
    INTEGER(i8) :: ipart, icell, npart_per_cell, n0, npart_left
    INTEGER(i8) :: cell_index, next_index, idx
    REAL(num) :: cell_frac_y, cy2
    REAL(num), DIMENSION(-1:1) :: gy
    REAL(num) :: temp_local, drift_local, npart_frac
    REAL(num) :: weight_local, x0
    TYPE(parameter_pack) :: parameters

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    IF (nproc > 1) THEN
      IF (SIZE(density) /= ny+2) THEN
        DEALLOCATE(density, temperature, drift)
        ALLOCATE(density(0:ny+1))
        ALLOCATE(temperature(0:ny+1, 1:3))
        ALLOCATE(drift(0:ny+1, 1:3))
        ALLOCATE(indices(ny))
        ALLOCATE(got_index(ny))
      ENDIF
    ENDIF

    errcode = c_err_none

    DO ispecies = 1, n_species
      CALL create_empty_partlist(append_list)
      npart_per_cell = FLOOR(species_list(ispecies)%npart_per_cell, KIND=i8)
      npart_frac = species_list(ispecies)%npart_per_cell - npart_per_cell
      npart_left = INT(ny * npart_frac)

      parameters%pack_ix = nx
      DO i = 1, 3
        DO iy = 0, ny+1
          parameters%pack_iy = iy
          temperature(iy,i) = evaluate_with_parameters( &
              species_list(ispecies)%temperature_function(i), &
              parameters, errcode)
          drift(iy,i) = evaluate_with_parameters( &
              species_list(ispecies)%drift_function(i), parameters, errcode)
        ENDDO
      ENDDO
      DO iy = 0, ny+1
        parameters%pack_iy = iy
        density(iy) = evaluate_with_parameters( &
            species_list(ispecies)%density_function, parameters, errcode)
        IF (density(iy) &
                > species_list(ispecies)%initial_conditions%density_max) THEN
          density(iy) = species_list(ispecies)%initial_conditions%density_max
        ENDIF
      ENDDO

      DO icell = 1, ny
        got_index(icell) = .FALSE.
        indices(icell) = icell
      ENDDO

      next_index = ny
      DO icell = 1, npart_left
        cell_index = INT(random() * (next_index - 1)) + 1
        idx = indices(cell_index)
        indices(cell_index:next_index-1) = indices(cell_index+1:next_index)
        got_index(idx) = .TRUE.
        next_index = next_index - 1
      ENDDO

      idx = 1
      DO icell = 1, ny
        IF (got_index(icell)) THEN
          indices(idx) = icell
          idx = idx + 1
        ENDIF
      ENDDO

      idx = 1
      cell_index = 0
      next_index = indices(idx)
      x0 = x_grid_max + 0.5_num * dx
      DO iy = 1, ny
        cell_index = cell_index + 1
        IF (cell_index == next_index) THEN
          n0 = 0
          idx = idx + 1
          next_index = indices(idx)
        ELSE
          n0 = 1
        ENDIF

        IF (density(iy) &
                < species_list(ispecies)%initial_conditions%density_min) THEN
          CYCLE
        ENDIF

        DO ipart = n0, npart_per_cell
          ! Place extra particle based on probability
          CALL create_particle(current)
          cell_frac_y = 0.5_num - random()
          current%part_pos(1) = x0 + random() * dx
          current%part_pos(2) = y(iy) - cell_frac_y * dy

          ! Always use the triangle particle weighting for simplicity
          cy2 = cell_frac_y**2
          gy(-1) = 0.5_num * (0.25_num + cy2 + cell_frac_y)
          gy( 0) = 0.75_num - cy2
          gy( 1) = 0.5_num * (0.25_num + cy2 - cell_frac_y)

          DO i = 1, 3
            temp_local = 0.0_num
            drift_local = 0.0_num
            DO isuby = -1, 1
              temp_local = temp_local + gy(isuby) * temperature(iy+isuby, i)
              drift_local = drift_local + gy(isuby) * drift(iy+isuby, i)
            ENDDO
            current%part_p(i) = momentum_from_temperature(&
                species_list(ispecies)%mass, temp_local, drift_local)
          ENDDO

          weight_local = 0.0_num
          DO isuby = -1, 1
            weight_local = weight_local + gy(isuby) * density(iy+isuby)
          ENDDO

          current%weight = weight_local * dx * dy &
              / species_list(ispecies)%npart_per_cell
#ifdef PARTICLE_DEBUG
          current%processor = rank
          current%processor_at_t0 = rank
#endif
          CALL add_particle_to_partlist(append_list, current)
        ENDDO
      ENDDO

      CALL append_partlist(species_list(ispecies)%attached_list, append_list)
    ENDDO

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
            DEALLOCATE(current)
          ENDIF
          current => next
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE remove_particles
#endif



  SUBROUTINE moving_window

#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: window_shift_real
    INTEGER :: window_shift_cells, errcode
#endif

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    IF (.NOT. window_started) THEN
      IF (time >= window_start_time) THEN
        bc_field(c_bd_x_min) = bc_x_min_after_move
        bc_field(c_bd_x_max) = bc_x_max_after_move
        CALL setup_boundaries
        IF (.NOT.ic_from_restart) window_shift_fraction = 0.0_num
        window_started = .TRUE.
      ENDIF
    ENDIF

    ! If we have a moving window then update the window position
    IF (window_started) THEN
      IF (use_window_stack) window_v_x = evaluate(window_v_x_stack, errcode)
      window_shift_fraction = window_shift_fraction + dt * window_v_x / dx
      window_shift_cells = FLOOR(window_shift_fraction)
      ! Allow for posibility of having jumped two cells at once
      IF (window_shift_cells > 0) THEN
        window_shift_real = REAL(window_shift_cells, num)
        IF (use_offset_grid) THEN
          window_shift(1) = window_shift(1) + window_shift_real * dx
        ENDIF
        CALL shift_window(window_shift_cells)
        CALL particle_bcs
        window_shift_fraction = window_shift_fraction - window_shift_real
      ENDIF
    ENDIF
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    ENDIF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE moving_window

END MODULE window
