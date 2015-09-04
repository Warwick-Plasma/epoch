MODULE window

  USE boundary
  USE partlist
  USE evaluator

  IMPLICIT NONE

  LOGICAL, SAVE :: window_started
  REAL(num), ALLOCATABLE :: density(:), temperature(:,:), drift(:,:)
  REAL(num), SAVE :: window_shift_fraction

CONTAINS

  SUBROUTINE initialise_window

#ifdef PER_SPECIES_WEIGHT
    INTEGER :: ierr
#endif

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    ALLOCATE(density(-2:ny+3))
    ALLOCATE(temperature(-2:ny+3, 1:3))
    ALLOCATE(drift(-2:ny+3, 1:3))
    window_started = .FALSE.
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    ENDIF
    CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
#endif

  END SUBROUTINE initialise_window



  SUBROUTINE deallocate_window

    IF (ALLOCATED(density)) DEALLOCATE(density)
    IF (ALLOCATED(temperature)) DEALLOCATE(temperature)
    IF (ALLOCATED(drift)) DEALLOCATE(drift)

  END SUBROUTINE deallocate_window



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE shift_window(window_shift_cells)

    INTEGER, INTENT(IN) :: window_shift_cells
    INTEGER :: iwindow

    ! Shift the window round one cell at a time.
    ! Inefficient, but it works
    DO iwindow = 1, window_shift_cells
      CALL insert_particles

      ! Shift the box around
      x_grid_min = x_grid_min + dx
      x_grid_max = x_grid_max + dx
      x_grid_mins = x_grid_mins + dx
      x_grid_maxs = x_grid_maxs + dx
      x_grid_min_local = x_grid_min_local + dx
      x_grid_max_local = x_grid_max_local + dx
      x_min = x_min + dx
      x_max = x_max + dx
      x_min_local = x_min_local + dx
      x_max_local = x_max_local + dx

      x = x + dx
      x_global = x_global + dx
      xb_global = xb_global + dx

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

    IF (x_max_boundary) THEN
      DO j = -2, ny+3
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
    INTEGER(i8) :: ipart, npart_per_cell, n0
    REAL(num) :: cell_y_r, cell_frac_y, cy2
    INTEGER :: cell_y
    REAL(num), DIMENSION(-1:1) :: gy
    REAL(num) :: temp_local, drift_local, npart_frac
    REAL(num) :: weight_local

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    IF (nproc > 1) THEN
      IF (SIZE(density) /= ny+6) THEN
        DEALLOCATE(density, temperature, drift)
        ALLOCATE(density(-2:ny+3))
        ALLOCATE(temperature(-2:ny+3, 1:3))
        ALLOCATE(drift(-2:ny+3, 1:3))
      ENDIF
    ENDIF

    errcode = c_err_none

    DO ispecies = 1, n_species
      CALL create_empty_partlist(append_list)
      npart_per_cell = FLOOR(species_list(ispecies)%npart_per_cell, KIND=i8)
      npart_frac = species_list(ispecies)%npart_per_cell - npart_per_cell
      IF (npart_frac > 0) THEN
        n0 = 0
      ELSE
        n0 = 1
      ENDIF

      DO i = 1, 3
        DO iy = -2, ny+3
          temperature(iy,i) = evaluate_at_point( &
              species_list(ispecies)%temperature_function(i), nx, iy, errcode)
          drift(iy,i) = evaluate_at_point( &
              species_list(ispecies)%drift_function(i), nx, iy, errcode)
        ENDDO
      ENDDO
      DO iy = -2, ny+3
        density(iy) = evaluate_at_point( &
            species_list(ispecies)%density_function, nx, iy, errcode)
        IF (density(iy) > initial_conditions(ispecies)%density_max) &
            density(iy) = initial_conditions(ispecies)%density_max
      ENDDO

      DO iy = 1, ny
        IF (density(iy) < initial_conditions(ispecies)%density_min) CYCLE
        DO ipart = n0, npart_per_cell
          ! Place extra particle based on probability
          IF (ipart == 0) THEN
            IF (npart_frac < random()) CYCLE
          ENDIF
          CALL create_particle(current)
          current%part_pos(1) = x_grid_max + dx + (random() - 0.5_num) * dx
          current%part_pos(2) = y(iy) + (random() - 0.5_num) * dy

          ! Always use the triangle particle weighting for simplicity
          cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
          cell_y = FLOOR(cell_y_r + 0.5_num)
          cell_frac_y = REAL(cell_y, num) - cell_y_r
          cell_y = cell_y + 1

          cy2 = cell_frac_y**2
          gy(-1) = 0.5_num * (0.25_num + cy2 + cell_frac_y)
          gy( 0) = 0.75_num - cy2
          gy( 1) = 0.5_num * (0.25_num + cy2 - cell_frac_y)

          DO i = 1, 3
            temp_local = 0.0_num
            drift_local = 0.0_num
            DO isuby = -1, 1
              temp_local = temp_local + gy(isuby) * temperature(cell_y+isuby, i)
              drift_local = drift_local + gy(isuby) * drift(cell_y+isuby, i)
            ENDDO
            current%part_p(i) = momentum_from_temperature(&
                species_list(ispecies)%mass, temp_local, drift_local)
          ENDDO

          weight_local = 0.0_num
          DO isuby = -1, 1
            weight_local = weight_local + gy(isuby) * dx * dy &
                / species_list(ispecies)%npart_per_cell * density(cell_y+isuby)
          ENDDO
          current%weight = weight_local
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
    INTEGER :: window_shift_cells
#else
    INTEGER :: ierr
#endif

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    IF (.NOT. window_started) THEN
      IF (time >= window_start_time) THEN
        bc_field(c_bd_x_min) = bc_x_min_after_move
        bc_field(c_bd_x_max) = bc_x_max_after_move
        CALL setup_particle_boundaries
        IF (.NOT.ic_from_restart) window_shift_fraction = 0.0_num
        window_started = .TRUE.
      ENDIF
    ENDIF

    ! If we have a moving window then update the window position
    IF (window_started) THEN
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
    CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
#endif

  END SUBROUTINE moving_window

END MODULE window
