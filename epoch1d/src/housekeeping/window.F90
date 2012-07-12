MODULE window

  USE mpi
  USE boundary
  USE partlist
  USE particle_temperature
  USE evaluator

  IMPLICIT NONE

  LOGICAL, SAVE :: window_started

CONTAINS

  SUBROUTINE initialise_window

#ifdef PER_PARTICLE_WEIGHT
    INTEGER :: ispecies
#endif

    IF (.NOT. move_window) RETURN

#ifdef PER_PARTICLE_WEIGHT
    DO ispecies = 1, n_species
      species_list(ispecies)%density = &
          initial_conditions(ispecies)%density(nx)
      species_list(ispecies)%temperature = &
          initial_conditions(ispecies)%temp(nx,:)
    ENDDO
    window_started = .FALSE.
#else
    IF (rank .EQ. 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    ENDIF
    CALL MPI_ABORT(comm, errcode, errcode)
#endif

  END SUBROUTINE initialise_window



  SUBROUTINE deallocate_window

  END SUBROUTINE deallocate_window



#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE shift_window(window_shift_cells)

    INTEGER, INTENT(IN) :: window_shift_cells
    INTEGER :: iwindow

    ! Shift the window round one cell at a time.
    ! Inefficient, but it works
    DO iwindow = 1, window_shift_cells
      CALL insert_particles

      ! Shift the box around
      x_min = x_min + dx
      x_max = x_max + dx
      x_mins = x_mins + dx
      x_maxs = x_maxs + dx
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

    CALL shift_field(ex)
    CALL shift_field(ey)
    CALL shift_field(ez)

    CALL shift_field(bx)
    CALL shift_field(by)
    CALL shift_field(bz)

    CALL shift_field(jx)
    CALL shift_field(jy)
    CALL shift_field(jz)

    IF (x_max_boundary) THEN
      ! Fix incoming field cell.
      ex(nx)   = ex_x_max
      ex(nx+1) = ex_x_max
      ey(nx+1) = ey_x_max
      ez(nx+1) = ez_x_max
      ex(nx-1) = 0.5_num * (ex(nx-2) + ex(nx))
      ey(nx)   = 0.5_num * (ey(nx-1) + ey(nx+1))
      ez(nx)   = 0.5_num * (ez(nx-1) + ez(nx+1))
      bx(nx+1) = bx_x_max
      by(nx)   = by_x_max
      bz(nx)   = bz_x_max
      bx(nx)   = 0.5_num * (bx(nx-1) + bx(nx+1))
      by(nx-1) = 0.5_num * (by(nx-2) + by(nx))
      bz(nx-1) = 0.5_num * (bz(nx-2) + bz(nx))
    ENDIF

  END SUBROUTINE shift_fields



  SUBROUTINE shift_field(field)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: field
    INTEGER :: i

    ! Shift field to the left by one cell
    DO i = -2, nx+2
      field(i) = field(i+1)
    ENDDO

    CALL field_bc(field)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    TYPE(particle_list) :: append_list
    INTEGER :: ispecies, i
    INTEGER(i8) :: ipart, npart_per_cell, n0
    REAL(num) :: temp_local, npart_frac

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    errcode = c_err_none

    DO ispecies = 1, n_species
      CALL create_empty_partlist(append_list)
      npart_per_cell = AINT(species_list(ispecies)%npart_per_cell, KIND=i8)
      npart_frac = species_list(ispecies)%npart_per_cell - npart_per_cell
      IF (npart_frac .GT. 0) THEN
        n0 = 0
      ELSE
        n0 = 1
      ENDIF

      DO i = 1, 3
        species_list(ispecies)%temperature(i) = evaluate_at_point( &
            species_list(ispecies)%temperature_function(i), nx, errcode)
      ENDDO
      species_list(ispecies)%density = evaluate_at_point( &
          species_list(ispecies)%density_function, nx, errcode)

      DO ipart = n0, npart_per_cell
        ! Place extra particle based on probability
        IF (ipart .EQ. 0) THEN
          IF (npart_frac .LT. random()) CYCLE
        ENDIF
        ALLOCATE(current)
        CALL init_particle(current)
        current%part_pos = x_max + dx + (random() - 0.5_num) * dx

        DO i = 1, 3
          temp_local = species_list(ispecies)%temperature(i)
          current%part_p(i) = momentum_from_temperature(&
              species_list(ispecies)%mass, temp_local, 0.0_num)
        ENDDO

        current%weight = dx * species_list(ispecies)%density &
            / species_list(ispecies)%npart_per_cell
#ifdef PARTICLE_DEBUG
        current%processor = rank
        current%processor_at_t0 = rank
#endif
        CALL add_particle_to_partlist(append_list, current)
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
          IF (current%part_pos .LT. x_min - 0.5_num * dx) THEN
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

#ifdef PER_PARTICLE_WEIGHT
    REAL(num), SAVE :: window_shift_fraction
    REAL(num) :: window_shift_real
    INTEGER :: window_shift_cells
#endif

    IF (.NOT. move_window) RETURN

#ifdef PER_PARTICLE_WEIGHT
    IF (.NOT. window_started) THEN
      IF (time .GE. window_start_time) THEN
        bc_field(c_bd_x_min) = bc_x_min_after_move
        bc_field(c_bd_x_max) = bc_x_max_after_move
        CALL setup_particle_boundaries
        window_shift_fraction = 0.0_num
        window_started = .TRUE.
      ENDIF
    ENDIF

    ! If we have a moving window then update the window position
    IF (window_started) THEN
      window_shift_fraction = window_shift_fraction + dt * window_v_x / dx
      window_shift_cells = FLOOR(window_shift_fraction)
      ! Allow for posibility of having jumped two cells at once
      IF (window_shift_cells .GT. 0) THEN
        window_shift_real = REAL(window_shift_cells, num)
        IF (use_offset_grid) THEN
          window_shift = window_shift + window_shift_real * dx
        ENDIF
        CALL shift_window(window_shift_cells)
        CALL particle_bcs
        window_shift_fraction = window_shift_fraction - window_shift_real
      ENDIF
    ENDIF
#else
    IF (rank .EQ. 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    ENDIF
    CALL MPI_ABORT(comm, errcode, errcode)
#endif

  END SUBROUTINE moving_window

END MODULE window
