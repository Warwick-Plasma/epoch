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

    INTEGER :: ispecies

    IF (.NOT. move_window) RETURN

#ifdef PER_PARTICLE_WEIGHT
    DO ispecies = 1, n_species
      ALLOCATE(species_list(ispecies)%density(-2:ny+3,-2:nz+3))
      ALLOCATE(species_list(ispecies)%temperature(-2:ny+3,-2:nz+3, 1:3))
      species_list(ispecies)%density = &
          initial_conditions(ispecies)%density(nx,:,:)
      species_list(ispecies)%temperature = &
          initial_conditions(ispecies)%temp(nx,:,:,:)
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

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      IF (ASSOCIATED(species_list(ispecies)%density)) &
          DEALLOCATE(species_list(ispecies)%density)
      IF (ASSOCIATED(species_list(ispecies)%temperature)) &
          DEALLOCATE(species_list(ispecies)%temperature)
    ENDDO

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

    INTEGER :: j, k

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
      DO k = -2, nz+3
        DO j = -2, ny+3
          ! Fix up incoming field cell. A future version will use
          ! equilibrium field, rather than zero.
          ex(nx+1,j,k) = 0.0_num
          ey(nx+1,j,k) = 0.0_num
          ez(nx+1,j,k) = 0.0_num
          ex(nx,j,k)   = 0.5_num * (ex(nx-1,j,k) + ex(nx+1,j,k))
          ey(nx,j,k)   = 0.5_num * (ey(nx-1,j,k) + ey(nx+1,j,k))
          ez(nx,j,k)   = 0.5_num * (ez(nx-1,j,k) + ez(nx+1,j,k))
          bx(nx+1,j,k) = 0.0_num
          by(nx+1,j,k) = 0.0_num
          bz(nx+1,j,k) = 0.0_num
          bx(nx,j,k)   = 0.5_num * (bx(nx-1,j,k) + bx(nx+1,j,k))
          by(nx,j,k)   = 0.5_num * (by(nx-1,j,k) + by(nx+1,j,k))
          bz(nx,j,k)   = 0.5_num * (bz(nx-1,j,k) + bz(nx+1,j,k))
        ENDDO
      ENDDO
    ENDIF

    IF (x_max_boundary) THEN
      DO k = -2, nz+3
        DO j = -2, ny+3
          ! Fix incoming field cell.
          ex(nx,j,k)   = ex_x_max(j,k)
          ex(nx+1,j,k) = ex_x_max(j,k)
          ey(nx+1,j,k) = ey_x_max(j,k)
          ez(nx+1,j,k) = ez_x_max(j,k)
          ex(nx-1,j,k) = 0.5_num * (ex(nx-2,j,k) + ex(nx,j,k))
          ey(nx,j,k)   = 0.5_num * (ey(nx-1,j,k) + ey(nx+1,j,k))
          ez(nx,j,k)   = 0.5_num * (ez(nx-1,j,k) + ez(nx+1,j,k))
          bx(nx+1,j,k) = bx_x_max(j,k)
          by(nx,j,k)   = by_x_max(j,k)
          bz(nx,j,k)   = bz_x_max(j,k)
          bx(nx,j,k)   = 0.5_num * (bx(nx-1,j,k) + bx(nx+1,j,k))
          by(nx-1,j,k) = 0.5_num * (by(nx-2,j,k) + by(nx,j,k))
          bz(nx-1,j,k) = 0.5_num * (bz(nx-2,j,k) + bz(nx,j,k))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE shift_fields



  SUBROUTINE shift_field(field)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: field
    INTEGER :: i, j, k

    ! Shift field to the left by one cell
    DO k = -2, nz+3
      DO j = -2, ny+3
        DO i = -2, nx+2
          field(i,j,k) = field(i+1,j,k)
        ENDDO
      ENDDO
    ENDDO

    CALL field_bc(field)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    TYPE(particle_list) :: append_list
    INTEGER :: ispecies, i, iy, iz, isuby, isubz
    INTEGER(KIND=8) :: ipart, npart_per_cell, n0
    REAL(num) :: cell_y_r, cell_frac_y, cy2
    REAL(num) :: cell_z_r, cell_frac_z, cz2
    INTEGER :: cell_y, cell_z
    REAL(num), DIMENSION(-1:1) :: gy, gz
    REAL(num) :: temp_local, npart_frac
    REAL(num) :: weight_local

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    errcode = c_err_none

    DO ispecies = 1, n_species
      CALL create_empty_partlist(append_list)
      npart_per_cell = AINT(species_list(ispecies)%npart_per_cell, KIND=8)
      npart_frac = species_list(ispecies)%npart_per_cell - npart_per_cell
      IF (npart_frac .GT. 0) THEN
        n0 = 0
      ELSE
        n0 = 1
      ENDIF

      DO i = 1, 3
        DO iz = -2, nz+3
          DO iy = -2, ny+3
            species_list(ispecies)%temperature(iy,iz,i) = evaluate_at_point( &
                species_list(ispecies)%temperature_function(i), nx, &
                iy, iz, errcode)
          ENDDO
        ENDDO
      ENDDO
      DO iz = -2, nz+3
        DO iy = -2, ny+3
          species_list(ispecies)%density(iy,iz) = evaluate_at_point( &
              species_list(ispecies)%density_function, nx, iy, iz, errcode)
        ENDDO
      ENDDO

      DO iz = 1, nz
        DO iy = 1, ny
          DO ipart = n0, npart_per_cell
            ! Place extra particle based on probability
            IF (ipart .EQ. 0) THEN
              IF (npart_frac .LT. random()) CYCLE
            ENDIF
            ALLOCATE(current)
            CALL init_particle(current)
            current%part_pos(1) = x_max + dx + (random() - 0.5_num) * dx
            current%part_pos(2) = y(iy) + (random() - 0.5_num) * dy
            current%part_pos(3) = z(iz) + (random() - 0.5_num) * dz

            ! Always use the triangle particle weighting for simplicity
            cell_y_r = (current%part_pos(2) - y_min_local) / dy
            cell_y = FLOOR(cell_y_r + 0.5_num)
            cell_frac_y = REAL(cell_y, num) - cell_y_r
            cell_y = cell_y + 1

            cell_z_r = (current%part_pos(3) - z_min_local) / dz
            cell_z = FLOOR(cell_z_r + 0.5_num)
            cell_frac_z = REAL(cell_z, num) - cell_z_r
            cell_z = cell_z + 1

            cy2 = cell_frac_y**2
            gy(-1) = 0.5_num * (0.25_num + cy2 + cell_frac_y)
            gy( 0) = 0.75_num - cy2
            gy( 1) = 0.5_num * (0.25_num + cy2 - cell_frac_y)

            cz2 = cell_frac_z**2
            gz(-1) = 0.5_num * (0.25_num + cz2 + cell_frac_z)
            gz( 0) = 0.75_num - cz2
            gz( 1) = 0.5_num * (0.25_num + cz2 - cell_frac_z)

            DO i = 1, 3
              temp_local = 0.0_num
              DO isubz = -1, 1
                DO isuby = -1, 1
                  temp_local = temp_local + gy(isuby) * gz(isubz) &
                      * species_list(ispecies) &
                          %temperature(cell_y+isuby, cell_z+isubz, i)
                ENDDO
              ENDDO
              current%part_p(i) = momentum_from_temperature(&
                  species_list(ispecies)%mass, temp_local, 0.0_num)
            ENDDO

            weight_local = 0.0_num
            DO isubz = -1, 1
              DO isuby = -1, 1
                weight_local = weight_local &
                    + gy(isuby) * gz(isubz) * dx * dy * dz &
                    / species_list(ispecies)%npart_per_cell &
                    * species_list(ispecies)%density(cell_y+isuby, cell_z+isubz)
              ENDDO
            ENDDO
            current%weight = weight_local
#ifdef PARTICLE_DEBUG
            current%processor = rank
            current%processor_at_t0 = rank
#endif
            CALL add_particle_to_partlist(append_list, current)
          ENDDO
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
          IF (current%part_pos(1) .LT. x_min - 0.5_num * dx) THEN
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

    REAL(num), SAVE :: window_shift_fraction
    REAL(num) :: window_shift_real
    INTEGER :: window_shift_cells

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
          window_shift(1) = window_shift(1) + window_shift_real * dx
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
