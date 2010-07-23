MODULE window

  USE helper
  USE boundary
  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE allocate_window

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      ALLOCATE(particle_species(ispecies)%density(-2:ny+2))
      ALLOCATE(particle_species(ispecies)%temperature(-2:ny+2, 1:3))
    ENDDO

  END SUBROUTINE allocate_window



  SUBROUTINE deallocate_window

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      IF (ASSOCIATED(particle_species(ispecies)%density)) &
          DEALLOCATE(particle_species(ispecies)%density)
      IF (ASSOCIATED(particle_species(ispecies)%temperature)) &
          DEALLOCATE(particle_species(ispecies)%temperature)
    ENDDO

  END SUBROUTINE deallocate_window



  SUBROUTINE shift_window

    INTEGER :: iwindow

    ! Shift the window round one cell at a time.
    ! Inefficient, but it works
    DO iwindow = 1, FLOOR(window_shift_fraction)
      CALL insert_particles
      ! Shift the box around
      x_mins = x_mins+dx
      x_min_local = x_min_local+dx
      x_min = x_min+dx
      x_global = x_global+dx
      x = x+dx

      x_maxs = x_maxs+dx
      x_max_local = x_max_local+dx
      x_max = x_max+dx
      CALL remove_particles

      ! Shift fields around
      CALL shift_field(ex)
      CALL shift_field(ey)
      CALL shift_field(ez)

      CALL shift_field(jx)
      CALL shift_field(jy)
      CALL shift_field(jz)

      CALL shift_field(bx)
      CALL shift_field(by)
      CALL shift_field(bz)
    ENDDO

  END SUBROUTINE shift_window



  SUBROUTINE shift_field(field)

    REAL(num), DIMENSION(-2:nx+3,-2:ny+3), INTENT(INOUT) :: field

    field(-2:nx+2,:) = field(-1:nx+3,:)
    CALL field_bc(field)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ipart, iy, i, isuby
    REAL(num) :: rand
    INTEGER :: clock, idum
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_y
    REAL(num), DIMENSION(-1:1) :: gy
    REAL(num) :: temp_local
#ifdef PER_PARTICLE_WEIGHT
    REAL(num) :: weight_local
#endif

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (coordinates(2) .EQ. nprocx-1) THEN
      clock = 9084263
      IF (use_random_seed) CALL SYSTEM_CLOCK(clock)
      idum = -(clock + rank)

      DO ispecies = 1, n_species
        DO iy = 1, ny
          DO ipart = 1, particle_species(ispecies)%npart_per_cell
            ALLOCATE(current)
            rand = random(idum) - 0.5_num
            current%part_pos(1) = x_max + dx + rand * dx
            rand = random(idum) - 0.5_num
            current%part_pos(2) = y(iy) + rand * dx

            cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
            cell_y = FLOOR(cell_y_r + 0.5_num)
            cell_frac_y = REAL(cell_y, num) - cell_y_r
            cell_y = cell_y + 1

            ! This uses the triangle shape function for particle weighting
            ! regardless of what is used by the rest of the code.
            ! All we are doing is injecting random noise so it shouldn't
            ! matter what the weighting is.
            gy(-1) = 0.5_num * (0.5_num + cell_frac_y)**2
            gy( 0) = 0.75_num - cell_frac_y**2
            gy( 1) = 0.5_num * (0.5_num - cell_frac_y)**2

            DO i = 1, 3
              temp_local = 0.0_num
              DO isuby = -1, +1
                temp_local = temp_local + gy(isuby) &
                    * particle_species(ispecies)%temperature(cell_y+isuby, i)
              ENDDO
              current%part_p(i) = &
                  momentum_from_temperature(particle_species(ispecies)%mass, &
                  temp_local, 0.0_num, idum)
            ENDDO

#ifdef PER_PARTICLE_WEIGHT
            weight_local = 0.0_num
            DO isuby = -1, +1
              weight_local = weight_local + gy(isuby) &
                  * particle_species(ispecies)%density(cell_y+isuby) &
                  / (REAL(particle_species(ispecies)%npart_per_cell, num) &
                  / (dx*dy))
            ENDDO
            current%weight = weight_local
#endif
#ifdef PARTICLE_DEBUG
            current%processor = rank
            current%processor_at_t0 = rank
#endif
            CALL add_particle_to_partlist(&
                particle_species(ispecies)%attached_list, current)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE insert_particles



  SUBROUTINE remove_particles

    TYPE(particle), POINTER :: current, next
    INTEGER :: ispecies

    IF (coordinates(2) .EQ. 0) THEN
      DO ispecies = 1, n_species
        current=>particle_species(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          next=>current%next
          IF (current%part_pos(1) .LT. x_min-0.5_num*dx) THEN
            CALL remove_particle_from_partlist(&
                particle_species(ispecies)%attached_list, current)
            DEALLOCATE(current)
          ENDIF
          current=>next
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE remove_particles

END MODULE window
