MODULE window

  USE helper
  USE boundary
  USE partlist

  IMPLICIT NONE

CONTAINS

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
      xb_global = xb_global+dx
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

    REAL(num), DIMENSION(-2:nx+3), INTENT(INOUT) :: field

    field(-2:nx+2) = field(-1:nx+3)
    CALL field_bc(field)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ipart, i
    REAL(num) :: rand
    REAL(num) :: temp_local
#ifdef PER_PARTICLE_WEIGHT
    REAL(num) :: weight_local
#endif

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (coordinates(1) .EQ. nprocx-1) THEN
      DO ispecies = 1, n_species
        DO ipart = 1, particle_species(ispecies)%npart_per_cell
          ALLOCATE(current)
          rand = random() - 0.5_num
          current%part_pos = x_max + dx + rand * dx

          DO i = 1, 3
            temp_local = particle_species(ispecies)%temperature(i)
            current%part_p(i) = &
                momentum_from_temperature(particle_species(ispecies)%mass, &
                temp_local, 0.0_num)
          ENDDO

#ifdef PER_PARTICLE_WEIGHT
          weight_local = particle_species(ispecies)%density &
                / (REAL(particle_species(ispecies)%npart_per_cell, num)/(dx))
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
    ENDIF

  END SUBROUTINE insert_particles



  SUBROUTINE remove_particles

    TYPE(particle), POINTER :: current, next
    INTEGER :: ispecies

    IF (coordinates(1) .EQ. 0) THEN
      DO ispecies = 1, n_species
        current=>particle_species(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          next=>current%next
          IF (current%part_pos .LT. x_min-0.5_num*dx) THEN
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
