MODULE split_particle

  USE boundary

  IMPLICIT NONE

  SAVE

  INTEGER(i8) :: npart_per_cell_min = 5
  LOGICAL :: use_split = .FALSE.

CONTAINS

  SUBROUTINE reorder_particles_to_grid

    INTEGER :: ispecies, ix
    INTEGER :: cell_x
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: local_count

    DO ispecies = 1, n_species
      local_count = species_list(ispecies)%attached_list%count
      CALL MPI_ALLREDUCE(local_count, species_list(ispecies)%global_count, &
          1, MPI_INTEGER8, MPI_SUM, comm, errcode)
      ALLOCATE(species_list(ispecies)%secondary_list(nx))
      DO ix = 1, nx
        CALL create_empty_partlist(&
            species_list(ispecies)%secondary_list(ix))
      ENDDO
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
        cell_x = FLOOR((current%part_pos - x_grid_min_local) / dx + 1.5_num)

        CALL remove_particle_from_partlist(&
            species_list(ispecies)%attached_list, current)
        CALL add_particle_to_partlist(&
            species_list(ispecies)%secondary_list(cell_x), current)
        current => next
      ENDDO
    ENDDO

  END SUBROUTINE reorder_particles_to_grid



  SUBROUTINE reattach_particles_to_mainlist

    INTEGER :: ispecies, ix

    DO ispecies = 1, n_species
      DO ix = 1, nx
        CALL append_partlist(species_list(ispecies)%attached_list, &
            species_list(ispecies)%secondary_list(ix))
      ENDDO
      DEALLOCATE(species_list(ispecies)%secondary_list)
    ENDDO

    CALL particle_bcs

  END SUBROUTINE reattach_particles_to_mainlist



  SUBROUTINE setup_split_particles

#ifdef PER_PARTICLE_WEIGHT
    INTEGER :: ispecies

    use_split = .FALSE.
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%split) THEN
        use_split = .TRUE.
        EXIT
      ENDIF
    ENDDO

    use_particle_lists = use_particle_lists .OR. use_split
    IF (use_split) need_random_state = .TRUE.
#endif

  END SUBROUTINE setup_split_particles



  SUBROUTINE split_particles

#ifdef PER_PARTICLE_WEIGHT
    INTEGER :: ispecies, ix
    INTEGER(i8) :: count
    TYPE(particle), POINTER :: current, new_particle
    TYPE(particle_list) :: append_list
    REAL(num) :: jitter_x

    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%split) CYCLE
      IF (species_list(ispecies)%npart_max > 0 &
          .AND. species_list(ispecies)%global_count &
          >= species_list(ispecies)%npart_max) CYCLE

      CALL create_empty_partlist(append_list)

      DO ix = 1, nx
        count = species_list(ispecies)%secondary_list(ix)%count
        IF (count > 0 .AND. count <= npart_per_cell_min) THEN
          current => species_list(ispecies)%secondary_list(ix)%head
          DO WHILE(ASSOCIATED(current) .AND. count <= npart_per_cell_min &
              .AND. current%weight >= 1.0_num)
            count = species_list(ispecies)%secondary_list(ix)%count
            jitter_x = (2 * random() - 1) * 0.25_num * dx
            current%weight = 0.5_num * current%weight
            ALLOCATE(new_particle)
            new_particle = current
#if PARTICLE_ID || PARTICLE_ID4
            new_particle%id = 0
#endif
            new_particle%part_pos = current%part_pos + jitter_x
            CALL add_particle_to_partlist(append_list, new_particle)
#ifdef PARTICLE_DEBUG
            ! If running with particle debugging, specify that this
            ! particle has been split
            new_particle%processor_at_t0 = -1
#endif
            NULLIFY(new_particle)

            current%part_pos = current%part_pos - jitter_x
            current => current%next
          ENDDO
        ENDIF
      ENDDO

      CALL append_partlist(species_list(ispecies)%attached_list, append_list)
    ENDDO
#endif

  END SUBROUTINE split_particles

END MODULE split_particle
