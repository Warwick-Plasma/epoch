MODULE split_particle

  USE shared_data
  USE helper

  IMPLICIT NONE

  SAVE

  INTEGER(KIND=8) :: npart_per_cell_min = 5

#ifdef SPLIT_PARTICLES_AFTER_PUSH
CONTAINS

  SUBROUTINE reorder_particles_to_grid

    INTEGER :: ispecies, ix
    INTEGER :: cell_x
    TYPE(particle), POINTER :: current, next
    INTEGER(KIND=8) :: local_count

    DO ispecies = 1, n_species
      local_count = particle_species(ispecies)%attached_list%count
      CALL MPI_ALLREDUCE(local_count, particle_species(ispecies)%global_count, &
          1, mpireal, MPI_SUM, comm, errcode)
      ALLOCATE(particle_species(ispecies)%secondary_list(0:nx+1))
      CALL create_empty_partlist(&
          particle_species(ispecies)%secondary_list(ix))
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next=>current%next
        cell_x = INT((current%part_pos-x_min_local)/dx)! +1
        CALL remove_particle_from_partlist(&
            particle_species(ispecies)%attached_list, current)
        CALL add_particle_to_partlist(&
            particle_species(ispecies)%secondary_list(cell_x), current)
        current=>next
      ENDDO
    ENDDO

  END SUBROUTINE reorder_particles_to_grid



  SUBROUTINE reattach_particles_to_mainlist

    INTEGER :: ispecies, ix

    DO ispecies = 1, n_species
      DO ix = 0, nx+1
        CALL append_partlist(particle_species(ispecies)%attached_list, &
            particle_species(ispecies)%secondary_list(ix))
      ENDDO
      DEALLOCATE(particle_species(ispecies)%secondary_list)
    ENDDO

    CALL particle_bcs

  END SUBROUTINE reattach_particles_to_mainlist



  SUBROUTINE split_particles

    INTEGER :: ispecies, ix
    INTEGER(KIND=8) :: count
    TYPE(particle), POINTER :: current, new_particle
    REAL(num) :: jitter_x

    DO ispecies = 1, n_species
      IF (.NOT. particle_species(ispecies)%split) CYCLE
      IF (particle_species(ispecies)%npart_max .GT. 0 &
          .AND. particle_species(ispecies)%global_count &
          .GE. particle_species(ispecies)%npart_max) CYCLE

      DO ix = 0, nx+1
        count = particle_species(ispecies)%secondary_list(ix)%count
        IF (count .GT. 0 .AND. count .LE. npart_per_cell_min) THEN
          current=>particle_species(ispecies)%secondary_list(ix)%head
          DO WHILE(ASSOCIATED(current) .AND. count .LE. npart_per_cell_min &
              .AND. current%weight .GE. 1.0_num)
            count = particle_species(ispecies)%secondary_list(ix)%count
            jitter_x = random()*dx/2.0_num - dx/4.0_num
            current%weight = current%weight/2.0_num
            ALLOCATE(new_particle)
            new_particle = current
            new_particle%part_pos = current%part_pos+jitter_x
            CALL add_particle_to_partlist(&
                particle_species(ispecies)%attached_list, new_particle)
#ifdef PARTICLE_DEBUG
            ! If running with particle debugging, specify that this
            ! particle has been split
            new_particle%processor_at_t0 = -1
#endif
            NULLIFY(new_particle)

            current%part_pos = current%part_pos-jitter_x
            current=>current%next
          ENDDO
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE split_particles
#endif

END MODULE split_particle
