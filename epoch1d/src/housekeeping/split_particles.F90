MODULE split_particle

  USE shared_data
  USE partlist
  USE helper
  USE boundary

  IMPLICIT NONE

  SAVE

  INTEGER(KIND=8) :: npart_per_cell_min = 5

#ifdef SPLIT_PARTICLES_AFTER_PUSH
CONTAINS

  SUBROUTINE reorder_particles_to_grid

    INTEGER :: ispecies, cell_x, cell_y
    TYPE(particle), POINTER :: current, next
    INTEGER(KIND=8) :: local_count

    DO ispecies = 1, n_species
      local_count = particle_species(ispecies)%attached_list%count
      CALL MPI_ALLREDUCE(local_count, particle_species(ispecies)%global_count, 1, mpireal, MPI_SUM, comm, errcode)
      ALLOCATE(particle_species(ispecies)%secondary_list(0:nx+1, 0:ny+1))
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          CALL create_empty_partlist(particle_species(ispecies)%secondary_list(ix, iy))
        ENDDO
      ENDDO
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next=>current%next
        cell_x = INT((current%part_pos(1)-x_start_local)/dx)! +1
        cell_y = INT((current%part_pos(2)-y_start_local)/dy)! +1
        CALL remove_particle_from_partlist(particle_species(ispecies)%attached_list, current)
        CALL add_particle_to_partlist(particle_species(ispecies)%secondary_list(cell_x, cell_y), current)
        current=>next
      ENDDO
    ENDDO

  END SUBROUTINE reorder_particles_to_grid



  SUBROUTINE reattach_particles_to_mainlist

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          CALL append_partlist(particle_species(ispecies)%attached_list, particle_species(ispecies)%secondary_list(ix, iy))
        ENDDO
      ENDDO
      DEALLOCATE(particle_species(ispecies)%secondary_list)
    ENDDO

    CALL particle_bcs

  END SUBROUTINE reattach_particles_to_mainlist



  SUBROUTINE split_particles

    INTEGER :: ispecies
    INTEGER(KIND=8) :: count
    TYPE(particle), POINTER :: current, NEW
    INTEGER :: clock, idum
    REAL(num) :: jitter_x, jitter_y

    ! Reseed random number generator
    CALL SYSTEM_CLOCK(clock)
    idum = -(clock+rank)

    DO ispecies = 1, n_species
      IF (.NOT. particle_species(ispecies)%split) CYCLE
      IF (particle_species(ispecies)%npart_max .GT. 0 .AND. particle_species(ispecies)%global_count .GE. particle_species(ispecies)%npart_max) CYCLE
      DO iy = 0, ny+1
        DO ix = 0, nx+1
          count = particle_species(ispecies)%secondary_list(ix, iy)%count
          IF (count .GT. 0 .AND. count .LE. npart_per_cell_min) THEN
            current=>particle_species(ispecies)%secondary_list(ix, iy)%head
            DO WHILE(ASSOCIATED(current) .AND. count .LE. npart_per_cell_min .AND. current%weight .GE. 1.0_num)
              count = particle_species(ispecies)%secondary_list(ix, iy)%count
              jitter_x = random(idum)*dx/2.0_num - dx/4.0_num
              jitter_y = random(idum)*dy/2.0_num - dy/4.0_num
              current%weight = current%weight/2.0_num
              ALLOCATE(NEW)
              NEW = current
              new%part_pos(1) = current%part_pos(1)+jitter_x
              new%part_pos(2) = current%part_pos(2)+jitter_y
              CALL add_particle_to_partlist(particle_species(ispecies)%attached_list, NEW)
#ifdef PART_DEBUG
              ! If running with particle debugging, specify that this particle has been split
              new%processor_at_t0 = -1
#endif
              NULLIFY(NEW)

              current%part_pos(1) = current%part_pos(1)-jitter_x
              current%part_pos(2) = current%part_pos(2)-jitter_y
              current=>current%next
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE split_particles
#endif

END MODULE split_particle
