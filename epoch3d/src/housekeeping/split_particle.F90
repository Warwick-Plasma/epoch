MODULE split_particle

  USE mpi
  USE boundary
  USE partlist
  USE random_generator

  IMPLICIT NONE

  SAVE

  INTEGER(KIND=8) :: npart_per_cell_min = 5
  LOGICAL :: split_none = .TRUE.

CONTAINS

  SUBROUTINE reorder_particles_to_grid

    INTEGER :: ispecies, ix, iy, iz
    INTEGER :: cell_x, cell_y, cell_z
    TYPE(particle), POINTER :: current, next
    INTEGER(KIND=8) :: local_count

    DO ispecies = 1, n_species
      local_count = species_list(ispecies)%attached_list%count
      CALL MPI_ALLREDUCE(local_count, species_list(ispecies)%global_count, &
          1, mpireal, MPI_SUM, comm, errcode)
      ALLOCATE(species_list(ispecies)%secondary_list(nx,ny,nz))
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            CALL create_empty_partlist(&
                species_list(ispecies)%secondary_list(ix,iy,iz))
          ENDDO
        ENDDO
      ENDDO
      current=>species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next=>current%next
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos(1) - x_min_local) / dx + 1.0_num)
        cell_y = FLOOR((current%part_pos(2) - y_min_local) / dy + 1.0_num)
        cell_z = FLOOR((current%part_pos(3) - z_min_local) / dz + 1.0_num)
#else
        cell_x = FLOOR((current%part_pos(1) - x_min_local) / dx + 1.5_num)
        cell_y = FLOOR((current%part_pos(2) - y_min_local) / dy + 1.5_num)
        cell_z = FLOOR((current%part_pos(3) - z_min_local) / dz + 1.5_num)
#endif
        CALL remove_particle_from_partlist(&
            species_list(ispecies)%attached_list, current)
        CALL add_particle_to_partlist(&
            species_list(ispecies)%secondary_list(cell_x,cell_y,cell_z), &
            current)
        current=>next
      ENDDO
    ENDDO

  END SUBROUTINE reorder_particles_to_grid



  SUBROUTINE reattach_particles_to_mainlist

    INTEGER :: ispecies, ix, iy, iz

    DO ispecies = 1, n_species
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            CALL append_partlist(species_list(ispecies)%attached_list, &
                species_list(ispecies)%secondary_list(ix,iy,iz))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(species_list(ispecies)%secondary_list)
    ENDDO

    CALL particle_bcs

  END SUBROUTINE reattach_particles_to_mainlist



  SUBROUTINE setup_split_particles

#ifdef PER_PARTICLE_WEIGHT
    INTEGER :: ispecies

    split_none = .TRUE.
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%split) THEN
        split_none = .FALSE.
        EXIT
      ENDIF
    ENDDO

    IF (.NOT.split_none) use_particle_lists = .TRUE.
#endif

  END SUBROUTINE setup_split_particles



  SUBROUTINE split_particles

#ifdef PER_PARTICLE_WEIGHT
    INTEGER :: ispecies, ix, iy, iz
    INTEGER(KIND=8) :: count
    TYPE(particle), POINTER :: current, new_particle
    TYPE(particle_list) :: append_list
    REAL(num) :: jitter_x, jitter_y, jitter_z

    IF (split_none) RETURN

    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%split) CYCLE
      IF (species_list(ispecies)%npart_max .GT. 0 &
          .AND. species_list(ispecies)%global_count &
          .GE. species_list(ispecies)%npart_max) CYCLE

      CALL create_empty_partlist(append_list)

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            count = species_list(ispecies)%secondary_list(ix,iy,iz)%count
            IF (count .GT. 0 .AND. count .LE. npart_per_cell_min) THEN
              current=>species_list(ispecies)%secondary_list(ix,iy,iz)%head
              DO WHILE(ASSOCIATED(current) .AND. count .LE. npart_per_cell_min &
                  .AND. current%weight .GE. 1.0_num)
                count = &
                    species_list(ispecies)%secondary_list(ix,iy,iz)%count
                jitter_x = (2 * random() - 1) * 0.25_num * dx
                jitter_y = (2 * random() - 1) * 0.25_num * dy
                jitter_z = (2 * random() - 1) * 0.25_num * dz
                current%weight = 0.5_num * current%weight
                ALLOCATE(new_particle)
                new_particle = current
                new_particle%part_pos(1) = current%part_pos(1) + jitter_x
                new_particle%part_pos(2) = current%part_pos(2) + jitter_y
                new_particle%part_pos(3) = current%part_pos(3) + jitter_z
                CALL add_particle_to_partlist(append_list, new_particle)
#ifdef PARTICLE_DEBUG
                ! If running with particle debugging, specify that this
                ! particle has been split
                new_particle%processor_at_t0 = -1
#endif
                NULLIFY(new_particle)

                current%part_pos(1) = current%part_pos(1) - jitter_x
                current%part_pos(2) = current%part_pos(2) - jitter_y
                current%part_pos(3) = current%part_pos(3) - jitter_z
                current=>current%next
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL append_partlist(species_list(ispecies)%attached_list, append_list)
    ENDDO
#endif

  END SUBROUTINE split_particles

END MODULE split_particle
