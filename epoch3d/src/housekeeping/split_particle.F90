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

MODULE split_particle

  USE boundary

  IMPLICIT NONE

  SAVE

  INTEGER(i8) :: npart_per_cell_min = 5
  LOGICAL :: use_split = .FALSE.

CONTAINS

  SUBROUTINE reorder_particles_to_grid

    INTEGER :: ispecies, ix, iy, iz
    INTEGER :: cell_x, cell_y, cell_z
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: local_count

    DO ispecies = 1, n_species
      local_count = species_list(ispecies)%attached_list%count
      CALL MPI_ALLREDUCE(local_count, species_list(ispecies)%global_count, &
          1, MPI_INTEGER8, MPI_SUM, comm, errcode)
      ALLOCATE(species_list(ispecies)%secondary_list(nx,ny,nz))
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            CALL create_empty_partlist(&
                species_list(ispecies)%secondary_list(ix,iy,iz))
          ENDDO
        ENDDO
      ENDDO
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
        cell_z = FLOOR((current%part_pos(3) - z_grid_min_local) / dz + 1.5_num)

        CALL remove_particle_from_partlist(&
            species_list(ispecies)%attached_list, current)
        CALL add_particle_to_partlist(&
            species_list(ispecies)%secondary_list(cell_x,cell_y,cell_z), &
            current)
        current => next
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

#ifndef PER_SPECIES_WEIGHT
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

#ifndef PER_SPECIES_WEIGHT
    INTEGER :: ispecies, ix, iy, iz
    INTEGER(i8) :: count
    TYPE(particle), POINTER :: current, new_particle
    TYPE(particle_list) :: append_list
    REAL(num) :: jitter_x, jitter_y, jitter_z

    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%split) CYCLE
      IF (species_list(ispecies)%npart_max > 0 &
          .AND. species_list(ispecies)%global_count &
          >= species_list(ispecies)%npart_max) CYCLE

      CALL create_empty_partlist(append_list)

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            count = species_list(ispecies)%secondary_list(ix,iy,iz)%count
            IF (count > 0 .AND. count <= npart_per_cell_min) THEN
              current => species_list(ispecies)%secondary_list(ix,iy,iz)%head
              DO WHILE(ASSOCIATED(current) .AND. count <= npart_per_cell_min &
                  .AND. current%weight >= 1.0_num)
                count = &
                    species_list(ispecies)%secondary_list(ix,iy,iz)%count
                jitter_x = (2 * random() - 1) * 0.25_num * dx
                jitter_y = (2 * random() - 1) * 0.25_num * dy
                jitter_z = (2 * random() - 1) * 0.25_num * dz
                current%weight = 0.5_num * current%weight
                ALLOCATE(new_particle)
                new_particle = current
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
                new_particle%id = 0
#endif
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
                current => current%next
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
