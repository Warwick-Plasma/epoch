! Copyright (C) 2009-2019 University of Warwick
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

    INTEGER :: ispecies, ix, iy
    INTEGER :: cell_x, cell_y
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: local_count
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    DO ispecies = 1, n_species
      local_count = species_list(ispecies)%attached_list%count
      CALL MPI_ALLREDUCE(local_count, species_list(ispecies)%global_count, &
          1, MPI_INTEGER8, MPI_SUM, comm, errcode)
      ALLOCATE(species_list(ispecies)%secondary_list(i0:nx+i1,i0:ny+i1))
      DO iy = i0, ny + i1
        DO ix = i0, nx + i1
          CALL create_empty_partlist(&
              species_list(ispecies)%secondary_list(ix,iy))
        END DO
      END DO
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next => current%next
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx) + 1
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy) + 1
#else
        cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
        cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
#endif

        CALL remove_particle_from_partlist(&
            species_list(ispecies)%attached_list, current)
        CALL add_particle_to_partlist(&
            species_list(ispecies)%secondary_list(cell_x,cell_y), current)
        current => next
      END DO
    END DO

  END SUBROUTINE reorder_particles_to_grid



  SUBROUTINE reattach_particles_to_mainlist

    INTEGER :: ispecies, ix, iy
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    DO ispecies = 1, n_species
      DO iy = i0, ny + i1
        DO ix = i0, nx + i1
          CALL append_partlist(species_list(ispecies)%attached_list, &
              species_list(ispecies)%secondary_list(ix,iy))
        END DO
      END DO
      DEALLOCATE(species_list(ispecies)%secondary_list)
    END DO

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
      END IF
    END DO

    use_particle_lists = use_particle_lists .OR. use_split
    IF (use_split) need_random_state = .TRUE.
#endif

  END SUBROUTINE setup_split_particles



  SUBROUTINE split_particles

#ifndef PER_SPECIES_WEIGHT
    INTEGER :: ispecies, ix, iy
    INTEGER(i8) :: count
    TYPE(particle), POINTER :: current, new_particle
    TYPE(particle_list) :: append_list
    REAL(num) :: jitter_x, jitter_y
    INTEGER :: i0, i1

    i0 = 1 - ng
    IF (use_field_ionisation) i0 = -ng
    i1 = 1 - i0

    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%split) CYCLE
      IF (species_list(ispecies)%npart_max > 0 &
          .AND. species_list(ispecies)%global_count &
          >= species_list(ispecies)%npart_max) CYCLE

      CALL create_empty_partlist(append_list)

      DO iy = i0, ny + i1
        DO ix = i0, nx + i1
          count = species_list(ispecies)%secondary_list(ix,iy)%count
          IF (count > 0 .AND. count <= npart_per_cell_min) THEN
            current => species_list(ispecies)%secondary_list(ix,iy)%head
            DO WHILE(ASSOCIATED(current) .AND. count <= npart_per_cell_min &
                .AND. current%weight >= 1.0_num)
              count = &
                  species_list(ispecies)%secondary_list(ix,iy)%count
              jitter_x = (2 * random() - 1) * 0.25_num * dx
              jitter_y = (2 * random() - 1) * 0.25_num * dy
              current%weight = 0.5_num * current%weight
              ALLOCATE(new_particle)
              new_particle = current
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
              new_particle%id = 0
#endif
              new_particle%part_pos(1) = current%part_pos(1) + jitter_x
              new_particle%part_pos(2) = current%part_pos(2) + jitter_y
              CALL add_particle_to_partlist(append_list, new_particle)
#ifdef PARTICLE_DEBUG
              ! If running with particle debugging, specify that this
              ! particle has been split
              new_particle%processor_at_t0 = -1
#endif
              NULLIFY(new_particle)

              current%part_pos(1) = current%part_pos(1) - jitter_x
              current%part_pos(2) = current%part_pos(2) - jitter_y
              current => current%next
            END DO
          END IF
        END DO
      END DO

      CALL append_partlist(species_list(ispecies)%attached_list, append_list)
    END DO
#endif

  END SUBROUTINE split_particles

END MODULE split_particle
