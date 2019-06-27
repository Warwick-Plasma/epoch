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

! Module to move particles between species based on energy
! written by M. G. Ramsay and N. J. Sircombe

MODULE particle_migration

  USE calc_df
  USE partlist
#ifdef PREFETCH
  USE prefetch
#endif

  IMPLICIT NONE

  PUBLIC :: initialise_migration, migrate_particles

CONTAINS

  SUBROUTINE migrate_particles(step)

    INTEGER, INTENT(IN) :: step
    INTEGER :: ispecies

    ! Is it time to migrate particles?
    IF (MOD(step, particle_migration_interval) /= 0) RETURN

    ! Update fluid energies & densities
    io_list = species_list
    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%migrate%fluid) CYCLE
      CALL update_fluid_energy(ispecies)
    END DO

    species_list(:)%migrate%done = .FALSE.

    DO ispecies = 1, n_species
      ! Migration must start with the most energetic end of the chain and work
      ! back
      CALL migration_chain(ispecies)
    END DO

  END SUBROUTINE migrate_particles



  RECURSIVE SUBROUTINE migration_chain(current_list)

    INTEGER, INTENT(IN) :: current_list
    INTEGER :: next_list

    ! Does this species need to be done?
    IF (.NOT. species_list(current_list)%migrate%this_species) RETURN
    IF (species_list(current_list)%migrate%done) RETURN

    ! If there's another species above this one in the chain, then do that one
    ! first
    IF (species_list(current_list)%migrate%promoteable) THEN
      next_list = species_list(current_list)%migrate%promote_to_species
      CALL migration_chain(next_list)
    END IF

    ! Do promotions and demotions on this species
    IF (species_list(current_list)%migrate%promoteable) &
        CALL promote_particles(current_list)
    IF (species_list(current_list)%migrate%demoteable) &
        CALL demote_particles(current_list)

    ! Flag this species as done
    species_list(current_list)%migrate%done = .TRUE.

  END SUBROUTINE migration_chain



  SUBROUTINE update_fluid_energy(fluid_list)

    INTEGER, INTENT(IN) :: fluid_list
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: tmp
    REAL(num), PARAMETER :: mean_steps = 0.25_num

    ALLOCATE(tmp(1-ng:nx+ng, 1-ng:ny+ng))

    CALL calc_temperature(tmp, fluid_list)

    species_list(fluid_list)%migrate%fluid_energy = mean_steps * tmp &
        + (1.0_num - mean_steps) &
        * species_list(fluid_list)%migrate%fluid_energy

    CALL calc_number_density(tmp, fluid_list)

    species_list(fluid_list)%migrate%fluid_density = mean_steps * tmp &
        + (1.0_num - mean_steps) &
        * species_list(fluid_list)%migrate%fluid_density

    DEALLOCATE(tmp)

  END SUBROUTINE update_fluid_energy



  SUBROUTINE promote_particles(from_list)

    INTEGER, INTENT(IN) :: from_list
    INTEGER :: to_list
    REAL(num) :: ke_multiplier, density_condition
    REAL(num) :: rsqrt_part_m
    REAL(num) :: part_ke, local_te, local_ne
    INTEGER(KIND=i8) :: ipart
    INTEGER :: ix, iy
    TYPE(particle), POINTER :: current, next
#include "particle_head.inc"

    to_list = species_list(from_list)%migrate%promote_to_species
    ke_multiplier = species_list(from_list)%migrate%promotion_energy_factor
    density_condition = species_list(from_list)%migrate%promotion_density
#ifndef PER_PARTICLE_CHARGE_MASS
    rsqrt_part_m = 1.0_num / SQRT(species_list(from_list)%mass)
#endif

    current => species_list(from_list)%attached_list%head
    DO ipart = 1, species_list(from_list)%attached_list%count
      next => current%next
#ifdef PREFETCH
      CALL prefetch_particle(next)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
      rsqrt_part_m = 1.0_num / SQRT(current%mass)
#endif
      part_ke = SUM((current%part_p(:) * rsqrt_part_m)**2)

#include "particle_to_grid.inc"

      local_te = 0.0_num
      local_ne = 0.0_num
      DO iy = sf_min, sf_max
      DO ix = sf_min, sf_max
        local_te = local_te + gx(ix) * gy(iy) &
            * species_list(from_list)%migrate%fluid_energy(cell_x+ix, &
            cell_y+iy)
        local_ne = local_ne + gx(ix) * gy(iy) &
            * species_list(from_list)%migrate%fluid_density(cell_x+ix, &
            cell_y+iy)
      END DO
      END DO

      IF (part_ke > ke_multiplier * 3.0_num * kb * local_te &
          .AND. local_ne < density_condition) &
          CALL swap_lists(from_list, to_list, current)

      current => next
    END DO

  END SUBROUTINE promote_particles



  SUBROUTINE demote_particles(from_list)

    INTEGER, INTENT(IN) :: from_list
    INTEGER :: to_list
    REAL(num) :: ke_multiplier, density_condition
    REAL(num) :: rsqrt_part_m
    REAL(num) :: part_ke, local_te, local_ne
    INTEGER(KIND=i8) :: ipart
    INTEGER :: ix, iy
    TYPE(particle), POINTER :: current, next
#include "particle_head.inc"

    to_list = species_list(from_list)%migrate%demote_to_species
    ke_multiplier = species_list(from_list)%migrate%demotion_energy_factor
    density_condition = species_list(from_list)%migrate%demotion_density
#ifndef PER_PARTICLE_CHARGE_MASS
    rsqrt_part_m = 1.0_num / SQRT(species_list(from_list)%mass)
#endif

    current => species_list(from_list)%attached_list%head
    DO ipart = 1, species_list(from_list)%attached_list%count
      next => current%next
#ifdef PREFETCH
      CALL prefetch_particle(next)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
      rsqrt_part_m = 1.0_num / SQRT(current%mass)
#endif
      part_ke = SUM((current%part_p(:) * rsqrt_part_m)**2)

#include "particle_to_grid.inc"

      local_te = 0.0_num
      local_ne = 0.0_num
      DO iy = sf_min, sf_max
      DO ix = sf_min, sf_max
        local_te = local_te + gx(ix) * gy(iy) &
            * species_list(to_list)%migrate%fluid_energy(cell_x+ix, &
            cell_y+iy)
        local_ne = local_ne + gx(ix) * gy(iy) &
            * species_list(to_list)%migrate%fluid_density(cell_x+ix, &
            cell_y+iy)
      END DO
      END DO

      IF (part_ke < ke_multiplier * 3.0_num * kb * local_te &
          .AND. local_ne >= density_condition) &
          CALL swap_lists(from_list, to_list, current)

      current => next
    END DO

  END SUBROUTINE demote_particles



  SUBROUTINE swap_lists(from_list, to_list, current)

    INTEGER, INTENT(IN) :: from_list, to_list
    TYPE(particle), POINTER :: current

    CALL remove_particle_from_partlist(species_list(from_list)%attached_list, &
        current)
    CALL add_particle_to_partlist(species_list(to_list)%attached_list, &
        current)

  END SUBROUTINE swap_lists



  SUBROUTINE initialise_migration

    INTEGER :: ispecies, ipromote, idemote, current_list, next_list

    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%migrate%this_species) CYCLE

      ipromote = species_list(ispecies)%migrate%promote_to_species
      idemote = species_list(ispecies)%migrate%demote_to_species

      IF ((ipromote >= 1) .AND. (ipromote <= n_species)) &
          species_list(ispecies)%migrate%promoteable = .TRUE.
      IF ((idemote >= 1) .AND. (idemote <= n_species)) &
          species_list(ispecies)%migrate%demoteable = .TRUE.

      IF ((.NOT. species_list(ispecies)%migrate%promoteable) &
          .AND. (.NOT. species_list(ispecies)%migrate%demoteable)) THEN
        species_list(ispecies)%migrate%this_species = .FALSE.
        IF (rank == 0) THEN
          WRITE(stdout, *) '*** WARNING ***'
          WRITE(stdout, *) 'No valid promotion or demotion species specified.'
          WRITE(stdout, *) 'Migration turned off for species ',  &
              TRIM(species_list(ispecies)%name)
        END IF
        CYCLE
      END IF

#ifndef PER_PARTICLE_CHARGE_MASS
      IF (species_list(ispecies)%migrate%promoteable) THEN
        IF (ABS(species_list(ispecies)%mass &
            - species_list(ipromote)%mass) > c_tiny &
            .OR. ABS(species_list(ispecies)%charge &
            - species_list(ipromote)%charge) > c_tiny) THEN
          species_list(ispecies)%migrate%promoteable = .FALSE.
          IF (rank == 0) THEN
            WRITE(stdout, *) '*** WARNING ***'
            WRITE(stdout, *) 'Attempting to promote between species with'
            WRITE(stdout, *) 'different charge and/or mass.'
            WRITE(stdout, *) 'Promotion turned off for species ',  &
                TRIM(species_list(ispecies)%name)
          END IF
        END IF
      END IF

      IF (species_list(ispecies)%migrate%demoteable) THEN
        IF (ABS(species_list(ispecies)%mass &
            - species_list(idemote)%mass) > c_tiny &
            .OR. ABS(species_list(ispecies)%charge &
            - species_list(idemote)%charge) > c_tiny) THEN
          species_list(ispecies)%migrate%demoteable = .FALSE.
          IF (rank == 0) THEN
            WRITE(stdout, *) '*** WARNING ***'
            WRITE(stdout, *) 'Attempting to demote between species with'
            WRITE(stdout, *) 'different charge and/or mass.'
            WRITE(stdout, *) 'Demotion turned off for species ',  &
                TRIM(species_list(ispecies)%name)
          END IF
        END IF
      END IF

      IF ((.NOT. species_list(ispecies)%migrate%promoteable) &
          .AND. (.NOT. species_list(ispecies)%migrate%demoteable)) THEN
        species_list(ispecies)%migrate%this_species = .FALSE.
        CYCLE
      END IF
#endif

      ! Fluids are 'background' species - species that can be promoted from
      ! and/or demoted to.
      IF (species_list(ispecies)%migrate%promoteable) &
          species_list(ispecies)%migrate%fluid = .TRUE.
      IF (species_list(ispecies)%migrate%demoteable) &
          species_list(idemote)%migrate%fluid = .TRUE.
    END DO

    ! Check for looped migration chains
    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%migrate%this_species) CYCLE

      current_list = ispecies
      DO WHILE(species_list(current_list)%migrate%promoteable)
        next_list = species_list(current_list)%migrate%promote_to_species
        IF (next_list == ispecies) THEN
          IF (rank == 0) THEN
            WRITE(stdout, *) '*** WARNING ***'
            WRITE(stdout, *) 'Looped promotion chain.'
            WRITE(stdout, *) 'Promotion from species ',  &
                TRIM(species_list(current_list)%name), ' to ', &
                TRIM(species_list(next_list)%name), ' disabled.'
          END IF
          species_list(current_list)%migrate%promote_to_species = 0
          species_list(current_list)%migrate%promoteable = .FALSE.
        ELSE
          current_list = next_list
        END IF
      END DO

      current_list = ispecies
      DO WHILE(species_list(current_list)%migrate%demoteable)
        next_list = species_list(current_list)%migrate%demote_to_species
        IF (next_list == ispecies) THEN
          IF (rank == 0) THEN
            WRITE(stdout, *) '*** WARNING ***'
            WRITE(stdout, *) 'Looped demotion chain.'
            WRITE(stdout, *) 'Demotion from species ',  &
                TRIM(species_list(current_list)%name), ' to ', &
                TRIM(species_list(next_list)%name), ' disabled.'
          END IF
          species_list(current_list)%migrate%demote_to_species = 0
          species_list(current_list)%migrate%demoteable = .FALSE.
        ELSE
          current_list = next_list
        END IF
      END DO

      IF ((.NOT. species_list(ispecies)%migrate%promoteable) &
          .AND. (.NOT. species_list(ispecies)%migrate%demoteable)) THEN
        species_list(ispecies)%migrate%this_species = .FALSE.
        IF (rank == 0) THEN
          WRITE(stdout, *) '*** WARNING ***'
          WRITE(stdout, *) 'No valid promotion or demotion species specified.'
          WRITE(stdout, *) 'Migration turned off for species ', &
              TRIM(species_list(ispecies)%name)
        END IF
      END IF
    END DO

    ! After all that, do we still need migration at all?
    IF (.NOT. ANY(species_list(:)%migrate%this_species)) THEN
      use_particle_migration = .FALSE.
      IF (rank == 0) THEN
        WRITE(stdout, *) '*** WARNING ***'
        WRITE(stdout, *) 'All particle migration has been disabled.'
      END IF

      RETURN
    END IF

    ! Need to keep time-averaged temperature and density fields for fluids.
    io_list = species_list
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%migrate%fluid) THEN
        ALLOCATE(species_list(ispecies)%migrate%fluid_energy(1-ng:nx+ng, &
            1-ng:ny+ng))
        CALL calc_temperature(species_list(ispecies)%migrate%fluid_energy, &
            ispecies)
        ALLOCATE(species_list(ispecies)%migrate%fluid_density(1-ng:nx+ng, &
            1-ng:ny+ng))
        CALL calc_number_density(species_list(ispecies)%migrate%fluid_density,&
            ispecies)
      END IF
    END DO

  END SUBROUTINE initialise_migration

END MODULE particle_migration
