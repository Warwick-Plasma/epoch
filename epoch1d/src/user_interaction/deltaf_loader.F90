! Copyright (C) 2016      Ben McMillan <B.F.McMillan@warwick.ac.uk>
! Copyright (C) 2016      Keith Bennett <K.Bennett@warwick.ac.uk>
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

MODULE deltaf_loader

  USE shared_data
  USE particles, ONLY: f0
  USE random_generator
  USE particle_temperature

  IMPLICIT NONE

  REAL(num), PARAMETER :: threshold = 0.5_num
  REAL(num), PARAMETER :: temp_fac = 4.0_num

CONTAINS

#ifdef DELTAF_METHOD
  ! Find local particle temperature: at the moment, just copied and
  ! pasted from particle_temperature.F90

  SUBROUTINE params_local(current, temperature, drift, temp_local, drift_local)

    TYPE(particle), POINTER, INTENT(IN) :: current
    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: temperature, drift
    REAL(num), INTENT(INOUT) :: temp_local, drift_local
    REAL(num) :: gf
    INTEGER :: ix, i

#include "particle_head.inc"

    ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

    temp_local = 0.0_num
    drift_local = 0.0_num
    DO ix = sf_min, sf_max
      i = cell_x + ix
      gf = gx(ix)
      temp_local = temp_local + gf * temperature(i)
      drift_local = drift_local + gf * drift(i)
    END DO

  END SUBROUTINE params_local



  SUBROUTINE density_local(current, density, dens_local)

    TYPE(particle), POINTER, INTENT(IN) :: current
    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: density
    REAL(num), INTENT(INOUT) :: dens_local
    INTEGER :: ix, i

#include "particle_head.inc"

#include "particle_to_grid.inc"

    dens_local = 0.0_num
    DO ix = sf_min, sf_max
      i = cell_x + ix

! Lines below commented out: should save density_map.
! Can we assume it is always 1 for deltaf cases?

!#ifdef PARTICLE_SHAPE_TOPHAT
!      IF (.NOT. density_map(i)) i = cell_x + 1 - ix
!#else
!      IF (.NOT. density_map(i)) THEN
!          i = cell_x + ix / 2
!#ifdef PARTICLE_SHAPE_BSPLINE3
!        IF (.NOT. density_map(i)) i = cell_x - ix / 2
!#endif
!      END IF
!#endif

      dens_local = dens_local + gx(ix) * density(i)
    END DO

  END SUBROUTINE density_local
#endif



  ! The initial particle distribution function f according to initial conditions.
  ! (temp_fac!=1 allows a scaled version for the tail-resolving method)

  SUBROUTINE f_dist(current, species, distribution, density, temp_fac)

    TYPE(particle_species), POINTER :: species
    TYPE(particle), POINTER :: current
    REAL(num), INTENT(OUT) :: distribution, density
    REAL(num), INTENT(IN)  :: temp_fac
    REAL(num), PARAMETER :: two_kb = 2.0_num * kb
    REAL(num) :: mass, two_kb_mass, two_pi_kb_mass
    REAL(num) :: f0_exponent, normalisation_term
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz

#ifdef PER_PARTICLE_CHARGE_MASS
    mass = current%mass
#else
    mass = species%mass
#endif
    two_kb_mass = two_kb * mass
    two_pi_kb_mass = pi * two_kb_mass
    CALL params_local(current, species%initial_conditions%temp(:,1), &
        species%initial_conditions%drift(:,1), Tx, driftx)
    CALL params_local(current, species%initial_conditions%temp(:,2), &
        species%initial_conditions%drift(:,2), Ty, drifty)
    CALL params_local(current, species%initial_conditions%temp(:,3), &
        species%initial_conditions%drift(:,3), Tz, driftz)
    CALL density_local(current, species%initial_conditions%density(:), &
        density)

    Tx = Tx * temp_fac
    Ty = Ty * temp_fac
    Tz = Tz * temp_fac
    ! To allow temperatures to be zero in y or z direction,
    f0_exponent = 0
    normalisation_term = 1
    IF (Tx /= 0) THEN
      f0_exponent = f0_exponent + (current%part_p(1) - driftx)**2 / Tx
      normalisation_term = normalisation_term * two_pi_kb_mass * Tx
    END IF
    IF (Ty /= 0) THEN
      f0_exponent = f0_exponent + (current%part_p(2) - drifty)**2 / Ty
      normalisation_term = normalisation_term * two_pi_kb_mass * Ty
    END IF
    IF (Tz /= 0) THEN
      f0_exponent = f0_exponent + (current%part_p(3) - driftz)**2 / Tz
      normalisation_term = normalisation_term * two_pi_kb_mass * Tz
    END IF
    f0_exponent = f0_exponent / two_kb_mass

    ! We want to calculate the distribution of markers.
    distribution = density * EXP(-f0_exponent)  / SQRT(normalisation_term)

  END SUBROUTINE f_dist



  ! Load a proportion of markers with one temperature, rest with a second
  ! temperature.
  SUBROUTINE twotemp_loadv

    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: ispecies
    REAL(num) :: tfac
    TYPE(particle_species), POINTER :: species
    LOGICAL :: twot
    INTEGER :: iend

    DO ispecies = 1, n_species
      species => species_list(ispecies)

      iend = LEN(TRIM(species%name))
      twot = (species%name(iend-1:iend) == '2t')

      IF (ABS(species%initial_conditions%density_back) <= c_tiny .OR. twot) &
          CYCLE

      partlist => species%attached_list
      current => partlist%head
      ipart = 0
      DO WHILE(ipart < partlist%count)
        IF (random() < threshold) THEN
          tfac = temp_fac
        ELSE
          tfac = 1.0_num
        END IF

        CALL set_single_particle_temperature(&
            species_list(ispecies)%initial_conditions%temp(:,1), c_dir_x, &
            species, current, &
            species_list(ispecies)%initial_conditions%drift(:,1), tfac)
        CALL set_single_particle_temperature(&
            species_list(ispecies)%initial_conditions%temp(:,2), c_dir_y, &
            species, current, &
            species_list(ispecies)%initial_conditions%drift(:,2), tfac)
        CALL set_single_particle_temperature(&
            species_list(ispecies)%initial_conditions%temp(:,3), c_dir_z, &
            species, current, &
            species_list(ispecies)%initial_conditions%drift(:,3), tfac)

        current => current%next
        ipart = ipart + 1
      END DO
    END DO

  END SUBROUTINE twotemp_loadv



  SUBROUTINE set_single_particle_temperature(temperature, direction, &
      part_species, current, drift, tfac)

    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_species), POINTER :: part_species
    TYPE(particle), POINTER, INTENT(INOUT) :: current
    REAL(num), DIMENSION(1-ng:), INTENT(IN) :: drift
    REAL(num), INTENT(IN) :: tfac
    REAL(num) :: mass, temp_local, drift_local, gf
    INTEGER :: ix, i

#include "particle_head.inc"

#ifdef PER_PARTICLE_CHARGE_MASS
    mass = current%mass
#else
    mass = part_species%mass
#endif

    ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

    temp_local = 0.0_num
    drift_local = 0.0_num
    DO ix = sf_min, sf_max
      i = cell_x + ix
      gf = gx(ix)
      temp_local = temp_local + gf * temperature(i)
      drift_local = drift_local + gf * drift(i)
    END DO
    temp_local = temp_local * tfac

    IF (direction == c_dir_x) THEN
      current%part_p(1) = &
          momentum_from_temperature(mass, temp_local, drift_local)
    END IF

    IF (direction == c_dir_y) THEN
      current%part_p(2) = &
          momentum_from_temperature(mass, temp_local, drift_local)
    END IF

    IF (direction == c_dir_z) THEN
      current%part_p(3) = &
          momentum_from_temperature(mass, temp_local, drift_local)
    END IF

  END SUBROUTINE set_single_particle_temperature



  SUBROUTINE deltaf_load

#ifdef DELTAF_METHOD
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz
    REAL(num) :: f0_exponent, distribution, mass, npart_per_cell, idx
    REAL(num) :: two_kb_mass, two_pi_kb_mass, two_pi_kb_mass3
    REAL(num) :: part_weight, normalisation_term
    REAL(num), PARAMETER :: two_kb = 2.0_num * kb
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    TYPE(particle_species), POINTER :: species
    INTEGER :: ipart, ispecies
#if DELTAF_DEBUG
    REAL(num) :: weight_back, f0_back
#endif

    ! f0 calculation: mainly, we need to calculate the phase space volumes.
    ! Calculate this based on the loading parameters. Easy to check
    ! that this is OK for a Maxwellian load by setting f0 = f0_back,
    ! and making sure the weights cancel.

    idx = 1.0_num / dx

    DO ispecies = 1, n_species
      species => species_list(ispecies)
      partlist => species%attached_list
      current => partlist%head

      IF (ABS(species%initial_conditions%density_back) <= c_tiny) CYCLE

      mass = species%mass
      two_kb_mass = two_kb * mass
      two_pi_kb_mass  = pi * two_kb_mass
      two_pi_kb_mass3 = two_pi_kb_mass**3
      part_weight = species_list(ispecies)%weight

      ipart = 0
      DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
        mass = current%mass
#else
        mass = species%mass
#endif
        two_kb_mass = two_kb * mass
        two_pi_kb_mass  = pi * two_kb_mass
        two_pi_kb_mass3 = two_pi_kb_mass**3

        CALL params_local(current, species%initial_conditions%temp(:,1), &
            species%initial_conditions%drift(:,1), Tx, driftx)
        CALL params_local(current, species%initial_conditions%temp(:,2), &
            species%initial_conditions%drift(:,2), Ty, drifty)
        CALL params_local(current, species%initial_conditions%temp(:,3), &
            species%initial_conditions%drift(:,3), Tz, driftz)

        ! To allow temperatures to be zero (for cases where deltaf not used or
        ! ignorable dir)
        f0_exponent = 0
        normalisation_term = 1
        IF (Tx /= 0) THEN
          f0_exponent = f0_exponent + (current%part_p(1) - driftx)**2 / Tx
          normalisation_term = normalisation_term * two_pi_kb_mass * Tx
        END IF
        IF (Ty /= 0) THEN
          f0_exponent = f0_exponent + (current%part_p(2) - drifty)**2 / Ty
          normalisation_term = normalisation_term * two_pi_kb_mass * Ty
        END IF
        IF (Tz /= 0) THEN
          f0_exponent = f0_exponent + (current%part_p(3) - driftz)**2 / Tz
          normalisation_term = normalisation_term * two_pi_kb_mass * Tz
        END IF
        f0_exponent = f0_exponent / two_kb_mass

        npart_per_cell = current%pvol

        ! We want to calculate the distribution of markers.
        distribution = EXP(-f0_exponent) * npart_per_cell * idx &
            / SQRT(normalisation_term)

        current%pvol = 1.0_num / distribution

#if DELTAF_DEBUG
        f0_back = f0(ispecies, mass, current%part_p)

        ! Checks for correct particle weight calculation.
        weight_back = f0_back * current%pvol
#ifndef PER_SPECIES_WEIGHT
        part_weight = current%weight
#endif
        WRITE(*,*) ipart, distribution, f0_exponent, npart_per_cell, &
            SQRT(two_pi_kb_mass3 * Tx * Ty * Tz), kb, mass
        WRITE(*,*) ipart, 'R', EXP(-f0_exponent), EXP(-f0_exponent) &
            * npart_per_cell / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
        WRITE(*,*) ipart, 'Q', distribution, f0_back, weight_back, part_weight
#endif

        current => current%next
        ipart = ipart + 1
      END DO
    END DO
#endif

  END SUBROUTINE deltaf_load



  SUBROUTINE deltaf_load_twotemp

#if defined(DELTAF_METHOD) && !defined(PER_SPECIES_WEIGHT)
    REAL(num) :: idx, f0_1, f0_2, distribution, npart_per_cell
    REAL(num) :: mass, density
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    TYPE(particle_species), POINTER :: species
    INTEGER :: ipart, ispecies
    INTEGER :: iend
    LOGICAL :: twot

    ! f0 calculation:
    ! 1) calculate the phase space volumes based on the loading.
    ! 2) set the initial particle weights appropriately.
    ! Check
    ! that this is OK for a Maxwellian load by setting f0 = f0_back,
    ! and making sure the weights cancel.

    idx = 1.0_num / dx

    DO ispecies = 1, n_species
      species => species_list(ispecies)
      partlist => species%attached_list
      current => partlist%head

      iend = LEN(TRIM(species%name))
      twot = (species%name(iend-1:iend) == '2t')

      ! Only perform operation if the species uses delta_f method
      IF (ABS(species%initial_conditions%density_back) <= c_tiny) CYCLE

      IF (rank == 0) THEN
        IF (twot) THEN
          WRITE(*,*) 'Delta-f load (tail resolving) of species ', &
              TRIM(species%name)
        ELSE
          WRITE(*,*) 'Delta-f load of species ', TRIM(species%name)
        END IF
      END IF
      mass = species%mass

      ipart = 0
      DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
        mass = current%mass
#endif
        npart_per_cell = current%pvol

        ! Calculate the desired initial particle distribution (f0_1)
        ! and the marker distribution (distribution)
        IF (twot) THEN
           CALL f_dist(current, species, f0_2, density, temp_fac)
           CALL f_dist(current, species, f0_1, density, 1.0_num)
           distribution = (npart_per_cell / density) * idx &
               * (threshold * f0_1 + (1.0_num - threshold) * f0_2)
        ELSE
           CALL f_dist(current, species, f0_1, density, 1.0_num)
           distribution = (npart_per_cell / density) * idx * f0_1
        END IF

        ! Phase space volume per marker is just the inverse of the marker
        ! distribution function.
        current%pvol = 1.0_num / distribution
        ! And the particle weight is the particle distribution
        ! divided by the marker distribution. This is already
        ! set correctly for a normal Maxwellian load.
        current%weight = f0_1 * current%pvol

        current => current%next
        ipart = ipart + 1
      END DO
    END DO
#endif

  END SUBROUTINE deltaf_load_twotemp

END MODULE deltaf_loader
