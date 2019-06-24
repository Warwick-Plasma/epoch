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

MODULE deltaf_loader

  USE shared_data

  IMPLICIT NONE

CONTAINS

#ifdef DELTAF_METHOD
  ! Find local particle temperature: at the moment, just copied and
  ! pasted from particle_temperature.F90

  SUBROUTINE params_local(current, temperature, drift, temp_local, drift_local)

    TYPE(particle), POINTER, INTENT(IN) :: current
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: temperature, drift
    REAL(num), INTENT(INOUT) :: temp_local, drift_local
    REAL(num) :: gf
    INTEGER :: ix, iy

#include "particle_head.inc"

    ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

    temp_local = 0.0_num
    drift_local = 0.0_num
    DO iy = sf_min, sf_max
    DO ix = sf_min, sf_max
      gf = gx(ix) * gy(iy)
      temp_local = temp_local + gf * temperature(cell_x+ix, cell_y+iy)
      drift_local = drift_local + gf * drift(cell_x+ix, cell_y+iy)
    END DO
    END DO

  END SUBROUTINE params_local
#endif



  SUBROUTINE deltaf_load(ispecies, species_temp, species_drift)

    REAL(num), DIMENSION(:,:,:), POINTER :: species_temp, species_drift
    INTEGER, INTENT(IN) :: ispecies
#ifdef DELTAF_METHOD
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz
    REAL(num) :: f0_exponent, distribution, mass, npart_per_cell, idx
    REAL(num) :: two_kb_mass, two_pi_kb_mass3, part_weight
    REAL(num), PARAMETER :: two_kb = 2.0_num * kb
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    TYPE(particle_species), POINTER :: species
    INTEGER :: ipart
#if DELTAF_DEBUG
    REAL(num) :: weight_back, f0_back
#endif

    ! f0 calculation: mainly, we need to calculate the phase space volumes.
    ! Calculate this based on the loading parameters. Easy to check
    ! that this is OK for a Maxwellian load by setting f0 = f0_back,
    ! and making sure the weights cancel.

    idx = 1.0_num / dx / dy

    species => species_list(ispecies)
    partlist => species%attached_list
    current => partlist%head

    mass = species%mass
    two_kb_mass = two_kb * mass
    two_pi_kb_mass3 = (pi * two_kb_mass)**3
    part_weight = species_list(ispecies)%weight

    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
      two_kb_mass = two_kb * mass
      two_pi_kb_mass3 = (pi * two_kb_mass)**3
#endif
      CALL params_local(current, species_temp(:,:,1), &
          species_drift(:,:,1), Tx, driftx)
      CALL params_local(current, species_temp(:,:,2), &
          species_drift(:,:,2), Ty, drifty)
      CALL params_local(current, species_temp(:,:,3), &
          species_drift(:,:,3), Tz, driftz)

      f0_exponent = ((current%part_p(1) - driftx)**2 / Tx &
                   + (current%part_p(2) - drifty)**2 / Ty &
                   + (current%part_p(3) - driftz)**2 / Tz) / two_kb_mass

      npart_per_cell = current%pvol

      ! We want to calculate the distribution of markers.
      distribution = EXP(-f0_exponent) * npart_per_cell * idx &
          / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
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
#endif

  END SUBROUTINE deltaf_load

END MODULE deltaf_loader
