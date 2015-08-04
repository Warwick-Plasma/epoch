! Copyright (C) 2010-2014 Keith Bennett <K.Bennett@warwick.ac.uk>
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

MODULE ic_module

  USE shared_data
  USE helper
  USE particles

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: manual_load

CONTAINS

#ifdef DELTAF_METHOD
  ! Find local particle temperature: at the moment, just copied and
  ! pasted from particle_temperature.F90

  SUBROUTINE params_local(current, temperature, drift, temp_local, drift_local)

    TYPE(particle), POINTER, INTENT(IN) :: current
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: temperature, drift
    REAL(num), INTENT(INOUT) :: temp_local, drift_local
    INTEGER :: ix, iy, iz

#include "particle_head.inc"

    ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

    temp_local = 0.0_num
    drift_local = 0.0_num
    DO iz = sf_min, sf_max
      DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          temp_local = temp_local + gx(ix) * gy(iy) * gz(iz) &
              * temperature(cell_x+ix, cell_y+iy, cell_z+iz)
          drift_local = drift_local + gx(ix) * gy(iy) * gz(iz) &
              * drift(cell_x+ix, cell_y+iy, cell_z+iz)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE params_local
#endif



  SUBROUTINE manual_load

#ifdef DELTAF_METHOD
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz
    REAL(num) :: f0_exponent, distribution, mass, npart_per_cell
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    TYPE(particle_species), POINTER :: species
    INTEGER :: ipart, ispecies
#if DELTAF_DEBUG
    REAL(num) :: part_weight, weight_back, f0_back
#endif

    ! f0 calculation: mainly, we need to calculate the phase space volumes.
    ! Calculate this based on the loading parameters. Easy to check
    ! that this is OK for a Maxwellian load by setting f0 = f0_back,
    ! and making sure the weights cancel.

    DO ispecies = 1, n_species
      species => species_list(ispecies)
      partlist => species%attached_list
      current => partlist%head
      ipart = 0
      DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
        mass = current%mass
#else
        mass = species%mass
#endif
        CALL params_local(current, initial_conditions(ispecies)%temp(:,:,:,1), &
            initial_conditions(ispecies)%drift(:,:,:,1), Tx, driftx)
        CALL params_local(current, initial_conditions(ispecies)%temp(:,:,:,2), &
            initial_conditions(ispecies)%drift(:,:,:,2), Ty, drifty)
        CALL params_local(current, initial_conditions(ispecies)%temp(:,:,:,3), &
            initial_conditions(ispecies)%drift(:,:,:,3), Tz, driftz)

        f0_exponent = ((current%part_p(1) - driftx)**2 / Tx &
                     + (current%part_p(2) - drifty)**2 / Ty &
                     + (current%part_p(3) - driftz)**2 / Tz) / (2 * kb * mass)

        npart_per_cell = current%pvol

        ! We want to calculate the distribution of markers.
        distribution = EXP(-f0_exponent) * (npart_per_cell / (dx * dy * dz)) &
            / SQRT((2 * pi * kb * mass)**3 * Tx * Ty * Tz)
        current%pvol = 1.0_num / distribution

#if DELTAF_DEBUG
        f0_back = f0(ispecies, mass, current%part_p(1), current%part_p(2), &
            current%part_p(3))

        ! Checks for correct particle weight calculation.
        weight_back = f0_back * current%pvol
#ifdef PER_SPECIES_WEIGHT
        part_weight = species_list(ispecies)%weight
#else
        part_weight = current%weight
#endif
        WRITE(*,*) ipart, distribution, f0_exponent, npart_per_cell, &
            SQRT((2 * pi * kb * mass)**3 * Tx * Ty * Tz), kb, mass
        WRITE(*,*) ipart, 'R', EXP(-f0_exponent), EXP(-f0_exponent) &
            * npart_per_cell / SQRT((2 * pi * kb * mass)**3 * Tx * Ty * Tz)
        WRITE(*,*) ipart, 'Q', distribution, f0_back, weight_back, part_weight
#endif

        current => current%next
        ipart = ipart + 1
      ENDDO
    ENDDO
#endif

  END SUBROUTINE manual_load

END MODULE ic_module
