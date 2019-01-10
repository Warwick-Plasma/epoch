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

MODULE particle_temperature

  USE shared_data
  USE random_generator

  IMPLICIT NONE

CONTAINS


  ! Subroutine to initialise a species with a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(species, temperature, drift)
    TYPE(particle_species), POINTER :: species
    REAL(num), DIMENSION(1-ng:, :), INTENT(IN) :: temperature
    REAL(num), DIMENSION(1-ng:, :), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    partlist => species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)

      CALL setup_particle_temperature_particle(current, temperature, drift, &
          species%mass)

      current => current%next
      ipart = ipart + 1
    END DO

  END SUBROUTINE setup_particle_temperature



  ! Subroutine to initialise a particle with a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature_particle(current, temperature, &
      drift, species_mass)
    TYPE(particle), POINTER :: current
    REAL(num), DIMENSION(1-ng:, :), INTENT(IN) :: temperature
    REAL(num), DIMENSION(1-ng:, :), INTENT(IN) :: drift
    REAL(num), INTENT(IN) :: species_mass
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    INTEGER(i8) :: ipart
    INTEGER :: ix, i
#include "particle_head.inc"

#ifdef PER_PARTICLE_CHARGE_MASS
    mass = current%mass
#else
    mass = species_mass
#endif

    ! Assume that temperature is cell centred
#include "particle_to_grid.inc"
    DO i = 1, 3
      temp_local = 0.0_num
      drift_local = 0.0_num
      DO ix = sf_min, sf_max
        temp_local = temp_local + gx(ix) * temperature(cell_x+ix, i)
        drift_local = drift_local + gx(ix) * drift(cell_x+ix, i)
      END DO

      current%part_p(i) = &
          momentum_from_temperature(mass, temp_local, drift_local)
    END DO

  END SUBROUTINE setup_particle_temperature_particle



  FUNCTION momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: momentum_from_temperature
    DOUBLE PRECISION :: stdev, mu

    stdev = DBLE(SQRT(temperature * kb * mass))
    mu = DBLE(drift)
    momentum_from_temperature = random_box_muller(stdev, mu)

  END FUNCTION momentum_from_temperature


  ! Function for generating momenta of thermal particles in a particular
  ! direction, e.g. the +x direction.
  ! These satisfy a Rayleigh distribution, formed by combining two
  ! normally-distributed (~N(0,sigma)) random variables as follows:
  ! Z = SQRT(X**2 + Y**2)
  FUNCTION flux_momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: flux_momentum_from_temperature
    REAL(num) :: mom1, mom2

    mom1 = momentum_from_temperature(mass, temperature, 0.0_num)
    mom2 = momentum_from_temperature(mass, temperature, 0.0_num)

    flux_momentum_from_temperature = SQRT(mom1**2 + mom2**2) + drift

  END FUNCTION flux_momentum_from_temperature


END MODULE particle_temperature
