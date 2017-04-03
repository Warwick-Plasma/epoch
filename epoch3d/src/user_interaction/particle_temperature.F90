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

  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_species, &
      drift)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_species), POINTER :: part_species
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: ix, iy, iz
#include "particle_head.inc"

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

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

      IF (direction == c_dir_x) current%part_p(1) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction == c_dir_y) current%part_p(2) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction == c_dir_z) current%part_p(3) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      current => current%next
      ipart = ipart + 1
    ENDDO

  END SUBROUTINE setup_particle_temperature



  FUNCTION momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: momentum_from_temperature

    REAL(num) :: stdev
    REAL(num) :: rand1, rand2, w
    REAL(num), SAVE :: val
    LOGICAL, SAVE :: cached = .FALSE.

    ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers
    ! The standard deviation (stdev) is related to temperature

    stdev = SQRT(temperature * kb * mass)

    IF (cached) THEN
      cached = .FALSE.
      momentum_from_temperature = val * stdev + drift
    ELSE
      cached = .TRUE.

      DO
        rand1 = random()
        rand2 = random()

        rand1 = 2.0_num * rand1 - 1.0_num
        rand2 = 2.0_num * rand2 - 1.0_num

        w = rand1**2 + rand2**2

        IF (w > c_tiny .AND. w < 1.0_num) EXIT
      ENDDO

      w = SQRT((-2.0_num * LOG(w)) / w)

      momentum_from_temperature = rand1 * w * stdev + drift
      val = rand2 * w
    ENDIF

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
