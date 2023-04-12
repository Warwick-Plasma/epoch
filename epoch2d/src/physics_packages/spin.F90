! Copyright (C) 2009-2023 University of Warwick
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

MODULE spin

#ifdef PARTICLE_SPIN

  USE constants
  USE shared_data
  USE random_generator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_particle_spin(species, new_particle)
    TYPE(particle_species), INTENT(IN) :: species
    TYPE(particle), POINTER :: new_particle

    SELECT CASE(species%spin_distribution)
      CASE(c_spin_uniform)
        CALL init_particle_spin_uniform(species, new_particle)
      CASE(c_spin_directed)
        CALL init_particle_spin_directed(species, new_particle)
    END SELECT

  END SUBROUTINE init_particle_spin


  SUBROUTINE init_particle_spin_directed(species, new_particle)
    TYPE(particle_species), INTENT(IN) :: species
    TYPE(particle), POINTER :: new_particle

    new_particle%spin = species%spin_orientation
  END SUBROUTINE init_particle_spin_directed

  SUBROUTINE init_particle_spin_uniform(species, new_particle)
    TYPE(particle_species), INTENT(IN) :: species
    TYPE(particle), POINTER :: new_particle

    REAL(num) :: theta, phi
    REAL(num) :: sx, sy, sz

    phi = 2.0_num * pi * random()
    theta = acos(1.0_num - 2.0_num * random())
    sx = sin(theta) * cos(phi)
    sy = sin(theta) * sin(phi)
    sz = cos(theta)

    new_particle%spin = (/ sx, sy, sz /)
  END SUBROUTINE init_particle_spin_uniform
  
#endif

END MODULE spin