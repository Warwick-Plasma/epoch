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
! #ifdef PARTICLE_SPIN    

  USE shared_data

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE init_particle_spin(species, new_particle)
    TYPE(particle_species), INTENT(INOUT) :: species
    TYPE(particle), INTENT(INOUT) :: new_particle

    new_particle%spin(1) = 2.0

  END SUBROUTINE init_particle_spin

! #endif
END MODULE spin