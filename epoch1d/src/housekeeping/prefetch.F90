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

MODULE prefetch

  USE shared_data

CONTAINS

  ! Uses Intel 12 prefetch call to move particle data into cache.
  SUBROUTINE prefetch_particle(p)

    TYPE(particle), INTENT(INOUT) :: p

#ifdef PREFETCH
#ifndef PER_SPECIES_WEIGHT
    CALL mm_prefetch(p%part_p(1)) !, 1)
    CALL mm_prefetch(p%weight) !, 1)
#else
ERROR - '-DPREFETCH' must not be used in conjunction with '-DPER_SPECIES_WEIGHT'
#endif
#endif

  END SUBROUTINE prefetch_particle

END MODULE prefetch
