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

MODULE shape_functions

  USE constants

  IMPLICIT NONE

CONTAINS

  SUBROUTINE particle_to_grid(cell_frac, gf)

    REAL(num), INTENT(IN) :: cell_frac
    REAL(num), DIMENSION(sf_min:sf_max), INTENT(OUT) :: gf
#ifndef PARTICLE_SHAPE_TOPHAT
    REAL(num) :: cf2
#endif
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third

    cf2 = cell_frac**2
    gf(-2) = fac1 * (0.5_num + cell_frac)**4
    gf(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac &
        + cf2 * (1.5_num - cell_frac - cf2))
    gf( 0) = 0.25_num * (fac3 + cf2 * (cf2 - 2.5_num))
    gf( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac &
        + cf2 * (1.5_num + cell_frac - cf2))
    gf( 2) = fac1 * (0.5_num - cell_frac)**4
#elif  PARTICLE_SHAPE_TOPHAT
    gf( 0) = 0.5_num + cell_frac
    gf( 1) = 0.5_num - cell_frac
#else
    cf2 = cell_frac**2
    gf(-1) = 0.5_num * (0.25_num + cf2 + cell_frac)
    gf( 0) = 0.75_num - cf2
    gf( 1) = 0.5_num * (0.25_num + cf2 - cell_frac)
#endif

  END SUBROUTINE particle_to_grid



  SUBROUTINE grid_to_particle(cell_frac, gf)

    REAL(num), INTENT(IN) :: cell_frac
    REAL(num), DIMENSION(sf_min:sf_max), INTENT(OUT) :: gf

    CALL particle_to_grid(cell_frac, gf)

  END SUBROUTINE grid_to_particle

END MODULE shape_functions
