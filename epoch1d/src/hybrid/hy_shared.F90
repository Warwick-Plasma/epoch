! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
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
!
!-------------------------------------------------------------------------------
!
! hy_shared.F90
!
! This module contains subroutines which are general enough to not belong to any
! one hybrid module.

MODULE hy_shared
#ifdef HYBRID

  USE boundary
  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE hy_grid_centred_var_at_particle(part_x, part_var, grid_var)

  ! Calculates the value of a grid-centred variable "part_var" stored in the
  ! grid "grid_var", averaged over the particle shape for a particle at part_x

  REAL(num), INTENT(IN) :: part_x
  REAL(num), INTENT(IN) :: grid_var(1-ng:nx+ng)
  REAL(num), INTENT(OUT) :: part_var
  INTEGER :: cell_x1
  REAL(num) :: cell_x_r
  REAL(num) :: cell_frac_x
  REAL(num), DIMENSION(sf_min:sf_max) :: gx
#ifdef PARTICLE_SHAPE_BSPLINE3
  REAL(num) :: cf2
  REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
  REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
  REAL(num) :: cf2
  REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

  ! The following method is lifted from photons.F90 (field_at_particle), for
  ! the cell-centered fields, taking into account the various particle shapes
#ifdef PARTICLE_SHAPE_TOPHAT
  cell_x_r = part_x / dx - 0.5_num
#else
  cell_x_r = part_x / dx
#endif
  cell_x1 = FLOOR(cell_x_r + 0.5_num)
  cell_frac_x = REAL(cell_x1, num) - cell_x_r
  cell_x1 = cell_x1 + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
  part_var = &
        gx(-2) * grid_var(cell_x1-2) &
      + gx(-1) * grid_var(cell_x1-1) &
      + gx( 0) * grid_var(cell_x1  ) &
      + gx( 1) * grid_var(cell_x1+1) &
      + gx( 2) * grid_var(cell_x1+2)
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
  part_var = &
        gx(0) * grid_var(cell_x1  ) &
      + gx(1) * grid_var(cell_x1+1)
#else
#include "triangle/gx.inc"
  part_var = &
        gx(-1) * grid_var(cell_x1-1) &
      + gx( 0) * grid_var(cell_x1  ) &
      + gx( 1) * grid_var(cell_x1+1)
#endif
  part_var = fac*part_var

END SUBROUTINE hy_grid_centred_var_at_particle

#endif
END MODULE hy_shared
