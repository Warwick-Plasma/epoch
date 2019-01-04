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

MODULE particle_ring_beam

  USE particle_temperature
  USE shared_data
  USE random_generator

  IMPLICIT NONE

CONTAINS

  ! Subroutine to initialise a ring beam particle distribution
  ! Assumes linear interpolation of temperatures between cells
  ! Create a beam distribution along z, and ring distribution around x and y
  ! then rotate to point along magnetic field
  SUBROUTINE setup_particle_ring_beam(species, temperatures, drifts)
    TYPE(particle_species), POINTER :: species
    REAL(num), DIMENSION(1-ng:, :), INTENT(IN) :: temperatures
    REAL(num), DIMENSION(1-ng:, :), INTENT(IN) :: drifts
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, gyrophase, pperp, ppara, part_x
    REAL(num), DIMENSION(1:3) :: temps_local, drifts_local
    TYPE(particle), POINTER :: current
    REAL(num), PARAMETER :: golden_angle = pi * (3.0_num - SQRT(5.0_num))
    INTEGER(i8) :: ipart
    INTEGER :: ix, i
    REAL(num), DIMENSION(1:3, 1:3) :: rotation_matrix
    REAL(num), DIMENSION(1:3) :: aligned_momentum
#include "particle_head.inc"

    partlist => species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = species%mass
#endif

      ! Assume that temperatures is cell centred
#include "particle_to_grid.inc"

      temps_local = 0.0_num
      drifts_local = 0.0_num
      DO i = 1, 3
        DO ix = sf_min, sf_max
          temps_local(i) = temps_local(i) + gx(ix) * temperatures(cell_x+ix, i)
          drifts_local(i) = drifts_local(i) + gx(ix) * drifts(cell_x+ix, i)
        END DO
      END DO

      part_x  = current%part_pos - x_grid_min_local

      pperp = momentum_from_temperature(mass, temps_local(1), drifts_local(3))
      ppara = momentum_from_temperature(mass, temps_local(3), drifts_local(3))
      gyrophase = ipart * golden_angle

      CALL setup_rotation_matrix(current, rotation_matrix)

      aligned_momentum(1) = pperp * COS(gyrophase)
      aligned_momentum(2) = pperp * SIN(gyrophase)
      aligned_momentum(3) = ppara
      
      current%part_p(1) = DOT_PRODUCT(rotation_matrix(1, :), aligned_momentum)
      current%part_p(2) = DOT_PRODUCT(rotation_matrix(2, :), aligned_momentum)
      current%part_p(3) = DOT_PRODUCT(rotation_matrix(3, :), aligned_momentum)

      current => current%next
      ipart = ipart + 1
    END DO

  END SUBROUTINE setup_particle_ring_beam

  SUBROUTINE setup_rotation_matrix(part, r)

    Type(Particle), POINTER, INTENT(INOUT) :: part
    REAL(num), DIMENSION(1:3,1:3), INTENT(INOUT) :: r
    REAL(num) :: th, ux, uy, uz, ax, ay, az, fx, fy, fz, det, u2, part_x
#include "fields_at_particle_variable_declarations.inc"

    part_x = part%part_pos

#include "fields_at_particle_implementation.inc"

    fx = bx_part
    fy = by_part
    fz = bz_part
    ! normalise
    fx = fx / SQRT(fx**2 + fy**2 + fz**2)
    fy = fy / SQRT(fx**2 + fy**2 + fz**2)
    fz = fz / SQRT(fx**2 + fy**2 + fz**2)

    ax = 1.0_num
    ay = 0.0_num
    az = 0.0_num

    ux = ay * fz - az * fy
    uy = az * fx - ax * fz
    uz = ax * fy - ay * fx

    th = ACOS(ax * fx + ay * fy + az * fz)
    
    r = 0.0_num
    r(1,1) = 1.0_num
    r(2,2) = 1.0_num
    r(3,3) = 1.0_num
    IF ((ax*fx + ay*fy + az*fz) .EQ. -1.0_num) THEN
      r = -r
    ENDIF

    u2 = SQRT(ux**2 + uy**2 + uz**2)
    IF (u2 .GT. 0.0_num) THEN
      ux = ux/u2
      uy = uy/u2
      uz = uz/u2
      r(1,1) = ux**2 + (1 - ux**2) * COS(th)
      r(1,2) = ux * uy * (1 - COS(th)) - uz * SIN(th)
      r(1,3) = ux * uz * (1 - COS(th)) + uy * SIN(th)
      r(2,1) = ux * uy * (1 - COS(th)) + uz * SIN(th)
      r(2,2) = uy**2 + (1 - uy**2) * COS(th)
      r(2,3) = uy * uz * (1 - COS(th)) - ux * SIN(th)
      r(3,1) = ux * uz * (1 - COS(th)) - uy * SIN(th)
      r(3,2) = uy * uz * (1 - COS(th)) + ux * SIN(th)
      r(3,3) = uz**2 + (1 - uz**2) * COS(th)
    END IF

    det = r(1,1)*(r(2,2)*r(3,3) - r(2,3)*r(3,2)) &
        + r(1,2)*(r(2,3)*r(3,1) - r(2,1)*r(3,3)) &
        + r(1,3)*(r(2,1)*r(3,2) - r(2,2)*r(3,1))
    r = r / det

  END SUBROUTINE setup_rotation_matrix



END MODULE particle_ring_beam
