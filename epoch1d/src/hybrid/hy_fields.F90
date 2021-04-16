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
!------------------------------------------------------------------------------
!
! hy_fields.F90
!
! This module houses the hybrid field solver, and the special zero-curl boundary
! condition for treating fields beyond the simulation window

MODULE hy_fields
#ifdef HYBRID

  USE hy_shared

  IMPLICIT NONE

  REAL(num), PRIVATE :: half_dt_by_dx
  REAL(num), PRIVATE :: i_mu_dx

CONTAINS

  SUBROUTINE setup_hybrid_fields

    ! Pre-calculate repeated variables

    half_dt_by_dx = 0.5_num * dt / dx

    i_mu_dx = 1.0_num / (mu0 * dx)

  END SUBROUTINE setup_hybrid_fields



  SUBROUTINE run_hybrid_fields

    ! This is the subroutine called from the hybrid PIC loop. We update B using
    ! the current E, sort out the boundary conditions, and recalculate E using
    ! the new B

    CALL half_B_step
    CALL hybrid_B_bc
    CALL calculate_E
    CALL hybrid_E_bc

  END SUBROUTINE run_hybrid_fields



  SUBROUTINE half_B_step

    ! This subroutine performs a half-step in the magnetic field, assuming a
    ! constant electric field and using:
    !
    ! dB/dt = -curl(E)
    !
    ! We calculate into 1 ghost cell to allow non-zero curls across simulation
    ! boundaries

    INTEGER :: ix
    REAL(num) :: dt_dey_dx, dt_dez_dx

    ! Update B by half a timestep
    DO ix = 0, nx+1

      ! Calculate gradients
      dt_dey_dx = (ey(ix+1) - ey(ix)) * half_dt_by_dx
      dt_dez_dx = (ez(ix+1) - ez(ix)) * half_dt_by_dx

      ! Update fields
      by(ix) = by(ix) + dt_dez_dx
      bz(ix) = bz(ix) - dt_dey_dx
    END DO

  END SUBROUTINE half_B_step



  SUBROUTINE hybrid_B_bc

    ! This script updates the ghost cells on the local processor, by either
    ! passing over those from the neighbouring processor, or by ensuring the
    ! curls beyond the boundary are zero

    INTEGER :: i

    ! Pass neighbouring ghost cells
    CALL field_bc(bx, ng)
    CALL field_bc(by, ng)
    CALL field_bc(bz, ng)

    ! Special cases for boundary processors
    DO i = 1, 2*c_ndims
      CALL field_zero_curl(bx, i)
      CALL field_zero_curl(by, i)
      CALL field_zero_curl(bz, i)
    END DO

  END SUBROUTINE hybrid_B_bc



  SUBROUTINE calculate_E

    ! Calculates the electric field for the current values of global variables
    ! bx, by, bz, j, and resistivity, using the equation:
    !
    ! E = resistivity * (curl(B)/mu_0 - J)
    !
    ! We have precalculated i_mu_dx = 1/(mu_0*dx)

    ! Note that B is staggered from E, whereas J is evaulated at the same point
    ! as E. Resistivity is a cell centred variable. We calculate into 1 ghost
    ! cell to allow non-zero curls across simulation boundaries

    INTEGER :: ix
    REAL(num) :: dby_dx, dbz_dx
    REAL(num) :: res_cell, res_x

    DO ix = 0, nx+1

      ! Calculate derivatives
      dby_dx = (by(ix) - by(ix-1)) * i_mu_dx
      dbz_dx = (bz(ix) - bz(ix-1)) * i_mu_dx

      ! Calculate average resistivity at the electric field position
      res_cell = resistivity(ix)
      res_x = 0.5_num * (resistivity(ix+1) + res_cell)

      ! Calculate the return currents
      jbx(ix) = -jx(ix)
      jby(ix) = - dbz_dx - jy(ix)
      jbz(ix) = dby_dx - jz(ix)

      ! Update fields
      ex(ix) = res_x * jbx(ix)
      ey(ix) = res_cell * jby(ix)
      ez(ix) = res_cell * jbz(ix)
    END DO

  END SUBROUTINE calculate_E



  SUBROUTINE hybrid_E_bc

    ! This script updates the ghost cells on the local processor, by either
    ! passing over those from the neighbouring processor, or by ensuring the
    ! curls beyond the boundary are zero

    INTEGER :: i

    ! Pass neighbouring ghost cells
    CALL field_bc(ex, ng)
    CALL field_bc(ey, ng)
    CALL field_bc(ez, ng)

    ! Special cases for boundary processors
    DO i = 1, 2*c_ndims
      CALL field_zero_curl(ex, i)
      CALL field_zero_curl(ey, i)
      CALL field_zero_curl(ez, i)
    END DO

  END SUBROUTINE hybrid_E_bc



  SUBROUTINE field_zero_curl(field, boundary)

    ! If we are on the simulation boundary, let the ghost cells match those of
    ! the first ghost cell (e.g. by(nx+2:nx+ng) = by(nx+1)). This forces
    ! the gradient across the first two ghost cells to be zero

    INTEGER, INTENT(IN) :: boundary
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, nn

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      DO i = 1, ng-1
        field(-i) = field(0)
      END DO
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      DO i = 1, ng
        field(nx+i) = field(nx)
      END DO
    END IF

  END SUBROUTINE field_zero_curl

#endif
END MODULE hy_fields
