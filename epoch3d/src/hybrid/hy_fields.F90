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

  REAL(num), PRIVATE :: half_dt_by_dx, half_dt_by_dy, half_dt_by_dz
  REAL(num), PRIVATE :: i_mu_dx, i_mu_dy, i_mu_dz

CONTAINS

  SUBROUTINE setup_hybrid_fields

    ! Pre-calculate repeated variables

    half_dt_by_dx = 0.5_num * dt / dx
    half_dt_by_dy = 0.5_num * dt / dy
    half_dt_by_dz = 0.5_num * dt / dz

    i_mu_dx = 1.0_num / (mu0 * dx)
    i_mu_dy = 1.0_num / (mu0 * dy)
    i_mu_dz = 1.0_num / (mu0 * dz)

  END SUBROUTINE setup_hybrid_fields



  SUBROUTINE run_hybrid_fields

    ! This is the subroutine called from the hybrid PIC loop. We update B using
    ! the current E, sort out the boundary conditions, and recalculate E using
    ! the new B

    CALL half_b_step
    CALL hybrid_b_bc
    CALL calculate_e
    CALL hybrid_e_bc

  END SUBROUTINE run_hybrid_fields



  SUBROUTINE half_b_step

    ! This subroutine performs a half-step in the magnetic field, assuming a
    ! constant electric field and using:
    !
    ! dB/dt = -curl(E)
    !
    ! We calculate into 1 ghost cell to allow non-zero curls across simulation
    ! boundaries

    INTEGER :: ix, iy, iz
    REAL(num) :: dt_dey_dx, dt_dez_dx, dt_dex_dy, dt_dez_dy
    REAL(num) :: dt_dex_dz, dt_dey_dz

    ! Update B by half a timestep
    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1

          ! Calculate gradients
          dt_dey_dx = (ey(ix+1,iy,iz) - ey(ix,iy,iz)) * half_dt_by_dx
          dt_dez_dx = (ez(ix+1,iy,iz) - ez(ix,iy,iz)) * half_dt_by_dx

          dt_dex_dy = (ex(ix,iy+1,iz) - ex(ix,iy,iz)) * half_dt_by_dy
          dt_dez_dy = (ez(ix,iy+1,iz) - ez(ix,iy,iz)) * half_dt_by_dy

          dt_dex_dz = (ex(ix,iy,iz+1) - ex(ix,iy,iz)) * half_dt_by_dz
          dt_dey_dz = (ey(ix,iy,iz+1) - ey(ix,iy,iz)) * half_dt_by_dz

          ! Update fields
          bx(ix, iy, iz) = bx(ix, iy, iz) - dt_dez_dy + dt_dey_dz
          by(ix, iy, iz) = by(ix, iy, iz) - dt_dex_dz + dt_dez_dx
          bz(ix, iy, iz) = bz(ix, iy, iz) - dt_dey_dx + dt_dex_dy
        END DO
      END DO
    END DO

  END SUBROUTINE half_b_step



  SUBROUTINE hybrid_b_bc

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

  END SUBROUTINE hybrid_b_bc



  SUBROUTINE calculate_e

    ! Calculates the electric field for the current values of global variables
    ! bx, by, bz, j, and resistivity, using the equation:
    !
    ! E = resistivity * (curl(B)/mu_0 - J)
    !
    ! We have precalculated i_mu_dx = 1/(mu_0*dx)

    ! Note that B is staggered from E, whereas J is evaulated at the same point
    ! as E. Resistivity is a cell centred variable. We calculate into 1 ghost
    ! cell to allow non-zero curls across simulation boundaries

    INTEGER :: ix, iy, iz
    REAL(num) :: dby_dx, dbz_dx, dbx_dy, dbz_dy, dbx_dz, dby_dz
    REAL(num) :: res_cell, res_x, res_y, res_z

    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1

          ! Calculate derivatives
          dby_dx = (by(ix,iy,iz) - by(ix-1,iy,iz)) * i_mu_dx
          dbz_dx = (bz(ix,iy,iz) - bz(ix-1,iy,iz)) * i_mu_dx
          dbx_dy = (bx(ix,iy,iz) - bx(ix,iy-1,iz)) * i_mu_dy
          dbz_dy = (bz(ix,iy,iz) - bz(ix,iy-1,iz)) * i_mu_dy
          dbx_dz = (bx(ix,iy,iz) - bx(ix,iy,iz-1)) * i_mu_dz
          dby_dz = (by(ix,iy,iz) - by(ix,iy,iz-1)) * i_mu_dz

          ! Calculate average resistivity at the electric field position
          res_cell = resistivity(ix,iy,iz)
          res_x = 0.5_num * (resistivity(ix+1,iy,iz) + res_cell)
          res_y = 0.5_num * (resistivity(ix,iy+1,iz) + res_cell)
          res_z = 0.5_num * (resistivity(ix,iy,iz+1) + res_cell)

          ! Calculate the return currents
          jbx(ix,iy,iz) = dbz_dy - dby_dz - jx(ix,iy,iz)
          jby(ix,iy,iz) = dbx_dz - dbz_dx - jy(ix,iy,iz)
          jbz(ix,iy,iz) = dby_dx - dbx_dy - jz(ix,iy,iz)

          ! Update fields
          ex(ix,iy,iz) = res_x * jbx(ix,iy,iz)
          ey(ix,iy,iz) = res_y * jby(ix,iy,iz)
          ez(ix,iy,iz) = res_z * jbz(ix,iy,iz)
        END DO
      END DO
    END DO

  END SUBROUTINE calculate_e



  SUBROUTINE hybrid_e_bc

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

  END SUBROUTINE hybrid_e_bc



  SUBROUTINE field_zero_curl(field, boundary)

    ! If we are on the simulation boundary, let the ghost cells match those of
    ! the first ghost cell (e.g. by(nx+2:nx+ng,1,1) = by(nx+1,1,1)). This forces
    ! the gradient across the first two ghost cells to be zero

    INTEGER, INTENT(IN) :: boundary
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, nn

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      DO i = 1, ng-1
        field(-i,:,:) = field(0,:,:)
      END DO
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      DO i = 1, ng
        field(nx+i,:,:) = field(nx,:,:)
      END DO

    ELSE IF (boundary == c_bd_y_min .AND. y_min_boundary) THEN
      DO i = 1, ng-1
        field(:,-i,:) = field(:,0,:)
      END DO
    ELSE IF (boundary == c_bd_y_max .AND. y_max_boundary) THEN
      DO i = 1, ng
        field(:,ny+i,:) = field(:,ny,:)
      END DO

    ELSE IF (boundary == c_bd_z_min .AND. z_min_boundary) THEN
      DO i = 1, ng-1
        field(:,:,-i) = field(:,:,0)
      END DO
    ELSE IF (boundary == c_bd_z_max .AND. z_max_boundary) THEN
      DO i = 1, ng
        field(:,:,nz+i) = field(:,:,nz)
      END DO
    END IF

  END SUBROUTINE field_zero_curl

#endif
END MODULE hy_fields
