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
! hy_heating.F90
!
! This module contains the scripts related to changing the background
! temperatures.

MODULE hy_heating
#ifdef HYBRID

  USE hy_shared

  IMPLICIT NONE

  REAL(num), PRIVATE, ALLOCATABLE :: eff_heat_capacity(:,:,:)
  REAL(num), PRIVATE, ALLOCATABLE :: ion_heat_const(:,:,:)
  REAL(num), PRIVATE, ALLOCATABLE :: ohmic_heat_const(:,:,:)

  REAL(num), PRIVATE, PARAMETER :: c_hy_equil = &
      2.0_num/3.0_num*(2.0_num*pi*kb)**(-1.5_num) * q0**4 * SQRT(m0)/epsilon0**2

  REAL(num), PRIVATE :: idx, idy, idz

CONTAINS

  SUBROUTINE setup_heating

    ! Precalculate useful variables for speed

    idx = 1.0_num/dx
    idy = 1.0_num/dy
    idz = 1.0_num/dz

  END SUBROUTINE setup_heating



  SUBROUTINE get_heat_capacity

    ! Calculates heat capacity as described by Davies, et al, (2002). Phys. Rev.
    ! E, 65(2), 026407, for each element of solid_array
    !
    ! C = 0.3 + 1.2*T'*(2.2 + T')/(1.1 + T')^2
    !
    ! Where T' = Z^(-4/3) * Te * kB / q0, as Te in this code is measured in K
    !
    ! We also calculate the effective heat capacity Sum(ne/C) for compound
    ! targets, where Sum refers to the sum over all solids in the cell
    !
    ! Finally, this routine pre-calculates heating constants for different
    ! heating processes

    REAL(num), ALLOCATABLE :: T_prime(:,:,:), i_kb_sumne2(:,:,:)
    INTEGER :: ix, iy, iz, i_sol

    ! Temporary varible to store T' array values
    ALLOCATE(T_prime(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    ! Cycle through solids and calculate C values
    DO i_sol = 1, solid_count
      T_prime = solid_array(i_sol)%Z_prime * hy_Te
      ALLOCATE(solid_array(i_sol)%heat_capacity(1-ng:nx+ng,1-ng:ny+ng, &
          1-ng:nz+ng))
      solid_array(i_sol)%heat_capacity = 0.3_num + 1.2_num * T_prime &
          * (2.2_num + T_prime)/(1.1_num + T_prime)**2
    END DO

    ! Deallocate temporary variable T_prime
    DEALLOCATE(T_prime)

    ! Get effective heat capacity
    ALLOCATE(eff_heat_capacity(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    eff_heat_capacity = 0.0_num
    DO i_sol = 1, solid_count
      eff_heat_capacity = eff_heat_capacity &
          + solid_array(i_sol)%el_density / solid_array(i_sol)%heat_capacity
    END DO

    ! Get 1/kb/Sum(ne)² values
    ALLOCATE(i_kb_sumne2(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    i_kb_sumne2 = 0.0_num
    DO i_sol = 1, solid_count
      i_kb_sumne2 = i_kb_sumne2 + solid_array(i_sol)%el_density
    END DO
    i_kb_sumne2 = 1.0_num/(kb*i_kb_sumne2**2)

    ! Precalculate the ion_heat_const array if ionisation loss is running
    IF (use_hybrid_collisions) THEN
      ALLOCATE(ion_heat_const(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ion_heat_const = i_kb_sumne2 / (dx * dy * dz)
    END IF

    ! Precalculate the ohmic_heat_const array if Ohmic heating is running
    IF (use_ohmic) THEN
      ALLOCATE(ohmic_heat_const(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ohmic_heat_const = dt * i_kb_sumne2
    END IF

    DEALLOCATE(i_kb_sumne2)

  END SUBROUTINE get_heat_capacity



  SUBROUTINE ionisation_heating(dE, part_x, part_y, part_z)

    ! The current particle has deposited energy dE into the solid, at position
    ! (part_x, part_y, part_z), due to ionisation energy loss. By finding the
    ! heat capacity at that point, this is converted to a temperature increase.
    !
    ! The temperature rise of a compound solid will be approximated to:
    !
    ! dT = (Energy change per unit vol.) * Sum(ne/C) / [Sum(ne)]² / kB
    !
    ! ion_heat_const = 1 / kB / sum(ne²) / (cell vol.)
    ! eff_heat_capacity = sum(ne/C)

    REAL(num), INTENT(IN) :: dE, part_x, part_y, part_z
    INTEGER :: ix, iy, iz
    REAL(num) :: part_C_eff, part_heat_const, delta_Te

    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, &
        part_heat_const, ion_heat_const)
    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, &
        part_C_eff, eff_heat_capacity)

    ! Calculate the temperature increase, and add this to Te
    delta_Te = dE * part_C_eff * part_heat_const

    ! Write temperature change to the grid (to current cell only, ignores shape)
#ifdef PARTICLE_SHAPE_TOPHAT
    ix = FLOOR(part_x * idx) + 1
    iy = FLOOR(part_y * idy) + 1
    iz = FLOOR(part_z * idz) + 1
#else
    ix = FLOOR(part_x * idx + 0.5_num) + 1
    iy = FLOOR(part_y * idy + 0.5_num) + 1
    iz = FLOOR(part_z * idz + 0.5_num) + 1
#endif
    hy_Te(ix,iy,iz) = hy_Te(ix,iy,iz) + delta_Te

  END SUBROUTINE ionisation_heating



  SUBROUTINE ohmic_heating

    ! Calculates the Ohmic heating of the simulation grid, as described by
    ! Davies, et al, (2002). Phys. Rev. E, 65(2), 026407. At this point in the
    ! simulation, we are at the end of a timestep, with E evaluated in the
    ! middle of this timestep. Assuming this is the average electric field over
    ! the timestep, we have a power disspiation of I²R, which is equivalent to a
    ! power density of J².(resistivity)
    !
    !  dT = (Power density) * dt  * Sum(ne/C) / [Sum(ne)]² / kB
    !
    ! ohmic_heat_const = dt / kB / Sum(ne)²
    ! eff_heat_capacity = sum(ne/C)

    REAL(num) :: j2
    INTEGER :: ix, iy, iz, i_sol

    ! Loop over all grid points to find temperature change
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx

          ! Tb is a cell-centred variable, but J has stagger - need to average
          j2 = 0.25_num * ((jbx(ix,iy,iz) + jbx(ix-1,iy,iz))**2 &
              + (jby(ix,iy,iz) + jby(ix,iy-1,iz))**2 &
              + (jbz(ix,iy,iz) + jbz(ix,iy,iz-1))**2)

          ! Calculate Ohmic heating, avoiding the 0/0 NaN
          IF (eff_heat_capacity(ix,iy,iz) > 0.0_num) THEN
            hy_Te(ix,iy,iz) = hy_Te(ix,iy,iz) &
                + j2*resistivity(ix,iy,iz)*ohmic_heat_const(ix,iy,iz) &
                * eff_heat_capacity(ix,iy,iz)
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE ohmic_heating



  SUBROUTINE clear_heat_capacity

    ! De-allocate the heat capacity in each solid. This array is recalculated
    ! each time-step anyway, so there's no need to save and resize it in the
    ! load balancer. Same applies for the effective heat capacity

    INTEGER :: i_sol

    DO i_sol = 1, solid_count
      DEALLOCATE(solid_array(i_sol)%heat_capacity)
    END DO
    DEALLOCATE(eff_heat_capacity)
    IF (use_hybrid_collisions) DEALLOCATE(ion_heat_const)
    IF (use_ohmic) DEALLOCATE(ohmic_heat_const)

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hy_Te, ng)

  END SUBROUTINE clear_heat_capacity



  SUBROUTINE thermal_equilibration

    ! Calculates transfer of thermal energy between electrons and ions using the
    ! thermal equilibration rate, from Spitzer (1962) - "Physics of Fully
    ! Ionized Gases", Second Edition, number 3, Interscience Tracts on Physics &
    ! Astronomy. Available for free on
    ! https://openlibrary.org/books/OL5856372M/Physics_of_fully_ionized_gases
    ! (last accessed 05/July/2020)
    !
    ! The thermal rate of change is taken to be the ratio of e-ion temperature
    ! difference to equilibration time (eq. 5-30), where equilibration time is
    ! also given (eq. 5-31). Note that these equations have been converted to SI
    ! for use in EPOCH
    !
    ! For reference:
    ! c_hy_equil = 2/3 * (2*pi*kb)**-3/2 * q0**4 * SQRT(m0)/eps0**2

    INTEGER :: ix, iy, iz
    REAL(num) :: eq_Zst, eq_cou_log, eq_A, eq_Te, eq_Ti, eq_ni, eq_mi
    REAL(num) :: temp_diff, eq_term, dTe_fac, dTi_fac

    DO iz = 1,nz
      DO iy = 1,ny
        DO ix = 1,nx

          ! Copy out grid variables for clarity
          eq_Zst = ion_charge(ix,iy,iz)
          eq_cou_log = ion_cou_log(ix,iy,iz)
          eq_A = ion_A(ix,iy,iz)
          eq_Te = hy_Te(ix,iy,iz)
          eq_Ti = hy_Ti(ix,iy,iz)
          eq_ni = ion_ni(ix,iy,iz)
          eq_mi = eq_A * amu

          temp_diff = eq_Ti - eq_Te
          eq_term = c_hy_equil * eq_ni * eq_cou_log &
              * SQRT(eq_mi/(eq_Te*eq_mi + eq_Ti*m0)**3)

          ! To prevent over-shooting, cap the temperature change at half the
          ! temperature difference
          dTe_fac = MIN(dt * eq_term * eq_Zst**2, 0.5_num)
          dTi_fac = MIN(dt * eq_term * eq_Zst**3, 0.5_num)

          ! Apply temperature change
          hy_Te(ix,iy,iz) = eq_Te + dTe_fac * temp_diff
          hy_Ti(ix,iy,iz) = eq_Ti - dTi_fac * temp_diff

        END DO
      END DO
    END DO

    CALL field_bc(hy_Te, ng)
    CALL field_bc(hy_Ti, ng)

  END SUBROUTINE thermal_equilibration

#endif
END MODULE hy_heating
