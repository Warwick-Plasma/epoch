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

  REAL(num), PRIVATE, ALLOCATABLE :: ohmic_heat_const(:,:)

  REAL(num), PRIVATE, ALLOCATABLE :: eff_heat_capacity(:,:)

  REAL(num), PRIVATE, PARAMETER :: c_hy_equil = &
      2.0_num/3.0_num*(2.0_num*pi*kb)**(-1.5_num) * q0**4 * SQRT(m0)/epsilon0**2

CONTAINS

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

    REAL(num), ALLOCATABLE :: t_prime(:,:), i_kb_sumne2(:,:)
    INTEGER :: i_sol

    ! Temporary varible to store T' array values
    ALLOCATE(t_prime(1-ng:nx+ng,1-ng:ny+ng))

    ! Cycle through solids and calculate C values
    DO i_sol = 1, solid_count
      t_prime = solid_array(i_sol)%z_prime * hy_te
      ALLOCATE(solid_array(i_sol)%heat_capacity(1-ng:nx+ng,1-ng:ny+ng))
      solid_array(i_sol)%heat_capacity = 0.3_num + 1.2_num * t_prime &
          * (2.2_num + t_prime)/(1.1_num + t_prime)**2
    END DO

    ! Deallocate temporary variable t_prime
    DEALLOCATE(t_prime)

    ! Get effective heat capacity
    ALLOCATE(eff_heat_capacity(1-ng:nx+ng,1-ng:ny+ng))
    eff_heat_capacity = 0.0_num
    DO i_sol = 1, solid_count
      eff_heat_capacity = eff_heat_capacity &
          + solid_array(i_sol)%el_density / solid_array(i_sol)%heat_capacity
    END DO

    ! Get 1/kb/Sum(ne)**2 values
    ALLOCATE(i_kb_sumne2(1-ng:nx+ng,1-ng:ny+ng))
    i_kb_sumne2 = 0.0_num
    DO i_sol = 1, solid_count
      i_kb_sumne2 = i_kb_sumne2 + solid_array(i_sol)%el_density
    END DO
    i_kb_sumne2 = 1.0_num/(kb*i_kb_sumne2**2)

    ! Precalculate the ohmic_heat_const array if Ohmic heating is running
    IF (use_ohmic) THEN
      ALLOCATE(ohmic_heat_const(1-ng:nx+ng,1-ng:ny+ng))
      ohmic_heat_const = dt * i_kb_sumne2
    END IF

    DEALLOCATE(i_kb_sumne2)

  END SUBROUTINE get_heat_capacity



  SUBROUTINE ohmic_heating

    ! Calculates the Ohmic heating of the simulation grid, as described by
    ! Davies, et al, (2002). Phys. Rev. E, 65(2), 026407. At this point in the
    ! simulation, we are at the end of a timestep, with E evaluated in the
    ! middle of this timestep. Assuming this is the average electric field over
    ! the timestep, we have a power disspiation of R*I**2, which is equivalent
    ! to a power density of J**2.(resistivity)
    !
    !  dT = (Power density) * dt  * Sum(ne/C) / [Sum(ne)]**2 / kB
    !
    ! ohmic_heat_const = dt / kB / Sum(ne)**2
    ! eff_heat_capacity = sum(ne/C)

    REAL(num) :: j2
    INTEGER :: ix, iy, i_sol

    ! Loop over all grid points to find temperature change
    DO iy = 1, ny
      DO ix = 1, nx

        ! Tb is a cell-centred variable, but J has stagger - need to average
        j2 = 0.25_num * ((jbx(ix,iy) + jbx(ix-1,iy))**2 &
            + (jby(ix,iy) + jby(ix,iy-1))**2) + jbz(ix,iy)**2

        ! Calculate Ohmic heating, avoiding the 0/0 NaN
        IF (eff_heat_capacity(ix,iy) > 0.0_num) THEN
          hy_te(ix,iy) = hy_te(ix,iy) &
              + j2*resistivity(ix,iy)*ohmic_heat_const(ix,iy) &
              * eff_heat_capacity(ix,iy)
        END IF
      END DO
    END DO

  END SUBROUTINE ohmic_heating



  SUBROUTINE clear_heat_capacity

    ! De-allocate the heat capacity in each solid. This array is recalculated
    ! each time-step anyway, so there's no need to save and resize it in the
    ! load balancer

    INTEGER :: i_sol

    DO i_sol = 1, solid_count
      DEALLOCATE(solid_array(i_sol)%heat_capacity)
    END DO
    DEALLOCATE(eff_heat_capacity)

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hy_te, ng)

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

    INTEGER :: ix, iy
    REAL(num) :: eq_zst, eq_cou_log, eq_a, eq_te, eq_ti, eq_ni, eq_mi
    REAL(num) :: temp_diff, eq_term, dte_fac, dti_fac

    DO iy = 1,ny
      DO ix = 1,nx

        ! Copy out grid variables for clarity
        eq_zst = ion_charge(ix,iy)
        eq_cou_log = ion_cou_log(ix,iy)
        eq_a = ion_a(ix,iy)
        eq_te = hy_te(ix,iy)
        eq_ti = hy_ti(ix,iy)
        eq_ni = ion_ni(ix,iy)
        eq_mi = eq_a * amu

        temp_diff = eq_ti - eq_te
        eq_term = c_hy_equil * eq_ni * eq_cou_log &
            * SQRT(eq_mi/(eq_te*eq_mi + eq_ti*m0)**3)

        ! To prevent over-shooting, cap the temperature change at half the
        ! temperature difference
        dte_fac = MIN(dt * eq_term * eq_zst**2, 0.5_num)
        dti_fac = MIN(dt * eq_term * eq_zst**3, 0.5_num)

        ! Apply temperature change
        hy_te(ix,iy) = eq_te + dte_fac * temp_diff
        hy_ti(ix,iy) = eq_ti - dti_fac * temp_diff

      END DO
    END DO

    CALL field_bc(hy_te, ng)
    CALL field_bc(hy_ti, ng)

  END SUBROUTINE thermal_equilibration

#endif
END MODULE hy_heating
