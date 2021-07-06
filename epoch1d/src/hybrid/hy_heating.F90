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

  REAL(num), PRIVATE, ALLOCATABLE :: ohmic_heat_const(:)

  REAL(num), PRIVATE, ALLOCATABLE :: eff_heat_capacity(:)

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

    REAL(num), ALLOCATABLE :: t_prime(:), i_kb_sumne2(:)
    INTEGER :: i_sol

    ! Temporary varible to store T' array values
    ALLOCATE(t_prime(1-ng:nx+ng))

    ! Cycle through solids and calculate C values
    DO i_sol = 1, solid_count
      t_prime = solid_array(i_sol)%z_prime * hy_te
      ALLOCATE(solid_array(i_sol)%heat_capacity(1-ng:nx+ng))
      solid_array(i_sol)%heat_capacity = 0.3_num + 1.2_num * t_prime &
          * (2.2_num + t_prime)/(1.1_num + t_prime)**2
    END DO

    ! Deallocate temporary variable t_prime
    DEALLOCATE(t_prime)

    ! Get effective heat capacity
    ALLOCATE(eff_heat_capacity(1-ng:nx+ng))
    eff_heat_capacity = 0.0_num
    DO i_sol = 1, solid_count
      eff_heat_capacity = eff_heat_capacity &
          + solid_array(i_sol)%el_density / solid_array(i_sol)%heat_capacity
    END DO

    ! Get 1/kb/Sum(ne)**2 values
    ALLOCATE(i_kb_sumne2(1-ng:nx+ng))
    i_kb_sumne2 = 0.0_num
    DO i_sol = 1, solid_count
      i_kb_sumne2 = i_kb_sumne2 + solid_array(i_sol)%el_density
    END DO
    i_kb_sumne2 = 1.0_num/(kb*i_kb_sumne2**2)

    ! Precalculate the ohmic_heat_const array if Ohmic heating is running
    IF (use_ohmic) THEN
      ALLOCATE(ohmic_heat_const(1-ng:nx+ng))
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
    INTEGER :: ix, i_sol

    ! Loop over all grid points to find temperature change
    DO ix = 1, nx

      ! Tb is a cell-centred variable, but J has stagger - need to average
      j2 = 0.25_num * (jbx(ix) + jbx(ix-1))**2  + jby(ix)**2 + jbz(ix)**2

      ! Calculate Ohmic heating, avoiding the 0/0 NaN
      IF (eff_heat_capacity(ix) > 0.0_num) THEN
        hy_te(ix) = hy_te(ix) + j2 * resistivity(ix) * ohmic_heat_const(ix) &
            * eff_heat_capacity(ix)
      END IF
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

#endif
END MODULE hy_heating
