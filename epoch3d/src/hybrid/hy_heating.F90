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

    REAL(num), ALLOCATABLE :: t_prime(:,:,:)
    INTEGER :: i_sol

    ! Temporary varible to store T' array values
    ALLOCATE(t_prime(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

    ! Cycle through solids and calculate C values
    DO i_sol = 1, solid_count
      t_prime = solid_array(i_sol)%z_prime * hy_te
      ALLOCATE(solid_array(i_sol)%heat_capacity(1-ng:nx+ng,1-ng:ny+ng, &
          1-ng:nz+ng))
      solid_array(i_sol)%heat_capacity = 0.3_num + 1.2_num * t_prime &
          * (2.2_num + t_prime)/(1.1_num + t_prime)**2
    END DO

    ! Deallocate temporary variable t_prime
    DEALLOCATE(t_prime)

    ! Get effective heat capacity
    ALLOCATE(eff_heat_capacity(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    eff_heat_capacity = 0.0_num
    DO i_sol = 1, solid_count
      eff_heat_capacity = eff_heat_capacity &
          + solid_array(i_sol)%el_density / solid_array(i_sol)%heat_capacity
    END DO

  END SUBROUTINE get_heat_capacity



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
