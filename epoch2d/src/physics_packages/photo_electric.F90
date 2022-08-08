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
!
! ------------------------------------------------------------------------------
!
! This is a simplified module for calculating the photo-electric attenuation of
! K-alpha photons in solid-density Cu, for different background electron
! temperatures.
!
! Data comes from AVI-opacity data, provided by Suxing Hu (Rochester University)


MODULE photo_electric
#ifdef HYBRID

  USE partlist
  USE calc_df
  USE setup
  USE hy_shared

  IMPLICIT NONE

  REAL(num), PARAMETER :: ni_cu_solid = 8.4912e28_num

CONTAINS

  SUBROUTINE run_photo_electric_Cu

    ! Cycle through all photon species, and apply photo-electric attenuation. We
    ! use I = I0 * exp(-kappa * c * dt), where c*dt is the photon distance in
    ! timestep dt, I0 is the initial intensity, and I is the final intensity.
    !
    ! kappa is an attenuation co-efficient, and a fit to AVI-opacity data has
    ! been provided by Suxing Hu. This only works for k-alpha photons in copper.
    !
    ! The drop in intensity is implemented by reducing the macro-photon weight
    ! each time-step.

    INTEGER :: ispecies, isol, z_temp
    TYPE(particle), POINTER :: photon
    REAL(num) :: part_x, part_y, part_ni, part_te
    REAL(num) :: kappa

#ifdef HYBRID
    ! Loop over all background solids
    DO isol = 1, solid_count

      ! Identify if the solid is copper
      z_temp = solid_array(isol)%z
      IF (.NOT. z_temp == 29) CYCLE

      ! Loop over all photon species which experience attenuation
      DO ispecies = 1, n_species
        IF (.NOT. species_list(ispecies)%species_type == c_species_id_photon) &
            CYCLE
        IF (.NOT. species_list(ispecies)%attenuate) CYCLE

        ! Loop over all photons in this species, and apply attenuation
        photon => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(photon))

          ! Find the ion number density of the copper at the photon position
          part_x = photon%part_pos(1) - x_grid_min_local
          part_y = photon%part_pos(2) - y_grid_min_local
          CALL hy_grid_centred_var_at_particle(part_x, part_y, part_ni, &
              solid_array(isol)%ion_density)

          ! Ignore photons outside the solid isol
          IF (part_ni < 1.0_num) THEN
            photon => photon%next
            CYCLE
          END IF

          ! Interpolate background electron temperature to photon position
          CALL hy_grid_centred_var_at_particle(part_x, part_y, part_te, hy_te)

          ! Calculate attenuation co-efficient
          kappa = kappa_copper(part_te, part_ni)

          ! Apply attenuation
          photon%weight = photon%weight * EXP(-kappa*c*dt)

          ! Point to the next photon in the particle list
          photon => photon%next
        END DO

      END DO
    END DO
#endif

  END SUBROUTINE run_photo_electric_Cu



  FUNCTION kappa_copper(te, ni)

    ! Calculates the copper attenuation co-efficient for K-alpha photons
    ! according to an AVI-opacity data fit provided by Suxing Hu (Rochester
    ! University).
    !
    ! This function returns kappa in units of 1/m.

    REAL(num), INTENT(IN) :: te, ni
    REAL(num) :: te_ev, x1, x2, x3, x4, x5, kappa
    REAL(num) :: kappa_copper

    ! Obtain electron temperature in units of eV
    te_ev = te * kb / q0

    ! Fit to obtain kappa at solid density copper
    IF (te_ev < 500.0_num) THEN
      x1 = 154.633275098830_num
      x2 = 0.684356560314504_num
      x3 = -1.129249310157342E-003_num
      x4 = -9.808822131732512E-007_num
      x5 = 1.505648924357310E-009_num
      kappa = 100.0_num * (x1 + x2*te_ev + x3*te_ev**2 + x4*te_ev**3 &
          + x5*te_ev**4)
    ELSE
      x1 = 475.299874660280_num
      x2 = -0.781298994640132_num
      x3 = 4.530821133169535E-004_num
      x4 = -8.910063459566747E-008_num
      kappa = 100.0_num * (x1 + x2*te_ev + x3*te_ev**2 + x4*te_ev**3)
    END IF

    ! Fit is for solid-density copper, correct for the density of the current
    ! solid
    kappa_copper = kappa * ni / ni_cu_solid

  END FUNCTION kappa_copper

#endif
END MODULE photo_electric
