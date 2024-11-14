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
!-------------------------------------------------------------------------------
!
! This module contains scripts used to calculate the decay of high energy
! energy photons into muon-antimuon pairs. These routines are called from
! bremsstrahlung.F90, and describe the optical depth update, mu-/mu+ energy
! sampling, and angular deflections of the pair particles.
!
! These algorithms have been adapted from Geant4, and details are provided in
! the Geant4 Physics Reference Manual. The model we are using is valid for
! photon energies between 1.5 MeV and 100 GeV.

MODULE brem_muon
#ifdef BREMSSTRAHLUNG

  USE partlist
  USE particles
  USE calc_df
  USE setup

  IMPLICIT NONE

  REAL(num), PARAMETER, PRIVATE :: g4_s = -0.88_num
  REAL(num), PARAMETER, PRIVATE :: classic_muon_radius = &
      q0**2 / (4.0_num * pi * epsilon0 * m_mu * c**2)
  REAL(num), PARAMETER, PRIVATE :: two_muon_rest = 2.0_num * m_mu * c**2
  REAL(num), PARAMETER, PRIVATE :: inv_gev = 1.0_num / (1.0e9_num * q0)
  REAL(num), PARAMETER, PRIVATE :: wm_fac = &
      1.0_num / (4.0_num * EXP(0.5_num) * m_mu * c**2 * inv_gev)

  REAL(num), ALLOCATABLE, PRIVATE :: sig_fac(:), g4_ec(:), w_sat_s(:)
  INTEGER, PRIVATE :: z2ind(100)

  ! Variables evaluted in cross-section calculation, and reused in pair-creation
  REAL(num), PRIVATE :: g4_dn, a_027

  INTEGER, PARAMETER, PRIVATE :: loop_break = 1000

  ! Kinematic variables for the muon-antimuon pair
  REAL(num), PRIVATE :: g4_gamma_minus, g4_gamma_plus
  REAL(num), PRIVATE :: g4_theta_minus, g4_theta_plus 
  REAL(num), PRIVATE :: g4_phi, g4_phi0

 CONTAINS

  SUBROUTINE init_brem_muon(size_z, z_values, z_to_index)

    ! The muon pair-production code uses functions which are only Z-dependent so
    ! it makes no sense to recalculate these each step if the Z materials
    ! present don't change over a simulation. Hence, we will pre-calculate some
    ! functions and store them in private arrays
    !
    ! Here, z_values is an array containing all the unique z values present in
    ! the simulation, z2index is an array such that z2index(z) returns the index
    ! in z_values corresponding to z, and size_z is the size of z_values

    INTEGER, INTENT(IN) :: z_values(:)
    INTEGER, INTENT(IN) :: z_to_index(100)
    INTEGER, INTENT(IN) :: size_z
    INTEGER :: iz
    REAL(num) :: z_ion, g4_b

    z2ind = z_to_index

    ! Allocate pre-calculated arrays, same indexing as z_values
    ALLOCATE(sig_fac(size_z), g4_ec(size_z), w_sat_s(size_z))

    ! Loop over all possible atomic numbers, calculate terms if this number is
    ! present
    DO iz = 1,size_z
        z_ion = REAL(z_values(iz), num)

        sig_fac(iz) = 28.0_num * alpha * z_ion**2 * classic_muon_radius**2 &
            / 9.0_num

        IF (z_values(iz) == 1) THEN 
          g4_b = 202.4_num 
        ELSE
          g4_b = 183_num 
        END IF
        g4_ec(iz) = -18.0_num + 4347.0_num / g4_b * z_ion**(1.0_num / 3.0_num)
        w_sat_s(iz) = (g4_b * z_ion**(1.0_num / 3.0_num) * 4.0_num &
            * EXP(0.5_num) * m_mu**2 / m0 * c**2 * inv_gev)**g4_s
    END DO

  END SUBROUTINE init_brem_muon



  SUBROUTINE brem_muon_update_depth

    ! Cycles through all photon species and updates the muon pair-production
    ! optical depths

    INTEGER :: iz, z_temp, ispecies
    REAL(num), ALLOCATABLE :: grid_num_density_ion(:)
    TYPE(particle), POINTER :: photon, next_photon
    REAL(num) :: part_x
    REAL(num) :: a_temp, cross_sec, part_ni, delta_opdep

    ALLOCATE(grid_num_density_ion(1-ng:nx+ng))

    ! Calculate the number density of each ion species
    DO iz = 1, n_species
      ! Identify if the charge is greater than 1
      z_temp = species_list(iz)%atomic_no
      IF (z_temp < 1 .OR. z_temp > 100) CYCLE

      a_temp = species_list(iz)%mass / amu

      CALL calc_number_density(grid_num_density_ion, iz)
      CALL field_bc(grid_num_density_ion, ng)

      ! Update the optical depth for each photon species
      DO ispecies = 1, n_species
        ! Only update optical_depth_bremsstrahlung for the photon species
        IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
          ! Cycle through all photons in this species
          photon => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(photon))
            ! Remember the next photon in the list, as the current photon could
            ! be removed
            next_photon => photon%next

            ! Photon energy is calculated on creation, and this energy never
            ! changes in the current version of the code (28/Nov/2023). Check it
            ! is high enough to pair produce.
            IF (photon%particle_energy < two_muon_rest) THEN
              photon => next_photon
              CYCLE
            END IF

            ! Calculate the cross section at this photon energy
            cross_sec = calc_muon_cross_sec(z_temp, a_temp, &
                photon%particle_energy)

            ! Get background number density at photon
            part_x = photon%part_pos - x_grid_min_local
            CALL grid_centred_var_at_particle_bm(part_x, part_ni, &
                grid_num_density_ion)

            ! Update the photon optical depth
            delta_opdep = part_ni * cross_sec * c * dt * boost_muons
            photon%optical_depth_brem_muon = &
                photon%optical_depth_brem_muon - delta_opdep

            ! If photon optical depth drops below 0, create an e-/e+ pair and
            ! remove the photon from the simulation
            IF (photon%optical_depth_brem_muon <= 0.0_num) THEN
              CALL set_pair_kinematics(z_temp, a_temp, photon)
              CALL generate_pair(photon, &
                  bethe_heitler_muon_species, &
                  bethe_heitler_antimuon_species, ispecies)
            END IF

            photon => next_photon
          END DO
        END IF
      END DO
    END DO

    DEALLOCATE(grid_num_density_ion)

  END SUBROUTINE brem_muon_update_depth



  FUNCTION calc_muon_cross_sec(z_int, a_ion, e_gamma)

    ! Calculates the muon pair-production cross section in the presence of an 
    ! ion species with atomic number z_int, mass number a_ion, for a photon of 
    ! energy e_gamma.
    !
    ! Parametrisation taken from Geant4 source code, which acts as a fit to
    ! cross section data. Variable names chosen to match Geant4 physics 
    ! reference manual, with g4_ prefix

    INTEGER, INTENT(IN) :: z_int
    REAL(num), INTENT(IN) :: a_ion, e_gamma
    REAL(num) :: e_gamma_gev
    REAL(num) :: g4_wm, g4_cf, g4_t, g4_eg
    REAL(num) :: calc_muon_cross_sec

    a_027 = a_ion**0.27_num 
    IF (z_int == 1) THEN 
      g4_dn = 1.49_num
    ELSE
      g4_dn = 1.54_num * a_027
    END IF

    e_gamma_gev = e_gamma * inv_gev

    ! Combine for the total cross section
    g4_wm = wm_fac / g4_dn
    g4_cf = 1.0_num &
        + 0.04_num * LOG(1.0_num + g4_ec(z2ind(z_int)) /  e_gamma_gev)
    g4_t = 1.479_num + 0.00799_num * g4_dn 
    g4_eg = (1.0_num - 4.0_num * m_mu * c**2 / e_gamma)**g4_t &
        * (w_sat_s(z2ind(z_int)) + e_gamma_gev**g4_s)**(1.0_num / g4_s)
        
    calc_muon_cross_sec = sig_fac(z2ind(z_int)) * LOG(1 + g4_wm * g4_cf * g4_eg)

  END FUNCTION calc_muon_cross_sec



  SUBROUTINE generate_pair(photon, imuon, ianti, iphoton)

    ! This subroutine is called when "photon" undergoes muon pair production. 
    ! Before this subroutine is called, the code runs set_pair_kinematics, which
    ! sets the module variables describing the energies and directions of the 
    ! muon pair. Macro-particles for mu- and mu+ are created here, and the 
    ! photon is removed from the simulation

    TYPE(particle), POINTER :: photon
    INTEGER, INTENT(IN) :: imuon, ianti, iphoton
    TYPE(particle), POINTER :: antimuon, muon
    REAL(num) :: energy_mu_p, momentum_mu_p, phi_p
    REAL(num) :: energy_mu_m, momentum_mu_m, phi_m
    REAl(num) :: norm, dir_x, dir_y, dir_z

    ! Create pair at photon position
    CALL create_particle(antimuon)
    CALL create_particle(muon)
    antimuon%part_pos = photon%part_pos
    muon%part_pos = photon%part_pos

    ! e- and e+ have the same weights as generating photon
    antimuon%weight = photon%weight / boost_muons
    muon%weight = photon%weight / boost_muons

    ! Calculate momentum magnitude going to each particle
    energy_mu_p = g4_gamma_plus * m_mu * c**2
    energy_mu_m = g4_gamma_minus * m_mu * c**2
    momentum_mu_p = SQRT(energy_mu_p**2 - (m_mu * c**2)**2)/c
    momentum_mu_m = SQRT(energy_mu_m**2 - (m_mu * c**2)**2)/c

    ! Calculate momentum direction of generating photon
    norm  = c / photon%particle_energy
    dir_x = photon%part_p(1) * norm
    dir_y = photon%part_p(2) * norm
    dir_z = photon%part_p(3) * norm

    ! Assign momenta to particles
    antimuon%part_p(1:3) = momentum_mu_p * (/ dir_x, dir_y, dir_z /)
    muon%part_p(1:3) = momentum_mu_m * (/ dir_x, dir_y, dir_z /)

    ! Rotate momenta for each particle
    phi_p = g4_phi0 + 0.5_num * g4_phi 
    phi_m = g4_phi0 - 0.5_num * g4_phi + pi
    CALL rotate_p(antimuon, COS(g4_theta_plus), phi_p, momentum_mu_p)
    CALL rotate_p(muon, COS(g4_theta_minus), phi_m, momentum_mu_m)

    ! Set particle energies
    antimuon%particle_energy = energy_mu_p
    muon%particle_energy = energy_mu_m

    ! Save pair to particle list
    CALL add_particle_to_partlist(species_list(ianti)%attached_list, antimuon)
    CALL add_particle_to_partlist(species_list(imuon)%attached_list, muon)

    ! Remove photon
    IF (random() < 1.0_num / boost_muons) THEN
      CALL remove_particle_from_partlist(species_list(iphoton)%attached_list, &
          photon)
      CALL destroy_particle(photon)
    ELSE
      photon%optical_depth_brem_muon = -LOG(random())
    END IF

  END SUBROUTINE generate_pair



  SUBROUTINE set_pair_kinematics(z_int, a_ion, photon)

    ! This subroutine sets module parameters to describe the energy and momentum
    ! direction of the mu-/mu+ pair. The method has been adapted from the
    ! Geant4 physics reference manual, where parameters t, psi and rho are
    ! sampled first, then the deflection angles are calculated.
    !
    ! Where the Geant4 source-code disagrees with the physics reference manual, 
    ! the code implementation has been chosen. Variable names have been chosen
    ! to match the Geant4 documentation or source-code where appropriate, and 
    ! are assigned the pre-fix g4_

    INTEGER, INTENT(IN) :: z_int
    REAL(num), INTENT(IN) :: a_ion 
    TYPE(particle), POINTER, INTENT(IN) :: photon
    INTEGER :: i, itheta
    REAL(num) :: e_photon, gamma_tot, inv_gamma, z_ion, z_13, inv_z_13
    REAL(num) :: delta_max, w_term_1, w_term_2, log_w_max, scaled_z13
    REAL(num) :: x_pm, a32_b1, a3_t2, a3t2_b1, f1_fac, inv_gamma_pmt, psi_fac
    REAL(num) :: g4_b, g4_b1, g4_x_min, g4_x_max, g4_w, g4_w_inf, g4_w_max
    REAL(num) :: g4_a3, g4_x_plus, g4_x_minus, g4_delta, g4_f1_max, g4_t, g4_f1 
    REAL(num) :: g4_f2, g4_f2_max, g4_psi, g4_c2, g4_rho, g4_rho_max2, g4_u

    ! Derived variables
    e_photon = SQRT(photon%part_p(1)**2 + photon%part_p(2)**2 &
        + photon%part_p(3)**2) * c
    gamma_tot = e_photon / (m_mu * c**2) 
    inv_gamma = 1.0_num / gamma_tot
    z_ion = REAL(z_int, num)
    z_13 = z_ion ** (1.0_num / 3.0_num)
    inv_z_13 = 1.0_num / z_13

    IF (z_ion == 1) THEN
      g4_b = 202.4_num
    ELSE
      g4_b = 183_num
    END IF

    ! Energy split
    g4_x_min = 0.5_num - sqrt(0.25_num - inv_gamma)
    g4_x_max = 0.5_num + sqrt(0.25_num - inv_gamma)

    ! Pre-calculate variables for x+, x- calculation
    g4_w_inf = g4_b * m_mu / (z_13 * g4_dn * m0)
    delta_max = 2.0_num * m_mu**2 * c**2 / e_photon
    w_term_1 = (g4_dn * EXP(0.5_num) - 2.0_num) / m_mu
    w_term_2 = g4_b * inv_z_13 * EXP(0.5_num) / m0
    g4_w_max = g4_w_inf * (1.0_num + delta_max * w_term_1) &
        / (1.0_num + delta_max * w_term_2)
    log_w_max = LOG(g4_w_max)

    ! Pre-calculate variable for t sampling
    g4_b1 = 1.0_num / (4.0_num * (0.138_num * a_027)**2)

    ! Pre-calculate variable for rho sampling
    scaled_z13 = m0 * z_13 / (183.0_num * m_mu)

    ! Keep sampling until we have physical deflection angles, or hit loop_break
    DO itheta = 0, loop_break

      ! Rejection sample x_plus and x_minus
      DO i = 0, loop_break
        g4_x_plus = x_min + random() * (g4_x_max - g4_x_min)
        g4_x_minus = 1.0_num - g4_x_plus
        x_pm = g4_x_plus * g4_x_minus

        g4_delta = 0.5_num * m_mu * inv_gamma / x_pm
        g4_w = g4_w_inf * (1.0_num + g4_delta * w_term_1) &
            / (1.0_num + g4_delta * w_term_2)

        ! Using condition in G4 source-code, not PRM
        IF ((1.0_num - 4.0_num * x_pm / 3.0_num) * LOG(g4_w) / log_w_max &
            > random()) THEN
          EXIT
        END IF
      END DO

      ! Sample auxilary t using G4 source-code method, not PRM
      g4_a3 = 0.5_num * inv_gamma / x_pm
      a32_b1 = g4_a3**2 + g4_b1
      g4_f1_max= -0.5_num * (1.0_num - x_pm) &
          * (2.0_num * g4_b1 + (a32_b1 + g4_a3**2) * LOG(g4_a3**2 / a32_b1)) &
          / g4_b1**3

      ! Rejection sample for t
      DO i = 0, loop_break
        g4_t = random()
        a3_t2 = (g4_a3 / g4_t)**2
        a3t2_b1 = a3_t2 + g4_b1
        f1_fac = 1.0_num - 2.0_num * x_pm &
            + 4.0_num * x_pm * g4_t * (1.0_num - g4_t)
        IF (ABS(g4_b1) < 0.0001_num * a3_t2) THEN
          ! Special case of a3_t2 = a3t2_b1 because of logarithm accuracy
          g4_f1 = f1_fac / (12.0_num * a3_t2**4)
        ELSE
          g4_f1 = - 0.5_num * f1_fac &
              * (2.0_num*g4_b1 + (a3t2_b1 + a3_t2)*LOG(a3_t2/a3t2_b1)) &
              / g4_b1**3
        END IF
        
        ! Check if f1 is outside the acceptable range
        IF (g4_f1 < 0.0_num .OR. g4_f1 > g4_f1_max) THEN
          CYCLE
        END IF

        ! Rejection test for sampled t
        IF (random() * g4_f1_max <= g4_f1) THEN
          EXIT
        END IF
      END DO

      ! Sample auxilary psi         
      g4_f2_max = 1.0_num - 2.0_num * x_pm &
          * (1.0_num - 4.0_num * g4_t * (1.0_num - g4_t))

      ! Generate psi by the rejection method
      DO i = 0, loop_break
        g4_psi = random() * 2.0_num * pi
        g4_f2 = f1_fac * (1.0_num + COS(2.0_num * g4_psi))

        ! Check if f2 is outside the acceptable range
        IF (g4_f2 < 0.0_num .OR. g4_f2 > g4_f2_max) THEN
          CYCLE
        END IF

        ! Successful f2 sample
        IF (random() * g4_f2_max <= g4_f2) THEN
          EXIT
        END IF
      END DO

      ! Sample auxilary rho by direct transformation
      inv_gamma_pmt = inv_gamma / (2.0_num * x_pm * g4_t)
      g4_c2 = 4.0_num * (inv_gamma_pmt**2 + scaled_z13**2)**2 / SQRT(x_pm)
      g4_rho_max2 = 1.9_num / a_027 * (1.0_num / g4_t - 1.0_num)
      g4_rho = (g4_c2 * ((1.0_num + g4_rho_max2**2 / g4_c2)**random() &
          - 1.0_num))**0.25_num

      ! Sample kinematic variables
      g4_gamma_plus = g4_x_plus * gamma_tot
      g4_gamma_minus = g4_x_minus * gamma_tot

      ! Temporary fix for sampling issue - can rarely get gamma < 1
      ! Should probably re-sample to do it correctly? 
      ! Is our sampling even right if a gamma < 1 for 464 MeV photon is allowed?
      IF (g4_gamma_plus < 1.0_num .OR. g4_gamma_minus < 1.0_num) THEN 
        g4_gamma_plus = 0.5_num * gamma_tot 
        g4_gamma_minus = 0.5_num * gamma_tot 
      END IF

      g4_u = SQRT(1.0_num / g4_t - 1.0_num)
      psi_fac = 0.5_num * g4_rho * COS(g4_psi)
      g4_theta_plus = (g4_u + psi_fac) / g4_gamma_plus 
      g4_theta_minus = (g4_u - psi_fac) / g4_gamma_minus
      g4_phi = g4_rho / g4_u * SIN(g4_psi)
      g4_phi0 = 2.0_num * pi * random()

      ! Check for physical theta parameters
      IF (ABS(g4_theta_plus) < pi .AND. ABS(g4_theta_minus) < pi) THEN 
        RETURN
      END IF
    END DO

  END SUBROUTINE set_pair_kinematics



  SUBROUTINE grid_centred_var_at_particle_bm(part_x, part_var, &
      grid_var)

    ! Calculates the value of a grid-centred variable part_var stored in the
    ! grid grid_var, averaged over the particle shape for a particle at position
    ! part_x

    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(IN) :: grid_var(1-ng)
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
#include "bspline3/part_var.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#include "tophat/part_var.inc"
#else
#include "triangle/gx.inc"
#include "triangle/part_var.inc"
#endif
    part_var = fac * part_var

  END SUBROUTINE grid_centred_var_at_particle_bm

#endif
END MODULE brem_muon
