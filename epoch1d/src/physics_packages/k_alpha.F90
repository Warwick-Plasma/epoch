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

MODULE k_alpha
#ifdef K_ALPHA

  USE partlist
  USE calc_df
  USE setup

  IMPLICIT NONE

  REAL(num), PARAMETER :: rbeb_const = 2.0_num * pi * a0**2 * alpha**4
  REAL(num), ALLOCATABLE :: rbeb_ke(:), rbeb_1s(:,:)
  REAL(num) :: binding_energy(0:28)

  INTEGER, PARAMETER :: sample_el_in = 100

  TYPE(interpolation_state), SAVE :: last_state

CONTAINS

  SUBROUTINE setup_k_alpha_module

    ! Confirm that the correct particles are present to actually use the
    ! k_alpha routines, and warn the user if they're not. Then call
    ! subroutines to initialise the k_alpha table array, and give all
    ! particles an optical depth for k_alpha emission (if appropriate)

    INTEGER :: i_sol, iu, io, ispecies
    LOGICAL :: found

    CALL calculate_cross_sections_rbeb

    ! Check if the correct background solid is present
#ifdef HYBRID
    IF (use_hybrid .AND. produce_k_alpha_photons) THEN

      ! Look for a copper solid
      found = .FALSE.
      DO i_sol = 1, solid_count
        IF (solid_array(i_sol)%z == 29) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO

      ! Print warning if no copper is detected
      IF (.NOT. found) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'No copper solid has been detected. K-alpha photons ', &
                'are currently only'
            WRITE(io,*) ' calculated for hot electrons moving through copper. '
          END DO
        END IF
      END IF
    END IF
#endif

    ! All electrons initially present in the simulation window will be assigned
    ! an optical depth, after passing which, K-alpha photons will be emitted
    IF (.NOT. ic_from_restart) THEN
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%species_type == c_species_id_electron) THEN
          CALL initialise_optical_depth(species_list(ispecies))
        END IF
      END DO
    END IF

  END SUBROUTINE setup_k_alpha_module



  SUBROUTINE calculate_cross_sections_rbeb

    ! This subroutine is called to calculate collisional ionisation cross
    ! sections of the 1s shell in copper atoms and ions. Calculated variables
    ! are saved to the following arrays:
    ! - rbeb_ke(i_ke), array of incident electron kinetic energies
    ! - rbeb_1s(i_ke, i_z), 2d array, where each (i_ke, i_z) pair gives the 1s
    !     cross section for a copper atom or ion, with charge state i_z
    !
    ! Method described in Kim et al, "Extension of the binary-encounter-dipole
    ! model to relativistic incident electrons" Phys. Rev. A 62(5) (2000).
    !
    ! This subroutine uses cross section (22) with the pre-factor of (37) from
    ! Kim (2000). This model works best for the inner shells of heavy nuclei,
    ! and is well suited for the 1s shell of Cu. The approximation B = U
    ! (binding energy = average orbital kinetic energy) has been made, due to a
    ! lack of U data

    INTEGER :: i_ke, i_z
    REAL(num) :: min_ke_ev, max_ke_ev, el_ke
    REAL(num) :: rbeb_t, rbeb_tp, rbeb_bp, rbeb_up
    REAL(num) :: beta_t2, beta_b2, beta_u2, beta2
    REAL(num) :: rbeb_pre, rbeb_fac, log_2bp, rbeb_1_middle, rbeb_1
    REAL(num) :: one_2tp, tp2_frac, ln_t_term, rbeb_2, rbeb_3, n_1s
    REAL(num) :: sig_add, sig_rbeb

    ! Binding energies of the copper 1s shell [eV]. Shell binding energies for
    ! atomic copper have been tabulated by Desclaux (1973), At. Data Nucl. Data
    ! Tables 12 (4). Ionic copper binding energies are deduced using the method
    ! of Carlson (1970) At. Data Nucl. Data Tables 2.
    ! - binding_energy(i_z) gives the 1s binding energy for charges state i_z
    binding_energy = (/9051.99_num, 9059.01_num, 9075.56_num, 9096.09_num, &
        9118.51_num, 9141.71_num, 9177.71_num, 9204.71_num, 9237.08_num, &
        9277.53_num, 9310.66_num, 9325.62_num, 9359.62_num, 9394.62_num, &
        9441.72_num, 9477.32_num, 9514.21_num, 9544.17_num, 9582.28_num, &
        9754.17_num, 9863.67_num, 9981.67_num, 10107.7_num, 10243.1_num,&
        10391.8_num, 10397.1_num, 10505.0_num, 11062.4_num, 11567.6_num/)

    ! Covert to SI units
    binding_energy = binding_energy * q0

    ! Range of table and number of sample points. In Kim (2000), agreement
    ! between RBEB and experiment is verified for some elements between ~10 eV
    ! to ~10 MeV
    min_ke_ev = binding_energy(0)/q0
    max_ke_ev = 1.0e9_num

    ! Electron kinetic energies [J]
    ALLOCATE(rbeb_ke(sample_el_in))
    ALLOCATE(rbeb_1s(sample_el_in, 0:28))

    ! Calculate electron energy and cross section at each table sample point
    DO i_ke = 1, sample_el_in

      ! Logarithmically space energy sample points between limits
      el_ke = min_ke_ev * q0 * &
         (max_ke_ev / min_ke_ev)**(REAL(i_ke-1,num)/REAL(sample_el_in-1,num))
      rbeb_ke(i_ke) = el_ke

      ! Caclulate cross section contribution from the 1s shell of copper
      DO i_z = 0, 28

        ! Only consider shells with binding energy less than the incident
        ! electron kinetic energy
        IF (el_ke < binding_energy(i_z)) CYCLE

        ! t
        rbeb_t = el_ke / binding_energy(i_z)

        ! t', b', u' (with approximation B = U)
        rbeb_tp = el_ke / (m0*c**2)
        rbeb_bp = binding_energy(i_z) / (m0*c**2)
        rbeb_up = rbeb_bp

        ! beta_t**2, beta_b**2, beta_u**2
        beta_t2 = 1.0_num - 1.0_num/(1.0_num + rbeb_tp)**2
        beta_b2 = 1.0_num - 1.0_num/(1.0_num + rbeb_bp)**2
        beta_u2 = 1.0_num - 1.0_num/(1.0_num + rbeb_up)**2

        ! Averaging pre-factor of (37) [Kim (2000)]
        beta2 = beta_t2 + beta_b2 + beta_u2
        rbeb_pre = 0.5_num * (1.0_num + beta2/beta_t2)

        ! Number of electrons in the 1s shell
        IF (i_z < 28) THEN
          n_1s = 2.0_num
        ELSE
          n_1s = 1.0_num
        END IF

        ! Calculate RBEB cross section terms
        rbeb_fac = rbeb_const * n_1s / (beta2 * rbeb_bp)

        log_2bp = LOG(2.0_num * rbeb_bp)
        rbeb_1_middle = LOG(beta_t2 / (1.0_num - beta_t2)) - beta_t2 - log_2bp
        rbeb_1 = 0.5_num * rbeb_1_middle * (1.0_num - 1.0_num / rbeb_t**2)

        one_2tp = 1.0_num + 2.0_num * rbeb_tp
        tp2_frac = 1.0_num / (1.0_num + 0.5_num * rbeb_tp)**2
        ln_t_term = LOG(rbeb_t) / (rbeb_t + 1.0_num) * one_2tp * tp2_frac
        rbeb_2 = 1.0_num - (1.0_num / rbeb_t) - ln_t_term

        rbeb_3 = rbeb_bp**2 * tp2_frac * (rbeb_t - 1.0_num) / 2.0_num

        ! Combine to RBEB cross section
        rbeb_1s(i_ke, i_z) = rbeb_pre * rbeb_fac * (rbeb_1 + rbeb_2 + rbeb_3)
      END DO
    END DO

  END SUBROUTINE calculate_cross_sections_rbeb



  FUNCTION check_k_alpha_variables()

    ! Function to check all required species are present

    INTEGER :: check_k_alpha_variables
    INTEGER :: io, iu
    INTEGER :: ispecies
    INTEGER :: first_electron = -1

    check_k_alpha_variables = c_err_none

    ! No special species required if k_alpha is turned off
    IF (.NOT.use_k_alpha) RETURN

    ! No special species required if we only do radiation reaction
    IF (.NOT.produce_k_alpha_photons) RETURN

    ! Identify if there exists any electron or positron species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .AND. first_electron == -1) THEN
        first_electron = ispecies
      END IF
    END DO

    ! Print warning if there is no electron species
    IF (first_electron < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No electron species specified.'
          WRITE(io,*) 'Specify using "identify:electron".'
          WRITE(io,*) 'k_alpha routines require at least one ' , &
              'species of electrons.'
        END DO
      END IF
      check_k_alpha_variables = c_err_missing_elements
      RETURN
    END IF

#ifdef PHOTONS
    ! photon_species can act as k_alpha_photon_species if no
    ! k_alpha species is defined
    IF (k_alpha_photon_species == -1 &
        .AND. .NOT. photon_species == -1 ) THEN
      k_alpha_photon_species = photon_species
    END IF
#endif

    ! Check if there exists a species to populate with k_alpha photons
    IF (k_alpha_photon_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No photon species specified. Specify using ', &
              '"identify:k_alpha"'
        END DO
      END IF
      check_k_alpha_variables = c_err_missing_elements
      RETURN
    END IF

  END FUNCTION check_k_alpha_variables



  SUBROUTINE initialise_optical_depth(current_species)

    ! Loops through all particles in a species and sets an optical depth for
    ! k_alpha emission

    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    ! Randomly sample the optical depth from an exponential distribution
    current => current_species%attached_list%head
    DO WHILE(ASSOCIATED(current))
      p_tau = random()
      current%optical_depth_k_alpha = -LOG(1.0_num - p_tau)
      current => current%next
    END DO

  END SUBROUTINE initialise_optical_depth



  SUBROUTINE k_alpha_update_optical_depth

    ! Updates the optical depth for electrons. This subroutine is responsible for
    ! also calling the function which calculates the optical depth change, and
    ! calling the generate_photon subroutine. This subroutine serves as the main
    ! interface to the k_alpha module for main-loop processes in
    ! epoch3d.F90

    ! Currently empty, K-alpha routines have no support for normal PIC codes,
    ! only hybrid-PIC.

  END SUBROUTINE k_alpha_update_optical_depth



#ifdef HYBRID
  SUBROUTINE hy_k_alpha_update_optical_depth

    ! Analogous to k_alpha_update_optical_depth, but tracks the passage of
    ! electrons through the hybrid solid type instead of ion macroparticles.
    ! This subroutine serves as the main interface to the k_alpha module for
    ! the hybrid PIC loop in hybrid.F90

    INTEGER :: isol, ispecies, iz, z_temp, z_st, part_type
    TYPE(particle), POINTER :: current
    REAL(num) :: part_ux, part_uy, part_uz, gamma_rel
    REAL(num) :: part_x, part_e, part_ke, part_ni, part_v
    REAL(num) :: part_te, part_z_st

    ! Consider the optical depth change due to each solid separately
    DO isol = 1, solid_count

      ! Identify if the solid is copper
      z_temp = solid_array(isol)%z
      IF (.NOT. z_temp == 29) CYCLE

      ! Update the optical depth for each electron species
      DO ispecies = 1, n_species
        part_type = species_list(ispecies)%species_type
        IF (part_type == c_species_id_electron) THEN
          ! Cycle through all particles in this species
          current => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(current))

            ! Get grid properties at electron position
            part_x = current%part_pos - x_grid_min_local

            ! Background number density
            CALL grid_centred_var_at_particle(part_x, part_ni, &
                solid_array(isol)%ion_density)

            ! Ignore electrons which are outside the copper target
            IF (part_ni < 1.0_num) THEN
              current => current%next
              CYCLE
            END IF

            ! Get electron kinetic energy
            part_ux = current%part_p(1) / mc0
            part_uy = current%part_p(2) / mc0
            part_uz = current%part_p(3) / mc0
            gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
            part_e = gamma_rel * m0 * c**2
            current%particle_energy = part_e
            part_ke = part_e - m0*c**2

            ! Don't update the optical depth if the particle hasn't moved
            IF (gamma_rel - 1.0_num < 1.0e-15_num) THEN
              current => current%next
              CYCLE
            END IF

            ! Get electron speed
            part_v = SQRT(current%part_p(1)**2 + current%part_p(2)**2 &
                + current%part_p(3)**2) * c**2 / part_e

            ! Determine the average ion charge state at this position
            CALL grid_centred_var_at_particle(part_x, part_z_st, ion_charge)
            z_st = NINT(part_z_st)

            ! Fully ionised copper cannot produce K-alpha emissions
            IF (z_st == 29) THEN
              current => current%next
              CYCLE
            END IF

            ! Background plasma temperature
            CALL grid_centred_var_at_particle(part_x, part_te, hy_te)

            ! Calculate the optical depth change for the current electron, due
            ! to moving through the current solid
            current%optical_depth_k_alpha = &
                current%optical_depth_k_alpha &
                - delta_optical_depth(part_v, part_ni, z_st, part_ke)

            ! If optical depth dropped below zero generate photon and reset
            ! optical depth
            IF (current%optical_depth_k_alpha <= 0.0_num) THEN
              CALL generate_photon(current, k_alpha_photon_species, part_te, &
                  z_st)
              current%optical_depth_k_alpha = reset_optical_depth()
            END IF
            current => current%next
          END DO
        END IF
      END DO
    END DO

  END SUBROUTINE hy_k_alpha_update_optical_depth
#endif



  FUNCTION delta_optical_depth(part_v, part_ni, z_st, part_ke)

    ! Calculate the change in optical depth during this timestep, given by:
    ! (ion number density) * (emission cross section) * (distance traversed)
    ! Here, cross sections are determined using the RBEB model for 1s
    ! ionisation
    !
    ! This function uses the electron speed (part_v), the background ion number
    ! density of the current solid at the particle positon (part_ni) and its
    ! ionisation state (z_st). The particle kinetic energy is also used
    ! (part_ke)

    REAL(num), INTENT(IN) :: part_v, part_ni, part_ke
    INTEGER, INTENT(IN) :: z_st
    REAL(num) :: cross_sec_val
    REAL(num) :: delta_optical_depth

    ! Interpolate cross section value for correct particle type
    cross_sec_val = find_value_from_table_1d_coll(part_ke, sample_el_in, &
        rbeb_ke, rbeb_1s(:,z_st), last_state)

    ! Optical depth update, optionally boosted by user factor k_alpha_weight
    delta_optical_depth = part_ni * cross_sec_val * part_v * dt / k_alpha_weight

  END FUNCTION delta_optical_depth



  FUNCTION reset_optical_depth()

    ! Draws a new random number for the exponentially distributed optical depths

    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth



  SUBROUTINE generate_photon(electron, iphoton, part_te, z_st)

    ! Generates a K-alpha photon, which moves in a direction randomly sampled
    ! from an isotropic distribution. Inputs include the generating electron
    ! particle, and the K-alpha photon species index, iphoton. The plasma
    ! plasma temperature at the electron position (part_te) is also required for
    ! energy sampling

    TYPE(particle), POINTER :: electron
    INTEGER, INTENT(IN) :: iphoton, z_st
    REAL(num), INTENT(IN) :: part_te
    INTEGER :: brem_index
    REAL(num) :: dir_x, dir_y, dir_z, mag_p
    REAL(num) :: rand_temp, photon_energy, photon_p, part_e
    REAL(num) :: el_e_final, p_mag_final
    REAL(num) :: scatter_cos_theta, scatter_phi
    TYPE(particle), POINTER :: new_photon
    LOGICAL :: add_photon

    ! Obtain electron direction (magnitude must be > 0 to prevent 1/0 issues)
    mag_p = MAX(SQRT(electron%part_p(1)**2 + electron%part_p(2)**2 &
        + electron%part_p(3)**2), c_tiny)
    dir_x = electron%part_p(1) / mag_p
    dir_y = electron%part_p(2) / mag_p
    dir_z = electron%part_p(3) / mag_p
    part_e = SQRT(mag_p**2 * c**2 + m0**2 * c**4)

    ! Calculate electron recoil (subtract weighted photon momentum)
    IF (use_k_alpha_recoil) THEN
      el_e_final = part_e - k_alpha_weight * binding_energy(z_st)
      p_mag_final = SQRT(el_e_final**2 - m0**2 * c**4) / c
      electron%part_p(1) = p_mag_final * dir_x
      electron%part_p(2) = p_mag_final * dir_y
      electron%part_p(3) = p_mag_final * dir_z
#if defined(HYBRID) || defined(PHOTONS) || defined(BREMSSTRAHLUNG)
      electron%particle_energy = el_e_final
#endif
    END IF

    ! The de-excitation may lead to the creation of an Auger electron instead of
    ! a K-alpha photon. At present, this probability is a user-defined variable,
    ! and no Auger electrons are added to the simulation window
    IF (random() < auger_frac) THEN
      RETURN
    END IF

    ! Determine photon energy
    photon_energy = calculate_photon_energy(part_te)

    ! Is this photon to be added into the simulation?
    add_photon = (photon_energy > photon_energy_min_k_alpha &
        .AND. produce_k_alpha_photons)

    ! Ensure photon_energy is a number we can handle at our precision
    IF (photon_energy < c_tiny) photon_energy = c_tiny

    ! Temporarily store photon momentum in the "new_photon" particle
    CALL create_particle(new_photon)
    photon_p = photon_energy / c
    new_photon%part_p(1) = dir_x * photon_p
    new_photon%part_p(2) = dir_y * photon_p
    new_photon%part_p(3) = dir_z * photon_p

    ! Rotate K-alpha photon momentum to form an isotropic distribution
    IF (use_brem_scatter .AND. add_photon) THEN
      scatter_cos_theta = 1.0_num - 2.0_num * random()
      scatter_phi = 2.0_num * pi * random()
      CALL rotate_p(new_photon, scatter_cos_theta, scatter_phi, photon_p)
    END IF

    ! This will only create photons that have energies above a user specified
    ! cutoff and if photon generation is turned on.
    IF (add_photon) THEN
      ! Put new photon at electron position
      new_photon%part_pos = electron%part_pos
#ifdef PHOTONS
      new_photon%optical_depth = reset_optical_depth()
#endif
      new_photon%particle_energy = photon_energy
      new_photon%weight = electron%weight * k_alpha_weight

      ! Add photon to its species
      CALL add_particle_to_partlist(species_list(iphoton)%attached_list, &
          new_photon)

    ! We are not adding this photon to the simulation, so remove the particle
    ELSE
      CALL destroy_particle(new_photon)
    END IF

  END SUBROUTINE generate_photon



  FUNCTION calculate_photon_energy(part_te)

    ! Samples a photon energy according to the fitting function of Suxing. In
    ! this model, the photon emission scales with plasma temperature

    REAL(num), INTENT(IN) :: part_te
    REAL(num) :: te_ev, e_ev, sig_ev, calculate_photon_energy

    ! Convert temperature to [eV]
    te_ev = part_te * kb / q0

    ! Randomly choose either K-alpha1 or K-alpha2, sample energy. The centrak
    ! line-out comes from the Suxing model, the line width comes from a
    ! user-defined polynomial function
    IF (random() < 2.0_num / 3.0_num) THEN
      ! K-alpha1, 2/3 probability
      e_ev = 8047.8_num + 9.304e-3_num - 0.2596_num * te_ev &
          + 8.725e-3_num * te_ev**2 - 8.444e-5_num * te_ev**3 &
          + 3.908e-7_num * te_ev**4 - 6.785e-10_num * te_ev**5

      ! Calculate standard deviation of Gaussian energy distribution
      sig_ev = sig_ka1_0 + sig_ka1_1 * te_ev + sig_ka1_2 * te_ev**2 &
          + sig_ka1_3 * te_ev**3 + sig_ka1_4 * te_ev**4

      ! Shift central line energy by a normally-distributed random number
      e_ev = e_ev + random_box_muller(sig_ev)
    ELSE
      ! K-alpha2, 1/3 probability
      e_ev = 8027.3_num - 8.486e-3_num - 0.2648_num * te_ev &
          + 8.837e-3_num * te_ev**2 - 8.335e-5_num * te_ev**3 &
          + 3.684e-7_num * te_ev**4 - 6.042e-10_num * te_ev**5

      ! Calculate standard deviation of Gaussian energy distribution
      sig_ev = sig_ka2_0 + sig_ka2_1 * te_ev + sig_ka2_2 * te_ev**2 &
          + sig_ka2_3 * te_ev**3 + sig_ka2_4 * te_ev**4

      ! Shift central line energy by a normally-distributed random number
      e_ev = e_ev + random_box_muller(sig_ev)
    END IF

    ! Convert back to [J]
    calculate_photon_energy = e_ev * q0

  END FUNCTION calculate_photon_energy



  SUBROUTINE shutdown_k_alpha_module

    ! Shuts down the k_alpha module
    DEALLOCATE(rbeb_ke, rbeb_1s)

  END SUBROUTINE shutdown_k_alpha_module



  FUNCTION find_value_from_table_1d_coll(x_in, nx, x, values, state)

    ! For a pair of arrays, x and values, of size nx, this function returns the
    ! interpolated value of "values" corresponding to x_in in the x array. This
    ! uses linear interpolation, unlike in photons.F90 which is logarithmic

    REAL(num) :: find_value_from_table_1d_coll
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(:), values(:)
    TYPE(interpolation_state), INTENT(INOUT) :: state
    REAL(num) :: fx, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.
    LOGICAL :: found

    ! Scan through x to find correct row of table
    i1 = state%ix1
    i2 = state%ix2
    found = .FALSE.
    IF (i1 /= i2) THEN
      xdif1 = x(i1) - x_in
      xdif2 = x(i2) - x_in
      IF (xdif1 * xdif2 < 0) THEN
        found = .TRUE.
      ELSE
        i1 = MIN(state%ix1+1,nx)
        i2 = MIN(state%ix2+1,nx)
        xdif1 = x(i1) - x_in
        xdif2 = x(i2) - x_in
        IF (xdif1 * xdif2 < 0) THEN
          found = .TRUE.
        ELSE
          i1 = MAX(state%ix1-1,1)
          i2 = MAX(state%ix2-1,1)
          xdif1 = x(i1) - x_in
          xdif2 = x(i2) - x_in
          IF (xdif1 * xdif2 < 0) THEN
            found = .TRUE.
          END IF
        END IF
      END IF
    END IF

    IF (.NOT.found) THEN
      ! Scan through x to find correct row of table
      i1 = 1
      i2 = nx
      xdif1 = x(i1) - x_in
      xdif2 = x(i2) - x_in
      IF (xdif1 * xdif2 < 0) THEN
        ! Use bisection to find the nearest cell
        DO
          im = (i1 + i2) / 2
          xdifm = x(im) - x_in
          IF (xdif1 * xdifm < 0) THEN
            i2 = im
          ELSE
            i1 = im
            xdif1 = xdifm
          END IF
          IF (i2 - i1 == 1) EXIT
        END DO
        found = .TRUE.
      END IF
    END IF

    IF (found) THEN
      ! Interpolate in x to find fraction between the cells
      fx = (x_in - x(i1)) / (x(i2) - x(i1))
      state%ix1 = i1
      state%ix2 = i2
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_1d_coll" in ', &
            'collision_ionise.F90 outside the range of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      ! Our x_in value falls outside of the x array - truncate the value
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
      state%ix1 = 1
      state%ix2 = 1
    END IF

    ! Corresponding number from value array, a fraction fx between i1 and i2
    value_interp = (1.0_num - fx) * values(i1) + fx * values(i2)
    find_value_from_table_1d_coll = value_interp
    state%x = x_in
    state%val1d = find_value_from_table_1d_coll

  END FUNCTION find_value_from_table_1d_coll



  SUBROUTINE grid_centred_var_at_particle(part_x, part_var, grid_var)

    ! Calculates the value of a grid-centred variable part_var stored in the
    ! grid grid_var, averaged over the particle shape for a particle at position
    ! (part_x)

    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(IN) :: grid_var(1-ng:)
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

  END SUBROUTINE grid_centred_var_at_particle

#endif
END MODULE k_alpha
