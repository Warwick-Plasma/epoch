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
! This module contains scripts used to calculate the nuclear trident process,
! where the collision between two charged particles creates an electron/positron
! pair. The theory comes from Bhabha (1935), and implementation details are
! given in Section II.D of Martinez (2019).
!
! Bhabha, H. J. (1935). Proc. Math. Phys. Eng. Sci., 152(877), 559-586.
! Martinez, B., et al (2019). Phys. Plasmas, 26(10).

MODULE brem_trident
#ifdef BREM_TRIDENT

  USE partlist
  USE particles
  USE calc_df
  USE setup

  IMPLICIT NONE

  REAL(num), ALLOCATABLE :: incident_energy(:), pair_energy_list(:)
  REAL(num), ALLOCATABLE :: pair_energies(:,:), cdf_pair_energies(:,:)
  REAL(num), ALLOCATABLE :: energy_split_table(:,:), cdf_energy_split_table(:,:)

  INTEGER :: sample_el_in = 100
  INTEGER :: sample_pairs = 300
  REAL(num) :: max_incident_ke = 10.0e9_num * q0

  TYPE(interpolation_state), SAVE :: last_state_pair
  TYPE(interpolation_state), SAVE :: last_state_positron

CONTAINS

  SUBROUTINE setup_brem_trident_tables

    ! Create look-up tables to calculate pair energies for a range of incident 
    ! electron energies

    INTEGER :: i_samp
    REAL(num) :: min_ke, max_ke, el_ke

    ALLOCATE(incident_energy(sample_el_in))
    ALLOCATE(pair_energy_list(sample_el_in))
    ALLOCATE(pair_energies(sample_el_in, sample_pairs))
    ALLOCATE(cdf_pair_energies(sample_el_in, sample_pairs))
    ALLOCATE(energy_split_table(sample_el_in, sample_pairs))
    ALLOCATE(cdf_energy_split_table(sample_el_in, sample_pairs))

    ! Minimum KE to make rest pair at rest, maximum to match brem tables 10 GeV
    min_ke = 2.0_num * m0 * c**2
    max_ke = max_incident_ke

    ! Logarithmically sample total energies between these limits
    DO i_samp = 1, sample_el_in
      el_ke = min_ke * &
          (max_ke / min_ke)**(REAL(i_samp-1,num)/REAL(sample_el_in-1,num))
      incident_energy(i_samp) = el_ke + m0 * c**2
    END DO

    ! Second sampling list, corresponding to the maximum pair energy from each
    ! incident_energy
    pair_energy_list = incident_energy - m0c2

    ! Populate lookup tables for pair creation
    CALL calculate_pair_energy_tables
    CALL calculate_energy_split_table_tables

  END SUBROUTINE setup_brem_trident_tables



  SUBROUTINE calculate_pair_energy_tables

    ! Sets the value of the 2D arrays pair_energies and cdf_pair_energies. For 
    ! each incident electron energy, pair_energies provides a line of possible 
    ! energies a trident-produced pair could have. For each of these 
    ! incident-pair energy combinations, there is a corresponding CDF value in 
    ! cdf_pair_energies, read from the differential cross section in Martinez 
    ! (2019) equation (30).

    INTEGER :: i_el, i_pair
    REAL(num) :: min_ke_pair, max_ke_pair, ke_pair 
    REAL(num) :: gamma_in, gamma_p
    REAL(num) :: mart_x, mart_x2_frac, mart_x2_log
    REAL(num) :: mart_c, mart_cz, mart_cr, mart_c_terms
    REAL(num) :: dcs_nonrel, dcs_rel, dcs, dcs_sum

    ! Loop over incident electron energies
    DO i_el = 1, sample_el_in
      min_ke_pair = 0.0_num
      max_ke_pair = (incident_energy(i_el) - m0 * c**2) - (2.0_num*m0*c**2)

      ! Pre-calculate terms for dcs calculation
      gamma_in = incident_energy(i_el) / (m0*c**2)
      dcs_sum = 0.0_num

      ! Martinez (2019) eq. (29)
      mart_x = 1.0_num / gamma_in
      mart_x2_frac = mart_x**2 / (1.0_num - mart_x**2)
      mart_x2_log = mart_x2_frac * LOG(1.0_num / mart_x**2)
      mart_c = 4.0_num * mart_x2_log - (4.0_num/3.0_num)*mart_x**2 &
          + (1.0_num/6.0_num)*mart_x**4
      mart_cz = 3.0_num * mart_x2_frac * (1.0_num - mart_x2_log) &
          - 2.6_num*mart_x**2 + 1.75_num*mart_x**4 - 0.9_num*mart_x**6 &
          + 0.2_num*mart_x**8
      mart_cr = -1.5_num * mart_x2_frac * (1.0_num - mart_x2_log) &
          + 0.8_num*mart_x**2 - 0.125_num*mart_x**4 - 0.05_num*mart_x**6 &
          + 0.025_num*mart_x**8
      mart_c_terms = (-161.0_num/60.0_num + mart_c + mart_cz + mart_cr)

      DO i_pair = 1, sample_pairs 
        ! Uniformly sample pair energies between the two limits
        ke_pair = min_ke_pair + (max_ke_pair - min_ke_pair) &
        * (REAL(i_pair-1,num) / REAL(sample_pairs-1,num)) 
        gamma_p = ke_pair / (m0*c**2) + 2.0_num
        pair_energies(i_el, i_pair) = gamma_p * m0 * c**2 

        ! Save summed differential cross section contribution at this gamma
        ! Martinez (2019) eq. (27) - (ignore (Z*re*alpha) shared prefactor)
        dcs_nonrel = (1.0_num/32.0_num) &
            * (LOG(gamma_in**2) + mart_c_terms) * (gamma_p - 2.0_num)**3
        ! Martinez (2019) eq. (28) - (ignore (Z*re*alpha) shared prefactor)
        dcs_rel = (56.0_num/(9.0_num*pi)) * LOG(gamma_p) &
            * LOG(gamma_in / gamma_p) / gamma_p
        ! Martinez (2019) eq. (30)
        dcs = dcs_nonrel * dcs_rel / (dcs_nonrel + dcs_rel)
        dcs_sum = dcs_sum + dcs

        ! Temporarily store the summed DCS in cdf_pair_energies
        cdf_pair_energies(i_el, i_pair) = dcs_sum
      END DO 

      ! Normalise cdf_pair_energies array to actually give CDF values 
      IF (dcs_sum > 0.0_num) THEN  
        cdf_pair_energies(i_el,:) = cdf_pair_energies(i_el,:) / dcs_sum
      ELSE
        ! Invalid dcs for all electron energies - no pairs at incident energy
        cdf_pair_energies(i_el,1) = 0.0_num
        cdf_pair_energies(i_el,2:sample_pairs) = 1.0_num
      END IF
    END DO

  END SUBROUTINE calculate_pair_energy_tables



  SUBROUTINE calculate_energy_split_table_tables

    ! Sets the value of the 2D array cdf_energy_split_table. For each possible 
    ! pair energy (electron + positron total energies), energy_split_table 
    ! provides a line of possible electron energies. The CDF of a pair-electron 
    ! having a particular energy is given in cdf_energy_split_table, where array 
    ! elements correspond to the matching pair/electron energy combo in 
    ! energy_split_table   

    INTEGER :: i_pair, i_el
    REAL(num) :: pair_energy
    REAL(num) :: min_ke_el, max_ke_el, ke_el
    REAL(num) :: gamma_pair, gamma_el, gamma_pos
    REAL(num) :: dcs, dcs_sum

    ! Loop over possible pair energies, using incident_energy as sample points
    DO i_pair = 1, sample_el_in
      pair_energy = pair_energy_list(i_pair)
      gamma_pair = pair_energy / (m0*c**2)

      min_ke_el = 0.0_num
      max_ke_el = pair_energy - 2.0_num*m0*c**2

      dcs_sum = 0.0_num
      DO i_el = 1, sample_pairs 
        ! Uniformly sample pair-electron energies between the two limits
        ke_el = min_ke_el + (max_ke_el - min_ke_el) &
            * (REAL(i_el-1,num) / REAL(sample_pairs-1,num))
        gamma_el = ke_el / (m0*c**2) + 1.0_num
        energy_split_table(i_pair, i_el) = gamma_el * m0 * c**2 

        ! Save summed differential cross section contribution at this electron
        ! energy
        gamma_pos = gamma_pair - gamma_el
        ! Martinez (2019) eq. (31)
        dcs = (gamma_pos**2 + gamma_el**2 &
            + 2.0_num*gamma_el*gamma_pos/3.0_num) &
            * LOG(gamma_el * gamma_pos / gamma_pair)
        dcs_sum = dcs_sum + MAX(dcs,0.0_num)

        ! Temporarily store the summed DCS in cdf_energy_split_table
        cdf_energy_split_table(i_pair, i_el) = dcs_sum
      END DO

      ! Normalise cdf_pair_energies array to actually give CDF values
      IF (dcs_sum > 0.0_num) THEN  
        cdf_energy_split_table(i_pair,:) = cdf_energy_split_table(i_pair,:) &
            / dcs_sum
      ELSE
        ! Invalid dcs for all electron energies - forbidden pair energy value
        cdf_energy_split_table(i_pair,1:sample_pairs/2) = 0.0_num
        cdf_energy_split_table(i_pair,sample_pairs/2+1:sample_pairs) = 1.0_num
      END IF
    END DO

  END SUBROUTINE calculate_energy_split_table_tables



  SUBROUTINE brem_trident_update_depth

    ! Cycles through all electron species and updates the nuclear trident 
    ! optical depths

    INTEGER :: iz, z_temp, ispecies
    REAL(num), ALLOCATABLE :: grid_num_density_ion(:,:,:)
    TYPE(particle), POINTER :: electron, next_electron
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: el_p2, el_energy, el_gamma, el_p, el_v
    REAL(num) :: cross_sec_fac, cross_sec, part_ni, delta_opdep

    ALLOCATE(grid_num_density_ion(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))

    ! Calculate the number density of each ion species
    DO iz = 1, n_species
      ! Identify if the charge is greater than 1
      z_temp = species_list(iz)%atomic_no
      IF (z_temp < 1 .OR. z_temp > 100) CYCLE

      cross_sec_fac = 5.22e-34_num * z_temp**2
      CALL calc_number_density(grid_num_density_ion, iz)
      CALL field_bc(grid_num_density_ion, ng)

      ! Update the optical depth for each electron species
      DO ispecies = 1, n_species
        ! Only update optical_depth_brem_tri for electrons
        IF (species_list(ispecies)%species_type == c_species_id_electron) THEN
          ! Cycle through all electrons in this species
          electron => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(electron))
            next_electron => electron%next

            el_p2 = electron%part_p(1)**2 + electron%part_p(2)**2 &
                + electron%part_p(3)**2
            el_energy = SQRT(el_p2*c**2 + m0c2**2)

            ! Ignore electrons too slow to undergo nuclear trident pair 
            ! production
            IF (el_energy < 3.0_num * m0c2) THEN
              electron => next_electron
              CYCLE
            END IF

            el_gamma = el_energy / m0c2 
            el_p = SQRT(el_p2)
            el_v = el_p / (el_gamma * m0)
            
            ! Nuclear trident cross section
            ! Martinez (2019) eq. (26)
            cross_sec = cross_sec_fac * LOG((el_gamma + 4.5_num)/6.89_num)**3

            ! Get background number density at incident electron position
            part_x = electron%part_pos(1) - x_grid_min_local
            part_y = electron%part_pos(2) - y_grid_min_local
            part_z = electron%part_pos(3) - z_grid_min_local
            CALL grid_centred_var_at_particle_bt(part_x, part_y, part_z, &
                part_ni, grid_num_density_ion)

            ! Update the incident electron optical depth
            delta_opdep = part_ni * cross_sec * el_v * dt
            electron%optical_depth_brem_tri = electron%optical_depth_brem_tri &
                - delta_opdep

            ! If electron optical depth drops below 0, create an e-/e+ pair,
            ! recoil the incident electron, and sample a new optical depth
            IF (electron%optical_depth_brem_tri <= 0.0_num) THEN
              CALL generate_pair(electron, el_energy, el_p, &
                  brem_trident_electron_species, brem_trident_positron_species)
                  electron%optical_depth_brem_tri = &
                      LOG(1.0_num / (1.0_num - random()))
            END IF

            electron => next_electron
          END DO
        END IF
      END DO
    END DO

    DEALLOCATE(grid_num_density_ion)

  END SUBROUTINE brem_trident_update_depth



  SUBROUTINE generate_pair(incident, el_in_energy, el_in_p, ielectron, &
      ipositron)

    ! This subroutine is called when an incident electron "incident" undergoes 
    ! nuclear trident pair production. Macro-particles for e- and e+ are 
    ! created, and the energies and directions are sampled using Martinez 
    ! (2019) equations

    TYPE(particle), POINTER :: incident
    INTEGER, INTENT(IN) :: ielectron, ipositron
    REAL(num), INTENT(IN) :: el_in_energy, el_in_p
    TYPE(particle), POINTER :: new_electron, new_positron
    REAL(num) :: rand_cdf, pair_energy, incident_p
    REAL(num) :: electron_e, electron_p
    REAL(num) :: positron_e, positron_p
    REAl(num) :: norm, dir_x, dir_y, dir_z

    ! Create pair at incident electron position
    CALL create_particle(new_electron)
    CALL create_particle(new_positron)
    new_electron%part_pos = incident%part_pos
    new_positron%part_pos = incident%part_pos

    ! e- and e+ have the same weights as generating incident e-
    new_electron%weight = incident%weight
    new_positron%weight = incident%weight

    ! Calculate the total energy of the pair
    rand_cdf = random()
    pair_energy = find_value_from_table_2d_bt(el_in_energy, rand_cdf, &
        sample_el_in, sample_pairs, &
        incident_energy, pair_energies, cdf_pair_energies, last_state_pair)

    ! Calculate the e+ energy
    rand_cdf = random()
    positron_e = find_value_from_table_2d_bt(pair_energy, rand_cdf, &
        sample_el_in, sample_pairs, pair_energy_list, energy_split_table, &
        cdf_energy_split_table, last_state_positron) 
    electron_e = pair_energy - positron_e

    ! Calculate momentum magnitude going to each particle
    electron_p = SQRT(electron_e**2 - m0c2**2)/c
    positron_p = SQRT(positron_e**2 - m0c2**2)/c

    ! Assume pair emitted in incident electron direction
    ! Proper treatment could be derived from Bhabha (27)
    norm  = 1.0_num / el_in_p
    dir_x = incident%part_p(1) * norm
    dir_y = incident%part_p(2) * norm
    dir_z = incident%part_p(3) * norm
    new_electron%part_p(1:3) = electron_p * (/ dir_x, dir_y, dir_z /)
    new_positron%part_p(1:3) = positron_p * (/ dir_x, dir_y, dir_z /)

    ! Set particle energies
    new_electron%particle_energy = electron_e
    new_positron%particle_energy = positron_e

    ! Save pair to particle list
    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

    ! Recoil incident electron
    incident_p = el_in_p - (electron_p + positron_p)
    incident%part_p(1:3) = incident_p * (/ dir_x, dir_y, dir_z /)

  END SUBROUTINE generate_pair



  FUNCTION find_value_from_table_1d_bt(x_in, nx, x, values, state)

    ! For a pair of arrays, x and values, of size nx, this function returns the
    ! interpolated value of "values" corresponding to x_in in the x array. This
    ! uses linear interpolation, unlike in photons.F90 which is logarithmic

    REAL(num) :: find_value_from_table_1d_bt
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
        PRINT*,'Argument to "find_value_from_table_1d_bt" in ', &
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
    find_value_from_table_1d_bt = value_interp
    state%x = x_in
    state%val1d = find_value_from_table_1d_bt

  END FUNCTION find_value_from_table_1d_bt



  FUNCTION find_value_from_table_2d_bt(x_in, p_value, nx, ny, x, y, p_table, &
      state)

    ! For each element of x, we have a 1D array of y values and a 1D array of P
    ! values, such that the 1D array x has corresponding 2D arrays y and
    ! p_table. The 2D arrays y and p_table are of equal size (nx,ny). This
    ! function interpolates in x_in first, creating an interpolated 1D array of
    ! y and p_table values. The second interpolation finds p_value in p_table,
    ! and the function returns the corresponding value in the interpolated 1D y
    ! array. Used for CDF look-ups.

    REAL(num) :: find_value_from_table_2d_bt
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(:), y(:,:), p_table(:,:)
    TYPE(interpolation_state), INTENT(INOUT) :: state
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, y_interp, xdif1, xdif2, xdifm
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
      ! Interpolate in x
      fx = (x_in - x(i1)) / (x(i2) - x(i1))
      state%ix1 = i1
      state%ix2 = i2
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_2d_bt" outside the ', &
            'range of the table.'
        PRINT*,'An incident electron kinetic energy exceeds tabulated values'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
      state%ix1 = 1
      state%ix2 = 1
    END IF

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = state%iy1
    i2 = state%iy2
    found = .FALSE.
    IF (i1 /= i2) THEN
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        found = .TRUE.
      ELSE
        i1 = MIN(state%iy1+1,ny)
        i2 = MIN(state%iy2+1,ny)
        xdif1 = p_table(ix,i1) - p_value
        xdif2 = p_table(ix,i2) - p_value
        IF (xdif1 * xdif2 < 0) THEN
          found = .TRUE.
        ELSE
          i1 = MAX(state%iy1-1,1)
          i2 = MAX(state%iy2-1,1)
          xdif1 = p_table(ix,i1) - p_value
          xdif2 = p_table(ix,i2) - p_value
          IF (xdif1 * xdif2 < 0) THEN
            found = .TRUE.
          END IF
        END IF
      END IF
    END IF

    IF (.NOT.found) THEN
      ! Scan through table row to find p_value
      i1 = 1
      i2 = ny
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        ! Use bisection to find the nearest cell
        DO
          im = (i1 + i2) / 2
          xdifm = p_table(ix,im) - p_value
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
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
      state%iy1 = i1
      state%iy2 = i2
    ELSE
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
      state%iy1 = 1
      state%iy2 = 1
    END IF

    y_lt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = state%iy1
    i2 = state%iy2
    found = .FALSE.
    IF (i1 /= i2) THEN
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        found = .TRUE.
      ELSE
        i1 = MIN(state%iy1+1,ny)
        i2 = MIN(state%iy2+1,ny)
        xdif1 = p_table(ix,i1) - p_value
        xdif2 = p_table(ix,i2) - p_value
        IF (xdif1 * xdif2 < 0) THEN
          found = .TRUE.
        ELSE
          i1 = MAX(state%iy1-1,1)
          i2 = MAX(state%iy2-1,1)
          xdif1 = p_table(ix,i1) - p_value
          xdif2 = p_table(ix,i2) - p_value
          IF (xdif1 * xdif2 < 0) THEN
            found = .TRUE.
          END IF
        END IF
      END IF
    END IF

    IF (.NOT.found) THEN
      ! Scan through table row to find p_value
      i1 = 1
      i2 = ny
      xdif1 = p_table(ix,i1) - p_value
      xdif2 = p_table(ix,i2) - p_value
      IF (xdif1 * xdif2 < 0) THEN
        ! Use bisection to find the nearest cell
        DO
          im = (i1 + i2) / 2
          xdifm = p_table(ix,im) - p_value
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
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
      state%iy1 = i1
      state%iy2 = i2
    ELSE
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
      state%iy1 = 1
      state%iy2 = 1
    END IF

    y_gt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ! Interpolate in x
    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table_2d_bt = y_interp
    state%x = x_in
    state%y = p_value
    state%val2d = find_value_from_table_2d_bt

  END FUNCTION find_value_from_table_2d_bt



  ! Calculates the value of a grid-centred variable part_var stored in the grid
  ! grid_var, averaged over the particle shape for a particle at position
  ! (part_x, part_y, part_z)

  SUBROUTINE grid_centred_var_at_particle_bt(part_x, part_y, part_z, part_var, &
      grid_var)

    REAL(num), INTENT(IN) :: part_x, part_y, part_z
    REAL(num), INTENT(IN) :: grid_var(1-ng:,1-ng:,1-ng:)
    REAL(num), INTENT(OUT) :: part_var
    INTEGER :: cell_x1, cell_y1, cell_z1
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
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
    cell_y_r = part_y / dy - 0.5_num
    cell_z_r = part_z / dz - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
    cell_z_r = part_z / dz
#endif
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1
    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1
    cell_z1 = FLOOR(cell_z_r + 0.5_num)
    cell_frac_z = REAL(cell_z1, num) - cell_z_r
    cell_z1 = cell_z1 + 1

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

  END SUBROUTINE grid_centred_var_at_particle_bt

#endif
END MODULE brem_trident
