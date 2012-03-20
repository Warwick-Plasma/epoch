#ifdef PHOTONS
MODULE photons

  USE shared_data
  USE random_generator
  USE partlist
  USE shape_functions

  IMPLICIT NONE

  SAVE
  REAL(num) :: part_pos_global, gamma_global, eta_global

CONTAINS

  SUBROUTINE setup_qed_module

    INTEGER :: ispecies

    ! Load the tables for the QED routines
    CALL setup_tables_qed
    DO ispecies = 1, n_species
      CALL initialise_optical_depth(species_list(ispecies))
    ENDDO

  END SUBROUTINE setup_qed_module



  SUBROUTINE shutdown_qed_module

    CALL deallocate_tables_qed

  END SUBROUTINE shutdown_qed_module



  FUNCTION check_qed_variables()

    INTEGER :: check_qed_variables
    INTEGER :: io, ispecies
    INTEGER :: first_electron = -1, first_positron = -1

    check_qed_variables = c_err_none

    IF (.NOT.use_qed) RETURN

    ! If you're only doing radiation reaction force then don't need any special
    ! species, so don't do any checking here
    IF (.NOT.produce_photons) RETURN

    ! Identify if there exists any electron species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type .EQ. c_species_id_electron &
          .AND. first_electron .EQ. -1) THEN
        first_electron = ispecies
      ENDIF
      IF (species_list(ispecies)%species_type .EQ. c_species_id_positron &
          .AND. first_positron .EQ. -1) THEN
        first_positron = ispecies
      ENDIF
    ENDDO

    IF (first_electron .LT. 0 .AND. first_positron .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No electron or positron species specified.'
          WRITE(io,*) 'Specify using "identify:electron" or "identify:positron"'
          WRITE(io,*) 'QED routines require at least one species of ', &
              'electrons or positrons.'
        ENDDO
      ENDIF
      check_qed_variables = c_err_missing_elements
      RETURN
    ENDIF

    IF (photon_species .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No photon species specified. Specify using ', &
              '"identify:photon"'
        ENDDO
      ENDIF
      check_qed_variables = c_err_missing_elements
      RETURN
    ENDIF

    ! If you're not producing pairs then you don't have to designate special
    ! electron or positron species so just return
    IF (.NOT.produce_pairs) RETURN

    IF (first_electron .LT. 0 .OR. first_positron .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'To use pair production routines need at least one ', &
              'positron species and one ', 'electron species. Specify ', &
              'using "identify:electron" or "identify:positron"'
        ENDDO
      ENDIF
      check_qed_variables = c_err_missing_elements
      RETURN
    ENDIF

    IF (breit_wheeler_positron_species .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No Breit-Wheeler positron species specified.'
          WRITE(io,*) 'Specify using "identify:bw_positron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_positron)%name), ' instead.'
        ENDDO
      ENDIF
      breit_wheeler_positron_species = first_positron
    ENDIF

#ifdef TRIDENT_PHOTONS
    IF (trident_positron_species .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No trident positron species specified.'
          WRITE(io,*) 'Specify using "identify:trident_positron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_positron)%name), ' instead.'
        ENDDO
      ENDIF
      trident_electron_species = first_positron
    ENDIF
#endif

    IF (breit_wheeler_electron_species .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No Breit-Wheeler electron species specified.'
          WRITE(io,*) 'Specify using "identify:bw_electron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_electron)%name), ' instead.'
        ENDDO
      ENDIF
      breit_wheeler_electron_species = first_electron
    ENDIF

#ifdef TRIDENT_PHOTONS
    IF (trident_electron_species .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'No trident electron species specified.'
          WRITE(io,*) 'Specify using "identify:trident_electron".'
          WRITE(io,*) 'Using species ', &
              TRIM(species_list(first_electron)%name), ' instead.'
        ENDDO
      ENDIF
      trident_electron_species = first_electron
    ENDIF
#endif
  END FUNCTION check_qed_variables



  SUBROUTINE setup_tables_qed

    ! Reads files epsilon.table, log_chi.table, energy_split.table
    ! and sets up appropriate tables

    REAL(num) :: etalog_min, etalog_max
    INTEGER :: ichi2, iepsilon, ieta, ichi, idummy

    OPEN(50, file='./src/physics_packages/TABLES/hsokolov.table', status='OLD')
    READ(50,*) n_sample_h
    ALLOCATE(log_hsokolov(2,n_sample_h))
    DO ieta = 1, n_sample_h
      READ(50,*) log_hsokolov(1:2,ieta)
    ENDDO
    CLOSE(50)

    OPEN(51, file='./src/physics_packages/TABLES/pairprod.table', status='OLD')
    READ(51,*) n_sample_t
    ALLOCATE(log_tpair_dummy(3,n_sample_t))
    ALLOCATE(log_tpair(2,n_sample_t))
    ALLOCATE(log_omegahat(2,n_sample_t))
    DO ichi = 1, n_sample_t
      READ(51,*) log_tpair_dummy(1:3,ichi)
    ENDDO
    CLOSE(51)

    DO ichi = 1, n_sample_t
      log_tpair(1,ichi) = log_tpair_dummy(1,ichi)
      log_tpair(2,ichi) = log_tpair_dummy(3,ichi)
      log_omegahat(1,ichi) = log_tpair_dummy(1,ichi)
      log_omegahat(2,ichi) = log_tpair_dummy(2,ichi)
    ENDDO

    OPEN(70, file='./src/physics_packages/TABLES/ksi_sokolov.table', &
        status='OLD')
    OPEN(80, file='./src/physics_packages/TABLES/chimin.table', status='OLD')

    READ(70,*) n_sample_eta, n_sample_chi, etalog_min, etalog_max

    ALLOCATE(chimin_table(n_sample_eta))
    ALLOCATE(log_eta(n_sample_eta))
    ALLOCATE(log_chi(n_sample_eta,n_sample_chi))
    ALLOCATE(p_photon_energy(n_sample_eta,n_sample_chi))

    READ(80,*) (chimin_table(ieta), ieta=1,n_sample_eta)

    DO ieta = 1, n_sample_eta
      log_eta(ieta) = etalog_min + (DBLE(ieta) - 1.0_num) &
          * (etalog_max - etalog_min) / (DBLE(n_sample_eta) - 1.0_num)
    ENDDO

    DO ieta = 1, n_sample_eta
      DO ichi = 1, n_sample_chi
        log_chi(ieta,ichi) = LOG10(chimin_table(ieta)) &
            + (DBLE(ichi) - 1.0_num) &
            * (LOG10(10.0_num**(log_eta(ieta)) / 2.0_num) &
            - LOG10(chimin_table(ieta))) / (DBLE(n_sample_chi) - 1.0_num)
      ENDDO
    ENDDO

    DO ieta = 1, n_sample_eta
      DO idummy = 1, ieta
        READ(70,*)
      ENDDO
      READ(70,*) (p_photon_energy(ieta,ichi), ichi=1,n_sample_chi)
      REWIND(70)
    ENDDO

    CLOSE(70)
    CLOSE(80)

    OPEN(10, file='./src/physics_packages/TABLES/log_chi2.table', status='OLD')
    READ(10,*) n_sample_chi2
    ALLOCATE(log_chi2(n_sample_chi2))
    DO ichi2 = 1, n_sample_chi2
      READ(10,*) log_chi2(ichi2)
    ENDDO
    CLOSE(10)

    OPEN(21, file='./src/physics_packages/TABLES/epsilon.table', status='OLD')
    READ(21,*) n_sample_epsilon
    ALLOCATE(epsilon_split(n_sample_epsilon))
    DO iepsilon = 1, n_sample_epsilon
      READ(21,*) epsilon_split(iepsilon)
    ENDDO
    CLOSE(21)

    ALLOCATE(p_energy(n_sample_chi2,n_sample_epsilon))
    OPEN(30, file='./src/physics_packages/TABLES/energy_split.table', &
        status='OLD')
    DO ichi2 = 1, n_sample_chi2
      DO iepsilon = 1, n_sample_epsilon
        READ(30,*) p_energy(ichi2,iepsilon)
      ENDDO
    ENDDO
    CLOSE(30)

  END SUBROUTINE setup_tables_qed



  SUBROUTINE deallocate_tables_qed

    DEALLOCATE(log_chi2)
    DEALLOCATE(epsilon_split)
    DEALLOCATE(p_energy)
    DEALLOCATE(log_hsokolov)
    DEALLOCATE(log_eta)
    DEALLOCATE(log_chi)
    DEALLOCATE(p_photon_energy)
    DEALLOCATE(log_tpair_dummy)
    DEALLOCATE(log_tpair)
    DEALLOCATE(log_omegahat)
    DEALLOCATE(chimin_table)

  END SUBROUTINE deallocate_tables_qed



  SUBROUTINE initialise_optical_depth(current_species)

    ! Resets optical depth (to random number) of all particles
    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    current => current_species%attached_list%head
    DO WHILE(ASSOCIATED(current))
      p_tau = random()
      current%optical_depth = LOG(1.0_num / (1.0_num - p_tau))

#ifdef TRIDENT_PHOTONS
      p_tau = random()
      current%optical_depth_tri = LOG(1.0_num / (1.0_num - p_tau))
#endif
      current => current%next
    ENDDO

  END SUBROUTINE initialise_optical_depth



  FUNCTION reset_optical_depth()

    ! Resets optical depth of particle
    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = LOG(1.0_num / (1.0_num - p_tau))

  END FUNCTION reset_optical_depth



  SUBROUTINE qed_update_optical_depth

    ! Updates the optical depth for electrons and photons
    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next_pt

    REAL(num) :: part_x
    REAL(num) :: part_px, part_py, part_pz
    REAL(num) :: part_vx, part_vy, part_vz
    REAL(num) :: eta, chi_val, part_e

    DO ispecies = 1, n_species

      ! First consider electrons
      IF (species_list(ispecies)%species_type .EQ. c_species_id_electron) THEN
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          ! Find eta at particle position
          part_x  = current%part_pos - x_min_local
          part_px = current%part_p(1)
          part_py = current%part_p(2)
          part_pz = current%part_p(3)

          eta = calculate_eta(part_x, part_px, part_py, part_pz)
          current%optical_depth = current%optical_depth &
              - delta_optical_depth(eta, part_px, part_py, part_pz)
#ifdef TRIDENT_PHOTONS
          current%optical_depth_tri = current%optical_depth_tri &
              - delta_optical_depth_tri(eta, part_px, part_py, part_pz)
#endif
          ! If optical depth dropped below zero generate photon...
          IF (current%optical_depth .LE. 0.0_num) THEN
            CALL generate_photon(current, photon_species, eta)
            ! ... and reset optical depth
            current%optical_depth = reset_optical_depth()
          ENDIF

#ifdef TRIDENT_PHOTONS
          IF (current%optical_depth_tri .LE. 0.0_num) THEN
            CALL generate_pair_tri(current, trident_electron_species, &
                trident_positron_species)
            ! ... and reset optical depth
            current%optical_depth_tri = reset_optical_depth()
          ENDIF
#endif
          current => current%next
        ENDDO
      ENDIF

      ! Next consider positrons
      IF (species_list(ispecies)%species_type .EQ. c_species_id_positron) THEN
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          ! Find eta at particle position
          part_x  = current%part_pos - x_min_local
          part_px = current%part_p(1)
          part_py = current%part_p(2)
          part_pz = current%part_p(3)

          eta = calculate_eta(part_x, part_px, part_py, part_pz)

          current%optical_depth = current%optical_depth &
              - delta_optical_depth(eta, part_px, part_py, part_pz)
#ifdef TRIDENT_PHOTONS
          current%optical_depth_tri = current%optical_depth_tri &
              - delta_optical_depth_tri(eta, part_px, part_py, part_pz)
#endif
          ! If optical depth dropped below zero generate photon...
          IF (current%optical_depth .LE. 0.0_num) THEN
            CALL generate_photon(current, photon_species, eta)
            ! ... and reset optical depth
            current%optical_depth = reset_optical_depth()
          ENDIF

#ifdef TRIDENT_PHOTONS
          IF (current%optical_depth_tri .LE. 0.0_num) THEN
            CALL generate_pair_tri(current, trident_electron_species, &
                trident_positron_species)
            ! ... and reset optical depth
            current%optical_depth_tri = reset_optical_depth()
          ENDIF
#endif
          current => current%next
        ENDDO
      ENDIF

      ! and finally photons
      IF (species_list(ispecies)%species_type .EQ. c_species_id_photon) THEN
        IF (produce_pairs) THEN
          current => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(current))
            ! Current may be deleted
            next_pt => current%next
            part_x  = current%part_pos - x_min_local
            part_vx = current%part_p(1)
            part_vy = current%part_p(2)
            part_vz = current%part_p(3)
            part_e = current%particle_energy
            chi_val = calculate_chi(part_x, part_vx, part_vy, &
                part_vz, part_e)

            current%optical_depth = current%optical_depth &
                - delta_optical_depth_photon(chi_val, part_e)
            ! If optical depth dropped below zero generate pair...
            IF (current%optical_depth .LE. 0.0_num) THEN
              CALL generate_pair(current, chi_val, photon_species, &
                  breit_wheeler_electron_species, &
                  breit_wheeler_positron_species)
            ENDIF
            current => next_pt
          ENDDO
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE qed_update_optical_depth



  FUNCTION delta_optical_depth(eta, part_px, part_py, part_pz)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth
    REAL(num), INTENT(IN) :: eta, part_px, part_py, part_pz
    REAL(num) :: hsokolov, gamma

    gamma = SQRT(1.0_num + (part_px**2 + part_py**2 + part_pz**2)/((m0*c)**2))

    hsokolov = find_value_from_table_1d(n_sample_h, log_hsokolov(1,:), &
        eta, log_hsokolov(2,:))

    delta_optical_depth = dt * eta * alpha_f * SQRT(3.0_num) * hsokolov &
        / (2.0_num * pi * tau_c * gamma)

  END FUNCTION delta_optical_depth



  FUNCTION delta_optical_depth_tri(eta, part_px, part_py, part_pz)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_tri
    REAL(num), INTENT(IN) :: eta, part_px, part_py, part_pz
    REAL(num) :: omegahat, gamma

    gamma = SQRT(1.0_num + (part_px**2 + part_py**2 + part_pz**2)/((m0*c)**2))

    omegahat = find_value_from_table_1d(n_sample_t, log_omegahat(1,:), &
        eta, log_omegahat(2,:))

    delta_optical_depth_tri = dt * eta * alpha_f**2 * 0.64_num * omegahat &
        / (2.0_num * pi * tau_c * gamma)

  END FUNCTION delta_optical_depth_tri



  FUNCTION delta_optical_depth_photon(chi_val, part_e)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_photon
    REAL(num), INTENT(IN) :: chi_val, part_e
    REAL(num) :: tpair

    tpair = find_value_from_table_1d(n_sample_t, log_tpair(1,:), &
        chi_val, log_tpair(2,:))

    delta_optical_depth_photon = (dt / tau_c) * alpha_f * (m0 * c**2 / part_e) &
        * chi_val * tpair

  END FUNCTION delta_optical_depth_photon



  FUNCTION calculate_eta(part_x, part_px, part_py, part_pz)

    REAL(num) :: calculate_eta
    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(IN) :: part_px, part_py, part_pz
    REAL(num) :: e_at_part(3), b_at_part(3)
    REAL(num) :: gamma, beta_x, beta_y, beta_z
    REAL(num) :: dir_x, dir_y, dir_z, flperp(3), i_e, tau0, roland_eta
    REAL(num) :: lambdac, coeff_eta, modp, modpclip, p_dot_e

    CALL field_at_particle(part_x, e_at_part, b_at_part)

    modp = part_px**2 + part_py**2 + part_pz**2
    modpclip = MAX(modp, c_non_zero)
    modp = SQRT(modp)
    modpclip = SQRT(modpclip)
    gamma = SQRT(1.0_num + modp**2 / (m0*c)**2)

    beta_x = part_px / (gamma * m0 * c)
    beta_y = part_py / (gamma * m0 * c)
    beta_z = part_pz / (gamma * m0 * c)
    dir_x = part_px / modpclip
    dir_y = part_py / modpclip
    dir_z = part_pz / modpclip

    lambdac = h_bar / m0 / c

    coeff_eta = SQRT(3.0_num * lambdac / (2.0_num * alpha_f * m0 * c**3))

    p_dot_e = (part_px * e_at_part(1) + part_py * e_at_part(2) &
        + part_pz * e_at_part(3)) / modpclip**2

    flperp(1) = q0 * (e_at_part(1) - p_dot_e * part_px &
        + c * (beta_y * b_at_part(3) - beta_z * b_at_part(2)))

    flperp(2) = q0 * (e_at_part(2) - p_dot_e * part_py &
        + c * (beta_z * b_at_part(1) - beta_x * b_at_part(3)))

    flperp(3) = q0 * (e_at_part(3) - p_dot_e * part_pz &
        + c * (beta_x * b_at_part(2) - beta_y * b_at_part(1)))

    ! Dipole emission intensity

    tau0 = q0**2 / (6.0_num * pi * epsilon0 * m0 * c**3)

    i_e = tau0 * gamma**2 * (flperp(1)**2 + flperp(2)**2 + flperp(3)**2 &
        + (q0 * m0 * c * (beta_x * e_at_part(1) + beta_y * e_at_part(2) &
        + beta_z * e_at_part(3)) / modpclip)**2) / m0

    roland_eta = coeff_eta * SQRT(i_e)

    ! Determine eta from fields
    calculate_eta = roland_eta

  END FUNCTION calculate_eta



  FUNCTION calculate_chi(part_x, part_vx, part_vy, part_vz, &
      part_e)

    REAL(num) :: calculate_chi
    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(IN) :: part_vx, part_vy, part_vz, part_e
    REAL(num) :: e_at_part(3), b_at_part(3), q(3)
    REAL(num) :: part_px, part_py, part_pz, part_e_n, dir_x, dir_y, dir_z
    REAL(num) :: e_dot_dir, calculate_chi_roland

    CALL field_at_particle(part_x, e_at_part, b_at_part)

    part_px = ((part_vx / c) * (part_e / c)) / (m0 * c)
    part_py = ((part_vy / c) * (part_e / c)) / (m0 * c)
    part_pz = ((part_vz / c) * (part_e / c)) / (m0 * c)
    part_e_n = part_e / (m0 * c**2)

    dir_x = part_vx / c
    dir_y = part_vy / c
    dir_z = part_vz / c

    e_dot_dir = e_at_part(1) * dir_x + e_at_part(2) * dir_y &
        + e_at_part(3) * dir_z

    q(1) = e_at_part(1) - e_dot_dir * dir_x &
        + c * (dir_y * b_at_part(3) - dir_z * b_at_part(2))
    q(2) = e_at_part(2) - e_dot_dir * dir_y &
        + c * (dir_z * b_at_part(1) - dir_x * b_at_part(3))
    q(3) = e_at_part(3) - e_dot_dir * dir_z &
        + c * (dir_x * b_at_part(2) - dir_y * b_at_part(1))

    calculate_chi_roland = 0.5_num * SQRT(q(1)**2 + q(2)**2 + q(3)**2) &
        * part_e_n / e_s

    ! Determine chi from fields
    calculate_chi = calculate_chi_roland

  END FUNCTION calculate_chi



  SUBROUTINE field_at_particle(part_x, e_at_part, b_at_part)

    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(OUT) :: e_at_part(3), b_at_part(3)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r

    ! The fraction of a cell between the particle position and cell boundary
    REAL(num) :: cell_frac_x

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min:sf_max) :: hx
    ! Temporary variables
    REAL(num) :: cf2
    INTEGER :: dcellx
    REAL(num) :: ex_part, ey_part, ez_part
    REAL(num) :: bx_part, by_part, bz_part

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! Work out the grid cell number for the particle.
    ! Not an integer in general.
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
#else
    cell_x_r = part_x / dx
#endif
    ! Round cell position to nearest cell
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    ! Calculate fraction of cell between nearest cell boundary and particle
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1

    ! Particle weight factors as described in the manual, page25
    ! These weight grid properties onto particles
    ! Also used to weight particle properties onto grid, used later
    ! to calculate J
#ifdef PARTICLE_SHAPE_BSPLINE3
    INCLUDE '../include/bspline3/gx.inc'
#elif  PARTICLE_SHAPE_TOPHAT
    INCLUDE '../include/tophat/gx.inc'
#else
    INCLUDE '../include/triangle/gx.inc'
#endif

    ! Now redo shifted by half a cell due to grid stagger.
    ! Use shifted version for ex in X, ey in Y, ez in Z
    ! And in Y&Z for bx, X&Z for by, X&Y for bz
    cell_x2 = FLOOR(cell_x_r)
    cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
    cell_x2 = cell_x2 + 1

    dcellx = 0
#ifdef PARTICLE_SHAPE_BSPLINE3
    INCLUDE '../include/bspline3/hx_dcell.inc'
#elif  PARTICLE_SHAPE_TOPHAT
    INCLUDE '../include/tophat/hx_dcell.inc'
#else
    INCLUDE '../include/triangle/hx_dcell.inc'
#endif

    ! These are the electric and magnetic fields interpolated to the
    ! particle position. They have been checked and are correct.
    ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
    INCLUDE '../include/bspline3/eb_part.inc'
#elif  PARTICLE_SHAPE_TOPHAT
    INCLUDE '../include/tophat/eb_part.inc'
#else
    INCLUDE '../include/triangle/eb_part.inc'
#endif

    ! update particle momenta using weighted fields
    ! ex_part etc are NOT fields at particle, but fac times
    ! field

    e_at_part(1) = fac * ex_part
    e_at_part(2) = fac * ey_part
    e_at_part(3) = fac * ez_part

    b_at_part(1) = fac * bx_part
    b_at_part(2) = fac * by_part
    b_at_part(3) = fac * bz_part

  END SUBROUTINE field_at_particle



  SUBROUTINE generate_photon(generating_electron, iphoton, eta)

    ! Generates a photon moving in same direction as electron
    ! (generates entirely new photon)

    ! Initialises photons with unit vector along generating electron
    ! momentum and speed c (STORED IN part_p!!!)

    TYPE(particle), POINTER :: generating_electron
    INTEGER, INTENT(IN) :: iphoton
    REAL(num), INTENT(IN) :: eta
    REAL(num) :: mult_x, mult_y, mult_z, mag_p, generating_gamma
    REAL(num) :: rand_temp
    TYPE(particle), POINTER :: new_photon

    ALLOCATE(new_photon)

    new_photon%part_pos = generating_electron%part_pos

    mag_p = MAX(SQRT(generating_electron%part_p(1)**2 &
        + generating_electron%part_p(2)**2 &
        + generating_electron%part_p(3)**2), c_non_zero)

    mult_x = generating_electron%part_p(1) / mag_p
    mult_y = generating_electron%part_p(2) / mag_p
    mult_z = generating_electron%part_p(3) / mag_p

    new_photon%part_p(1) = c * mult_x
    new_photon%part_p(2) = c * mult_y
    new_photon%part_p(3) = c * mult_z

    new_photon%optical_depth = reset_optical_depth()

    generating_gamma = SQRT(1.0_num + (mag_p / m0 / c)**2)

    ! Determine photon energy

    rand_temp = random()
    new_photon%particle_energy = calculate_photon_energy(rand_temp, eta, &
        generating_gamma)
    new_photon%weight = generating_electron%weight

    ! Calculate electron recoil

    mag_p = mag_p - new_photon%particle_energy / c

    generating_electron%part_p(1) = mult_x * mag_p
    generating_electron%part_p(2) = mult_y * mag_p
    generating_electron%part_p(3) = mult_z * mag_p

    ! This will only create photons that have energies above a user specified
    ! cutoff and if photon generation is turned on. E+/- recoil is always
    ! considered
    IF (new_photon%particle_energy .LE. photon_energy_min &
        .OR. .NOT.produce_photons) THEN
      DEALLOCATE(new_photon)
    ELSE
      CALL add_particle_to_partlist(species_list(iphoton)%attached_list, &
          new_photon)
    ENDIF

  END SUBROUTINE generate_photon



  FUNCTION calculate_photon_energy(rand_seed, eta, generating_gamma)

    REAL(num) :: calculate_photon_energy
    REAL(num), INTENT(IN) :: rand_seed, eta, generating_gamma
    REAL(num) :: chi_final

    chi_final = find_value_from_table_alt(rand_seed, eta, log_eta, log_chi, &
        p_photon_energy, n_sample_eta, n_sample_chi)

    calculate_photon_energy = (2.0_num * chi_final / eta) * generating_gamma &
        * m0 * c**2

  END FUNCTION calculate_photon_energy



  SUBROUTINE generate_pair(generating_photon, chi_val, iphoton, ielectron, &
      ipositron)

    ! Generates a pair moving in same direction as photon
    TYPE(particle), POINTER :: generating_photon
    REAL(num), INTENT(IN) :: chi_val
    INTEGER, INTENT(IN) :: iphoton, ielectron, ipositron
    REAL(num) :: mult_x, mult_y, mult_z, mag_p
    REAL(num) :: probability_split, epsilon_frac
    TYPE(particle), POINTER :: new_electron, new_positron

    ALLOCATE(new_electron)
    ALLOCATE(new_positron)

    new_electron%part_pos = generating_photon%part_pos
    new_positron%part_pos = generating_photon%part_pos

    mult_x = generating_photon%part_p(1) / c
    mult_y = generating_photon%part_p(2) / c
    mult_z = generating_photon%part_p(3) / c

    ! Determine how to split the energy amoung e-/e+
    ! IS CHI HERE SAME AS ROLAND'S? DEFINED BSinT/B_s

    probability_split = random()

    epsilon_frac = find_value_from_table(probability_split, chi_val, &
        log_chi2, epsilon_split, p_energy, n_sample_chi2, n_sample_epsilon)

    mag_p = MAX(generating_photon%particle_energy / c, c_non_zero)

    new_electron%part_p(1) = epsilon_frac * mag_p * mult_x
    new_electron%part_p(2) = epsilon_frac * mag_p * mult_y
    new_electron%part_p(3) = epsilon_frac * mag_p * mult_z

    new_positron%part_p(1) = (1.0_num - epsilon_frac) * mag_p * mult_x
    new_positron%part_p(2) = (1.0_num - epsilon_frac) * mag_p * mult_y
    new_positron%part_p(3) = (1.0_num - epsilon_frac) * mag_p * mult_z

    new_electron%optical_depth = reset_optical_depth()
    new_positron%optical_depth = reset_optical_depth()

#ifdef TRIDENT_PHOTONS
    new_electron%optical_depth_tri = reset_optical_depth()
    new_positron%optical_depth_tri = reset_optical_depth()
#endif

    new_electron%weight = generating_photon%weight
    new_positron%weight = generating_photon%weight

    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

    ! Remove photon
    CALL remove_particle_from_partlist(species_list(iphoton)%attached_list, &
        generating_photon)

    DEALLOCATE(generating_photon)

  END SUBROUTINE generate_pair



  SUBROUTINE generate_pair_tri(generating_electron, ielectron, ipositron)

    ! Generates a pair moving in same direction as photon
    TYPE(particle), POINTER :: generating_electron
    INTEGER, INTENT(IN) :: ielectron, ipositron
    TYPE(particle), POINTER :: new_electron, new_positron

    ALLOCATE(new_electron)
    ALLOCATE(new_positron)

    new_electron%part_pos = generating_electron%part_pos
    new_positron%part_pos = generating_electron%part_pos

    new_electron%part_p = 0.0_num
    new_positron%part_p = 0.0_num

    new_electron%optical_depth = reset_optical_depth()
    new_positron%optical_depth = reset_optical_depth()

#ifdef TRIDENT_PHOTONS
    new_electron%optical_depth_tri = reset_optical_depth()
    new_positron%optical_depth_tri = reset_optical_depth()
#endif

    new_electron%weight = generating_electron%weight
    new_positron%weight = generating_electron%weight

    CALL add_particle_to_partlist(species_list(ielectron)%attached_list, &
        new_electron)
    CALL add_particle_to_partlist(species_list(ipositron)%attached_list, &
        new_positron)

  END SUBROUTINE generate_pair_tri



  FUNCTION find_value_from_table_1d(nx, x, x_value, values)

    REAL(num) :: find_value_from_table_1d
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(nx), x_value, values(nx)
    LOGICAL :: l_found_element
    INTEGER :: ix, index_gt, index_lt
    REAL(num) :: fx, value_interp

    l_found_element = .FALSE.
    DO ix = 1, nx
      IF (x(ix) .GT. LOG10(x_value) .AND. .NOT.l_found_element) THEN
        index_gt = ix
        IF (index_gt .GT. 1) THEN
          index_lt = index_gt
        ELSE
          index_lt = ix - 1
        ENDIF
        l_found_element = .TRUE.
      ENDIF
    ENDDO

    IF (.NOT.l_found_element) THEN
      index_gt = nx
      index_lt = nx
    ENDIF

    ! Interpolate in x
    IF (index_lt .GT. index_gt) THEN
      fx = 0.0_num
    ELSE
      fx = (LOG10(x_value) - x(index_lt)) / (x(index_gt) - x(index_lt))
    ENDIF

    value_interp = (1.0_num - fx) * values(index_lt) + fx * values(index_gt)

    find_value_from_table_1d = 10.0_num**value_interp

  END FUNCTION find_value_from_table_1d



  FUNCTION find_value_from_table_alt(p, x_value, x, y, p_table, nx, ny)

    REAL(num) :: find_value_from_table_alt
    REAL(num), INTENT(IN) :: p, x_value, x(nx), y(nx,ny), p_table(nx,ny)
    INTEGER, INTENT(IN) :: nx, ny
    INTEGER :: ix, iy, index_gt, index_lt
    INTEGER :: index_lt_lt, index_lt_gt, index_gt_lt, index_gt_gt
    REAL(num) :: chi_interp
    REAL(num) :: f_eta, f_chi1, chi_low, chi_up
    LOGICAL :: l_found_element

    ! Scan through chi to find correct row of table

    l_found_element = .FALSE.
    DO ix = 1, nx
      IF (x(ix) .GT. LOG10(x_value) .AND. .NOT.l_found_element) THEN
        index_gt = ix
        IF (index_gt .GT. 1) THEN
          index_lt = index_gt
        ELSE
          index_lt = ix - 1
        ENDIF
        l_found_element = .TRUE.
      ENDIF
    ENDDO

    IF (.NOT.l_found_element) THEN
      index_gt = nx
      index_lt = nx
    ENDIF

    ! Scan through table row to find chi

    l_found_element = .FALSE.
    DO iy = 1, ny
      IF (p_table(index_lt,iy) .GT. p .AND. .NOT.l_found_element) THEN
        index_gt_lt = iy
        IF (index_gt_lt .GT. 1) THEN
          index_lt_lt = index_gt_lt
        ELSE
          index_lt_lt = iy - 1
        ENDIF
        l_found_element = .TRUE.
      ENDIF
    ENDDO

    IF (.NOT.l_found_element) THEN
      index_lt_lt = ny
      index_gt_lt = ny
    ENDIF

    IF (index_lt_lt .GT. index_gt_lt) THEN
      f_chi1 = 0.0_num
    ELSE
      f_chi1 = (p - p_table(index_lt,index_lt_lt)) &
          / (p_table(index_lt,index_gt_lt) - p_table(index_lt,index_lt_lt))
    ENDIF

    chi_low = (1.0_num - f_chi1) * y(index_lt,index_lt_lt) &
        + f_chi1 * y(index_lt,index_gt_lt)

    l_found_element = .FALSE.
    DO iy = 1, ny
      IF (p_table(index_gt,iy) .GT. p .AND. .NOT.l_found_element) THEN
        index_gt_gt = iy
        IF (index_gt_gt .GT. 1) THEN
          index_lt_gt = index_gt_gt
        ELSE
          index_lt_gt = iy - 1
        ENDIF
        l_found_element = .TRUE.
      ENDIF
    ENDDO

    IF (.NOT.l_found_element) THEN
      index_lt_gt = ny
      index_gt_gt = ny
    ENDIF

    IF (index_lt_gt .GT. index_gt_gt) THEN
      f_chi1 = 0.0_num
    ELSE
      f_chi1 = (p - p_table(index_gt,index_lt_gt)) &
          / (p_table(index_gt,index_gt_gt) - p_table(index_gt,index_lt_gt))
    ENDIF

    chi_up = (1.0_num - f_chi1) * y(index_gt,index_lt_gt) &
        + f_chi1 * y(index_gt,index_gt_gt)

    IF (index_lt .GT. index_gt) THEN
      f_eta = 0.0_num
    ELSE
      f_eta = (LOG10(x_value) - x(index_lt)) / (x(index_gt) - x(index_lt))
    ENDIF

    chi_interp = (1.0_num - f_eta) * chi_low + f_eta * chi_up

    find_value_from_table_alt = 10.0_num**chi_interp

  END FUNCTION find_value_from_table_alt



  FUNCTION find_value_from_table(p, x_value, x, y, p_table, nx, ny)

    REAL(num) :: find_value_from_table
    REAL(num), INTENT(IN) :: p, x_value, x(nx), y(ny), p_table(nx,ny)
    INTEGER, INTENT(IN) :: nx, ny
    INTEGER :: ix, iy, index_gt, index_lt
    REAL(num) :: f_epsilon, epsilon_interp
    REAL(num) :: f_chi
    LOGICAL :: l_found_element
    REAL(num) :: p_interp(ny)

    ! Scan through chi to find correct row of table

    l_found_element = .FALSE.
    DO ix = 1, nx
      IF (x(ix) .GT. LOG10(x_value) .AND. .NOT.l_found_element) THEN
        index_gt = ix
        IF (index_gt .GT. 1) THEN
          index_lt = index_gt
        ELSE
          index_lt = ix - 1
        ENDIF
        l_found_element = .TRUE.
      ENDIF
    ENDDO

    IF (.NOT.l_found_element) THEN
      index_gt = nx
      index_lt = nx
    ENDIF

    ! Interpolate p in chi

    IF (index_lt .GT. index_gt) THEN
      f_chi = 0.0_num
    ELSE
      f_chi = (LOG10(x_value) - x(index_lt)) / (x(index_gt) - x(index_lt))
    ENDIF

    DO iy = 1, ny
      p_interp(iy) = (1.0_num - f_chi) * p_table(index_lt,iy) &
          + f_chi * p_table(index_gt,iy)
    ENDDO

    ! Scan through table row to find epsilon

    l_found_element = .FALSE.
    DO iy = 1, ny
      IF (p_interp(iy) .GT. p .AND. .NOT.l_found_element) THEN
        index_gt = iy
        IF (index_gt .GT. 1) THEN
          index_lt = index_gt
        ELSE
          index_lt = iy - 1
        ENDIF
        l_found_element = .TRUE.
      ENDIF
    ENDDO

    IF (.NOT.l_found_element) THEN
      index_gt = ny
      index_lt = ny
    ENDIF

    IF (index_lt .GT. index_gt) THEN
      f_epsilon = 0.0_num
    ELSE
      f_epsilon = (p - p_interp(index_lt)) &
          / (p_interp(index_gt) - p_interp(index_lt))
    ENDIF

    epsilon_interp = (1.0_num - f_epsilon) * epsilon_split(index_lt) &
        + f_epsilon * epsilon_split(index_gt)

    find_value_from_table = epsilon_interp

  END FUNCTION find_value_from_table

END MODULE photons
#endif
