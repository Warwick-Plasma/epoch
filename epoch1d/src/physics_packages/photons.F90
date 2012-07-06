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

    REAL(num) :: etalog_min, etalog_max, etalog_dx, chi_min, chi_dx
    REAL(num), ALLOCATABLE :: realbuf(:)
    INTEGER :: i, n, ichi2, iepsilon, ieta, ichi, bufsize, intbuf(6)

    IF (rank .EQ. 0) THEN
      OPEN(unit=lu, file=TRIM(qed_table_location)//'/hsokolov.table', &
          status='OLD')
      READ(lu,*) n_sample_h
      ALLOCATE(log_hsokolov(n_sample_h,2))
      DO ieta = 1, n_sample_h
        READ(lu,*) log_hsokolov(ieta,1), log_hsokolov(ieta,2)
      ENDDO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/pairprod.table', &
          status='OLD')
      READ(lu,*) n_sample_t
      ALLOCATE(log_tpair(n_sample_t,2))
      ALLOCATE(log_omegahat(n_sample_t,2))
      DO ichi = 1, n_sample_t
        READ(lu,*) log_tpair(ichi,1), log_omegahat(ichi,2), log_tpair(ichi,2)
      ENDDO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/ksi_sokolov.table', &
          status='OLD')
      READ(lu,*) n_sample_eta, n_sample_chi, etalog_min, etalog_max
      ALLOCATE(p_photon_energy(n_sample_eta,n_sample_chi))
      DO ieta = 1, n_sample_eta
        READ(lu,*) (p_photon_energy(ieta,ichi), ichi=1,n_sample_chi)
      ENDDO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/chimin.table', &
          status='OLD')
      ALLOCATE(chimin_table(n_sample_eta))
      READ(lu,*) (chimin_table(ieta), ieta=1,n_sample_eta)
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/log_chi2.table', &
          status='OLD')
      READ(lu,*) n_sample_chi2
      ALLOCATE(log_chi2(n_sample_chi2))
      DO ichi2 = 1, n_sample_chi2
        READ(lu,*) log_chi2(ichi2)
      ENDDO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/epsilon.table', &
          status='OLD')
      READ(lu,*) n_sample_epsilon
      ALLOCATE(epsilon_split(n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        READ(lu,*) epsilon_split(iepsilon)
      ENDDO
      CLOSE(unit=lu)

      OPEN(unit=lu, file=TRIM(qed_table_location)//'/energy_split.table', &
          status='OLD')
      ALLOCATE(p_energy(n_sample_chi2,n_sample_epsilon))
      DO ichi2 = 1, n_sample_chi2
        DO iepsilon = 1, n_sample_epsilon
          READ(lu,*) p_energy(ichi2,iepsilon)
        ENDDO
      ENDDO
      CLOSE(unit=lu)

      ! Pack data for broadcasting to other processes

      intbuf(1) = n_sample_h
      intbuf(2) = n_sample_t
      intbuf(3) = n_sample_eta
      intbuf(4) = n_sample_chi
      intbuf(5) = n_sample_chi2
      intbuf(6) = n_sample_epsilon

      CALL MPI_BCAST(intbuf, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)

      bufsize = n_sample_h * 2
      bufsize = bufsize + n_sample_t * 3
      bufsize = bufsize + 2 + n_sample_eta * n_sample_chi
      bufsize = bufsize + n_sample_eta
      bufsize = bufsize + n_sample_chi2
      bufsize = bufsize + n_sample_epsilon
      bufsize = bufsize + n_sample_chi2 * n_sample_epsilon

      ALLOCATE(realbuf(bufsize))
      n = 1

      DO i = 1, 2
        DO ieta = 1, n_sample_h
          realbuf(n) = log_hsokolov(ieta,i)
          n = n + 1
        ENDDO
      ENDDO

      DO ichi = 1, n_sample_t
        realbuf(n) = log_tpair(ichi,1)
        n = n + 1
        realbuf(n) = log_tpair(ichi,2)
        n = n + 1
        realbuf(n) = log_omegahat(ichi,2)
        n = n + 1
      ENDDO

      realbuf(n) = etalog_min
      n = n + 1
      realbuf(n) = etalog_max
      n = n + 1

      DO ichi = 1, n_sample_chi
        DO ieta = 1, n_sample_eta
          realbuf(n) = p_photon_energy(ieta,ichi)
          n = n + 1
        ENDDO
      ENDDO

      DO ieta = 1, n_sample_eta
        realbuf(n) = chimin_table(ieta)
        n = n + 1
      ENDDO

      DO ichi2 = 1, n_sample_chi2
        realbuf(n) = log_chi2(ichi2)
        n = n + 1
      ENDDO

      DO iepsilon = 1, n_sample_epsilon
        realbuf(n) = epsilon_split(iepsilon)
        n = n + 1
      ENDDO

      DO iepsilon = 1, n_sample_epsilon
        DO ichi2 = 1, n_sample_chi2
          realbuf(n) = p_energy(ichi2,iepsilon)
          n = n + 1
        ENDDO
      ENDDO

      CALL MPI_BCAST(realbuf, bufsize, mpireal, 0, MPI_COMM_WORLD, errcode)

      DEALLOCATE(realbuf)
    ELSE
      ! Unpack data from rank zero

      CALL MPI_BCAST(intbuf, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)

      n_sample_h       = intbuf(1)
      n_sample_t       = intbuf(2)
      n_sample_eta     = intbuf(3)
      n_sample_chi     = intbuf(4)
      n_sample_chi2    = intbuf(5)
      n_sample_epsilon = intbuf(6)

      bufsize = n_sample_h * 2
      bufsize = bufsize + n_sample_t * 3
      bufsize = bufsize + 2 + n_sample_eta * n_sample_chi
      bufsize = bufsize + n_sample_eta
      bufsize = bufsize + n_sample_chi2
      bufsize = bufsize + n_sample_epsilon
      bufsize = bufsize + n_sample_chi2 * n_sample_epsilon

      ALLOCATE(realbuf(bufsize))
      n = 1

      CALL MPI_BCAST(realbuf, bufsize, mpireal, 0, MPI_COMM_WORLD, errcode)

      ALLOCATE(log_hsokolov(n_sample_h,2))
      DO i = 1, 2
        DO ieta = 1, n_sample_h
          log_hsokolov(ieta,i) = realbuf(n)
          n = n + 1
        ENDDO
      ENDDO

      ALLOCATE(log_tpair(n_sample_t,2))
      ALLOCATE(log_omegahat(n_sample_t,2))
      DO ichi = 1, n_sample_t
        log_tpair(ichi,1)    = realbuf(n)
        n = n + 1
        log_tpair(ichi,2)    = realbuf(n)
        n = n + 1
        log_omegahat(ichi,2) = realbuf(n)
        n = n + 1
      ENDDO

      etalog_min = realbuf(n)
      n = n + 1
      etalog_max = realbuf(n)
      n = n + 1

      ALLOCATE(p_photon_energy(n_sample_eta,n_sample_chi))
      DO ichi = 1, n_sample_chi
        DO ieta = 1, n_sample_eta
          p_photon_energy(ieta,ichi) = realbuf(n)
          n = n + 1
        ENDDO
      ENDDO

      ALLOCATE(chimin_table(n_sample_eta))
      DO ieta = 1, n_sample_eta
        chimin_table(ieta) = realbuf(n)
        n = n + 1
      ENDDO

      ALLOCATE(log_chi2(n_sample_chi2))
      DO ichi2 = 1, n_sample_chi2
        log_chi2(ichi2) = realbuf(n)
        n = n + 1
      ENDDO

      ALLOCATE(epsilon_split(n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        epsilon_split(iepsilon) = realbuf(n)
        n = n + 1
      ENDDO

      ALLOCATE(p_energy(n_sample_chi2,n_sample_epsilon))
      DO iepsilon = 1, n_sample_epsilon
        DO ichi2 = 1, n_sample_chi2
          p_energy(ichi2,iepsilon) = realbuf(n)
          n = n + 1
        ENDDO
      ENDDO

      DEALLOCATE(realbuf)
    ENDIF

    log_omegahat(:,1) = log_tpair(:,1)

    ALLOCATE(log_eta(n_sample_eta))
    ALLOCATE(log_chi(n_sample_eta,n_sample_chi))

    etalog_dx = (etalog_max - etalog_min) / REAL(n_sample_eta-1,num)
    DO ieta = 1, n_sample_eta
      log_eta(ieta) = etalog_min + REAL(ieta-1,num) * etalog_dx
      chi_min = LOG10(chimin_table(ieta))
      chi_dx  = (log_eta(ieta) - LOG10(2.0_num) - chi_min) &
          / REAL(n_sample_chi-1,num)
      DO ichi = 1, n_sample_chi
        log_chi(ieta,ichi) = chi_min + REAL(ichi-1,num) * chi_dx
      ENDDO
    ENDDO

    DEALLOCATE(chimin_table)

  END SUBROUTINE setup_tables_qed



  SUBROUTINE deallocate_tables_qed

    DEALLOCATE(log_chi2)
    DEALLOCATE(epsilon_split)
    DEALLOCATE(p_energy)
    DEALLOCATE(log_hsokolov)
    DEALLOCATE(log_eta)
    DEALLOCATE(log_chi)
    DEALLOCATE(p_photon_energy)
    DEALLOCATE(log_tpair)
    DEALLOCATE(log_omegahat)

  END SUBROUTINE deallocate_tables_qed



  SUBROUTINE initialise_optical_depth(current_species)

    ! Resets optical depth (to random number) of all particles
    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    current => current_species%attached_list%head
    DO WHILE(ASSOCIATED(current))
      p_tau = random()
      current%optical_depth = -LOG(1.0_num - p_tau)

#ifdef TRIDENT_PHOTONS
      p_tau = random()
      current%optical_depth_tri = -LOG(1.0_num - p_tau)
#endif
      current => current%next
    ENDDO

  END SUBROUTINE initialise_optical_depth



  FUNCTION reset_optical_depth()

    ! Resets optical depth of particle
    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth



  SUBROUTINE qed_update_optical_depth

    ! Updates the optical depth for electrons and photons
    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next_pt

    REAL(num) :: part_x
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: dir_x, dir_y, dir_z
    REAL(num) :: eta, chi_val, part_e, gamma

    DO ispecies = 1, n_species

      ! First consider electrons
      IF (species_list(ispecies)%species_type .EQ. c_species_id_electron) THEN
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          ! Find eta at particle position
          part_x  = current%part_pos - x_min_local
          part_ux = current%part_p(1) / mc0
          part_uy = current%part_p(2) / mc0
          part_uz = current%part_p(3) / mc0
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

          eta = calculate_eta(part_x, part_ux, part_uy, &
              part_uz, gamma)

          current%optical_depth = &
              current%optical_depth - delta_optical_depth(eta, gamma)
#ifdef TRIDENT_PHOTONS
          current%optical_depth_tri = &
              current%optical_depth_tri - delta_optical_depth_tri(eta, gamma)
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
          part_ux = current%part_p(1) / mc0
          part_uy = current%part_p(2) / mc0
          part_uz = current%part_p(3) / mc0
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

          eta = calculate_eta(part_x, part_ux, part_uy, &
              part_uz, gamma)

          current%optical_depth = &
              current%optical_depth - delta_optical_depth(eta, gamma)
#ifdef TRIDENT_PHOTONS
          current%optical_depth_tri = &
              current%optical_depth_tri - delta_optical_depth_tri(eta, gamma)
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
            dir_x = current%part_p(1) / c
            dir_y = current%part_p(2) / c
            dir_z = current%part_p(3) / c
            part_e  = current%particle_energy / m0 / c**2
            chi_val = calculate_chi(part_x, dir_x, dir_y, &
                dir_z, part_e)

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



  FUNCTION delta_optical_depth(eta, gamma)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth
    REAL(num), INTENT(IN) :: eta, gamma
    REAL(num) :: hsokolov

    hsokolov = find_value_from_table_1d(eta, n_sample_h, log_hsokolov(:,1), &
        log_hsokolov(:,2))

    delta_optical_depth = dt * eta * alpha_f * SQRT(3.0_num) * hsokolov &
        / (2.0_num * pi * tau_c * gamma)

  END FUNCTION delta_optical_depth



  FUNCTION delta_optical_depth_tri(eta, gamma)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_tri
    REAL(num), INTENT(IN) :: eta, gamma
    REAL(num) :: omegahat

    omegahat = find_value_from_table_1d(eta, n_sample_t, log_omegahat(:,1), &
        log_omegahat(:,2))

    delta_optical_depth_tri = dt * eta * alpha_f**2 * 0.64_num * omegahat &
        / (2.0_num * pi * tau_c * gamma)

  END FUNCTION delta_optical_depth_tri



  FUNCTION delta_optical_depth_photon(chi_val, part_e)

    ! Function that calcualtes the change to the optical depth
    REAL(num) :: delta_optical_depth_photon
    REAL(num), INTENT(IN) :: chi_val, part_e
    REAL(num) :: tpair

    tpair = find_value_from_table_1d(chi_val, n_sample_t, log_tpair(:,1), &
        log_tpair(:,2))

    delta_optical_depth_photon = dt / tau_c * alpha_f / part_e * chi_val * tpair

  END FUNCTION delta_optical_depth_photon



  FUNCTION calculate_eta(part_x, part_ux, part_uy, part_uz, &
      gamma)

    REAL(num) :: calculate_eta
    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(IN) :: part_ux, part_uy, part_uz, gamma
    REAL(num) :: e_at_part(3), b_at_part(3)
    REAL(num) :: beta_x, beta_y, beta_z
    REAL(num) :: flperp(3), i_e, tau0, roland_eta
    REAL(num) :: lambdac, coeff_eta, moduclip, moduclip2, u_dot_e

    CALL field_at_particle(part_x, e_at_part, b_at_part)

    moduclip2 = MAX(part_ux**2 + part_uy**2 + part_uz**2, c_non_zero)
    moduclip = SQRT(moduclip2)

    beta_x = part_ux / gamma
    beta_y = part_uy / gamma
    beta_z = part_uz / gamma

    lambdac = h_bar / mc0

    coeff_eta = SQRT(3.0_num * lambdac / (2.0_num * alpha_f * m0 * c**3))

    u_dot_e = (part_ux * e_at_part(1) + part_uy * e_at_part(2) &
        + part_uz * e_at_part(3)) / moduclip2

    flperp(1) = q0 * (e_at_part(1) - u_dot_e * part_ux &
        + c * (beta_y * b_at_part(3) - beta_z * b_at_part(2)))

    flperp(2) = q0 * (e_at_part(2) - u_dot_e * part_uy &
        + c * (beta_z * b_at_part(1) - beta_x * b_at_part(3)))

    flperp(3) = q0 * (e_at_part(3) - u_dot_e * part_uz &
        + c * (beta_x * b_at_part(2) - beta_y * b_at_part(1)))

    ! Dipole emission intensity

    tau0 = q0**2 / (6.0_num * pi * epsilon0 * m0 * c**3)

    i_e = tau0 * gamma**2 * (flperp(1)**2 + flperp(2)**2 + flperp(3)**2 &
        + (q0 * (beta_x * e_at_part(1) + beta_y * e_at_part(2) &
        + beta_z * e_at_part(3)) / moduclip)**2) / m0

    roland_eta = coeff_eta * SQRT(i_e)

    ! Determine eta from fields
    calculate_eta = roland_eta

  END FUNCTION calculate_eta



  FUNCTION calculate_chi(part_x, dir_x, dir_y, dir_z, part_e)

    REAL(num) :: calculate_chi
    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(IN) :: dir_x, dir_y, dir_z, part_e
    REAL(num) :: e_at_part(3), b_at_part(3), q(3)
    REAL(num) :: e_dot_dir, calculate_chi_roland

    CALL field_at_particle(part_x, e_at_part, b_at_part)

    e_dot_dir = e_at_part(1) * dir_x + e_at_part(2) * dir_y &
        + e_at_part(3) * dir_z

    q(1) = e_at_part(1) - e_dot_dir * dir_x &
        + c * (dir_y * b_at_part(3) - dir_z * b_at_part(2))
    q(2) = e_at_part(2) - e_dot_dir * dir_y &
        + c * (dir_z * b_at_part(1) - dir_x * b_at_part(3))
    q(3) = e_at_part(3) - e_dot_dir * dir_z &
        + c * (dir_x * b_at_part(2) - dir_y * b_at_part(1))

    calculate_chi_roland = 0.5_num * SQRT(q(1)**2 + q(2)**2 + q(3)**2) &
        * part_e / e_s

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
    REAL(num) :: dir_x, dir_y, dir_z, mag_p, generating_gamma
    REAL(num) :: rand_temp
    TYPE(particle), POINTER :: new_photon

    ALLOCATE(new_photon)

    new_photon%part_pos = generating_electron%part_pos

    mag_p = MAX(SQRT(generating_electron%part_p(1)**2 &
        + generating_electron%part_p(2)**2 &
        + generating_electron%part_p(3)**2), c_non_zero)

    dir_x = generating_electron%part_p(1) / mag_p
    dir_y = generating_electron%part_p(2) / mag_p
    dir_z = generating_electron%part_p(3) / mag_p

    new_photon%part_p(1) = c * dir_x
    new_photon%part_p(2) = c * dir_y
    new_photon%part_p(3) = c * dir_z

    new_photon%optical_depth = reset_optical_depth()

    generating_gamma = SQRT(1.0_num + (mag_p / m0 / c)**2)

    ! Determine photon energy

    rand_temp = random()
    new_photon%particle_energy = calculate_photon_energy(rand_temp, eta, &
        generating_gamma)
    new_photon%weight = generating_electron%weight

    ! Calculate electron recoil

    mag_p = mag_p - new_photon%particle_energy / c

    generating_electron%part_p(1) = dir_x * mag_p
    generating_electron%part_p(2) = dir_y * mag_p
    generating_electron%part_p(3) = dir_z * mag_p

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

    chi_final = find_value_from_table_alt(eta, rand_seed, &
        n_sample_eta, n_sample_chi, log_eta, log_chi, p_photon_energy)

    calculate_photon_energy = (2.0_num * chi_final / eta) * generating_gamma &
        * m0 * c**2

  END FUNCTION calculate_photon_energy



  SUBROUTINE generate_pair(generating_photon, chi_val, iphoton, ielectron, &
      ipositron)

    ! Generates a pair moving in same direction as photon
    TYPE(particle), POINTER :: generating_photon
    REAL(num), INTENT(IN) :: chi_val
    INTEGER, INTENT(IN) :: iphoton, ielectron, ipositron
    REAL(num) :: dir_x, dir_y, dir_z, mag_p
    REAL(num) :: probability_split, epsilon_frac
    TYPE(particle), POINTER :: new_electron, new_positron

    ALLOCATE(new_electron)
    ALLOCATE(new_positron)

    new_electron%part_pos = generating_photon%part_pos
    new_positron%part_pos = generating_photon%part_pos

    dir_x = generating_photon%part_p(1) / c
    dir_y = generating_photon%part_p(2) / c
    dir_z = generating_photon%part_p(3) / c

    ! Determine how to split the energy amoung e-/e+
    ! IS CHI HERE SAME AS ROLAND'S? DEFINED BSinT/B_s

    probability_split = random()

    epsilon_frac = find_value_from_table(chi_val, probability_split, &
        n_sample_chi2, n_sample_epsilon, log_chi2, epsilon_split, p_energy)

    mag_p = MAX(generating_photon%particle_energy / c, c_non_zero)

    new_electron%part_p(1) = epsilon_frac * mag_p * dir_x
    new_electron%part_p(2) = epsilon_frac * mag_p * dir_y
    new_electron%part_p(3) = epsilon_frac * mag_p * dir_z

    new_positron%part_p(1) = (1.0_num - epsilon_frac) * mag_p * dir_x
    new_positron%part_p(2) = (1.0_num - epsilon_frac) * mag_p * dir_y
    new_positron%part_p(3) = (1.0_num - epsilon_frac) * mag_p * dir_z

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



  FUNCTION find_value_from_table_1d(x_in, nx, x, values)

    REAL(num) :: find_value_from_table_1d
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(nx), values(nx)
    REAL(num) :: fx, x_value, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(MAX(x_in,c_non_zero))

    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 .LT. 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm .LT. 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        ENDIF
        IF (i2 - i1 .EQ. 1) EXIT
      ENDDO
      ! Interpolate in x
      fx = (x_value - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_1d" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      ENDIF
      IF (xdif1 .GE. 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      ENDIF
    ENDIF

    value_interp = (1.0_num - fx) * values(i1) + fx * values(i2)

    find_value_from_table_1d = 10.0_num**value_interp

  END FUNCTION find_value_from_table_1d



  FUNCTION find_value_from_table_alt(x_in, p_value, nx, ny, x, y, p_table)

    REAL(num) :: find_value_from_table_alt
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(nx,ny), p_table(nx,ny)
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, x_value, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(x_in)

    ! Scan through x to find correct row of table
    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 .LT. 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm .LT. 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        ENDIF
        IF (i2 - i1 .EQ. 1) EXIT
      ENDDO
      ! Interpolate in x
      fx = (x_value - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      ENDIF
      IF (xdif1 .GE. 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      ENDIF
    ENDIF

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 .LT. 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm .LT. 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        ENDIF
        IF (i2 - i1 .EQ. 1) EXIT
      ENDDO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      ENDIF
      IF (xdif1 .GE. 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      ENDIF
    ENDIF

    y_lt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 .LT. 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm .LT. 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        ENDIF
        IF (i2 - i1 .EQ. 1) EXIT
      ENDDO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      ENDIF
      IF (xdif1 .GE. 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      ENDIF
    ENDIF

    y_gt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ! Interpolate in x

    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table_alt = 10.0_num**y_interp

  END FUNCTION find_value_from_table_alt



  FUNCTION find_value_from_table(x_in, p_value, nx, ny, x, y, p_table)

    REAL(num) :: find_value_from_table
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(ny), p_table(nx,ny)
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, x_value, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = LOG10(x_in)

    ! Scan through x to find correct row of table
    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 .LT. 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm .LT. 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        ENDIF
        IF (i2 - i1 .EQ. 1) EXIT
      ENDDO
      ! Interpolate in x
      fx = (x_value - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      ENDIF
      IF (xdif1 .GE. 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      ENDIF
    ENDIF

    index_lt = i1
    index_gt = i2

    ix = index_lt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 .LT. 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm .LT. 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        ENDIF
        IF (i2 - i1 .EQ. 1) EXIT
      ENDDO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      ENDIF
      IF (xdif1 .GE. 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      ENDIF
    ENDIF

    y_lt = (1.0_num - fp) * y(i1) + fp * y(i2)

    ix = index_gt
    ! Scan through table row to find p_value
    i1 = 1
    i2 = ny
    xdif1 = p_table(ix,i1) - p_value
    xdif2 = p_table(ix,i2) - p_value
    IF (xdif1 * xdif2 .LT. 0) THEN
      ! Use bisection to find the nearest cell
      DO
        im = (i1 + i2) / 2
        xdifm = p_table(ix,im) - p_value
        IF (xdif1 * xdifm .LT. 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        ENDIF
        IF (i2 - i1 .EQ. 1) EXIT
      ENDDO
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank .EQ. 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      ENDIF
      IF (xdif1 .GE. 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      ENDIF
    ENDIF

    y_gt = (1.0_num - fp) * y(i1) + fp * y(i2)

    ! Interpolate in x

    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table = y_interp

  END FUNCTION find_value_from_table

END MODULE photons
#endif
