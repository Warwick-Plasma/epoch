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

MODULE bremsstrahlung
#ifdef BREMSSTRAHLUNG

  USE partlist
  USE calc_df

  IMPLICIT NONE

CONTAINS

  ! Confirm that the correct particles are present to actually use the
  ! bremsstrahlung routines, and warn the user if they're not. Then call
  ! subroutines to initialise the bremsstrahlung table array, and give all
  ! particles an optical depth for bremsstrahlung emission (if appropriate)
  SUBROUTINE setup_bremsstrahlung_module

    INTEGER :: ispecies, iu, io
    LOGICAL :: found

    IF (produce_bremsstrahlung_photons .AND. .NOT. ic_from_restart) THEN
      ! Look for electron and photon species. This is a secondary sanity check -
      ! full check appears in check_bremsstrahlung_variables
      found = .FALSE.
      DO ispecies = 1, n_species
        IF (species_list(ispecies)%count <= 0) CYCLE
        IF (species_list(ispecies)%species_type == c_species_id_electron &
            .OR. ispecies == bremsstrahlung_photon_species ) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO

     ! Print warning if the correct species aren't present
      IF (.NOT. found) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Electron and photon species are either ', &
                'unspecified or contain no'
            WRITE(io,*) 'particles. Bremsstrahlung routines will do nothing.'
          END DO
        END IF
      END IF
    END IF

    ! Load the tables for the bremsstrahlung routines
    CALL setup_tables_bremsstrahlung

    IF (.NOT. ic_from_restart) THEN
      DO ispecies = 1, n_species
        CALL initialise_optical_depth(species_list(ispecies))
      END DO
    END IF

  END SUBROUTINE setup_bremsstrahlung_module



  ! Function to ensure all required species are present
  FUNCTION check_bremsstrahlung_variables()
    INTEGER :: check_bremsstrahlung_variables
    INTEGER :: io, iu
    INTEGER :: ispecies
    INTEGER :: first_electron = -1

    check_bremsstrahlung_variables = c_err_none

    ! No special species required if bremsstrahlung is turned off
    IF (.NOT.use_bremsstrahlung) RETURN

    ! No special species required if we only do radiation reaction
    IF (.NOT.produce_bremsstrahlung_photons) RETURN

    ! Identify if there exists any electron species
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
          WRITE(io,*) 'Bremsstrahlung routines require at least one species' , &
              ' of electrons.'
        END DO
      END IF
      check_bremsstrahlung_variables = c_err_missing_elements
      RETURN
    END IF

#ifdef PHOTONS
    ! photon_species can act as bremsstrahlung_photon_species if no
    ! bremsstrahlung species is defined
    IF (bremsstrahlung_photon_species == -1 &
        .AND. .NOT. photon_species == -1 ) THEN
      bremsstrahlung_photon_species = photon_species
    END IF
#endif

    ! Check if there exists a species to populate with bremsstrahlung photons
    IF (bremsstrahlung_photon_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No photon species specified. Specify using ', &
              '"identify:brem_photon"'
        END DO
      END IF
      check_bremsstrahlung_variables = c_err_missing_elements
      RETURN
    END IF

  END FUNCTION check_bremsstrahlung_variables



  ! Subroutine to generate tables of bremsstrahlung cross sections as a function
  ! of incident electron energy, T. Also computes a table of cumlative
  ! probabilities for each photon energy, k, and each value of T. All tables are
  ! stored in an array (brem_array) as the variable type brem_tables. Each
  ! element of brem_array correpsonds to a different target atomic number.
  SUBROUTINE setup_tables_bremsstrahlung

    INTEGER :: Z_temp, iu, io
    INTEGER :: i_species, iZ, jZ
    INTEGER, ALLOCATABLE :: Z_flags(:)
    CHARACTER(LEN=3) :: Z_string
    INTEGER :: size_k, size_T
    INTEGER :: i, j
    INTEGER :: buf_index, buf_size
    INTEGER, ALLOCATABLE :: int_buf(:)
    REAL(num), ALLOCATABLE :: real_buf(:)

    ! Do all analysis on rank 0, then MPI_BROADCAST to other ranks
    IF (rank == 0) THEN

      ! For each unique atomic number, Z, let Z_flags(Z) = 1
      ALLOCATE(Z_flags(100))
      Z_flags(:) = 0
      DO i_species = 1, n_species
        Z_temp = species_list(i_species)%atomic_no
        IF (Z_temp > 0 .AND. Z_temp < 101) THEN
          Z_flags(Z_temp) = 1;

        ! We only have tables up to Z=100
        ELSE IF (Z_temp > 100) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Species with atomic numbers over 100 cannot be ', &
                'modelled using bremsstrahlung libraries, and will be ignored'
          END DO
        END IF
      END DO

      ! Create brem_array to store tables for each atomic number, Z. The value
      ! of Z_to_index(Z) is the brem_array index for tables of atomic number Z.
      ! Z_values contains the Z values present.
      size_brem_array = SUM(Z_flags)
      ALLOCATE(brem_array(size_brem_array))
      ALLOCATE(Z_values(size_brem_array))
      ALLOCATE(Z_to_index(100))
      iZ = 1
      DO jZ = 1 ,100
        IF (Z_flags(jZ) == 1) THEN
          Z_to_index(jZ) = iZ
          Z_values(iZ) = jZ
          iZ = iZ + 1
        ELSE
          Z_to_index(jZ) = 0
        END IF
      END DO
      DEALLOCATE(Z_flags)

      ! Obtain tables for each brem_array value
      DO iZ = 1,size_brem_array
        Z_temp = Z_values(iZ)

        ! Open table for given Z value
        IF (Z_temp < 10) THEN
          WRITE(Z_string, '(I1)') Z_temp
        ELSE IF (Z_temp < 100) THEN
          WRITE(Z_string, '(I2)') Z_temp
        ELSE IF (Z_temp == 100) THEN
          WRITE(Z_string, '(I3)') Z_temp
        ELSE
          CYCLE
        END IF
        OPEN(unit = lu, file = TRIM(bremsstrahlung_table_location) // '/br' &
            // Z_string, status = 'OLD')


        ! Extract data from the table file
        ! File format: (number of k values) (number of e- energy values)
        !              Line of electron eneries [J]
        !              Line of corresponding cross sections [m²]
        !              Table of k values [J]
        !              Table of CDF values for photon emission
        ! In the table, rows correspond to constant electron energy, and columns
        ! at different photon energies, k. Data taken from Geant4 scaled
        ! differential cross section libraries, and processed to include a cut
        ! off k for each electron energy. This ignores a fraction (1e-9) of the
        ! radiated energy produced by the lowest energy photons. This is for
        ! consistency with photons.F90.
        READ(lu,*) size_k, size_T
        ALLOCATE(brem_array(iZ)%cross_section(size_T))
        ALLOCATE(brem_array(iZ)%k_table(size_T,size_k))
        ALLOCATE(brem_array(iZ)%cdf_table(size_T, size_k))
        ALLOCATE(brem_array(iZ)%E_table(size_T))
        brem_array(iZ)%size_T = size_T
        brem_array(iZ)%size_k = size_k

        ! Read electron energies and cross sections
        READ(lu,*) brem_array(iZ)%E_table(1:size_T)
        READ(lu,*) brem_array(iZ)%cross_section(1:size_T)

        ! Read table of k values
        DO i = 1, size_T
            READ(lu,*) brem_array(iZ)%k_table(i,1:size_k)
        END DO

        ! Read table of cdf values
        DO i = 1, size_T
            READ(lu,*) brem_array(iZ)%cdf_table(i,1:size_k)
        END DO

        CLOSE(unit = lu)

      END DO

      ! Pack data for broadcast to other ranks
      ! How many brem_tables do we need to load in brem_array?
      CALL MPI_BCAST(size_brem_array, 1, MPI_INTEGER, 0, comm, errcode)

      ! Pack Z values and array sizes into int_buf
      buf_size = 3*size_brem_array + 100
      ALLOCATE(int_buf(buf_size))
      buf_index = 1
      DO i = 1, size_brem_array
        int_buf(buf_index) = Z_values(i)
        buf_index = buf_index + 1
      END DO
      DO i = 1, 100
        int_buf(buf_index) = Z_to_index(i)
        buf_index = buf_index + 1
      END DO
      DO i = 1, size_brem_array
        int_buf(buf_index) = brem_array(i)%size_k
        buf_index = buf_index + 1
        int_buf(buf_index) = brem_array(i)%size_T
        buf_index = buf_index + 1
      END DO
      CALL MPI_BCAST(int_buf, buf_size, MPI_INTEGER, 0, comm, errcode)

      ! Pack array data into real_buf and broadcast
      buf_size = 0
      DO i = 1, size_brem_array
        buf_size = brem_array(i)%size_T*(2*brem_array(i)%size_k + 2) + buf_size
      END DO
      ALLOCATE(real_buf(buf_size))
      buf_index = 1
      DO i = 1, size_brem_array
        size_k = brem_array(i)%size_k
        DO j = 1, brem_array(i)%size_T
          real_buf(buf_index:buf_index+size_k-1) = brem_array(i)%cdf_table(j,:)
          buf_index = buf_index + size_k
          real_buf(buf_index:buf_index+size_k-1) = brem_array(i)%k_table(j,:)
          buf_index = buf_index + size_k
          real_buf(buf_index) = brem_array(i)%cross_section(j)
          buf_index = buf_index + 1
          real_buf(buf_index) = brem_array(i)%E_table(j)
          buf_index = buf_index + 1
        END DO
      END DO
      CALL MPI_BCAST(real_buf, buf_size, mpireal, 0, comm, errcode)

    ! All other ranks
    ELSE

      ! First call to get number of Z values present in simulation
      CALL MPI_BCAST(size_brem_array, 1, MPI_INTEGER, 0, comm, errcode)

      ALLOCATE(brem_array(size_brem_array))

      ! Second call to get sizes to allocate arrays, and integer Z values
      buf_size = 3*size_brem_array + 100
      ALLOCATE(int_buf(buf_size))
      CALL MPI_BCAST(int_buf, buf_size, MPI_INTEGER, 0, comm, errcode)

      ! Get Z arrays
      ALLOCATE(Z_values(size_brem_array))
      ALLOCATE(Z_to_index(100))
      buf_index = 1
      DO i = 1, size_brem_array
        Z_values(i) = int_buf(buf_index)
        buf_index = buf_index + 1
      END DO
      DO i = 1, 100
        Z_to_index(i) = int_buf(buf_index)
        buf_index = buf_index + 1
      END DO

      ! Get brem table sizes
      DO i = 1,size_brem_array
        brem_array(i)%size_k = int_buf(buf_index)
        buf_index = buf_index + 1
        brem_array(i)%size_T = int_buf(buf_index)
        buf_index = buf_index + 1
      END DO

      ! Allocate brem tables
      DO i = 1, size_brem_array
        size_T = brem_array(i)%size_T
        size_k = brem_array(i)%size_k
        ALLOCATE(brem_array(i)%cdf_table(size_T,size_k))
        ALLOCATE(brem_array(i)%k_table(size_T,size_k))
        ALLOCATE(brem_array(i)%cross_section(size_T))
        ALLOCATE(brem_array(i)%E_table(size_T))
      END DO

      ! Final MPI call to get data for tables
      buf_size = 0
      DO i = 1,size_brem_array
        buf_size = brem_array(i)%size_T*(2*brem_array(i)%size_k + 2) + buf_size
      END DO
      ALLOCATE(real_buf(buf_size))
      CALL MPI_BCAST(real_buf, buf_size, mpireal, 0, comm, errcode)

      ! Load real_buf into brem_array
      buf_index = 1
      DO i = 1, size_brem_array
        size_k = brem_array(i)%size_k
        DO j = 1, brem_array(i)%size_T
          brem_array(i)%cdf_table(j,:) = real_buf(buf_index:buf_index+size_k-1)
          buf_index = buf_index + size_k
          brem_array(i)%k_table(j,:) = real_buf(buf_index:buf_index+size_k-1)
          buf_index = buf_index + size_k
          brem_array(i)%cross_section(j) = real_buf(buf_index)
          buf_index = buf_index + 1
          brem_array(i)%E_table(j) = real_buf(buf_index)
          buf_index = buf_index + 1
        END DO
      END DO

    END IF

    DEALLOCATE(int_buf)
    DEALLOCATE(real_buf)

  END SUBROUTINE setup_tables_bremsstrahlung



  ! Loops through all particles in a species and sets a bremsstrahlung optical
  ! depth
  SUBROUTINE initialise_optical_depth(current_species)

    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    ! Randomly sample the optical depth from an exponential distribution
    current => current_species%attached_list%head
    DO WHILE(ASSOCIATED(current))
      p_tau = random()
      current%optical_depth_bremsstrahlung = -LOG(1.0_num - p_tau)
      current => current%next
    END DO

  END SUBROUTINE initialise_optical_depth



  ! Deallocates the shared_data created in setup_tables_bremsstrahlung
  SUBROUTINE deallocate_tables_bremsstrahlung

    INTEGER :: i

    DO i = 1,size_brem_array
      DEALLOCATE(brem_array(i)%k_table)
      DEALLOCATE(brem_array(i)%cdf_table)
      DEALLOCATE(brem_array(i)%cross_section)
      DEALLOCATE(brem_array(i)%E_table)
    END DO

    DEALLOCATE(brem_array)

  END SUBROUTINE deallocate_tables_bremsstrahlung



  ! Updates the optical depth for electrons. This subroutine is responsible for
  ! also calling the function which calculates the optical depth change, and
  ! calling the generate_photon subroutine. This subroutine serves as the main
  ! interface to the bremsstrahlung module for main-loop processes in
  ! epoch3d.F90
  SUBROUTINE bremsstrahlung_update_optical_depth

    INTEGER :: ispecies, iZ, Z_temp, iArray, jArray, kArray
    TYPE(particle), POINTER :: current
    REAL(num), ALLOCATABLE :: grid_num_density_electron_temp(:,:,:)
    REAL(num), ALLOCATABLE :: grid_num_density_electron(:,:,:)
    REAL(num), ALLOCATABLE :: grid_temperature_electron_temp(:,:,:)
    REAL(num), ALLOCATABLE :: grid_temperature_electron(:,:,:)
    REAL(num), ALLOCATABLE :: grid_root_temp_over_num(:,:,:)
    REAL(num) :: grid_num_density_ion(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng)
    REAL(num) :: part_ux, part_uy, part_uz, gamma_rel
    REAL(num) :: part_x, part_y, part_z, part_E, part_ni, part_v
    REAL(num) :: part_root_Te_over_ne, plasma_factor

    ! Calculate electron number density and temperature, summed over all
    ! electron species
    IF (use_plasma_screening) THEN

      ALLOCATE(grid_num_density_electron_temp(1-ng:nx+ng, 1-ng:ny+ng, &
          1-ng:nz+ng))
      ALLOCATE(grid_num_density_electron(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      ALLOCATE(grid_temperature_electron_temp(1-ng:nx+ng, 1-ng:ny+ng, &
          1-ng:nz+ng))
      ALLOCATE(grid_temperature_electron(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      ALLOCATE(grid_root_temp_over_num(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
      grid_num_density_electron = 0.0_num
      grid_temperature_electron = 0.0_num

      DO ispecies = 1, n_species
        IF (species_list(ispecies)%species_type == c_species_id_electron) THEN
          CALL calc_number_density(grid_num_density_electron_temp, ispecies)
          CALL calc_temperature(grid_temperature_electron_temp, ispecies)

          ! Sum the densities, perform a weighted mean to find mean temperature
          grid_num_density_electron = grid_num_density_electron &
              + grid_num_density_electron_temp
          grid_temperature_electron = grid_temperature_electron &
              + grid_temperature_electron_temp*grid_num_density_electron_temp
        END IF
      END DO
      grid_temperature_electron = grid_temperature_electron &
          / grid_num_density_electron

      ! Create a grid of sqrt(Te/ne) values
      DO kArray = 1-ng, nz+ng
        DO jArray = 1-ng, ny+ng
          DO iArray = 1-ng, nx+ng
            IF (grid_num_density_electron(iArray, jArray, kArray) < &
                1.0e-10_num) THEN
              grid_root_temp_over_num(iArray, jArray, kArray) = 0.0_num
            ELSE IF (grid_temperature_electron(iArray, jArray, kArray) < &
                1.0e-10_num) THEN
              grid_root_temp_over_num(iArray, jArray, kArray) = 0.0_num
            ELSE
              grid_root_temp_over_num(iArray, jArray, kArray) = &
                  SQRT(grid_temperature_electron(iArray, jArray, kArray) &
                  / grid_num_density_electron(iArray, jArray, kArray))
            END IF
          END DO
        END DO
      END DO

      CALL field_bc(grid_root_temp_over_num, ng)

      DEALLOCATE(grid_num_density_electron_temp)
      DEALLOCATE(grid_temperature_electron_temp)
      DEALLOCATE(grid_num_density_electron)
      DEALLOCATE(grid_temperature_electron)

    ELSE

      ! This will neglect thermal contributions to the bremsstrahlung cross
      ! section
      plasma_factor = 1.0_num

    END IF

    ! Calculate the number density of each ion species
    DO iZ = 1, n_species

      ! Identify if the charge is greater than 1
      Z_temp = species_list(iZ)%atomic_no

      IF (Z_temp < 1 .OR. Z_temp > 100) THEN
        CYCLE
      ENDIF
      CALL calc_number_density(grid_num_density_ion,iZ)
      CALL field_bc(grid_num_density_ion, ng)

      ! Update the optical depth for each electron species
      DO ispecies = 1, n_species

        ! Only update optical_depth_bremsstrahlung for the electron species
        IF (species_list(ispecies)%species_type == c_species_id_electron) THEN

          ! Cycle through all electrons in this species
          current => species_list(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(current))

            ! Get electron energy
            part_ux = current%part_p(1) / mc0
            part_uy = current%part_p(2) / mc0
            part_uz = current%part_p(3) / mc0
            gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
            part_E = gamma_rel*m0*c**2
            current%particle_energy = part_E

            ! Don't update the optical depth if the particle hasn't moved
            IF (gamma_rel - 1.0_num < 1.0e-15_num) THEN
              current => current%next
              CYCLE
            END IF

            ! Get electron speed
            part_v = SQRT(current%part_p(1)**2 + current%part_p(2)**2 &
                + current%part_p(3)**2) * c**2 / part_E

            ! Get number density at electron
            part_x = current%part_pos(1) - x_grid_min_local
            part_y = current%part_pos(2) - y_grid_min_local
            part_z = current%part_pos(3) - z_grid_min_local
            CALL grid_centred_var_at_particle(part_x, part_y, part_z, part_ni,&
                iZ, grid_num_density_ion)

            ! Update the optical depth for the screening option chosen
            IF (use_plasma_screening) THEN

              ! Obtain extra parameters needed for plasma screening model
              CALL grid_centred_var_at_particle(part_x, part_y, part_z, &
                  part_root_Te_over_ne, iZ, grid_root_temp_over_num)

              plasma_factor = get_plasma_factor( &
                  NINT(species_list(iZ)%charge/q0), Z_temp, &
                  part_root_Te_over_ne)

            END IF

            current%optical_depth_bremsstrahlung = &
                current%optical_depth_bremsstrahlung &
                - delta_optical_depth(Z_temp, part_E, part_v, part_ni, &
                plasma_factor)

            ! If optical depth dropped below zero generate photon and reset
            ! optical depth
            IF (current%optical_depth_bremsstrahlung <= 0.0_num) THEN
              CALL generate_photon(current, Z_temp, &
                  bremsstrahlung_photon_species)
              current%optical_depth_bremsstrahlung = reset_optical_depth()
            END IF

            current => current%next

          END DO
        END IF
      END DO
    END DO

    IF (use_plasma_screening) THEN
      DEALLOCATE(grid_root_temp_over_num)
    END IF

  END SUBROUTINE bremsstrahlung_update_optical_depth



  ! Calculate the change in optical depth during this timestep, given by:
  ! (ion number density) * (emission cross section) * (distance traversed)
  ! Here, cross sections are determined from Geant4 look-up tables
  FUNCTION delta_optical_depth(Z, part_E, part_v, part_ni, plasma_factor)

    INTEGER, INTENT(in) :: Z
    REAL(num), INTENT(in) :: part_E, part_v, part_ni, plasma_factor
    INTEGER :: brem_index
    REAL(num) :: cross_sec_val
    REAL(num) :: delta_optical_depth

    brem_index = Z_to_index(Z)

    cross_sec_val = find_value_from_table_1d(part_E, &
        brem_array(brem_index)%size_T, brem_array(brem_index)%E_table, &
        brem_array(brem_index)%cross_section) * plasma_factor

    delta_optical_depth = part_ni * cross_sec_val * part_v * dt &
        / photon_weight

  END FUNCTION delta_optical_depth



  ! Calculate the enhancement of the cross section due to plasma screening
  ! Z: Charge of ion species
  ! A: Atomic number of ion species
  ! part_root_Te_over_ne: Number density of background electron species at the
  !                       electron divided by the temperature of this species,
  !                       evaluated at the electron position.
  FUNCTION get_plasma_factor(Z, A, part_root_Te_over_ne)

    INTEGER, INTENT(in) :: A, Z
    REAL(num), INTENT(in) :: part_root_Te_over_ne
    REAL(num) :: A_third, term1, term2
    REAL(num) :: get_plasma_factor

    A_third = A**(1.0_num/3.0_num)
    term1 = LOG(plasma_screen_const_1/A_third)
    term2 = LOG(plasma_screen_const_2 * part_root_Te_over_ne * A_third)
    get_plasma_factor = 1.0_num &
       + ((REAL(Z,num)/REAL(A,num))**2 * term2 / term1)
    get_plasma_factor = MAX(1.0_num, get_plasma_factor)

  END FUNCTION get_plasma_factor



  ! Draws a new random number for the exponentially distributed optical depths
  FUNCTION reset_optical_depth()

    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth



  ! Generates a photon moving in same direction as the electron, calculating
  ! electron recoil if appropriate
  SUBROUTINE generate_photon(electron, Z, iphoton)

    TYPE(particle), POINTER :: electron
    INTEGER, INTENT(IN) :: iphoton, Z
    INTEGER :: brem_index
    REAL(num) :: dir_x, dir_y, dir_z, mag_p
    REAL(num) :: rand_temp, photon_energy, part_E
    TYPE(particle), POINTER :: new_photon

    ! Obtain electron direction (magnitude must be > 0 to prevent 1/0 issues)
    part_E = electron%particle_energy
    mag_p = MAX(SQRT(electron%part_p(1)**2 + electron%part_p(2)**2 &
        + electron%part_p(3)**2), c_tiny)
    dir_x = electron%part_p(1) / mag_p
    dir_y = electron%part_p(2) / mag_p
    dir_z = electron%part_p(3) / mag_p

    ! Determine photon energy
    rand_temp = random()
    brem_index = Z_to_index(Z)
    photon_energy = find_value_from_table_alt(part_E, rand_temp, &
        brem_array(brem_index)%size_T, brem_array(brem_index)%size_k, &
        brem_array(brem_index)%E_table, brem_array(brem_index)%k_table, &
        brem_array(brem_index)%cdf_table)

    ! Calculate electron recoil
    IF (use_bremsstrahlung_recoil) THEN
      mag_p = mag_p - photon_weight * photon_energy / c
      electron%part_p(1) = dir_x * mag_p
      electron%part_p(2) = dir_y * mag_p
      electron%part_p(3) = dir_z * mag_p
      electron%particle_energy = electron%particle_energy - photon_energy
    END IF

    ! This will only create photons that have energies above a user specified
    ! cutoff and if photon generation is turned on.
    IF (photon_energy > photon_energy_min_bremsstrahlung .AND. &
        produce_bremsstrahlung_photons) THEN

      ! Ensure photon_energy is a number we can handle at our precision
      IF (photon_energy < c_tiny) photon_energy = c_tiny

      ! Create new photon at the electron position, in the electron direction
      CALL create_particle(new_photon)
      new_photon%part_pos = electron%part_pos
      new_photon%part_p(1) = dir_x * photon_energy / c
      new_photon%part_p(2) = dir_y * photon_energy / c
      new_photon%part_p(3) = dir_z * photon_energy / c
#ifdef PHOTONS
      new_photon%optical_depth = reset_optical_depth()
#endif
      new_photon%particle_energy = photon_energy
      new_photon%weight = electron%weight * photon_weight

      CALL add_particle_to_partlist(species_list(iphoton)%attached_list, &
          new_photon)

    END IF

  END SUBROUTINE generate_photon



  ! Calculates the value of a grid-centred variable part_var stored in the grid
  ! grid_var, averaged over the particle shape for a particle at position
  ! (part_x, part_y, part_z) and of species current_species
  SUBROUTINE grid_centred_var_at_particle(part_x, part_y, part_z, part_var, &
      current_species, grid_var)

    REAL(num), INTENT(in) :: part_x, part_y, part_z
    INTEGER, INTENT(in) :: current_species
    REAL(num), INTENT(in) :: grid_var(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng)
    REAL(num), INTENT(out) :: part_var
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
    part_var = fac*part_var

  END SUBROUTINE grid_centred_var_at_particle



  ! For a pair of arrays, x and values, of size nx, this function returns the
  ! interpolated value of "values" corresponding to x_in in the x array. This
  ! uses linear interpolation, unlike in photons.F90 which is logarithmic
  FUNCTION find_value_from_table_1d(x_in, nx, x, values)

    REAL(num) :: find_value_from_table_1d
    REAL(num), INTENT(IN) :: x_in
    INTEGER, INTENT(IN) :: nx
    REAL(num), INTENT(IN) :: x(nx), values(nx)
    REAL(num) :: fx, x_value, value_interp, xdif1, xdif2, xdifm
    INTEGER :: i1, i2, im
    LOGICAL, SAVE :: warning = .TRUE.

    x_value = x_in

    ! Use bisection to find the nearest cells i1, i2
    i1 = 1
    i2 = nx
    xdif1 = x(i1) - x_value
    xdif2 = x(i2) - x_value
    IF (xdif1 * xdif2 < 0) THEN
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_value
        IF (xdif1 * xdifm < 0) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO

      ! Interpolate in x to find fraction between the cells
      fx = (x_value - x(i1)) / (x(i2) - x(i1))

    ! Our x_in value falls outside of the x array - truncate the value
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_1d" in bremsstrahlung.F90', &
            ' outside the range of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
    END IF

    ! Corresponding number from value array, a fraction fx between i1 and i2
    value_interp = (1.0_num - fx) * values(i1) + fx * values(i2)
    find_value_from_table_1d = value_interp

  END FUNCTION find_value_from_table_1d



  ! For each element of x, we have a 1D array of y values and a 1D array of P
  ! values, such that the 1D array x has a corresponding 2D array of y and P.
  ! The P values represent the cumulative probability of finding a particular y
  ! value for a given x. This function returns a y value for a given probability
  ! p_value, at a particular x value x_in, and is used in the code to obtain
  ! photon energies (y) for a given incident electron energy (x).
  FUNCTION find_value_from_table_alt(x_in, p_value, nx, ny, x, y, p_table)

    REAL(num) :: find_value_from_table_alt
    REAL(num), INTENT(IN) :: x_in, p_value
    INTEGER, INTENT(IN) :: nx, ny
    REAL(num), INTENT(IN) :: x(nx), y(nx,ny), p_table(nx,ny)
    INTEGER :: ix, index_lt, index_gt, i1, i2, im
    REAL(num) :: fx, fp, y_lt, y_gt, y_interp, xdif1, xdif2, xdifm
    LOGICAL, SAVE :: warning = .TRUE.

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
      ! Interpolate in x
      fx = (x_in - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fx = 0.0_num
      ELSE
        fx = 1.0_num
      END IF
    END IF

    index_lt = i1
    index_gt = i2

    ix = index_lt
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
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_lt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ix = index_gt
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
      ! Interpolate in x
      fp = (p_value - p_table(ix,i1)) / (p_table(ix,i2) - p_table(ix,i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'Argument to "find_value_from_table_alt" outside the range ', &
            'of the table.'
        PRINT*,'Using truncated value. No more warnings will be issued.'
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0) THEN
        fp = 0.0_num
      ELSE
        fp = 1.0_num
      END IF
    END IF

    y_gt = (1.0_num - fp) * y(ix,i1) + fp * y(ix,i2)

    ! Interpolate in x
    y_interp = (1.0_num - fx) * y_lt + fx * y_gt

    find_value_from_table_alt = y_interp

  END FUNCTION find_value_from_table_alt



  ! Shuts down the bremsstrahlung module - to keep consistency with photons.F90
  SUBROUTINE shutdown_bremsstrahlung_module

    CALL deallocate_tables_bremsstrahlung

  END SUBROUTINE

#endif
END MODULE bremsstrahlung
