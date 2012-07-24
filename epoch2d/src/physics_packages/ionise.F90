MODULE ionise

  USE calc_df
  USE split_particle
  USE shape_functions
  USE random_generator
  USE numerics

  IMPLICIT NONE

  REAL(num), PARAMETER :: ionisation_exponent = -1.0 / 3.0_num
  REAL(num), PARAMETER :: bessel_constant = SQRT(8.0_num / pi)

  REAL(num), DIMENSION(:), ALLOCATABLE :: released_mass_fraction
  REAL(num), DIMENSION(:), ALLOCATABLE :: effective_n_exponent
  REAL(num), DIMENSION(:), ALLOCATABLE :: adk_scaling
  REAL(num), DIMENSION(:), ALLOCATABLE :: ionisation_constant
  REAL(num), DIMENSION(:), ALLOCATABLE :: smallest_e_mag
  REAL(num), DIMENSION(:), ALLOCATABLE :: bsi_constant
  REAL(num), DIMENSION(:), ALLOCATABLE :: bsi_scaling
  REAL(num), DIMENSION(:), ALLOCATABLE :: bsi_threshold
  REAL(num), DIMENSION(:), ALLOCATABLE :: adk_maximum
  REAL(num), DIMENSION(:), ALLOCATABLE :: adk_bsi_cap
  REAL(num), DIMENSION(:), ALLOCATABLE :: keldysh
  REAL(num), DIMENSION(:), ALLOCATABLE :: multi_constant
  REAL(num), DIMENSION(:), ALLOCATABLE :: k_photons_energy
  REAL(num), DIMENSION(:), ALLOCATABLE :: k_photons_exponent
  REAL(num), DIMENSION(:), ALLOCATABLE :: adk_multiphoton_cap

CONTAINS

  ! Do most of the calculations on a per-species basis once at start of
  ! simulation to save on computational time

  SUBROUTINE initialise_ionisation

    INTEGER :: i, bessel_error

    ALLOCATE(released_mass_fraction(n_species), &
        effective_n_exponent(n_species), adk_scaling(n_species), &
        ionisation_constant(n_species), smallest_e_mag(n_species))
    IF (use_bsi) ALLOCATE(bsi_constant(n_species), bsi_scaling(n_species), &
        bsi_threshold(n_species), adk_maximum(n_species), &
        adk_bsi_cap(n_species))
    IF (use_multiphoton) ALLOCATE(keldysh(n_species), &
        multi_constant(n_species), k_photons_energy(n_species), &
        k_photons_exponent(n_species), adk_multiphoton_cap(n_species))

    DO i = 1, n_species
      IF (species_list(i)%ionise) THEN
        ! Mass fraction of released electron
        released_mass_fraction(i) = &
            species_list(species_list(i)%release_species)%mass &
            / species_list(i)%mass
        ! Effective principle quantum number
        effective_n_exponent(i) = &
            (species_list(species_list(i)%ionise_to_species)%charge / ev) &
            / SQRT(2.0_num * species_list(i)%ionisation_energy / hartree)
        ! Electric field strength scaling in ADK
        adk_scaling(i) = 2.0 * SQRT((2.0_num &
            * species_list(i)%ionisation_energy / hartree)**3)
        ! Constant in ADK equation
        ionisation_constant(i) = SQRT(6.0_num / pi) &
            * species_list(i)%ionisation_energy / hartree * 2.0_num**(2.0_num &
            * effective_n_exponent(i)) / (effective_n_exponent(i) &
            * gamma_fn(2.0_num * effective_n_exponent(i)))

        IF (use_bsi) THEN
          ! Constant in BSI equation
          bsi_constant(i) = (species_list(i)%ionisation_energy / hartree) &
              * SQRT(2.0_num * species_list(i)%ionisation_energy / hartree) &
              / (pi * species_list(species_list(i)%ionise_to_species)%charge &
              / ev)
          ! Electric field strength scaling in BSI
          bsi_scaling(i) = (species_list(i)%ionisation_energy / hartree)**2 &
              / (4.0_num &
              * species_list(species_list(i)%ionise_to_species)%charge / ev)
          ! Transition to BSI electric field strength
          bsi_threshold(i) = &
              (species_list(species_list(i)%ionise_to_species)%charge / ev)**3 &
              / (2.0_num * effective_n_exponent(i))**4
          ! Turning point in ADK rate
          adk_maximum(i) = MAX(adk_scaling(i) / (3.0_num * (2.0_num &
              * effective_n_exponent(i) + species_list(i)%l - 1.5_num)), &
              bsi_threshold(i))
        ENDIF

        IF (use_multiphoton) THEN
          ! Number of photons required for ionisation in multiphoton
          k_photons_exponent(i) = &
              REAL(FLOOR((species_list(i)%ionisation_energy / hartree) &
              / (laser_x_min%omega*atomic_time)) + 1, num)
          multi_constant(i) = factorial(k_photons_exponent(i))

          ! If K! is too large then the multiphoton ionisation rate is zero
          IF (multi_constant(i) .LT. SQRT(HUGE(0.0_num))) THEN
            multi_constant(i) = c * (atomic_time / a0) * multi_constant(i)**2 &
                * REAL(species_list(i)%n**5,num) &
                * (laser_x_min%omega * atomic_time)**((10.0_num &
                * k_photons_exponent(i) - 1.0_num) / 3.0_num) &
                * SQRT(k_photons_exponent(i)) &
                * (2.0_num * k_photons_exponent(i) + 1.0_num)
          ELSE
            multi_constant(i) = 0.0_num
          ENDIF

          ! Constant in multiphoton equations, calculated like this to trap any
          ! floating underflow
          IF (multi_constant(i) .NE. 0.0_num) multi_constant(i) = 4.8_num &
              * (1.3_num * c * (atomic_time / a0) / (8.0_num * pi &
              * laser_x_min%omega * atomic_time))**k_photons_exponent(i) &
              / multi_constant(i)
          ! Energy in K photons
          k_photons_energy(i) = &
              k_photons_exponent(i) * h_bar * laser_x_min%omega
          ! Constant used in multiphoton equation
          k_photons_exponent(i) = 4.0_num * k_photons_exponent(i) - 2.0_num
        ENDIF

        ! Constant used in ADK equation
        effective_n_exponent(i) = 2.0_num * effective_n_exponent(i) - 1.5_num
        IF (use_bsi) THEN
          ! Threshold ADK ionisation rate at transition to BSI regime
          adk_bsi_cap(i) = ionisation_constant(i) * (adk_scaling(i) &
              / bsi_threshold(i))**effective_n_exponent(i) &
              * EXP(ionisation_exponent * adk_scaling(i) / bsi_threshold(i)) &
              * (bessel_constant * SQRT(adk_scaling(i) / bsi_threshold(i)) &
              * EXP(adk_scaling(i) / bsi_threshold(i)) * RKBESL(adk_scaling(i) &
              / bsi_threshold(i), 0.5_num, species_list(i)%l + 1, 1, &
              bessel_error) - 1.0_num)
        ENDIF

        IF (use_multiphoton) THEN
          ! The transition between multiphoton and tunnelling is controlled by
          ! this parameter, if E > keldysh then we are in the tunnelling regime.
          ! Where the ionisation energy is exceptionally large, the smallest E
          ! we can use in the tunnelling equations can be larger than the E at
          ! which we define the multiphoton to tunnelling transition. In this
          ! case we ensure a monotonic increasing ionisation rate by forcing the
          ! transition to occur at the smallest usable E instead.
          keldysh(i) = MAX(adk_scaling(i) / (0.99472065388909858_num &
              * c_largest_exp), laser_x_min%omega * atomic_time * SQRT(2.0_num &
              * species_list(i)%ionisation_energy / hartree) / 0.5)

          ! Threshold ADK ionisation rate at transition to multiphoton regime
          ! Required to ensure the rate is always monotonic increasing
          adk_multiphoton_cap(i) = ionisation_constant(i) * (adk_scaling(i) &
              / keldysh(i))**effective_n_exponent(i) * EXP(ionisation_exponent &
              * adk_scaling(i) / keldysh(i)) * (bessel_constant &
              * SQRT(adk_scaling(i) / keldysh(i)) * EXP(adk_scaling(i) &
              / keldysh(i)) * RKBESL(adk_scaling(i) / keldysh(i), 0.5_num, &
              species_list(i)%l + 1, 1, bessel_error) - 1.0_num)

          ! If multiphoton ionisation rate will always be zero for ion, then set
          ! the smallest electric field strength for which the rate will be
          ! calculated to the largest machine number. Otherwise calculate the
          ! smallest electric field strength for which multiphoton can be
          ! calculated
          IF (multi_constant(i) .EQ. 0.0_num) THEN
            smallest_e_mag(i) = c_largest_number
          ELSE
            smallest_e_mag(i) = (TINY(0.0_num) &
            / MIN(multi_constant(i), 1.0_num))**(1.0 / k_photons_exponent(i))
          ENDIF
        ELSE
          ! The smallest E we can use in the ADK equation is numerically
          ! determined and so may not be portable; if over/underflow errors
          ! start occurring this might require a smarter solution
          smallest_e_mag(i) = adk_scaling(i) / (0.99472065388909858_num &
              * c_largest_exp)
        ENDIF
      ELSE
        released_mass_fraction(i) = 0.0_num
        effective_n_exponent(i) = 0.0_num
        adk_scaling(i) = 0.0_num
        ionisation_constant(i) = 0.0_num
        smallest_e_mag(i) = 0.0_num
        IF (use_bsi) THEN
          bsi_constant(i) = 0.0_num
          bsi_scaling(i) = 0.0_num
          bsi_threshold(i) = 0.0_num
          adk_maximum(i) = 0.0_num
          adk_bsi_cap(i) = 0.0_num
        ENDIF
        IF (use_multiphoton) THEN
          keldysh(i) = 0.0_num
          multi_constant(i) = 0.0_num
          k_photons_energy(i) = 0.0_num
          k_photons_exponent(i) = 0.0_num
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE initialise_ionisation



  ! Redirect to correct routine; the only difference between these routines is
  ! in the rate calculated but it saves doing this logic for each ionisation

  SUBROUTINE ionise_particles

    ! Barrier suppression model from "Molecular dissociative ionisation using a
    ! classical over-the-barrier approach" Posthumus et al
    ! IOP Conference Series 1997 Vol. 154 Pg. 298-307
    IF (use_bsi) THEN
      ! Multiphoton model from "Multiphoton Processes in Atoms" Delone et al
      ! Springer-Verlag 2000 Vol. 13 Pg. 97
      IF (use_multiphoton) THEN
        CALL multiphoton_tunnelling_bsi
      ELSE
        CALL tunnelling_bsi
      ENDIF
    ELSE
      IF (use_multiphoton) THEN
        CALL multiphoton_tunnelling
      ELSE
        CALL tunnelling
        ! Tunnelling model from "Tunnel ionization of complex atoms and of
        ! atomic ions in an alternating electromagnetic field" Ammosov et al
        ! Soviet Physics JETP 1986 Vol. 64 Pg. 1191-1194
      ENDIF
    ENDIF

  END SUBROUTINE ionise_particles



  SUBROUTINE multiphoton_tunnelling_bsi

    INTEGER :: i, current_state, bessel_error
    INTEGER :: ix, cell_x1, cell_x2, dcellx
    INTEGER :: iy, cell_y1, cell_y2, dcelly
    REAL(num) :: rate, ex_part, ey_part, ez_part, e_part_mag, time_left, sample
    REAL(num) :: cf2, j_ion(3)
    REAL(num) :: gx(sf_min:sf_max), hx(sf_min:sf_max)
    REAL(num) :: gy(sf_min:sf_max), hy(sf_min:sf_max)
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x, idx
    REAL(num) :: part_y, part_y2, cell_y_r, cell_frac_y, idy
    LOGICAL :: multiphoton_ionised

    TYPE(particle), POINTER :: current, new, next
    TYPE(particle_list) :: ionised_list(n_species)

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    idx = 1.0_num / dx
    idy = 1.0_num / dy

    ! Stores ionised species until close of ionisation run. Main purpose of this
    ! method is to ensure proper statistics (i.e. prevent ionisation rate being
    ! calculated for full dt when particle has already ionised in time step)
    ! and to stop electric field at particle being calculated more than once
    DO i = 1, n_species
      CALL create_empty_partlist(ionised_list(i))
    ENDDO

    ! Ionise a species at a time
    DO i = 1, n_species
      ! Skip particle if it cannot be ionised
      IF ( .NOT. species_list(i)%ionise) CYCLE
      ! Start with first particle in the list
      current => species_list(i)%attached_list%head

      ! Try to ionise every particle of the species
      DO WHILE(ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        dcellx = 0
        dcelly = 0
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
        ! This can be done with electric field smoothing but hasn't been
        ! necessary since the statistics were changed away from using number
        ! densities.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#else
#include "triangle/e_part.inc"
#endif

        ! Electric field strength in atomic units
        e_part_mag = SQRT(ex_part**2 + ey_part**2 + ez_part**2) &
            / atomic_electric_field
        next => current%next
        ! Need to keep track of what ionisation state current particle has
        ! moved to
        current_state = i
        time_left = dt / atomic_time
        j_ion = 0.0_num

        ! This cycles through every ionisation level for the particle until it
        ! is no longer ionising in the field
        DO WHILE(time_left .GT. 0.0_num &
            .AND. species_list(current_state)%ionise)
          ! If the electron is ionised in the multiphoton regime, the resultant
          ! velocity is different so we need to track this
          multiphoton_ionised = .FALSE.
          ! If we're past the maximum of the ADK rate then we're definitely in
          ! the BSI region...
          IF (e_part_mag .GT. adk_maximum(current_state)) THEN
            rate = bsi_constant(current_state) * (1.0_num &
                - bsi_scaling(current_state) / e_part_mag) &
                + adk_bsi_cap(current_state)
          ! ... otherwise we need to check if we're in the tunnelling or
          ! multiphoton regime...
          ELSEIF (e_part_mag .GT. keldysh(current_state)) THEN
            ! Calculate ADK ionisation rate
            rate = ionisation_constant(current_state) &
                * (adk_scaling(current_state) &
                / e_part_mag)**effective_n_exponent(current_state) &
                * EXP(ionisation_exponent * adk_scaling(current_state) &
                / e_part_mag) * (bessel_constant &
                * SQRT(adk_scaling(current_state) / e_part_mag) &
                * EXP(adk_scaling(current_state) / e_part_mag) &
                * RKBESL(adk_scaling(current_state) / e_part_mag, 0.5_num, &
                species_list(current_state)%l + 1, 1, bessel_error) - 1.0_num)
            ! If we're in the BSI regime, choose the smallest of either the ADK
            ! or BSI rate. This allows a smooth transition to BSI
            IF (e_part_mag .GT. bsi_threshold(current_state)) &
                rate = MIN(bsi_constant(current_state) &
                * (1.0_num - bsi_scaling(current_state) / e_part_mag) &
                + adk_bsi_cap(current_state), rate)
          ! If we're in the multiphoton regime, make sure the electric field
          ! strength is larger than the minimum value for multiphoton
          ELSEIF (e_part_mag .GT. smallest_e_mag(current_state)) THEN
            rate = MIN(adk_multiphoton_cap(current_state), &
                multi_constant(current_state) &
                * e_part_mag**k_photons_exponent(current_state))
            multiphoton_ionised = .TRUE.
          ELSE
            ! If we got here then the electric field strength was too small for
            ! any ionisation
            EXIT
          ENDIF

          sample = random()
          ! Calculate probability of ionisation using a cumulative distribution
          ! function modelling ionisation in a field as an exponential decay
          IF (sample .LT. 1.0_num - exp(-1.0_num * rate * time_left)) THEN
            IF (species_list(current_state)%release_species .GT. 0) THEN
              ALLOCATE(new)
              ! Create electron for release
#ifdef PER_PARTICLE_WEIGHT
              new%weight = current%weight
#endif
              new%part_pos = current%part_pos
              ! Electron is released without acceleration so simply use momentum
              ! conservation to split the particle
              new%part_p = current%part_p &
                  * released_mass_fraction(current_state)
              current%part_p = current%part_p - new%part_p
              ! Using multiphoton ionisation, the additional energy from photons
              ! accelerates the electron in the direction of the electric field.
              ! This is an approximation as the ejection angle ranges widely
              ! with a maxima at theta = 0 with respect to the field
              IF (multiphoton_ionised) new%part_p = new%part_p + SQRT(2.0_num &
                  * m0 * (k_photons_energy(current_state) &
                  - species_list(current_state)%ionisation_energy)) &
                  * (/ ex_part, ey_part, ez_part /) / (e_part_mag &
                  * atomic_electric_field)
#ifdef PER_PARTICLE_CHARGE_MASS
              new%charge = species_list( &
                  species_list(current_state)%release_species)%charge
              new%mass = species_list( &
                  species_list(current_state)%release_species)%mass
              current%charge = current%charge - new%charge
              current%mass = current%mass - new%mass
#endif
#ifdef PARTICLE_DEBUG
              new%processor = rank
              new%processor_at_t0 = rank
#endif
              ! Put electron into particle lists
              CALL add_particle_to_partlist(species_list(species_list( &
                  current_state)%release_species)%attached_list, new)
            ENDIF
            ! Calculates the time of ionisation using inverse sampling, and
            ! subtracts it from the time step. Ensures diminishing time for
            ! successive ionisations
            time_left = time_left + log(1.0_num - sample) / rate
            ! Current correction as proposed from Mulser et al 1998, true from
            ! ejection energy <e_j> << m_e*c**2, i.e. sub-relativistic ejection
            ! velocity. This shall be true for all laser gamma factors, as BSI
            ! techniques release electron at rest in zero field approximation,
            ! and BSI encompasses both tunneling and over-barrier ionisation
            ! rates. Multiphoton rates will only be used for low intensity
            ! lasers
            IF (multiphoton_ionised) THEN
              j_ion = j_ion + k_photons_energy(current_state)
            ELSE
              j_ion = j_ion + species_list(current_state)%ionisation_energy
            ENDIF
            current_state = species_list(current_state)%ionise_to_species
          ELSE
            time_left = 0.0_num
          ENDIF
        ENDDO

        ! Finally the ion is moved to the ionised list following multiple
        ! ionisation, and current correction is applied
        IF (current_state .NE. i) THEN
          CALL remove_particle_from_partlist(species_list(i)%attached_list, &
              current)
          CALL add_particle_to_partlist(ionised_list(current_state), current)

          j_ion = fac * j_ion * current%weight &
              * (/ ex_part, ey_part, ez_part /) &
              / (dt * (atomic_electric_field * e_part_mag)**2)

          IF (j_ion(1) .NE. 0.0_num .OR. j_ion(2) .NE. 0.0_num .OR. &
              j_ion(3) .NE. 0.0_num) THEN
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              jx(cell_x2+ix, cell_y1+iy) = &
                  jx(cell_x2+ix, cell_y1+iy) + hx(ix) * gy(iy) * j_ion(1)
              jy(cell_x1+ix, cell_y2+iy) = &
                  jy(cell_x1+ix, cell_y2+iy) + gx(ix) * hy(iy) * j_ion(2)
              jz(cell_x1+ix, cell_y1+iy) = &
                  jz(cell_x1+ix, cell_y1+iy) + gx(ix) * gy(iy) * j_ion(3)
            ENDDO
            ENDDO
          ENDIF
        ENDIF
        current => next
      ENDDO
    ENDDO

    ! Clean up procedure; put ionised ions back into the correct particle lists
    DO i = 1, n_species
      CALL append_partlist(species_list(i)%attached_list, ionised_list(i))
    ENDDO
    ! Put ionised particles back into partlists

  END SUBROUTINE multiphoton_tunnelling_bsi



  SUBROUTINE multiphoton_tunnelling

    INTEGER :: i, current_state, bessel_error
    INTEGER :: ix, cell_x1, cell_x2, dcellx
    INTEGER :: iy, cell_y1, cell_y2, dcelly
    REAL(num) :: rate, ex_part, ey_part, ez_part, e_part_mag, time_left, sample
    REAL(num) :: cf2, j_ion(3)
    REAL(num) :: gx(sf_min:sf_max), hx(sf_min:sf_max)
    REAL(num) :: gy(sf_min:sf_max), hy(sf_min:sf_max)
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x, idx
    REAL(num) :: part_y, part_y2, cell_y_r, cell_frac_y, idy
    LOGICAL :: multiphoton_ionised

    TYPE(particle), POINTER :: current, new, next
    TYPE(particle_list) :: ionised_list(n_species)

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    idx = 1.0_num / dx
    idy = 1.0_num / dy

    ! Stores ionised species until close of ionisation run. Main purpose of this
    ! method is to ensure proper statistics (i.e. prevent ionisation rate being
    ! calculated for full dt when particle has already ionised in time step)
    ! and to stop electric field at particle being calculated more than once
    DO i = 1, n_species
      CALL create_empty_partlist(ionised_list(i))
    ENDDO

    ! Ionise a species at a time
    DO i = 1, n_species
      ! Skip particle if it cannot be ionised
      IF ( .NOT. species_list(i)%ionise) CYCLE
      ! Start with first particle in the list
      current => species_list(i)%attached_list%head

      ! Try to ionise every particle of the species
      DO WHILE(ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        dcellx = 0
        dcelly = 0
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
        ! This can be done with electric field smoothing but hasn't been
        ! necessary since the statistics were changed away from using number
        ! densities.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#else
#include "triangle/e_part.inc"
#endif

        ! Electric field strength in atomic units
        e_part_mag = SQRT(ex_part**2 + ey_part**2 + ez_part**2) &
            / atomic_electric_field
        next => current%next
        ! Need to keep track of what ionisation state current particle has
        ! moved to
        current_state = i
        time_left = dt / atomic_time
        j_ion = 0.0_num

        ! This cycles through every ionisation level for the particle until it
        ! is no longer ionising in the field
        DO WHILE(time_left .GT. 0.0_num &
            .AND. species_list(current_state)%ionise)
          ! If the electron is ionised in the multiphoton regime, the resultant
          ! velocity is different so we need to track this
          multiphoton_ionised = .FALSE.
          ! Check if we're in the tunnelling or multiphoton regime...
          IF (e_part_mag .GT. keldysh(current_state)) THEN
            ! Calculate ADK ionisation rate
            rate = ionisation_constant(current_state) &
                * (adk_scaling(current_state) &
                / e_part_mag)**effective_n_exponent(current_state) &
                * EXP(ionisation_exponent * adk_scaling(current_state) &
                / e_part_mag) * (bessel_constant &
                * SQRT(adk_scaling(current_state) / e_part_mag) &
                * EXP(adk_scaling(current_state) / e_part_mag) &
                * RKBESL(adk_scaling(current_state) / e_part_mag, 0.5_num, &
                species_list(current_state)%l + 1, 1, bessel_error) - 1.0_num)
          ! If we're in the multiphoton regime, make sure the electric field
          ! strength is larger than the minimum value for multiphoton
          ELSEIF (e_part_mag .GT. smallest_e_mag(current_state)) THEN
            rate = MIN(adk_multiphoton_cap(current_state), &
                multi_constant(current_state) &
                * e_part_mag**k_photons_exponent(current_state))
            multiphoton_ionised = .TRUE.
          ELSE
            ! If we got here then the electric field strength was too small for
            ! any ionisation
            EXIT
          ENDIF

          sample = random()
          ! Calculate probability of ionisation using a cumulative distribution
          ! function modelling ionisation in a field as an exponential decay
          IF (sample .LT. 1.0_num - exp(-1.0_num * rate * time_left)) THEN
            IF (species_list(current_state)%release_species .GT. 0) THEN
              ALLOCATE(new)
              ! Create electron for release
#ifdef PER_PARTICLE_WEIGHT
              new%weight = current%weight
#endif
              new%part_pos = current%part_pos
              ! Electron is released without acceleration so simply use momentum
              ! conservation to split the particle
              new%part_p = current%part_p &
                  * released_mass_fraction(current_state)
              current%part_p = current%part_p - new%part_p
              ! Using multiphoton ionisation, the additional energy from photons
              ! accelerates the electron in the direction of the electric field.
              ! This is an approximation as the ejection angle ranges widely
              ! with a maxima at theta = 0 with respect to the field
              IF (multiphoton_ionised) new%part_p = new%part_p + SQRT(2.0_num &
                  * m0 * (k_photons_energy(current_state) &
                  - species_list(current_state)%ionisation_energy)) &
                  * (/ ex_part, ey_part, ez_part /) / (e_part_mag &
                  * atomic_electric_field)
#ifdef PER_PARTICLE_CHARGE_MASS
              new%charge = species_list( &
                  species_list(current_state)%release_species)%charge
              new%mass = species_list( &
                  species_list(current_state)%release_species)%mass
              current%charge = current%charge - new%charge
              current%mass = current%mass - new%mass
#endif
#ifdef PARTICLE_DEBUG
              new%processor = rank
              new%processor_at_t0 = rank
#endif
              ! Put electron into particle lists
              CALL add_particle_to_partlist(species_list(species_list( &
                  current_state)%release_species)%attached_list, new)
            ENDIF
            ! Calculates the time of ionisation using inverse sampling, and
            ! subtracts it from the time step. Ensures diminishing time for
            ! successive ionisations
            time_left = time_left + log(1.0_num - sample) / rate
            ! Current correction as proposed from Mulser et al 1998, true from
            ! ejection energy <e_j> << m_e*c**2, i.e. sub-relativistic ejection
            ! velocity. This shall be true for all laser gamma factors, as BSI
            ! techniques release electron at rest in zero field approximation,
            ! and BSI encompasses both tunneling and over-barrier ionisation
            ! rates. Multiphoton rates will only be used for low intensity
            ! lasers
            IF (multiphoton_ionised) THEN
              j_ion = j_ion + k_photons_energy(current_state)
            ELSE
              j_ion = j_ion + species_list(current_state)%ionisation_energy
            ENDIF
            current_state = species_list(current_state)%ionise_to_species
          ELSE
            time_left = 0.0_num
          ENDIF
        ENDDO

        ! Finally the ion is moved to the ionised list following multiple
        ! ionisation, and current correction is applied
        IF (current_state .NE. i) THEN
          CALL remove_particle_from_partlist(species_list(i)%attached_list, &
              current)
          CALL add_particle_to_partlist(ionised_list(current_state), current)

          j_ion = fac * j_ion * current%weight &
              * (/ ex_part, ey_part, ez_part /) &
              / (dt * (atomic_electric_field * e_part_mag)**2)

          IF (j_ion(1) .NE. 0.0_num .OR. j_ion(2) .NE. 0.0_num .OR. &
              j_ion(3) .NE. 0.0_num) THEN
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              jx(cell_x2+ix, cell_y1+iy) = &
                  jx(cell_x2+ix, cell_y1+iy) + hx(ix) * gy(iy) * j_ion(1)
              jy(cell_x1+ix, cell_y2+iy) = &
                  jy(cell_x1+ix, cell_y2+iy) + gx(ix) * hy(iy) * j_ion(2)
              jz(cell_x1+ix, cell_y1+iy) = &
                  jz(cell_x1+ix, cell_y1+iy) + gx(ix) * gy(iy) * j_ion(3)
            ENDDO
            ENDDO
          ENDIF
        ENDIF
        current => next
      ENDDO
    ENDDO

    ! Clean up procedure; put ionised ions back into the correct particle lists
    DO i = 1, n_species
      CALL append_partlist(species_list(i)%attached_list, ionised_list(i))
    ENDDO
    ! Put ionised particles back into partlists

  END SUBROUTINE multiphoton_tunnelling



  SUBROUTINE tunnelling_bsi

    INTEGER :: i, current_state, bessel_error
    INTEGER :: ix, cell_x1, cell_x2, dcellx
    INTEGER :: iy, cell_y1, cell_y2, dcelly
    REAL(num) :: rate, ex_part, ey_part, ez_part, e_part_mag, time_left, sample
    REAL(num) :: cf2, j_ion(3)
    REAL(num) :: gx(sf_min:sf_max), hx(sf_min:sf_max)
    REAL(num) :: gy(sf_min:sf_max), hy(sf_min:sf_max)
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x, idx
    REAL(num) :: part_y, part_y2, cell_y_r, cell_frac_y, idy

    TYPE(particle), POINTER :: current, new, next
    TYPE(particle_list) :: ionised_list(n_species)

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    idx = 1.0_num / dx
    idy = 1.0_num / dy

    ! Stores ionised species until close of ionisation run. Main purpose of this
    ! method is to ensure proper statistics (i.e. prevent ionisation rate being
    ! calculated for full dt when particle has already ionised in time step)
    ! and to stop electric field at particle being calculated more than once
    DO i = 1, n_species
      CALL create_empty_partlist(ionised_list(i))
    ENDDO

    ! Ionise a species at a time
    DO i = 1, n_species
      ! Skip particle if it cannot be ionised
      IF ( .NOT. species_list(i)%ionise) CYCLE
      ! Start with first particle in the list
      current => species_list(i)%attached_list%head

      ! Try to ionise every particle of the species
      DO WHILE(ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        dcellx = 0
        dcelly = 0
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
        ! This can be done with electric field smoothing but hasn't been
        ! necessary since the statistics were changed away from using number
        ! densities.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#else
#include "triangle/e_part.inc"
#endif

        ! Electric field strength in atomic units
        e_part_mag = SQRT(ex_part**2 + ey_part**2 + ez_part**2) &
            / atomic_electric_field
        next => current%next
        ! Need to keep track of what ionisation state current particle has
        ! moved to
        current_state = i
        time_left = dt / atomic_time
        j_ion = 0.0_num

        ! This cycles through every ionisation level for the particle until it
        ! is no longer ionising in the field
        DO WHILE(time_left .GT. 0.0_num &
            .AND. species_list(current_state)%ionise)
          ! If we're past the maximum of the ADK rate then we're definitely in
          ! the BSI region...
          IF (e_part_mag .GT. adk_maximum(current_state)) THEN
            rate = bsi_constant(current_state) * (1.0_num &
                - bsi_scaling(current_state) / e_part_mag) &
                + adk_bsi_cap(current_state)
          ! ... otherwise we check if the electric field is large enough for
          ! tunnelling ionisation...
          ELSEIF (e_part_mag .GT. smallest_e_mag(current_state)) THEN
            ! Calculate ADK ionisation rate
            rate = ionisation_constant(current_state) &
                * (adk_scaling(current_state) &
                / e_part_mag)**effective_n_exponent(current_state) &
                * EXP(ionisation_exponent * adk_scaling(current_state) &
                / e_part_mag) * (bessel_constant &
                * SQRT(adk_scaling(current_state) / e_part_mag) &
                * EXP(adk_scaling(current_state) / e_part_mag) &
                * RKBESL(adk_scaling(current_state) / e_part_mag, 0.5_num, &
                species_list(current_state)%l + 1, 1, bessel_error) - 1.0_num)
            ! ... but we could still be in the lower part of the in the BSI
            ! regime, choose the smallest of either the ADK or BSI rate.
            ! This allows a smooth transition to BSI
            IF (e_part_mag .GT. bsi_threshold(current_state)) &
                rate = MIN(bsi_constant(current_state) &
                * (1.0_num - bsi_scaling(current_state) / e_part_mag) &
                + adk_bsi_cap(current_state), rate)
          ELSE
            ! If we got here then the electric field strength was too small for
            ! any ionisation
            EXIT
          ENDIF

          sample = random()
          ! Calculate probability of ionisation using a cumulative distribution
          ! function modelling ionisation in a field as an exponential decay
          IF (sample .LT. 1.0_num - exp(-1.0_num * rate * time_left)) THEN
            IF (species_list(current_state)%release_species .GT. 0) THEN
              ALLOCATE(new)
              ! Create electron for release
#ifdef PER_PARTICLE_WEIGHT
              new%weight = current%weight
#endif
              new%part_pos = current%part_pos
              ! Electron is released without acceleration so simply use momentum
              ! conservation to split the particle
              new%part_p = current%part_p &
                  * released_mass_fraction(current_state)
              current%part_p = current%part_p - new%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
              new%charge = species_list( &
                  species_list(current_state)%release_species)%charge
              new%mass = species_list( &
                  species_list(current_state)%release_species)%mass
              current%charge = current%charge - new%charge
              current%mass = current%mass - new%mass
#endif
#ifdef PARTICLE_DEBUG
              new%processor = rank
              new%processor_at_t0 = rank
#endif
              ! Put electron into particle lists
              CALL add_particle_to_partlist(species_list(species_list( &
                  current_state)%release_species)%attached_list, new)
            ENDIF
            ! Calculates the time of ionisation using inverse sampling, and
            ! subtracts it from the time step. Ensures diminishing time for
            ! successive ionisations
            time_left = time_left + log(1.0_num - sample) / rate
            ! Current correction as proposed from Mulser et al 1998, true from
            ! ejection energy <e_j> << m_e*c**2, i.e. sub-relativistic ejection
            ! velocity. This shall be true for all laser gamma factors, as BSI
            ! techniques release electron at rest in zero field approximation,
            ! and BSI encompasses both tunneling and over-barrier ionisation
            ! rates
            j_ion = j_ion + species_list(current_state)%ionisation_energy
            current_state = species_list(current_state)%ionise_to_species
          ELSE
            time_left = 0.0_num
          ENDIF
        ENDDO

        ! Finally the ion is moved to the ionised list following multiple
        ! ionisation, and current correction is applied
        IF (current_state .NE. i) THEN
          CALL remove_particle_from_partlist(species_list(i)%attached_list, &
              current)
          CALL add_particle_to_partlist(ionised_list(current_state), current)

          j_ion = fac * j_ion * current%weight &
              * (/ ex_part, ey_part, ez_part /) &
              / (dt * (atomic_electric_field * e_part_mag)**2)

          IF (j_ion(1) .NE. 0.0_num .OR. j_ion(2) .NE. 0.0_num .OR. &
              j_ion(3) .NE. 0.0_num) THEN
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              jx(cell_x2+ix, cell_y1+iy) = &
                  jx(cell_x2+ix, cell_y1+iy) + hx(ix) * gy(iy) * j_ion(1)
              jy(cell_x1+ix, cell_y2+iy) = &
                  jy(cell_x1+ix, cell_y2+iy) + gx(ix) * hy(iy) * j_ion(2)
              jz(cell_x1+ix, cell_y1+iy) = &
                  jz(cell_x1+ix, cell_y1+iy) + gx(ix) * gy(iy) * j_ion(3)
            ENDDO
            ENDDO
          ENDIF
        ENDIF
        current => next
      ENDDO
    ENDDO

    ! Clean up procedure; put ionised ions back into the correct particle lists
    DO i = 1, n_species
      CALL append_partlist(species_list(i)%attached_list, ionised_list(i))
    ENDDO
    ! Put ionised particles back into partlists

  END SUBROUTINE tunnelling_bsi



  SUBROUTINE tunnelling

    INTEGER :: i, current_state, bessel_error
    INTEGER :: ix, cell_x1, cell_x2, dcellx
    INTEGER :: iy, cell_y1, cell_y2, dcelly
    REAL(num) :: rate, ex_part, ey_part, ez_part, e_part_mag, time_left, sample
    REAL(num) :: cf2, j_ion(3)
    REAL(num) :: gx(sf_min:sf_max), hx(sf_min:sf_max)
    REAL(num) :: gy(sf_min:sf_max), hy(sf_min:sf_max)
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x, idx
    REAL(num) :: part_y, part_y2, cell_y_r, cell_frac_y, idy

    TYPE(particle), POINTER :: current, new, next
    TYPE(particle_list) :: ionised_list(n_species)

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    idx = 1.0_num / dx
    idy = 1.0_num / dy

    ! Stores ionised species until close of ionisation run. Main purpose of this
    ! method is to ensure proper statistics (i.e. prevent ionisation rate being
    ! calculated for full dt when particle has already ionised in time step)
    ! and to stop electric field at particle being calculated more than once
    DO i = 1, n_species
      CALL create_empty_partlist(ionised_list(i))
    ENDDO

    ! Ionise a species at a time
    DO i = 1, n_species
      ! Skip particle if it cannot be ionised
      IF ( .NOT. species_list(i)%ionise) CYCLE
      ! Start with first particle in the list
      current => species_list(i)%attached_list%head

      ! Try to ionise every particle of the species
      DO WHILE(ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        dcellx = 0
        dcelly = 0
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
        ! This can be done with electric field smoothing but hasn't been
        ! necessary since the statistics were changed away from using number
        ! densities.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#else
#include "triangle/e_part.inc"
#endif

        ! Electric field strength in atomic units
        e_part_mag = SQRT(ex_part**2 + ey_part**2 + ez_part**2) &
            / atomic_electric_field
        next => current%next
        ! Need to keep track of what ionisation state current particle has
        ! moved to
        current_state = i
        time_left = dt / atomic_time
        j_ion = 0.0_num

        ! This cycles through every ionisation level for the particle until it
        ! is no longer ionising in the field
        DO WHILE(time_left .GT. 0.0_num &
            .AND. species_list(current_state)%ionise)
          ! Check if the electric field is strong enough for tunnelling
          ! ionisation
          IF (e_part_mag .GT. smallest_e_mag(current_state)) THEN
            rate = ionisation_constant(current_state) &
                * (adk_scaling(current_state) &
                / e_part_mag)**effective_n_exponent(current_state) &
                * EXP(ionisation_exponent * adk_scaling(current_state) &
                / e_part_mag) * (bessel_constant &
                * SQRT(adk_scaling(current_state) / e_part_mag) &
                * EXP(adk_scaling(current_state) / e_part_mag) &
                * RKBESL(adk_scaling(current_state) / e_part_mag, 0.5_num, &
                species_list(current_state)%l + 1, 1, bessel_error) - 1.0_num)
          ELSE
            ! If we got here then the electric field strength was too small for
            ! any ionisation
            EXIT
          ENDIF

          sample = random()
          ! Calculate probability of ionisation using a cumulative distribution
          ! function modelling ionisation in a field as an exponential decay
          IF (sample .LT. 1.0_num - exp(-1.0_num * rate * time_left)) THEN
            IF (species_list(current_state)%release_species .GT. 0) THEN
              ALLOCATE(new)
              ! Create electron for release
#ifdef PER_PARTICLE_WEIGHT
              new%weight = current%weight
#endif
              new%part_pos = current%part_pos
              ! Electron is released without acceleration so simply use momentum
              ! conservation to split the particle
              new%part_p = current%part_p &
                  * released_mass_fraction(current_state)
              current%part_p = current%part_p - new%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
              new%charge = species_list( &
                  species_list(current_state)%release_species)%charge
              new%mass = species_list( &
                  species_list(current_state)%release_species)%mass
              current%charge = current%charge - new%charge
              current%mass = current%mass - new%mass
#endif
#ifdef PARTICLE_DEBUG
              new%processor = rank
              new%processor_at_t0 = rank
#endif
              ! Put electron into particle lists
              CALL add_particle_to_partlist(species_list(species_list( &
                  current_state)%release_species)%attached_list, new)
            ENDIF
            ! Calculates the time of ionisation using inverse sampling, and
            ! subtracts it from the time step. Ensures diminishing time for
            ! successive ionisations
            time_left = time_left + log(1.0_num - sample) / rate
            ! Current correction as proposed from Mulser et al 1998, true from
            ! ejection energy <e_j> << m_e*c**2, i.e. sub-relativistic ejection
            ! velocity. This shall be true for all laser gamma factors, as BSI
            ! techniques release electron at rest in zero field approximation,
            ! and BSI encompasses both tunneling and over-barrier ionisation
            ! rates
            j_ion = j_ion + species_list(current_state)%ionisation_energy
            current_state = species_list(current_state)%ionise_to_species
          ELSE
            time_left = 0.0_num
          ENDIF
        ENDDO

        ! Finally the ion is moved to the ionised list following multiple
        ! ionisation, and current correction is applied
        IF (current_state .NE. i) THEN
          CALL remove_particle_from_partlist(species_list(i)%attached_list, &
              current)
          CALL add_particle_to_partlist(ionised_list(current_state), current)

          j_ion = fac * j_ion * current%weight &
              * (/ ex_part, ey_part, ez_part /) &
              / (dt * (atomic_electric_field * e_part_mag)**2)

          IF (j_ion(1) .NE. 0.0_num .OR. j_ion(2) .NE. 0.0_num .OR. &
              j_ion(3) .NE. 0.0_num) THEN
            DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              jx(cell_x2+ix, cell_y1+iy) = &
                  jx(cell_x2+ix, cell_y1+iy) + hx(ix) * gy(iy) * j_ion(1)
              jy(cell_x1+ix, cell_y2+iy) = &
                  jy(cell_x1+ix, cell_y2+iy) + gx(ix) * hy(iy) * j_ion(2)
              jz(cell_x1+ix, cell_y1+iy) = &
                  jz(cell_x1+ix, cell_y1+iy) + gx(ix) * gy(iy) * j_ion(3)
            ENDDO
            ENDDO
          ENDIF
        ENDIF
        current => next
      ENDDO
    ENDDO

    ! Clean up procedure; put ionised ions back into the correct particle lists
    DO i = 1, n_species
      CALL append_partlist(species_list(i)%attached_list, ionised_list(i))
    ENDDO
    ! Put ionised particles back into partlists

  END SUBROUTINE tunnelling

END MODULE ionise
