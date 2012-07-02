MODULE particles

  USE boundary
  USE shape_functions
  USE current_smooth
  USE prefetch

  IMPLICIT NONE

CONTAINS

  SUBROUTINE push_particles

    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_x3

    ! Xi (space factor see page 38 in manual)
    ! The code now uses gx and hx instead of xi0 and xi1

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
    REAL(num), DIMENSION(sf0-2:sf1+1) :: jxh
    REAL(num), DIMENSION(sf0-1:sf1+1) :: jyh
    REAL(num), DIMENSION(sf0-1:sf1+1) :: jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, part_mc, ipart_mc, part_weight

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy, part_mc2
#endif

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: uxp, uxm, uyp, uym, uzp, uzm
    REAL(num) :: tau, taux, tauy, tauz, taux2, tauy2, tauz2

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio, ccmratio

    ! Used by J update
    INTEGER :: xmin, xmax
    REAL(num) :: wx, wy

    ! Temporary variables
    REAL(num) :: idx
    REAL(num) :: idtf, idxf
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: fcx, fcy, fjx, fjy, fjz
    REAL(num) :: root, fac, dtfac, gamma, cf2
    REAL(num) :: delta_x, part_vy, part_vz
    INTEGER :: ispecies, ix, dcellx
    INTEGER(i8) :: ipart
#ifdef PARTICLE_PROBES
    LOGICAL :: probes_for_species
#endif
#ifdef TRACER_PARTICLES
    LOGICAL :: not_tracer_species
#endif

    TYPE(particle), POINTER :: current, next

#ifdef PREFETCH
    CALL prefetch_particle(species_list(1)%attached_list%head)
#endif

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    jxh = 0.0_num
    jyh = 0.0_num
    jzh = 0.0_num

    gx = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    ! particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    fac = 1.0_num / 24.0_num
#elif  PARTICLE_SHAPE_TOPHAT
    fac = 1.0_num
#else
    fac = 0.5_num
#endif
    dtfac = 0.5_num * dt * fac

    idtf = idt * fac
    idxf = idx * fac

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
#ifdef PHOTONS
      IF (species_list(ispecies)%species_type .EQ. c_species_id_photon) THEN
        IF (photon_dynamics) CALL push_photons(ispecies)
        CYCLE
      ENDIF
#endif
#ifdef PARTICLE_PROBES
      current_probe => species_list(ispecies)%attached_probes
      probes_for_species = ASSOCIATED(current_probe)
#endif
#ifdef TRACER_PARTICLES
      not_tracer_species = .NOT. species_list(ispecies)%tracer
#endif

#ifndef PER_PARTICLE_WEIGHT
      part_weight = species_list(ispecies)%weight
      fcx = idtf * part_weight
      fcy = idxf * part_weight
#endif
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q   = species_list(ispecies)%charge
      part_mc  = c * species_list(ispecies)%mass
      ipart_mc = 1.0_num / part_mc
      cmratio  = part_q * dtfac * ipart_mc
      ccmratio = c * cmratio
#ifdef PARTICLE_PROBES
      part_mc2 = c * part_mc
#endif
#endif
      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next=>current%next
#ifdef PREFETCH
        CALL prefetch_particle(next)
#endif
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
        fcx = idtf * part_weight
        fcy = idxf * part_weight
#endif
#ifdef PARTICLE_PROBES
        init_part_x = current%part_pos
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q   = current%charge
        part_mc  = c * current%mass
        ipart_mc = 1.0_num / part_mc
        cmratio  = part_q * dtfac * ipart_mc
        ccmratio = c * cmratio
#ifdef PARTICLE_PROBES
        part_mc2 = c * part_mc
#endif
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos - x_min_local
        part_ux = current%part_p(1) * ipart_mc
        part_uy = current%part_p(2) * ipart_mc
        part_uz = current%part_p(3) * ipart_mc

        ! Calculate v(t) from p(t)
        ! See PSC manual page (25-27)
        root = dtco2 / SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)

        ! Move particles to half timestep position to first order
        part_x = part_x + part_ux * root

        ! Work out the grid cell number for the particle.
        ! Not an integer in general.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
#else
        cell_x_r = part_x * idx
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
        INCLUDE 'include/bspline3/gx.inc'
#elif  PARTICLE_SHAPE_TOPHAT
        INCLUDE 'include/tophat/gx.inc'
#else
        INCLUDE 'include/triangle/gx.inc'
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        dcellx = 0
#ifdef PARTICLE_SHAPE_BSPLINE3
        INCLUDE 'include/bspline3/hx_dcell.inc'
#elif  PARTICLE_SHAPE_TOPHAT
        INCLUDE 'include/tophat/hx_dcell.inc'
#else
        INCLUDE 'include/triangle/hx_dcell.inc'
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
        INCLUDE 'include/bspline3/eb_part.inc'
#elif  PARTICLE_SHAPE_TOPHAT
        INCLUDE 'include/tophat/eb_part.inc'
#else
        INCLUDE 'include/triangle/eb_part.inc'
#endif

        ! update particle momenta using weighted fields
        uxm = part_ux + cmratio * ex_part
        uym = part_uy + cmratio * ey_part
        uzm = part_uz + cmratio * ez_part

        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
        root = ccmratio / SQRT(uxm**2 + uym**2 + uzm**2 + 1.0_num)

        taux = bx_part * root
        tauy = by_part * root
        tauz = bz_part * root

        taux2 = taux**2
        tauy2 = tauy**2
        tauz2 = tauz**2

        tau = 1.0_num / (1.0_num + taux2 + tauy2 + tauz2)

        uxp = ((1.0_num + taux2 - tauy2 - tauz2) * uxm &
            + 2.0_num * ((taux * tauy + tauz) * uym &
            + (taux * tauz - tauy) * uzm)) * tau
        uyp = ((1.0_num - taux2 + tauy2 - tauz2) * uym &
            + 2.0_num * ((tauy * tauz + taux) * uzm &
            + (tauy * taux - tauz) * uxm)) * tau
        uzp = ((1.0_num - taux2 - tauy2 + tauz2) * uzm &
            + 2.0_num * ((tauz * taux + tauy) * uxm &
            + (tauz * tauy - taux) * uym)) * tau

        ! Rotation over, go to full timestep
        part_ux = uxp + cmratio * ex_part
        part_uy = uyp + cmratio * ey_part
        part_uz = uzp + cmratio * ez_part

        ! Calculate particle velocity from particle momentum
        gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        root = c / gamma

        delta_x = part_ux * root * dto2
        part_vy = part_uy * root
        part_vz = part_uz * root

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + delta_x

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = part_x + x_min_local
        current%part_p   = part_mc * (/ part_ux, part_uy, part_uz /)

#ifdef PARTICLE_PROBES
        final_part_x = current%part_pos
#endif
        ! Original code calculates densities of electrons, ions and neutrals
        ! here. This has been removed to reduce memory footprint

        ! If the code is compiled with tracer particle support then put in an
        ! IF statement so that the current is not calculated for this species
#ifdef TRACER_PARTICLES
        IF (not_tracer_species) THEN
#endif
          ! Now advance to t+1.5dt to calculate current. This is detailed in
          ! the manual between pages 37 and 41. The version coded up looks
          ! completely different to that in the manual, but is equivalent.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x = part_x + delta_x

#ifdef PARTICLE_SHAPE_TOPHAT
          cell_x_r = part_x * idx - 0.5_num
#else
          cell_x_r = part_x * idx
#endif
          cell_x3 = FLOOR(cell_x_r + 0.5_num)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r
          cell_x3 = cell_x3 + 1

          hx = 0.0_num

          dcellx = cell_x3 - cell_x1
#ifdef PARTICLE_SHAPE_BSPLINE3
          INCLUDE 'include/bspline3/hx_dcell.inc'
#elif  PARTICLE_SHAPE_TOPHAT
          INCLUDE 'include/tophat/hx_dcell.inc'
#else
          INCLUDE 'include/triangle/hx_dcell.inc'
#endif

          ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
          ! the current update much simpler
          hx = hx - gx

          ! Remember that due to CFL condition particle can never cross more
          ! than one gridcell in one timestep

          xmin = sf_min + (dcellx - 1) / 2
          xmax = sf_max + (dcellx + 1) / 2

          ! Set these to zero due to diffential inside loop
          jxh = 0.0_num

          fjx = fcx * part_q
          fjy = fcy * part_q * part_vy
          fjz = fcy * part_q * part_vz

          DO ix = xmin, xmax
            wx =  hx(ix)
            wy =  gx(ix) + 0.5_num * hx(ix)

            ! This is the bit that actually solves d(rho)/dt = -div(J)
            jxh(ix) = jxh(ix-1) - fjx * wx
            jyh(ix) = fjy * wy
            jzh(ix) = fjz * wy

            jx(cell_x1+ix) = jx(cell_x1+ix) + jxh(ix)
            jy(cell_x1+ix) = jy(cell_x1+ix) + jyh(ix)
            jz(cell_x1+ix) = jz(cell_x1+ix) + jzh(ix)
          ENDDO
#ifdef TRACER_PARTICLES
        ENDIF
#endif
#ifdef PARTICLE_PROBES
        IF (probes_for_species) THEN
          ! Compare the current particle with the parameters of any probes in
          ! the system. These particles are copied into a separate part of the
          ! output file.

          current_probe => species_list(ispecies)%attached_probes

          ! Cycle through probes
          DO WHILE(ASSOCIATED(current_probe))
            ! Note that this is the energy of a single REAL particle in the
            ! pseudoparticle, NOT the energy of the pseudoparticle
            probe_energy = (gamma - 1.0_num) * part_mc2

            ! Unidirectional probe
            IF (probe_energy .GT. current_probe%ek_min) THEN
              IF (probe_energy .LT. current_probe%ek_max) THEN

                d_init  = current_probe%normal &
                    * (current_probe%point - init_part_x)
                d_final = current_probe%normal &
                    * (current_probe%point - final_part_x)
                IF (d_final .LT. 0.0_num .AND. d_init .GE. 0.0_num) THEN
                  ! this particle is wanted so copy it to the list associated
                  ! with this probe
                  ALLOCATE(particle_copy)
                  particle_copy = current
                  CALL add_particle_to_partlist(&
                      current_probe%sampled_particles, particle_copy)
                  NULLIFY(particle_copy)
                ENDIF

              ENDIF
            ENDIF
            current_probe => current_probe%next
          ENDDO
        ENDIF
#endif
        current=>next
      ENDDO
    ENDDO

    CALL current_bcs
    CALL particle_bcs

    IF (smooth_currents) CALL smooth_current()

  END SUBROUTINE push_particles

#ifdef PHOTONS
  SUBROUTINE push_photons(ispecies)

    ! Very simple photon pusher
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: delta_x
    INTEGER,INTENT(IN) :: ispecies
    TYPE(particle), POINTER :: current

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy
    LOGICAL :: probes_for_species
#endif

#ifdef PARTICLE_PROBES
    current_probe => species_list(ispecies)%attached_probes
    probes_for_species = ASSOCIATED(current_probe)
#endif

    ! set current to point to head of list
    current => species_list(ispecies)%attached_list%head
    ! loop over photons
    DO WHILE(ASSOCIATED(current))
      delta_x = current%part_p(1) * dt
#ifdef PARTICLE_PROBES
      init_part_x = current%part_pos
#endif
      current%part_pos = current%part_pos + delta_x
#ifdef PARTICLE_PROBES
      final_part_x = current%part_pos
#endif

#ifdef PARTICLE_PROBES
      IF (probes_for_species) THEN
        ! Compare the current particle with the parameters of any probes in
        ! the system. These particles are copied into a separate part of the
        ! output file.

        current_probe => species_list(ispecies)%attached_probes

        ! Cycle through probes
        DO WHILE(ASSOCIATED(current_probe))
          ! Note that this is the energy of a single REAL particle in the
          ! pseudoparticle, NOT the energy of the pseudoparticle
          probe_energy = current%particle_energy

          ! Unidirectional probe
          IF (probe_energy .GT. current_probe%ek_min) THEN
            IF (probe_energy .LT. current_probe%ek_max) THEN

              d_init  = current_probe%normal &
                  * (current_probe%point - init_part_x)
              d_final = current_probe%normal &
                  * (current_probe%point - final_part_x)
              IF (d_final .LT. 0.0_num .AND. d_init .GE. 0.0_num) THEN
                ! this particle is wanted so copy it to the list associated
                ! with this probe
                ALLOCATE(particle_copy)
                particle_copy = current
                CALL add_particle_to_partlist(&
                    current_probe%sampled_particles, particle_copy)
                NULLIFY(particle_copy)
              ENDIF

            ENDIF
          ENDIF
          current_probe => current_probe%next
        ENDDO
      ENDIF
#endif

      current => current%next
    ENDDO

  END SUBROUTINE push_photons
#endif

END MODULE particles
