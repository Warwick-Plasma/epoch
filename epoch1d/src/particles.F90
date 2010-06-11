MODULE particles

  USE boundary
  USE shape_functions
  !USE current_smooth

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

    INTEGER, PARAMETER :: sf = sf_order

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    REAL(num), DIMENSION(-sf-2:sf+1) :: jxh
    REAL(num), DIMENSION(-sf-1:sf+1) :: jyh
    REAL(num), DIMENSION(-sf-1:sf+1) :: jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x
    REAL(num) :: part_px, part_py, part_pz
    REAL(num) :: part_q, part_mc, part_weight

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: probe_energy
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
    REAL(num), DIMENSION(-sf-1:sf+1) :: gx

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(-sf-1:sf+1) :: hx

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: pxp, pxm, pyp, pym, pzp, pzm
    REAL(num) :: tau, taux, tauy, tauz

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    ! Used by J update
    INTEGER :: xmin, xmax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: idx
    REAL(num) :: idxf, idtf
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: fcx, fcy, fjx, fjy, fjz
    REAL(num) :: root, mean, fac, dtfac, gamma_mass_c, cf2
    REAL(num) :: delta_x, part_vy, part_vz
    INTEGER :: ispecies, ix, dcellx
    INTEGER(KIND=8) :: ipart

    TYPE(particle), POINTER :: current, next

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

    part_weight = weight
    fcx = idtf * part_weight
    fcy = idxf * part_weight

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGEMASS
      part_q  = particle_species(ispecies)%charge
      part_mc = c * particle_species(ispecies)%mass
#endif
      !DEC$ VECTOR ALWAYS
      DO ipart = 1, particle_species(ispecies)%attached_list%count
        next=>current%next
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
        fcx = idtf * part_weight
        fcy = idxf * part_weight
#endif
#ifdef PARTICLE_PROBES
        init_part_x = current%part_pos
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos - x_min_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        ! Use a lookup table for charge and mass to save memory
        ! No reason not to do this (I think), check properly later
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
        part_mc = c * current%mass
#endif

        ! Calculate v(t+0.5dt) from p(t)
        ! See PSC manual page (25-27)
        root = dtco2 / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)

        ! Move particles to half timestep position to first order
        part_x = part_x + part_px * root

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
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
        cmratio = part_q * dtfac
        pxm = part_px + cmratio * ex_part
        pym = part_py + cmratio * ey_part
        pzm = part_pz + cmratio * ez_part

        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
        root = c * cmratio / SQRT(part_mc**2 + pxm**2 + pym**2 + pzm**2)

        taux = bx_part * root
        tauy = by_part * root
        tauz = bz_part * root

        tau = 1.0_num / (1.0_num + taux**2 + tauy**2 + tauz**2)

        pxp = ((1.0_num + taux**2 - tauy**2 - tauz**2) * pxm &
            + 2.0_num * ((taux * tauy + tauz) * pym &
            + (taux * tauz - tauy) * pzm)) * tau
        pyp = ((1.0_num - taux**2 + tauy**2 - tauz**2) * pym &
            + 2.0_num * ((tauy * tauz + taux) * pzm &
            + (tauy * taux - tauz) * pxm)) * tau
        pzp = ((1.0_num - taux**2 - tauy**2 + tauz**2) * pzm &
            + 2.0_num * ((tauz * taux + tauy) * pxm &
            + (tauz * tauy - taux) * pym)) * tau

        ! Rotation over, go to full timestep
        part_px = pxp + cmratio * ex_part
        part_py = pyp + cmratio * ey_part
        part_pz = pzp + cmratio * ez_part

        ! Calculate particle velocity from particle momentum
        gamma_mass_c = SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)
        root = c / gamma_mass_c

        delta_x = part_px * root * dto2
        part_vy = part_py * root
        part_vz = part_pz * root

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + delta_x

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = part_x + x_min_local
        current%part_p   = (/ part_px, part_py, part_pz /)

#ifdef PARTICLE_PROBES
        final_part_x = current%part_pos
#endif
        ! Original code calculates densities of electrons, ions and neutrals
        ! here. This has been removed to reduce memory footprint

        ! If the code is compiled with tracer particle support then put in an
        ! IF statement so that the current is not calculated for this species
#ifdef TRACER_PARTICLES
        IF (.NOT. particle_species(ispecies)%tracer) THEN
#endif
          ! Now advance to t+1.5dt to calculate current. This is detailed in
          ! the manual between pages 37 and 41. The version coded up looks
          ! completely different to that in the manual, but is equivalent.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x = part_x + delta_x

#ifdef PARTICLE_SHAPE_TOPHAT
          cell_x_r = part_x / dx - 0.5_num
#else
          cell_x_r = part_x / dx
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

          xmin = -sf_order + (dcellx - 1) / 2
          xmax =  sf_order + (dcellx + 1) / 2

          ! Set these to zero due to diffential inside loop
          jxh = 0.0_num

          fjx = fcx * part_q
          fjy = fcy * part_q * part_vy
          fjz = fcy * part_q * part_vz

          DO ix = xmin, xmax
            wx = hx(ix)
            wy = gx(ix) + 0.5_num * hx(ix)
            wz = gx(ix) + 0.5_num * hx(ix)

            ! This is the bit that actually solves d(rho)/dt = -div(J)
            jxh(ix) = jxh(ix-1) - fjx * wx
            jyh(ix) = fjy * wy
            jzh(ix) = fjz * wz

            jx(cell_x1+ix) = &
                jx(cell_x1+ix) + jxh(ix)
            jy(cell_x1+ix) = &
                jy(cell_x1+ix) + jyh(ix)
            jz(cell_x1+ix) = &
                jz(cell_x1+ix) + jzh(ix)
          ENDDO
#ifdef TRACER_PARTICLES
        ENDIF
#endif
#ifdef PARTICLE_PROBES
        ! Compare the current particle with the parameters of any probes in the
        ! system. These particles are copied into a separate part of the output
        ! file.

        current_probe=>particle_species(ispecies)%attached_probes

        ! Cycle through probes
        DO WHILE(ASSOCIATED(current_probe))
          ! Note that this is the energy of a single REAL particle in the
          ! pseudoparticle, NOT the energy of the pseudoparticle
          probe_energy = c * (gamma_mass_c - part_mc)

          ! right energy? (in J)
          IF (probe_energy .GT. current_probe%ek_min) THEN
            IF ((probe_energy .LT. current_probe%ek_max) &
                .OR. (current_probe%ek_max .LT. 0.0_num)) THEN

              IF (current_probe%left_to_right) THEN
                IF (init_part_x .LT. current_probe%probe_point &
                    .AND. final_part_x .GT. current_probe%probe_point) THEN
                  ! this particle is wanted so copy it to the list associated
                  ! with this probe
                  ALLOCATE(particle_copy)
                  particle_copy = current
                  CALL add_particle_to_partlist(&
                      current_probe%sampled_particles, particle_copy)
                  NULLIFY(particle_copy)
                ENDIF
              ELSE
                IF (init_part_x .GT. current_probe%probe_point &
                    .AND. final_part_x .LT. current_probe%probe_point) THEN
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
          ENDIF
          current_probe=>current_probe%next
        ENDDO
#endif
        current=>next
      ENDDO
    ENDDO

    ! domain is decomposed. Just add currents at edges
    CALL processor_summation_bcs(jx)
    CALL field_bc(jx)
    CALL processor_summation_bcs(jy)
    CALL field_bc(jy)
    CALL processor_summation_bcs(jz)
    CALL field_bc(jz)

    CALL particle_bcs

  END SUBROUTINE push_particles

END MODULE particles
