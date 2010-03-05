MODULE particles

  USE boundary
  USE shape_functions
  !USE current_smooth

  IMPLICIT NONE

CONTAINS

  SUBROUTINE push_particles

    IMPLICIT NONE

    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_x3

    ! Xi (space factor see page 38 in manual)
    REAL(num), ALLOCATABLE, DIMENSION(:) :: xi0x
    REAL(num), ALLOCATABLE, DIMENSION(:) :: xi1x

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    REAL(num), ALLOCATABLE, DIMENSION(:) :: jxh, jyh, jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_px, part_py, part_pz, part_q, part_m
    REAL(num) :: root, part_vx, part_vy, part_vz, part_weight

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
    REAL(num), DIMENSION(-2:2) :: gx

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * hmx + F(j) * h0x + F(j+1) * hpx
    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(-2:2) :: hx

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+ and P- from Page27 of manual
    REAL(num) :: pxp, pxm, pyp, pym, pzp, pzm

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    ! Tau variables from Page27 of manual
    REAL(num) :: tau, taux, tauy, tauz

    ! Used by J update
    INTEGER :: xmin, xmax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: mean
    INTEGER :: ispecies

    TYPE(particle), POINTER :: current, next

    ALLOCATE(xi0x(-3:3))
    ALLOCATE(xi1x(-3:3))

    ALLOCATE(jxh(-4:3))
    ALLOCATE(jyh(-3:3))
    ALLOCATE(jzh(-3:3))

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    jxh = 0.0_num
    jyh = 0.0_num
    jzh = 0.0_num

    ekbar_sum = 0.0_num
    ct = 0.0_num

    part_weight = weight

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
      DO ipart = 1, particle_species(ispecies)%attached_list%count
        next=>current%next
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
#endif
        ! Set the weighting functions to zero for each new particle
        xi0x = 0.0_num
        xi1x = 0.0_num

        ! Copy the particle properties out for speed
        part_x  = current%part_pos - x_start_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        ! Use a lookup table for charge and mass to SAVE memory
        ! No reason not to do this (I think), check properly later
#ifdef PER_PARTICLE_CHARGEMASS
        part_q = current%charge
        part_m = current%mass
#else
        part_q  = particle_species(ispecies)%charge
        part_m  = particle_species(ispecies)%mass
#endif

#ifdef PARTICLE_PROBES
        init_part_x = current%part_pos
#endif

        ! Calculate v(t+0.5dt) from p(t)
        ! See PSC manual page (25-27)
        root = 1.0_num &
            / SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
        part_vx = part_px * root
        part_vy = part_py * root
        part_vz = part_pz * root

        ! Move particles to half timestep position to first order
        part_x = part_x + part_vx*dt/2.0_num

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
        cell_x_r = part_x/dx
        ! Round cell position to nearest cell
        cell_x1 = NINT(cell_x_r)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1+1

        ! These are now the weighting factors correct for field weighting
        CALL grid_to_particle(cell_frac_x, gx)

        ! particle weighting factors in 1D
        ! These wieght particle properties onto grid
        ! This is used later to calculate J
        CALL particle_to_grid(cell_frac_x, xi0x(-2:2))

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x_r = part_x/dx - 0.5_num
        cell_x2  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r
        cell_x2 = cell_x2+1

        ! Grid weighting factors in 3D (3D analogue of equation 4.77 page 25
        ! of manual). These weight grid properties onto particles
        CALL grid_to_particle(cell_frac_x, hx)

        ex_part = 0.0_num
        ey_part = 0.0_num
        ez_part = 0.0_num
        bx_part = 0.0_num
        by_part = 0.0_num
        bz_part = 0.0_num

        DO ix = -sf_order, sf_order
          ex_part = ex_part + hx(ix)*ex(cell_x2+ix)
          ey_part = ey_part + gx(ix)*ey(cell_x1+ix)
          ez_part = ez_part + gx(ix)*ez(cell_x1+ix)

          bx_part = bx_part + gx(ix)*bx(cell_x1+ix)
          by_part = by_part + hx(ix)*by(cell_x2+ix)
          bz_part = bz_part + hx(ix)*bz(cell_x2+ix)
        ENDDO

        ! update particle momenta using weighted fields
        cmratio = part_q * 0.5_num * dt
        pxm = part_px + cmratio * ex_part
        pym = part_py + cmratio * ey_part
        pzm = part_pz + cmratio * ez_part

        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon

        root = cmratio / SQRT(part_m**2 + (pxm**2 + pym**2 + pzm**2)/c**2)
        taux = bx_part * root
        tauy = by_part * root
        tauz = bz_part * root

        tau = 1.0_num / (1.0_num + taux**2 + tauy**2 + tauz**2)
        pxp = ((1.0_num+taux*taux-tauy*tauy-tauz*tauz)*pxm &
            + (2.0_num*taux*tauy+2.0_num*tauz)*pym &
            + (2.0_num*taux*tauz-2.0_num*tauy)*pzm)*tau
        pyp = ((2.0_num*taux*tauy-2.0_num*tauz)*pxm &
            + (1.0_num-taux*taux+tauy*tauy-tauz*tauz)*pym &
            + (2.0_num*tauy*tauz+2.0_num*taux)*pzm)*tau
        pzp = ((2.0_num*taux*tauz+2.0_num*tauy)*pxm &
            + (2.0_num*tauy*tauz-2.0_num*taux)*pym &
            + (1.0_num-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau

        ! Rotation over, go to full timestep
        part_px = pxp + cmratio * ex_part
        part_py = pyp + cmratio * ey_part
        part_pz = pzp + cmratio * ez_part

        ! Calculate particle velocity from particle momentum
        root = 1.0_num &
            / SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
        part_vx = part_px * root
        part_vy = part_py * root

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + part_vx * dt/2.0_num

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = part_x + x_start_local
        current%part_p  (1) = part_px
        current%part_p  (2) = part_py
        current%part_p  (3) = part_pz

#ifdef PARTICLE_PROBES
        final_part_x = current%part_pos
#endif

        ! If the code is compiled with tracer particle support then put in an
        ! If statement so that the current is not calculated for this species
#ifdef TRACER_PARTICLES
        IF (.NOT. particle_species(ispecies)%tracer) THEN
#endif

          ! Now advance to t+1.5dt to calculate current. This is detailed in
          ! the manual between pages 37 and 41. The version coded up looks
          ! completely different to that in the manual, but is equivalent.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x = part_x + part_vx * dt/2.0_num

          cell_x_r = part_x / dx
          cell_x3  = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r
          cell_x3 = cell_x3+1

          CALL particle_to_grid(cell_frac_x, &
              xi1x(cell_x3-cell_x1-2:cell_x3-cell_x1+2))

          ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
          ! the current update much simpler
          xi1x = xi1x - xi0x

          ! Remember that due to CFL condition particle can never cross more
          ! than one gridcell in one timestep

          IF (cell_x3 .EQ. cell_x1) THEN
            ! particle is still in same cell at t+1.5dt as at t+0.5dt
            xmin = -sf_order
            xmax = +sf_order
          ELSE IF (cell_x3 .EQ. cell_x1 - 1) THEN
            ! particle has moved one cell to left
            xmin = -sf_order-1
            xmax = +sf_order
          ELSE ! IF (cell_x3 .EQ. cell_x1 + 1) THEN
            ! particle has moved one cell to right
            xmin = -sf_order
            xmax = +sf_order+1
          ENDIF

          ! Set these to zero due to diffential inside loop
          jxh = 0.0_num
          jyh = 0.0_num
          jzh = 0.0_num

          DO ix = xmin, xmax
            wx = xi1x(ix)
            wy = xi0x(ix) + 0.5_num * xi1x(ix)
            wz = xi0x(ix) + 0.5_num * xi1x(ix)

            ! This is the bit that actually solves d(rho)/dt = -div(J)
            jxh(ix) = jxh(ix-1) - part_q * wx * 1.0_num/dt * part_weight
            jyh(ix) = part_q * part_vy * wy  * part_weight/dx
            jzh(ix) = part_q * part_vz * wz  * part_weight/dx

            jx(cell_x1+ix) = jx(cell_x1+ix) +jxh(ix)
            jy(cell_x1+ix) = jy(cell_x1+ix) +jyh(ix)
            jz(cell_x1+ix) = jz(cell_x1+ix) +jzh(ix)
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
          probe_energy = &
              (SQRT(1.0_num + (part_px**2 + part_py**2 + part_pz**2) &
              / (part_m * c)**2) - 1.0_num) * (part_m * c**2)

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
          current_probe=> current_probe%next
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

    DO ispecies = 1, n_species
      CALL processor_summation_bcs(ekbar_sum(:,ispecies))
      CALL field_bc(ekbar_sum(:,ispecies))
      CALL processor_summation_bcs(ct(:,ispecies))
      CALL field_bc(ct(:,ispecies))
    ENDDO

    ! Calculate the mean kinetic energy for each species in space
    ekbar = 0.0_num
    DO ispecies = 1, n_species
      DO ix = 1, nx
        mean = ekbar_sum(ix, ispecies)/MAX(ct(ix, ispecies), c_non_zero)
        ekbar(ix, ispecies) = mean
      ENDDO
    ENDDO

    DEALLOCATE(xi0x)
    DEALLOCATE(xi1x)
    DEALLOCATE(jxh)
    DEALLOCATE(jyh)
    DEALLOCATE(jzh)

    CALL particle_bcs

  END SUBROUTINE push_particles

END MODULE particles
