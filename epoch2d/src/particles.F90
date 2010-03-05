MODULE particles

  USE boundary
  USE shape_functions
  USE current_smooth

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
    INTEGER :: cell_x1, cell_x2, cell_x3, cell_y1, cell_y2, cell_y3

    ! Xi (space factor see page 38 in manual)
    REAL(num), ALLOCATABLE, DIMENSION(:) :: xi0x, xi0y
    REAL(num), ALLOCATABLE, DIMENSION(:) :: xi1x, xi1y
    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    REAL(num), ALLOCATABLE, DIMENSION(:,:) :: jxh, jyh, jzh

    ! Timestep used in calculating current update
    ! dt_j = dt for relativistic code, not for non-relativistic
    REAL(num) :: dt_j

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_px, part_py, part_pz, part_q, part_m
    REAL(num) :: root, part_vx, part_vy, part_vz, part_weight

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, init_part_y, final_part_x, final_part_y
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: probe_x1, probe_y1, probe_x2, probe_y2, probe_a, probe_b
    REAL(num) :: probe_temp, probe_energy
#endif

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy

    ! particle weight factors as described in the manual (FIXREF)
    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(-2:2) :: hx, hy

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+ and P- from Boris1970
    REAL(num) :: pxp, pxm, pyp, pym, pzp, pzm

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    ! Tau variables from Boris1970
    REAL(num) :: tau, taux, tauy, tauz

    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: mean
    INTEGER :: ispecies

    TYPE(particle), POINTER :: current, next

    ALLOCATE(xi0x(-3:3), xi0y(-3:3))
    ALLOCATE(xi1x(-3:3), xi1y(-3:3))

    ALLOCATE(jxh(-4:3, -3:3))
    ALLOCATE(jyh(-3:4, -4:3))
    ALLOCATE(jzh(-3:3, -3:3))

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    jxh = 0.0_num
    jyh = 0.0_num
    jzh = 0.0_num

    ekbar_sum = 0.0_num
    ct = 0.0_num

    part_weight = weight

#ifdef NEWTONIAN
    max_part_v = 0.0_num
#endif

    dt_j = dt

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
        xi0y = 0.0_num
        xi1y = 0.0_num

        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local
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
        init_part_x = current%part_pos(1)
        init_part_y = current%part_pos(2)
#endif

        ! Calculate v(t+0.5dt) from p(t)
        ! See PSC manual page (25-27)
#ifdef NEWTONIAN
        root = 1.0_num / part_m
#else
        root = 1.0_num &
            / SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
#endif
        part_vx = part_px * root
        part_vy = part_py * root
        part_vz = part_pz * root

        ! Move particles to half timestep position to first order
        part_x = part_x + part_vx*dt/2.0_num
        part_y = part_y + part_vy*dt/2.0_num

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
        cell_x_r = part_x/dx
        ! Round cell position to nearest cell
        cell_x1 = NINT(cell_x_r)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1+1

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
        cell_y_r = part_y/dy
        ! Round cell position to nearest cell
        cell_y1 = NINT(cell_y_r)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1+1

        ! particle weight factors as described in the manual (FIXREF)
        ! These weight grid properties onto particles
        CALL grid_to_particle(cell_frac_x, gx)
        CALL grid_to_particle(cell_frac_y, gy)

        ! particle weight factors as described in the manual (FIXREF)
        ! These wieght particle properties onto grid
        ! This is used later to calculate J

        CALL particle_to_grid(cell_frac_x, xi0x(-2:2))
        CALL particle_to_grid(cell_frac_y, xi0y(-2:2))

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x_r = part_x/dx - 0.5_num
        cell_x2  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r
        cell_x2 = cell_x2+1

        cell_y_r = part_y/dy - 0.5_num
        cell_y2  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r
        cell_y2 = cell_y2+1

        CALL grid_to_particle(cell_frac_x, hx)
        CALL grid_to_particle(cell_frac_y, hy)

        ! These are the electric an magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.

        ex_part = 0.0_num
        ey_part = 0.0_num
        ez_part = 0.0_num
        bx_part = 0.0_num
        by_part = 0.0_num
        bz_part = 0.0_num

        DO ix = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            ex_part = ex_part + hx(ix)*gy(iy)*ex(cell_x2+ix, cell_y1+iy)
            ey_part = ey_part + gx(ix)*hy(iy)*ey(cell_x1+ix, cell_y2+iy)
            ez_part = ez_part + gx(ix)*gy(iy)*ez(cell_x1+ix, cell_y1+iy)

            bx_part = bx_part + gx(ix)*hy(iy)*bx(cell_x1+ix, cell_y2+iy)
            by_part = by_part + hx(ix)*gy(iy)*by(cell_x2+ix, cell_y1+iy)
            bz_part = bz_part + hx(ix)*hy(iy)*bz(cell_x2+ix, cell_y2+iy)
          ENDDO
        ENDDO

        ! update particle momenta using weighted fields
        cmratio = part_q * 0.5_num * dt
        pxm = part_px + cmratio * ex_part
        pym = part_py + cmratio * ey_part
        pzm = part_pz + cmratio * ez_part

        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
#ifdef NEWTONIAN
        root = cmratio / part_m
#else
        root = cmratio / SQRT(part_m**2 + (pxm**2 + pym**2 + pzm**2)/c**2)
#endif
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
#ifdef NEWTONIAN
        root = 1.0_num / part_m
#else
        root = 1.0_num &
            / SQRT(part_m**2 + (part_px**2 + part_py**2 + part_pz**2)/c**2)
#endif
        part_vx = part_px * root
        part_vy = part_py * root
        part_vz = part_pz * root

#ifdef NEWTONIAN
        IF (ABS(part_vx) .GT. max_part_v) max_part_v = ABS(part_vx)
        IF (ABS(part_vy) .GT. max_part_v) max_part_v = ABS(part_vy)
#endif

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + part_vx * dt/2.0_num
        part_y = part_y + part_vy * dt/2.0_num

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos(1) = part_x + x_start_local
        current%part_pos(2) = part_y + y_start_local
        current%part_p  (1) = part_px
        current%part_p  (2) = part_py
        current%part_p  (3) = part_pz

#ifdef PARTICLE_PROBES
        final_part_x = current%part_pos(1)
        final_part_y = current%part_pos(2)
#endif

        ! If the code is compiled with tracer particle support then put in an
        ! If statement so that the current is not calculated for this species
#ifdef TRACER_PARTICLES
        IF (.NOT. particle_species(ispecies)%tracer) THEN
#endif

#ifdef NEWTONIAN
          dt_j = 0.05_num * MIN(MIN(dx, dy)/MAX(ABS(part_vx), ABS(part_vy)), dt)
#endif

          ! Now advance to t+1.5dt to calculate current.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x = part_x + part_vx * dt_j/2.0_num
          part_y = part_y + part_vy * dt_j/2.0_num

          cell_x_r = part_x / dx
          cell_x3  = NINT(cell_x_r)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r
          cell_x3 = cell_x3+1

          cell_y_r = part_y / dy
          cell_y3  = NINT(cell_y_r)
          cell_frac_y = REAL(cell_y3, num) - cell_y_r
          cell_y3 = cell_y3+1

          CALL particle_to_grid(cell_frac_x, &
              xi1x(cell_x3-cell_x1-2:cell_x3-cell_x1+2))
          CALL particle_to_grid(cell_frac_y, &
              xi1y(cell_y3-cell_y1-2:cell_y3-cell_y1+2))

          ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
          ! the current update much simpler
          xi1x = xi1x - xi0x
          xi1y = xi1y - xi0y

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
            xmax = sf_order+1
          ENDIF

          IF (cell_y3 .EQ. cell_y1) THEN
            ! particle is still in same cell at t+1.5dt as at t+0.5dt
            ymin = -sf_order
            ymax = +sf_order
          ELSE IF (cell_y3 .EQ. cell_y1 - 1) THEN
            ! particle has moved one cell to left
            ymin = -sf_order-1
            ymax = +sf_order
          ELSE ! IF (cell_y3 .EQ. cell_y1 + 1) THEN
            ! particle has moved one cell to right
            ymin = -sf_order
            ymax = +sf_order+1
          ENDIF

          ! Set these to zero due to diffential inside loop
          jxh = 0.0_num
          jyh = 0.0_num
          jzh = 0.0_num

          DO iy = ymin, ymax
            DO ix = xmin, xmax
              wx = xi1x(ix) * (xi0y(iy) + 0.5_num * xi1y(iy))
              wy = xi1y(iy) * (xi0x(ix) + 0.5_num * xi1x(ix))
              wz = xi0x(ix) * xi0y(iy) + 0.5_num*xi1x(ix)*xi0y(iy) &
                  + 0.5_num*xi0x(ix)*xi1y(iy) &
                  + 1.0_num/3.0_num * xi1x(ix) * xi1y(iy)

              ! This is the bit that actually solves d(rho)/dt = -div(J)
              jxh(ix, iy) = &
                  jxh(ix-1, iy) - part_q * wx * 1.0_num/dt_j * part_weight/dy
              jyh(ix, iy) = &
                  jyh(ix, iy-1) - part_q * wy * 1.0_num/dt_j * part_weight/dx
              jzh(ix, iy) = part_q * part_vz * wz * part_weight/(dx*dy)

              jx(cell_x1+ix, cell_y1+iy) = &
                  jx(cell_x1+ix, cell_y1+iy) +jxh(ix, iy)
              jy(cell_x1+ix, cell_y1+iy) = &
                  jy(cell_x1+ix, cell_y1+iy) +jyh(ix, iy)
              jz(cell_x1+ix, cell_y1+iy) = &
                  jz(cell_x1+ix, cell_y1+iy) +jzh(ix, iy)
            ENDDO
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

              ! probe lines are defined by two points. particles crossing the
              ! line are recorded assumes transit from left -> right. Put the
              ! top at the bottom and bottom at the top to catch
              ! right-> left moving particles.
              probe_x1 = current_probe%vertex_bottom(1)
              probe_y1 = current_probe%vertex_bottom(2)
              probe_x2 = current_probe%vertex_top(1)
              probe_y2 = current_probe%vertex_top(2)

              probe_a = (probe_y1 * probe_x2 - probe_x1 * probe_y2) &
                  / (probe_y1 - probe_y2)
              probe_b = (probe_x1 - probe_x2) / (probe_y1 - probe_y2)

              ! direction?
              IF (probe_y2 .LT. probe_y1) THEN
                probe_temp = init_part_x
                init_part_x = final_part_x
                final_part_x = probe_temp
                probe_temp = init_part_y
                init_part_y = final_part_y
                final_part_y = probe_temp
              ENDIF

              IF ((final_part_y .LT. MAX(probe_y1, probe_y2)) &
                  .AND. (final_part_y .GT. MIN(probe_y1, probe_y2))) THEN
                IF (init_part_x .LE. probe_a + probe_b * init_part_y) THEN
                  IF (final_part_x .GT. probe_a + probe_b * final_part_y) THEN

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
      CALL processor_summation_bcs(ekbar_sum(:,:,ispecies))
      CALL field_bc(ekbar_sum(:,:,ispecies))
      CALL processor_summation_bcs(ct(:,:,ispecies))
      CALL field_bc(ct(:,:,ispecies))
    ENDDO

    ! Calculate the mean kinetic energy for each species in space
    ekbar = 0.0_num
    DO ispecies = 1, n_species
      DO iy = 1, ny
        DO ix = 1, nx
          mean = ekbar_sum(ix, iy, ispecies) &
              / MAX(ct(ix, iy, ispecies), c_non_zero)
          ekbar(ix, iy, ispecies) = mean
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(xi0x)
    DEALLOCATE(xi1x)
    DEALLOCATE(xi0y)
    DEALLOCATE(xi1y)
    DEALLOCATE(jxh)
    DEALLOCATE(jyh)
    DEALLOCATE(jzh)

    CALL particle_bcs

    IF (smooth_currents) CALL smooth_current()

  END SUBROUTINE push_particles

END MODULE particles
