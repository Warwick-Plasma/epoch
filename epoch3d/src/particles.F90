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
    INTEGER :: cell_y1, cell_y2, cell_y3
    INTEGER :: cell_z1, cell_z2, cell_z3

    ! Xi (space factor see page 38 in manual)
    ! The code now uses gx and hx instead of xi0 and xi1

    INTEGER, PARAMETER :: sf = sf_order

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    REAL(num), DIMENSION(-sf-2:sf+1,-sf-1:sf+1,-sf-1:sf+1) :: jxh
    REAL(num), DIMENSION(-sf-1:sf+1,-sf-2:sf+1,-sf-1:sf+1) :: jyh
    REAL(num), DIMENSION(-sf-1:sf+1,-sf-1:sf+1,-sf-2:sf+1) :: jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: part_px, part_py, part_pz
    REAL(num) :: part_vx, part_vy, part_vz
    REAL(num) :: part_q, part_mc, part_weight

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifdef PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    REAL(num) :: init_part_z, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy
#endif

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(-sf-1:sf+1) :: gx, gy, gz

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(-sf-1:sf+1) :: hx, hy, hz

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: pxp, pxm, pyp, pym, pzp, pzm
    REAL(num) :: tau, taux, tauy, tauz

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio

    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: idx, idy, idz
    REAL(num) :: idtxy, idtyz, idtxz
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: fcx, fcy, fcz, fjx, fjy, fjz
    REAL(num) :: root, mean, fac, dtfac, third
    INTEGER :: ispecies, dcellx, dcelly, dcellz
    INTEGER(KIND=8) :: ipart

    TYPE(particle), POINTER :: current, next

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    jxh = 0.0_num
    jyh = 0.0_num
    jzh = 0.0_num

    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz
    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    third = 1.0_num / 3.0_num
    ! particle weighting multiplication factor
#ifdef SPLINE_FOUR
    fac = 1.0_num / 24.0_num
#else
    fac = 0.5_num
#endif
    dtfac = 0.5_num * dt * fac**3

    idtyz = idt * idy * idz * fac**3
    idtxz = idt * idx * idz * fac**3
    idtxy = idt * idx * idy * fac**3

    part_weight = weight
    fcx = idtyz * part_weight
    fcy = idtxz * part_weight
    fcz = idtxy * part_weight

    DO ispecies = 1, n_species
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGEMASS
      part_q  = particle_species(ispecies)%charge
      part_mc = c * particle_species(ispecies)%mass
#endif
      ! -- this option needs more testing -- DEC$ IVDEP
      !DEC$ VECTOR ALWAYS
      !DEC$ NOPREFETCH current
      DO ipart = 1, particle_species(ispecies)%attached_list%count
        next=>current%next
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
        fcx = idtyz * part_weight
        fcy = idtxz * part_weight
        fcz = idtxy * part_weight
#endif
#ifdef PARTICLE_PROBES
        init_part_x = current%part_pos(1)
        init_part_y = current%part_pos(2)
        init_part_z = current%part_pos(3)
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        ! Use a lookup table for charge and mass to SAVE memory
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
        part_y = part_y + part_py * root
        part_z = part_z + part_pz * root

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
        cell_x_r = part_x / dx
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
        cell_y_r = part_y / dy
        ! Round cell position to nearest cell
        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
        cell_z_r = part_z / dz
        ! Round cell position to nearest cell
        cell_z1 = FLOOR(cell_z_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_z = REAL(cell_z1, num) - cell_z_r
        cell_z1 = cell_z1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
#ifdef SPLINE_FOUR
        gx(-2) = (1.5_num - cell_frac_x)**4
        gx(-1) = 4.0_num * (1.1875_num + cell_frac_x * (cell_frac_x &
            * (1.5_num + cell_frac_x - cell_frac_x**2) - 2.75_num))
        gx( 0) = 6.0_num * (115/48 + cell_frac_x**2 &
            * (cell_frac_x**2 - 2.5_num))
        gx( 1) = 4.0_num * (1.1875_num + cell_frac_x * (cell_frac_x &
            * (1.5_num - cell_frac_x - cell_frac_x**2) + 2.75_num))
        gx( 2) = (1.5_num + cell_frac_x)**4

        gy(-2) = (1.5_num - cell_frac_y)**4
        gy(-1) = 4.0_num * (1.1875_num + cell_frac_y * (cell_frac_y &
            * (1.5_num + cell_frac_y - cell_frac_y**2) - 2.75_num))
        gy( 0) = 6.0_num * (115/48 + cell_frac_y**2 &
            * (cell_frac_y**2 - 2.5_num))
        gy( 1) = 4.0_num * (1.1875_num + cell_frac_y * (cell_frac_y &
            * (1.5_num - cell_frac_y - cell_frac_y**2) + 2.75_num))
        gy( 2) = (1.5_num + cell_frac_y)**4

        gz(-2) = (1.5_num - cell_frac_z)**4
        gz(-1) = 4.0_num * (1.1875_num + cell_frac_z * (cell_frac_z &
            * (1.5_num + cell_frac_z - cell_frac_z**2) - 2.75_num))
        gz( 0) = 6.0_num * (115/48 + cell_frac_z**2 &
            * (cell_frac_z**2 - 2.5_num))
        gz( 1) = 4.0_num * (1.1875_num + cell_frac_z * (cell_frac_z &
            * (1.5_num - cell_frac_z - cell_frac_z**2) + 2.75_num))
        gz( 2) = (1.5_num + cell_frac_z)**4
#else
        gx(-1) = (0.5_num + cell_frac_x)**2
        gx( 0) =  1.5_num - 2.0_num * cell_frac_x**2
        gx( 1) = (0.5_num - cell_frac_x)**2

        gy(-1) = (0.5_num + cell_frac_y)**2
        gy( 0) =  1.5_num - 2.0_num * cell_frac_y**2
        gy( 1) = (0.5_num - cell_frac_y)**2

        gz(-1) = (0.5_num + cell_frac_z)**2
        gz( 0) =  1.5_num - 2.0_num * cell_frac_z**2
        gz( 1) = (0.5_num - cell_frac_z)**2
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

        cell_z2 = FLOOR(cell_z_r)
        cell_frac_z = REAL(cell_z2, num) - cell_z_r + 0.5_num
        cell_z2 = cell_z2 + 1

#ifdef SPLINE_FOUR
        hx(-2) = (1.5_num - cell_frac_x)**4
        hx(-1) = 4.0_num * (1.1875_num + cell_frac_x * (cell_frac_x &
            * (1.5_num + cell_frac_x - cell_frac_x**2) - 2.75_num))
        hx( 0) = 6.0_num * (115/48 + cell_frac_x**2 &
            * (cell_frac_x**2 - 2.5_num))
        hx( 1) = 4.0_num * (1.1875_num + cell_frac_x * (cell_frac_x &
            * (1.5_num - cell_frac_x - cell_frac_x**2) + 2.75_num))
        hx( 2) = (1.5_num + cell_frac_x)**4

        hy(-2) = (1.5_num - cell_frac_y)**4
        hy(-1) = 4.0_num * (1.1875_num + cell_frac_y * (cell_frac_y &
            * (1.5_num + cell_frac_y - cell_frac_y**2) - 2.75_num))
        hy( 0) = 6.0_num * (115/48 + cell_frac_y**2 &
            * (cell_frac_y**2 - 2.5_num))
        hy( 1) = 4.0_num * (1.1875_num + cell_frac_y * (cell_frac_y &
            * (1.5_num - cell_frac_y - cell_frac_y**2) + 2.75_num))
        hy( 2) = (1.5_num + cell_frac_y)**4

        hz(-2) = (1.5_num - cell_frac_z)**4
        hz(-1) = 4.0_num * (1.1875_num + cell_frac_z * (cell_frac_z &
            * (1.5_num + cell_frac_z - cell_frac_z**2) - 2.75_num))
        hz( 0) = 6.0_num * (115/48 + cell_frac_z**2 &
            * (cell_frac_z**2 - 2.5_num))
        hz( 1) = 4.0_num * (1.1875_num + cell_frac_z * (cell_frac_z &
            * (1.5_num - cell_frac_z - cell_frac_z**2) + 2.75_num))
        hz( 2) = (1.5_num + cell_frac_z)**4
#else
        hx(-1) = (0.5_num + cell_frac_x)**2
        hx( 0) =  1.5_num - 2.0_num * cell_frac_x**2
        hx( 1) = (0.5_num - cell_frac_x)**2

        hy(-1) = (0.5_num + cell_frac_y)**2
        hy( 0) =  1.5_num - 2.0_num * cell_frac_y**2
        hy( 1) = (0.5_num - cell_frac_y)**2

        hz(-1) = (0.5_num + cell_frac_z)**2
        hz( 0) =  1.5_num - 2.0_num * cell_frac_z**2
        hz( 1) = (0.5_num - cell_frac_z)**2
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
#ifdef SPLINE_FOUR
        ex_part = &
              gz(-2) * (gy(-2) * (hx(-2)*ex(cell_x2-2, cell_y1-2, cell_z1-2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-2, cell_z1-2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-2, cell_z1-2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-2, cell_z1-2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-2, cell_z1-2)) &
            +           gy(-1) * (hx(-2)*ex(cell_x2-2, cell_y1-1, cell_z1-2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1-2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1-2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1-2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-1, cell_z1-2)) &
            +           gy( 0) * (hx(-2)*ex(cell_x2-2, cell_y1  , cell_z1-2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1-2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1-2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1-2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1  , cell_z1-2)) &
            +           gy( 1) * (hx(-2)*ex(cell_x2-2, cell_y1+1, cell_z1-2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1-2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1-2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1-2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+1, cell_z1-2)) &
            +           gy( 2) * (hx(-2)*ex(cell_x2-2, cell_y1+2, cell_z1-2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+2, cell_z1-2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+2, cell_z1-2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+2, cell_z1-2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+2, cell_z1-2))) &
            + gz(-1) * (gy(-2) * (hx(-2)*ex(cell_x2-2, cell_y1-2, cell_z1-1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-2, cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-2, cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-2, cell_z1-1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-2, cell_z1-1)) &
            +           gy(-1) * (hx(-2)*ex(cell_x2-2, cell_y1-1, cell_z1-1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1-1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-1, cell_z1-1)) &
            +           gy( 0) * (hx(-2)*ex(cell_x2-2, cell_y1  , cell_z1-1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1-1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1  , cell_z1-1)) &
            +           gy( 1) * (hx(-2)*ex(cell_x2-2, cell_y1+1, cell_z1-1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1-1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+1, cell_z1-1)) &
            +           gy( 2) * (hx(-2)*ex(cell_x2-2, cell_y1+2, cell_z1-1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+2, cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+2, cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+2, cell_z1-1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+2, cell_z1-1))) &
            + gz( 0) * (gy(-2) * (hx(-2)*ex(cell_x2-2, cell_y1-2, cell_z1  ) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-2, cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-2, cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-2, cell_z1  ) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-2, cell_z1  )) &
            +           gy(-1) * (hx(-2)*ex(cell_x2-2, cell_y1-1, cell_z1  ) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1  ) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-1, cell_z1  )) &
            +           gy( 0) * (hx(-2)*ex(cell_x2-2, cell_y1  , cell_z1  ) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1  ) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1  , cell_z1  )) &
            +           gy( 1) * (hx(-2)*ex(cell_x2-2, cell_y1+1, cell_z1  ) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1  ) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+1, cell_z1  )) &
            +           gy( 2) * (hx(-2)*ex(cell_x2-2, cell_y1+2, cell_z1  ) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+2, cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+2, cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+2, cell_z1  ) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+2, cell_z1  ))) &
            + gz( 1) * (gy(-2) * (hx(-2)*ex(cell_x2-2, cell_y1-2, cell_z1+1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-2, cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-2, cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-2, cell_z1+1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-2, cell_z1+1)) &
            +           gy(-1) * (hx(-2)*ex(cell_x2-2, cell_y1-1, cell_z1+1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1+1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-1, cell_z1+1)) &
            +           gy( 0) * (hx(-2)*ex(cell_x2-2, cell_y1  , cell_z1+1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1+1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1  , cell_z1+1)) &
            +           gy( 1) * (hx(-2)*ex(cell_x2-2, cell_y1+1, cell_z1+1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1+1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+1, cell_z1+1)) &
            +           gy( 2) * (hx(-2)*ex(cell_x2-2, cell_y1+2, cell_z1+1) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+2, cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+2, cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+2, cell_z1+1) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+2, cell_z1+1))) &
            + gz( 2) * (gy(-2) * (hx(-2)*ex(cell_x2-2, cell_y1-2, cell_z1+2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-2, cell_z1+2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-2, cell_z1+2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-2, cell_z1+2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-2, cell_z1+2)) &
            +           gy(-1) * (hx(-2)*ex(cell_x2-2, cell_y1-1, cell_z1+2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1+2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1+2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1+2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1-1, cell_z1+2)) &
            +           gy( 0) * (hx(-2)*ex(cell_x2-2, cell_y1  , cell_z1+2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1+2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1+2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1+2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1  , cell_z1+2)) &
            +           gy( 1) * (hx(-2)*ex(cell_x2-2, cell_y1+1, cell_z1+2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1+2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1+2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1+2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+1, cell_z1+2)) &
            +           gy( 2) * (hx(-2)*ex(cell_x2-2, cell_y1+2, cell_z1+2) &
            +                     hx(-1)*ex(cell_x2-1, cell_y1+2, cell_z1+2) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+2, cell_z1+2) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+2, cell_z1+2) &
            +                     hx( 2)*ex(cell_x2+2, cell_y1+2, cell_z1+2)))

        ey_part = &
              gz(-2) * (hy(-2) * (gx(-2)*ey(cell_x1-2, cell_y2-2, cell_z1-2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-2, cell_z1-2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-2, cell_z1-2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-2, cell_z1-2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-2, cell_z1-2)) &
            +           hy(-1) * (gx(-2)*ey(cell_x1-2, cell_y2-1, cell_z1-2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1-2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1-2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1-2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-1, cell_z1-2)) &
            +           hy( 0) * (gx(-2)*ey(cell_x1-2, cell_y2  , cell_z1-2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1-2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1-2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1-2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2  , cell_z1-2)) &
            +           hy( 1) * (gx(-2)*ey(cell_x1-2, cell_y2+1, cell_z1-2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1-2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1-2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1-2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+1, cell_z1-2)) &
            +           hy( 2) * (gx(-2)*ey(cell_x1-2, cell_y2+2, cell_z1-2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+2, cell_z1-2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+2, cell_z1-2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+2, cell_z1-2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+2, cell_z1-2))) &
            + gz(-1) * (hy(-2) * (gx(-2)*ey(cell_x1-2, cell_y2-2, cell_z1-1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-2, cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-2, cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-2, cell_z1-1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-2, cell_z1-1)) &
            +           hy(-1) * (gx(-2)*ey(cell_x1-2, cell_y2-1, cell_z1-1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1-1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-1, cell_z1-1)) &
            +           hy( 0) * (gx(-2)*ey(cell_x1-2, cell_y2  , cell_z1-1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1-1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2  , cell_z1-1)) &
            +           hy( 1) * (gx(-2)*ey(cell_x1-2, cell_y2+1, cell_z1-1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1-1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+1, cell_z1-1)) &
            +           hy( 2) * (gx(-2)*ey(cell_x1-2, cell_y2+2, cell_z1-1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+2, cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+2, cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+2, cell_z1-1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+2, cell_z1-1))) &
            + gz( 0) * (hy(-2) * (gx(-2)*ey(cell_x1-2, cell_y2-2, cell_z1  ) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-2, cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-2, cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-2, cell_z1  ) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-2, cell_z1  )) &
            +           hy(-1) * (gx(-2)*ey(cell_x1-2, cell_y2-1, cell_z1  ) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1  ) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-1, cell_z1  )) &
            +           hy( 0) * (gx(-2)*ey(cell_x1-2, cell_y2  , cell_z1  ) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1  ) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2  , cell_z1  )) &
            +           hy( 1) * (gx(-2)*ey(cell_x1-2, cell_y2+1, cell_z1  ) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1  ) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+1, cell_z1  )) &
            +           hy( 2) * (gx(-2)*ey(cell_x1-2, cell_y2+2, cell_z1  ) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+2, cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+2, cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+2, cell_z1  ) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+2, cell_z1  ))) &
            + gz( 1) * (hy(-2) * (gx(-2)*ey(cell_x1-2, cell_y2-2, cell_z1+1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-2, cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-2, cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-2, cell_z1+1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-2, cell_z1+1)) &
            +           hy(-1) * (gx(-2)*ey(cell_x1-2, cell_y2-1, cell_z1+1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1+1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-1, cell_z1+1)) &
            +           hy( 0) * (gx(-2)*ey(cell_x1-2, cell_y2  , cell_z1+1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1+1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2  , cell_z1+1)) &
            +           hy( 1) * (gx(-2)*ey(cell_x1-2, cell_y2+1, cell_z1+1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1+1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+1, cell_z1+1)) &
            +           hy( 2) * (gx(-2)*ey(cell_x1-2, cell_y2+2, cell_z1+1) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+2, cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+2, cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+2, cell_z1+1) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+2, cell_z1+1))) &
            + gz( 2) * (hy(-2) * (gx(-2)*ey(cell_x1-2, cell_y2-2, cell_z1+2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-2, cell_z1+2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-2, cell_z1+2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-2, cell_z1+2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-2, cell_z1+2)) &
            +           hy(-1) * (gx(-2)*ey(cell_x1-2, cell_y2-1, cell_z1+2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1+2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1+2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1+2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2-1, cell_z1+2)) &
            +           hy( 0) * (gx(-2)*ey(cell_x1-2, cell_y2  , cell_z1+2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1+2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1+2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1+2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2  , cell_z1+2)) &
            +           hy( 1) * (gx(-2)*ey(cell_x1-2, cell_y2+1, cell_z1+2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1+2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1+2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1+2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+1, cell_z1+2)) &
            +           hy( 2) * (gx(-2)*ey(cell_x1-2, cell_y2+2, cell_z1+2) &
            +                     gx(-1)*ey(cell_x1-1, cell_y2+2, cell_z1+2) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+2, cell_z1+2) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+2, cell_z1+2) &
            +                     gx( 2)*ey(cell_x1+2, cell_y2+2, cell_z1+2)))

        ez_part = &
              hz(-2) * (gy(-2) * (gx(-2)*ez(cell_x1-2, cell_y1-2, cell_z2-2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-2, cell_z2-2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-2, cell_z2-2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-2, cell_z2-2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-2, cell_z2-2)) &
            +           gy(-1) * (gx(-2)*ez(cell_x1-2, cell_y1-1, cell_z2-2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2-2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2-2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2-2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-1, cell_z2-2)) &
            +           gy( 0) * (gx(-2)*ez(cell_x1-2, cell_y1  , cell_z2-2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2-2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2-2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2-2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1  , cell_z2-2)) &
            +           gy( 1) * (gx(-2)*ez(cell_x1-2, cell_y1+1, cell_z2-2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2-2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2-2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2-2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+1, cell_z2-2)) &
            +           gy( 2) * (gx(-2)*ez(cell_x1-2, cell_y1+2, cell_z2-2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+2, cell_z2-2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+2, cell_z2-2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+2, cell_z2-2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+2, cell_z2-2))) &
            + hz(-1) * (gy(-2) * (gx(-2)*ez(cell_x1-2, cell_y1-2, cell_z2-1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-2, cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-2, cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-2, cell_z2-1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-2, cell_z2-1)) &
            +           gy(-1) * (gx(-2)*ez(cell_x1-2, cell_y1-1, cell_z2-1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2-1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-1, cell_z2-1)) &
            +           gy( 0) * (gx(-2)*ez(cell_x1-2, cell_y1  , cell_z2-1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2-1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1  , cell_z2-1)) &
            +           gy( 1) * (gx(-2)*ez(cell_x1-2, cell_y1+1, cell_z2-1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2-1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+1, cell_z2-1)) &
            +           gy( 2) * (gx(-2)*ez(cell_x1-2, cell_y1+2, cell_z2-1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+2, cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+2, cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+2, cell_z2-1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+2, cell_z2-1))) &
            + hz( 0) * (gy(-2) * (gx(-2)*ez(cell_x1-2, cell_y1-2, cell_z2  ) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-2, cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-2, cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-2, cell_z2  ) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-2, cell_z2  )) &
            +           gy(-1) * (gx(-2)*ez(cell_x1-2, cell_y1-1, cell_z2  ) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2  ) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-1, cell_z2  )) &
            +           gy( 0) * (gx(-2)*ez(cell_x1-2, cell_y1  , cell_z2  ) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2  ) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1  , cell_z2  )) &
            +           gy( 1) * (gx(-2)*ez(cell_x1-2, cell_y1+1, cell_z2  ) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2  ) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+1, cell_z2  )) &
            +           gy( 2) * (gx(-2)*ez(cell_x1-2, cell_y1+2, cell_z2  ) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+2, cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+2, cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+2, cell_z2  ) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+2, cell_z2  ))) &
            + hz( 1) * (gy(-2) * (gx(-2)*ez(cell_x1-2, cell_y1-2, cell_z2+1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-2, cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-2, cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-2, cell_z2+1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-2, cell_z2+1)) &
            +           gy(-1) * (gx(-2)*ez(cell_x1-2, cell_y1-1, cell_z2+1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2+1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-1, cell_z2+1)) &
            +           gy( 0) * (gx(-2)*ez(cell_x1-2, cell_y1  , cell_z2+1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2+1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1  , cell_z2+1)) &
            +           gy( 1) * (gx(-2)*ez(cell_x1-2, cell_y1+1, cell_z2+1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2+1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+1, cell_z2+1)) &
            +           gy( 2) * (gx(-2)*ez(cell_x1-2, cell_y1+2, cell_z2+1) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+2, cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+2, cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+2, cell_z2+1) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+2, cell_z2+1))) &
            + hz( 2) * (gy(-2) * (gx(-2)*ez(cell_x1-2, cell_y1-2, cell_z2+2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-2, cell_z2+2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-2, cell_z2+2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-2, cell_z2+2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-2, cell_z2+2)) &
            +           gy(-1) * (gx(-2)*ez(cell_x1-2, cell_y1-1, cell_z2+2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2+2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2+2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2+2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1-1, cell_z2+2)) &
            +           gy( 0) * (gx(-2)*ez(cell_x1-2, cell_y1  , cell_z2+2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2+2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2+2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2+2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1  , cell_z2+2)) &
            +           gy( 1) * (gx(-2)*ez(cell_x1-2, cell_y1+1, cell_z2+2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2+2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2+2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2+2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+1, cell_z2+2)) &
            +           gy( 2) * (gx(-2)*ez(cell_x1-2, cell_y1+2, cell_z2+2) &
            +                     gx(-1)*ez(cell_x1-1, cell_y1+2, cell_z2+2) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+2, cell_z2+2) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+2, cell_z2+2) &
            +                     gx( 2)*ez(cell_x1+2, cell_y1+2, cell_z2+2)))

        bx_part = &
              hz(-2) * (hy(-2) * (gx(-2)*bx(cell_x1-2, cell_y2-2, cell_z2-2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-2, cell_z2-2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-2, cell_z2-2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-2, cell_z2-2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-2, cell_z2-2)) &
            +           hy(-1) * (gx(-2)*bx(cell_x1-2, cell_y2-1, cell_z2-2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2-2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2-2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2-2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-1, cell_z2-2)) &
            +           hy( 0) * (gx(-2)*bx(cell_x1-2, cell_y2  , cell_z2-2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2-2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2-2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2-2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2  , cell_z2-2)) &
            +           hy( 1) * (gx(-2)*bx(cell_x1-2, cell_y2+1, cell_z2-2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2-2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2-2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2-2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+1, cell_z2-2)) &
            +           hy( 2) * (gx(-2)*bx(cell_x1-2, cell_y2+2, cell_z2-2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+2, cell_z2-2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+2, cell_z2-2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+2, cell_z2-2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+2, cell_z2-2))) &
            + hz(-1) * (hy(-2) * (gx(-2)*bx(cell_x1-2, cell_y2-2, cell_z2-1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-2, cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-2, cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-2, cell_z2-1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-2, cell_z2-1)) &
            +           hy(-1) * (gx(-2)*bx(cell_x1-2, cell_y2-1, cell_z2-1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2-1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-1, cell_z2-1)) &
            +           hy( 0) * (gx(-2)*bx(cell_x1-2, cell_y2  , cell_z2-1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2-1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2  , cell_z2-1)) &
            +           hy( 1) * (gx(-2)*bx(cell_x1-2, cell_y2+1, cell_z2-1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2-1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+1, cell_z2-1)) &
            +           hy( 2) * (gx(-2)*bx(cell_x1-2, cell_y2+2, cell_z2-1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+2, cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+2, cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+2, cell_z2-1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+2, cell_z2-1))) &
            + hz( 0) * (hy(-2) * (gx(-2)*bx(cell_x1-2, cell_y2-2, cell_z2  ) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-2, cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-2, cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-2, cell_z2  ) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-2, cell_z2  )) &
            +           hy(-1) * (gx(-2)*bx(cell_x1-2, cell_y2-1, cell_z2  ) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2  ) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-1, cell_z2  )) &
            +           hy( 0) * (gx(-2)*bx(cell_x1-2, cell_y2  , cell_z2  ) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2  ) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2  , cell_z2  )) &
            +           hy( 1) * (gx(-2)*bx(cell_x1-2, cell_y2+1, cell_z2  ) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2  ) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+1, cell_z2  )) &
            +           hy( 2) * (gx(-2)*bx(cell_x1-2, cell_y2+2, cell_z2  ) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+2, cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+2, cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+2, cell_z2  ) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+2, cell_z2  ))) &
            + hz( 1) * (hy(-2) * (gx(-2)*bx(cell_x1-2, cell_y2-2, cell_z2+1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-2, cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-2, cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-2, cell_z2+1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-2, cell_z2+1)) &
            +           hy(-1) * (gx(-2)*bx(cell_x1-2, cell_y2-1, cell_z2+1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2+1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-1, cell_z2+1)) &
            +           hy( 0) * (gx(-2)*bx(cell_x1-2, cell_y2  , cell_z2+1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2+1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2  , cell_z2+1)) &
            +           hy( 1) * (gx(-2)*bx(cell_x1-2, cell_y2+1, cell_z2+1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2+1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+1, cell_z2+1)) &
            +           hy( 2) * (gx(-2)*bx(cell_x1-2, cell_y2+2, cell_z2+1) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+2, cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+2, cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+2, cell_z2+1) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+2, cell_z2+1))) &
            + hz( 2) * (hy(-2) * (gx(-2)*bx(cell_x1-2, cell_y2-2, cell_z2+2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-2, cell_z2+2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-2, cell_z2+2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-2, cell_z2+2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-2, cell_z2+2)) &
            +           hy(-1) * (gx(-2)*bx(cell_x1-2, cell_y2-1, cell_z2+2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2+2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2+2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2+2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2-1, cell_z2+2)) &
            +           hy( 0) * (gx(-2)*bx(cell_x1-2, cell_y2  , cell_z2+2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2+2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2+2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2+2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2  , cell_z2+2)) &
            +           hy( 1) * (gx(-2)*bx(cell_x1-2, cell_y2+1, cell_z2+2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2+2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2+2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2+2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+1, cell_z2+2)) &
            +           hy( 2) * (gx(-2)*bx(cell_x1-2, cell_y2+2, cell_z2+2) &
            +                     gx(-1)*bx(cell_x1-1, cell_y2+2, cell_z2+2) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+2, cell_z2+2) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+2, cell_z2+2) &
            +                     gx( 2)*bx(cell_x1+2, cell_y2+2, cell_z2+2)))

        by_part = &
              hz(-2) * (gy(-2) * (hx(-2)*by(cell_x2-2, cell_y1-2, cell_z2-2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-2, cell_z2-2) &
            +                     hx( 0)*by(cell_x2  , cell_y1-2, cell_z2-2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-2, cell_z2-2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-2, cell_z2-2)) &
            +           gy(-1) * (hx(-2)*by(cell_x2-2, cell_y1-1, cell_z2-2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2-2) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2-2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2-2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-1, cell_z2-2)) &
            +           gy( 0) * (hx(-2)*by(cell_x2-2, cell_y1  , cell_z2-2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1  , cell_z2-2) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2-2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2-2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1  , cell_z2-2)) &
            +           gy( 1) * (hx(-2)*by(cell_x2-2, cell_y1+1, cell_z2-2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2-2) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2-2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2-2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+1, cell_z2-2)) &
            +           gy( 2) * (hx(-2)*by(cell_x2-2, cell_y1+2, cell_z2-2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+2, cell_z2-2) &
            +                     hx( 0)*by(cell_x2  , cell_y1+2, cell_z2-2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+2, cell_z2-2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+2, cell_z2-2))) &
            + hz(-1) * (gy(-2) * (hx(-2)*by(cell_x2-2, cell_y1-2, cell_z2-1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-2, cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1-2, cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-2, cell_z2-1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-2, cell_z2-1)) &
            +           gy(-1) * (hx(-2)*by(cell_x2-2, cell_y1-1, cell_z2-1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2-1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-1, cell_z2-1)) &
            +           gy( 0) * (hx(-2)*by(cell_x2-2, cell_y1  , cell_z2-1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1  , cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2-1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1  , cell_z2-1)) &
            +           gy( 1) * (hx(-2)*by(cell_x2-2, cell_y1+1, cell_z2-1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2-1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+1, cell_z2-1)) &
            +           gy( 2) * (hx(-2)*by(cell_x2-2, cell_y1+2, cell_z2-1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+2, cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1+2, cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+2, cell_z2-1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+2, cell_z2-1))) &
            + hz( 0) * (gy(-2) * (hx(-2)*by(cell_x2-2, cell_y1-2, cell_z2  ) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-2, cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1-2, cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-2, cell_z2  ) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-2, cell_z2  )) &
            +           gy(-1) * (hx(-2)*by(cell_x2-2, cell_y1-1, cell_z2  ) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2  ) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-1, cell_z2  )) &
            +           gy( 0) * (hx(-2)*by(cell_x2-2, cell_y1  , cell_z2  ) &
            +                     hx(-1)*by(cell_x2-1, cell_y1  , cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2  ) &
            +                     hx( 2)*by(cell_x2+2, cell_y1  , cell_z2  )) &
            +           gy( 1) * (hx(-2)*by(cell_x2-2, cell_y1+1, cell_z2  ) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2  ) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+1, cell_z2  )) &
            +           gy( 2) * (hx(-2)*by(cell_x2-2, cell_y1+2, cell_z2  ) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+2, cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1+2, cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+2, cell_z2  ) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+2, cell_z2  ))) &
            + hz( 1) * (gy(-2) * (hx(-2)*by(cell_x2-2, cell_y1-2, cell_z2+1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-2, cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1-2, cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-2, cell_z2+1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-2, cell_z2+1)) &
            +           gy(-1) * (hx(-2)*by(cell_x2-2, cell_y1-1, cell_z2+1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2+1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-1, cell_z2+1)) &
            +           gy( 0) * (hx(-2)*by(cell_x2-2, cell_y1  , cell_z2+1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1  , cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2+1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1  , cell_z2+1)) &
            +           gy( 1) * (hx(-2)*by(cell_x2-2, cell_y1+1, cell_z2+1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2+1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+1, cell_z2+1)) &
            +           gy( 2) * (hx(-2)*by(cell_x2-2, cell_y1+2, cell_z2+1) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+2, cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1+2, cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+2, cell_z2+1) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+2, cell_z2+1))) &
            + hz( 2) * (gy(-2) * (hx(-2)*by(cell_x2-2, cell_y1-2, cell_z2+2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-2, cell_z2+2) &
            +                     hx( 0)*by(cell_x2  , cell_y1-2, cell_z2+2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-2, cell_z2+2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-2, cell_z2+2)) &
            +           gy(-1) * (hx(-2)*by(cell_x2-2, cell_y1-1, cell_z2+2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2+2) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2+2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2+2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1-1, cell_z2+2)) &
            +           gy( 0) * (hx(-2)*by(cell_x2-2, cell_y1  , cell_z2+2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1  , cell_z2+2) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2+2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2+2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1  , cell_z2+2)) &
            +           gy( 1) * (hx(-2)*by(cell_x2-2, cell_y1+1, cell_z2+2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2+2) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2+2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2+2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+1, cell_z2+2)) &
            +           gy( 2) * (hx(-2)*by(cell_x2-2, cell_y1+2, cell_z2+2) &
            +                     hx(-1)*by(cell_x2-1, cell_y1+2, cell_z2+2) &
            +                     hx( 0)*by(cell_x2  , cell_y1+2, cell_z2+2) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+2, cell_z2+2) &
            +                     hx( 2)*by(cell_x2+2, cell_y1+2, cell_z2+2)))

        bz_part = &
              gz(-2) * (hy(-2) * (hx(-2)*bz(cell_x2-2, cell_y2-2, cell_z1-2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-2, cell_z1-2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-2, cell_z1-2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-2, cell_z1-2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-2, cell_z1-2)) &
            +           hy(-1) * (hx(-2)*bz(cell_x2-2, cell_y2-1, cell_z1-2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1-2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1-2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1-2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-1, cell_z1-2)) &
            +           hy( 0) * (hx(-2)*bz(cell_x2-2, cell_y2  , cell_z1-2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1-2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1-2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1-2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2  , cell_z1-2)) &
            +           hy( 1) * (hx(-2)*bz(cell_x2-2, cell_y2+1, cell_z1-2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1-2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1-2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1-2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+1, cell_z1-2)) &
            +           hy( 2) * (hx(-2)*bz(cell_x2-2, cell_y2+2, cell_z1-2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+2, cell_z1-2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+2, cell_z1-2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+2, cell_z1-2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+2, cell_z1-2))) &
            + gz(-1) * (hy(-2) * (hx(-2)*bz(cell_x2-2, cell_y2-2, cell_z1-1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-2, cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-2, cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-2, cell_z1-1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-2, cell_z1-1)) &
            +           hy(-1) * (hx(-2)*bz(cell_x2-2, cell_y2-1, cell_z1-1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1-1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-1, cell_z1-1)) &
            +           hy( 0) * (hx(-2)*bz(cell_x2-2, cell_y2  , cell_z1-1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1-1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2  , cell_z1-1)) &
            +           hy( 1) * (hx(-2)*bz(cell_x2-2, cell_y2+1, cell_z1-1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1-1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+1, cell_z1-1)) &
            +           hy( 2) * (hx(-2)*bz(cell_x2-2, cell_y2+2, cell_z1-1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+2, cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+2, cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+2, cell_z1-1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+2, cell_z1-1))) &
            + gz( 0) * (hy(-2) * (hx(-2)*bz(cell_x2-2, cell_y2-2, cell_z1  ) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-2, cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-2, cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-2, cell_z1  ) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-2, cell_z1  )) &
            +           hy(-1) * (hx(-2)*bz(cell_x2-2, cell_y2-1, cell_z1  ) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1  ) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-1, cell_z1  )) &
            +           hy( 0) * (hx(-2)*bz(cell_x2-2, cell_y2  , cell_z1  ) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1  ) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2  , cell_z1  )) &
            +           hy( 1) * (hx(-2)*bz(cell_x2-2, cell_y2+1, cell_z1  ) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1  ) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+1, cell_z1  )) &
            +           hy( 2) * (hx(-2)*bz(cell_x2-2, cell_y2+2, cell_z1  ) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+2, cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+2, cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+2, cell_z1  ) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+2, cell_z1  ))) &
            + gz( 1) * (hy(-2) * (hx(-2)*bz(cell_x2-2, cell_y2-2, cell_z1+1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-2, cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-2, cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-2, cell_z1+1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-2, cell_z1+1)) &
            +           hy(-1) * (hx(-2)*bz(cell_x2-2, cell_y2-1, cell_z1+1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1+1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-1, cell_z1+1)) &
            +           hy( 0) * (hx(-2)*bz(cell_x2-2, cell_y2  , cell_z1+1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1+1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2  , cell_z1+1)) &
            +           hy( 1) * (hx(-2)*bz(cell_x2-2, cell_y2+1, cell_z1+1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1+1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+1, cell_z1+1)) &
            +           hy( 2) * (hx(-2)*bz(cell_x2-2, cell_y2+2, cell_z1+1) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+2, cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+2, cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+2, cell_z1+1) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+2, cell_z1+1))) &
            + gz( 2) * (hy(-2) * (hx(-2)*bz(cell_x2-2, cell_y2-2, cell_z1+2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-2, cell_z1+2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-2, cell_z1+2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-2, cell_z1+2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-2, cell_z1+2)) &
            +           hy(-1) * (hx(-2)*bz(cell_x2-2, cell_y2-1, cell_z1+2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1+2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1+2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1+2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2-1, cell_z1+2)) &
            +           hy( 0) * (hx(-2)*bz(cell_x2-2, cell_y2  , cell_z1+2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1+2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1+2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1+2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2  , cell_z1+2)) &
            +           hy( 1) * (hx(-2)*bz(cell_x2-2, cell_y2+1, cell_z1+2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1+2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1+2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1+2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+1, cell_z1+2)) &
            +           hy( 2) * (hx(-2)*bz(cell_x2-2, cell_y2+2, cell_z1+2) &
            +                     hx(-1)*bz(cell_x2-1, cell_y2+2, cell_z1+2) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+2, cell_z1+2) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+2, cell_z1+2) &
            +                     hx( 2)*bz(cell_x2+2, cell_y2+2, cell_z1+2)))
#else
        ex_part = &
              gz(-1) * (gy(-1) * (hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1-1)) &
            +           gy( 0) * (hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1-1)) &
            +           gy( 1) * (hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1-1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1-1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1-1))) &
            + gz( 0) * (gy(-1) * (hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1  )) &
            +           gy( 0) * (hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1  )) &
            +           gy( 1) * (hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1  ) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1  ) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1  ))) &
            + gz( 1) * (gy(-1) * (hx(-1)*ex(cell_x2-1, cell_y1-1, cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1-1, cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1-1, cell_z1+1)) &
            +           gy( 0) * (hx(-1)*ex(cell_x2-1, cell_y1  , cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1  , cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1  , cell_z1+1)) &
            +           gy( 1) * (hx(-1)*ex(cell_x2-1, cell_y1+1, cell_z1+1) &
            +                     hx( 0)*ex(cell_x2  , cell_y1+1, cell_z1+1) &
            +                     hx( 1)*ex(cell_x2+1, cell_y1+1, cell_z1+1)))

        ey_part = &
              gz(-1) * (hy(-1) * (gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1-1)) &
            +           hy( 0) * (gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1-1)) &
            +           hy( 1) * (gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1-1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1-1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1-1))) &
            + gz( 0) * (hy(-1) * (gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1  )) &
            +           hy( 0) * (gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1  )) &
            +           hy( 1) * (gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1  ) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1  ) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1  ))) &
            + gz( 1) * (hy(-1) * (gx(-1)*ey(cell_x1-1, cell_y2-1, cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2-1, cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2-1, cell_z1+1)) &
            +           hy( 0) * (gx(-1)*ey(cell_x1-1, cell_y2  , cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2  , cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2  , cell_z1+1)) &
            +           hy( 1) * (gx(-1)*ey(cell_x1-1, cell_y2+1, cell_z1+1) &
            +                     gx( 0)*ey(cell_x1  , cell_y2+1, cell_z1+1) &
            +                     gx( 1)*ey(cell_x1+1, cell_y2+1, cell_z1+1)))

        ez_part = &
              hz(-1) * (gy(-1) * (gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2-1)) &
            +           gy( 0) * (gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2-1)) &
            +           gy( 1) * (gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2-1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2-1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2-1))) &
            + hz( 0) * (gy(-1) * (gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2  )) &
            +           gy( 0) * (gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2  )) &
            +           gy( 1) * (gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2  ) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2  ) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2  ))) &
            + hz( 1) * (gy(-1) * (gx(-1)*ez(cell_x1-1, cell_y1-1, cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1-1, cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1-1, cell_z2+1)) &
            +           gy( 0) * (gx(-1)*ez(cell_x1-1, cell_y1  , cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1  , cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1  , cell_z2+1)) &
            +           gy( 1) * (gx(-1)*ez(cell_x1-1, cell_y1+1, cell_z2+1) &
            +                     gx( 0)*ez(cell_x1  , cell_y1+1, cell_z2+1) &
            +                     gx( 1)*ez(cell_x1+1, cell_y1+1, cell_z2+1)))

        bx_part = &
              hz(-1) * (hy(-1) * (gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2-1)) &
            +           hy( 0) * (gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2-1)) &
            +           hy( 1) * (gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2-1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2-1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2-1))) &
            + hz( 0) * (hy(-1) * (gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2  )) &
            +           hy( 0) * (gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2  )) &
            +           hy( 1) * (gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2  ) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2  ) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2  ))) &
            + hz( 1) * (hy(-1) * (gx(-1)*bx(cell_x1-1, cell_y2-1, cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2-1, cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2-1, cell_z2+1)) &
            +           hy( 0) * (gx(-1)*bx(cell_x1-1, cell_y2  , cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2  , cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2  , cell_z2+1)) &
            +           hy( 1) * (gx(-1)*bx(cell_x1-1, cell_y2+1, cell_z2+1) &
            +                     gx( 0)*bx(cell_x1  , cell_y2+1, cell_z2+1) &
            +                     gx( 1)*bx(cell_x1+1, cell_y2+1, cell_z2+1)))

        by_part = &
              hz(-1) * (gy(-1) * (hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2-1)) &
            +           gy( 0) * (hx(-1)*by(cell_x2-1, cell_y1  , cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2-1)) &
            +           gy( 1) * (hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2-1) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2-1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2-1))) &
            + hz( 0) * (gy(-1) * (hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2  )) &
            +           gy( 0) * (hx(-1)*by(cell_x2-1, cell_y1  , cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2  )) &
            +           gy( 1) * (hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2  ) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2  ) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2  ))) &
            + hz( 1) * (gy(-1) * (hx(-1)*by(cell_x2-1, cell_y1-1, cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1-1, cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1-1, cell_z2+1)) &
            +           gy( 0) * (hx(-1)*by(cell_x2-1, cell_y1  , cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1  , cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1  , cell_z2+1)) &
            +           gy( 1) * (hx(-1)*by(cell_x2-1, cell_y1+1, cell_z2+1) &
            +                     hx( 0)*by(cell_x2  , cell_y1+1, cell_z2+1) &
            +                     hx( 1)*by(cell_x2+1, cell_y1+1, cell_z2+1)))

        bz_part = &
              gz(-1) * (hy(-1) * (hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1-1)) &
            +           hy( 0) * (hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1-1)) &
            +           hy( 1) * (hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1-1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1-1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1-1))) &
            + gz( 0) * (hy(-1) * (hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1  )) &
            +           hy( 0) * (hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1  )) &
            +           hy( 1) * (hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1  ) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1  ) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1  ))) &
            + gz( 1) * (hy(-1) * (hx(-1)*bz(cell_x2-1, cell_y2-1, cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2-1, cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2-1, cell_z1+1)) &
            +           hy( 0) * (hx(-1)*bz(cell_x2-1, cell_y2  , cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2  , cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2  , cell_z1+1)) &
            +           hy( 1) * (hx(-1)*bz(cell_x2-1, cell_y2+1, cell_z1+1) &
            +                     hx( 0)*bz(cell_x2  , cell_y2+1, cell_z1+1) &
            +                     hx( 1)*bz(cell_x2+1, cell_y2+1, cell_z1+1)))
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
        root = c / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)

        part_vx = part_px * root
        part_vy = part_py * root
        part_vz = part_pz * root

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + part_vx * dto2
        part_y = part_y + part_vy * dto2
        part_z = part_z + part_vz * dto2

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = (/ part_x + x_min_local, &
            part_y + y_min_local, part_z + z_min_local /)
        current%part_p   = (/ part_px, part_py, part_pz /)

#ifdef PARTICLE_PROBES
        final_part_x = current%part_pos(1)
        final_part_y = current%part_pos(2)
        final_part_z = current%part_pos(3)
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
          part_x = part_x + part_vx * dto2
          part_y = part_y + part_vy * dto2
          part_z = part_z + part_vz * dto2

          cell_x_r = part_x / dx
          cell_x3 = FLOOR(cell_x_r + 0.5_num)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r
          cell_x3 = cell_x3 + 1

          cell_y_r = part_y / dy
          cell_y3 = FLOOR(cell_y_r + 0.5_num)
          cell_frac_y = REAL(cell_y3, num) - cell_y_r
          cell_y3 = cell_y3 + 1

          cell_z_r = part_z / dz
          cell_z3 = FLOOR(cell_z_r + 0.5_num)
          cell_frac_z = REAL(cell_z3, num) - cell_z_r
          cell_z3 = cell_z3 + 1

          hx = 0.0_num
          hy = 0.0_num
          hz = 0.0_num

          dcellx = cell_x3 - cell_x1
          dcelly = cell_y3 - cell_y1
          dcellz = cell_z3 - cell_z1

#ifdef SPLINE_FOUR
          hx(dcellx-2) = (1.5_num - cell_frac_x)**4
          hx(dcellx-1) = 4.0_num * (1.1875_num + cell_frac_x * (cell_frac_x &
              * (1.5_num + cell_frac_x - cell_frac_x**2) - 2.75_num))
          hx(dcellx  ) = 6.0_num * (115/48 + cell_frac_x**2 &
              * (cell_frac_x**2 - 2.5_num))
          hx(dcellx+1) = 4.0_num * (1.1875_num + cell_frac_x * (cell_frac_x &
              * (1.5_num - cell_frac_x - cell_frac_x**2) + 2.75_num))
          hx(dcellx+2) = (1.5_num + cell_frac_x)**4

          hy(dcelly-2) = (1.5_num - cell_frac_y)**4
          hy(dcelly-1) = 4.0_num * (1.1875_num + cell_frac_y * (cell_frac_y &
              * (1.5_num + cell_frac_y - cell_frac_y**2) - 2.75_num))
          hy(dcelly  ) = 6.0_num * (115/48 + cell_frac_y**2 &
              * (cell_frac_y**2 - 2.5_num))
          hy(dcelly+1) = 4.0_num * (1.1875_num + cell_frac_y * (cell_frac_y &
              * (1.5_num - cell_frac_y - cell_frac_y**2) + 2.75_num))
          hy(dcelly+2) = (1.5_num + cell_frac_y)**4

          hz(dcellz-2) = (1.5_num - cell_frac_z)**4
          hz(dcellz-1) = 4.0_num * (1.1875_num + cell_frac_z * (cell_frac_z &
              * (1.5_num + cell_frac_z - cell_frac_z**2) - 2.75_num))
          hz(dcellz  ) = 6.0_num * (115/48 + cell_frac_z**2 &
              * (cell_frac_z**2 - 2.5_num))
          hz(dcellz+1) = 4.0_num * (1.1875_num + cell_frac_z * (cell_frac_z &
              * (1.5_num - cell_frac_z - cell_frac_z**2) + 2.75_num))
          hz(dcellz+2) = (1.5_num + cell_frac_z)**4
#else
          hx(dcellx-1) = (0.5_num + cell_frac_x)**2
          hx(dcellx  ) =  1.5_num - 2.0_num * ABS(cell_frac_x)**2
          hx(dcellx+1) = (0.5_num - cell_frac_x)**2

          hy(dcelly-1) = (0.5_num + cell_frac_y)**2
          hy(dcelly  ) =  1.5_num - 2.0_num * ABS(cell_frac_y)**2
          hy(dcelly+1) = (0.5_num - cell_frac_y)**2

          hz(dcellz-1) = (0.5_num + cell_frac_z)**2
          hz(dcellz  ) =  1.5_num - 2.0_num * ABS(cell_frac_z)**2
          hz(dcellz+1) = (0.5_num - cell_frac_z)**2
#endif

          ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
          ! the current update much simpler
          hx = hx - gx
          hy = hy - gy
          hz = hz - gz

          ! Remember that due to CFL condition particle can never cross more
          ! than one gridcell in one timestep

          xmin = -sf_order + (dcellx - 1) / 2
          xmax =  sf_order + (dcellx + 1) / 2

          ymin = -sf_order + (dcelly - 1) / 2
          ymax =  sf_order + (dcelly + 1) / 2

          zmin = -sf_order + (dcellz - 1) / 2
          zmax =  sf_order + (dcellz + 1) / 2

          ! Set these to zero due to diffential inside loop
          jxh = 0.0_num
          jyh = 0.0_num
          jzh = 0.0_num

          fjx = fcx * part_q
          fjy = fcy * part_q
          fjz = fcz * part_q

          DO iz = zmin, zmax
            DO iy = ymin, ymax
              DO ix = xmin, xmax
                wx =  hx(ix) * (gy(iy) * (gz(iz) + 0.5_num * hz(iz)) &
                    + hy(iy) * (third  *  hz(iz) + 0.5_num * gz(iz)))
                wy =  hy(iy) * (gx(ix) * (gz(iz) + 0.5_num * hz(iz)) &
                    + hx(ix) * (third  *  hz(iz) + 0.5_num * gz(iz)))
                wz =  hz(iz) * (gx(ix) * (gy(iy) + 0.5_num * hy(iy)) &
                    + hx(ix) * (third  *  hy(iy) + 0.5_num * gy(iy)))

                ! This is the bit that actually solves d(rho)/dt = -div(J)
                jxh(ix, iy, iz) = jxh(ix-1, iy, iz) - fjx * wx
                jyh(ix, iy, iz) = jyh(ix, iy-1, iz) - fjy * wy
                jzh(ix, iy, iz) = jzh(ix, iy, iz-1) - fjz * wz

                jx(cell_x1+ix, cell_y1+iy, cell_z1+iz) = &
                    jx(cell_x1+ix, cell_y1+iy, cell_z1+iz) + jxh(ix, iy, iz)
                jy(cell_x1+ix, cell_y1+iy, cell_z1+iz) = &
                    jy(cell_x1+ix, cell_y1+iy, cell_z1+iz) + jyh(ix, iy, iz)
                jz(cell_x1+ix, cell_y1+iy, cell_z1+iz) = &
                    jz(cell_x1+ix, cell_y1+iy, cell_z1+iz) + jzh(ix, iy, iz)

              ENDDO
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
              c * (SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2) &
              - part_mc)

          ! right energy? (in J)
          IF (probe_energy .GT. current_probe%ek_min) THEN
            IF ((probe_energy .LT. current_probe%ek_max) &
                .OR. (current_probe%ek_max .LT. 0.0_num)) THEN

              d_init = SUM(current_probe%normal &
                  * (current_probe%corner(1,:) &
                  - (/init_part_x, init_part_y, init_part_z/)))
              d_final = SUM(current_probe%normal &
                  * (current_probe%corner(1,:) &
                  - (/final_part_x, final_part_y, final_part_z/)))
              IF (SIGN(1.0_num, d_init)*SIGN(1.0_num, d_final) &
                  .LE. 0.0_num) THEN
                ! this particle is wanted so copy it to the list associated
                ! with this probe
                ALLOCATE(particle_copy)
                particle_copy = current
                CALL add_particle_to_partlist(current_probe%sampled_particles, &
                    particle_copy)
                NULLIFY(particle_copy)
              ENDIF

            ENDIF
          ENDIF
          current_probe=>current_probe%next
        ENDDO
#endif
        current=>next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(jx)
    CALL processor_summation_bcs(jy)
    CALL processor_summation_bcs(jz)

    CALL particle_bcs

  END SUBROUTINE push_particles

END MODULE particles
