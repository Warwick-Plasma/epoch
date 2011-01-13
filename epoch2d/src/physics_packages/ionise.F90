MODULE ionise

  USE calc_df
  USE helper

  IMPLICIT NONE

CONTAINS

  SUBROUTINE ionise_particles

#ifdef PARTICLE_IONISE
    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next, new_part
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x
    REAL(num) :: part_y, part_y2, cell_y_r, cell_frac_y
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, hx, gy, hy
    REAL(num) :: ex_part, ey_part, ez_part, e_part2
    REAL(num) :: number_density_part, ndp_low, ndp_high
    INTEGER :: cell_x1, cell_x2, ix
    INTEGER :: cell_y1, cell_y2, iy
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: number_density, nd_low, nd_high
    REAL(num) :: lambda_db, e_photon, t_eff, saha_rhs, ion_frac, rand
    REAL(num) :: fac, tfac, lfac, cf2
    INTEGER :: idum, next_species
    INTEGER, PARAMETER :: dcellx = 0, dcelly = 0

    idum = -1445
    rand = random(idum)
#ifdef PARTICLE_SHAPE_BSPLINE3
    fac = 1.0_num / 24.0_num
#elif  PARTICLE_SHAPE_TOPHAT
    fac = 1.0_num
#else
    fac = 0.5_num
#endif
    tfac = epsilon0 * fac**2 / kb / 3.0_num
    lfac = h_planck / SQRT(2.0_num * pi * m0 * kb)

    ALLOCATE(number_density(-2:nx+3,-2:ny+3))
    ALLOCATE(nd_low (-2:nx+3,-2:ny+3))
    ALLOCATE(nd_high(-2:nx+3,-2:ny+3))

    DO ispecies = 1, n_species
      IF (.NOT. particle_species(ispecies)%ionise) CYCLE
      CALL calc_number_density(nd_low, ispecies)
      CALL calc_number_density(nd_high, &
          particle_species(ispecies)%ionise_to_species)
      number_density = nd_low+nd_high
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next=>current%next
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local

        ! Work out number of grid cells in the particle is
        ! Not in general an integer
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        ! Round cell position to nearest cell
        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! These are now the weighting factors correct for field weighting
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

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        ! Grid weighting factors in 3D (3D analogue of equation 4.77 page 25
        ! of manual)
        ! These weight grid properties onto particles
#ifdef PARTICLE_SHAPE_BSPLINE3
        INCLUDE '../include/bspline3/hx_dcell.inc'
#elif  PARTICLE_SHAPE_TOPHAT
        INCLUDE '../include/tophat/hx_dcell.inc'
#else
        INCLUDE '../include/triangle/hx_dcell.inc'
#endif

        ex_part = 0.0_num
        ey_part = 0.0_num
        ez_part = 0.0_num
        number_density_part = 0.0_num
        ndp_low = 0.0_num
        ndp_high = 0.0_num
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            ex_part = ex_part + hx(ix) * gy(iy) * ex(cell_x2+ix, cell_y1+iy)
            ey_part = ey_part + gx(ix) * hy(iy) * ey(cell_x1+ix, cell_y2+iy)
            ez_part = ez_part + gx(ix) * gy(iy) * ez(cell_x1+ix, cell_y1+iy)
            number_density_part = number_density_part &
                + gx(ix) * gy(iy) * number_density(cell_x1+ix, cell_y1+iy)
            ndp_low  = ndp_low &
                + gx(ix) * gy(iy) * nd_low(cell_x1+ix, cell_y1+iy)
            ndp_high = ndp_high &
                + gx(ix) * gy(iy) * nd_high(cell_x1+ix, cell_y1+iy)
          ENDDO
        ENDDO
        e_part2 = ex_part**2 + ey_part**2 + ez_part**2
        number_density_part = fac * number_density_part

        ! This is a first attempt at using the 1 level Saha equation to
        ! calculate an ionisation fraction. This isn't really a very good model!

        ! e_photon = 0.5_num * epsilon0 * fac**2 * e_part2 * dx * dy
        ! t_eff = 2.0_num / 3.0_num * e_photon &
        !    / (kb * number_density_part * dx * dy)
        t_eff = tfac * e_part2 / number_density_part
        IF (t_eff .GT. 1.0e-6_num) THEN
          lambda_db = lfac / SQRT(t_eff)
          saha_rhs = &
              EXP(-particle_species(ispecies)%ionisation_energy / kb / t_eff) &
              / lambda_db**3
          ion_frac = -saha_rhs &
              + SQRT(saha_rhs**2 + 2.0_num * number_density_part * saha_rhs)
          ion_frac = ion_frac / number_density_part
          IF (ion_frac .GT. 1.0_num) ion_frac = 1.0_num
        ELSE
          ion_frac = 0.0_num
        ENDIF

        rand = random(idum)
        ! After all that, we now know the target ionisation fraction, so
        ! subtract the current fraction and ionise
        IF (rand .LT. (ion_frac-ndp_high/MAX(ndp_low, c_non_zero))) THEN
          CALL remove_particle_from_partlist(&
              particle_species(ispecies)%attached_list, current)
          next_species = particle_species(ispecies)%ionise_to_species
          CALL add_particle_to_partlist(&
              particle_species(next_species)%attached_list, current)
          next_species = particle_species(ispecies)%release_species
          IF (next_species .GT. 0) THEN
            ALLOCATE(new_part)
#ifdef PER_PARTICLE_WEIGHT
            new_part%weight = current%weight
#endif
            new_part%part_pos = current%part_pos
            new_part%part_p = current%part_p
#ifdef PER_PARTICLE_CHARGE_MASS
            new_part%charge = particle_species(next_species)%charge
            new_part%mass = particle_species(next_species)%mass
#endif
#ifdef PARTICLE_DEBUG
            new_part%processor = rank
            new_part%processor_at_t0 = rank
#endif
            CALL add_particle_to_partlist(&
                particle_species(next_species)%attached_list, new_part)
          ENDIF
        ENDIF

        current=>next
      ENDDO

    ENDDO

    DEALLOCATE(number_density)
    DEALLOCATE(nd_low, nd_high)
#endif

  END SUBROUTINE ionise_particles

END MODULE ionise
