MODULE ionise

  USE shared_data
  USE partlist
  USE calc_df
  USE helper
  IMPLICIT NONE

CONTAINS

  SUBROUTINE ionise_particles

#ifdef PART_IONISE
    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next, new_part
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x
    REAL(num) :: part_y, part_y2, cell_y_r, cell_frac_y
    REAL(num), DIMENSION(-1:1) :: hx, hy, gx, gy
    REAL(num) :: ex_part, ey_part, ez_part, e_part
    REAL(num) :: number_density_part, ndp_low, ndp_high
    INTEGER :: cell_x1, cell_y1, cell_x2, cell_y2, ix, iy, next_species
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: number_density, nd_low, nd_high
    REAL(num) :: lambda_db, e_photon, t_eff, saha_rhs, ion_frac, rand
    INTEGER :: idum

    idum = -1445
    rand = random(idum)

    ALLOCATE(number_density(-2:nx+3, -2:ny+3))
    ALLOCATE(nd_low(-2:nx+3, -2:ny+3), nd_high(-2:nx+3, -2:ny+3))
    DO ispecies = 1, n_species
      IF (.NOT. particle_species(ispecies)%ionise) CYCLE
      CALL calc_number_density(nd_low, ispecies)
      CALL calc_number_density(nd_high, particle_species(ispecies)%ionise_to_species)
      number_density = nd_low+nd_high
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next=>current%next
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local

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

        ! These are now the weighting factors correct for field weighting
        gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
        gx( 0) = 0.75_num - cell_frac_x**2
        gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

        gy(-1) = 0.5_num * (0.5_num + cell_frac_y)**2
        gy( 0) = 0.75_num - cell_frac_y**2
        gy( 1) = 0.5_num * (0.5_num - cell_frac_y)**2

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

        ! Grid weighting factors in 3D (3D analogue of equation 4.77 page 25 of manual)
        ! These weight grid properties onto particles
        hx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
        hx( 0) = 0.75_num - cell_frac_x**2
        hx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

        hy(-1) = 0.5_num * (0.5_num + cell_frac_y)**2
        hy( 0) = 0.75_num - cell_frac_y**2
        hy( 1) = 0.5_num * (0.5_num - cell_frac_y)**2

        ex_part = 0.0_num
        ey_part = 0.0_num
        ez_part = 0.0_num
        number_density_part = 0.0_num
        ndp_low = 0.0_num
        ndp_high = 0.0_num
        DO ix = -1, 1
          DO iy = -1, 1
            ex_part = ex_part+hx(ix)*gy(iy)*ex(cell_x2+ix, cell_y1+iy)
            ey_part = ey_part+gx(ix)*hy(iy)*ex(cell_x1+ix, cell_y2+iy)
            ez_part = ez_part+gx(ix)*hy(iy)*ex(cell_x1+ix, cell_y1+iy)
            number_density_part = gx(ix)*gy(iy)* number_density(cell_x1+ix, cell_y1+iy)
            ndp_low = gx(ix)*gy(iy)*nd_low(cell_x1+ix, cell_y1+iy)
            ndp_high = gx(ix)*gy(iy)*nd_high(cell_x1+ix, cell_y1+iy)

          ENDDO
        ENDDO
        e_part = SQRT(ex_part**2+ey_part**2+ez_part**2)

        ! This is a first attempt at using the 1 level Saha equation to calculate an
        ! Ionisation fraction. This isn't really a very good model!

        e_photon = 0.5_num * epsilon0 * (e_part)**2 * dx * dy
        t_eff = 2.0_num/3.0_num*e_photon/(kb*number_density_part*dx*dy)
        IF (t_eff .GT. 1.0e-6_num) THEN
          lambda_db = SQRT(h_planck**2/(2.0_num*pi*m0*kb*t_eff))
          saha_rhs = 2.0_num/lambda_db**3 * EXP(-particle_species(ispecies)%ionisation_energy/(kb*t_eff))
          ion_frac = 0.5_num * (-saha_rhs + SQRT(saha_rhs**2+4.0_num*number_density_part*saha_rhs))
          ion_frac = ion_frac/number_density_part
        ELSE
          ion_frac = 0.0_num
        ENDIF
        IF (ion_frac .GT. 1.0_num) ion_frac = 1.0_num

        rand = random(idum)
        ! After all that, we now know the target ionisation fraction, so subtract the current fraction and ionise
        IF (rand .LT. (ion_frac-ndp_high/MAX(ndp_low, c_non_zero))) THEN
          CALL remove_particle_from_partlist(particle_species(ispecies)%attached_list, current)
          next_species = particle_species(ispecies)%ionise_to_species
          CALL add_particle_to_partlist(particle_species(next_species)%attached_list, current)
          next_species = particle_species(ispecies)%release_species
          IF (next_species .GT. 0) THEN
            ALLOCATE(new_part)
#ifdef PER_PARTICLE_WEIGHT
            new_part%weight = current%weight
#endif
            new_part%part_pos = current%part_pos
            new_part%part_p = current%part_p
#ifdef PER_PARTICLE_CHARGEMASS
            new_part%charge = particle_species(next_species)%charge
            new_part%mass = particle_species(next_species)%mass
#endif
#ifdef PART_DEBUG
            new_part%processor = rank
            new_part%processor_at_t0 = rank
#endif
            CALL add_particle_to_partlist(particle_species(next_species)%attached_list, new_part)
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
