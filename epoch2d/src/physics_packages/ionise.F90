MODULE ionise

  USE calc_df
  USE split_particle

  IMPLICIT NONE

CONTAINS

  SUBROUTINE ionise_particles

#ifdef PARTICLE_IONISE
    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next, new_part
    TYPE(particle_list), ALLOCATABLE :: append_lists(:)
    REAL(num) :: part_x, part_x2, cell_x_r, cell_frac_x
    REAL(num) :: part_y, part_y2, cell_y_r, cell_frac_y
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, hx, gy, hy
    REAL(num) :: ex_part, ey_part, ez_part, e_part
    REAL(num) :: number_at_part, delta_en, erat
    INTEGER :: cell_x1, cell_x2, ix
    INTEGER :: cell_y1, cell_y2, iy
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: number_density
    REAL(num) :: rand, mom_frac, ion_rate, ion_prob
    ! Specific to Keldysh equation
    REAL(num) :: e_norm, ion_en, prefactor
    REAL(num) :: ion_fac, fac, cf2
    REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
    INTEGER :: next_species
    INTEGER, PARAMETER :: dcellx = 0, dcelly = 0

    rand = random()
#ifdef PARTICLE_SHAPE_BSPLINE3
    fac = 1.0_num / 24.0_num
#elif  PARTICLE_SHAPE_TOPHAT
    fac = 1.0_num
#else
    fac = 0.5_num
#endif

    ! This prefactor is specific to the Keldysh model of ionisation
    prefactor = SQRT(6.0_num * pi) / (2.0_num)**(5.0_num/4.0_num) / 2.42e-17_num

    ALLOCATE(number_density(-2:nx+3,-2:ny+3))

    ALLOCATE(append_lists(n_species))
    DO ispecies = 1, n_species
      CALL create_empty_partlist(append_lists(ispecies))
    ENDDO

    DO ispecies = 1, n_species
      IF (.NOT. species_list(ispecies)%ionise) CYCLE
      mom_frac = &
          species_list(species_list(ispecies)%ionise_to_species)%mass
      mom_frac = mom_frac / (mom_frac + species_list(ispecies)%mass)
      CALL calc_number_density(number_density, ispecies)

      current=>species_list(ispecies)%attached_list%head
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
        number_at_part = 0.0_num
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            ex_part = ex_part + hx(ix) * gy(iy) * ex(cell_x2+ix, cell_y1+iy)
            ey_part = ey_part + gx(ix) * hy(iy) * ey(cell_x1+ix, cell_y2+iy)
            ez_part = ez_part + gx(ix) * gy(iy) * ez(cell_x1+ix, cell_y1+iy)
            number_at_part = number_at_part + dx * dy &
                * gx(ix) * gy(iy) * number_density(cell_x1+ix, cell_y1+iy)
          ENDDO
        ENDDO
        e_part = fac * SQRT(ex_part**2 + ey_part**2 + ez_part**2)
        number_at_part = fac**2 * number_at_part

        ! In this section need to calculate the ionisation rate
        ! And store it in the variable 'ion_rate'.
        ! The rest of the routine is general

        ! The code is supplied with a
        ! simple implementation of Keldysh equation
        ! Keldysh (Zh. Eksp. Teor. Fiz. 47, 1945 (1964))
        ! In English (Sov. Phys. JETP 20, 1307 (1965))

        e_norm = e_part / 5.1401e11_num
        ion_en = species_list(ispecies)%ionisation_energy &
            / (27.211_num * ev)
        ion_fac = e_norm / (2.0_num * ion_en)**1.5_num
        ion_rate = prefactor * ion_en * SQRT(ion_fac) &
            * EXP(-two_thirds / ion_fac)

        ! End of ionisation rate calculation

        ! Convert an ionisation rate to a probability
        ion_prob = dt * ion_rate * current%weight / number_at_part
        rand = random()
        ! After all that, we now know the target ionisation probability
        ! So convert
        IF (rand .LT. ion_prob) THEN
          !delta_en = current%weight &
          !    * species_list(ispecies)%ionisation_energy
          !erat = SQRT(1.0_num - 2.0_num * delta_en &
          !    / (e_part**2 * epsilon0 * dx * dy))
          !DO iy = -1, 1
          !  DO ix = -1, 1
          !    ex(cell_x2+ix, cell_y1+iy) = &
          !        erat * hx(ix) * gy(iy) * ex(cell_x2+ix, cell_y1+iy)
          !    ey(cell_x1+ix, cell_y2+iy) = &
          !        erat * gx(ix) * hy(iy) * ey(cell_x1+ix, cell_y2+iy)
          !    ez(cell_x1+ix, cell_y1+iy) = &
          !        erat * gx(ix) * gy(iy) * ez(cell_x1+ix, cell_y1+iy)
          !  ENDDO
          !ENDDO
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, current)
          next_species = species_list(ispecies)%ionise_to_species
          CALL add_particle_to_partlist(&
              species_list(next_species)%attached_list, current)
          next_species = species_list(ispecies)%release_species
          IF (next_species .GT. 0) THEN
            ALLOCATE(new_part)
            CALL init_particle(new_part)
#ifdef PER_PARTICLE_WEIGHT
            new_part%weight = current%weight
#endif
            new_part%part_pos = current%part_pos
            new_part%part_p = current%part_p * mom_frac
            current%part_p = current%part_p * (1.0_num - mom_frac)
#ifdef PER_PARTICLE_CHARGE_MASS
            new_part%charge = species_list(next_species)%charge
            new_part%mass = species_list(next_species)%mass
#endif
#ifdef PARTICLE_DEBUG
            new_part%processor = rank
            new_part%processor_at_t0 = rank
#endif
            CALL add_particle_to_partlist(append_lists(next_species), new_part)
          ENDIF
        ENDIF

        current=>next
      ENDDO
    ENDDO

    DO ispecies = 1, n_species
      CALL append_partlist(species_list(ispecies)%attached_list, &
          append_lists(ispecies))
    ENDDO

    DEALLOCATE(number_density, append_lists)
#endif

  END SUBROUTINE ionise_particles

END MODULE ionise
