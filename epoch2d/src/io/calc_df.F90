MODULE calc_df

  USE boundary
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_mass_density(data_array, current_species)

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    REAL(num) :: fac, idx
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_m  = species_list(ispecies)%mass
#endif
#ifndef PER_PARTICLE_WEIGHT
      fac = species_list(ispecies)%weight * idx
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m  = current%mass
#endif
#ifdef PER_PARTICLE_WEIGHT
        fac = current%weight * idx
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        wdata = part_m * fac
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_ekbar(data_array, current_species)

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_px, part_py, part_pz, part_mc
    ! The weight of a particle
    REAL(num) :: l_weight, l_weightc
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: data_array
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: ct
    INTEGER, INTENT(IN) :: current_species
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    ALLOCATE(ct(-2:nx+3,-2:ny+3))
    data_array = 0.0_num
    ct = 0.0_num

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc = c * species_list(ispecies)%mass
#endif
#ifndef PER_PARTICLE_WEIGHT
      l_weight = species_list(ispecies)%weight
      l_weightc = c * l_weight
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
        l_weightc = c * l_weight
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        wdata = (SQRT(part_px**2 + part_py**2 + part_pz**2 + part_mc**2) &
          - part_mc) * l_weightc
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
            ct(cell_x+ix, cell_y+iy) = &
                ct(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * l_weight
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    CALL processor_summation_bcs(ct)

    data_array = data_array / MAX(ct, c_non_zero)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(ct)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_charge_density(data_array, current_species)

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    REAL(num) :: fac, idx
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    idx = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q  = species_list(ispecies)%charge
#endif
#ifndef PER_PARTICLE_WEIGHT
      fac = species_list(ispecies)%weight * idx
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
#endif
#ifdef PER_PARTICLE_WEIGHT
        fac = current%weight * idx
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        wdata = part_q * fac
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_charge_density



  SUBROUTINE calc_number_density(data_array, current_species)

    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    REAL(num) :: idx
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    data_array = 0.0_num

    idx   = 1.0_num / dx / dy

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_WEIGHT
      wdata = species_list(ispecies)%weight * idx
#endif
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_WEIGHT
        wdata = current%weight * idx
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_number_density



  SUBROUTINE calc_temperature(sigma, current_species)

    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: l_weight
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    REAL(num) :: gf
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    ALLOCATE(meanx(-2:nx+3,-2:ny+3))
    ALLOCATE(meany(-2:nx+3,-2:ny+3))
    ALLOCATE(meanz(-2:nx+3,-2:ny+3))
    ALLOCATE(part_count(-2:nx+3,-2:ny+3))
    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    DO ispecies = spec_start, spec_end
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      sqrt_part_m  = SQRT(species_list(ispecies)%mass)
#endif
#ifndef PER_PARTICLE_WEIGHT
      l_weight = species_list(ispecies)%weight
#endif
      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            gf = gx(ix) * gy(iy) * l_weight
            meanx(cell_x+ix, cell_y+iy) = &
                meanx(cell_x+ix, cell_y+iy) + gf * part_pmx
            meany(cell_x+ix, cell_y+iy) = &
                meany(cell_x+ix, cell_y+iy) + gf * part_pmy
            meanz(cell_x+ix, cell_y+iy) = &
                meanz(cell_x+ix, cell_y+iy) + gf * part_pmz
            part_count(cell_x+ix, cell_y+iy) = &
                part_count(cell_x+ix, cell_y+iy) + gf
          ENDDO
        ENDDO
        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(meanx)
    CALL processor_summation_bcs(meany)
    CALL processor_summation_bcs(meanz)
    CALL processor_summation_bcs(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      sqrt_part_m  = SQRT(species_list(ispecies)%mass)
#endif
      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            gf = gx(ix) * gy(iy)
            sigma(cell_x+ix, cell_y+iy) = &
                sigma(cell_x+ix, cell_y+iy) + gf &
                * ((part_pmx - meanx(cell_x+ix, cell_y+iy))**2 &
                + (part_pmy - meany(cell_x+ix, cell_y+iy))**2 &
                + (part_pmz - meanz(cell_x+ix, cell_y+iy))**2)
            part_count(cell_x+ix, cell_y+iy) = &
                part_count(cell_x+ix, cell_y+iy) + gf
          ENDDO
        ENDDO
        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(sigma)
    CALL processor_summation_bcs(part_count)

    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / 2.0_num

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num), DIMENSION(-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    INTERFACE
      FUNCTION evaluator(a_particle, species_eval)
        USE shared_data
        TYPE(particle), POINTER :: a_particle
        INTEGER, INTENT(IN) :: species_eval
        REAL(num) :: evaluator
      END FUNCTION evaluator
    END INTERFACE

    data_array = 0.0_num

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>species_list(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx
        cell_y_r = (current%part_pos(2) - y_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        wdata = evaluator(current, ispecies)
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_on_grid_with_evaluator

END MODULE calc_df
