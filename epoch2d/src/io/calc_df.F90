MODULE calc_df

  USE boundary
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_mass_density(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
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
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
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

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: l_weight
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    ALLOCATE(wt(-2:nx+3,-2:ny+3))
    data_array = 0.0_num
    wt = 0.0_num

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc = c * species_list(ispecies)%mass
#ifndef PER_PARTICLE_WEIGHT
      l_weight = species_list(ispecies)%weight
#endif
      fac = part_mc * l_weight * c
#else
#ifndef PER_PARTICLE_WEIGHT
      l_weight = species_list(ispecies)%weight
      fac = part_mc * l_weight * c
#endif
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif
        fac = part_mc * l_weight * c
#else
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
        fac = part_mc * l_weight * c
#endif
#endif
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc

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

        gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        wdata = (gamma - 1.0_num) * fac
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
            wt(cell_x+ix, cell_y+iy) = &
                wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * l_weight
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    CALL processor_summation_bcs(wt)

    data_array = data_array / MAX(wt, c_non_zero)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_ekflux(data_array, current_species, direction)

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: l_weight
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma, ek, part_flux, xfac, yfac, zfac
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    ALLOCATE(wt(-2:nx+3,-2:ny+3))
    data_array = 0.0_num
    wt = 0.0_num

    xfac = c * dy
    yfac = c * dx
    zfac = c * dx * dy

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
      current=>species_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc = c * species_list(ispecies)%mass
#ifndef PER_PARTICLE_WEIGHT
      l_weight = species_list(ispecies)%weight
#endif
      fac = part_mc * l_weight * c
#else
#ifndef PER_PARTICLE_WEIGHT
      l_weight = species_list(ispecies)%weight
      fac = part_mc * l_weight * c
#endif
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_mc = c * current%mass
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif
        fac = part_mc * l_weight * c
#else
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
        fac = part_mc * l_weight * c
#endif
#endif
        part_ux = current%part_p(1) / part_mc
        part_uy = current%part_p(2) / part_mc
        part_uz = current%part_p(3) / part_mc

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

        gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        ek = (gamma - 1.0_num) * fac

        SELECT CASE(direction)
        CASE(-c_dir_x)
          ! negative flux in x
          part_flux = xfac * part_ux / gamma
          wdata = -ek * MIN(part_flux, 0.0_num)
        CASE( c_dir_x)
          ! positive flux in x
          part_flux = xfac * part_ux / gamma
          wdata =  ek * MAX(part_flux, 0.0_num)
        CASE(-c_dir_y)
          ! negative flux in y
          part_flux = yfac * part_uy / gamma
          wdata = -ek * MIN(part_flux, 0.0_num)
        CASE( c_dir_y)
          ! positive flux in y
          part_flux = yfac * part_uy / gamma
          wdata =  ek * MAX(part_flux, 0.0_num)
        CASE(-c_dir_z)
          ! negative flux in z
          part_flux = zfac * part_uz / gamma
          wdata = -ek * MIN(part_flux, 0.0_num)
        CASE( c_dir_z)
          ! positive flux in z
          part_flux = zfac * part_uz / gamma
          wdata =  ek * MAX(part_flux, 0.0_num)
        END SELECT

        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            data_array(cell_x+ix, cell_y+iy) = &
                data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * wdata
            wt(cell_x+ix, cell_y+iy) = &
                wt(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * l_weight
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    CALL processor_summation_bcs(wt)

    data_array = data_array / MAX(wt, c_non_zero)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekflux



  SUBROUTINE calc_poynt_flux(data_array, current_species, direction)

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    INTEGER :: ix, iy
    REAL(num) :: ex_cc, ey_cc, ez_cc, bx_cc, by_cc, bz_cc

    SELECT CASE(direction)
    CASE(c_dir_x)
      DO iy = 1, ny
        DO ix = 1, nx
          ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy))
          ez_cc = ez(ix, iy)
          by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy))
          bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                           +  bz(ix-1, iy  ) + bz(ix, iy  ))
          data_array(ix,iy) = (ey_cc * bz_cc - ez_cc * by_cc) / mu0
        ENDDO
      ENDDO
    CASE(c_dir_y)
      DO iy = 1, ny
        DO ix = 1, nx
          ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy))
          ez_cc = ez(ix, iy)
          bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy))
          bz_cc = 0.25_num * (bz(ix-1, iy-1) + bz(ix, iy-1) &
                           +  bz(ix-1, iy  ) + bz(ix, iy  ))
          data_array(ix,iy) = (ez_cc * bx_cc - ex_cc * bz_cc) / mu0
        ENDDO
      ENDDO
    CASE(c_dir_z)
      DO iy = 1, ny
        DO ix = 1, nx
          ex_cc = 0.5_num  * (ex(ix-1, iy  ) + ex(ix, iy))
          ey_cc = 0.5_num  * (ey(ix  , iy-1) + ey(ix, iy))
          bx_cc = 0.5_num  * (bx(ix  , iy-1) + bx(ix, iy))
          by_cc = 0.5_num  * (by(ix-1, iy  ) + by(ix, iy))
          data_array(ix,iy) = (ex_cc * by_cc - ey_cc * bx_cc) / mu0
        ENDDO
      ENDDO
    END SELECT

  END SUBROUTINE calc_poynt_flux



  SUBROUTINE calc_charge_density(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
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
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
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

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
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
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
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

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: l_weight
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
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
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
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
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
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
            sigma(cell_x+ix, cell_y+iy) = sigma(cell_x+ix, cell_y+iy) + gf &
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

    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / c_ndims

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    INTEGER :: ispecies, ix, iy, spec_start, spec_end
    TYPE(particle), POINTER :: current
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
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
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
