MODULE calc_df

  USE boundary
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_mass_density(data_array, current_species)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    REAL(num) :: fac, idx

    data_array = 0.0_num

    idx = 1.0_num / dx / dy / dz

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_m  = particle_species(ispecies)%mass
#endif
#ifndef PER_PARTICLE_WEIGHT
      fac = particle_species(ispecies)%weight * idx
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
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.0_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.0_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.0_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.5_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.5_num
#endif
        cell_x = FLOOR(cell_x_r)
        cell_y = FLOOR(cell_y_r)
        cell_z = FLOOR(cell_z_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r + 0.5_num
        cell_frac_y = REAL(cell_y, num) - cell_y_r + 0.5_num
        cell_frac_z = REAL(cell_z, num) - cell_z_r + 0.5_num

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = part_m * fac
        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
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

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_px, part_py, part_pz, part_mc

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! The weight of a particle
    REAL(num) :: l_weight, l_weightc

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ct
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end

    ALLOCATE(ct(-2:nx+3,-2:ny+3,-2:nz+3))
    data_array = 0.0_num
    ct = 0.0_num

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc = c * particle_species(ispecies)%mass
#endif
#ifndef PER_PARTICLE_WEIGHT
      l_weight = particle_species(ispecies)%weight
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
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.0_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.0_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.0_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.5_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.5_num
#endif
        cell_x = FLOOR(cell_x_r)
        cell_y = FLOOR(cell_y_r)
        cell_z = FLOOR(cell_z_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r + 0.5_num
        cell_frac_y = REAL(cell_y, num) - cell_y_r + 0.5_num
        cell_frac_z = REAL(cell_z, num) - cell_z_r + 0.5_num

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = (SQRT(part_px**2 + part_py**2 + part_pz**2 + part_mc**2) &
          - part_mc) * l_weightc
        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
              ct(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  ct(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * l_weight
            ENDDO
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

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    REAL(num) :: fac, idx

    data_array = 0.0_num

    idx = 1.0_num / dx / dy / dz

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q  = particle_species(ispecies)%charge
#endif
#ifndef PER_PARTICLE_WEIGHT
      fac = particle_species(ispecies)%weight * idx
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
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.0_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.0_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.0_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.5_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.5_num
#endif
        cell_x = FLOOR(cell_x_r)
        cell_y = FLOOR(cell_y_r)
        cell_z = FLOOR(cell_z_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r + 0.5_num
        cell_frac_y = REAL(cell_y, num) - cell_y_r + 0.5_num
        cell_frac_z = REAL(cell_z, num) - cell_z_r + 0.5_num

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = part_q * fac
        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
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

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    REAL(num) :: idx

    data_array = 0.0_num

    idx   = 1.0_num / dx / dy / dz

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_WEIGHT
      wdata = particle_species(ispecies)%weight * idx
#endif
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_WEIGHT
        wdata = current%weight * idx
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.0_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.0_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.0_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.5_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.5_num
#endif
        cell_x = FLOOR(cell_x_r)
        cell_y = FLOOR(cell_y_r)
        cell_z = FLOOR(cell_z_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r + 0.5_num
        cell_frac_y = REAL(cell_y, num) - cell_y_r + 0.5_num
        cell_frac_z = REAL(cell_z, num) - cell_z_r + 0.5_num

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: sigma
    INTEGER, INTENT(IN) :: current_species

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! The weight of a particle
    REAL(num) :: l_weight

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    REAL(num) :: gf

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    ALLOCATE(meanx(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(meany(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(meanz(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(part_count(-2:nx+3,-2:ny+3,-2:nz+3))
    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      sqrt_part_m  = SQRT(particle_species(ispecies)%mass)
#endif
#ifndef PER_PARTICLE_WEIGHT
      l_weight = particle_species(ispecies)%weight
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
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.0_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.0_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.0_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.5_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.5_num
#endif
        cell_x = FLOOR(cell_x_r)
        cell_y = FLOOR(cell_y_r)
        cell_z = FLOOR(cell_z_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r + 0.5_num
        cell_frac_y = REAL(cell_y, num) - cell_y_r + 0.5_num
        cell_frac_z = REAL(cell_z, num) - cell_z_r + 0.5_num

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * gz(iz) * l_weight
              meanx(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  meanx(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmx
              meany(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  meany(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmy
              meanz(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  meanz(cell_x+ix, cell_y+iy, cell_z+iz) + gf * part_pmz
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  part_count(cell_x+ix, cell_y+iy, cell_z+iz) + gf
            ENDDO
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
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      sqrt_part_m  = SQRT(particle_species(ispecies)%mass)
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
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.0_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.0_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.0_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.5_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.5_num
#endif
        cell_x = FLOOR(cell_x_r)
        cell_y = FLOOR(cell_y_r)
        cell_z = FLOOR(cell_z_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r + 0.5_num
        cell_frac_y = REAL(cell_y, num) - cell_y_r + 0.5_num
        cell_frac_z = REAL(cell_z, num) - cell_z_r + 0.5_num

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              gf = gx(ix) * gy(iy) * gz(iz)
              sigma(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  sigma(cell_x+ix, cell_y+iy, cell_z+iz) + gf &
                  * ((part_pmx - meanx(cell_x+ix, cell_y+iy, cell_z+iz))**2 &
                  + (part_pmy - meany(cell_x+ix, cell_y+iy, cell_z+iz))**2 &
                  + (part_pmz - meanz(cell_x+ix, cell_y+iy, cell_z+iz))**2)
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  part_count(cell_x+ix, cell_y+iy, cell_z+iz) + gf
            ENDDO
          ENDDO
        ENDDO
        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(sigma)
    CALL processor_summation_bcs(part_count)

    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / 3.0_num

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end

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
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.0_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.0_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.0_num
#else
        cell_x_r = (current%part_pos(1) - x_min_local) / dx + 1.5_num
        cell_y_r = (current%part_pos(2) - y_min_local) / dy + 1.5_num
        cell_z_r = (current%part_pos(3) - z_min_local) / dz + 1.5_num
#endif
        cell_x = FLOOR(cell_x_r)
        cell_y = FLOOR(cell_y_r)
        cell_z = FLOOR(cell_z_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r + 0.5_num
        cell_frac_y = REAL(cell_y, num) - cell_y_r + 0.5_num
        cell_frac_z = REAL(cell_z, num) - cell_z_r + 0.5_num

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = evaluator(current, ispecies)
        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
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
