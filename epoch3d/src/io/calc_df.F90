MODULE calc_df

  USE boundary
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_mass_density(data_array, current_species)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z, part_m

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
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    REAL(num) :: fac, idxyz

    data_array = 0.0_num

    l_weight = weight
    idxyz = 1.0_num / dx / dy / dz
    fac   = weight  / dx / dy / dz

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGEMASS
      part_m  = particle_species(ispecies)%mass
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
#ifdef PER_PARTICLE_CHARGEMASS
        part_m  = current%mass
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
        fac = l_weight * idxyz
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
        cell_z_r = part_z / dz - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
        cell_z_r = part_z / dz
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y + 1

        cell_z = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = part_m * fac
        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
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
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_ekbar(data_array, current_species)

    ! Contains the integer cell position of the particle in x,y,z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z, part_px, part_py, part_pz, part_m

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
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
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

    l_weight = weight

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGEMASS
      part_m  = particle_species(ispecies)%mass
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
        part_m  = current%mass
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
        cell_z_r = part_z / dz - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
        cell_z_r = part_z / dz
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y + 1

        cell_z = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = (SQRT(part_px**2 + part_py**2 + part_pz**2 + (part_m*c)**2) &
          - part_m * c) * c * l_weight
        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              data_array(cell_x+ix,cell_y+iy,cell_z+iz) = &
                  data_array(cell_x+ix,cell_y+iy,cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
              ct(cell_x+ix,cell_y+iy,cell_z+iz) = &
                  ct(cell_x+ix,cell_y+iy,cell_z+iz) &
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
    CALL field_zero_gradient(data_array, .TRUE.)

    DEALLOCATE(ct)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_charge_density(data_array, current_species)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z, part_q

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
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    REAL(num) :: fac, idxyz

    data_array = 0.0_num

    l_weight = weight
    idxyz = 1.0_num / dx / dy / dz
    fac   = weight  / dx / dy / dz

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGEMASS
      part_q  = particle_species(ispecies)%charge
#endif
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
        fac = l_weight * idxyz
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
        cell_z_r = part_z / dz - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
        cell_z_r = part_z / dz
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y + 1

        cell_z = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = part_q * fac
        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
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
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_charge_density



  SUBROUTINE calc_number_density(data_array, current_species)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z

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
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    REAL(num) :: idxyz

    data_array = 0.0_num

    l_weight = weight
    idxyz = 1.0_num / dx / dy / dz
    wdata = weight  / dx / dy / dz

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local

#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
        wdata = l_weight * idxyz
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
        cell_z_r = part_z / dz - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
        cell_z_r = part_z / dz
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y + 1

        cell_z = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
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
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_number_density



  SUBROUTINE calc_temperature(data_array, current_species)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z, part_px, part_py, part_pz, part_m

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
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: wdata1, wdata2

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE ::  part_count, mass, sigma, mean
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end

    data_array = 0.0_num

    l_weight = weight

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    ALLOCATE(mean(-2:nx+3,-2:ny+3,-2:nz+3), part_count(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(mass(-2:nx+3,-2:ny+3,-2:nz+3), sigma(-2:nx+3,-2:ny+3,-2:nz+3))
    data_array = 0.0_num
    mean = 0.0_num
    part_count = 0.0_num
    mass = 0.0_num
    sigma = 0.0_num

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGEMASS
      part_m  = particle_species(ispecies)%mass
#endif
      DO WHILE(ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
        part_m  = current%mass
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
        cell_z_r = part_z / dz - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
        cell_z_r = part_z / dz
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y + 1

        cell_z = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata1 = SQRT(part_px**2 + part_py**2 + part_pz**2) * l_weight
        wdata2 = part_m * l_weight
        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              mean(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  mean(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata1
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  part_count(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * l_weight
              mass(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  mass(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata2
            ENDDO
          ENDDO
        ENDDO
        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(mean)
    CALL field_zero_gradient(mean, .TRUE.)
    CALL processor_summation_bcs(part_count)
    CALL field_zero_gradient(part_count, .TRUE.)
    CALL processor_summation_bcs(mass)
    CALL field_zero_gradient(mass, .TRUE.)

    mean = mean / part_count
    mass = mass / part_count

    WHERE (part_count .LT. 1.0_num) mean = 0.0_num
    WHERE (part_count .LT. 1.0_num) mass = 0.0_num

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
        cell_z_r = part_z / dz - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
        cell_z_r = part_z / dz
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y + 1

        cell_z = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata1 = SQRT(part_px**2 + part_py**2 + part_pz**2)
        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              sigma(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  sigma(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) &
                  * (wdata1 - mean(cell_x+ix, cell_y+iy, cell_z+iz))**2
              part_count(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  part_count(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz)
            ENDDO
          ENDDO
        ENDDO
        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(sigma)
    CALL field_zero_gradient(sigma, .TRUE.)
    CALL processor_summation_bcs(part_count)
    CALL field_zero_gradient(part_count, .TRUE.)

    sigma = sigma / MAX(part_count, 0.1_num)
    WHERE (part_count .LT. 1.0) sigma = 0.0_num

    data_array = sigma / (kb * MAX(mass, m0))

    DEALLOCATE(part_count, mean, mass, sigma)

    CALL field_bc(data_array)
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y, cell_z

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
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
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x / dx - 0.5_num
        cell_y_r = part_y / dy - 0.5_num
        cell_z_r = part_z / dz - 0.5_num
#else
        cell_x_r = part_x / dx
        cell_y_r = part_y / dy
        cell_z_r = part_z / dz
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x + 1

        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y + 1

        cell_z = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z + 1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        wdata = evaluator(current, ispecies)
        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
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
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_on_grid_with_evaluator

END MODULE calc_df
