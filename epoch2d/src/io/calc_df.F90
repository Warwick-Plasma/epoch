MODULE calc_df

  USE shared_data
  USE boundary
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_mass_density(data_array, current_species)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_px, part_py, part_pz, part_q, part_m

    ! Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! The weight of a particle
    REAL(num) :: l_weight
    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy
    ! The data to be weighted onto the grid
    REAL(num) :: DATA

    REAL(num), DIMENSION(-2:, -2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, spec_start, spec_end

    data_array = 0.0_num

    l_weight = weight

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
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
        part_m  = current%mass
#else
        part_q  = particle_species(ispecies)%charge
        part_m  = particle_species(ispecies)%mass
#endif

#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

        cell_x_r = part_x / dx
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            DATA = part_m * l_weight / (dx*dy)
            data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_ekbar(data_array, current_species)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_px, part_py, part_pz, part_q, part_m

    ! Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! The weight of a particle
    REAL(num) :: l_weight

    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy
    ! The data to be weighted onto the grid
    REAL(num) :: DATA

    REAL(num), DIMENSION(-2:, -2:), INTENT(INOUT) :: data_array
    REAL(num), DIMENSION(:, :), ALLOCATABLE :: ct
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, spec_start, spec_end

    ALLOCATE(ct(-2:nx+3, -2:ny+3))
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
      DO WHILE (ASSOCIATED(current))

        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
        part_m  = current%mass
#else
        part_q  = particle_species(ispecies)%charge
        part_m  = particle_species(ispecies)%mass
#endif

#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

        cell_x_r = part_x / dx
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        IF (cell_y .LT. -2) PRINT *, cell_y, rank, part_y

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            DATA = SQRT(((part_px*l_weight)**2+(part_py*l_weight)**2+(part_pz*l_weight)**2)*c**2 + (part_m*l_weight)**2*c**4) - (part_m*l_weight)*c**2
            data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
            ct(cell_x+ix, cell_y+iy) = ct(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * l_weight
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
    INTEGER :: cell_x, cell_y

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_px, part_py, part_pz, part_q, part_m

    ! Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! The weight of a particle
    REAL(num) :: l_weight

    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy
    ! The data to be weighted onto the grid
    REAL(num) :: DATA

    REAL(num), DIMENSION(-2:, -2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, spec_start, spec_end

    data_array = 0.0_num

    l_weight = weight

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
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
        part_m  = current%mass
#else
        part_q  = particle_species(ispecies)%charge
        part_m  = particle_species(ispecies)%mass
#endif

#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

        cell_x_r = part_x / dx
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            DATA = part_q * l_weight / (dx*dy)
            data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
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
    INTEGER :: cell_x, cell_y

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_px, part_py, part_pz, part_q, part_m

    ! Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! The weight of a particle
    REAL(num) :: l_weight

    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy
    ! The data to be weighted onto the grid
    REAL(num) :: DATA

    REAL(num), DIMENSION(-2:, -2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, spec_start, spec_end

    data_array = 0.0_num

    l_weight = weight

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
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
        part_m  = current%mass
#else
        part_q  = particle_species(ispecies)%charge
        part_m  = particle_species(ispecies)%mass
#endif

#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

        cell_x_r = part_x / dx
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            DATA = 1.0_num * l_weight / (dx*dy)
            data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
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
    INTEGER :: cell_x, cell_y

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_px, part_py, part_pz, part_q, part_m

    ! Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! The weight of a particle
    REAL(num) :: l_weight

    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy
    ! The data to be weighted onto the grid
    REAL(num) :: DATA

    REAL(num), DIMENSION(-2:, -2:), INTENT(INOUT) :: data_array
    REAL(num), DIMENSION(:, :), ALLOCATABLE ::  part_count, mass, sigma, mean
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, spec_start, spec_end

    data_array = 0.0_num

    l_weight = weight

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

!!$    ALLOCATE(p_max(-2:nx+3, -2:ny+3), p_min(-2:nx+3, -2:ny+3))
    ALLOCATE(mean(-2:nx+3, -2:ny+3), part_count(-2:nx+3, -2:ny+3))
    ALLOCATE(mass(-2:nx+3, -2:ny+3), sigma(-2:nx+3, -2:ny+3))
    data_array = 0.0_num
    mean = 0.0_num
    part_count = 0.0_num
    mass = 0.0_num
    sigma = 0.0_num

!!$    p_min = 1.0e13_num
!!$    p_max = -1.0e13_num

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
        part_m  = current%mass
#else
        part_q  = particle_species(ispecies)%charge
        part_m  = particle_species(ispecies)%mass
#endif
#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

        cell_x_r = part_x / dx !
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy !
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            DATA = SQRT(part_px**2+part_py**2+part_pz**2) * l_weight
!!$                IF (data .GE. p_min(cell_x+ix, cell_y+iy) .AND. data .LE. p_max(cell_x+ix, cell_y+iy)) THEN
            mean(cell_x+ix, cell_y+iy) = mean(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
            DATA = l_weight
            part_count(cell_x+ix, cell_y+iy) = part_count(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
            DATA = part_m * l_weight
            mass(cell_x+ix, cell_y+iy) = mass(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
!!$                ENDIF
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

    mean = mean/part_count
    mass = mass/part_count

    WHERE (part_count .LT. 1.0_num) mean = 0.0_num
    WHERE (part_count .LT. 1.0_num) mass = 0.0_num

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)

        cell_x_r = part_x / dx !
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy !
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            DATA = SQRT(part_px**2+part_py**2+part_pz**2)
!!$                IF (data .GE. p_min(cell_x+ix, cell_y+iy) .AND. data .LE. p_max(cell_x+ix, cell_y+iy)) THEN
            sigma(cell_x+ix, cell_y+iy) = sigma(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * (DATA-mean(cell_x+ix, cell_y+iy))**2
            part_count(cell_x+ix, cell_y+iy) = part_count(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy)
!!$                ENDIF
          ENDDO
        ENDDO
        current=>current%next
      ENDDO
    ENDDO
    CALL processor_summation_bcs(sigma)
    CALL field_zero_gradient(sigma, .TRUE.)
    CALL processor_summation_bcs(part_count)
    CALL field_zero_gradient(part_count, .TRUE.)
    sigma = sigma/MAX(part_count, 0.1_num)
    WHERE (part_count .LT. 1.0) sigma = 0.0_num

    data_array = sigma/(kb * MAX(mass, m0))

    DEALLOCATE(part_count, mean, mass, sigma)
!!$    DEALLOCATE(p_min, p_max)

    CALL field_bc(data_array)
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x, cell_y

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y

    ! Contains the floating point version of the cell number (never actually used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! The weight of a particle
    REAL(num) :: l_weight

    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy
    ! The data to be weighted onto the grid
    REAL(num) :: DATA

    REAL(num), DIMENSION(-2:, -2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, spec_start, spec_end

    INTERFACE
      FUNCTION evaluator(a_particle, species_eval)
        USE shared_data
        TYPE(particle), POINTER :: a_particle
        INTEGER, INTENT(IN) :: species_eval
        REAL(num) :: evaluator
      END FUNCTION evaluator
    END INTERFACE

    data_array = 0.0_num

    l_weight = weight

    spec_start = current_species
    spec_end = current_species

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
    ENDIF

    DO ispecies = spec_start, spec_end
      current=>particle_species(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))

        part_x  = current%part_pos(1) - x_start_local
        part_y  = current%part_pos(2) - y_start_local

        cell_x_r = part_x / dx !
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy !
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)

        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            DATA = evaluator(current, ispecies)
            data_array(cell_x+ix, cell_y+iy) = data_array(cell_x+ix, cell_y+iy) + gx(ix) * gy(iy) * DATA
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)
    CALL field_zero_gradient(data_array, .TRUE.)

  END SUBROUTINE calc_on_grid_with_evaluator

END MODULE calc_df
