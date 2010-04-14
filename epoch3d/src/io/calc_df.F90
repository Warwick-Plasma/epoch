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

    ! particle weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: data

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
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
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
#ifdef PER_PARTICLE_CHARGEMASS
        part_m  = current%mass
#else
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

        cell_z_r = part_z / dz
        cell_z  = NINT(cell_z_r)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              data = part_m * l_weight / (dx*dy*dz)
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * data
            ENDDO
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)

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

    ! Particle Weight factors as described in the manual (FIXREF)
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz

    ! The data to be weighted onto the grid
    REAL(num) :: data

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ct
    INTEGER, INTENT(IN) :: current_species

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies, spec_start, spec_end

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
#else
        part_m  = particle_species(ispecies)%mass
#endif

#ifdef PER_PARTICLE_WEIGHT
        l_weight = current%weight
#endif

        cell_x_r = part_x / dx
        cell_x = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy
        cell_y = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        cell_z_r = part_z / dz
        cell_z = NINT(cell_z_r)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              data = SQRT(((part_px*l_weight)**2 + (part_py*l_weight)**2 &
                  + (part_pz*l_weight)**2)*c**2 + (part_m*l_weight)**2*c**4) &
                  - (part_m*l_weight)*c**2
              data_array(cell_x+ix,cell_y+iy,cell_z+iz) = &
                  data_array(cell_x+ix,cell_y+iy,cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * data
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
    REAL(num) :: data

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
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
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
#ifdef PER_PARTICLE_CHARGEMASS
        part_q  = current%charge
#else
        part_q  = particle_species(ispecies)%charge
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

        cell_z_r = part_z / dz
        cell_z  = NINT(cell_z_r)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              data = part_q * l_weight / (dx*dy*dz)
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * data
            ENDDO
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)

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
    REAL(num) :: data

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
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
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local

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

        cell_z_r = part_z / dz
        cell_z  = NINT(cell_z_r)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              data = l_weight / (dx*dy*dz)
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * data
            ENDDO
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)

  END SUBROUTINE calc_number_density



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

    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    ! The data to be weighted onto the grid
    REAL(num) :: data

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: data_array
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

        cell_x_r = part_x / dx
        cell_x  = NINT(cell_x_r)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_x = cell_x+1

        cell_y_r = part_y / dy
        cell_y  = NINT(cell_y_r)
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_y = cell_y+1

        cell_z_r = part_z / dz
        cell_z  = NINT(cell_z_r)
        cell_frac_z = REAL(cell_z, num) - cell_z_r
        cell_z = cell_z+1

        CALL particle_to_grid(cell_frac_x, gx)
        CALL particle_to_grid(cell_frac_y, gy)
        CALL particle_to_grid(cell_frac_z, gz)

        DO iz = -sf_order, sf_order
          DO iy = -sf_order, sf_order
            DO ix = -sf_order, sf_order
              data = evaluator(current, ispecies)
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * data
            ENDDO
          ENDDO
        ENDDO

        current=>current%next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array)

  END SUBROUTINE calc_on_grid_with_evaluator

END MODULE calc_df
