MODULE particle_temperature

  USE random_generator
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_species, &
      drift)

    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_species), POINTER :: part_species
    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    TYPE(particle), POINTER :: current
    INTEGER(KIND=8) :: ipart
    INTEGER :: ix, iy
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    INTEGER :: cell_x, cell_y

    partlist=>part_species%attached_list
    current=>partlist%head
    ipart = 0
    DO WHILE(ipart .LT. partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

      ! Assume that temperature is cell centred
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

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          temp_local = temp_local &
              + gx(ix) * gy(iy) * temperature(cell_x+ix, cell_y+iy)
          drift_local = drift_local &
              + gx(ix) * gy(iy) * drift(cell_x+ix, cell_y+iy)
        ENDDO
      ENDDO

      IF (direction .EQ. c_dir_x) current%part_p(1) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction .EQ. c_dir_y) current%part_p(2) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction .EQ. c_dir_z) current%part_p(3) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      current=>current%next
      ipart = ipart + 1
    ENDDO

  END SUBROUTINE setup_particle_temperature



  FUNCTION momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: momentum_from_temperature

    REAL(num) :: stdev
    REAL(num) :: rand1, rand2, w
    REAL(num), SAVE :: val
    LOGICAL, SAVE :: cached = .FALSE.

    ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers
    ! The standard deviation (stdev) is related to temperature

    stdev = SQRT(temperature * kb * mass)

    IF (cached) THEN
      cached = .FALSE.
      momentum_from_temperature = val * stdev + drift
    ELSE
      cached = .TRUE.

      DO
        rand1 = random()
        rand2 = random()

        rand1 = 2.0_num * rand1 - 1.0_num
        rand2 = 2.0_num * rand2 - 1.0_num

        w = rand1**2 + rand2**2

        IF (w .LT. 1.0_num) EXIT
      ENDDO

      w = SQRT((-2.0_num * LOG(w)) / w)

      momentum_from_temperature = rand1 * w * stdev + drift
      val = rand2 * w
    ENDIF

  END FUNCTION momentum_from_temperature

END MODULE particle_temperature
