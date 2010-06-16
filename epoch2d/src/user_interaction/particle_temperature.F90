MODULE particle_temperature

  USE shape_functions

  IMPLICIT NONE

  REAL(num), SAVE :: max_rand

CONTAINS

  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_family, &
      drift, idum)

    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_family), POINTER :: part_family
    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: drift
    INTEGER, INTENT(INOUT) :: idum
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num), DIMENSION(-2:2) :: gx, gy
    TYPE(particle), POINTER :: current
    INTEGER :: cell_x, cell_y
    INTEGER(KIND=8) :: ipart
    INTEGER :: ix, iy

    partlist=>part_family%attached_list
    current=>partlist%head
    ipart = 0
    DO WHILE(ipart .LT. partlist%count)
#ifdef PER_PARTICLE_CHARGEMASS
      mass = current%mass
#else
      mass = part_family%mass
#endif

      ! Assume that temperature is cell centred
#ifdef PARTICLE_SHAPE_TOPHAT
      cell_x_r = (current%part_pos(1) - x_min_local) / dx - 1.0_num
      cell_y_r = (current%part_pos(2) - y_min_local) / dy - 1.0_num
#else
      cell_x_r = (current%part_pos(1) - x_min_local) / dx - 0.5_num
      cell_y_r = (current%part_pos(2) - y_min_local) / dy - 0.5_num
#endif
      cell_x = FLOOR(cell_x_r + 0.5_num)
      cell_frac_x = REAL(cell_x, num) - cell_x_r
      cell_x = cell_x + 1

      cell_y = FLOOR(cell_y_r + 0.5_num)
      cell_frac_y = REAL(cell_y, num) - cell_y_r
      cell_y = cell_y + 1

      CALL grid_to_particle(cell_frac_x, gx)
      CALL grid_to_particle(cell_frac_y, gy)

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO iy = -sf_order, sf_order
        DO ix = -sf_order, sf_order
          temp_local = temp_local &
              + gx(ix) * gy(iy) * temperature(cell_x+ix, cell_y+iy)
          drift_local = drift_local &
              + gx(ix) * gy(iy) * drift(cell_x+ix, cell_y+iy)
        ENDDO
      ENDDO

      IF (IAND(direction, c_dir_x) .NE. 0) current%part_p(1) = &
          momentum_from_temperature(mass, temp_local, drift_local, idum)

      IF (IAND(direction, c_dir_y) .NE. 0) current%part_p(2) = &
          momentum_from_temperature(mass, temp_local, drift_local, idum)

      IF (IAND(direction, c_dir_z) .NE. 0) current%part_p(3) = &
          momentum_from_temperature(mass, temp_local, drift_local, idum)

      current=>current%next
      ipart = ipart + 1
    ENDDO

  END SUBROUTINE setup_particle_temperature



  FUNCTION momentum_from_temperature(mass, temperature, drift, idum)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: momentum_from_temperature

    REAL(num) :: stdev
    REAL(num) :: rand1, rand2, w

    ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers
    ! The standard deviation (stdev) is related to temperature

    stdev = SQRT(temperature * kb * mass)

    DO
      rand1 = random(idum)
      rand2 = random(idum)

      rand1 = 2.0_num * rand1 - 1.0_num
      rand2 = 2.0_num * rand2 - 1.0_num

      w = rand1**2 + rand2**2

      IF (w .LT. 1.0_num) EXIT
    ENDDO

    w = SQRT((-2.0_num * LOG(w)) / w)

    momentum_from_temperature = rand1 * w * stdev + drift

  END FUNCTION momentum_from_temperature



  FUNCTION random(idum)

    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: random
    INTEGER, PARAMETER :: ia = 16807, im = 2147483647, iq = 127773
    INTEGER, PARAMETER :: ir = 2836, mask = 123459876
    REAL(dbl), PARAMETER :: am = 1.0_dbl / 2147483647.0_dbl

    INTEGER :: k

    idum = IEOR(idum, mask)
    k = idum / iq

    idum = ia*(idum-k*iq)-ir*k
    IF (idum .LT. 0) THEN
      idum = idum+im
    ENDIF

    random = am*idum
    idum = IEOR(idum, mask)

    IF (random .GT. max_rand) max_rand = random

  END FUNCTION random

END MODULE particle_temperature
