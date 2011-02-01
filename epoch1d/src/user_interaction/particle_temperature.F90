MODULE particle_temperature

  USE shape_functions

  IMPLICIT NONE

  REAL(num), SAVE :: max_rand

CONTAINS

  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_family, &
      drift, idum)

    REAL(num), DIMENSION(-2:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_family), POINTER :: part_family
    REAL(num), DIMENSION(-2:), INTENT(IN) :: drift
    INTEGER, INTENT(INOUT) :: idum
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num), DIMENSION(sf_min:sf_max) :: gx
    TYPE(particle), POINTER :: current
    INTEGER :: cell_x
    INTEGER(KIND=8) :: ipart
    INTEGER :: ix

    partlist=>part_family%attached_list
    current=>partlist%head
    ipart = 0
    DO WHILE(ipart .LT. partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_family%mass
#endif

      ! Assume that temperature is cell centred
#ifdef PARTICLE_SHAPE_TOPHAT
      cell_x_r = (current%part_pos - x_min_local) / dx - 1.0_num
#else
      cell_x_r = (current%part_pos - x_min_local) / dx - 0.5_num
#endif
      cell_x = FLOOR(cell_x_r + 0.5_num)
      cell_frac_x = REAL(cell_x, num) - cell_x_r
      cell_x = cell_x + 1

      CALL grid_to_particle(cell_frac_x, gx)

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO ix = sf_min, sf_max
        temp_local = temp_local + gx(ix) * temperature(cell_x+ix)
        drift_local = drift_local + gx(ix) * drift(cell_x+ix)
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
    REAL(num), SAVE :: val
    LOGICAL, SAVE :: cached = .FALSE.

    ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers
    ! The standard deviation (stdev) is related to temperature

    IF (cached) THEN
      cached = .FALSE.
      momentum_from_temperature = val
    ELSE
      cached = .TRUE.
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
      val = rand2 * w * stdev + drift
    ENDIF

  END FUNCTION momentum_from_temperature



  FUNCTION old_random(idum)

    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: old_random
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

    old_random = am*idum
    idum = IEOR(idum, mask)

    IF (old_random .GT. max_rand) max_rand = old_random

  END FUNCTION old_random



  FUNCTION random(idum)

    INTEGER :: idum
    REAL :: random
    INTEGER, PARAMETER :: mbig = 1000000000, mseed = 161803398, mz = 0
    REAL(dbl), PARAMETER :: fac = 1.d0/mbig
    INTEGER :: i, ii, k, mj, mk
    INTEGER, SAVE :: iff = 0, ma(55)
    INTEGER, SAVE :: inext, inextp

    IF (idum .LT. 0 .OR. iff .EQ. 0) THEN
      iff = 1
      mj = ABS(mseed - ABS(idum))
      ma(55) = mj
      mk = 1
      DO i = 1,54
        ii = mod(21 * i, 55)
        ma(ii) = mk
        mk = mj - mk
        IF (mk .LT. mz) mk = mk + mbig
        mj = ma(ii)
      ENDDO

      DO k = 1,4
        DO i = 1,55
          ma(i) = ma(i) - ma(1 + MOD(i + 30, 55))
          IF (ma(i) .LT. mz) ma(i) = ma(i) + mbig
        ENDDO
      ENDDO

      inext = 0
      inextp = 31
      idum = 1
    ENDIF

    inext = inext + 1
    IF (inext .EQ. 56) inext = 1
    inextp = inextp + 1
    IF (inextp .EQ. 56) inextp = 1
    mj = ma(inext) - ma(inextp)
    IF (mj .LT. mz) mj = mj + mbig
    ma(inext) = mj

    random = mj * fac

    IF (random .GT. max_rand) max_rand = random

  END FUNCTION random

END MODULE particle_temperature
