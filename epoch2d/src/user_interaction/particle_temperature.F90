! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE particle_temperature

  USE constants
  USE random_generator
  USE evaluator

  IMPLICIT NONE

  REAL(num), PARAMETER, PRIVATE :: max_average_its = 20.0_num

CONTAINS

  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_species, &
      drift)

    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_species), POINTER :: part_species
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: ix, iy
#include "particle_head.inc"

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

      ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO iy = sf_min, sf_max
        DO ix = sf_min, sf_max
          temp_local = temp_local &
              + gx(ix) * gy(iy) * temperature(cell_x+ix, cell_y+iy)
          drift_local = drift_local &
              + gx(ix) * gy(iy) * drift(cell_x+ix, cell_y+iy)
        END DO
      END DO

      IF (direction == c_dir_x) current%part_p(1) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction == c_dir_y) current%part_p(2) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      IF (direction == c_dir_z) current%part_p(3) = &
          momentum_from_temperature(mass, temp_local, drift_local)

      current => current%next
      ipart = ipart + 1
    END DO

  END SUBROUTINE setup_particle_temperature



  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature_relativistic(temperature, &
      part_species, drift)

    REAL(num), DIMENSION(1-ng:,1-ng:,:), INTENT(IN) :: temperature
    TYPE(particle_species), POINTER :: part_species
    REAL(num), DIMENSION(1-ng:,1-ng:,:), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass
    REAL(num), DIMENSION(c_ndirs) ::  temp_local, drift_local
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: ix, iy, idir
#include "particle_head.inc"

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

      ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO idir = 1, c_ndirs
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            temp_local(idir) = temp_local(idir) &
                + gx(ix) * gy(iy) * temperature(cell_x+ix, cell_y+iy, idir)
            drift_local(idir) = drift_local(idir) &
                + gx(ix) * gy(iy) * drift(cell_x+ix, cell_y+iy, idir)
          END DO
        END DO
      END DO

      current%part_p = momentum_from_temperature_relativistic(mass, &
          temp_local, part_species%fractional_tail_cutoff, drift_local)

      current => current%next
      ipart = ipart + 1
    END DO

  END SUBROUTINE setup_particle_temperature_relativistic



  ! Subroutine to initialise an arbitrary distribution function
  SUBROUTINE setup_particle_dist_fn(part_species, drift)

    TYPE(particle_species), POINTER :: part_species
    REAL(num), DIMENSION(1-ng:,1-ng:,:), INTENT(IN) :: drift
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass
    REAL(num), DIMENSION(c_ndirs) ::  drift_local
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart, it, ipart_global, it_global
    REAL(num) :: average_its
    INTEGER :: ix, iy, idir, err
    TYPE(parameter_pack) :: parameters
    REAL(num), DIMENSION(c_ndirs, 2) :: ranges
    CHARACTER(LEN=25) :: string
    INTEGER :: iu, io
#include "particle_head.inc"

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    it = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

      ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

      drift_local = 0.0_num
      DO idir = 1, c_ndirs
        DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
            drift_local(idir) = drift_local(idir) &
                + gx(ix) * gy(iy) * drift(cell_x+ix, cell_y+iy, idir)
          END DO
        END DO
      END DO

      parameters%use_grid_position = .FALSE.
      parameters%pack_ix = cell_x
      parameters%pack_pos = current%part_pos
      err = 0
      CALL evaluate_with_parameters_to_array(part_species%dist_fn_range(1), &
          parameters, 2, ranges(1,:), err)
      CALL evaluate_with_parameters_to_array(part_species%dist_fn_range(2), &
          parameters, 2, ranges(2,:), err)
      CALL evaluate_with_parameters_to_array(part_species%dist_fn_range(3), &
          parameters, 2, ranges(3,:), err)

      CALL sample_from_deck_expression(current, part_species%dist_fn, &
          parameters, ranges, mass, drift_local, it)

      current => current%next
      ipart = ipart + 1
    END DO

    CALL MPI_REDUCE(ipart, ipart_global, 1, MPI_INTEGER8, MPI_SUM, 0, comm, &
        errcode)
    CALL MPI_REDUCE(it, it_global, 1, MPI_INTEGER8, MPI_SUM, 0, comm, &
        errcode)
    average_its = REAL(it_global, num) / MAX(REAL(ipart_global, num), c_tiny)

    IF (rank == 0 .AND. average_its >= max_average_its) THEN
      WRITE(string,'(F8.1)') average_its
      DO iu = 1, nio_units ! Print to stdout and to file
        io = ios_units(iu)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'Setup distribution for particles of species ', '"' &
            // TRIM(part_species%name) // '"'
        WRITE(io,*) 'taking ' // TRIM(ADJUSTL(string)) &
            // ' iteratations per particle on average.'
        WRITE(io,*) 'Possibly try smaller momentum range'
      END DO
    END IF

  END SUBROUTINE setup_particle_dist_fn



  FUNCTION momentum_from_temperature_relativistic(mass, temperature, cutoff, &
      drift)

    REAL(num), INTENT(IN) :: mass
    REAL(num), DIMENSION(c_ndirs), INTENT(IN) :: temperature
    REAL(num), INTENT(IN) :: cutoff
    REAL(num), DIMENSION(c_ndirs), INTENT(IN) :: drift
    REAL(num), DIMENSION(c_ndirs) :: momentum_from_temperature_relativistic

    ! Three parameters for calculating the range of momenta
    ! Includes different combinations of physical constants
    REAL(num), PARAMETER :: param1 = -3.07236e-40_num
    REAL(num), PARAMETER :: param2 = 2.35985e-80_num
    ! c^2/kb
    REAL(num), PARAMETER :: c2_k = 6.509658203714208e39_num
    REAL(num) :: rand, probability
    REAL(num), DIMENSION(c_ndirs) :: momentum, mmc
    REAL(num) :: mod_momentum, mass_c, pfac, drift_2
    REAL(num) :: momentum1_2, momentum2_2, momentum3_2
    REAL(num) :: temp, p_max_x, p_max_y, p_max_z, p_max
    REAL(num) :: temp_norm1, temp_norm2, temp_norm3
    REAL(num) :: temp_fac1, temp_fac2, temp_fac3
    INTEGER :: dof
    REAL(num) :: inter1, inter2, inter3
    REAL(num) :: gamma_before, gamma_after, gamma_drift, gamma_drift_fac
    LOGICAL :: no_drift

    dof = COUNT(temperature > c_tiny)

    ! If there are no degrees of freedom them sampling is unnecessary
    ! plasma is cold
    IF (dof == 0) THEN
      momentum_from_temperature_relativistic = 0.0_num
      RETURN
    END IF

    temp = SUM(temperature) / REAL(dof, num)
    temp_norm1 = temperature(1) / temp
    temp_norm2 = temperature(2) / temp
    temp_norm3 = temperature(3) / temp
    temp_fac1 = 1.0_num / MAX(temp_norm1, c_tiny)
    temp_fac2 = 1.0_num / MAX(temp_norm2, c_tiny)
    temp_fac3 = 1.0_num / MAX(temp_norm3, c_tiny)
    mass_c = mass * c

    p_max = SQRT(param1 * mass * temp * LOG(cutoff) &
        + param2 * temp**2 * LOG(cutoff)**2) / mass

    p_max_x = p_max * SQRT(temp_norm1)
    p_max_y = p_max * SQRT(temp_norm2)
    p_max_z = p_max * SQRT(temp_norm3)

    pfac = -c2_k * mass / temp

    drift_2 = DOT_PRODUCT(drift, drift)
    IF (drift_2 < c_tiny) THEN
      no_drift = .TRUE.
    ELSE
      no_drift = .FALSE.
      gamma_drift = SQRT(1.0_num + drift_2 / mass_c**2)
      gamma_drift_fac = 0.5_num / gamma_drift
    END IF

    ! Loop around until a momentum is accepted for this particle
    DO
      ! Generate random x and y momenta between p_min and p_max
      momentum(1) = random() * 2.0_num * p_max_x - p_max_x
      momentum(2) = random() * 2.0_num * p_max_y - p_max_y
      momentum(3) = random() * 2.0_num * p_max_z - p_max_z

      momentum1_2 = momentum(1)**2
      momentum2_2 = momentum(2)**2
      momentum3_2 = momentum(3)**2

      mod_momentum = SQRT(momentum1_2 + momentum2_2 + momentum3_2)

      ! From that value, have to generate the probability that a particle
      ! with that momentum should be accepted.
      ! This is just the particle distribution function scaled to have
      ! a maximum of 1 (or lower).
      ! In general you will have to work this out yourself
      inter1 = momentum1_2 * temp_fac1
      inter2 = momentum2_2 * temp_fac2
      inter3 = momentum3_2 * temp_fac3
      probability = EXP(pfac * (SQRT(1.0_num + inter1 + inter2 + inter3) &
                                - 1.0_num))

      ! Once you know your probability you just generate a random number
      ! between 0 and 1 and if the generated number is less than the
      ! probability then accept the particle and exit this loop.
      rand = random()
      IF (rand > probability) CYCLE

      mmc = momentum * mass_c
      IF (no_drift) EXIT

      CALL drift_lorentz_transform(mmc, mass_c, drift, &
          gamma_before, gamma_after, gamma_drift)

      rand = random()
      IF (rand < gamma_drift_fac * (gamma_after / gamma_before)) EXIT
    END DO

    momentum_from_temperature_relativistic = mmc

  END FUNCTION momentum_from_temperature_relativistic



  ! Subroutine takes a rest frame momentum and a drift momentum and Lorentz
  ! transforms the momentum subject to the specified drift.
  SUBROUTINE drift_lorentz_transform(p, mass_c, drift, &
      gamma_before, gamma_after, gamma_drift)

    REAL(num), DIMENSION(c_ndirs), INTENT(INOUT) :: p
    REAL(num), INTENT(IN) :: mass_c
    REAL(num), DIMENSION(c_ndirs), INTENT(IN) :: drift
    REAL(num), INTENT(OUT) :: gamma_before, gamma_after
    REAL(num), INTENT(IN) :: gamma_drift
    REAL(num), DIMENSION(c_ndirs) :: p_mc, beta
    REAL(num), DIMENSION(c_ndirs+1) :: p4_in
    REAL(num), DIMENSION(c_ndirs,c_ndirs+1) :: boost_tensor
    REAL(num) :: e_prime, imc, gamma_m1_beta2, p2
    INTEGER :: i, j

    imc = 1.0_num / mass_c
    p_mc = p * imc
    gamma_before = SQRT(1.0_num + DOT_PRODUCT(p_mc, p_mc))
    e_prime = gamma_before * mass_c

    beta = -drift * imc / gamma_drift ! Lorentz beta vector

    gamma_m1_beta2 = (gamma_drift - 1.0_num) / DOT_PRODUCT(beta, beta)

    boost_tensor(1,1) = -beta(1) * gamma_drift
    boost_tensor(2,1) = -beta(2) * gamma_drift
    boost_tensor(3,1) = -beta(3) * gamma_drift

    boost_tensor(1,2) = 1.0_num + gamma_m1_beta2 * beta(1)**2
    boost_tensor(2,2) = gamma_m1_beta2 * beta(1) * beta(2)
    boost_tensor(3,2) = gamma_m1_beta2 * beta(1) * beta(3)

    boost_tensor(1,3) = gamma_m1_beta2 * beta(1) * beta(2)
    boost_tensor(2,3) = 1.0_num + gamma_m1_beta2 * beta(2)**2
    boost_tensor(3,3) = gamma_m1_beta2 * beta(2) * beta(3)

    boost_tensor(1,4) = gamma_m1_beta2 * beta(1) * beta(3)
    boost_tensor(2,4) = gamma_m1_beta2 * beta(2) * beta(3)
    boost_tensor(3,4) = 1.0_num + gamma_m1_beta2 * beta(3)**2

    p4_in = [e_prime, p(1), p(2), p(3)]
    p2 = 0.0_num

    DO i = 1, 3
      p(i) = 0.0_num
      DO j = 1, 4
        p(i) = p(i) + p4_in(j) * boost_tensor(i,j)
      END DO
      p2 = p2 + (p(i) * imc)**2
    END DO

    gamma_after = SQRT(1.0_num + p2)

  END SUBROUTINE drift_lorentz_transform



  FUNCTION momentum_from_temperature(mass, temperature, drift)

    REAL(num), INTENT(IN) :: mass, temperature, drift
    REAL(num) :: momentum_from_temperature
    DOUBLE PRECISION :: stdev, mu

    stdev = DBLE(SQRT(temperature * kb * mass))
    mu = DBLE(drift)
    momentum_from_temperature = random_box_muller(stdev, mu)

  END FUNCTION momentum_from_temperature



  ! Function for generating momenta of thermal particles in a particular
  ! direction, e.g. the +x direction.
  ! Samples distribution v f(v - drift) where f is Maxwellian
  ! including for drift == 0
  ! Drift is signed, and direction should be +/- 1 giving the sign of wanted
  ! velocity

  FUNCTION flux_momentum_from_temperature(mass, temperature, drift, direction)

    REAL(num), INTENT(IN) :: mass, temperature, drift, direction
    REAL(num) :: flux_momentum_from_temperature, mom1, mom2
    REAL(num) :: v, vt2, vt, vd, vmin, vmax, vrange, vexp
    REAL(num) :: norm, fac
    INTEGER :: i

    IF (ABS(drift) > c_tiny) THEN
      vt2 = kb * temperature / mass
      IF (vt2 < c_tiny) THEN
        flux_momentum_from_temperature = drift
        RETURN
      END IF

      vt = SQRT(vt2)
      vrange = 3.0_num * vt
      vd = direction * drift / mass
      vmax = vd + vrange

      IF (vmax < c_tiny) THEN
        flux_momentum_from_temperature = 0.0_num
        RETURN
      END IF

      vmin = MAX(vd - vrange, 0.0_num)
      vrange = vmax - vmin

      fac = -0.5_num / vt2
      vexp = 0.5_num * (vd + SQRT(vd**2 + 4.0_num * vt2))
      norm = 1.0_num / (vexp * EXP(fac * (vexp - vd)**2))

      DO i = 1, 1000
        v = vmin + random() * vrange

        IF (random() < norm * v * EXP(fac * (v - vd)**2)) THEN
          flux_momentum_from_temperature = direction * v * mass
          RETURN
        END IF
      END DO
    ELSE
      ! These satisfy a Rayleigh distribution, formed by combining two
      ! normally-distributed (~N(0,sigma)) random variables as follows:
      ! Z = SQRT(X**2 + Y**2)

      mom1 = momentum_from_temperature(mass, temperature, 0.0_num)
      mom2 = momentum_from_temperature(mass, temperature, 0.0_num)

      flux_momentum_from_temperature = direction * SQRT(mom1**2 + mom2**2)
    END IF

  END FUNCTION flux_momentum_from_temperature



  ! Function to take a deck expression and sample until it returns a value

  SUBROUTINE sample_from_deck_expression(part, stack, parameters, &
      ranges, mass, drift, it_r)

    TYPE(particle), INTENT(INOUT) :: part
    TYPE(primitive_stack), INTENT(INOUT) :: stack
    TYPE(parameter_pack), INTENT(INOUT) :: parameters
    REAL(num), DIMENSION(c_ndirs,2), INTENT(IN) :: ranges
    REAL(num), INTENT(IN) :: mass
    REAL(num), DIMENSION(c_ndirs) , INTENT(IN) :: drift
    INTEGER(i8), INTENT(INOUT), OPTIONAL :: it_r
    REAL(num) :: rand, probability, mass_c, drift_2
    REAL(num) :: gamma_before, gamma_after, gamma_drift, gamma_drift_fac
    REAL(num) :: range_diff1, range_diff2, range_diff3
    INTEGER :: err
    INTEGER(i8) :: it
    LOGICAL :: no_drift

    err = c_err_none
    IF (PRESENT(it_r)) it = it_r

    drift_2 = DOT_PRODUCT(drift, drift)
    IF (drift_2 < c_tiny) THEN
      no_drift = .TRUE.
    ELSE
      no_drift = .FALSE.
      mass_c = mass * c
      gamma_drift = SQRT(1.0_num + drift_2 / mass_c**2)
      gamma_drift_fac = 0.5_num / gamma_drift
    END IF

    range_diff1 = ranges(1,2) - ranges(1,1)
    range_diff2 = ranges(2,2) - ranges(2,1)
    range_diff3 = ranges(3,2) - ranges(3,1)

    DO
      ! These lines are setting global variables that are later used by
      ! the deck parser
      parameters%pack_p(1) = random() * range_diff1 + ranges(1,1)
      parameters%pack_p(2) = random() * range_diff2 + ranges(2,1)
      parameters%pack_p(3) = random() * range_diff3 + ranges(3,1)

      ! pack spatial information has already been set before calling
      probability = evaluate_with_parameters(stack, parameters, err)
      IF (err /= c_err_none .AND. rank == 0) THEN
        PRINT*, 'Unable to evaluate distribution function'
        CALL abort_code(c_err_bad_setup)
      END IF

      it = it + 1

      rand = random()
      IF (rand > probability) CYCLE

      IF (no_drift) EXIT

      CALL drift_lorentz_transform(parameters%pack_p, mass_c, drift, &
          gamma_before, gamma_after, gamma_drift)

      rand = random()
      IF (rand < gamma_drift_fac * (gamma_after / gamma_before)) EXIT
    END DO

    IF (PRESENT(it_r)) it_r = it

    part%part_p = parameters%pack_p

  END SUBROUTINE sample_from_deck_expression

END MODULE particle_temperature
