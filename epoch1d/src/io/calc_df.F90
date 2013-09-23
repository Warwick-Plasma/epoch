MODULE calc_df

  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_boundary(data_array)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array

    CALL processor_summation_bcs(data_array, ng)

    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic .AND. x_min_boundary) THEN
      data_array(1) = data_array(1) + data_array( 0)
      data_array(2) = data_array(2) + data_array(-1)
      data_array(3) = data_array(3) + data_array(-2)
    ENDIF
    IF (bc_particle(c_bd_x_max) .NE. c_bc_periodic .AND. x_max_boundary) THEN
      data_array(nx-2) = data_array(nx-2) + data_array(nx+3)
      data_array(nx-1) = data_array(nx-1) + data_array(nx+2)
      data_array(nx  ) = data_array(nx  ) + data_array(nx+1)
    ENDIF

  END SUBROUTINE calc_boundary



  SUBROUTINE calc_mass_density(data_array, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_m = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type .EQ. c_species_id_photon) CYCLE
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_m  = io_list(ispecies)%mass
      fac = io_list(ispecies)%weight * idx
      wdata = part_m * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_m  = current%mass
#ifdef PER_PARTICLE_WEIGHT
        fac = current%weight * idx
#endif
        wdata = part_m * fac
#else
#ifdef PER_PARTICLE_WEIGHT
        fac = current%weight * idx
        wdata = part_m * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * wdata
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_mass_density



  SUBROUTINE calc_ekbar(data_array, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: l_weight
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma
    REAL(num), DIMENSION(:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(-2:nx+3))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    l_weight = 1.0_num

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      l_weight = io_list(ispecies)%weight
      fac = part_mc * l_weight * c

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

#ifdef PHOTONS
        IF (io_list(ispecies)%species_type .NE. c_species_id_photon) THEN
#endif
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          wdata = (gamma - 1.0_num) * fac
#ifdef PHOTONS
        ELSE
          wdata = current%particle_energy * l_weight
        ENDIF
#endif

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * wdata
          wt(cell_x+ix) = wt(cell_x+ix) + gx(ix) * l_weight
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_non_zero)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekbar



  SUBROUTINE calc_ekflux(data_array, current_species, direction)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: l_weight
    ! The data to be weighted onto the grid
    REAL(num) :: wdata = 0.0_num
    REAL(num) :: fac, gamma, ek, part_flux, xfac, yfac, zfac
    REAL(num), DIMENSION(:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(-2:nx+3))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    l_weight = 1.0_num

    xfac = c
    yfac = c * dx
    zfac = c * dx

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_mc = c * io_list(ispecies)%mass
      l_weight = io_list(ispecies)%weight
      fac = part_mc * l_weight * c

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

#ifdef PHOTONS
        IF (io_list(ispecies)%species_type .NE. c_species_id_photon) THEN
#endif
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          ek = (gamma - 1.0_num) * fac
#ifdef PHOTONS
        ELSE
          fac = c / current%particle_energy
          part_ux = current%part_p(1) * fac
          part_uy = current%part_p(2) * fac
          part_uz = current%part_p(3) * fac
          ek = current%particle_energy * l_weight
          gamma = 1.0_num
        ENDIF
#endif

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

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * wdata
          wt(cell_x+ix) = wt(cell_x+ix) + gx(ix) * l_weight
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    CALL calc_boundary(wt)

    data_array = data_array / MAX(wt, c_non_zero)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

    DEALLOCATE(wt)

  END SUBROUTINE calc_ekflux



  SUBROUTINE calc_poynt_flux(data_array, direction)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: direction
    INTEGER :: ix
    REAL(num) :: ex_cc, ey_cc, ez_cc, bx_cc, by_cc, bz_cc

    SELECT CASE(direction)
    CASE(c_dir_x)
      DO ix = 1, nx
        ey_cc = ey(ix)
        ez_cc = ez(ix)
        by_cc = 0.5_num * (by(ix-1) + by(ix))
        bz_cc = 0.5_num * (bz(ix-1) + bz(ix))
        data_array(ix) = (ey_cc * bz_cc - ez_cc * by_cc) / mu0
      ENDDO
    CASE(c_dir_y)
      DO ix = 1, nx
        ex_cc = 0.5_num * (ex(ix-1) + ex(ix))
        ez_cc = ez(ix)
        bx_cc = bx(ix)
        bz_cc = 0.5_num * (bz(ix-1) + bz(ix))
        data_array(ix) = (ez_cc * bx_cc - ex_cc * bz_cc) / mu0
      ENDDO
    CASE(c_dir_z)
      DO ix = 1, nx
        ex_cc = 0.5_num * (ex(ix-1) + ex(ix))
        ey_cc = ey(ix)
        bx_cc = bx(ix)
        by_cc = 0.5_num * (by(ix-1) + by(ix))
        data_array(ix) = (ex_cc * by_cc - ey_cc * bx_cc) / mu0
      ENDDO
    END SELECT

  END SUBROUTINE calc_poynt_flux



  SUBROUTINE calc_charge_density(data_array, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type .EQ. c_species_id_photon) CYCLE
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      part_q  = io_list(ispecies)%charge
      fac = io_list(ispecies)%weight * idx
      wdata = part_q * fac

      DO WHILE (ASSOCIATED(current))
        ! Copy the particle properties out for speed
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
#ifdef PER_PARTICLE_WEIGHT
        fac = current%weight * idx
#endif
        wdata = part_q * fac
#else
#ifdef PER_PARTICLE_WEIGHT
        fac = current%weight * idx
        wdata = part_q * fac
#endif
#endif

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * wdata
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_charge_density



  SUBROUTINE calc_number_density(data_array, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: ispecies, ix, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num

    idx   = 1.0_num / dx

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      wdata = io_list(ispecies)%weight * idx

      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_WEIGHT
        wdata = current%weight * idx
#endif

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * wdata
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_number_density



  SUBROUTINE calc_temperature(sigma, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: l_weight
    REAL(num) :: gf
    REAL(num), DIMENSION(:), ALLOCATABLE :: part_count, meanx, meany, meanz
    INTEGER :: ispecies, ix, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    ALLOCATE(meanx(-2:nx+3))
    ALLOCATE(meany(-2:nx+3))
    ALLOCATE(meanz(-2:nx+3))
    ALLOCATE(part_count(-2:nx+3))
    meanx = 0.0_num
    meany = 0.0_num
    meanz = 0.0_num
    part_count = 0.0_num
    sigma = 0.0_num

    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
      l_weight = io_list(ispecies)%weight

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

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          gf = gx(ix) * l_weight
          meanx(cell_x+ix) = meanx(cell_x+ix) + gf * part_pmx
          meany(cell_x+ix) = meany(cell_x+ix) + gf * part_pmy
          meanz(cell_x+ix) = meanz(cell_x+ix) + gf * part_pmz
          part_count(cell_x+ix) = part_count(cell_x+ix) + gf
        ENDDO
        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(meanx)
    CALL calc_boundary(meany)
    CALL calc_boundary(meanz)
    CALL calc_boundary(part_count)

    part_count = MAX(part_count, 1.e-6_num)

    meanx = meanx / part_count
    meany = meany / part_count
    meanz = meanz / part_count

    part_count = 0.0_num
    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)

      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          gf = gx(ix)
          sigma(cell_x+ix) = sigma(cell_x+ix) + gf &
              * ((part_pmx - meanx(cell_x+ix))**2 &
              + (part_pmy - meany(cell_x+ix))**2 &
              + (part_pmz - meanz(cell_x+ix))**2)
          part_count(cell_x+ix) = part_count(cell_x+ix) + gf
        ENDDO
        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! 3/2 kT = <p^2>/(2m)
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / 3.0_num

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    INTEGER :: ispecies, ix, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

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
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
      DO WHILE (ASSOCIATED(current))
#include "particle_to_grid.inc"

        wdata = evaluator(current, ispecies)
        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * wdata
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_on_grid_with_evaluator



  SUBROUTINE calc_per_species_current(data_array, current_species, direction)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction

    REAL(num) :: part_px, part_py, part_pz
    REAL(num) :: part_q, part_mc, part_weight
    REAL(num) :: part_j = 0.0_num
    REAL(num) :: idxyz, root
    INTEGER :: ispecies, spec_start, spec_end, ix
    INTEGER(i8) :: ipart
    LOGICAL :: spec_sum

    TYPE (particle), POINTER :: current, next
#include "particle_head.inc"

    data_array = 0.0_num

    idxyz = 1.0_num / dx

    spec_start = current_species
    spec_end = current_species
    spec_sum = .FALSE.

    IF (current_species .LE. 0) THEN
      spec_start = 1
      spec_end = n_species
      spec_sum = .TRUE.
    ENDIF

    DO ispecies = spec_start, spec_end
      IF (io_list(ispecies)%species_type .EQ. c_species_id_photon) CYCLE
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head

      part_weight = io_list(ispecies)%weight
      part_q  = io_list(ispecies)%charge
      part_mc = c * io_list(ispecies)%mass

      DO ipart = 1, io_list(ispecies)%attached_list%count
        next => current%next
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q  = current%charge
        part_mc = c * current%mass
#endif

        ! Copy the particle properties out for speed
        part_px = current%part_p(1)
        part_py = current%part_p(2)
        part_pz = current%part_p(3)
        root = c / SQRT(part_mc**2 + part_px**2 + part_py**2 + part_pz**2)
        SELECT CASE (direction)
          CASE(c_dir_x)
            part_j = part_q * part_px * root * part_weight * idxyz
          CASE(c_dir_y)
            part_j = part_q * part_py * root * part_weight * idxyz
          CASE(c_dir_z)
            part_j = part_q * part_pz * root * part_weight * idxyz
        END SELECT

#include "particle_to_grid.inc"

        DO ix = sf_min, sf_max
          data_array(cell_x+ix) = data_array(cell_x+ix) + gx(ix) * part_j
        ENDDO
        current => next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array, ng, direction)

  END SUBROUTINE calc_per_species_current



  SUBROUTINE calc_per_species_jx(data_array, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_x)

  END SUBROUTINE calc_per_species_jx



  SUBROUTINE calc_per_species_jy(data_array, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_y)

  END SUBROUTINE calc_per_species_jy



  SUBROUTINE calc_per_species_jz(data_array, current_species)

    REAL(num), DIMENSION(-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_z)

  END SUBROUTINE calc_per_species_jz



  SUBROUTINE calc_total_energy_sum

    REAL(num) :: particle_energy, field_energy
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_mc, l_weight, fac, gamma
    REAL(num) :: sum_out(2), sum_in(2)
    REAL(num), PARAMETER :: c2 = c**2
    INTEGER :: ispecies, i
    TYPE(particle), POINTER :: current

    particle_energy = 0.0_num

    ! Sum over all particles to calculate total kinetic energy
    DO ispecies = 1, n_species
#ifdef TRACER_PARTICLES
      IF (species_list(ispecies)%tracer) CYCLE
#endif
      current => species_list(ispecies)%attached_list%head
      part_mc = c * species_list(ispecies)%mass
      l_weight = species_list(ispecies)%weight
      fac = part_mc * l_weight * c

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

#ifdef PHOTONS
        IF (species_list(ispecies)%species_type .NE. c_species_id_photon) THEN
#endif
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          particle_energy = particle_energy + (gamma - 1.0_num) * fac
#ifdef PHOTONS
        ELSE
          particle_energy = particle_energy + current%particle_energy * l_weight
        ENDIF
#endif

        current => current%next
      ENDDO
    ENDDO

    ! EM field energy
    field_energy = 0.0_num
    DO i = 1, nx
      field_energy = field_energy + ex(i)**2 + ey(i)**2 &
          + ez(i)**2 + c2 * (bx(i)**2 + by(i)**2 + bz(i)**2)
    ENDDO
    field_energy = 0.5_num * epsilon0 * field_energy * dx

    sum_out(1) = particle_energy
    sum_out(2) = field_energy
    CALL MPI_ALLREDUCE(sum_out, sum_in, 2, mpireal, MPI_SUM, comm, errcode)
    total_particle_energy = sum_in(1)
    total_field_energy = sum_in(2)

  END SUBROUTINE calc_total_energy_sum

END MODULE calc_df
