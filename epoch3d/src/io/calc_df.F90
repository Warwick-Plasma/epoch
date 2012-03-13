MODULE calc_df

  USE boundary
  USE shape_functions

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calc_boundary(data_array)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER :: i, j, k

    CALL processor_summation_bcs(data_array)

    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic .AND. x_min_boundary) THEN
      DO k = -2, nz+3
        DO j = -2, ny+3
          data_array(1,j,k) = data_array(1,j,k) + data_array( 0,j,k)
          data_array(2,j,k) = data_array(2,j,k) + data_array(-1,j,k)
          data_array(3,j,k) = data_array(3,j,k) + data_array(-2,j,k)
        ENDDO
      ENDDO
    ENDIF
    IF (bc_particle(c_bd_x_max) .NE. c_bc_periodic .AND. x_max_boundary) THEN
      DO k = -2, nz+3
        DO j = -2, ny+3
          data_array(nx-2,j,k) = data_array(nx-2,j,k) + data_array(nx+3,j,k)
          data_array(nx-1,j,k) = data_array(nx-1,j,k) + data_array(nx+2,j,k)
          data_array(nx  ,j,k) = data_array(nx  ,j,k) + data_array(nx+1,j,k)
        ENDDO
      ENDDO
    ENDIF

    IF (bc_particle(c_bd_y_min) .NE. c_bc_periodic .AND. y_min_boundary) THEN
      DO k = -2, nz+3
        DO i = -2, nx+3
          data_array(i,1,k) = data_array(i,1,k) + data_array(i, 0,k)
          data_array(i,2,k) = data_array(i,2,k) + data_array(i,-1,k)
          data_array(i,3,k) = data_array(i,3,k) + data_array(i,-2,k)
        ENDDO
      ENDDO
    ENDIF
    IF (bc_particle(c_bd_y_max) .NE. c_bc_periodic .AND. y_max_boundary) THEN
      DO k = -2, nz+3
        DO i = -2, nx+3
          data_array(i,ny-2,k) = data_array(i,ny-2,k) + data_array(i,ny+3,k)
          data_array(i,ny-1,k) = data_array(i,ny-1,k) + data_array(i,ny+2,k)
          data_array(i,ny  ,k) = data_array(i,ny  ,k) + data_array(i,ny+1,k)
        ENDDO
      ENDDO
    ENDIF

    IF (bc_particle(c_bd_z_min) .NE. c_bc_periodic .AND. z_min_boundary) THEN
      DO j = -2, ny+3
        DO i = -2, nx+3
          data_array(i,j,1) = data_array(i,j,1) + data_array(i,j, 0)
          data_array(i,j,2) = data_array(i,j,2) + data_array(i,j,-1)
          data_array(i,j,3) = data_array(i,j,3) + data_array(i,j,-2)
        ENDDO
      ENDDO
    ENDIF
    IF (bc_particle(c_bd_z_max) .NE. c_bc_periodic .AND. z_max_boundary) THEN
      DO j = -2, ny+3
        DO i = -2, nx+3
          data_array(i,j,nz-2) = data_array(i,j,nz-2) + data_array(i,j,nz+3)
          data_array(i,j,nz-1) = data_array(i,j,nz-1) + data_array(i,j,nz+2)
          data_array(i,j,nz  ) = data_array(i,j,nz  ) + data_array(i,j,nz+1)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE calc_boundary



  SUBROUTINE calc_mass_density(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_m
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_m = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy / dz

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
#ifndef PER_PARTICLE_CHARGE_MASS
      part_m  = io_list(ispecies)%mass
#ifndef PER_PARTICLE_WEIGHT
      fac = io_list(ispecies)%weight * idx
#endif
      wdata = part_m * fac
#else
#ifndef PER_PARTICLE_WEIGHT
      fac = io_list(ispecies)%weight * idx
      wdata = part_m * fac
#endif
#endif
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

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
          ENDDO
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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: l_weight
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(-2:nx+3,-2:ny+3,-2:nz+3))
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
#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc = c * io_list(ispecies)%mass
#ifndef PER_PARTICLE_WEIGHT
      l_weight = io_list(ispecies)%weight
#endif
      fac = part_mc * l_weight * c
#else
#ifndef PER_PARTICLE_WEIGHT
      l_weight = io_list(ispecies)%weight
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
#ifdef PHOTONS
        IF (io_list(ispecies)%species_type .NE. c_species_id_photon) THEN
#endif
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
#ifdef PHOTONS
        ELSE
          part_ux = current%part_p(1) / c
          part_uy = current%part_p(2) / c
          part_uz = current%part_p(3) / c
        ENDIF
#endif

#include "particle_to_grid.inc"

#ifdef PHOTONS
        IF (io_list(ispecies)%species_type .NE. c_species_id_photon) THEN
#endif
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          wdata = (gamma - 1.0_num) * fac
#ifdef PHOTONS
        ELSE
          wdata = current%particle_energy * l_weight
        ENDIF
#endif
        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
              wt(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  wt(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * l_weight
            ENDDO
          ENDDO
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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_ux, part_uy, part_uz, part_mc
    ! The weight of a particle
    REAL(num) :: l_weight
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, gamma, ek, part_flux, xfac, yfac, zfac
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: wt
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    ALLOCATE(wt(-2:nx+3,-2:ny+3,-2:nz+3))
    data_array = 0.0_num
    wt = 0.0_num
    part_mc  = 1.0_num
    l_weight = 1.0_num

    xfac = c * dy * dz
    yfac = c * dx * dz
    zfac = c * dx * dy

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
#ifndef PER_PARTICLE_CHARGE_MASS
      part_mc = c * io_list(ispecies)%mass
#ifndef PER_PARTICLE_WEIGHT
      l_weight = io_list(ispecies)%weight
#endif
      fac = part_mc * l_weight * c
#else
#ifndef PER_PARTICLE_WEIGHT
      l_weight = io_list(ispecies)%weight
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
#ifdef PHOTONS
        IF (io_list(ispecies)%species_type .NE. c_species_id_photon) THEN
#endif
          part_ux = current%part_p(1) / part_mc
          part_uy = current%part_p(2) / part_mc
          part_uz = current%part_p(3) / part_mc
#ifdef PHOTONS
        ELSE
          part_ux = current%part_p(1) / c
          part_uy = current%part_p(2) / c
          part_uz = current%part_p(3) / c
        ENDIF
#endif

#include "particle_to_grid.inc"

#ifdef PHOTONS
        IF (io_list(ispecies)%species_type .NE. c_species_id_photon) THEN
#endif
          gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
          ek = (gamma - 1.0_num) * fac
#ifdef PHOTONS
        ELSE
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

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
              wt(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  wt(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * l_weight
            ENDDO
          ENDDO
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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: direction
    INTEGER :: ix, iy, iz
    REAL(num) :: ex_cc, ey_cc, ez_cc, bx_cc, by_cc, bz_cc

    SELECT CASE(direction)
    CASE(c_dir_x)
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ey_cc = 0.5_num  * (ey(ix  , iy-1, iz  ) + ey(ix, iy, iz))
            ez_cc = 0.5_num  * (ez(ix  , iy  , iz-1) + ez(ix, iy, iz))
            by_cc = 0.25_num * (by(ix-1, iy  , iz-1) + by(ix  , iy  , iz-1) &
                             +  by(ix-1, iy  , iz  ) + by(ix  , iy  , iz  ))
            bz_cc = 0.25_num * (bz(ix-1, iy-1, iz  ) + bz(ix  , iy-1, iz  ) &
                             +  bz(ix-1, iy  , iz  ) + bz(ix  , iy  , iz  ))
            data_array(ix,iy,iz) = (ey_cc * bz_cc - ez_cc * by_cc) / mu0
          ENDDO
        ENDDO
      ENDDO
    CASE(c_dir_y)
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ex_cc = 0.5_num  * (ex(ix-1, iy  , iz  ) + ex(ix, iy, iz))
            ez_cc = 0.5_num  * (ez(ix  , iy  , iz-1) + ez(ix, iy, iz))
            bx_cc = 0.25_num * (bx(ix  , iy-1, iz-1) + bx(ix  , iy  , iz-1) &
                             +  bx(ix  , iy-1, iz  ) + bx(ix  , iy  , iz  ))
            bz_cc = 0.25_num * (bz(ix-1, iy-1, iz  ) + bz(ix  , iy-1, iz  ) &
                             +  bz(ix-1, iy  , iz  ) + bz(ix  , iy  , iz  ))
            data_array(ix,iy,iz) = (ez_cc * bx_cc - ex_cc * bz_cc) / mu0
          ENDDO
        ENDDO
      ENDDO
    CASE(c_dir_z)
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ex_cc = 0.5_num  * (ex(ix-1, iy  , iz  ) + ex(ix, iy, iz))
            ey_cc = 0.5_num  * (ey(ix  , iy-1, iz  ) + ey(ix, iy, iz))
            bx_cc = 0.25_num * (bx(ix  , iy-1, iz-1) + bx(ix  , iy  , iz-1) &
                             +  bx(ix  , iy-1, iz  ) + bx(ix  , iy  , iz  ))
            by_cc = 0.25_num * (by(ix-1, iy  , iz-1) + by(ix  , iy  , iz-1) &
                             +  by(ix-1, iy  , iz  ) + by(ix  , iy  , iz  ))
            data_array(ix,iy,iz) = (ex_cc * by_cc - ey_cc * bx_cc) / mu0
          ENDDO
        ENDDO
      ENDDO
    END SELECT

  END SUBROUTINE calc_poynt_flux



  SUBROUTINE calc_charge_density(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_q
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: fac, idx
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num
    part_q = 0.0_num
    fac = 0.0_num

    idx = 1.0_num / dx / dy / dz

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
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q  = io_list(ispecies)%charge
#ifndef PER_PARTICLE_WEIGHT
      fac = io_list(ispecies)%weight * idx
#endif
      wdata = part_q * fac
#else
#ifndef PER_PARTICLE_WEIGHT
      fac = io_list(ispecies)%weight * idx
      wdata = part_q * fac
#endif
#endif
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

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
          ENDDO
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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    REAL(num) :: idx
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
    TYPE(particle), POINTER :: current
    LOGICAL :: spec_sum
#include "particle_head.inc"

    data_array = 0.0_num

    idx   = 1.0_num / dx / dy / dz

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
#ifndef PER_PARTICLE_WEIGHT
      wdata = io_list(ispecies)%weight * idx
#endif
      DO WHILE (ASSOCIATED(current))
#ifdef PER_PARTICLE_WEIGHT
        wdata = current%weight * idx
#endif

#include "particle_to_grid.inc"

        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
          ENDDO
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

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: sigma
    INTEGER, INTENT(IN) :: current_species
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_pmx, part_pmy, part_pmz, sqrt_part_m
    ! The weight of a particle
    REAL(num) :: l_weight
    REAL(num) :: gf
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: part_count, meanx, meany, meanz
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
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
#ifdef TRACER_PARTICLES
      IF (spec_sum .AND. io_list(ispecies)%tracer) CYCLE
#endif
      current => io_list(ispecies)%attached_list%head
#ifndef PER_PARTICLE_CHARGE_MASS
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
#endif
#ifndef PER_PARTICLE_WEIGHT
      l_weight = io_list(ispecies)%weight
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

#include "particle_to_grid.inc"

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
#ifndef PER_PARTICLE_CHARGE_MASS
      sqrt_part_m  = SQRT(io_list(ispecies)%mass)
#endif
      DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
        sqrt_part_m  = SQRT(current%mass)
#endif
        ! Copy the particle properties out for speed
        part_pmx = current%part_p(1) / sqrt_part_m
        part_pmy = current%part_p(2) / sqrt_part_m
        part_pmz = current%part_p(3) / sqrt_part_m

#include "particle_to_grid.inc"

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
        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(sigma)
    CALL calc_boundary(part_count)

    ! N/2 kT = <p^2>/(2m), where N is the number of degrees of freedom
    sigma = sigma / MAX(part_count, 1.e-6_num) / kb / REAL(dof)

    DEALLOCATE(part_count, meanx, meany, meanz)

  END SUBROUTINE calc_temperature



  SUBROUTINE calc_on_grid_with_evaluator(data_array, current_species, evaluator)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species
    ! The data to be weighted onto the grid
    REAL(num) :: wdata
    INTEGER :: ispecies, ix, iy, iz, spec_start, spec_end
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
        DO iz = sf_min, sf_max
          DO iy = sf_min, sf_max
            DO ix = sf_min, sf_max
              data_array(cell_x+ix, cell_y+iy, cell_z+iz) = &
                  data_array(cell_x+ix, cell_y+iy, cell_z+iz) &
                  + gx(ix) * gy(iy) * gz(iz) * wdata
            ENDDO
          ENDDO
        ENDDO

        current => current%next
      ENDDO
    ENDDO

    CALL calc_boundary(data_array)
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(data_array, c_stagger_centre, ix)
    ENDDO

  END SUBROUTINE calc_on_grid_with_evaluator



  ! This subroutine calculates 'per species' currents, in the same way as in the
  ! particle push, but without looping over species index. Hot electron current
  ! and ion current are needed for CKD scheme.

  SUBROUTINE calc_per_species_current(data_array, current_species, direction)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species, direction

    INTEGER :: cell_x1, cell_x3
    INTEGER :: cell_y1, cell_y3
    INTEGER :: cell_z1, cell_z3
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
    REAL(num), DIMENSION(sf0-2:sf1+1,sf0-1:sf1+1,sf0-1:sf1+1) :: jxh
    REAL(num), DIMENSION(sf0-1:sf1+1,sf0-2:sf1+1,sf0-1:sf1+1) :: jyh
    REAL(num), DIMENSION(sf0-1:sf1+1,sf0-1:sf1+1,sf0-2:sf1+1) :: jzh
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, ipart_mc, part_weight
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy, gz
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy, hz
    INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL(num) :: wx, wy, wz
    REAL(num) :: idx, idy, idz
    REAL(num) :: idtyz, idtxz, idtxy
    REAL(num) :: idt, dtc
    REAL(num) :: fcx, fcy, fcz, fjx, fjy, fjz
    REAL(num) :: root, fac, third, gamma, cf2
    REAL(num) :: delta_x, delta_y, delta_z
    INTEGER :: ispecies, ix, iy, iz, dcellx, dcelly, dcellz
    INTEGER(KIND=8) :: ipart
    TYPE(particle), POINTER :: current, next
    INTEGER :: spec_start, spec_end
    LOGICAL :: spec_sum

    data_array = 0.0_num
    gx = 0.0_num
    gy = 0.0_num
    gz = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idz = 1.0_num / dz
    idt = 1.0_num / dt
    dtc = c * dt
    third = 1.0_num / 3.0_num
    ! particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    fac = 1.0_num / 24.0_num
#elif  PARTICLE_SHAPE_TOPHAT
    fac = 1.0_num
#else
    fac = 0.5_num
#endif

    idtyz = idt * idy * idz * fac**3
    idtxz = idt * idx * idz * fac**3
    idtxy = idt * idx * idy * fac**3

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

#ifndef PER_PARTICLE_WEIGHT
      part_weight = io_list(ispecies)%weight
      fcx = idtyz * part_weight
      fcy = idtxz * part_weight
      fcz = idtxy * part_weight
#endif
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q   = io_list(ispecies)%charge
      ipart_mc = 1.0_num / c / io_list(ispecies)%mass
#endif
      !DEC$ VECTOR ALWAYS
      DO ipart = 1, io_list(ispecies)%attached_list%count
        next => current%next
#ifdef PER_PARTICLE_WEIGHT
        part_weight = current%weight
        fcx = idtyz * part_weight
        fcy = idtxz * part_weight
        fcz = idtxy * part_weight
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q   = current%charge
        ipart_mc = 1.0_num / c / current%mass
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_min_local
        part_y  = current%part_pos(2) - y_min_local
        part_z  = current%part_pos(3) - z_min_local
        part_ux = current%part_p(1) * ipart_mc
        part_uy = current%part_p(2) * ipart_mc
        part_uz = current%part_p(3) * ipart_mc

        ! Work out the grid cell number for the particle.
        ! Not an integer in general.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
        cell_z_r = part_z * idz - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
        cell_z_r = part_z * idz
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        cell_z1 = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z1, num) - cell_z_r
        cell_z1 = cell_z1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
#ifdef PARTICLE_SHAPE_BSPLINE3
        INCLUDE '../include/bspline3/gx.inc'
#elif  PARTICLE_SHAPE_TOPHAT
        INCLUDE '../include/tophat/gx.inc'
#else
        INCLUDE '../include/triangle/gx.inc'
#endif

        ! Calculate particle velocity from particle momentum
        gamma = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        root = dtc / gamma

        delta_x = part_ux * root
        delta_y = part_uy * root
        delta_z = part_uz * root

        ! Move particles to end of time step
        part_x = part_x + delta_x
        part_y = part_y + delta_y
        part_z = part_z + delta_z

#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
        cell_z_r = part_z * idz - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
        cell_z_r = part_z * idz
#endif
        cell_x3 = FLOOR(cell_x_r + 0.5_num)
        cell_frac_x = REAL(cell_x3, num) - cell_x_r
        cell_x3 = cell_x3 + 1

        cell_y3 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y3, num) - cell_y_r
        cell_y3 = cell_y3 + 1

        cell_z3 = FLOOR(cell_z_r + 0.5_num)
        cell_frac_z = REAL(cell_z3, num) - cell_z_r
        cell_z3 = cell_z3 + 1

        hx = 0.0_num
        hy = 0.0_num
        hz = 0.0_num

        dcellx = cell_x3 - cell_x1
        dcelly = cell_y3 - cell_y1
        dcellz = cell_z3 - cell_z1
#ifdef PARTICLE_SHAPE_BSPLINE3
        INCLUDE '../include/bspline3/hx_dcell.inc'
#elif  PARTICLE_SHAPE_TOPHAT
        INCLUDE '../include/tophat/hx_dcell.inc'
#else
        INCLUDE '../include/triangle/hx_dcell.inc'
#endif

        ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
        ! the current update much simpler
        hx = hx - gx
        hy = hy - gy
        hz = hz - gz

        ! Remember that due to CFL condition particle can never cross more
        ! than one gridcell in one timestep

        xmin = sf_min + (dcellx - 1) / 2
        xmax = sf_max + (dcellx + 1) / 2

        ymin = sf_min + (dcelly - 1) / 2
        ymax = sf_max + (dcelly + 1) / 2

        zmin = sf_min + (dcellz - 1) / 2
        zmax = sf_max + (dcellz + 1) / 2

        ! This is the bit that actually solves d(rho)/dt = -div(J)
        SELECT CASE (direction)
          CASE(c_dir_x)
            ! Set this to zero due to diffential inside loop
            jxh = 0.0_num
            fjx = fcx * part_q
            DO iz = zmin, zmax
              DO iy = ymin, ymax
                DO ix = xmin, xmax
                  wx =  hx(ix) * (gy(iy) * (gz(iz) + 0.5_num * hz(iz)) &
                      + hy(iy) * (third  *  hz(iz) + 0.5_num * gz(iz)))
                  jxh(ix, iy, iz) = jxh(ix-1, iy, iz) - fjx * wx
                  data_array(cell_x1+ix, cell_y1+iy, cell_z1+iz) = &
                      data_array(cell_x1+ix, cell_y1+iy, cell_z1+iz) &
                          + jxh(ix, iy, iz)
                ENDDO
              ENDDO
            ENDDO
          CASE(c_dir_y)
            ! Set this to zero due to diffential inside loop
            jyh = 0.0_num
            fjy = fcy * part_q
            DO iz = zmin, zmax
              DO iy = ymin, ymax
                DO ix = xmin, xmax
                  wy =  hy(iy) * (gx(ix) * (gz(iz) + 0.5_num * hz(iz)) &
                      + hx(ix) * (third  *  hz(iz) + 0.5_num * gz(iz)))
                  jyh(ix, iy, iz) = jyh(ix, iy-1, iz) - fjy * wy
                  data_array(cell_x1+ix, cell_y1+iy, cell_z1+iz) = &
                      data_array(cell_x1+ix, cell_y1+iy, cell_z1+iz) &
                          + jyh(ix, iy, iz)
                ENDDO
              ENDDO
            ENDDO
          CASE(c_dir_z)
            ! Set this to zero due to diffential inside loop
            jzh = 0.0_num
            fjz = fcz * part_q
            DO iz = zmin, zmax
              DO iy = ymin, ymax
                DO ix = xmin, xmax
                  wz =  hz(iz) * (gx(ix) * (gy(iy) + 0.5_num * hy(iy)) &
                      + hx(ix) * (third  *  hy(iy) + 0.5_num * gy(iy)))
                  jzh(ix, iy, iz) = jzh(ix, iy, iz-1) - fjz * wz
                  data_array(cell_x1+ix, cell_y1+iy, cell_z1+iz) = &
                      data_array(cell_x1+ix, cell_y1+iy, cell_z1+iz) &
                          + jzh(ix, iy, iz)
                ENDDO
              ENDDO
            ENDDO
        END SELECT
        current => next
      ENDDO
    ENDDO

    CALL processor_summation_bcs(data_array, direction)

  END SUBROUTINE calc_per_species_current



  SUBROUTINE calc_per_species_jx(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_x)

  END SUBROUTINE calc_per_species_jx



  SUBROUTINE calc_per_species_jy(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_y)

  END SUBROUTINE calc_per_species_jy



  SUBROUTINE calc_per_species_jz(data_array, current_species)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
    INTEGER, INTENT(IN) :: current_species

    CALL calc_per_species_current(data_array, current_species, c_dir_z)

  END SUBROUTINE calc_per_species_jz

END MODULE calc_df
