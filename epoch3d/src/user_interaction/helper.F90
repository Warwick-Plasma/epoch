! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE helper

  USE boundary
  USE strings
  USE partlist
  USE calc_df

  IMPLICIT NONE

CONTAINS

  SUBROUTINE set_thermal_bcs

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species

    DO ispecies = 1, n_species
      species => species_list(ispecies)

      ! Set temperature at boundary for thermal bcs.

      IF (bc_particle(c_bd_x_min) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_min(:,:,:) = &
            species_list(ispecies)%initial_conditions%temp(1,:,:,:)
      ENDIF
      IF (bc_particle(c_bd_x_max) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_max(:,:,:) = &
            species_list(ispecies)%initial_conditions%temp(nx,:,:,:)
      ENDIF
      IF (bc_particle(c_bd_y_min) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_y_min(:,:,:) = &
            species_list(ispecies)%initial_conditions%temp(:,1,:,:)
      ENDIF
      IF (bc_particle(c_bd_y_max) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_y_max(:,:,:) = &
            species_list(ispecies)%initial_conditions%temp(:,ny,:,:)
      ENDIF
      IF (bc_particle(c_bd_z_min) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_z_min(:,:,:) = &
            species_list(ispecies)%initial_conditions%temp(:,:,1,:)
      ENDIF
      IF (bc_particle(c_bd_z_max) == c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_z_max(:,:,:) = &
            species_list(ispecies)%initial_conditions%temp(:,:,nz,:)
      ENDIF
    ENDDO

  END SUBROUTINE set_thermal_bcs



  SUBROUTINE auto_load

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species

    CALL set_thermal_bcs

    DO ispecies = 1, n_species
      species => species_list(ispecies)

#ifndef PER_SPECIES_WEIGHT
      CALL setup_particle_density(&
          species_list(ispecies)%initial_conditions%density, species, &
          species_list(ispecies)%initial_conditions%density_min, &
          species_list(ispecies)%initial_conditions%density_max)
#else
      CALL non_uniform_load_particles(&
          species_list(ispecies)%initial_conditions%density, species, &
          species_list(ispecies)%initial_conditions%density_min, &
          species_list(ispecies)%initial_conditions%density_max)
#endif
      CALL setup_particle_temperature(&
          species_list(ispecies)%initial_conditions%temp(:,:,:,1), c_dir_x, &
          species, species_list(ispecies)%initial_conditions%drift(:,:,:,1))
      CALL setup_particle_temperature(&
          species_list(ispecies)%initial_conditions%temp(:,:,:,2), c_dir_y, &
          species, species_list(ispecies)%initial_conditions%drift(:,:,:,2))
      CALL setup_particle_temperature(&
          species_list(ispecies)%initial_conditions%temp(:,:,:,3), c_dir_z, &
          species, species_list(ispecies)%initial_conditions%drift(:,:,:,3))
    ENDDO

    IF (rank == 0) THEN
      DO ispecies = 1, n_species
        species => species_list(ispecies)
        IF (species%count < 0) THEN
          WRITE(*,*) 'No particles specified for species ', &
              '"' // TRIM(species%name) // '"'
#ifndef NO_IO
          WRITE(stat_unit,*) 'No particles specified for species ', &
              '"' // TRIM(species%name) // '"'
#endif
          species%count = 0
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE auto_load



  SUBROUTINE allocate_ic

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      ALLOCATE(species_list(ispecies)%initial_conditions&
          %density(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(species_list(ispecies)%initial_conditions&
          %temp(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:3))
      ALLOCATE(species_list(ispecies)%initial_conditions&
          %drift(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,1:3))

      species_list(ispecies)%initial_conditions%density = 1.0_num
      species_list(ispecies)%initial_conditions%temp = 0.0_num
      species_list(ispecies)%initial_conditions%drift = 0.0_num
      species_list(ispecies)%initial_conditions%density_min = EPSILON(1.0_num)
      species_list(ispecies)%initial_conditions%density_max = HUGE(1.0_num)
    ENDDO

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      DEALLOCATE(species_list(ispecies)%initial_conditions%density)
      DEALLOCATE(species_list(ispecies)%initial_conditions%temp)
      DEALLOCATE(species_list(ispecies)%initial_conditions%drift)
    ENDDO

  END SUBROUTINE deallocate_ic



#ifdef PER_SPECIES_WEIGHT
  SUBROUTINE non_uniform_load_particles(density, species, density_min, &
      density_max)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: density
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(INOUT) :: density_min, density_max
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_per_cell
    REAL(num) :: density_total, density_total_global, density_average
    REAL(num) :: npart_per_cell_average
    INTEGER(i8) :: npart_this_proc_new, ipart, npart_this_species
    INTEGER :: ix, iy, iz
    CHARACTER(LEN=15) :: string
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next

    partlist => species%attached_list

    num_valid_cells_local = 0
    density_total = 0.0_num

    DO iz = 1-ng, nz+ng
    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng
      IF (density(ix,iy,iz) > density_max) density(ix,iy,iz) = density_max
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      IF (density(ix,iy,iz) >= density_min) THEN
        num_valid_cells_local = num_valid_cells_local + 1
        density_total = density_total + density(ix,iy,iz)
      ENDIF
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    CALL MPI_ALLREDUCE(num_valid_cells_local, num_valid_cells_global, 1, &
        MPI_INTEGER8, MPI_SUM, comm, errcode)

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell_average = FLOOR(species%npart_per_cell, num)
    ELSE
      npart_per_cell_average = REAL(species%count, num) &
          / REAL(num_valid_cells_global, num)
    ENDIF

    IF (npart_per_cell_average <= 0) RETURN

    CALL MPI_ALLREDUCE(density_total, density_total_global, 1, mpireal, &
        MPI_SUM, comm, errcode)
    density_average = density_total_global / REAL(num_valid_cells_global, num)

    npart_this_proc_new = 0
    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix, iy, iz) / density_average &
          * npart_per_cell_average)
      npart_this_proc_new = npart_this_proc_new + npart_per_cell
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)

    ! Randomly place npart_per_cell particles into each valid cell
    current => partlist%head
    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix, iy, iz) / density_average &
          * npart_per_cell_average)

      ipart = 0
      DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
#ifdef PER_PARTICLE_CHARGE_MASS
        ! Even if particles have per particle charge and mass, assume
        ! that initially they all have the same charge and mass (user
        ! can easily over_ride)
        current%charge = species%charge
        current%mass = species%mass
#endif
        current%part_pos(1) = x(ix) + (random() - 0.5_num) * dx
        current%part_pos(2) = y(iy) + (random() - 0.5_num) * dy
        current%part_pos(3) = z(iz) + (random() - 0.5_num) * dz

        ipart = ipart + 1
        current => current%next
      ENDDO
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current => next
    ENDDO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species
    species%weight = density_total_global * dx * dy * dz / npart_this_species

    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#ifndef NO_IO
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#endif
    ENDIF

    CALL particle_bcs

  END SUBROUTINE non_uniform_load_particles
#endif



  ! This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE load_particles(species, load_list)

    TYPE(particle_species), POINTER :: species
    LOGICAL, DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: load_list
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: valid_cell_list
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: ipart, npart_per_cell, num_int, num_total, idx
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_this_species, num_new_particles, npart_left
    INTEGER(i8), ALLOCATABLE :: num_valid_cells_all(:), num_idx(:)
    REAL(num) :: valid_cell_frac, num_real, f0, f1
    REAL(num), ALLOCATABLE :: num_frac(:)
    INTEGER(i8) :: cell_x
    INTEGER(i8) :: cell_y
    INTEGER(i8) :: cell_z
    INTEGER(i8) :: i, ipos
    INTEGER :: ix, iy, iz
    CHARACTER(LEN=15) :: string
    LOGICAL :: sweep

    npart_this_species = species%count
    IF (npart_this_species <= 0) RETURN

    num_valid_cells_local = 0
    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      IF (load_list(ix, iy, iz)) &
          num_valid_cells_local = num_valid_cells_local + 1
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
      num_new_particles = &
          FLOOR(species%npart_per_cell * num_valid_cells_local, KIND=i8)
    ELSE
      ALLOCATE(num_valid_cells_all(nproc), num_idx(nproc), num_frac(nproc))

      ! Calculate global number of particles per cell
      CALL MPI_ALLGATHER(num_valid_cells_local, 1, MPI_INTEGER8, &
          num_valid_cells_all, 1, MPI_INTEGER8, comm, errcode)

      num_valid_cells_global = 0
      DO i = 1,nproc
        num_valid_cells_global = num_valid_cells_global + num_valid_cells_all(i)
      ENDDO

      IF (num_valid_cells_global == 0) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Intial condition settings mean that there are no cells ' &
              // 'where particles may'
          WRITE(*,*) 'validly be placed for species "' // TRIM(species%name) &
              // '". ', 'Code will now terminate.'
          CALL abort_code(c_err_bad_setup)
        ENDIF
      ENDIF

      valid_cell_frac = REAL(num_valid_cells_local, num) &
          / REAL(num_valid_cells_global, num)
      num_real = npart_this_species * valid_cell_frac
      num_new_particles = FLOOR(num_real, KIND=i8)

      ! Work out which processors get the remaining fractional numbers
      ! of particles

      ! Get a list of the fractional part on each processor, along with
      ! the total
      num_total = 0
      DO i = 1,nproc
        valid_cell_frac = REAL(num_valid_cells_all(i), num) &
            / REAL(num_valid_cells_global, num)
        num_real = npart_this_species * valid_cell_frac
        num_int = FLOOR(num_real, KIND=i8)
        num_frac(i) = num_real - num_int
        num_idx (i) = i - 1
        num_total = num_total + num_int
      ENDDO
      num_total = npart_this_species - num_total

      IF (num_total > 0) THEN
        ! Sort the list of fractions into decreasing order using bubble sort
        sweep = .TRUE.
        DO WHILE(sweep)
          sweep = .FALSE.
          f0 = num_frac(1)
          DO i = 2,nproc
            f1 = num_frac(i)
            IF (f1 > f0) THEN
              num_frac(i-1) = f1
              num_frac(i) = f0
              f1 = f0
              idx = num_idx(i-1)
              num_idx(i-1) = num_idx(i)
              num_idx(i) = idx
              sweep = .TRUE.
            ENDIF
            f0 = f1
          ENDDO
        ENDDO

        ! Accumulate fractional particles until they have all been accounted
        ! for. If any of them have been assigned to the current processor,
        ! add them and exit the loop.

        DO i = 1,nproc
          IF (num_idx(i) == rank) THEN
            num_new_particles = num_new_particles + 1
            EXIT
          ENDIF
          num_total = num_total - 1
          IF (num_total <= 0) EXIT
        ENDDO
      ENDIF

      DEALLOCATE(num_valid_cells_all, num_idx, num_frac)

      species%npart_per_cell = &
          REAL(npart_this_species,num) / num_valid_cells_global
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
    ENDIF

    partlist => species%attached_list

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    ! Randomly place npart_per_cell particles into each valid cell
    npart_left = num_new_particles
    current => partlist%head
    IF (npart_per_cell > 0) THEN

      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        IF (.NOT. load_list(ix, iy, iz)) CYCLE

        ipart = 0
        DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
#ifdef PER_PARTICLE_CHARGE_MASS
          ! Even if particles have per particle charge and mass, assume
          ! that initially they all have the same charge and mass (user
          ! can easily over_ride)
          current%charge = species%charge
          current%mass = species%mass
#endif
          current%part_pos(1) = x(ix) + (random() - 0.5_num) * dx
          current%part_pos(2) = y(iy) + (random() - 0.5_num) * dy
          current%part_pos(3) = z(iz) + (random() - 0.5_num) * dz

          ipart = ipart + 1
          current => current%next

          ! One particle sucessfully placed
          npart_left = npart_left - 1
        ENDDO
      ENDDO ! ix
      ENDDO ! iy
      ENDDO ! iz

    ENDIF

    ! When num_new_particles does not equal
    ! npart_per_cell * num_valid_cells_local there will be particles left
    ! over that didn't get placed.
    ! The following loop randomly place remaining particles into valid cells.
    IF (npart_left > 0) THEN
      ALLOCATE(valid_cell_list(num_valid_cells_local))

      ipos = 0
      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        IF (load_list(ix,iy,iz)) THEN
          ipos = ipos + 1
          valid_cell_list(ipos) = ix - 1 + nx * (iy - 1 + ny * (iz - 1))
        ENDIF
      ENDDO ! ix
      ENDDO ! iy
      ENDDO ! iz

      DO i = 1, npart_left
        ipos = INT(random() * (num_valid_cells_local - 1)) + 1
        ipos = valid_cell_list(ipos)

        cell_z = ipos / (nx * ny)
        ipos = ipos - (nx * ny) * cell_z
        cell_z = cell_z + 1

        cell_y = ipos / nx
        ipos = ipos - nx * cell_y
        cell_y = cell_y + 1

        cell_x = ipos + 1

        current%part_pos(1) = x(cell_x) + (random() - 0.5_num) * dx
        current%part_pos(2) = y(cell_y) + (random() - 0.5_num) * dy
        current%part_pos(3) = z(cell_z) + (random() - 0.5_num) * dz

        current => current%next
      ENDDO

      DEALLOCATE(valid_cell_list)
    ENDIF

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current => next
    ENDDO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species

    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#ifndef NO_IO
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
#endif
    ENDIF

    CALL particle_bcs

  END SUBROUTINE load_particles



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE setup_particle_density(density_in, species, density_min, &
      density_max)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(IN) :: density_in
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(IN) :: density_min, density_max
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: npart_in_cell
    REAL(num) :: wdata
    TYPE(particle_list), POINTER :: partlist
    INTEGER :: ix, iy, iz, i, j, k, isubx, isuby, isubz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: density
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: density_map
#include "particle_head.inc"

    ALLOCATE(density(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(density_map(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    density = density_in
    density_map = .FALSE.

    CALL field_bc(density, ng)

    DO iz = 1-ng, nz+ng
    DO iy = 1-ng, ny+ng
    DO ix = 1-ng, nx+ng
      IF (density(ix,iy,iz) > density_max) density(ix,iy,iz) = density_max
      IF (density(ix,iy,iz) >= density_min) THEN
        density_map(ix,iy,iz) = .TRUE.
      ELSE
        density(ix,iy,iz) = 0.0_num
      ENDIF
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    ! Uniformly load particles in space
    CALL load_particles(species, density_map)

    ALLOCATE(npart_in_cell(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    npart_in_cell = 0

    partlist => species%attached_list
    ! If using per particle weighing then use the weight function to match the
    ! uniform pseudoparticle density to the real particle density
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
      IF (.NOT. ASSOCIATED(current)) PRINT *, 'Bad Particle'

#include "particle_to_grid.inc"

      ! Calculate density at the particle position
      wdata = 0.0_num
      DO isubz = sf_min, sf_max
        i = cell_x
        j = cell_y
        k = cell_z + isubz
#ifdef PARTICLE_SHAPE_TOPHAT
        IF (.NOT. density_map(i,j,k)) k = cell_z + 1 - isubz
#else
        IF (.NOT. density_map(i,j,k)) THEN
          k = cell_z + isubz / 2
#ifdef PARTICLE_SHAPE_BSPLINE3
          IF (.NOT. density_map(i,j,k)) k = cell_z - isubz / 2
#endif
        ENDIF
#endif
        DO isuby = sf_min, sf_max
          i = cell_x
          j = cell_y + isuby
#ifdef PARTICLE_SHAPE_TOPHAT
          IF (.NOT. density_map(i,j,k)) j = cell_y + 1 - isuby
#else
          IF (.NOT. density_map(i,j,k)) THEN
            j = cell_y + isuby / 2
#ifdef PARTICLE_SHAPE_BSPLINE3
            IF (.NOT. density_map(i,j,k)) j = cell_y - isuby / 2
#endif
          ENDIF
#endif
          DO isubx = sf_min, sf_max
            i = cell_x + isubx
#ifdef PARTICLE_SHAPE_TOPHAT
            IF (.NOT. density_map(i,j,k)) i = cell_x + 1 - isubx
#else
            IF (.NOT. density_map(i,j,k)) THEN
              i = cell_x + isubx / 2
#ifdef PARTICLE_SHAPE_BSPLINE3
              IF (.NOT. density_map(i,j,k)) i = cell_x - isubx / 2
#endif
            ENDIF
#endif
            wdata = wdata + gx(isubx) * gy(isuby) * gz(isubz) * density(i,j,k)
          ENDDO ! isubx
        ENDDO ! isuby
      ENDDO ! isubz

      current%weight = wdata
      npart_in_cell(cell_x,cell_y,cell_z) = &
          npart_in_cell(cell_x,cell_y,cell_z) + 1

      current => current%next
      ipart = ipart + 1
    ENDDO
    DEALLOCATE(density_map)
    DEALLOCATE(density)

    wdata = dx * dy * dz

    partlist => species%attached_list
    ! Second loop renormalises particle weights
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PARTICLE_SHAPE_TOPHAT
      cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx) + 1
      cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy) + 1
      cell_z = FLOOR((current%part_pos(3) - z_grid_min_local) / dz) + 1
#else
      cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
      cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
      cell_z = FLOOR((current%part_pos(3) - z_grid_min_local) / dz + 1.5_num)
#endif

      current%weight = current%weight * wdata &
          / npart_in_cell(cell_x,cell_y,cell_z)

      current => current%next
      ipart = ipart + 1
    ENDDO

    DEALLOCATE(npart_in_cell)

  END SUBROUTINE setup_particle_density
#endif



  FUNCTION sample_dist_function(axis, dist_fn)

    REAL(num), DIMENSION(:), INTENT(IN) :: axis, dist_fn
    REAL(num), DIMENSION(:), ALLOCATABLE :: cdf
    REAL(num) :: position, d_cdf
    INTEGER :: n_points, ipoint, start, endpoint, current
    REAL(num) :: sample_dist_function

    n_points = SIZE(dist_fn)
    ALLOCATE(cdf(n_points))

    cdf(1) = dist_fn(1)
    DO ipoint = 2, n_points
      cdf(ipoint) = cdf(ipoint-1) + dist_fn(ipoint)
    ENDDO

    cdf = cdf / cdf(n_points)

    position = random()
    sample_dist_function = 0.0_num

    start = 1
    endpoint = n_points
    current = (start + endpoint) / 2

    DO current = 1, n_points-1
      IF (cdf(current) <= position .AND. cdf(current+1) >= position) THEN
        d_cdf = cdf(current+1) - cdf(current)
        sample_dist_function = (axis(current) * (position - cdf(current)) &
            + axis(current+1) * (cdf(current+1) - position)) / d_cdf
        EXIT
      ENDIF
    ENDDO

    DEALLOCATE(cdf)

  END FUNCTION sample_dist_function

END MODULE helper
