MODULE helper

  USE mpi
  USE boundary
  USE random_generator
  USE strings

  IMPLICIT NONE

CONTAINS

  SUBROUTINE auto_load

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species

    DO ispecies = 1, n_species
      species=>species_list(ispecies)

      ! Set temperature at boundary for thermal bcs.

      IF (bc_particle(c_bd_x_min) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_min(-2:ny+3,-2:nz+3,1:3) = &
            initial_conditions(ispecies)%temp(1,-2:ny+3,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_max(-2:ny+3,-2:nz+3,1:3) = &
            initial_conditions(ispecies)%temp(nx,-2:ny+3,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_y_min) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_y_min(-2:nx+3,-2:nz+3,1:3) = &
            initial_conditions(ispecies)%temp(-2:nx+3,1,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_y_max) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_y_max(-2:nx+3,-2:nz+3,1:3) = &
            initial_conditions(ispecies)%temp(-2:nx+3,ny,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_z_min) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_z_min(-2:nx+3,-2:ny+3,1:3) = &
            initial_conditions(ispecies)%temp(-2:nx+3,-2:ny+3,1,1:3)
      ENDIF
      IF (bc_particle(c_bd_z_max) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_z_max(-2:nx+3,-2:ny+3,1:3) = &
            initial_conditions(ispecies)%temp(-2:nx+3,-2:ny+3,nz,1:3)
      ENDIF

#ifdef PER_PARTICLE_WEIGHT
      CALL setup_particle_density(initial_conditions(ispecies)%density, &
          species, initial_conditions(ispecies)%density_min, &
          initial_conditions(ispecies)%density_max)
#else
      CALL non_uniform_load_particles(initial_conditions(ispecies)%density, &
          species, initial_conditions(ispecies)%density_min, &
          initial_conditions(ispecies)%density_max)
#endif
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,:,:,1), c_dir_x, species, &
          initial_conditions(ispecies)%drift(:,:,:,1))
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,:,:,2), c_dir_y, species, &
          initial_conditions(ispecies)%drift(:,:,:,2))
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,:,:,3), c_dir_z, species, &
          initial_conditions(ispecies)%drift(:,:,:,3))
    ENDDO

  END SUBROUTINE auto_load



  SUBROUTINE allocate_ic

    INTEGER :: ispecies

    ALLOCATE(initial_conditions(1:n_species))
    DO ispecies = 1, n_species
      ALLOCATE(initial_conditions(ispecies)%density(-2:nx+3,-2:ny+3,-2:nz+3))
      ALLOCATE(initial_conditions(ispecies)%temp (-2:nx+3,-2:ny+3,-2:nz+3,1:3))
      ALLOCATE(initial_conditions(ispecies)%drift(-2:nx+3,-2:ny+3,-2:nz+3,1:3))

      initial_conditions(ispecies)%density = 1.0_num
      initial_conditions(ispecies)%temp = 0.0_num
      initial_conditions(ispecies)%drift = 0.0_num
      initial_conditions(ispecies)%density_min = EPSILON(1.0_num)
      initial_conditions(ispecies)%density_max = HUGE(1.0_num)
    ENDDO

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies, ix, iy, iz
    REAL(num) :: min_dt, omega, k_max

    min_dt = 1000000.0_num
    k_max = 2.0_num * pi / MIN(dx, dy)

    ! Identify the plasma frequency (Bohm-Gross dispersion relation)
    ! Note that this doesn't get strongly relativistic plasmas right
    DO ispecies = 1, n_species
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            omega = &
                SQRT((initial_conditions(ispecies)%density(ix, iy, iz)*q0**2) &
                / (species_list(ispecies)%mass * epsilon0) &
                + 6.0_num * k_max**2 * kb &
                * MAXVAL(initial_conditions(ispecies)%temp(ix, iy, iz,:)) &
                / (species_list(ispecies)%mass))
            IF (2.0_num * pi / omega .LT. min_dt) min_dt = 2.0_num * pi / omega
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(min_dt, dt_plasma_frequency, 1, mpireal, MPI_MIN, &
        comm, errcode)
    ! Must resolve plasma frequency
    dt_plasma_frequency = dt_plasma_frequency / 2.0_num

    DO ispecies = 1, n_species
      DEALLOCATE(initial_conditions(ispecies)%density)
      DEALLOCATE(initial_conditions(ispecies)%temp)
      DEALLOCATE(initial_conditions(ispecies)%drift)
    ENDDO
    DEALLOCATE(initial_conditions)

  END SUBROUTINE deallocate_ic



  SUBROUTINE non_uniform_load_particles(density, species, density_min, &
      density_max)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: density
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(INOUT) :: density_min, density_max
#ifndef PER_PARTICLE_WEIGHT
    INTEGER(KIND=8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(KIND=8) :: npart_per_cell
    REAL(num) :: density_total, density_total_global, density_average
    REAL(num) :: npart_per_cell_average
    INTEGER(KIND=8) :: npart_this_proc_new, ipart, npart_this_species
    INTEGER :: ix, iy, iz
    CHARACTER(LEN=15) :: string
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next

    partlist=>species%attached_list

    num_valid_cells_local = 0
    density_total = 0.0_num

    DO iz = -2, nz+3
      DO iy = -2, ny+3
        DO ix = -2, nx+3
          IF (density(ix,iy,iz) .GT. density_max) &
              density(ix,iy,iz) = density_max
        ENDDO
      ENDDO
    ENDDO

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          IF (density(ix,iy,iz) .GE. density_min) THEN
            num_valid_cells_local = num_valid_cells_local + 1
            density_total = density_total + density(ix,iy,iz)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(num_valid_cells_local, num_valid_cells_global, 1, &
        MPI_INTEGER8, MPI_SUM, comm, errcode)

    IF (species%npart_per_cell .GE. 0) THEN
      npart_per_cell_average = AINT(species%npart_per_cell, num)
    ELSE
      npart_per_cell_average = REAL(species%count, num) &
          / REAL(num_valid_cells_global, num)
    ENDIF

    IF (npart_per_cell_average .LE. 0) RETURN

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
        ENDDO
      ENDDO
    ENDDO

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)

    ! Randomly place npart_per_cell particles into each valid cell
    current=>partlist%head
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          npart_per_cell = NINT(density(ix, iy, iz) / density_average &
              * npart_per_cell_average)

          ipart = 0
          DO WHILE(ASSOCIATED(current) .AND. ipart .LT. npart_per_cell)
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
            current=>current%next
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next=>current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current=>next
    ENDDO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species
    species%weight = density_total_global * dx * dy * dz / npart_this_species

    IF (rank .EQ. 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
    ENDIF

    CALL particle_bcs
#else
    IF (rank .EQ. 0) THEN
      WRITE(*,*) 'non_uniform_load_particles() only available when using', &
          ' per species weighting'
    ENDIF
    CALL MPI_ABORT(comm, errcode, errcode)
#endif

  END SUBROUTINE non_uniform_load_particles



  ! This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE load_particles(species, load_list)

    TYPE(particle_species), POINTER :: species
    LOGICAL, DIMENSION(-2:,-2:,-2:), INTENT(IN) :: load_list
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: valid_cell_list
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER(KIND=8) :: ipart, npart_per_cell, num_int, num_total, idx
    INTEGER(KIND=8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(KIND=8) :: npart_this_species, num_new_particles, npart_left
    INTEGER(KIND=8), ALLOCATABLE :: num_valid_cells_all(:), num_idx(:)
    REAL(num) :: valid_cell_frac, num_real, f0, f1
    REAL(num), ALLOCATABLE :: num_frac(:)
    INTEGER :: cell_x
    INTEGER :: cell_y
    INTEGER :: cell_z
    INTEGER(KIND=8) :: i, ipos
    INTEGER :: ierr, ix, iy, iz
    CHARACTER(LEN=15) :: string
    LOGICAL :: sweep

    npart_this_species = species%count
    IF (npart_this_species .LT. 0) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*,*) 'Unable to continue, species ', &
          TRIM(species%name), ' has not had a number of particles set'
      ENDIF
      CALL MPI_ABORT(comm, errcode, ierr)
    ELSE IF (npart_this_species .EQ. 0) THEN
      RETURN
    ENDIF

    num_valid_cells_local = 0
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          IF (load_list(ix, iy, iz)) &
              num_valid_cells_local = num_valid_cells_local + 1
        ENDDO
      ENDDO
    ENDDO

    IF (species%npart_per_cell .GE. 0) THEN
      npart_per_cell = AINT(species%npart_per_cell, KIND=8)
      num_new_particles = &
          AINT(species%npart_per_cell * num_valid_cells_local, KIND=8)
    ELSE
      ALLOCATE(num_valid_cells_all(nproc), num_idx(nproc), num_frac(nproc))

      ! Calculate global number of particles per cell
      CALL MPI_ALLGATHER(num_valid_cells_local, 1, MPI_INTEGER8, &
          num_valid_cells_all, 1, MPI_INTEGER8, comm, errcode)

      num_valid_cells_global = 0
      DO i = 1,nproc
        num_valid_cells_global = num_valid_cells_global + num_valid_cells_all(i)
      ENDDO

      IF (num_valid_cells_global .EQ. 0) THEN
        IF (rank .EQ. 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Intial condition settings mean that there are no cells ' &
              // 'where particles may'
          WRITE(*,*) 'validly be placed for at least one species. Code will ' &
              // 'now terminate.'
          CALL MPI_ABORT(comm, errcode, ierr)
        ENDIF
      ENDIF

      valid_cell_frac = REAL(num_valid_cells_local, num) &
          / REAL(num_valid_cells_global, num)
      num_real = npart_this_species * valid_cell_frac
      num_new_particles = AINT(num_real, KIND=8)

      ! Work out which processors get the remaining fractional numbers
      ! of particles

      ! Get a list of the fractional part on each processor, along with
      ! the total
      num_total = 0
      DO i = 1,nproc
        valid_cell_frac = REAL(num_valid_cells_all(i), num) &
            / REAL(num_valid_cells_global, num)
        num_real = npart_this_species * valid_cell_frac
        num_int = AINT(num_real, KIND=8)
        num_frac(i) = num_real - num_int
        num_idx (i) = i
        num_total = num_total + num_int
      ENDDO
      num_total = npart_this_species - num_total

      IF (num_total .GT. 0) THEN
        ! Sort the list of fractions into decreasing order using bubble sort
        sweep = .TRUE.
        DO WHILE(sweep)
          sweep = .FALSE.
          f0 = num_frac(1)
          DO i = 2,nproc
            f1 = num_frac(i)
            IF (f1 .GT. f0) THEN
              num_frac(i-1) = f1
              num_frac(i) = f0
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

        DO i = 1,nproc-1
          num_int = AINT(num_frac(i) / (num_frac(i+1)+EPSILON(1.0_num)), KIND=8)
          IF (num_idx(i) .EQ. rank) THEN
            num_new_particles = num_new_particles + num_int
            EXIT
          ENDIF
          num_total = num_total - num_int
          IF (num_total .LE. 0) EXIT
        ENDDO
      ENDIF

      DEALLOCATE(num_valid_cells_all, num_idx, num_frac)

      species%npart_per_cell = &
          REAL(npart_this_species,num) / num_valid_cells_global
      npart_per_cell = AINT(species%npart_per_cell, KIND=8)
    ENDIF

    partlist=>species%attached_list

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    ! Randomly place npart_per_cell particles into each valid cell
    npart_left = num_new_particles
    current=>partlist%head
    IF (npart_per_cell .GT. 0) THEN

      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            IF (.NOT. load_list(ix, iy, iz)) CYCLE

            ipart = 0
            DO WHILE(ASSOCIATED(current) .AND. ipart .LT. npart_per_cell)
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
              current=>current%next

              ! One particle sucessfully placed
              npart_left = npart_left - 1
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    ! When num_new_particles does not equal
    ! npart_per_cell * num_valid_cells_local there will be particles left
    ! over that didn't get placed.
    ! The following loop randomly place remaining particles into valid cells.
    IF (npart_left .GT. 0) THEN
      ALLOCATE(valid_cell_list(num_valid_cells_local))

      ipos = 0
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            IF (load_list(ix,iy,iz)) THEN
              ipos = ipos + 1
              valid_cell_list(ipos) = ix - 1 + nx * (iy - 1 + ny * (iz - 1))
            ENDIF
          ENDDO
        ENDDO
      ENDDO

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

        current=>current%next
      ENDDO

      DEALLOCATE(valid_cell_list)
    ENDIF

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next=>current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current=>next
    ENDDO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species

    IF (rank .EQ. 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
    ENDIF

    CALL particle_bcs

  END SUBROUTINE load_particles



  SUBROUTINE setup_particle_density(density_in, species, density_min, &
      density_max)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: density_in
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(IN) :: density_min, density_max
#ifdef PER_PARTICLE_WEIGHT
    REAL(num) :: weight_local
    TYPE(particle), POINTER :: current
    INTEGER(KIND=8) :: ipart
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: weight_fn
    REAL(num) :: wdata
    TYPE(particle_list), POINTER :: partlist
    INTEGER :: ix, iy, iz, i, j, k, isubx, isuby, isubz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: density
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: density_map
#include "particle_head.inc"

    ALLOCATE(density(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(density_map(-2:nx+3,-2:ny+3,-2:nz+3))
    density = density_in
    density_map = .FALSE.

    CALL field_bc(density)

    DO iz = -2, nz+3
      DO iy = -2, ny+3
        DO ix = -2, nx+3
          IF (density(ix,iy,iz) .GT. density_max) &
              density(ix,iy,iz) = density_max
        ENDDO
      ENDDO
    ENDDO

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          IF (density(ix,iy,iz) .GE. density_min) density_map(ix,iy,iz) = .TRUE.
        ENDDO
      ENDDO
    ENDDO

    ! Uniformly load particles in space
    CALL load_particles(species, density_map)

    ALLOCATE(weight_fn(-2:nx+3,-2:ny+3,-2:nz+3))
    CALL MPI_BARRIER(comm, errcode)
    weight_fn = 0.0_num

    partlist=>species%attached_list
    ! If using per particle weighing then use the weight function to match the
    ! uniform pseudoparticle density to the real particle density
    current=>partlist%head
    ipart = 0
    ! First loop converts number density into weight function
    DO WHILE(ipart .LT. partlist%count)
      IF (.NOT. ASSOCIATED(current)) PRINT *, "Bad Particle"

#include "particle_to_grid.inc"

      DO isubz = sf_min, sf_max
        i = cell_x
        j = cell_y
        k = cell_z + isubz
#ifdef PARTICLE_SHAPE_TOPHAT
        IF (.NOT. density_map(i,j,k)) k = cell_z + 1 - isubz
#else
        IF (.NOT. density_map(i,j,k)) k = cell_z - isubz / 2
#endif
        DO isuby = sf_min, sf_max
          i = cell_x
          j = cell_y + isuby
#ifdef PARTICLE_SHAPE_TOPHAT
          IF (.NOT. density_map(i,j,k)) j = cell_y + 1 - isuby
#else
          IF (.NOT. density_map(i,j,k)) j = cell_y - isuby / 2
#endif
          DO isubx = sf_min, sf_max
            i = cell_x + isubx
#ifdef PARTICLE_SHAPE_TOPHAT
            IF (.NOT. density_map(i,j,k)) i = cell_x + 1 - isubx
#else
            IF (.NOT. density_map(i,j,k)) i = cell_x - isubx / 2
#endif
            weight_fn(i,j,k) = weight_fn(i,j,k) &
                + gx(isubx) * gy(isuby) * gz(isubz)
          ENDDO
        ENDDO
      ENDDO

      current=>current%next
      ipart = ipart + 1
    ENDDO
    DEALLOCATE(density_map)

    CALL processor_summation_bcs(weight_fn)
    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic) THEN
      IF (x_min_boundary) weight_fn(0 ,:,:) = weight_fn(1   ,:,:)
      IF (x_max_boundary) weight_fn(nx,:,:) = weight_fn(nx-1,:,:)
    ENDIF
    IF (bc_particle(c_bd_y_min) .NE. c_bc_periodic) THEN
      IF (y_min_boundary) weight_fn(:,0 ,:) = weight_fn(:,1   ,:)
      IF (y_max_boundary) weight_fn(:,ny,:) = weight_fn(:,ny-1,:)
    ENDIF
    IF (bc_particle(c_bd_z_min) .NE. c_bc_periodic) THEN
      IF (z_min_boundary) weight_fn(:,:,0 ) = weight_fn(:,:,1   )
      IF (z_max_boundary) weight_fn(:,:,nz) = weight_fn(:,:,nz-1)
    ENDIF
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(weight_fn, c_stagger_centre, ix)
    ENDDO

    wdata = dx * dy * dz
    DO iz = -2, nz+3
      DO iy = -2, ny+3
        DO ix = -2, nx+3
          IF (weight_fn(ix, iy, iz) .GT. 0.0_num) THEN
            weight_fn(ix, iy, iz) = &
                wdata * density(ix, iy, iz) / weight_fn(ix, iy, iz)
          ELSE
            weight_fn(ix, iy, iz) = 0.0_num
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(density)

    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic) THEN
      IF (x_min_boundary) weight_fn(0 ,:,:) = weight_fn(1   ,:,:)
      IF (x_max_boundary) weight_fn(nx,:,:) = weight_fn(nx-1,:,:)
    ENDIF
    IF (bc_particle(c_bd_y_min) .NE. c_bc_periodic) THEN
      IF (y_min_boundary) weight_fn(:,0 ,:) = weight_fn(:,1   ,:)
      IF (y_max_boundary) weight_fn(:,ny,:) = weight_fn(:,ny-1,:)
    ENDIF
    IF (bc_particle(c_bd_z_min) .NE. c_bc_periodic) THEN
      IF (z_min_boundary) weight_fn(:,:,0 ) = weight_fn(:,:,1   )
      IF (z_max_boundary) weight_fn(:,:,nz) = weight_fn(:,:,nz-1)
    ENDIF
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(weight_fn, c_stagger_centre, ix)
    ENDDO

    partlist=>species%attached_list
    ! Second loop actually assigns weights to particles
    ! Again assumes linear interpolation
    current=>partlist%head
    ipart = 0
    DO WHILE(ipart .LT. partlist%count)
#include "particle_to_grid.inc"

      weight_local = 0.0_num
      DO isubz = sf_min, sf_max
        DO isuby = sf_min, sf_max
          DO isubx = sf_min, sf_max
            weight_local = weight_local + gx(isubx) * gy(isuby) * gz(isubz) &
                * weight_fn(cell_x+isubx, cell_y+isuby, cell_z+isubz)
          ENDDO
        ENDDO
      ENDDO
      current%weight = weight_local

      current=>current%next
      ipart = ipart + 1
    ENDDO

    DEALLOCATE(weight_fn)
#else
    IF (rank .EQ. 0) THEN
      WRITE(*,*) 'setup_particle_density() only available when using', &
          ' per particle weighting'
    ENDIF
    CALL MPI_ABORT(comm, errcode, errcode)
#endif

  END SUBROUTINE setup_particle_density



  FUNCTION sample_dist_function(axis, dist_fn)

    REAL(num), DIMENSION(:), INTENT(IN) :: axis, dist_fn
    REAL(num), DIMENSION(:), ALLOCATABLE :: cdf
    REAL(num) :: position, d_cdf
    INTEGER :: n_points, ipoint, start, endpoint, current
    REAL(num) :: sample_dist_function

    n_points = SIZE(dist_fn)
    ALLOCATE(cdf(1:n_points))
    DO ipoint = 1, n_points
      cdf(ipoint) = SUM(dist_fn(1:ipoint))
    ENDDO
    cdf = cdf / SUM(dist_fn)

    position = random()
    sample_dist_function = 0.0_num

    start = 1
    endpoint = n_points
    current = (start+endpoint) / 2

    DO current = 1, n_points-1
      IF (cdf(current) .LE. position .AND. cdf(current+1) .GE. position) THEN
        d_cdf = cdf(current+1)-cdf(current)
        sample_dist_function = &
            (axis(current)*(position-cdf(current)) / d_cdf &
            + axis(current+1)*(cdf(current+1)-position) / d_cdf)
        EXIT
      ENDIF
    ENDDO

    DEALLOCATE(cdf)

  END FUNCTION sample_dist_function

END MODULE helper
