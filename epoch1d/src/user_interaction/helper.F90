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
      species => species_list(ispecies)

      ! Set temperature at boundary for thermal bcs.

      IF (bc_particle(c_bd_x_min) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_min(1:3) = &
            initial_conditions(ispecies)%temp(1,1:3)
      ENDIF
      IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
        species_list(ispecies)%ext_temp_x_max(1:3) = &
            initial_conditions(ispecies)%temp(nx,1:3)
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
          initial_conditions(ispecies)%temp(:,1), c_dir_x, species, &
          initial_conditions(ispecies)%drift(:,1))
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,2), c_dir_y, species, &
          initial_conditions(ispecies)%drift(:,2))
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,3), c_dir_z, species, &
          initial_conditions(ispecies)%drift(:,3))
    ENDDO

    IF (rank .EQ. 0) THEN
      DO ispecies = 1, n_species
        species => species_list(ispecies)
        IF (species%count .LT. 0) THEN
          WRITE(*,*) 'No particles specified for species ', TRIM(species%name)
          WRITE(stat_unit,*) &
              'No particles specified for species ', TRIM(species%name)
          species%count = 0
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE auto_load



  SUBROUTINE allocate_ic

    INTEGER :: ispecies

    ALLOCATE(initial_conditions(1:n_species))
    DO ispecies = 1, n_species
      ALLOCATE(initial_conditions(ispecies)%density(-2:nx+3))
      ALLOCATE(initial_conditions(ispecies)%temp (-2:nx+3,1:3))
      ALLOCATE(initial_conditions(ispecies)%drift(-2:nx+3,1:3))

      initial_conditions(ispecies)%density = 1.0_num
      initial_conditions(ispecies)%temp = 0.0_num
      initial_conditions(ispecies)%drift = 0.0_num
      initial_conditions(ispecies)%density_min = EPSILON(1.0_num)
      initial_conditions(ispecies)%density_max = HUGE(1.0_num)
    ENDDO

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies

    DO ispecies = 1, n_species
      DEALLOCATE(initial_conditions(ispecies)%density)
      DEALLOCATE(initial_conditions(ispecies)%temp)
      DEALLOCATE(initial_conditions(ispecies)%drift)
    ENDDO
    IF (.NOT. move_window) DEALLOCATE(initial_conditions)

  END SUBROUTINE deallocate_ic



#ifndef PER_PARTICLE_WEIGHT
  SUBROUTINE non_uniform_load_particles(density, species, density_min, &
      density_max)

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: density
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(INOUT) :: density_min, density_max
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_per_cell
    REAL(num) :: density_total, density_total_global, density_average
    REAL(num) :: npart_per_cell_average
    INTEGER(i8) :: npart_this_proc_new, ipart, npart_this_species
    INTEGER :: ix
    CHARACTER(LEN=15) :: string
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next

    partlist => species%attached_list

    num_valid_cells_local = 0
    density_total = 0.0_num

    DO ix = -2, nx+3
      IF (density(ix) .GT. density_max) density(ix) = density_max
    ENDDO ! ix

    DO ix = 1, nx
      IF (density(ix) .GE. density_min) THEN
        num_valid_cells_local = num_valid_cells_local + 1
        density_total = density_total + density(ix)
      ENDIF
    ENDDO ! ix

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
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix) / density_average &
          * npart_per_cell_average)
      npart_this_proc_new = npart_this_proc_new + npart_per_cell
    ENDDO ! ix

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)

    ! Randomly place npart_per_cell particles into each valid cell
    current => partlist%head
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix) / density_average &
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
        current%part_pos = x(ix) + (random() - 0.5_num) * dx

        ipart = ipart + 1
        current => current%next
      ENDDO
    ENDDO ! ix

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
    species%weight = density_total_global * dx / npart_this_species

    IF (rank .EQ. 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
    ENDIF

    CALL particle_bcs

  END SUBROUTINE non_uniform_load_particles
#endif



  ! This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE load_particles(species, load_list)

    TYPE(particle_species), POINTER :: species
    LOGICAL, DIMENSION(-2:), INTENT(IN) :: load_list
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
    INTEGER(i8) :: i, ipos
    INTEGER :: ierr, ix
    CHARACTER(LEN=15) :: string
    LOGICAL :: sweep

    npart_this_species = species%count
    IF (npart_this_species .LE. 0) RETURN

    num_valid_cells_local = 0
    DO ix = 1, nx
      IF (load_list(ix)) num_valid_cells_local = num_valid_cells_local + 1
    ENDDO ! ix

    IF (species%npart_per_cell .GE. 0) THEN
      npart_per_cell = AINT(species%npart_per_cell, KIND=i8)
      num_new_particles = &
          AINT(species%npart_per_cell * num_valid_cells_local, KIND=i8)
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
          WRITE(*,*) 'validly be placed for species "' // TRIM(species%name) &
              // '". ', 'Code will now terminate.'
          CALL MPI_ABORT(comm, errcode, ierr)
        ENDIF
      ENDIF

      valid_cell_frac = REAL(num_valid_cells_local, num) &
          / REAL(num_valid_cells_global, num)
      num_real = npart_this_species * valid_cell_frac
      num_new_particles = AINT(num_real, KIND=i8)

      ! Work out which processors get the remaining fractional numbers
      ! of particles

      ! Get a list of the fractional part on each processor, along with
      ! the total
      num_total = 0
      DO i = 1,nproc
        valid_cell_frac = REAL(num_valid_cells_all(i), num) &
            / REAL(num_valid_cells_global, num)
        num_real = npart_this_species * valid_cell_frac
        num_int = AINT(num_real, KIND=i8)
        num_frac(i) = num_real - num_int
        num_idx (i) = i - 1
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
          IF (num_idx(i) .EQ. rank) THEN
            num_new_particles = num_new_particles + 1
            EXIT
          ENDIF
          num_total = num_total - 1
          IF (num_total .LE. 0) EXIT
        ENDDO
      ENDIF

      DEALLOCATE(num_valid_cells_all, num_idx, num_frac)

      species%npart_per_cell = &
          REAL(npart_this_species,num) / num_valid_cells_global
      npart_per_cell = AINT(species%npart_per_cell, KIND=i8)
    ENDIF

    partlist => species%attached_list

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    ! Randomly place npart_per_cell particles into each valid cell
    npart_left = num_new_particles
    current => partlist%head
    IF (npart_per_cell .GT. 0) THEN

      DO ix = 1, nx
        IF (.NOT. load_list(ix)) CYCLE

        ipart = 0
        DO WHILE(ASSOCIATED(current) .AND. ipart .LT. npart_per_cell)
#ifdef PER_PARTICLE_CHARGE_MASS
          ! Even if particles have per particle charge and mass, assume
          ! that initially they all have the same charge and mass (user
          ! can easily over_ride)
          current%charge = species%charge
          current%mass = species%mass
#endif
          current%part_pos = x(ix) + (random() - 0.5_num) * dx

          ipart = ipart + 1
          current => current%next

          ! One particle sucessfully placed
          npart_left = npart_left - 1
        ENDDO
      ENDDO ! ix

    ENDIF

    ! When num_new_particles does not equal
    ! npart_per_cell * num_valid_cells_local there will be particles left
    ! over that didn't get placed.
    ! The following loop randomly place remaining particles into valid cells.
    IF (npart_left .GT. 0) THEN
      ALLOCATE(valid_cell_list(num_valid_cells_local))

      ipos = 0
      DO ix = 1, nx
        IF (load_list(ix)) THEN
          ipos = ipos + 1
          valid_cell_list(ipos) = ix - 1
        ENDIF
      ENDDO ! ix

      DO i = 1, npart_left
        ipos = INT(random() * (num_valid_cells_local - 1)) + 1
        ipos = valid_cell_list(ipos)

        cell_x = ipos + 1

        current%part_pos = x(cell_x) + (random() - 0.5_num) * dx

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

    IF (rank .EQ. 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
      WRITE(stat_unit,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', TRIM(species%name)
    ENDIF

    CALL particle_bcs

  END SUBROUTINE load_particles



#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE setup_particle_density(density_in, species, density_min, &
      density_max)

    REAL(num), DIMENSION(-2:), INTENT(IN) :: density_in
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(IN) :: density_min, density_max
    REAL(num) :: weight_local
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    REAL(num), DIMENSION(:), ALLOCATABLE :: weight_fn
    REAL(num) :: wdata
    TYPE(particle_list), POINTER :: partlist
    INTEGER :: ix, i, isubx
    REAL(num), DIMENSION(:), ALLOCATABLE :: density
    LOGICAL, DIMENSION(:), ALLOCATABLE :: density_map
#include "particle_head.inc"

    ALLOCATE(density(-2:nx+3))
    ALLOCATE(density_map(-2:nx+3))
    density = density_in
    density_map = .FALSE.

    CALL field_bc(density, ng)

    DO ix = -2, nx+3
      IF (density(ix) .GT. density_max) density(ix) = density_max
    ENDDO ! ix

    DO ix = 1, nx
      IF (density(ix) .GE. density_min) density_map(ix) = .TRUE.
    ENDDO ! ix

    ! Uniformly load particles in space
    CALL load_particles(species, density_map)

    ALLOCATE(weight_fn(-2:nx+3))
    CALL MPI_BARRIER(comm, errcode)
    weight_fn = 0.0_num

    partlist => species%attached_list
    ! If using per particle weighing then use the weight function to match the
    ! uniform pseudoparticle density to the real particle density
    current => partlist%head
    ipart = 0
    ! First loop converts number density into weight function
    DO WHILE(ipart .LT. partlist%count)
      IF (.NOT. ASSOCIATED(current)) PRINT *, 'Bad Particle'

#include "particle_to_grid.inc"

      DO isubx = sf_min, sf_max
        i = cell_x + isubx
#ifdef PARTICLE_SHAPE_TOPHAT
        IF (.NOT. density_map(i)) i = cell_x + 1 - isubx
#else
        IF (.NOT. density_map(i)) i = cell_x - isubx / 2
#endif
        weight_fn(i) = weight_fn(i) + gx(isubx)
      ENDDO

      current => current%next
      ipart = ipart + 1
    ENDDO
    DEALLOCATE(density_map)

    CALL processor_summation_bcs(weight_fn, ng)
    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic) THEN
      IF (x_min_boundary) weight_fn(0   ) = weight_fn(1 )
      IF (x_max_boundary) weight_fn(nx+1) = weight_fn(nx)
    ENDIF
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(weight_fn, c_stagger_centre, ix)
    ENDDO

    wdata = dx
    DO ix = -2, nx+3
      IF (weight_fn(ix) .GT. 0.0_num) THEN
        weight_fn(ix) = wdata * density(ix) / weight_fn(ix)
      ELSE
        weight_fn(ix) = 0.0_num
      ENDIF
    ENDDO ! ix
    DEALLOCATE(density)

    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic) THEN
      IF (x_min_boundary) weight_fn(0   ) = weight_fn(1 )
      IF (x_max_boundary) weight_fn(nx+1) = weight_fn(nx)
    ENDIF
    DO ix = 1, 2*c_ndims
      CALL field_zero_gradient(weight_fn, c_stagger_centre, ix)
    ENDDO

    partlist => species%attached_list
    ! Second loop actually assigns weights to particles
    ! Again assumes linear interpolation
    current => partlist%head
    ipart = 0
    DO WHILE(ipart .LT. partlist%count)
#include "particle_to_grid.inc"

      weight_local = 0.0_num
      DO isubx = sf_min, sf_max
        weight_local = weight_local + gx(isubx) * weight_fn(cell_x+isubx)
      ENDDO ! isubx
      current%weight = weight_local

      current => current%next
      ipart = ipart + 1
    ENDDO

    DEALLOCATE(weight_fn)

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
      IF (cdf(current) .LE. position .AND. cdf(current+1) .GE. position) THEN
        d_cdf = cdf(current+1) - cdf(current)
        sample_dist_function = (axis(current) * (position - cdf(current)) &
            + axis(current+1) * (cdf(current+1) - position)) / d_cdf
        EXIT
      ENDIF
    ENDDO

    DEALLOCATE(cdf)

  END FUNCTION sample_dist_function

END MODULE helper
