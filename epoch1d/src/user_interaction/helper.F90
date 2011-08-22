MODULE helper

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
      IF (move_window) THEN
        species_list(ispecies)%density = &
            initial_conditions(ispecies)%density(nx)
        species_list(ispecies)%temperature = &
            initial_conditions(ispecies)%temp(nx,:)
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

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies, ix
    REAL(num) :: min_dt, omega, k_max

    min_dt = 1000000.0_num
    k_max = 2.0_num * pi / dx

    ! Identify the plasma frequency (Bohm-Gross dispersion relation)
    ! Note that this doesn't get strongly relativistic plasmas right
    DO ispecies = 1, n_species
      DO ix = 1, nx
        omega = SQRT((initial_conditions(ispecies)%density(ix) * q0**2) &
            / (species_list(ispecies)%mass * epsilon0) &
            + 6.0_num * k_max**2 * kb &
            * MAXVAL(initial_conditions(ispecies)%temp(ix,:)) &
            / (species_list(ispecies)%mass))
        IF (2.0_num * pi / omega .LT. min_dt) min_dt = 2.0_num * pi / omega
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

    REAL(num), DIMENSION(-2:), INTENT(INOUT) :: density
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(INOUT) :: density_min, density_max
#ifndef PER_PARTICLE_WEIGHT
    INTEGER(KIND=8) :: num_valid_cells, num_valid_cells_global
    INTEGER(KIND=8) :: npart_per_cell
    REAL(num) :: npart_per_cell_average
    REAL(num) :: density_total, density_total_global, density_average
    INTEGER(KIND=8) :: npart_this_proc_new, ipart, npart_this_species
    INTEGER :: ix
    CHARACTER(LEN=15) :: string
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next

    partlist=>species%attached_list

    num_valid_cells = 0
    density_total = 0.0_num

    DO ix = -2, nx+3
      IF (density(ix) .GT. density_max) density(ix) = density_max
    ENDDO

    DO ix = 1, nx
      IF (density(ix) .GE. density_min) THEN
        num_valid_cells = num_valid_cells + 1
        density_total = density_total + density(ix)
      ENDIF
    ENDDO

    CALL MPI_ALLREDUCE(num_valid_cells, num_valid_cells_global, 1, &
        MPI_INTEGER8, MPI_MAX, comm, errcode)
    npart_per_cell_average = REAL(species%count, num) &
        / REAL(num_valid_cells_global, num)

    IF (npart_per_cell_average .LE. 0) RETURN

    CALL MPI_ALLREDUCE(density_total, density_total_global, 1, mpireal, &
        MPI_SUM, comm, errcode)
    density_average = density_total_global / REAL(num_valid_cells_global, num)

    npart_this_proc_new = 0
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix) / density_average &
          * npart_per_cell_average)
      npart_this_proc_new = npart_this_proc_new + npart_per_cell
    ENDDO

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)

    ! Randomly place npart_per_cell particles into each valid cell
    current=>partlist%head
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
        current=>current%next
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
    species%weight = density_total_global * dx / npart_this_species

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
    LOGICAL, DIMENSION(-2:), INTENT(IN) :: load_list
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: valid_cell_list
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER(KIND=8) :: ipart, npart_per_cell
    INTEGER(KIND=8) :: num_valid_cells, num_valid_cells_local
    INTEGER(KIND=8) :: npart_this_species, num_new_particles, npart_left
    REAL(num) :: valid_cell_frac
    INTEGER :: cell_x
    INTEGER(KIND=8) :: i, ipos
    INTEGER :: ierr, ix
    CHARACTER(LEN=15) :: string

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
    DO ix = 1, nx
      IF (load_list(ix)) num_valid_cells_local = num_valid_cells_local + 1
    ENDDO

    ! Calculate global number of particles per cell
    CALL MPI_ALLREDUCE(num_valid_cells_local, num_valid_cells, 1, &
        MPI_INTEGER8, MPI_SUM, comm, errcode)

    IF (num_valid_cells .EQ. 0) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*,*) '*** ERROR ***'
        WRITE(*,*) 'Intial condition settings mean that there are no cells ' &
            // 'where particles may'
        WRITE(*,*) 'validly be placed for at least one species. Code will ' &
            // 'now terminate.'
        CALL MPI_ABORT(comm, errcode, ierr)
      ENDIF
    ENDIF

    partlist=>species%attached_list

    valid_cell_frac = &
        REAL(num_valid_cells_local, num) / REAL(num_valid_cells, num)
    num_new_particles = NINT(npart_this_species*valid_cell_frac, KIND=8)

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    npart_per_cell = npart_this_species / num_valid_cells
    species%npart_per_cell = npart_per_cell
    IF (species%npart_per_cell .EQ. 0) species%npart_per_cell = 1

    ! Randomly place npart_per_cell particles into each valid cell
    npart_left = num_new_particles
    current=>partlist%head
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
          current=>current%next

          ! One particle sucessfully placed
          npart_left = npart_left - 1
        ENDDO
      ENDDO

    ENDIF

    ! When num_new_particles does not equal npart_per_cell * num_valid_cells
    ! there will be particles left over that didn't get placed.
    ! The following loop randomly place remaining particles into valid cells.
    IF (npart_left .GT. 0) THEN
      ALLOCATE(valid_cell_list(num_valid_cells_local))

      ipos = 0
      DO ix = 1, nx
        IF (load_list(ix)) THEN
          ipos = ipos + 1
          valid_cell_list(ipos) = ix - 1
        ENDIF
      ENDDO

      DO i = 1, npart_left
        ipos = INT(random() * (num_valid_cells_local - 1)) + 1
        ipos = valid_cell_list(ipos)

        cell_x = ipos + 1

        current%part_pos = x(cell_x) + (random() - 0.5_num) * dx

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

    REAL(num), DIMENSION(-2:), INTENT(IN) :: density_in
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(IN) :: density_min, density_max
#ifdef PER_PARTICLE_WEIGHT
    REAL(num) :: weight_local
    TYPE(particle), POINTER :: current
    INTEGER(KIND=8) :: ipart
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

    CALL field_bc(density)

    DO ix = -2, nx+3
      IF (density(ix) .GT. density_max) density(ix) = density_max
    ENDDO

    DO ix = 1, nx
      IF (density(ix) .GE. density_min) density_map(ix) = .TRUE.
    ENDDO

    ! Uniformly load particles in space
    CALL load_particles(species, density_map)

    ALLOCATE(weight_fn(-2:nx+3))
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

      DO isubx = sf_min, sf_max
        i = cell_x + isubx
#ifdef PARTICLE_SHAPE_TOPHAT
        IF (.NOT. density_map(i)) i = cell_x + 1 - isubx
#else
        IF (.NOT. density_map(i)) i = cell_x - isubx / 2
#endif
        weight_fn(i) = weight_fn(i) + gx(isubx)
      ENDDO

      current=>current%next
      ipart = ipart + 1
    ENDDO
    DEALLOCATE(density_map)

    CALL processor_summation_bcs(weight_fn)
    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic) THEN
      IF (x_min_boundary) weight_fn(0 ) = weight_fn(1   )
      IF (x_max_boundary) weight_fn(nx) = weight_fn(nx-1)
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
    ENDDO
    DEALLOCATE(density)

    IF (bc_particle(c_bd_x_min) .NE. c_bc_periodic) THEN
      IF (x_min_boundary) weight_fn(0 ) = weight_fn(1   )
      IF (x_max_boundary) weight_fn(nx) = weight_fn(nx-1)
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
      DO isubx = sf_min, sf_max
        weight_local = weight_local + gx(isubx) * weight_fn(cell_x+isubx)
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
