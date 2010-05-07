MODULE helper

  USE boundary
  USE shape_functions
  USE strings

  IMPLICIT NONE

  SAVE
  REAL(num) :: max_rand

CONTAINS

  SUBROUTINE auto_load

    INTEGER :: ispecies
    INTEGER :: clock, idum
    TYPE(particle_family), POINTER :: part_family

    CALL SYSTEM_CLOCK(clock)
    idum = -(clock+rank+1)
    DO ispecies = 1, n_species
      part_family=>particle_species(ispecies)
      IF (move_window) THEN
        particle_species(ispecies)%density = &
            initial_conditions(ispecies)%rho(nx,:,:)
        particle_species(ispecies)%temperature = &
            initial_conditions(ispecies)%temp(nx,:,:,:)
      ENDIF
#ifdef PER_PARTICLE_WEIGHT
      CALL setup_particle_density(initial_conditions(ispecies)%rho, &
          part_family, initial_conditions(ispecies)%minrho, &
          initial_conditions(ispecies)%maxrho, idum)
#else
      CALL non_uniform_load_particles(initial_conditions(ispecies)%rho, &
          part_family, initial_conditions(ispecies)%minrho, &
          initial_conditions(ispecies)%maxrho, idum)
#endif
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,:,:,1), c_dir_x, part_family, &
          initial_conditions(ispecies)%drift(:,:,:,1), idum)
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,:,:,2), c_dir_y, part_family, &
          initial_conditions(ispecies)%drift(:,:,:,2), idum)
      CALL setup_particle_temperature(&
          initial_conditions(ispecies)%temp(:,:,:,3), c_dir_z, part_family, &
          initial_conditions(ispecies)%drift(:,:,:,3), idum)
    ENDDO

  END SUBROUTINE auto_load



  SUBROUTINE allocate_ic

    INTEGER :: ispecies

    ALLOCATE(initial_conditions(1:n_species))
    DO ispecies = 1, n_species
      ALLOCATE(initial_conditions(ispecies)%rho  (-2:nx+3,-2:ny+3,-2:nz+3))
      ALLOCATE(initial_conditions(ispecies)%temp (-2:nx+3,-2:ny+3,-2:nz+3,1:3))
      ALLOCATE(initial_conditions(ispecies)%drift(-2:nx+3,-2:ny+3,-2:nz+3,1:3))

      initial_conditions(ispecies)%rho = 1.0_num
      initial_conditions(ispecies)%temp = 0.0_num
      initial_conditions(ispecies)%minrho = 0.0_num
      initial_conditions(ispecies)%maxrho = 0.0_num
      initial_conditions(ispecies)%drift = 0.0_num
    ENDDO

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies, ix, iy
    REAL(num) :: min_dt, omega, k_max

    min_dt = 1000000.0_num
    k_max = 2.0_num * pi / MIN(dx, dy)

    ! Identify the plasma frequency (Bohm-Gross dispersion relation)
    ! Note that this doesn't get strongly relativistic plasmas right
    DO ispecies = 1, n_species
      DO iy = 1, ny
        DO ix = 1, nx
          omega = SQRT((initial_conditions(ispecies)%rho(ix, iy, iz) * q0**2) &
              / (particle_species(ispecies)%mass * epsilon0) &
              + 3.0_num * k_max**2 * kb &
              * MAXVAL(initial_conditions(ispecies)%temp(ix, iy, iz,:)) &
              / (particle_species(ispecies)%mass))
          IF (2.0_num * pi/omega .LT. min_dt) min_dt = 2.0_num * pi /omega
        ENDDO
      ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(min_dt, dt_plasma_frequency, 1, mpireal, MPI_MIN, &
        comm, errcode)
    ! Must resolve plasma frequency
    dt_plasma_frequency = dt_plasma_frequency/2.0_num

    DO ispecies = 1, n_species
      DEALLOCATE(initial_conditions(ispecies)%rho)
      DEALLOCATE(initial_conditions(ispecies)%temp)
      DEALLOCATE(initial_conditions(ispecies)%drift)
    ENDDO
    DEALLOCATE(initial_conditions)

  END SUBROUTINE deallocate_ic



  SUBROUTINE non_uniform_load_particles(density, species_list, minrho, &
      maxrho, idum)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: density
    TYPE(particle_family), POINTER :: species_list
    REAL(num), INTENT(INOUT) :: minrho, maxrho
    INTEGER, INTENT(INOUT) :: idum

    INTEGER(KIND=8) :: num_valid_cells, num_valid_cells_global
    INTEGER(KIND=8) :: npart_per_cell_average
    INTEGER(KIND=8) :: npart_per_cell
    REAL(num) :: density_total, density_total_global, density_average
    INTEGER(KIND=8) :: npart_this_proc_new, ipart, npart_this_species
    REAL(num) :: rpos

    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next

    partlist=>species_list%attached_list

    num_valid_cells = 0
    density_total = 0.0_num
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          IF (density(ix, iy, iz) .GE. minrho) THEN
            num_valid_cells = num_valid_cells+1
            density_total = density_total + density(ix, iy, iz)
          ENDIF
          IF (density(ix, iy, iz) .GT. maxrho .AND. maxrho .GT. 0.0_num) &
              density(ix, iy, iz) = maxrho
        ENDDO
      ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(num_valid_cells, num_valid_cells_global, 1, &
        MPI_INTEGER8, MPI_MAX, comm, errcode)
    npart_per_cell_average = species_list%count/num_valid_cells_global
    IF (npart_per_cell_average .EQ. 0) npart_per_cell_average = 1

    CALL MPI_ALLREDUCE(density_total, density_total_global, 1, mpireal, &
        MPI_SUM, comm, errcode)
    density_average = density_total_global/REAL(num_valid_cells_global, num)

    ! Assume that a cell with the average density has the average number of
    ! particles per cell. Now calculate the new minimum density
    minrho = density_average/REAL(npart_per_cell_average, num)
    ! Set the particle weight
    weight = minrho * dx * dy *dz

    ! Recalculate the number of valid cells and the summed density
    num_valid_cells = 0
    density_total = 0.0_num
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
          IF (density(ix, iy, iz) .GE. minrho) THEN
            num_valid_cells = num_valid_cells+1
            density_total = density_total + density(ix, iy, iz)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    npart_this_proc_new = &
        INT(density_total / density_average * REAL(npart_per_cell_average, num))

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)
    current=>partlist%head
    DO iz = 1, nz
      DO ix = 1, nx
        DO iy = 1, ny
          ipart = 0
          npart_per_cell = INT(density(ix, iy, iz) / density_average &
              * REAL(npart_per_cell_average, num))
          DO WHILE(ASSOCIATED(current) .AND. ipart .LT. npart_per_cell)
#ifdef PER_PARTICLE_CHARGEMASS
            ! Even if particles have per particle charge and mass, assume
            ! that initially they all have the same charge and mass (user
            ! can easily over_ride)
            current%charge = species_list%charge
            current%mass = species_list%mass
#endif
            rpos = random(idum)-0.5_num
            rpos = (rpos*dx)+x(ix)
            current%part_pos(1) = rpos
            rpos = random(idum)-0.5_num
            rpos = (rpos*dy)+y(iy)
            current%part_pos(2) = rpos
            rpos = random(idum)-0.5_num
            rpos = (rpos*dz)+z(iz)
            current%part_pos(3) = rpos
            ipart = ipart+1
            current=>current%next
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO WHILE(ASSOCIATED(current))
      next=>current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current=>next
    ENDDO
    CALL MPI_REDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, 0, comm, errcode)
    species_list%count = npart_this_species
    IF (rank .EQ. 0) THEN
      PRINT *, "Loaded", npart_this_species, "particles of species ", &
          TRIM(species_list%name)
      WRITE(20, *) "Loaded", npart_this_species, "particles of species ", &
          TRIM(species_list%name)
    ENDIF
    CALL particle_bcs

  END SUBROUTINE non_uniform_load_particles



  ! This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE load_particles(species_list, load_list, idum)

    INTEGER, INTENT(INOUT) :: idum
    TYPE(particle_family), POINTER :: species_list
    LOGICAL, DIMENSION(-2:,-2:,-2:), INTENT(IN) :: load_list
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER(KIND=8) :: ipart, npart_per_cell
    INTEGER(KIND=8) :: num_valid_cells, num_valid_cells_local
    INTEGER(KIND=8) :: npart_this_species, num_new_particles, npart_left
    REAL(num) :: valid_cell_frac
    REAL(dbl) :: rpos
    INTEGER :: lower_x, upper_x, cell_x
    INTEGER :: lower_y, upper_y, cell_y
    INTEGER :: lower_z, upper_z, cell_z
    REAL(num) :: cell_x_r
    REAL(num) :: cell_y_r
    REAL(num) :: cell_z_r
    INTEGER(KIND=8) :: i
    INTEGER :: j, ierr
    CHARACTER(LEN=15) :: string

    lower_x = 1
    lower_y = 1
    lower_z = 1
    upper_x = nx
    upper_y = ny
    upper_z = nz

    IF (coordinates(3) .EQ. nprocx-1) upper_x = nx+1
    IF (coordinates(3) .EQ. 0) lower_x = 0
    IF (coordinates(2) .EQ. nprocy-1) upper_y = ny+1
    IF (coordinates(2) .EQ. 0) lower_y = 0
    IF (coordinates(1) .EQ. nprocz-1) upper_z = nz+1
    IF (coordinates(1) .EQ. 0) lower_z = 0

    partlist=>species_list%attached_list

    npart_this_species = species_list%count
    IF (npart_this_species .LT. 0) THEN
      IF (rank .EQ. 0) PRINT *, "Unable to continue, species ", &
          TRIM(species_list%name), " has not had a number of particles set"
      CALL MPI_ABORT(comm, errcode, ierr)
    ENDIF
    IF (npart_this_species .EQ. 0) RETURN
    num_valid_cells_local = 0
    DO iz = lower_z, upper_z
      DO iy = lower_y, upper_y
        DO ix = lower_x, upper_x
          IF (load_list(ix,iy,iz)) &
              num_valid_cells_local = num_valid_cells_local+1
        ENDDO
      ENDDO
    ENDDO

    ! Calculate global number of particles per cell
    CALL MPI_ALLREDUCE(num_valid_cells_local, num_valid_cells, 1, &
        MPI_INTEGER8, MPI_SUM, comm, errcode)

    IF (num_valid_cells .EQ. 0) THEN
      IF (rank .EQ. 0) THEN
        WRITE(*,*) '***ERROR***'
        WRITE(*,*) 'Intial condition settings mean that there are no cells ' &
            // 'where particles may'
        WRITE(*,*) 'validly be placed for at least one species. Code will ' &
            // 'now terminate.'
        CALL MPI_ABORT(comm, errcode, ierr)
      ENDIF
    ENDIF

    valid_cell_frac = &
        REAL(num_valid_cells_local, num)/REAL(num_valid_cells, num)
    num_new_particles = INT(npart_this_species*valid_cell_frac, KIND=8)
    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    npart_per_cell = npart_this_species/num_valid_cells
    species_list%npart_per_cell = npart_per_cell
    IF (species_list%npart_per_cell .EQ. 0) species_list%npart_per_cell = 1

    ipart = 0
    ix = 1
    iy = 1

    npart_left = npart_this_species
    current=>partlist%head
    IF (npart_per_cell .GT. 0) THEN

      DO iz = lower_z, upper_z
        DO iy = lower_y, upper_y
          DO ix = lower_x, upper_x
            ipart = 0
            IF (load_list(ix, iy, iz)) THEN
              DO WHILE(ASSOCIATED(current) .AND. ipart .LT. npart_per_cell)
#ifdef PER_PARTICLE_CHARGEMASS
                ! Even if particles have per particle charge and mass, assume
                ! that initially they all have the same charge and mass (user
                ! can easily over_ride)
                current%charge = species_list%charge
                current%mass = species_list%mass
#endif
                rpos = random(idum)-0.5_num
                rpos = (rpos*dx)+x(ix)
                current%part_pos(1) = rpos
                rpos = random(idum)-0.5_num
                rpos = (rpos*dy)+y(iy)
                current%part_pos(2) = rpos
                rpos = random(idum)-0.5_num
                rpos = (rpos*dz)+z(iz)
                current%part_pos(3) = rpos
                ipart = ipart+1
                current=>current%next
                ! One particle sucessfully placed
                npart_left = npart_left-1
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    DO i = 1, npart_left*4
      IF (.NOT. ASSOCIATED(current)) EXIT
      DO j = 1, 200
        rpos = random(idum)*(x(nx)-x(1) + dx) - dx/2.0_num
        current%part_pos(1) = rpos+x(1)
        rpos = random(idum)*(y(ny)-y(1) + dy) - dy/2.0_num
        current%part_pos(2) = rpos+y(1)
        rpos = random(idum)*(z(nz)-z(1) + dz) - dz/2.0_num
        current%part_pos(3) = rpos+z(1)

        cell_x_r = (current%part_pos(1)-x_min_local)/dx - 0.5_num
        cell_x = NINT(cell_x_r)
        cell_x = cell_x+1

        cell_y_r = (current%part_pos(2)-y_min_local)/dy - 0.5_num
        cell_y = NINT(cell_y_r)
        cell_y = cell_y+1

        cell_z_r = (current%part_pos(3)-z_min_local)/dz - 0.5_num
        cell_z = NINT(cell_z_r)
        cell_z = cell_z+1

        IF (load_list(cell_x, cell_y, cell_z)) THEN
          EXIT
        ENDIF
      ENDDO
      current=>current%next
    ENDDO

    DO WHILE(ASSOCIATED(current))
      next=>current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current=>next
    ENDDO
    CALL MPI_REDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, 0, comm, errcode)
    species_list%count = npart_this_species
    IF (rank .EQ. 0) THEN
      CALL integer_as_string(npart_this_species, string)
      PRINT *, "Loaded ", TRIM(ADJUSTL(string)), " particles of species ", &
          TRIM(species_list%name)
      WRITE(20, *) "Loaded ", TRIM(ADJUSTL(string)), " particles of species ", &
          TRIM(species_list%name)
    ENDIF
    CALL particle_bcs

  END SUBROUTINE load_particles



  ! Subroutine to initialise a thermal particle distribution
  ! Assumes linear interpolation of temperature between cells
  SUBROUTINE setup_particle_temperature(temperature, direction, part_family, &
      drift, idum)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: temperature
    INTEGER, INTENT(IN) :: direction
    TYPE(particle_family), POINTER :: part_family
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: drift
    INTEGER, INTENT(INOUT) :: idum
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local, drift_local
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cell_z_r, cell_frac_z
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    TYPE(particle), POINTER :: current
    INTEGER :: cell_x, cell_y, cell_z
    INTEGER(KIND=8) :: ipart

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
      cell_x_r = (current%part_pos(1)-x_min_local)/dx - 0.5_num
      cell_x = NINT(cell_x_r)
      cell_frac_x = REAL(cell_x, num) - cell_x_r
      cell_x = cell_x+1

      cell_y_r = (current%part_pos(2)-y_min_local)/dy - 0.5_num
      cell_y = NINT(cell_y_r)
      cell_frac_y = REAL(cell_y, num) - cell_y_r
      cell_y = cell_y+1

      cell_z_r = (current%part_pos(3)-z_min_local)/dz - 0.5_num
      cell_z = NINT(cell_z_r)
      cell_frac_z = REAL(cell_z, num) - cell_z_r
      cell_z = cell_z+1

      CALL grid_to_particle(cell_frac_x, gx)
      CALL grid_to_particle(cell_frac_y, gy)
      CALL grid_to_particle(cell_frac_z, gz)

      temp_local = 0.0_num
      drift_local = 0.0_num
      DO iz = -sf_order, sf_order
        DO iy = -sf_order, sf_order
          DO ix = -sf_order, sf_order
            temp_local = temp_local + gx(ix) * gy(iy) * gz(iz) &
                * temperature(cell_x+ix,cell_y+iy,cell_z+iz)
            drift_local = drift_local + gx(ix) * gy(iy) * gz(iz) &
                * drift(cell_x+ix,cell_y+iy,cell_z+iz)
          ENDDO
        ENDDO
      ENDDO

      IF (IAND(direction, c_dir_x) .NE. 0) current%part_p(1) = &
          momentum_from_temperature(mass, temp_local, idum) + drift_local

      IF (IAND(direction, c_dir_y) .NE. 0) current%part_p(2) = &
          momentum_from_temperature(mass, temp_local, idum) + drift_local

      IF (IAND(direction, c_dir_z) .NE. 0) current%part_p(3) = &
          momentum_from_temperature(mass, temp_local, idum) + drift_local

      current=>current%next
      ipart = ipart+1
    ENDDO

  END SUBROUTINE setup_particle_temperature



  SUBROUTINE setup_particle_density(density_in, part_family, min_density, &
      max_density, idum)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: density_in
    TYPE(particle_family), POINTER :: part_family
    REAL(num), INTENT(IN) :: min_density, max_density
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: weight_local
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cell_z_r, cell_frac_z
    TYPE(particle), POINTER :: current
    INTEGER :: cell_x, cell_y, cell_z
    INTEGER(KIND=8) :: ipart
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: weight_fn, temp
    REAL(num), DIMENSION(-2:2) :: gx, gy, gz
    REAL(num) :: data
    TYPE(particle_list), POINTER :: partlist
    INTEGER :: isubx, isuby, isubz, ierr
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: density
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: density_map

#ifdef PER_PARTICLE_WEIGHT
    ALLOCATE(density(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(density_map(-2:nx+3, -2:ny+3, -2:nz+3))
    density = density_in

    CALL field_bc(density)

    density_map = .FALSE.
    DO iz = -2, nz+3
      DO iy = -2, ny+3
        DO ix = -2, nx+3
          IF (density(ix, iy, iz) .GT. min_density) THEN
            density_map(ix, iy, iz) = .TRUE.
          ENDIF
          IF (density(ix, iy, iz) .GT. max_density &
              .AND. max_density .GT. 0.0_num) THEN
            density(ix, iy, iz) = max_density
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! Uniformly load particles in space
    CALL load_particles(part_family, density_map, idum)
    DEALLOCATE(density_map)

    ALLOCATE(weight_fn(-2:nx+3, -2:ny+3, -2:nz+3))
    ALLOCATE(temp(-2:nx+3, -2:ny+3, -2:nz+3))
    CALL MPI_BARRIER(comm, errcode)
    weight_fn = 0.0_num
    temp = 0.0_num

    partlist=>part_family%attached_list
    ! If using per particle weighing then use the weight function to match the
    ! uniform pseudoparticle density to the real particle density
    current=>partlist%head
    ipart = 0
    ! First loop converts number density into weight function
    DO WHILE(ipart .LT. partlist%count)
      IF (.NOT. ASSOCIATED(current)) PRINT *, "Bad Particle"
      cell_x_r = (current%part_pos(1)-x_min_local) / dx - 0.5_num
      cell_x = NINT(cell_x_r)
      cell_frac_x = REAL(cell_x, num) - cell_x_r
      cell_x = cell_x+1

      cell_y_r = (current%part_pos(2)-y_min_local) / dy - 0.5_num
      cell_y = NINT(cell_y_r)
      cell_frac_y = REAL(cell_y, num) - cell_y_r
      cell_y = cell_y+1

      cell_z_r = (current%part_pos(3)-z_min_local) / dz - 0.5_num
      cell_z = NINT(cell_z_r)
      cell_frac_z = REAL(cell_z, num) - cell_z_r
      cell_z = cell_z+1

      CALL particle_to_grid(cell_frac_x, gx)
      CALL particle_to_grid(cell_frac_y, gy)
      CALL particle_to_grid(cell_frac_z, gz)

      data = 1.0_num/(dx*dy*dz) ! Simply want to count particles per metre^2
      DO isubz = -sf_order, sf_order
        DO isuby = -sf_order, sf_order
          DO isubx = -sf_order, sf_order
            weight_fn(cell_x+isubx, cell_y+isuby, cell_z+isubz) = &
                weight_fn(cell_x+isubx, cell_y+isuby, cell_z+isubz) &
                + gx(isubx) * gy(isuby) * gz(isubz) * data
          ENDDO
        ENDDO
      ENDDO

      current=>current%next
      ipart = ipart+1
    ENDDO

    CALL processor_summation_bcs(weight_fn)
    IF (proc_x_min .EQ. MPI_PROC_NULL) weight_fn(0, :,:) = weight_fn(1,   :,:)
    IF (proc_x_max .EQ. MPI_PROC_NULL) weight_fn(nx,:,:) = weight_fn(nx-1,:,:)
    IF (proc_y_min .EQ. MPI_PROC_NULL) weight_fn(:,0, :) = weight_fn(:,1,   :)
    IF (proc_y_max .EQ. MPI_PROC_NULL) weight_fn(:,ny,:) = weight_fn(:,ny-1,:)
    IF (proc_z_min .EQ. MPI_PROC_NULL) weight_fn(:,:,0 ) = weight_fn(:,:,1   )
    IF (proc_z_max .EQ. MPI_PROC_NULL) weight_fn(:,:,nz) = weight_fn(:,:,nz-1)
    CALL field_zero_gradient(weight_fn, .TRUE.)

    DO iz = -2, nz+2
      DO iy = -2, ny+2
        DO ix = -2, nx+2
          IF (weight_fn(ix, iy, iz) .GT. 0.0_num) THEN
            weight_fn(ix, iy, iz) = density(ix, iy, iz)/weight_fn(ix, iy, iz)
          ELSE
            weight_fn(ix, iy, iz) = 0.0_num
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    IF (proc_x_min .EQ. MPI_PROC_NULL) weight_fn(0, :,:) = weight_fn(1,   :,:)
    IF (proc_x_max .EQ. MPI_PROC_NULL) weight_fn(nx,:,:) = weight_fn(nx-1,:,:)
    IF (proc_y_min .EQ. MPI_PROC_NULL) weight_fn(:,0, :) = weight_fn(:,1,   :)
    IF (proc_y_max .EQ. MPI_PROC_NULL) weight_fn(:,ny,:) = weight_fn(:,ny-1,:)
    IF (proc_z_min .EQ. MPI_PROC_NULL) weight_fn(:,:,0 ) = weight_fn(:,:,1   )
    IF (proc_z_max .EQ. MPI_PROC_NULL) weight_fn(:,:,nz) = weight_fn(:,:,nz-1)
    CALL field_zero_gradient(weight_fn, .TRUE.)

    partlist=>part_family%attached_list
    ! Second loop actually assigns weights to particles
    ! Again assumes linear interpolation
    current=>partlist%head
    ipart = 0
    DO WHILE(ipart .LT. partlist%count)
      cell_x_r = (current%part_pos(1)-x_min_local) / dx ! - 0.5_num
      cell_x = NINT(cell_x_r)
      cell_frac_x = REAL(cell_x, num) - cell_x_r
      cell_x = cell_x+1

      cell_y_r = (current%part_pos(2)-y_min_local) / dy ! - 0.5_num
      cell_y = NINT(cell_y_r)
      cell_frac_y = REAL(cell_y, num) - cell_y_r
      cell_y = cell_y+1

      cell_z_r = (current%part_pos(3)-z_min_local) / dz ! - 0.5_num
      cell_z = NINT(cell_z_r)
      cell_frac_z = REAL(cell_z, num) - cell_z_r
      cell_z = cell_z+1

      CALL grid_to_particle(cell_frac_x, gx)
      CALL grid_to_particle(cell_frac_y, gy)
      CALL grid_to_particle(cell_frac_z, gz)

      weight_local = 0.0_num
      DO isubz = -sf_order, +sf_order
        DO isuby = -sf_order, +sf_order
          DO isubx = -sf_order, +sf_order
            weight_local = &
                weight_local + gx(isubx)*gy(isuby)*gz(isubz) &
                * weight_fn(cell_x+isubx, cell_y+isuby, cell_z+isubz)
          ENDDO
        ENDDO
      ENDDO
      current%weight = weight_local
      current=>current%next
      ipart = ipart+1
    ENDDO

    DEALLOCATE(weight_fn)
    DEALLOCATE(density)
#else
    IF (rank .EQ. 0) &
        PRINT *, "Autoloader only available when using per particle weighting"
    CALL MPI_ABORT(comm, errcode, ierr)
#endif

  END SUBROUTINE setup_particle_density



  FUNCTION momentum_from_temperature(mass, temperature, idum)

    REAL(num), INTENT(IN) :: mass, temperature
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: momentum_from_temperature

    REAL(num) :: stdev
    REAL(num) :: rand1, rand2, w

    ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers
    ! The standard deviation (stdev) is related to temperature

    stdev = SQRT(2.0_num*temperature*kb*mass)

    DO
      rand1 = random(idum)
      rand2 = random(idum)

      rand1 = 2.0_num*rand1 - 1.0_num
      rand2 = 2.0_num*rand2 - 1.0_num

      w = rand1**2 + rand2**2

      IF (w .LT. 1.0_num) EXIT
    ENDDO

    w = SQRT((-2.0_num * LOG(w) )/w)

    momentum_from_temperature = rand1 * w * stdev

  END FUNCTION momentum_from_temperature



  FUNCTION sample_dist_function(axis, dist_fn, idum)

    REAL(num), DIMENSION(:), INTENT(IN) :: axis, dist_fn
    INTEGER, INTENT(INOUT) :: idum
    REAL(num), DIMENSION(:), ALLOCATABLE :: cdf
    REAL(num) :: position, d_cdf
    INTEGER :: n_points, ipoint, start, endpoint, current
    REAL(num) :: sample_dist_function

    n_points = SIZE(dist_fn)
    ALLOCATE(cdf(1:n_points))
    DO ipoint = 1, n_points
      cdf(ipoint) = SUM(dist_fn(1:ipoint))
    ENDDO
    cdf = cdf/SUM(dist_fn)

    position = random(idum)
    sample_dist_function = 0.0_num

    start = 1
    endpoint = n_points
    current = (start+endpoint)/2

    DO current = 1, n_points-1
      IF (cdf(current) .LE. position .AND. cdf(current+1) .GE. position) THEN
        d_cdf = cdf(current+1)-cdf(current)
        sample_dist_function = &
            (axis(current)*(position-cdf(current))/d_cdf &
            + axis(current+1)*(cdf(current+1)-position)/d_cdf)
        EXIT
      ENDIF
    ENDDO

    DEALLOCATE(cdf)

  END FUNCTION sample_dist_function



  FUNCTION random(idum)

    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: random
    INTEGER, PARAMETER :: ia = 16807, im = 2147483647, iq = 127773
    INTEGER, PARAMETER :: ir = 2836, mask = 123459876
    REAL(dbl), PARAMETER :: am = 1.0_dbl/2147483647.0_dbl

    INTEGER :: k

    idum = IEOR(idum, mask)
    k = idum/iq

    idum = ia*(idum-k*iq)-ir*k
    IF (idum .LT. 0) THEN
      idum = idum+im
    ENDIF

    random = am*idum
    idum = IEOR(idum, mask)

    IF (random .GT. max_rand) max_rand = random

  END FUNCTION random

END MODULE helper
