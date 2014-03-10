MODULE setup

  USE encoded_source
  USE fields
  USE mpi_subtype_control
  USE version_data
  USE welcome
  USE split_particle
  USE shunt
  USE laser
  USE window
  USE timer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control, minimal_init, restart_data
  PUBLIC :: open_files, close_files, flush_stat_file
  PUBLIC :: setup_species, after_deck_last, set_dt
  PUBLIC :: read_cpu_split

  TYPE(particle), POINTER, SAVE :: iterator_list
  CHARACTER(LEN=11+data_dir_max_length), SAVE :: stat_file

CONTAINS

  SUBROUTINE minimal_init

    INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
    INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
    INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(r=3000)
    INTEGER :: ierr

    IF (num .EQ. r4) THEN
      ! Should use MPI_SIZEOF() but this breaks on scalimpi
      realsize = 4
      mpireal = MPI_REAL4
    ELSE IF (num .EQ. r8) THEN
      realsize = 8
      mpireal = MPI_REAL8
    ELSE
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Cannot determine size of real'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    ENDIF

    dt_plasma_frequency = 0.0_num
    dt_multiplier = 0.95_num
    stdout_frequency = 0
    cpml_thickness = 6
    cpml_kappa_max = 20.0_num
    cpml_a_max = 0.15_num
    cpml_sigma_max = 0.7_num
    cpml_x_min_offset = 0
    cpml_x_max_offset = 0

    window_shift = 0.0_num
    npart_global = -1
    smooth_currents = .FALSE.
    use_balance = .FALSE.
    use_random_seed = .FALSE.
    use_offset_grid = .FALSE.
    use_particle_lists = .FALSE.
    use_multiphoton = .TRUE.
    use_bsi = .TRUE.
    need_random_state = .FALSE.
    force_first_to_be_restartable = .FALSE.
    force_final_to_be_restartable = .FALSE.
    full_dump_every = -1
    restart_dump_every = -1
    nsteps = -1
    t_end = HUGE(1.0_num)
    particles_max_id = 0

    NULLIFY(laser_x_min)
    NULLIFY(laser_x_max)

    NULLIFY(dist_fns)
    NULLIFY(io_block_list)

    run_date = get_unix_time()

    CALL set_field_order(2)

    CALL timer_init

    ! This array is true if a field component is staggered in the
    ! given direction.
    stagger = .FALSE.
    stagger(c_dir_x,c_stagger_ex) = .TRUE.
    stagger(c_dir_y,c_stagger_ey) = .TRUE.
    stagger(c_dir_z,c_stagger_ez) = .TRUE.

    stagger(c_dir_x:c_dir_z,c_stagger_bx) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_by) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_bz) = .TRUE.
    stagger(c_dir_x,c_stagger_bx) = .FALSE.
    stagger(c_dir_y,c_stagger_by) = .FALSE.
    stagger(c_dir_z,c_stagger_bz) = .FALSE.

    ALLOCATE(x(1))
    x = 0.0_num

  END SUBROUTINE minimal_init



  SUBROUTINE after_control

    CALL setup_grid
    CALL set_initial_values

  END SUBROUTINE after_control



  SUBROUTINE setup_grid

    INTEGER :: iproc, ix
    REAL(num) :: xb_min

    length_x = x_max - x_min
    dx = length_x / REAL(nx_global-2*cpml_thickness, num)
    x_grid_min = x_min - dx * cpml_thickness
    x_grid_max = x_max + dx * cpml_thickness

    ! Shift grid to cell centres.
    ! At some point the grid may be redefined to be node centred.

    xb_min = x_grid_min
    x_grid_min = x_grid_min + dx / 2.0_num
    x_grid_max = x_grid_max - dx / 2.0_num

    ! Setup global grid
    DO ix = -2, nx_global + 3
      x_global(ix) = x_grid_min + (ix - 1) * dx
    ENDDO
    DO ix = 1, nx_global + 1
      xb_global(ix) = xb_min + (ix - 1) * dx
      xb_offset_global(ix) = xb_global(ix)
    ENDDO

    DO iproc = 0, nprocx-1
      x_grid_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO

    x_grid_min_local = x_grid_mins(x_coords)
    x_grid_max_local = x_grid_maxs(x_coords)

    x_min_local = x_grid_min_local + (cpml_x_min_offset - 0.5_num) * dx
    x_max_local = x_grid_max_local - (cpml_x_max_offset - 0.5_num) * dx

    ! Setup local grid
    DO ix = -2, nx + 3
      x(ix) = x_global(nx_global_min+ix-1)
    ENDDO

  END SUBROUTINE setup_grid



  SUBROUTINE after_deck_last

    INTEGER :: i

    CALL setup_data_averaging
    CALL setup_split_particles
    CALL setup_field_boundaries

    IF (cpml_boundaries) THEN
      CALL allocate_cpml_fields
      CALL set_cpml_helpers(nx, nx_global_min, nx_global_max)
    ELSE
      cpml_thickness = 0
      cpml_kappa_max = 1.0_num
      cpml_a_max = 0.0_num
      cpml_sigma_max = 0.0_num
      DO i = 1, n_io_blocks
        io_block_list(i)%dumpmask(c_dump_cpml_psi_eyx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_ezx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_byx) = c_io_never
        io_block_list(i)%dumpmask(c_dump_cpml_psi_bzx) = c_io_never
      ENDDO
    ENDIF

  END SUBROUTINE after_deck_last



  SUBROUTINE setup_data_averaging()

    INTEGER :: io, ib, nspec_local, nstep_average, mask
    REAL(num) :: dt_average
    TYPE(averaged_data_block), POINTER :: avg

    IF (.NOT. any_average) RETURN

    DO io = 1, n_io_blocks
      IF (io_block_list(io)%any_average) THEN
        dt_min_average = t_end
      ELSE
        dt_min_average = -1.0_num
      ENDIF
    ENDDO

    DO io = 1, num_vars_to_dump
      ib = averaged_var_block(io)
      IF (ib .NE. 0) THEN
        dt_average  = io_block_list(ib)%dt_average
        nstep_average = io_block_list(ib)%nstep_average

        IF (nstep_average .GT. 0 .AND. dt_average .GT. 0) THEN
          io_block_list(ib)%dt_min_average = &
              MIN(io_block_list(ib)%dt_min_average, &
              dt_average / REAL(nstep_average, num))
        ENDIF

        mask = io_block_list(ib)%dumpmask(io)
        nspec_local = 0
        IF (IAND(mask, c_io_no_sum) .EQ. 0) &
            nspec_local = 1
        IF (IAND(mask, c_io_species) .NE. 0) &
            nspec_local = nspec_local + n_species

        IF (nspec_local .LE. 0) CYCLE

        avg => io_block_list(ib)%averaged_data(io)
        IF (avg%dump_single) THEN
          ALLOCATE(avg%r4array(-2:nx+3, nspec_local))
          avg%r4array = 0.0_num
        ELSE
          ALLOCATE(avg%array(-2:nx+3, nspec_local))
          avg%array = 0.0_num
        ENDIF

        avg%species_sum = 0
        avg%n_species = 0
        IF (IAND(mask, c_io_no_sum) .EQ. 0) &
            avg%species_sum = 1
        IF (IAND(mask, c_io_species) .NE. 0) &
            avg%n_species = n_species
        avg%real_time = 0.0_num
        avg%started = .FALSE.
      ENDIF
    ENDDO

  END SUBROUTINE setup_data_averaging



  SUBROUTINE setup_species

    INTEGER :: ispecies

    ALLOCATE(species_list(n_species))
    ALLOCATE(io_list_data(n_species))
    io_list => species_list

    DO ispecies = 1, n_species
      species_list(ispecies)%name = blank
      species_list(ispecies)%mass = -1.0_num
      species_list(ispecies)%charge = 0.0_num
      species_list(ispecies)%weight = 1.0_num
      species_list(ispecies)%dumpmask = c_io_always
      species_list(ispecies)%count = -1
      species_list(ispecies)%id = 0
      species_list(ispecies)%npart_per_cell = -1
      species_list(ispecies)%split = .FALSE.
      species_list(ispecies)%npart_max = 0
      species_list(ispecies)%global_count = 0
      species_list(ispecies)%species_type = c_species_id_generic
      species_list(ispecies)%immobile = .FALSE.
      NULLIFY(species_list(ispecies)%next)
      NULLIFY(species_list(ispecies)%prev)
      NULLIFY(species_list(ispecies)%ext_temp_x_min)
      NULLIFY(species_list(ispecies)%ext_temp_x_max)
      NULLIFY(species_list(ispecies)%secondary_list)
    ENDDO

    DO ispecies = 1, n_species
      CALL initialise_stack(species_list(ispecies)%density_function)
      CALL initialise_stack(species_list(ispecies)%temperature_function(1))
      CALL initialise_stack(species_list(ispecies)%temperature_function(2))
      CALL initialise_stack(species_list(ispecies)%temperature_function(3))
      CALL initialise_stack(species_list(ispecies)%drift_function(1))
      CALL initialise_stack(species_list(ispecies)%drift_function(2))
      CALL initialise_stack(species_list(ispecies)%drift_function(3))
      CALL set_stack_zero(species_list(ispecies)%density_function)
      CALL set_stack_zero(species_list(ispecies)%temperature_function(1))
      CALL set_stack_zero(species_list(ispecies)%temperature_function(2))
      CALL set_stack_zero(species_list(ispecies)%temperature_function(3))
      CALL set_stack_zero(species_list(ispecies)%drift_function(1))
      CALL set_stack_zero(species_list(ispecies)%drift_function(2))
      CALL set_stack_zero(species_list(ispecies)%drift_function(3))
      species_list(ispecies)%electron = .FALSE.
      species_list(ispecies)%ionise = .FALSE.
      species_list(ispecies)%ionise_to_species = -1
      species_list(ispecies)%release_species = -1
      species_list(ispecies)%n = 0
      species_list(ispecies)%l = 0
      species_list(ispecies)%ionisation_energy = HUGE(0.0_num)
      species_list(ispecies)%migrate%this_species = .FALSE.
      species_list(ispecies)%migrate%fluid = .FALSE.
      species_list(ispecies)%migrate%done = .FALSE.
      species_list(ispecies)%migrate%promoteable = .FALSE.
      species_list(ispecies)%migrate%demoteable = .FALSE.
      species_list(ispecies)%migrate%promote_to_species = 0
      species_list(ispecies)%migrate%demote_to_species = 0
      species_list(ispecies)%migrate%promotion_energy_factor = 1.0_num
      species_list(ispecies)%migrate%demotion_energy_factor = 1.0_num
      species_list(ispecies)%migrate%promotion_density = HUGE(1.0_num)
      species_list(ispecies)%migrate%demotion_density = 0.0_num
#ifdef TRACER_PARTICLES
      species_list(ispecies)%tracer = .FALSE.
#endif
#ifdef PARTICLE_PROBES
      NULLIFY(species_list(ispecies)%attached_probes)
#endif
    ENDDO

  END SUBROUTINE setup_species



  SUBROUTINE setup_field_boundaries

    INTEGER :: nx0, nx1

    nx0 = 1
    nx1 = nx
    IF (bc_field(c_bd_x_min) .EQ. c_bc_cpml_laser) nx0 = cpml_x_min_laser_idx-1
    IF (bc_field(c_bd_x_max) .EQ. c_bc_cpml_laser) nx1 = cpml_x_max_laser_idx+1

    ex_x_min = 0.5_num * (ex(nx0) + ex(nx0-1))
    ey_x_min = ey(nx0)
    ez_x_min = ez(nx0)
    ex_x_max = 0.5_num * (ex(nx1) + ex(nx1-1))
    ey_x_max = ey(nx1)
    ez_x_max = ez(nx1)

    bx_x_min = bx(nx0)
    by_x_min = 0.5_num * (by(nx0) + by(nx0-1))
    bz_x_min = 0.5_num * (bz(nx0) + bz(nx0-1))
    bx_x_max = bx(nx1)
    by_x_max = 0.5_num * (by(nx1) + by(nx1-1))
    bz_x_max = 0.5_num * (bz(nx1) + bz(nx1-1))

  END SUBROUTINE setup_field_boundaries



  SUBROUTINE open_files

#ifdef NO_IO
    RETURN
#else
    INTEGER :: errcode, ierr
    LOGICAL :: exists

    IF (rank .EQ. 0) THEN
      WRITE(stat_file, '(a, ''/epoch1d.dat'')') TRIM(data_dir)
      IF (ic_from_restart) THEN
        INQUIRE(file=stat_file, exist=exists)
        IF (exists) THEN
          OPEN(unit=stat_unit, status='OLD', position='APPEND', &
              file=stat_file, iostat=errcode)
        ELSE
          OPEN(unit=stat_unit, status='NEW', file=stat_file, iostat=errcode)
        ENDIF
      ELSE
        OPEN(unit=stat_unit, status='REPLACE', file=stat_file, iostat=errcode)
      ENDIF
      IF (errcode .NE. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Cannot create "epoch1d.dat" output file. The most common ' &
            // 'cause of this problem '
        PRINT*, 'is that the ouput directory does not exist'
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        STOP
      ENDIF
      IF (ic_from_restart) THEN
        WRITE(stat_unit,*)
        WRITE(stat_unit,*) 'Restarting from ', TRIM(restart_filename)
        WRITE(stat_unit,*) ascii_header
      ELSE
        WRITE(stat_unit,*) ascii_header
        WRITE(stat_unit,*)
      ENDIF
    ENDIF
#endif

  END SUBROUTINE open_files



  SUBROUTINE flush_stat_file

#ifdef NO_IO
    RETURN
#else
    INTEGER :: errcode

    IF (rank .EQ. 0) THEN
      CLOSE(unit=stat_unit)
      OPEN(unit=stat_unit, status='OLD', position='APPEND', &
          file=stat_file, iostat=errcode)
    ENDIF
#endif

  END SUBROUTINE flush_stat_file



  SUBROUTINE close_files

#ifdef NO_IO
    RETURN
#else
    IF (rank .EQ. 0) CLOSE(unit=stat_unit)
#endif

  END SUBROUTINE close_files



  SUBROUTINE set_initial_values

    INTEGER :: seed

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    ! Set up random number seed
    seed = 7842432
    IF (use_random_seed) CALL SYSTEM_CLOCK(seed)
    seed = seed + rank

    CALL random_init(seed)

  END SUBROUTINE set_initial_values



  SUBROUTINE set_plasma_frequency_dt

    INTEGER :: ispecies, ix
    REAL(num) :: min_dt, omega2, omega, k_max

    IF (ic_from_restart) RETURN

    min_dt = 1000000.0_num
    k_max = 2.0_num * pi / dx

    ! Identify the plasma frequency (Bohm-Gross dispersion relation)
    ! Note that this doesn't get strongly relativistic plasmas right
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type .NE. c_species_id_photon) THEN
        DO ix = 1, nx
          omega2 = initial_conditions(ispecies)%density(ix) * q0**2 &
              / species_list(ispecies)%mass / epsilon0 &
              + 3.0_num * k_max**2 * kb &
              * MAXVAL(initial_conditions(ispecies)%temp(ix,:)) &
              / species_list(ispecies)%mass
          IF (omega2 .LE. c_tiny) CYCLE
          omega = SQRT(omega2)
          IF (2.0_num * pi / omega .LT. min_dt) min_dt = 2.0_num * pi / omega
        ENDDO ! ix
      ENDIF
    ENDDO

    CALL MPI_ALLREDUCE(min_dt, dt_plasma_frequency, 1, mpireal, MPI_MIN, &
        comm, errcode)
    ! Must resolve plasma frequency
    dt_plasma_frequency = dt_plasma_frequency / 2.0_num

  END SUBROUTINE set_plasma_frequency_dt



  SUBROUTINE set_dt        ! sets CFL limited step

    INTEGER :: io

    CALL set_plasma_frequency_dt
    CALL set_laser_dt

    dt = cfl * dx / c
    IF (dt_plasma_frequency .GT. c_tiny) dt = MIN(dt, dt_plasma_frequency)
    IF (dt_laser .GT. c_tiny) dt = MIN(dt, dt_laser)
    dt = dt_multiplier * dt

    IF (.NOT. any_average) RETURN

    DO io = 1, n_io_blocks
      IF (.NOT. io_block_list(io)%any_average) CYCLE

      io_block_list(io)%average_time = MAX(io_block_list(io)%dt_average, &
          dt * io_block_list(io)%nstep_average)

      IF (io_block_list(io)%dt_min_average .GT. 0 &
          .AND. io_block_list(io)%dt_min_average .LT. dt) THEN
        IF (rank .EQ. 0) THEN
          PRINT*,'*** WARNING ***'
          PRINT*,'Time step is too small to satisfy "nstep_average"'
          PRINT*,'Averaging will occur over fewer time steps than specified'
          PRINT*,'Set "dt_multiplier" less than ', &
              dt_multiplier * io_block_list(io)%dt_min_average / dt, &
              ' to fix this'
        ENDIF
        io_block_list(io)%dt_min_average = -1
      ENDIF
    ENDDO

  END SUBROUTINE set_dt



  SUBROUTINE find_species_by_blockid(specname, species_number)

    CHARACTER(LEN=*), INTENT(IN) :: specname
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: ispecies, i1, i2

    CALL strip_species_blockid(specname, i1, i2)

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(specname(i1:i2), species_list(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE find_species_by_blockid



  SUBROUTINE find_species_by_id(species_id, species_number)

    CHARACTER(LEN=*), INTENT(IN) :: species_id
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: ispecies

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(species_id, species_list(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE find_species_by_id



  SUBROUTINE copy_string(str_in, str_out)

    CHARACTER(LEN=*), INTENT(IN) :: str_in
    CHARACTER(LEN=*), INTENT(INOUT) :: str_out
    INTEGER :: len1, len2, olen, i

    len1 = LEN_TRIM(str_in)
    len2 = LEN(str_out)
    olen = MIN(len1,len2)
    IF (olen .GT. 0) THEN
      str_out(1:olen) = str_in(1:olen)
      DO i = olen+1,len2
        str_out(i:i) = ' '
      ENDDO
    ELSE
      DO i = 1,len2
        str_out(i:i) = ' '
      ENDDO
    ENDIF

  END SUBROUTINE copy_string



  SUBROUTINE find_species_by_id_or_blockid(species_id, block_id, species_number)

    CHARACTER(LEN=*), INTENT(INOUT) :: species_id
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: i1, i2

    IF (str_cmp(species_id, '__unknown__')) THEN
      CALL strip_species_blockid(block_id, i1, i2)
      IF (i2 .LE. LEN_TRIM(species_id)) THEN
        CALL copy_string(block_id(i1:i2), species_id)
      ELSE
        CALL copy_string('__unknown__', species_id)
      ENDIF
    ENDIF

    CALL find_species_by_id(species_id, species_number)

  END SUBROUTINE find_species_by_id_or_blockid



  SUBROUTINE strip_species_blockid(name, i1, i2)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: i1, i2
    INTEGER :: ii

    i1 = 2
    i2 = LEN_TRIM(name)
    DO ii = 1, i2
      IF (name(ii:ii) .EQ. '/') RETURN
      i1 = i1 + 1
    ENDDO
    i1 = 1

  END SUBROUTINE strip_species_blockid



  SUBROUTINE restart_data(step)

    INTEGER, INTENT(OUT) :: step
    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, str1
    CHARACTER(LEN=c_id_length) :: species_id
    CHARACTER(LEN=c_max_string_length) :: name, len_string
    INTEGER :: blocktype, datatype, code_io_version, string_len, ispecies
    INTEGER :: ierr, i, is, iblock, nblocks, ndims, geometry
    INTEGER(i8) :: npart, npart_local
    INTEGER(i8), ALLOCATABLE :: nparts(:), npart_locals(:), npart_proc(:)
    INTEGER, DIMENSION(4) :: dims
    INTEGER, ALLOCATABLE :: random_states_per_proc(:)
    REAL(num), DIMENSION(2*c_ndims) :: extents
    LOGICAL :: restart_flag, got_full
    LOGICAL, ALLOCATABLE :: species_found(:)
    TYPE(sdf_file_handle) :: sdf_handle
    TYPE(particle_species), POINTER :: species
    INTEGER, POINTER :: species_subtypes(:)

    got_full = .FALSE.
    npart_global = 0
    step = -1

    IF (rank .EQ. 0) THEN
      PRINT*,'Attempting to restart from file: ',TRIM(full_restart_filename)
    ENDIF

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    ! Reset io_block parameters
    DO i = 1, n_io_blocks
      IF (io_block_list(i)%dt_snapshot .GT. 0.0_num) THEN
        io_block_list(i)%time_prev = time
      ELSE
        io_block_list(i)%time_prev = 0.0_num
      ENDIF
      IF (io_block_list(i)%nstep_snapshot .GT. 0) THEN
        io_block_list(i)%nstep_prev = step
      ELSE
        io_block_list(i)%nstep_prev = 0
      ENDIF
      IF (ASSOCIATED(io_block_list(i)%dump_at_nsteps)) THEN
        DO is = 1, SIZE(io_block_list(i)%dump_at_nsteps)
          IF (step .GE. io_block_list(i)%dump_at_nsteps(is)) THEN
            io_block_list(i)%dump_at_nsteps(is) = HUGE(1)
          ENDIF
        ENDDO
      ENDIF
      IF (ASSOCIATED(io_block_list(i)%dump_at_times)) THEN
        DO is = 1, SIZE(io_block_list(i)%dump_at_times)
          IF (time .GE. io_block_list(i)%dump_at_times(is)) THEN
            io_block_list(i)%dump_at_times(is) = HUGE(1.0_num)
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    IF (.NOT. restart_flag) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file is not a restart dump. Unable to continue.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    ENDIF

    IF (.NOT.str_cmp(code_name, 'Epoch1d')) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Epoch1d. Unable to ', &
            'continue.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    ENDIF

    IF (string_len > c_max_string_length) THEN
      IF (rank .EQ. 0) THEN
        CALL integer_as_string(string_len, len_string)
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file string lengths are too large to read.'
        PRINT*, 'Please increase the size of "c_max_string_length" in ', &
            'shared_data.F90 to ','be at least ',TRIM(len_string)
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    ENDIF

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    IF (rank .EQ. 0) THEN
      PRINT*, 'Loading snapshot for time', time
      CALL create_ascii_header
#ifndef NO_IO
      WRITE(stat_unit,*) ascii_header
      WRITE(stat_unit,*)
#endif
    ENDIF

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    dt_from_restart = 0.0_num

    IF (rank .EQ. 0) PRINT*, 'Input file contains', nblocks, 'blocks'

    CALL sdf_read_blocklist(sdf_handle)

    ALLOCATE(nparts(n_species))
    ALLOCATE(npart_locals(n_species))
    ALLOCATE(species_found(n_species))
    nparts = 0
    npart_locals = 0
    species_found = .FALSE.

    ! Scan file for particle species and allocate storage
    ! Both particle grids and CPU split blocks are interrogated. CPU split
    ! blocks take precedence if found and "use_exact_restart" is set.
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype .EQ. c_blocktype_cpu_split) THEN
        IF (.NOT.use_exact_restart .OR. datatype .NE. c_datatype_integer8) CYCLE
      ELSE IF (blocktype .NE. c_blocktype_point_mesh) THEN
        CYCLE
      ENDIF

      CALL sdf_read_point_mesh_info(sdf_handle, npart, geometry, species_id)
      CALL find_species_by_id_or_blockid(species_id, block_id, ispecies)
      IF (ispecies .EQ. 0) THEN
        IF (rank .EQ. 0) THEN
          IF (.NOT. str_cmp(species_id(1:6), 'subset')) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Particle species "', TRIM(species_id), &
                '" from restart dump ', 'not found in input deck. Ignoring.'
          ENDIF
        ENDIF
        CYCLE
      ENDIF

      IF (species_found(ispecies)) CYCLE

      species => species_list(ispecies)

      IF (blocktype .EQ. c_blocktype_cpu_split) THEN
        CALL sdf_read_cpu_split_info(sdf_handle, dims, geometry)

        ALLOCATE(npart_proc(dims(1)))
        CALL sdf_read_srl_cpu_split(sdf_handle, npart_proc)

        npart = 0
        DO i = 1,dims(1)
          npart = npart + npart_proc(i)
        ENDDO
        npart_local = npart_proc(rank+1)
        DEALLOCATE(npart_proc)

        npart_locals(ispecies) = npart_local
        nparts(ispecies) = npart
        species_found(ispecies) = .TRUE.
      ELSE IF (blocktype .EQ. c_blocktype_point_mesh) THEN
        CALL sdf_read_point_mesh_info(sdf_handle, npart)

        nparts(ispecies) = npart

        npart_local = npart / nproc
        IF (npart_local * nproc .NE. npart) THEN
          IF (rank .LT. npart - npart_local * nproc) &
              npart_local = npart_local + 1
        ENDIF

        npart_locals(ispecies) = npart_local
      ENDIF
    ENDDO

    DEALLOCATE(species_found)

    ! Do the species allocation
    DO ispecies = 1, n_species
      species => species_list(ispecies)
      npart_local = npart_locals(ispecies)

      CALL create_allocated_partlist(species%attached_list, npart_local)

      npart_global = npart_global + nparts(ispecies)
      species%count = nparts(ispecies)
    ENDDO

    DEALLOCATE(nparts, npart_locals)

    CALL create_subtypes_for_load(species_subtypes)

    CALL sdf_seek_start(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_array)
        IF (use_exact_restart .AND. need_random_state &
            .AND. str_cmp(block_id, 'random_states')) THEN
          ALLOCATE(random_states_per_proc(4*nproc))
          CALL sdf_read_srl(sdf_handle, random_states_per_proc)
          CALL set_random_state(random_states_per_proc(4*rank+1:4*(rank+1)))
          DEALLOCATE(random_states_per_proc)
        ELSE IF (str_cmp(block_id, 'file_numbers')) THEN
          CALL sdf_read_array_info(sdf_handle, dims)
          IF (ndims .NE. 1 .OR. dims(1) .NE. SIZE(file_numbers)) THEN
            IF (rank .EQ. 0) THEN
              PRINT*, '*** WARNING ***'
              PRINT*, 'Output file numbers do not agree. Ignoring.'
            ENDIF
          ELSE
            CALL sdf_read_srl(sdf_handle, file_numbers)
          ENDIF
        ENDIF
      CASE(c_blocktype_constant)
        IF (str_cmp(block_id, 'dt_plasma_frequency')) THEN
          CALL sdf_read_srl(sdf_handle, dt_plasma_frequency)
        ELSE IF (str_cmp(block_id, 'dt')) THEN
          CALL sdf_read_srl(sdf_handle, dt_from_restart)
        ELSE IF (str_cmp(block_id, 'window_shift_fraction')) THEN
          CALL sdf_read_srl(sdf_handle, window_shift_fraction)
        ELSE IF (block_id(1:7) .EQ. 'weight/') THEN
          CALL find_species_by_blockid(block_id, ispecies)
          IF (ispecies .EQ. 0) CYClE
          CALL sdf_read_srl(sdf_handle, species_list(ispecies)%weight)
        ELSE IF (block_id(1:5) .EQ. 'nppc/') THEN
          CALL find_species_by_blockid(block_id, ispecies)
          IF (ispecies .EQ. 0) CYClE
          CALL sdf_read_srl(sdf_handle, species_list(ispecies)%npart_per_cell)
        ELSE IF (block_id(1:10) .EQ. 'time_prev/') THEN
          DO i = 1, n_io_blocks
            IF (str_cmp(block_id(11:), io_block_list(i)%name)) THEN
              CALL sdf_read_srl(sdf_handle, io_block_list(i)%time_prev)
              EXIT
            ENDIF
          ENDDO
        ELSE IF (block_id(1:11) .EQ. 'nstep_prev/') THEN
          DO i = 1, n_io_blocks
            IF (str_cmp(block_id(12:), io_block_list(i)%name)) THEN
              CALL sdf_read_srl(sdf_handle, io_block_list(i)%nstep_prev)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      CASE(c_blocktype_plain_mesh)
        IF (.NOT.got_full) THEN
          IF (str_cmp(block_id, 'grid') &
              .OR. str_cmp(block_id, 'grid_full')) THEN
            CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims, extents)
            x_min = extents(1)
            x_max = extents(c_ndims+1)
            IF (str_cmp(block_id, 'grid_full')) got_full = .TRUE.
          ENDIF
        ENDIF
      CASE(c_blocktype_point_mesh)
        CALL sdf_read_point_mesh_info(sdf_handle, npart, geometry, species_id)

        CALL find_species_by_id_or_blockid(species_id, block_id, ispecies)
        IF (ispecies .EQ. 0) CYClE
        species => species_list(ispecies)

        npart_local = species%attached_list%count
        iterator_list => species%attached_list%head

        CALL sdf_read_point_mesh(sdf_handle, npart_local, &
            species_subtypes(ispecies), it_part)

      CASE(c_blocktype_plain_variable)
        CALL sdf_read_plain_variable_info(sdf_handle, dims, str1, mesh_id)

        IF (.NOT.str_cmp(mesh_id, 'grid')) CYCLE

        IF (dims(1) .NE. nx_global) THEN
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Number of gridpoints in restart dump does not match', &
                ' the input deck.'
            CALL integer_as_string(nx_global, str1)
            PRINT*, 'Input deck grid: ', TRIM(str1)
            CALL integer_as_string(dims(1), str1)
            PRINT*, 'Restart dump grid: ', TRIM(str1)
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
          STOP
        ENDIF

        IF (str_cmp(block_id, 'ex')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ex, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'ey')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ey, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'ez')) THEN
          CALL sdf_read_plain_variable(sdf_handle, ez, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'bx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, bx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'by')) THEN
          CALL sdf_read_plain_variable(sdf_handle, by, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'bz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, bz, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'jx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jx, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'jy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jy, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'jz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jz, &
              subtype_field, subarray_field_big)

        ELSE IF (str_cmp(block_id, 'cpml_psi_eyx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_eyx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_ezx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_ezx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_byx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_byx, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'cpml_psi_bzx')) THEN
          CALL sdf_read_plain_variable(sdf_handle, cpml_psi_bzx, &
              subtype_field, subarray_field)

        ENDIF

      CASE(c_blocktype_point_variable)
        CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id, &
            str1, species_id)

        CALL find_species_by_id_or_blockid(species_id, mesh_id, ispecies)
        IF (ispecies .EQ. 0) CYClE
        species => species_list(ispecies)

        IF (npart .NE. species%count) THEN
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Malformed restart dump. Number of particle variables', &
                ' does not match grid.'
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
          STOP
        ENDIF

        iterator_list => species%attached_list%head
        npart_local = species%attached_list%count

        IF (block_id(1:3) .EQ. 'px/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_px)

        ELSE IF (block_id(1:3) .EQ. 'py/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_py)

        ELSE IF (block_id(1:3) .EQ. 'pz/') THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_pz)

        ELSE IF (block_id(1:3) .EQ. 'id/') THEN
#if PARTICLE_ID || PARTICLE_ID4
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_id)
#else
          IF (rank .EQ. 0) THEN
            PRINT*, '*** WARNING ***'
            PRINT*, 'Discarding particle IDs.'
            PRINT*, 'To use, please recompile with the -DPARTICLE_ID option.'
          ENDIF
#endif

        ELSE IF (block_id(1:7) .EQ. 'weight/') THEN
#ifdef PER_PARTICLE_WEIGHT
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_weight)
#else
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with per particle weight.'
            PRINT*, 'Please recompile with the -DPER_PARTICLE_WEIGHT option.'
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
          STOP
#endif

        ELSE IF (block_id(1:14) .EQ. 'optical depth/') THEN
#ifdef PHOTONS
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_optical_depth)
#else
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with optical depths.'
            PRINT*, 'Please recompile with the -DPHOTONS option.'
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
          STOP
#endif

        ELSE IF (block_id(1:11) .EQ. 'qed energy/') THEN
#ifdef PHOTONS
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_qed_energy)
#else
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with QED energies.'
            PRINT*, 'Please recompile with the -DPHOTONS option.'
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
          STOP
#endif

        ELSE IF (block_id(1:14) .EQ. 'trident depth/') THEN
#ifdef TRIDENT_PHOTONS
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              species_subtypes(ispecies), it_optical_depth_trident)
#else
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Cannot load dump file with Trident optical depths.'
            PRINT*, 'Please recompile with the -DTRIDENT_PHOTONS option.'
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
          STOP
#endif
        ENDIF
      END SELECT
    ENDDO

    CALL sdf_close(sdf_handle)
    CALL free_subtypes_for_load(species_subtypes)

    CALL setup_grid

  END SUBROUTINE restart_data



  SUBROUTINE read_cpu_split

    CHARACTER(LEN=c_id_length) :: code_name, block_id
    CHARACTER(LEN=c_max_string_length) :: name
    INTEGER :: ierr, step, code_io_version, string_len, nblocks, ndims
    INTEGER :: blocktype, datatype, geometry, iblock, npx
    INTEGER, DIMENSION(4) :: dims
    LOGICAL :: restart_flag
    TYPE(sdf_file_handle) :: sdf_handle

    CALL sdf_open(sdf_handle, full_restart_filename, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    IF (.NOT. restart_flag) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file is not a restart dump. Unable to continue.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    ENDIF

    IF (.NOT.str_cmp(code_name, 'Epoch1d')) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Epoch1d. Unable to ', &
            'continue.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    ENDIF

    nblocks = sdf_read_nblocks(sdf_handle)

    CALL sdf_read_blocklist(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype .EQ. c_blocktype_cpu_split &
          .AND. datatype .NE. c_datatype_integer8) THEN
        CALL sdf_read_cpu_split_info(sdf_handle, dims, geometry)
        npx = dims(1) + 1
        IF (npx .EQ. nproc) THEN
          nprocx = npx
          ALLOCATE(old_x_max(nprocx))
          CALL sdf_read_srl_cpu_split(sdf_handle, old_x_max)
        ELSE
          IF (rank .EQ. 0) THEN
            PRINT*, '*** WARNING ***'
            PRINT'('' SDF restart file was generated using'', &
                & i4,'' CPUs.'')', npx
            PRINT*, 'Ignoring "use_exact_restart" flag.'
          ENDIF
          use_exact_restart = .FALSE.
        ENDIF
        EXIT
      ENDIF
    ENDDO

    CALL sdf_close(sdf_handle)

  END SUBROUTINE read_cpu_split



  FUNCTION it_part(array, npart_this_it, start, direction, param)

    REAL(num) :: it_part
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction
    INTEGER, INTENT(IN), OPTIONAL :: param

    INTEGER(i8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur => iterator_list

    DO ipart = 1, npart_this_it
      cur%part_pos = array(ipart)
      cur => cur%next
    ENDDO

    it_part = 0

  END FUNCTION it_part



  FUNCTION it_px(array, npart_this_it, start, param)

    REAL(num) :: it_px
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(1) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_px = 0

  END FUNCTION it_px



  FUNCTION it_py(array, npart_this_it, start, param)

    REAL(num) :: it_py
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(2) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_py = 0

  END FUNCTION it_py



  FUNCTION it_pz(array, npart_this_it, start, param)

    REAL(num) :: it_pz
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(3) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_pz = 0

  END FUNCTION it_pz



#ifdef PER_PARTICLE_WEIGHT
  FUNCTION it_weight(array, npart_this_it, start, param)

    REAL(num) :: it_weight
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%weight = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_weight = 0

  END FUNCTION it_weight
#endif



#if PARTICLE_ID || PARTICLE_ID4
  FUNCTION it_id(array, npart_this_it, start, param)

    REAL(num) :: it_id
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
#ifdef PARTICLE_ID4
      iterator_list%id = NINT(array(ipart))
#else
      iterator_list%id = NINT(array(ipart),i8)
#endif
      iterator_list => iterator_list%next
    ENDDO

    it_id = 0

  END FUNCTION it_id
#endif



#ifdef PHOTONS
  FUNCTION it_optical_depth(array, npart_this_it, start, param)

    REAL(num) :: it_optical_depth
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%optical_depth = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_optical_depth = 0

  END FUNCTION it_optical_depth



  FUNCTION it_qed_energy(array, npart_this_it, start, param)

    REAL(num) :: it_qed_energy
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%particle_energy = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_qed_energy = 0

  END FUNCTION it_qed_energy



#ifdef TRIDENT_PHOTONS
  FUNCTION it_optical_depth_trident(array, npart_this_it, start, param)

    REAL(num) :: it_optical_depth_trident
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN), OPTIONAL :: param
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%optical_depth_tri = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_optical_depth_trident = 0

  END FUNCTION it_optical_depth_trident
#endif
#endif

END MODULE setup
