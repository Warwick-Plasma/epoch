MODULE setup

  USE mpi
  USE sdf
  USE encoded_source
  USE fields
  USE mpi_subtype_control
  USE partlist
  USE shared_data
  USE strings
  USE version_data
  USE welcome
  USE random_generator
  USE split_particle

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control, minimal_init, restart_data
  PUBLIC :: open_files, close_files, flush_stat_file
  PUBLIC :: setup_species, after_deck_last

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

    dumpmask = 0
    comm = MPI_COMM_NULL

    dt_plasma_frequency = 0.0_num
    dt_multiplier = 0.95_num
    stdout_frequency = 0
    cpml_thickness = 6
    cpml_kappa_max = 20.0_num
    cpml_a_max = 0.15_num
    cpml_sigma_max = 0.7_num

    window_shift = 0.0_num
    npart_global = -1
    smooth_currents = .FALSE.
    dlb = .FALSE.
    use_random_seed = .FALSE.
    use_offset_grid = .FALSE.
    use_particle_lists = .FALSE.
    force_final_to_be_restartable = .FALSE.
    dump_source_code = .TRUE.
    dump_input_decks = .TRUE.
    full_dump_every = -1
    restart_dump_every = -1
    dt_snapshot = -1.0_num
    nstep_snapshot = -1
    nsteps = -1
    t_end = HUGE(1.0_num)
    particles_max_id = 0

    NULLIFY(laser_x_min)
    NULLIFY(laser_x_max)

    NULLIFY(dist_fns)

    run_date = get_unix_time()

    CALL set_field_order(2)

    CALL init_source_code()

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

  END SUBROUTINE minimal_init



  SUBROUTINE after_control

    INTEGER :: iproc, ix
    REAL(num) :: xb_min

    length_x = x_max - x_min
    dx = length_x / REAL(nx_global-2*cpml_thickness, num)
    x_min = x_min - dx * cpml_thickness
    x_max = x_max + dx * cpml_thickness
    length_x = x_max - x_min

    ! Shift grid to cell centres.
    ! At some point the grid may be redefined to be node centred.

    xb_min = x_min
    x_min = x_min + dx / 2.0_num
    x_max = x_max - dx / 2.0_num

    ! Setup global grid
    DO ix = -2, nx_global + 3
      x_global(ix) = x_min + (ix - 1) * dx
    ENDDO
    DO ix = 1, nx_global + 1
      xb_global(ix) = xb_min + (ix - 1) * dx
      xb_offset_global(ix) = xb_global(ix)
    ENDDO

    DO iproc = 0, nprocx-1
      x_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO

    x_min_local = x_mins(x_coords)
    x_max_local = x_maxs(x_coords)

    nx_global_min = cell_x_min(x_coords+1)
    nx_global_max = cell_x_max(x_coords+1)

    ! Setup local grid
    DO ix = -2, nx + 3
      x(ix) = x_global(nx_global_min+ix-1)
    ENDDO

    CALL set_initial_values

  END SUBROUTINE after_control



  SUBROUTINE after_deck_last

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
      dumpmask(c_dump_cpml_psi_eyx) = 0
      dumpmask(c_dump_cpml_psi_ezx) = 0
      dumpmask(c_dump_cpml_psi_byx) = 0
      dumpmask(c_dump_cpml_psi_bzx) = 0
    ENDIF

  END SUBROUTINE after_deck_last



  SUBROUTINE setup_data_averaging()

    INTEGER :: io, nspec_local

    dt_min_average = -1.0_num
    IF (.NOT. any_average) RETURN

    DO io = 1, num_vars_to_dump
      IF (IAND(dumpmask(io), c_io_averaged) .NE. 0) THEN
        nspec_local = 0
        IF (IAND(dumpmask(io), c_io_no_sum) .EQ. 0) &
            nspec_local = 1
        IF (IAND(dumpmask(io), c_io_species) .NE. 0) &
            nspec_local = nspec_local + n_species

        IF (nspec_local .LE. 0) CYCLE

        ALLOCATE(averaged_data(io)%array(-2:nx+3,nspec_local))
        averaged_data(io)%array = 0.0_num
        averaged_data(io)%real_time = c_non_zero
      ENDIF
    ENDDO

    IF (nstep_average .GT. 0 .AND. dt_average .GT. 0) THEN
      dt_min_average = dt_average / REAL(nstep_average, num)
    ENDIF

  END SUBROUTINE setup_data_averaging



  SUBROUTINE setup_species

    INTEGER :: ispecies

    ALLOCATE(species_list(n_species))
    ALLOCATE(io_list_data(n_species))
    ALLOCATE(particle_file_lengths(n_species))
    ALLOCATE(particle_file_offsets(n_species))
    io_list => species_list

    DO ispecies = 1, n_species
      species_list(ispecies)%name = blank
      species_list(ispecies)%mass = -1.0_num
      species_list(ispecies)%charge = 0.0_num
      species_list(ispecies)%dumpmask = c_io_always
      species_list(ispecies)%count = -1
      species_list(ispecies)%id = 0
      species_list(ispecies)%npart_per_cell = -1
      species_list(ispecies)%density = 0.0_num
      species_list(ispecies)%temperature = 0.0_num
      species_list(ispecies)%split = .FALSE.
      species_list(ispecies)%npart_max = 0
      species_list(ispecies)%global_count = 0
      NULLIFY(species_list(ispecies)%next)
      NULLIFY(species_list(ispecies)%prev)
      NULLIFY(species_list(ispecies)%ext_temp_x_min)
      NULLIFY(species_list(ispecies)%ext_temp_x_max)
      NULLIFY(species_list(ispecies)%secondary_list)
#ifdef PARTICLE_IONISE
      species_list(ispecies)%ionise = .FALSE.
      species_list(ispecies)%ionise_to_species = -1
      species_list(ispecies)%release_species = -1
      species_list(ispecies)%ionisation_energy = 0.0_num
#endif
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

    CHARACTER(LEN=16) :: string
    INTEGER :: errcode, ierr
    LOGICAL :: exists

    IF (rank .EQ. 0) THEN
#ifdef NO_IO
      stat_file = '/dev/null'
#else
      WRITE(stat_file, '(a, ''/epoch1d.dat'')') TRIM(data_dir)
#endif
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
        CALL MPI_ABORT(comm, errcode, ierr)
        STOP
      ENDIF
      IF (ic_from_restart) THEN
        CALL integer_as_string(restart_snapshot, string)
        WRITE(stat_unit,*)
        WRITE(stat_unit,*) 'Restarting from ', TRIM(string)
        WRITE(stat_unit,*) ascii_header
      ELSE
        WRITE(stat_unit,*) ascii_header
        WRITE(stat_unit,*)
      ENDIF
    ENDIF

  END SUBROUTINE open_files



  SUBROUTINE flush_stat_file

    INTEGER :: errcode

    IF (rank .EQ. 0) THEN
      CLOSE(unit=stat_unit)
      OPEN(unit=stat_unit, status='OLD', position='APPEND', &
          file=stat_file, iostat=errcode)
    ENDIF

  END SUBROUTINE flush_stat_file



  SUBROUTINE close_files

    IF (rank .EQ. 0) CLOSE(unit=stat_unit)

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



  SUBROUTINE find_species_by_name(specname, species_number)

    CHARACTER(LEN=*), INTENT(IN) :: specname
    INTEGER, INTENT(OUT) :: species_number
    INTEGER :: ispecies, ierr, i1, i2

    CALL strip_species_name(specname, i1, i2)

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(specname(i1:i2), species_list(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE find_species_by_name



  SUBROUTINE strip_species_name(name, i1, i2)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: i1, i2
    INTEGER :: ii

    i2 = LEN_TRIM(name)
    DO ii = i2, 1, -1
      IF (name(ii:ii) .EQ. '/') RETURN
      i1 = ii
    ENDDO

  END SUBROUTINE strip_species_name



  SUBROUTINE restart_data(step)

    INTEGER, INTENT(OUT) :: step
    CHARACTER(LEN=20+data_dir_max_length) :: filename
    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, str1
    CHARACTER(LEN=c_max_string_length) :: name, len_string
    INTEGER :: blocktype, datatype, code_io_version, string_len, ispecies
    INTEGER :: ierr, i1, i2, iblock, nblocks, ndims, found_species
    INTEGER(KIND=8) :: npart, npart_local
    INTEGER, DIMENSION(4) :: dims
    LOGICAL :: restart_flag
    TYPE(sdf_file_handle) :: sdf_handle
    TYPE(particle_species), POINTER :: species
    TYPE(particle_list), POINTER :: partlist
    INTEGER, POINTER :: species_subtypes(:)

    npart_global = 0
    step = -1

    ! Create the filename for the last snapshot
    WRITE(filename, '(a, ''/'', i4.4, ''.sdf'')') TRIM(data_dir), &
        restart_snapshot
    CALL sdf_open(sdf_handle, filename, rank, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, step, time, code_name, code_io_version, &
        string_len, restart_flag)

    IF (dt_snapshot .GT. 0.0_num) THEN
      time_next = time + dt_snapshot
    ELSE
      time_next = time
    ENDIF
    IF (nstep_snapshot .GT. 0) THEN
      nstep_next = step + nstep_snapshot
    ELSE
      nstep_next = step
    ENDIF

    IF (.NOT. restart_flag) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file is not a restart dump. Unable to continue.'
      ENDIF
      CALL MPI_ABORT(comm, errcode, ierr)
      STOP
    ENDIF

    IF (.NOT.str_cmp(code_name, 'Epoch1d')) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Epoch1d. Unable to ', &
            'continue.'
      ENDIF
      CALL MPI_ABORT(comm, errcode, ierr)
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
      CALL MPI_ABORT(comm, errcode, ierr)
      STOP
    ENDIF

    IF (rank .EQ. 0) PRINT*, 'Loading snapshot for time', time

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    CALL create_ascii_header
    WRITE(stat_unit,*) ascii_header
    WRITE(stat_unit,*)

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    IF (rank .EQ. 0) PRINT*, 'Input file contains', nblocks, 'blocks'

    CALL sdf_read_blocklist(sdf_handle)

    ! Scan file for particle species and allocate storage
    found_species = 0
    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      IF (blocktype .EQ. c_blocktype_point_mesh) THEN
        CALL sdf_read_point_mesh_info(sdf_handle, npart)

        CALL find_species_by_name(name, ispecies)
        IF (ispecies .EQ. 0) THEN
          IF (rank .EQ. 0) THEN
            CALL strip_species_name(name, i1, i2)
            PRINT*, '*** WARNING ***'
            PRINT*, 'Particle species "', name(i1:i2), '" from restart dump ', &
                'not found in input deck. Ignoring.'
          ENDIF
          CYCLE
        ENDIF
        species => species_list(ispecies)

        IF (ASSOCIATED(species%attached_list%head)) THEN
          IF (rank .EQ. 0) THEN
            CALL strip_species_name(name, i1, i2)
            PRINT*, '*** ERROR ***'
            PRINT*, 'Duplicate meshes for species "', name(i1:i2),'"'
          ENDIF
          CALL MPI_ABORT(comm, errcode, ierr)
          STOP
        ENDIF

        npart_local = npart / nproc
        IF (npart_local * nproc .NE. npart) THEN
          IF (rank .LT. npart - npart_local * nproc) &
              npart_local = npart_local + 1
        ENDIF

        CALL create_allocated_partlist(species%attached_list, npart_local)

        npart_global = npart_global + npart
        species%count = npart
        found_species = found_species + 1
      ENDIF
    ENDDO

    IF (found_species .NE. n_species) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Number of species in restart dump does not match input.deck'
      ENDIF
      CALL MPI_ABORT(comm, errcode, ierr)
      STOP
    ENDIF

    CALL create_subtypes_for_load(species_subtypes)

    CALL sdf_seek_start(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_constant)
        IF (str_cmp(block_id, 'dt_plasma_frequency')) THEN
          CALL sdf_read_srl(sdf_handle, dt_plasma_frequency)
        ELSE IF (block_id(1:7) .EQ. 'weight/') THEN
          CALL find_species_by_name(block_id, ispecies)
          IF (ispecies .EQ. 0) CYClE
          CALL sdf_read_srl(sdf_handle, species_list(ispecies)%weight)
        ENDIF
      !CASE(c_blocktype_plain_mesh)
        !CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims)

        !CALL sdf_read_srl_plain_mesh(sdf_handle, x, y)
      CASE(c_blocktype_point_mesh)
        CALL sdf_read_point_mesh_info(sdf_handle, npart)

        CALL find_species_by_name(name, ispecies)
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
          CALL MPI_ABORT(comm, errcode, ierr)
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
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'jy')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jy, &
              subtype_field, subarray_field)

        ELSE IF (str_cmp(block_id, 'jz')) THEN
          CALL sdf_read_plain_variable(sdf_handle, jz, &
              subtype_field, subarray_field)

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
        CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id)

        CALL find_species_by_name(mesh_id, ispecies)
        IF (ispecies .EQ. 0) CYClE
        species => species_list(ispecies)

        IF (npart .NE. species%count) THEN
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Malformed restart dump. Number of particle variables', &
                ' does not match grid.'
          ENDIF
          CALL MPI_ABORT(comm, errcode, ierr)
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
          CALL MPI_ABORT(comm, errcode, ierr)
          STOP
#endif
        ENDIF
      END SELECT
    ENDDO

    CALL sdf_close(sdf_handle)
    CALL free_subtypes_for_load(species_subtypes)

  END SUBROUTINE restart_data



  FUNCTION it_part(array, npart_this_it, start, direction)

    REAL(num) :: it_part
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction

    INTEGER(KIND=8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur => iterator_list

    DO ipart = 1, npart_this_it
      cur%part_pos = array(ipart)
      cur => cur%next
    ENDDO

    it_part = 0

  END FUNCTION it_part



  FUNCTION it_px(array, npart_this_it, start)

    REAL(num) :: it_px
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(1) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_px = 0

  END FUNCTION it_px



  FUNCTION it_py(array, npart_this_it, start)

    REAL(num) :: it_py
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(2) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_py = 0

  END FUNCTION it_py



  FUNCTION it_pz(array, npart_this_it, start)

    REAL(num) :: it_pz
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%part_p(3) = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_pz = 0

  END FUNCTION it_pz



#ifdef PER_PARTICLE_WEIGHT
  FUNCTION it_weight(array, npart_this_it, start)

    REAL(num) :: it_weight
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
      iterator_list%weight = array(ipart)
      iterator_list => iterator_list%next
    ENDDO

    it_weight = 0

  END FUNCTION it_weight
#endif



#if PARTICLE_ID || PARTICLE_ID4
  FUNCTION it_id(array, npart_this_it, start)

    REAL(num) :: it_id
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER :: ipart

    DO ipart = 1, npart_this_it
#ifdef PARTICLE_ID4
      iterator_list%id = NINT(array(ipart))
#else
      iterator_list%id = NINT(array(ipart),8)
#endif
      iterator_list => iterator_list%next
    ENDDO

    it_id = 0

  END FUNCTION it_id
#endif

END MODULE setup
