MODULE setup

  USE sdf
  USE encoded_source
  USE fields
  USE mpi_subtype_control
  USE partlist
  USE shared_data
  USE strings
  USE version_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control, minimal_init, restart_data
  PUBLIC :: open_files, close_files
  PUBLIC :: setup_species

  TYPE(particle), POINTER, SAVE :: iterator_list

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
    dt_multiplier = 1.0_num
    stdout_frequency = 0

    window_shift = 0.0_num
    npart_global = -1
    dlb = .FALSE.
    use_random_seed = .FALSE.
    use_offset_grid = .FALSE.
    force_final_to_be_restartable = .FALSE.
    full_dump_every = -1
    restart_dump_every = -1

    NULLIFY(laser_x_min)
    NULLIFY(laser_x_max)
    NULLIFY(laser_y_max)
    NULLIFY(laser_y_min)
    NULLIFY(laser_z_max)
    NULLIFY(laser_z_min)

    NULLIFY(dist_fns)

    run_date = get_unix_time()

    CALL set_field_order(2)

    CALL init_source_code()

    stagger = 0
    stagger(1,c_stagger_ex) = 1
    stagger(2,c_stagger_ey) = 1
    stagger(3,c_stagger_ez) = 1

    stagger(1:3,c_stagger_bx) = 1
    stagger(1:3,c_stagger_by) = 1
    stagger(1:3,c_stagger_bz) = 1
    stagger(1,c_stagger_bx) = 0
    stagger(2,c_stagger_by) = 0
    stagger(3,c_stagger_bz) = 0

  END SUBROUTINE minimal_init



  SUBROUTINE after_control

    INTEGER :: iproc, ix, iy, iz
    REAL(num) :: xb_min, yb_min, zb_min

    length_x = x_max - x_min
    length_y = y_max - y_min
    length_z = z_max - z_min
    dx = length_x / REAL(nx_global, num)
    dy = length_y / REAL(ny_global, num)
    dz = length_z / REAL(nz_global, num)

    ! Shift grid to cell centres.
    ! At some point the grid may be redefined to be node centred.

    xb_min = x_min
    yb_min = y_min
    zb_min = z_min
    x_min = x_min + dx / 2.0_num
    x_max = x_max - dx / 2.0_num
    y_min = y_min + dy / 2.0_num
    y_max = y_max - dy / 2.0_num
    z_min = z_min + dz / 2.0_num
    z_max = z_max - dz / 2.0_num
    length_x = x_max - x_min
    length_y = y_max - y_min
    length_z = z_max - z_min

    ! Setup global grid
    DO ix = -2, nx_global + 3
      x_global(ix) = x_min + (ix - 1) * dx
    ENDDO
    DO ix = 1, nx_global + 1
      xb_global(ix) = xb_min + (ix - 1) * dx
      xb_offset_global(ix) = xb_global(ix)
    ENDDO
    DO iy = -2, ny_global + 3
      y_global(iy) = y_min + (iy - 1) * dy
    ENDDO
    DO iy = 1, ny_global + 1
      yb_global(iy) = yb_min + (iy - 1) * dy
      yb_offset_global(iy) = yb_global(iy)
    ENDDO
    DO iz = -2, nz_global + 3
      z_global(iz) = z_min + (iz - 1) * dz
    ENDDO
    DO iz = 1, nz_global + 1
      zb_global(iz) = zb_min + (iz - 1) * dz
      zb_offset_global(iz) = zb_global(iz)
    ENDDO

    DO iproc = 0, nprocx-1
      x_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocy-1
      y_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocz-1
      z_mins(iproc) = z_global(cell_z_min(iproc+1))
      z_maxs(iproc) = z_global(cell_z_max(iproc+1))
    ENDDO

    x_min_local = x_mins(coordinates(c_ndims))
    x_max_local = x_maxs(coordinates(c_ndims))
    y_min_local = y_mins(coordinates(c_ndims-1))
    y_max_local = y_maxs(coordinates(c_ndims-1))
    z_min_local = z_mins(coordinates(c_ndims-2))
    z_max_local = z_maxs(coordinates(c_ndims-2))

    nx_global_min = cell_x_min(coordinates(c_ndims)+1)
    nx_global_max = cell_x_max(coordinates(c_ndims)+1)
    ny_global_min = cell_y_min(coordinates(c_ndims-1)+1)
    ny_global_max = cell_y_max(coordinates(c_ndims-1)+1)
    nz_global_min = cell_z_min(coordinates(c_ndims-2)+1)
    nz_global_max = cell_z_max(coordinates(c_ndims-2)+1)

    ! Setup local grid
    DO ix = -2, nx + 3
      x(ix) = x_global(nx_global_min+ix-1)
    ENDDO
    DO iy = -2, ny + 3
      y(iy) = y_global(ny_global_min+iy-1)
    ENDDO
    DO iz = -2, nz + 3
      z(iz) = z_global(nz_global_min+iz-1)
    ENDDO

    CALL set_initial_values

    CALL setup_data_averaging

  END SUBROUTINE after_control



  SUBROUTINE setup_data_averaging()

    INTEGER :: ioutput, n_species_local
    REAL(num) :: min_av_time

    dt_min_average = -1.0_num
    IF (.NOT. any_average) RETURN

    min_av_time = t_end
    DO ioutput = 1, num_vars_to_dump
      IF (IAND(dumpmask(ioutput), c_io_averaged) .NE. 0) THEN
        averaged_data(ioutput)%average_over_real_time = average_time
        min_av_time = &
            MIN(min_av_time, averaged_data(ioutput)%average_over_real_time)
        n_species_local = 1
        IF (IAND(dumpmask(ioutput), c_io_species) .NE. 0) &
            n_species_local = n_species + 1
        ALLOCATE(averaged_data(ioutput)%array(-2:nx+3,-2:ny+3,-2:nz+3, &
            n_species_local))
        averaged_data(ioutput)%array = 0.0_num
        averaged_data(ioutput)%real_time_after_average = 0.0_num
      ELSE
        dumpmask(ioutput) = IOR(dumpmask(ioutput), c_io_snapshot)
      ENDIF
    ENDDO

    IF (min_cycles_per_average .GT. 0) &
        dt_min_average = min_av_time / REAL(min_cycles_per_average, num)

  END SUBROUTINE setup_data_averaging



  SUBROUTINE setup_species

    INTEGER :: ispecies

    ALLOCATE(particle_species(n_species))
    ALLOCATE(particle_file_lengths(n_species))
    ALLOCATE(particle_file_offsets(n_species))

    DO ispecies = 1, n_species
      particle_species(ispecies)%name = blank
      particle_species(ispecies)%mass = -1.0_num
      particle_species(ispecies)%charge = 0.0_num
      particle_species(ispecies)%dump = .TRUE.
      particle_species(ispecies)%count = -1
      particle_species(ispecies)%id = 0
      particle_species(ispecies)%npart_per_cell = 0
      NULLIFY(particle_species(ispecies)%density)
      NULLIFY(particle_species(ispecies)%temperature)
      NULLIFY(particle_species(ispecies)%next)
      NULLIFY(particle_species(ispecies)%prev)
#ifdef SPLIT_PARTICLES_AFTER_PUSH
      particle_species(ispecies)%split = .FALSE.
      particle_species(ispecies)%npart_max = 0
      particle_species(ispecies)%global_count = 0
      NULLIFY(particle_species(ispecies)%secondary_list)
#endif
#ifdef PARTICLE_IONISE
      particle_species(ispecies)%ionise = .FALSE.
      particle_species(ispecies)%ionise_to_species = -1
      particle_species(ispecies)%release_species = -1
      particle_species(ispecies)%critical_field = 0.0_num
      particle_species(ispecies)%ionisation_energy = 0.0_num
#endif
#ifdef TRACER_PARTICLES
      particle_species(ispecies)%tracer = .FALSE.
#endif
#ifdef PARTICLE_PROBES
      NULLIFY(particle_species(ispecies)%attached_probes)
#endif
    ENDDO

  END SUBROUTINE setup_species



  SUBROUTINE open_files

    CHARACTER(LEN=11+data_dir_max_length) :: file2
    CHARACTER(LEN=16) :: string
    INTEGER :: errcode, ierr
    LOGICAL :: exists

    IF (rank .EQ. 0) THEN
      WRITE(file2, '(a, "/epoch3d.dat")') TRIM(data_dir)
      IF (ic_from_restart) THEN
        INQUIRE(file=file2, exist=exists)
        IF (exists) THEN
          OPEN(unit=20, status='OLD', access='APPEND', file=file2, &
              iostat=errcode)
        ELSE
          OPEN(unit=20, status='NEW', file=file2, iostat=errcode)
        ENDIF
      ELSE
        OPEN(unit=20, status='REPLACE', file=file2, iostat=errcode)
      ENDIF
      IF (errcode .NE. 0) THEN
        PRINT*, '***ERROR***'
        PRINT*, 'Cannot create "epoch3d.dat" output file. The most common ' &
            // 'cause of this problem '
        PRINT*, 'is that the ouput directory does not exist'
        CALL MPI_ABORT(comm, errcode, ierr)
        STOP
      ENDIF
      IF (ic_from_restart) THEN
        CALL integer_as_string(restart_snapshot, string)
        WRITE(20,*)
        WRITE(20,*) 'Restarting from ', TRIM(string)
      ENDIF
      WRITE(20,*) ascii_header
      WRITE(20,*)
    ENDIF

  END SUBROUTINE open_files



  SUBROUTINE close_files

    CLOSE(unit=20)

  END SUBROUTINE close_files



  SUBROUTINE set_initial_values

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

  END SUBROUTINE set_initial_values



  SUBROUTINE find_species_by_name(specname, species)

    CHARACTER(LEN=*), INTENT(IN) :: specname
    TYPE(particle_family), POINTER :: species
    INTEGER :: ispecies, species_number, errcode, ierr

    species_number = 0
    DO ispecies = 1,n_species
      IF (str_cmp(specname,particle_species(ispecies)%name)) THEN
        species_number = ispecies
        EXIT
      ENDIF
    ENDDO

    IF (species_number .EQ. 0) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Particle species "', TRIM(specname),'" from restart dump ', &
            'not found in input deck.'
      ENDIF
      CALL MPI_ABORT(comm, errcode, ierr)
      STOP
    ENDIF

    species => particle_species(species_number)

  END SUBROUTINE find_species_by_name



  SUBROUTINE restart_data(snap)

    INTEGER, INTENT(OUT) :: snap
    CHARACTER(LEN=20+data_dir_max_length) :: filename
    CHARACTER(LEN=c_id_length) :: code_name, block_id, mesh_id, str1, str2, str3
    CHARACTER(LEN=c_max_string_length) :: name
    INTEGER :: geometry, blocktype, datatype, code_io_version
    INTEGER :: ierr, ii, i1, i2, iblock, nblocks, ndims
    INTEGER(KIND=8) :: npart, npart_local
    INTEGER, DIMENSION(4) :: dims
    LOGICAL :: constant_weight, restart_flag
    TYPE(sdf_file_handle) :: sdf_handle
    TYPE(particle), POINTER :: current
    TYPE(particle_family), POINTER :: species

    npart_global = 0
    constant_weight = .FALSE.
    snap = -1

    ! Create the filename for the last snapshot
    WRITE(filename, '(a, "/", i4.4, ".sdf")') TRIM(data_dir), restart_snapshot
    CALL sdf_open(sdf_handle, filename, rank, comm, c_sdf_read)

    CALL sdf_read_header(sdf_handle, snap, time, code_name, code_io_version, &
        restart_flag)

    IF (.NOT. restart_flag) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF file is not a restart dump. Unable to continue.'
      ENDIF
      CALL MPI_ABORT(comm, errcode, ierr)
      STOP
    ENDIF

    IF (str_cmp(code_name, 'Eden3d')) THEN
      IF (rank .EQ. 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'SDF restart file was not generated by Eden3d. Unable to ', &
            'continue.'
      ENDIF
      CALL MPI_ABORT(comm, errcode, ierr)
      STOP
    ENDIF

    IF (rank .EQ. 0) PRINT*, 'Loading snapshot for time', time

    nblocks = sdf_read_nblocks(sdf_handle)
    jobid = sdf_read_jobid(sdf_handle)

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    IF (rank .EQ. 0) PRINT*, 'Input file contains', nblocks, 'blocks'

    CALL sdf_read_blocklist(sdf_handle)

    DO iblock = 1, nblocks
      CALL sdf_read_next_block_header(sdf_handle, block_id, name, blocktype, &
          ndims, datatype)

      SELECT CASE(blocktype)
      CASE(c_blocktype_constant)
        IF (str_cmp(block_id, 'weight')) THEN
          CALL sdf_read_srl(sdf_handle, weight)
          constant_weight = .TRUE.
        ENDIF
      !CASE(c_blocktype_plain_mesh)
        !CALL sdf_read_plain_mesh_info(sdf_handle, geometry, dims)

        !CALL sdf_read_srl_plain_mesh(sdf_handle, x, y, z)
      CASE(c_blocktype_point_mesh)
        CALL sdf_read_point_mesh_info(sdf_handle, npart)

        i2 = LEN_TRIM(name)
        DO ii = 1,i2
          i1 = ii
          IF (name(ii:ii) .EQ. '/') EXIT
        ENDDO
        i1 = i1 + 1

        CALL find_species_by_name(name(i1:i2), species)

        IF (ASSOCIATED(species%attached_list%head)) THEN
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Duplicate meshes for species "', TRIM(name(i1:i2)),'"'
          ENDIF
          CALL MPI_ABORT(comm, errcode, ierr)
          STOP
        ENDIF

        npart_local = npart / nproc
        IF (npart_local * nproc .NE. npart) THEN
          IF (rank .EQ. 0) THEN
            PRINT*, 'Cannot evenly subdivide particles over', nproc, &
                'processors. Trying to fix'
          ENDIF
          IF (rank .LT. npart - npart_local * nproc) &
              npart_local = npart_local + 1
        ENDIF

        CALL create_subtypes_for_load(npart_local)
        CALL create_allocated_partlist(species%attached_list, npart_local)

        iterator_list => species%attached_list%head
        npart_global = npart_global + npart
        species%count = npart

        CALL sdf_read_point_mesh(sdf_handle, npart_local, &
            subtype_particle_var, it_part)

      CASE(c_blocktype_plain_variable)
        CALL sdf_read_plain_variable_info(sdf_handle, dims)

        IF (dims(1) .NE. nx_global .OR. dims(2) .NE. ny_global &
            .OR. dims(3) .NE. nz_global) THEN
          IF (rank .EQ. 0) THEN
            PRINT*, '*** ERROR ***'
            PRINT*, 'Number of gridpoints in restart dump does not match', &
                ' the input deck.'
            CALL integer_as_string(nx_global, str1)
            CALL integer_as_string(ny_global, str2)
            CALL integer_as_string(nz_global, str3)
            PRINT*, 'Input deck grid: ', TRIM(str1), ',', TRIM(str2), &
                ',', TRIM(str3)
            CALL integer_as_string(dims(1), str1)
            CALL integer_as_string(dims(2), str2)
            CALL integer_as_string(dims(3), str3)
            PRINT*, 'Restart dump grid: ',TRIM(str1), ',', TRIM(str2), &
                ',', TRIM(str3)
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

        ENDIF

      CASE(c_blocktype_point_variable)
        CALL sdf_read_point_variable_info(sdf_handle, npart, mesh_id)

        i2 = LEN_TRIM(mesh_id)
        DO ii = 1,i2
          i1 = ii
          IF (mesh_id(ii:ii) .EQ. '/') EXIT
        ENDDO
        i1 = i1 + 1

        CALL find_species_by_name(mesh_id(i1:i2), species)

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
        npart = species%attached_list%count

        IF (str_cmp(block_id, 'px')) THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              subtype_particle_var, it_px)

        ELSE IF (str_cmp(block_id, 'py')) THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              subtype_particle_var, it_py)

        ELSE IF (str_cmp(block_id, 'pz')) THEN
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              subtype_particle_var, it_pz)

        ELSE IF (str_cmp(block_id, 'weight')) THEN
#ifdef PER_PARTICLE_WEIGHT
          CALL sdf_read_point_variable(sdf_handle, npart_local, &
              subtype_particle_var, it_weight)
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
    CALL free_subtypes_for_load()

#ifdef PER_PARTICLE_WEIGHT
    IF (constant_weight) THEN
      DO ii = 1,n_species
        current => particle_species(ii)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          current%weight = weight
          current => current%next
        ENDDO
      ENDDO
    ENDIF
#endif

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
      cur%part_pos(direction) = array(ipart)
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

END MODULE setup
