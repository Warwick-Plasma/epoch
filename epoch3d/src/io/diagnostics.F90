MODULE diagnostics

  USE calc_df
  USE sdf
  USE deck
  USE dist_fn
  USE encoded_source
  USE iterators
  USE probes
  USE version_data
  USE setup
  USE random_generator
  USE strings

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_routines

  TYPE(sdf_file_handle) :: sdf_handle
  INTEGER(i8), ALLOCATABLE :: species_offset(:)
  INTEGER(i8), ALLOCATABLE :: ejected_offset(:)
  LOGICAL :: reset_ejected, done_species_offset_init, done_subset_init
  INTEGER :: isubset
  INTEGER, DIMENSION(num_vars_to_dump) :: iomask
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: iodumpmask

CONTAINS

  SUBROUTINE output_routines(step)   ! step = step index

    INTEGER, INTENT(INOUT) :: step
    LOGICAL :: print_arrays, first_call, last_call
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename, filename_desc
    CHARACTER(LEN=c_max_string_length) :: dump_type
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array
    INTEGER :: code, i, io, random_state(4)
    INTEGER, ALLOCATABLE :: random_states_per_proc(:)
    INTEGER, DIMENSION(c_ndims) :: dims
    LOGICAL :: restart_flag, convert
    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species

    CHARACTER(LEN=1), DIMENSION(3) :: dim_tags = (/'x', 'y', 'z'/)
    CHARACTER(LEN=5), DIMENSION(6) :: dir_tags = &
        (/'x_max', 'y_max', 'z_max', 'x_min', 'y_min', 'z_min'/)
    INTEGER, DIMENSION(6) :: fluxdir = &
        (/c_dir_x, c_dir_y, c_dir_z, -c_dir_x, -c_dir_y, -c_dir_z/)

#ifdef NO_IO
    RETURN
#endif

    IF (rank .EQ. 0 .AND. stdout_frequency .GT. 0 &
        .AND. MOD(step, stdout_frequency) .EQ. 0) THEN
      WRITE(*, '(''Time'', g20.12, '' and iteration'', i7, '' after '', &
          & f8.1, '' seconds'')') time, step, MPI_WTIME() - walltime_start
    ENDIF

    CALL io_test(step, print_arrays, first_call, last_call)

    IF (.NOT.print_arrays) RETURN

    dims = (/nx_global, ny_global, nz_global/)

    ! Allows a maximum of 10^999 output dumps, should be enough for anyone
    ! (feel free to laugh when this isn't the case)
    WRITE(filename_desc, '(''(a, "/", i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
        n_zeros, n_zeros
    WRITE(filename, filename_desc) TRIM(data_dir), output_file

    ! Always dump the variables with the 'Every' attribute
    code = c_io_always

    ! Only dump variables with the 'FULL' attributre on full dump intervals
    IF (MOD(output_file, full_dump_every) .EQ. 0 &
        .AND. full_dump_every .GT. -1) code = IOR(code, c_io_full)
    IF (MOD(output_file, restart_dump_every) .EQ. 0 &
        .AND. restart_dump_every .GT. -1) code = IOR(code, c_io_restartable)
    IF (first_call .AND. force_first_to_be_restartable) &
        code = IOR(code, c_io_restartable)
    IF (last_call .AND. force_final_to_be_restartable) &
        code = IOR(code, c_io_restartable)

    CALL create_subtypes(code)

    ! Set a restart_flag to pass to the file header
    IF (IAND(code, c_io_restartable) .NE. 0) THEN
      restart_flag = .TRUE.
    ELSE
      restart_flag = .FALSE.
    ENDIF
    reset_ejected = .FALSE.

    ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))

    ! open the file
    CALL sdf_open(sdf_handle, filename, comm, c_sdf_write)
    CALL sdf_write_header(sdf_handle, 'Epoch3d', 1, step, time, restart_flag, &
        jobid)
    CALL sdf_write_run_info(sdf_handle, c_version, c_revision, c_minor_rev, &
        c_commit_id, sha1sum, c_compile_machine, c_compile_flags, defines, &
        c_compile_date, run_date)
    CALL sdf_write_cpu_split(sdf_handle, 'cpu_rank', 'CPUs/Original rank', &
        cell_x_max, cell_y_max, cell_z_max)

    IF (restart_flag) THEN
      CALL sdf_write_srl(sdf_handle, 'dt_plasma_frequency', 'Time increment', &
          dt_plasma_frequency)

      DO io = 1, n_io_blocks
        CALL sdf_write_srl(sdf_handle, &
            'time_prev/'//TRIM(io_block_list(io)%name), &
            'time_prev/'//TRIM(io_block_list(io)%name), &
            io_block_list(io)%time_prev)
        CALL sdf_write_srl(sdf_handle, &
            'nstep_prev/'//TRIM(io_block_list(io)%name), &
            'nstep_prev/'//TRIM(io_block_list(io)%name), &
            io_block_list(io)%nstep_prev)
      ENDDO

      DO ispecies = 1, n_species
        species => io_list(ispecies)
        CALL sdf_write_srl(sdf_handle, 'nppc/' // TRIM(species%name), &
            'Particles/Particles Per Cell/' // TRIM(species%name), &
            species%npart_per_cell)
      ENDDO

      IF (need_random_state) THEN
        CALL get_random_state(random_state)
        ALLOCATE(random_states_per_proc(4*nproc))
        CALL MPI_GATHER(random_state, 4, MPI_INTEGER, &
            random_states_per_proc, 4, MPI_INTEGER, 0, comm, errcode)
        CALL sdf_write_srl(sdf_handle, 'random_states', 'Random States', &
            random_states_per_proc)
        DEALLOCATE(random_states_per_proc)
      ENDIF
    ENDIF

    iomask = iodumpmask(1,:)

    ! Write the cartesian mesh
    IF (IAND(iomask(c_dump_grid), code) .NE. 0) THEN
      convert = (IAND(iomask(c_dump_grid), c_io_dump_single) .NE. 0 &
          .AND. (IAND(code,c_io_restartable) .EQ. 0 &
          .OR. IAND(iomask(c_dump_grid), c_io_restartable) .EQ. 0))
      IF (.NOT. use_offset_grid) THEN
        CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
            xb_global, yb_global, zb_global, convert)
      ELSE
        CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid', 'Grid/Grid', &
            xb_offset_global, yb_offset_global, zb_offset_global, convert)
        CALL sdf_write_srl_plain_mesh(sdf_handle, 'grid_full', &
            'Grid/Grid_Full', xb_global, yb_global, zb_global, convert)
      ENDIF
    ENDIF

    CALL write_field(c_dump_ex, code, 'ex', 'Electric Field/Ex', 'V/m', &
        c_stagger_ex, ex)
    CALL write_field(c_dump_ey, code, 'ey', 'Electric Field/Ey', 'V/m', &
        c_stagger_ey, ey)
    CALL write_field(c_dump_ez, code, 'ez', 'Electric Field/Ez', 'V/m', &
        c_stagger_ez, ez)

    CALL write_field(c_dump_bx, code, 'bx', 'Magnetic Field/Bx', 'T', &
        c_stagger_bx, bx)
    CALL write_field(c_dump_by, code, 'by', 'Magnetic Field/By', 'T', &
        c_stagger_by, by)
    CALL write_field(c_dump_bz, code, 'bz', 'Magnetic Field/Bz', 'T', &
        c_stagger_bz, bz)

    CALL write_field(c_dump_jx, code, 'jx', 'Current/Jx', 'A/m^2', &
        c_stagger_jx, jx)
    CALL write_field(c_dump_jy, code, 'jy', 'Current/Jy', 'A/m^2', &
        c_stagger_jy, jy)
    CALL write_field(c_dump_jz, code, 'jz', 'Current/Jz', 'A/m^2', &
        c_stagger_jz, jz)

    IF (cpml_boundaries) THEN
      CALL write_field(c_dump_cpml_psi_eyx, code, 'cpml_psi_eyx', 'CPML/Ey_x', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_eyx)
      CALL write_field(c_dump_cpml_psi_ezx, code, 'cpml_psi_ezx', 'CPML/Ez_x', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_ezx)
      CALL write_field(c_dump_cpml_psi_byx, code, 'cpml_psi_byx', 'CPML/By_x', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_byx)
      CALL write_field(c_dump_cpml_psi_bzx, code, 'cpml_psi_bzx', 'CPML/Bz_x', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_bzx)

      CALL write_field(c_dump_cpml_psi_exy, code, 'cpml_psi_exy', 'CPML/Ex_y', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_exy)
      CALL write_field(c_dump_cpml_psi_ezy, code, 'cpml_psi_ezy', 'CPML/Ez_y', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_ezy)
      CALL write_field(c_dump_cpml_psi_bxy, code, 'cpml_psi_bxy', 'CPML/Bx_y', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_bxy)
      CALL write_field(c_dump_cpml_psi_bzy, code, 'cpml_psi_bzy', 'CPML/Bz_y', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_bzy)

      CALL write_field(c_dump_cpml_psi_exz, code, 'cpml_psi_exz', 'CPML/Ex_z', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_exz)
      CALL write_field(c_dump_cpml_psi_eyz, code, 'cpml_psi_eyz', 'CPML/Ey_z', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_eyz)
      CALL write_field(c_dump_cpml_psi_bxz, code, 'cpml_psi_bxz', 'CPML/Bx_z', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_bxz)
      CALL write_field(c_dump_cpml_psi_byz, code, 'cpml_psi_byz', 'CPML/By_z', &
          'A/m^2', c_stagger_cell_centre, cpml_psi_byz)
    ENDIF

    IF (n_subsets .GT. 0) THEN
      DO i = 1, n_species
        CALL create_empty_partlist(io_list_data(i)%attached_list)
      ENDDO
    ENDIF

    DO isubset = 1, n_subsets + 1
      done_species_offset_init = .FALSE.
      done_subset_init = .FALSE.
      IF (isubset .GT. 1) io_list => io_list_data
      iomask = iodumpmask(isubset,:)

      CALL write_particle_grid(code)

#ifdef PER_PARTICLE_WEIGHT
      CALL write_particle_variable(c_dump_part_weight, code, 'Weight', '', &
          iterate_weight)
#else
      IF (IAND(iomask(c_dump_part_weight), code) .NE. 0) THEN
        CALL build_species_subset

        DO ispecies = 1, n_species
          species => io_list(ispecies)
          IF (IAND(species%dumpmask, code) .NE. 0 &
              .OR. IAND(code, c_io_restartable) .NE. 0) THEN
            CALL sdf_write_srl(sdf_handle, 'weight/' // TRIM(species%name), &
                'Particles/Weight/' // TRIM(species%name), species%weight)
          ENDIF
        ENDDO
      ENDIF
#endif
      CALL write_particle_variable(c_dump_part_px, code, 'Px', 'kg.m/s', &
          iterate_px)
      CALL write_particle_variable(c_dump_part_py, code, 'Py', 'kg.m/s', &
          iterate_py)
      CALL write_particle_variable(c_dump_part_pz, code, 'Pz', 'kg.m/s', &
          iterate_pz)

      CALL write_particle_variable(c_dump_part_vx, code, 'Vx', 'm/s', &
          iterate_vx)
      CALL write_particle_variable(c_dump_part_vy, code, 'Vy', 'm/s', &
          iterate_vy)
      CALL write_particle_variable(c_dump_part_vz, code, 'Vz', 'm/s', &
          iterate_vz)

      CALL write_particle_variable(c_dump_part_charge, code, 'Q', 'C', &
          iterate_charge)
      CALL write_particle_variable(c_dump_part_mass, code, 'Mass', 'kg', &
          iterate_mass)
      CALL write_particle_variable(c_dump_part_ek, code, 'Ek', 'J', &
          iterate_ek)
#ifdef PARTICLE_DEBUG
      CALL write_particle_variable(c_dump_part_grid, code, 'Processor', &
          '', iterate_processor)
      CALL write_particle_variable(c_dump_part_grid, code, 'Processor_at_t0', &
          '', iterate_processor0)
#endif
#if PARTICLE_ID || PARTICLE_ID4
      CALL write_particle_variable(c_dump_part_id, code, 'ID', '#', &
          iterate_id)
#endif

      ! These are derived variables from the particles
      CALL write_nspecies_field(c_dump_ekbar, code, 'ekbar', &
          'EkBar', 'J', c_stagger_cell_centre, &
          calc_ekbar, array)

      CALL write_nspecies_field(c_dump_mass_density, code, 'mass_density', &
          'Mass_Density', 'kg/m^3', c_stagger_cell_centre, &
          calc_mass_density, array)

      CALL write_nspecies_field(c_dump_charge_density, code, 'charge_density', &
          'Charge_Density', 'C/m^3', c_stagger_cell_centre, &
          calc_charge_density, array)

      CALL write_nspecies_field(c_dump_number_density, code, 'number_density', &
          'Number_Density', '1/m^3', c_stagger_cell_centre, &
          calc_number_density, array)

      CALL write_nspecies_field(c_dump_temperature, code, 'temperature', &
          'Temperature', 'K', c_stagger_cell_centre, &
          calc_temperature, array)

      CALL write_nspecies_field(c_dump_jx, code, 'jx', &
          'Jx', 'A/m^2', c_stagger_cell_centre, &
          calc_per_species_jx, array)

      CALL write_nspecies_field(c_dump_jy, code, 'jy', &
          'Jy', 'A/m^2', c_stagger_cell_centre, &
          calc_per_species_jy, array)

      CALL write_nspecies_field(c_dump_jz, code, 'jz', &
          'Jz', 'A/m^2', c_stagger_cell_centre, &
          calc_per_species_jz, array)

      CALL write_nspecies_flux(c_dump_ekflux, code, 'ekflux', &
          'EkFlux', 'W/m^2', c_stagger_cell_centre, &
          calc_ekflux, array, fluxdir, dir_tags)

      CALL write_field_flux(c_dump_poynt_flux, code, 'poynt_flux', &
          'Poynting Flux', 'W/m^2', c_stagger_cell_centre, &
          calc_poynt_flux, array, fluxdir(1:3), dim_tags)

      IF (isubset .NE. 1) THEN
        DO i = 1, n_species
          CALL append_partlist(species_list(i)%attached_list, &
              io_list(i)%attached_list)
        ENDDO
        DO i = 1, n_species
          CALL create_empty_partlist(io_list(i)%attached_list)
        ENDDO
      ENDIF
    ENDDO

    io_list => species_list
    iomask = iodumpmask(1,:)

    IF (IAND(iomask(c_dump_dist_fns), code) .NE. 0) THEN
      CALL write_dist_fns(sdf_handle, code)
    ENDIF

#ifdef PARTICLE_PROBES
    IF (IAND(iomask(c_dump_probes), code) .NE. 0) THEN
      CALL write_probes(sdf_handle, code)
    ENDIF
#endif

    IF (restart_flag) THEN
      IF (dump_input_decks) CALL write_input_decks(sdf_handle)
      IF (dump_source_code .AND. SIZE(source_code) .GT. 0) &
          CALL sdf_write_source_code(sdf_handle, 'code', &
              'base64_packed_source_code', source_code, last_line, 0)
    ENDIF

    IF (IAND(iomask(c_dump_absorption), code) .NE. 0) THEN
      CALL MPI_ALLREDUCE(laser_absorb_local, laser_absorbed, 1, mpireal, &
          MPI_SUM, comm, errcode)
      CALL MPI_ALLREDUCE(laser_inject_local, laser_injected, 1, mpireal, &
          MPI_SUM, comm, errcode)
      IF (laser_injected .GT. 0.0_num) THEN
        laser_absorbed = laser_absorbed / laser_injected
      ELSE
        laser_absorbed = 0.0_num
      ENDIF
      CALL sdf_write_srl(sdf_handle, 'abs_frac', 'Absorption/Abs_frac', &
          laser_absorbed)
      CALL sdf_write_srl(sdf_handle, 'laser_enTotal', &
          'Absorption/Laser_enTotal', laser_injected)
    ENDIF

    ! close the file
    CALL sdf_close(sdf_handle)

    IF (rank .EQ. 0) THEN
      DO io = 1, n_io_blocks
        IF (io_block_list(io)%dump) THEN
          dump_type = TRIM(io_block_list(io)%name)
          CALL append_filename(dump_type, output_file)
        ENDIF
      ENDDO
      IF (IAND(code, c_io_restartable) .NE. 0) THEN
        dump_type = 'restart'
        CALL append_filename(dump_type, output_file)
      ENDIF
      IF (IAND(code, c_io_full) .NE. 0) THEN
        dump_type = 'full'
        CALL append_filename(dump_type, output_file)
      ENDIF
      WRITE(stat_unit, '(''Wrote '', a7, '' dump number'', i5, '' at time'', &
          & g20.12, '' and iteration'', i7)') dump_type, output_file, time, step
      CALL flush_stat_file()
    ENDIF

    output_file = output_file + 1

    DEALLOCATE(array)
    IF (ALLOCATED(species_offset)) DEALLOCATE(species_offset)
    IF (ALLOCATED(ejected_offset)) DEALLOCATE(ejected_offset)
    CALL free_subtypes()

    IF (reset_ejected) THEN
      DO i = 1, n_species
        CALL destroy_partlist(ejected_list(i)%attached_list)
      ENDDO
    ENDIF

  END SUBROUTINE output_routines



  SUBROUTINE append_filename(listname, output_file)
    ! This routine updates a list of each output type (eg. full, dump, normal)
    ! that can be passed to VisIt to filter the output.

    CHARACTER(LEN=*), INTENT(IN) :: listname
    INTEGER, INTENT(IN) :: output_file
    CHARACTER(LEN=data_dir_max_length+64) :: listfile, filename_desc
    INTEGER :: ierr
    LOGICAL :: exists

    WRITE(filename_desc, '(''(i'', i3.3, ''.'', i3.3, '', ".sdf")'')') &
        n_zeros, n_zeros
    listfile = TRIM(data_dir) // '/' // TRIM(listname) // '.visit'

    INQUIRE(file=listfile, exist=exists)
    IF (exists) THEN
      OPEN(unit=lu, status='OLD', position='APPEND', file=listfile, iostat=ierr)
    ELSE
      OPEN(unit=lu, status='NEW', file=listfile, iostat=errcode)
    ENDIF

    WRITE(lu,filename_desc) output_file
    CLOSE(lu)

  END SUBROUTINE append_filename



  SUBROUTINE io_test(step, print_arrays, first_call, last_call)

    INTEGER, INTENT(IN) :: step
    LOGICAL, INTENT(OUT) :: print_arrays, first_call, last_call
    INTEGER :: id, io, is, nstep_next
    REAL(num) :: t0, t1, time_first, av_time_first
    LOGICAL, SAVE :: first = .TRUE.

    IF (.NOT.ALLOCATED(iodumpmask)) &
        ALLOCATE(iodumpmask(n_subsets+1,num_vars_to_dump))

    print_arrays = .FALSE.
    first_call = first
    iomask = c_io_never
    iodumpmask = c_io_never

    IF ((time .GE. t_end .OR. step .EQ. nsteps) .AND. dump_last) THEN
      last_call = .TRUE.
      print_arrays = .TRUE.
    ELSE
      last_call = .FALSE.
    ENDIF

    DO io = 1, n_io_blocks
      IF (last_call) THEN
        io_block_list(io)%dump = .TRUE.
      ELSE
        io_block_list(io)%dump = .FALSE.
      ENDIF
      IF (first .AND. dump_first) io_block_list(io)%dump = .TRUE.

      ! Work out the time that the next dump will occur based on the
      ! current timestep
      t0 = HUGE(1.0_num)
      t1 = HUGE(1.0_num)
      IF (io_block_list(io)%dt_snapshot .GE. 0.0_num) &
          t0 = io_block_list(io)%time_prev + io_block_list(io)%dt_snapshot
      IF (io_block_list(io)%nstep_snapshot .GE. 0) THEN
        nstep_next = io_block_list(io)%nstep_prev &
            + io_block_list(io)%nstep_snapshot
        t1 = time + dt * (nstep_next - step)
      ENDIF

      IF (t0 .LT. t1) THEN
        ! Next I/O dump based on dt_snapshot
        time_first = t0
        IF (io_block_list(io)%dt_snapshot .GT. 0 .AND. time .GE. t0) THEN
          io_block_list(io)%time_prev = time
          io_block_list(io)%dump = .TRUE.
        ENDIF
      ELSE
        ! Next I/O dump based on nstep_snapshot
        time_first = t1
        IF (io_block_list(io)%nstep_snapshot .GT. 0 &
            .AND. step .GE. nstep_next) THEN
          io_block_list(io)%nstep_prev = step
          io_block_list(io)%dump = .TRUE.
        ENDIF
      ENDIF

      IF (io_block_list(io)%dump) THEN
        print_arrays = .TRUE.
        iomask = IOR(iomask, io_block_list(io)%dumpmask)
        IF (n_subsets .NE. 0) THEN
          DO is = 1, n_subsets
            iodumpmask(1+is,:) = &
                IOR(iodumpmask(1+is,:), subset_list(is)%dumpmask(io,:))
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    IF (dt * nsteps .LT. time_first) THEN
      av_time_first = dt * nsteps
    ELSE
      av_time_first = time_first
    ENDIF

    DO io = 1, n_io_blocks
      IF (.NOT. io_block_list(io)%any_average) CYCLE

      IF (time .GE. av_time_first - io_block_list(io)%average_time) THEN
        DO id = 1, num_vars_to_dump
          IF (IAND(io_block_list(io)%dumpmask(id), c_io_averaged) .NE. 0) THEN
            CALL average_field(id, io_block_list(io)%averaged_data(id))
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    IF (first) THEN
      first = .FALSE.
      IF (.NOT.dump_first) print_arrays = .FALSE.
    ENDIF

    iodumpmask(1,:) = iomask

  END SUBROUTINE io_test



  SUBROUTINE average_field(ioutput, avg)

    INTEGER, INTENT(IN) :: ioutput
    TYPE(averaged_data_block) :: avg
    INTEGER :: n_species_local, ispecies, species_sum
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: array

    avg%real_time = avg%real_time + dt
    avg%started = .TRUE.

    species_sum = 0
    n_species_local = 0
    IF (IAND(iomask(ioutput), c_io_no_sum) .EQ. 0) species_sum = 1
    IF (IAND(iomask(ioutput), c_io_species) .NE. 0) n_species_local = n_species

    n_species_local = n_species_local + species_sum

    IF (n_species_local .LE. 0) RETURN

    IF (avg%dump_single) THEN
      SELECT CASE(ioutput)
      CASE(c_dump_ex)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) + REAL(ex * dt, r4)
      CASE(c_dump_ey)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) + REAL(ey * dt, r4)
      CASE(c_dump_ez)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) + REAL(ez * dt, r4)
      CASE(c_dump_bx)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) + REAL(bx * dt, r4)
      CASE(c_dump_by)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) + REAL(by * dt, r4)
      CASE(c_dump_bz)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) + REAL(bz * dt, r4)
      CASE(c_dump_jx)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) &
            + REAL(jx(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng) * dt, r4)
      CASE(c_dump_jy)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) &
            + REAL(jy(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng) * dt, r4)
      CASE(c_dump_jz)
        avg%r4array(:,:,:,1) = avg%r4array(:,:,:,1) &
            + REAL(jz(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng) * dt, r4)
      CASE(c_dump_ekbar)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_ekbar(array, ispecies-species_sum)
          avg%r4array(:,:,:,ispecies) = avg%r4array(:,:,:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_mass_density)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_mass_density(array, ispecies-species_sum)
          avg%r4array(:,:,:,ispecies) = avg%r4array(:,:,:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_charge_density)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_charge_density(array, ispecies-species_sum)
          avg%r4array(:,:,:,ispecies) = avg%r4array(:,:,:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_number_density)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_number_density(array, ispecies-species_sum)
          avg%r4array(:,:,:,ispecies) = avg%r4array(:,:,:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_temperature)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-species_sum)
          avg%r4array(:,:,:,ispecies) = avg%r4array(:,:,:,ispecies) &
              + REAL(array * dt, r4)
        ENDDO
        DEALLOCATE(array)
      END SELECT
    ELSE
      SELECT CASE(ioutput)
      CASE(c_dump_ex)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) + ex * dt
      CASE(c_dump_ey)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) + ey * dt
      CASE(c_dump_ez)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) + ez * dt
      CASE(c_dump_bx)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) + bx * dt
      CASE(c_dump_by)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) + by * dt
      CASE(c_dump_bz)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) + bz * dt
      CASE(c_dump_jx)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) &
            + jx(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng) * dt
      CASE(c_dump_jy)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) &
            + jy(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng) * dt
      CASE(c_dump_jz)
        avg%array(:,:,:,1) = avg%array(:,:,:,1) &
            + jz(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng) * dt
      CASE(c_dump_ekbar)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_ekbar(array, ispecies-species_sum)
          avg%array(:,:,:,ispecies) = avg%array(:,:,:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_mass_density)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_mass_density(array, ispecies-species_sum)
          avg%array(:,:,:,ispecies) = avg%array(:,:,:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_charge_density)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_charge_density(array, ispecies-species_sum)
          avg%array(:,:,:,ispecies) = avg%array(:,:,:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_number_density)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_number_density(array, ispecies-species_sum)
          avg%array(:,:,:,ispecies) = avg%array(:,:,:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      CASE(c_dump_temperature)
        ALLOCATE(array(-2:nx+3,-2:ny+3,-2:nz+3))
        DO ispecies = 1, n_species_local
          CALL calc_temperature(array, ispecies-species_sum)
          avg%array(:,:,:,ispecies) = avg%array(:,:,:,ispecies) + array * dt
        ENDDO
        DEALLOCATE(array)
      END SELECT
    ENDIF

  END SUBROUTINE average_field



  SUBROUTINE energy_account()

  END SUBROUTINE energy_account



  SUBROUTINE write_field(id, code, block_id, name, units, stagger, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:,:,:), INTENT(IN) :: array
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: should_dump, subtype, subarray, io
    LOGICAL :: convert
    TYPE(averaged_data_block), POINTER :: avg

    IF (IAND(iomask(id), code) .EQ. 0) RETURN

    dims = (/nx_global, ny_global, nz_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))
    convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
        .AND. (IAND(code,c_io_restartable) .EQ. 0 &
        .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))

    IF (convert) THEN
      subtype  = subtype_field_r4
      IF (id .EQ. c_dump_jx .OR. id .EQ. c_dump_jy .OR. id .EQ. c_dump_jz) THEN
        subarray = subarray_field_big_r4
      ELSE
        subarray = subarray_field_r4
      ENDIF
    ELSE
      subtype  = subtype_field
      IF (id .EQ. c_dump_jx .OR. id .EQ. c_dump_jy .OR. id .EQ. c_dump_jz) THEN
        subarray = subarray_field_big
      ELSE
        subarray = subarray_field
      ENDIF
    ENDIF

    IF (IAND(iomask(id), should_dump) .NE. 0) THEN
      CALL sdf_write_plain_variable(sdf_handle, TRIM(block_id), &
          TRIM(name), TRIM(units), dims, stagger, 'grid', array, &
          subtype, subarray, convert)
    ENDIF

    DO io = 1, n_io_blocks
      IF (io_block_list(io)%dump) THEN
        avg => io_block_list(io)%averaged_data(id)
        IF (IAND(iomask(id), c_io_averaged) .NE. 0 .AND. avg%started) THEN
          IF (avg%dump_single) THEN
            avg%r4array = avg%r4array / REAL(avg%real_time, r4)

            CALL sdf_write_plain_variable(sdf_handle, &
                TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
                TRIM(units), dims, stagger, 'grid', &
                avg%r4array(:,:,:,1), subtype_field_r4, subarray_field_r4)

            avg%r4array = 0.0_num
          ELSE
            avg%array = avg%array / avg%real_time

            CALL sdf_write_plain_variable(sdf_handle, &
                TRIM(block_id) // '_averaged', TRIM(name) // '_averaged', &
                TRIM(units), dims, stagger, 'grid', &
                avg%array(:,:,:,1), subtype_field, subarray_field)

            avg%array = 0.0_num
          ENDIF

          avg%real_time = 0.0_num
          avg%started = .FALSE.
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE write_field



  SUBROUTINE write_nspecies_field(id, code, block_id, name, units, stagger, &
      func, array)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: array
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: should_dump, subtype, subarray, ispecies, species_sum
    INTEGER :: len1, len2, len3, len4, len5, io
    CHARACTER(LEN=c_id_length) :: temp_block_id
    CHARACTER(LEN=c_max_string_length) :: temp_name, len_string
    LOGICAL :: convert
    TYPE(averaged_data_block), POINTER :: avg

    INTERFACE
      SUBROUTINE func(data_array, current_species)
        USE constants
        REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: current_species
      END SUBROUTINE func
    END INTERFACE

    IF (IAND(iomask(id), code) .EQ. 0) RETURN

    dims = (/nx_global, ny_global, nz_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))
    convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
        .AND. (IAND(code,c_io_restartable) .EQ. 0 &
        .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))

    IF (convert) THEN
      subtype  = subtype_field_r4
      subarray = subarray_field_r4
    ELSE
      subtype  = subtype_field
      subarray = subarray_field
    ENDIF

    IF (IAND(iomask(id), should_dump) .NE. 0) THEN
      IF (IAND(iomask(id), c_io_no_sum) .EQ. 0 &
          .AND. IAND(iomask(id), c_io_field) .EQ. 0) THEN
        CALL build_species_subset

        IF (isubset .EQ. 1) THEN
          temp_block_id = TRIM(block_id)
          temp_name = 'Derived/' // TRIM(name)
        ELSE
          temp_block_id = TRIM(block_id) // '/s_' // &
              TRIM(subset_list(isubset-1)%name)
          temp_name = 'Derived/' // TRIM(name) // '/Subset_' // &
              TRIM(subset_list(isubset-1)%name)
        ENDIF
        CALL func(array, 0)
        CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
            TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
            subtype, subarray, convert)
      ENDIF

      IF (IAND(iomask(id), c_io_species) .NE. 0) THEN
        CALL build_species_subset

        len1 = LEN_TRIM(block_id) + 1
        len2 = LEN_TRIM(name) + 9
        DO ispecies = 1, n_species
          len3 = LEN_TRIM(io_list(ispecies)%name)
          len4 = len3
          len5 = len3
          IF ((len1 + len3) > c_id_length) len4 = c_id_length - len1
          IF ((len2 + len3) > c_max_string_length) THEN
            len5 = c_max_string_length - len2
            IF (rank .EQ. 0) THEN
              CALL integer_as_string((len2+len3), len_string)
              PRINT*, '*** WARNING ***'
              PRINT*, 'Output block name ','Derived/' // TRIM(name) // '/' &
                  // TRIM(io_list(ispecies)%name),' is truncated.'
              PRINT*, 'Either shorten the species name or increase ', &
                  'the size of "c_max_string_length" ', &
                  'to at least ',TRIM(len_string)
            ENDIF
          ENDIF

          temp_block_id = TRIM(block_id) // '/' // &
              TRIM(io_list(ispecies)%name(1:len4))
          temp_name = 'Derived/' // TRIM(name) // '/' // &
              TRIM(io_list(ispecies)%name(1:len5))
          CALL func(array, ispecies)
          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
              subtype, subarray, convert)
        ENDDO
      ENDIF
    ENDIF

    IF (isubset .NE. 1) RETURN

    ! Write averaged data
    DO io = 1, n_io_blocks
      IF (io_block_list(io)%dump) THEN
        avg => io_block_list(io)%averaged_data(id)
        IF (IAND(iomask(id), c_io_averaged) .NE. 0 .AND. avg%started) THEN
          IF (avg%dump_single) THEN
            avg%r4array = avg%r4array / REAL(avg%real_time, r4)

            species_sum = 0
            IF (IAND(iomask(id), c_io_no_sum) .EQ. 0 &
                .AND. IAND(iomask(id), c_io_field) .EQ. 0) THEN
              species_sum = 1
              CALL sdf_write_plain_variable(sdf_handle, &
                  TRIM(block_id) // '_averaged', &
                  'Derived/' // TRIM(name) // '_averaged', &
                  TRIM(units), dims, stagger, 'grid', &
                  avg%r4array(:,:,:,1), subtype_field_r4, subarray_field_r4)
            ENDIF

            IF (IAND(iomask(id), c_io_species) .NE. 0) THEN
              len1 = LEN_TRIM(block_id) + 10
              len2 = LEN_TRIM(name) + 18

              DO ispecies = 1, n_species
                len3 = LEN_TRIM(io_list(ispecies)%name)
                len4 = len3
                len5 = len3
                IF ((len1 + len3) > c_id_length) len4 = c_id_length - len1
                IF ((len2 + len3) > c_max_string_length) THEN
                  len5 = c_max_string_length - len2
                  IF (rank .EQ. 0) THEN
                    CALL integer_as_string((len2+len3), len_string)
                    PRINT*, '*** WARNING ***'
                    PRINT*, 'Output block name ','Derived/' // TRIM(name) // &
                        '_averaged/' // TRIM(io_list(ispecies)%name), &
                        ' is truncated.'
                    PRINT*, 'Either shorten the species name or increase ', &
                        'the size of "c_max_string_length" ', &
                        'to at least ',TRIM(len_string)
                  ENDIF
                ENDIF

                temp_block_id = TRIM(block_id) // '_averaged/' // &
                    TRIM(io_list(ispecies)%name(1:len4))
                temp_name = 'Derived/' // TRIM(name) // '_averaged/' // &
                    TRIM(io_list(ispecies)%name(1:len5))
                CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                    TRIM(temp_name), TRIM(units), dims, stagger, 'grid', &
                    avg%r4array(:,:,:,ispecies+species_sum), &
                    subtype_field_r4, subarray_field_r4)
              ENDDO
            ENDIF

            avg%r4array = 0.0_num
          ELSE
            avg%array = avg%array / avg%real_time

            species_sum = 0
            IF (IAND(iomask(id), c_io_no_sum) .EQ. 0 &
                .AND. IAND(iomask(id), c_io_field) .EQ. 0) THEN
              species_sum = 1
              CALL sdf_write_plain_variable(sdf_handle, &
                  TRIM(block_id) // '_averaged', &
                  'Derived/' // TRIM(name) // '_averaged', &
                  TRIM(units), dims, stagger, 'grid', &
                  avg%array(:,:,:,1), subtype_field, subarray_field)
            ENDIF

            IF (IAND(iomask(id), c_io_species) .NE. 0) THEN
              len1 = LEN_TRIM(block_id) + 10
              len2 = LEN_TRIM(name) + 18

              DO ispecies = 1, n_species
                len3 = LEN_TRIM(io_list(ispecies)%name)
                len4 = len3
                len5 = len3
                IF ((len1 + len3) > c_id_length) len4 = c_id_length - len1
                IF ((len2 + len3) > c_max_string_length) THEN
                  len5 = c_max_string_length - len2
                  IF (rank .EQ. 0) THEN
                    CALL integer_as_string((len2+len3), len_string)
                    PRINT*, '*** WARNING ***'
                    PRINT*, 'Output block name ','Derived/' // TRIM(name) // &
                        '_averaged/' // TRIM(io_list(ispecies)%name), &
                        ' is truncated.'
                    PRINT*, 'Either shorten the species name or increase ', &
                        'the size of "c_max_string_length" ', &
                        'to at least ',TRIM(len_string)
                  ENDIF
                ENDIF

                temp_block_id = TRIM(block_id) // '_averaged/' // &
                    TRIM(io_list(ispecies)%name(1:len4))
                temp_name = 'Derived/' // TRIM(name) // '_averaged/' // &
                    TRIM(io_list(ispecies)%name(1:len5))
                CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                    TRIM(temp_name), TRIM(units), dims, stagger, 'grid', &
                    avg%array(:,:,:,ispecies+species_sum), &
                    subtype_field, subarray_field)
              ENDDO
            ENDIF

            avg%array = 0.0_num
          ENDIF

          avg%real_time = 0.0_num
          avg%started = .FALSE.
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE write_nspecies_field



  SUBROUTINE write_field_flux(id, code, block_id, name, units, stagger, &
      func, array, fluxdir, dir_tags)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: array
    INTEGER, DIMENSION(:), INTENT(IN) :: fluxdir
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: dir_tags
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: should_dump, subtype, subarray, ndirs, idir
    CHARACTER(LEN=c_id_length) :: temp_block_id
    CHARACTER(LEN=c_max_string_length) :: temp_name
    LOGICAL :: convert

    INTERFACE
      SUBROUTINE func(data_array, direction)
        USE constants
        REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: direction
      END SUBROUTINE func
    END INTERFACE

    IF (IAND(iomask(id), code) .EQ. 0) RETURN

    ndirs = SIZE(fluxdir)
    dims = (/nx_global, ny_global, nz_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))
    convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
        .AND. (IAND(code,c_io_restartable) .EQ. 0 &
        .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))

    IF (convert) THEN
      subtype  = subtype_field_r4
      subarray = subarray_field_r4
    ELSE
      subtype  = subtype_field
      subarray = subarray_field
    ENDIF

    IF (IAND(iomask(id), should_dump) .NE. 0) THEN
      DO idir = 1, ndirs
        temp_block_id = TRIM(block_id) // '/' // &
            TRIM(dir_tags(idir))
        temp_name = 'Derived/' // TRIM(name) // '/' // &
            TRIM(dir_tags(idir))
        CALL func(array, fluxdir(idir))
        CALL sdf_write_plain_variable(sdf_handle, &
            TRIM(ADJUSTL(temp_block_id)), TRIM(ADJUSTL(temp_name)), &
            TRIM(units), dims, stagger, 'grid', &
            array, subtype, subarray, convert)
      ENDDO
    ENDIF

    ! Flux variables not currently averaged

  END SUBROUTINE write_field_flux



  SUBROUTINE write_nspecies_flux(id, code, block_id, name, units, stagger, &
      func, array, fluxdir, dir_tags)

    INTEGER, INTENT(IN) :: id, code
    CHARACTER(LEN=*), INTENT(IN) :: block_id, name, units
    INTEGER, INTENT(IN) :: stagger
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: array
    INTEGER, DIMENSION(:), INTENT(IN) :: fluxdir
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: dir_tags
    INTEGER, DIMENSION(c_ndims) :: dims
    INTEGER :: should_dump, subtype, subarray, ispecies, ndirs, idir
    INTEGER :: len1, len2, len3, len4, len5
    CHARACTER(LEN=c_id_length) :: temp_block_id
    CHARACTER(LEN=c_max_string_length) :: temp_name, len_string
    LOGICAL :: convert

    INTERFACE
      SUBROUTINE func(data_array, current_species, direction)
        USE constants
        REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: data_array
        INTEGER, INTENT(IN) :: current_species, direction
      END SUBROUTINE func
    END INTERFACE

    IF (IAND(iomask(id), code) .EQ. 0) RETURN

    ndirs = SIZE(fluxdir)
    dims = (/nx_global, ny_global, nz_global/)

    ! Want the code to output unaveraged data if either restarting or
    ! requested by the user
    should_dump = IOR(c_io_snapshot, IAND(code,c_io_restartable))
    should_dump = IOR(should_dump, NOT(c_io_averaged))
    convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
        .AND. (IAND(code,c_io_restartable) .EQ. 0 &
        .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))

    IF (convert) THEN
      subtype  = subtype_field_r4
      subarray = subarray_field_r4
    ELSE
      subtype  = subtype_field
      subarray = subarray_field
    ENDIF

    IF (IAND(iomask(id), should_dump) .NE. 0) THEN
      IF (IAND(iomask(id), c_io_no_sum) .EQ. 0 &
          .AND. IAND(iomask(id), c_io_field) .EQ. 0) THEN
        CALL build_species_subset

        DO idir = 1, ndirs
          temp_block_id = TRIM(block_id) // '/' // &
              TRIM(dir_tags(idir))
          temp_name = 'Derived/' // TRIM(name) // '/' // &
              TRIM(dir_tags(idir))
          CALL func(array, 0, fluxdir(idir))
          CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
              TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
              subtype, subarray, convert)
        ENDDO
      ENDIF

      IF (IAND(iomask(id), c_io_species) .NE. 0) THEN
        CALL build_species_subset

        len1 = LEN_TRIM(block_id) + LEN_TRIM(dir_tags(1)) + 2
        len2 = LEN_TRIM(name) + LEN_TRIM(dir_tags(1)) + 10

        DO ispecies = 1, n_species
          len3 = LEN_TRIM(io_list(ispecies)%name)
          len4 = len3
          len5 = len3
          IF ((len1 + len3) > c_id_length) len4 = c_id_length - len1
          IF ((len2 + len3) > c_max_string_length) THEN
            len5 = c_max_string_length - len2
            IF (rank .EQ. 0) THEN
              CALL integer_as_string((len2+len3), len_string)
              PRINT*, '*** WARNING ***'
              PRINT*, 'Output block name ','Derived/' // TRIM(name) // '_' &
                  // TRIM(io_list(ispecies)%name) // '/' &
                  // TRIM(dir_tags(1)),' is truncated.'
              PRINT*, 'Either shorten the species name or increase ', &
                  'the size of "c_max_string_length" ', &
                  'to at least ',TRIM(len_string)
            ENDIF
          ENDIF

          DO idir = 1, ndirs
            temp_block_id = TRIM(block_id) // '_' // &
                TRIM(dir_tags(idir)) // '/' // &
                TRIM(io_list(ispecies)%name(1:len4))
            temp_name = 'Derived/' // TRIM(name) // '_' // &
                TRIM(dir_tags(idir)) // '/' // &
                TRIM(io_list(ispecies)%name(1:len5))
            CALL func(array, ispecies, fluxdir(idir))
            CALL sdf_write_plain_variable(sdf_handle, TRIM(temp_block_id), &
                TRIM(temp_name), TRIM(units), dims, stagger, 'grid', array, &
                subtype, subarray, convert)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    ! Flux variables not currently averaged

  END SUBROUTINE write_nspecies_flux



  SUBROUTINE build_species_subset

    INTEGER :: i, l
    TYPE(particle), POINTER :: current, next
    LOGICAL :: use_particle
    REAL(num) :: gamma, random_num, part_mc

    IF (done_subset_init) RETURN
    done_subset_init = .TRUE.

    IF (isubset .EQ. 1) THEN
      io_list => species_list
      RETURN
    ENDIF

    io_list => io_list_data

    l = isubset - 1
    DO i = 1, n_species
      io_list(i) = species_list(i)
      io_list(i)%count = 0
      io_list(i)%name = 'subset_' // TRIM(subset_list(l)%name) // '/' // &
          TRIM(species_list(i)%name)
      CALL create_empty_partlist(io_list(i)%attached_list)

      IF (.NOT. subset_list(l)%use_species(i)) THEN
        io_list(i)%dumpmask = c_io_never
        CYCLE
      ENDIF

      part_mc = c * species_list(i)%mass

      current => species_list(i)%attached_list%head
      DO WHILE (ASSOCIATED(current))
        next => current%next
        use_particle = .TRUE.
        IF (subset_list(l)%use_gamma) THEN
#ifdef PER_PARTICLE_CHARGE_MASS
          part_mc = c * current%mass
#endif
          gamma = SQRT(SUM((current%part_p / part_mc)**2) + 1.0_num)
          IF (subset_list(l)%use_gamma_min &
              .AND. gamma .LT. subset_list(l)%gamma_min) use_particle = .FALSE.
          IF (subset_list(l)%use_gamma_max &
              .AND. gamma .GT. subset_list(l)%gamma_max) use_particle = .FALSE.
        ENDIF

        IF (subset_list(l)%use_x_min &
            .AND. current%part_pos(1) .LT. subset_list(l)%x_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_x_max &
            .AND. current%part_pos(1) .GT. subset_list(l)%x_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_y_min &
            .AND. current%part_pos(2) .LT. subset_list(l)%y_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_y_max &
            .AND. current%part_pos(2) .GT. subset_list(l)%y_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_z_min &
            .AND. current%part_pos(3) .LT. subset_list(l)%z_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_z_max &
            .AND. current%part_pos(3) .GT. subset_list(l)%z_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_px_min &
            .AND. current%part_p(1) .LT. subset_list(l)%px_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_px_max &
            .AND. current%part_p(1) .GT. subset_list(l)%px_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_py_min &
            .AND. current%part_p(2) .LT. subset_list(l)%py_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_py_max &
            .AND. current%part_p(2) .GT. subset_list(l)%py_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_pz_min &
            .AND. current%part_p(3) .LT. subset_list(l)%pz_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_pz_max &
            .AND. current%part_p(3) .GT. subset_list(l)%pz_max) &
                use_particle = .FALSE.

#ifdef PER_PARTICLE_WEIGHT
        IF (subset_list(l)%use_weight_min &
            .AND. current%weight .LT. subset_list(l)%weight_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_weight_max &
            .AND. current%weight .GT. subset_list(l)%weight_max) &
                use_particle = .FALSE.

#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        IF (subset_list(l)%use_charge_min &
            .AND. current%charge .LT. subset_list(l)%charge_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_charge_max &
            .AND. current%charge .GT. subset_list(l)%charge_max) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_mass_min &
            .AND. current%mass .LT. subset_list(l)%mass_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_mass_max &
            .AND. current%mass .GT. subset_list(l)%mass_max) &
                use_particle = .FALSE.

#endif
#if PARTICLE_ID || PARTICLE_ID4
        IF (subset_list(l)%use_id_min &
            .AND. current%id .LT. subset_list(l)%id_min) &
                use_particle = .FALSE.

        IF (subset_list(l)%use_id_max &
            .AND. current%id .LT. subset_list(l)%id_max) &
                use_particle = .FALSE.
#endif

        IF (subset_list(l)%use_random) THEN
          random_num = random()
          IF (random_num .GT. subset_list(l)%random_fraction) &
              use_particle = .FALSE.
        ENDIF

        IF (use_particle) THEN
          ! Move particle to io_list
          CALL remove_particle_from_partlist(species_list(i)%attached_list, &
              current)
          CALL add_particle_to_partlist(io_list(i)%attached_list, current)
        ENDIF
        current => next
      ENDDO
    ENDDO

  END SUBROUTINE build_species_subset



  SUBROUTINE species_offset_init

    INTEGER(i8) :: species_count
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: npart_species_per_proc
    INTEGER :: i, ispecies
    TYPE(particle_species), POINTER :: spec

    IF (done_species_offset_init) RETURN
    done_species_offset_init = .TRUE.

    IF (.NOT.ALLOCATED(species_offset)) ALLOCATE(species_offset(n_species))

    CALL build_species_subset

    ALLOCATE(npart_species_per_proc(nproc))
    species_offset = 0

    npart_global = 0
    DO ispecies = 1, n_species
      spec => io_list(ispecies)

      CALL MPI_ALLGATHER(spec%attached_list%count, 1, MPI_INTEGER8, &
          npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
      species_count = 0
      DO i = 1, nproc
        IF (rank .EQ. i-1) species_offset(ispecies) = species_count
        species_count = species_count + npart_species_per_proc(i)
      ENDDO
      spec%count = species_count
      npart_global = npart_global + species_count

      CALL sdf_write_cpu_split(sdf_handle, 'cpu/' // TRIM(spec%name), &
          'CPU split/' // TRIM(spec%name), npart_species_per_proc)
    ENDDO

    IF (track_ejected_particles &
        .AND. .NOT.ALLOCATED(ejected_offset)) THEN
      ALLOCATE(ejected_offset(n_species))
      ejected_offset = 0

      DO ispecies = 1, n_species
        spec => ejected_list(ispecies)

        CALL MPI_ALLGATHER(spec%attached_list%count, 1, MPI_INTEGER8, &
            npart_species_per_proc, 1, MPI_INTEGER8, comm, errcode)
        species_count = 0
        DO i = 1, nproc
          IF (rank .EQ. i-1) ejected_offset(ispecies) = species_count
          species_count = species_count + npart_species_per_proc(i)
        ENDDO
        spec%count = species_count

        CALL sdf_write_cpu_split(sdf_handle, 'cpu/' // TRIM(spec%name), &
            'CPU split/' // TRIM(spec%name), npart_species_per_proc)
      ENDDO
    ENDIF

    DEALLOCATE(npart_species_per_proc)

  END SUBROUTINE species_offset_init



  SUBROUTINE write_particle_grid(code)

    INTEGER, INTENT(IN) :: code
    INTEGER :: ispecies, id
    LOGICAL :: convert

    id = c_dump_part_grid
    IF (IAND(iomask(id), code) .NE. 0) THEN
      CALL build_species_subset

      convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
          .AND. (IAND(code,c_io_restartable) .EQ. 0 &
          .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) .NE. 0 &
            .OR. IAND(code, c_io_restartable) .NE. 0) THEN
          CALL species_offset_init()
          IF (npart_global .EQ. 0) RETURN

          CALL sdf_write_point_mesh(sdf_handle, &
              'grid/' // TRIM(current_species%name), &
              'Grid/Particles/' // TRIM(current_species%name), &
              io_list(ispecies)%count, c_dimension_3d, &
              iterate_particles, species_offset(ispecies), convert)
        ENDIF
      ENDDO
    ENDIF

    IF (isubset .NE. 1) RETURN

    id = c_dump_ejected_particles
    IF (IAND(iomask(id), code) .NE. 0) THEN
      convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
          .AND. (IAND(code,c_io_restartable) .EQ. 0 &
          .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global .EQ. 0) RETURN

        current_species => ejected_list(ispecies)
        CALL sdf_write_point_mesh(sdf_handle, &
            'grid/' // TRIM(current_species%name), &
            'Grid/Particles/' // TRIM(current_species%name), &
            ejected_list(ispecies)%count, c_dimension_3d, &
            iterate_particles, ejected_offset(ispecies), convert)
      ENDDO
    ENDIF

  END SUBROUTINE write_particle_grid



  SUBROUTINE write_particle_variable(id_in, code, name, units, iterator)

    INTEGER, INTENT(IN) :: id_in, code
    CHARACTER(LEN=*), INTENT(IN) :: name, units
    CHARACTER(LEN=c_id_length) :: temp_block_id

    INTERFACE
      FUNCTION iterator(array, npart_it, start)
        USE constants
        REAL(num) :: iterator
        REAL(num), DIMENSION(:), INTENT(OUT) :: array
        INTEGER, INTENT(INOUT) :: npart_it
        LOGICAL, INTENT(IN) :: start
      END FUNCTION iterator
    END INTERFACE

    INTEGER :: ispecies, id
    LOGICAL :: convert, found

    id = id_in
    IF (IAND(iomask(id), code) .NE. 0) THEN
      CALL build_species_subset

      convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
          .AND. (IAND(code,c_io_restartable) .EQ. 0 &
          .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))

      DO ispecies = 1, n_species
        current_species => io_list(ispecies)

        IF (IAND(current_species%dumpmask, code) .NE. 0 &
            .OR. IAND(code, c_io_restartable) .NE. 0) THEN
          CALL species_offset_init()
          IF (npart_global .EQ. 0) RETURN

          found = sdf_get_block_id(sdf_handle, &
              'grid/' // TRIM(current_species%name), temp_block_id)
          CALL sdf_write_point_variable(sdf_handle, &
              lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
              'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
              TRIM(units), io_list(ispecies)%count, temp_block_id, &
              iterator, species_offset(ispecies), convert)
        ENDIF
      ENDDO
    ENDIF

    id = c_dump_ejected_particles
    IF (IAND(iomask(id), code) .NE. 0) THEN
      convert = (IAND(iomask(id), c_io_dump_single) .NE. 0 &
          .AND. (IAND(code,c_io_restartable) .EQ. 0 &
          .OR. IAND(iomask(id), c_io_restartable) .EQ. 0))
      reset_ejected = .TRUE.

      DO ispecies = 1, n_species
        CALL species_offset_init()
        IF (npart_global .EQ. 0) RETURN

        current_species => ejected_list(ispecies)
        found = sdf_get_block_id(sdf_handle, &
            'grid/' // TRIM(current_species%name), temp_block_id)
        CALL sdf_write_point_variable(sdf_handle, &
            lowercase(TRIM(name) // '/' // TRIM(current_species%name)), &
            'Particles/' // TRIM(name) // '/' // TRIM(current_species%name), &
            TRIM(units), ejected_list(ispecies)%count, temp_block_id, &
            iterator, ejected_offset(ispecies), convert)
      ENDDO
    ENDIF

  END SUBROUTINE write_particle_variable



  FUNCTION lowercase(string_in) RESULT(string_out)

    CHARACTER(LEN=*), PARAMETER :: lwr = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER(LEN=*), PARAMETER :: upr = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(LEN=*), INTENT(IN) :: string_in
    CHARACTER(LEN=LEN(string_in)) :: string_out
    INTEGER :: i, idx

    string_out = string_in

    DO i = 1, LEN(string_out)
      idx = INDEX(upr, string_out(i:i))
      IF (idx .NE. 0) string_out(i:i) = lwr(idx:idx)
    ENDDO

  END FUNCTION lowercase

END MODULE diagnostics
