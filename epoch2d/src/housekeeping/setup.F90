MODULE setup

  USE input_cartesian
  USE input_particle
  USE iocontrol
  USE strings
  USE partlist
  USE mpi_subtype_control
  USE fields
  USE version_data
  USE encoded_source

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control, minimal_init, restart_data
  PUBLIC :: open_files, close_files
  PUBLIC :: setup_species

  SAVE
  TYPE(particle_list) :: main_root
  INTEGER, DIMENSION(:), ALLOCATABLE :: species_id

CONTAINS

  SUBROUTINE minimal_init

    REAL(num) :: dummy = 1.0_num

    IF (num .EQ. 4) mpireal = MPI_REAL
    dumpmask = 0
    comm = MPI_COMM_NULL
    CALL MPI_SIZEOF(dummy, realsize, errcode)

    dt_plasma_frequency = 0.0_num
    stdout_frequency = 0

    window_shift = 0.0_num
    npart_global = -1
    smooth_currents = .FALSE.
    dlb = .FALSE.

    NULLIFY(laser_x_min)
    NULLIFY(laser_x_max)
    NULLIFY(laser_y_max)
    NULLIFY(laser_y_min)

    NULLIFY(dist_fns)

    run_date = get_unix_time()

    CALL set_field_order(2)

    CALL init_source_code()

  END SUBROUTINE minimal_init



  SUBROUTINE after_control

    INTEGER :: iproc, ix, iy, ix_min, iy_min

    length_x = x_max - x_min
    length_y = y_max - y_min
    dx = length_x / REAL(nx_global-1, num)
    dy = length_y / REAL(ny_global-1, num)

    ! Setup global grid
    DO ix = -1, nx_global + 2
      x_global(ix) = x_min + (ix - 1) * dx
      x_offset_global(ix) = x_global(ix)
    ENDDO
    DO iy = -1, ny_global + 2
      y_global(iy) = y_min + (iy - 1) * dy
      y_offset_global(iy) = y_global(iy)
    ENDDO

    DO iproc = 0, nprocx-1
      x_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocy-1
      y_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO

    x_min_local = x_mins(coordinates(2))
    x_max_local = x_maxs(coordinates(2))
    y_min_local = y_mins(coordinates(1))
    y_max_local = y_maxs(coordinates(1))

    ix_min = cell_x_min(coordinates(2)+1)
    iy_min = cell_y_min(coordinates(1)+1)

    ! Setup local grid
    DO ix = -1, nx+2
      x(ix) = x_global(ix_min+ix-1)
    ENDDO
    DO iy = -1, ny+2
      y(iy) = y_global(iy_min+iy-1)
    ENDDO

    CALL set_initial_values

    CALL setup_data_averaging

  END SUBROUTINE after_control



  SUBROUTINE setup_data_averaging()

    INTEGER :: ioutput

    DO ioutput = 1, num_vars_to_dump
      IF (IAND(dumpmask(ioutput), c_io_averaged) .NE. 0) THEN
        ALLOCATE(averaged_data(ioutput)%data(-2:nx+3, -2:ny+3))
      ENDIF
    ENDDO

  END SUBROUTINE setup_data_averaging



  SUBROUTINE setup_species

    INTEGER :: ispecies

    ALLOCATE(particle_species(n_species))
    ALLOCATE(particle_file_lengths(n_species))
    ALLOCATE(particle_file_offsets(n_species))

    DO ispecies = 1, n_species
      particle_species(ispecies)%name = blank
      particle_species(ispecies)%mass = -1.0_num
      particle_species(ispecies)%dump = .TRUE.
      particle_species(ispecies)%count = -1
#ifdef SPLIT_PARTICLES_AFTER_PUSH
      particle_species(ispecies)%split = .FALSE.
      particle_species(ispecies)%npart_max = 0.0_num
#endif
#ifdef PARTICLE_IONISE
      particle_species(ispecies)%ionise = .FALSE.
      particle_species(ispecies)%ionise_to_species = -1
      particle_species(ispecies)%release_species = -1
      particle_species(ispecies)%critical_field = 0.0_num
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
      WRITE(file2, '(a, "/epoch2d.dat")') TRIM(data_dir)
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
        PRINT *, '***ERROR***'
        PRINT *, 'Cannot create "epoch2d.dat" output file. The most common ' &
            // 'cause of this problem '
        PRINT *, 'is that the ouput directory does not exist'
        CALL MPI_ABORT(comm, errcode, ierr)
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



  SUBROUTINE restart_data(snap)

    INTEGER, INTENT(OUT) :: snap
    CHARACTER(LEN=20+data_dir_max_length) :: filename
    CHARACTER(LEN=max_string_len) :: name, class, mesh_name, mesh_class
    INTEGER :: block_type, nd
    INTEGER :: sof
    INTEGER(KIND=8) :: npart_l
    REAL(num), DIMENSION(2) :: extents, stagger
    INTEGER, DIMENSION(1) :: dims
    REAL(KIND=8) :: time_d
    INTEGER :: coord_type, ierr, nblocks
    TYPE(particle), POINTER :: current, next
    LOGICAL :: constant_weight
    INTEGER(KIND=8) :: npart, ipart, ix

    npart_global = 0
    constant_weight = .FALSE.
    snap = -1

    ! Create the filename for the last snapshot
    WRITE(filename, '(a, "/", i4.4, ".cfd")') TRIM(data_dir), restart_snapshot
    CALL cfd_open(filename, rank, comm, MPI_MODE_RDONLY, snap, time_d)
    IF (snap .GE. 0) THEN
      time = time_d
      IF (rank .EQ. 0) PRINT *, "Loading snapshot for time", time
    ENDIF
    ! open the file
    nblocks = cfd_get_nblocks()
    jobid = cfd_get_jobid()

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    IF (rank .EQ. 0) PRINT *, "Input file contains", nblocks, "blocks"
    DO ix = 1, nblocks
      CALL cfd_get_next_block_info_all(name, class, block_type)
      ! IF (rank .EQ. 0) PRINT *, "Loading block", ix, name, block_type
      IF (block_type .EQ. c_type_snapshot .AND. snap .LT. 0) THEN
        CALL cfd_get_snapshot(time_d, snap)
        time = time_d
        IF (rank .EQ. 0) PRINT *, "Loading snapshot for time", time
      ENDIF

      SELECT CASE(block_type)
      CASE(c_type_mesh_variable)
        CALL cfd_get_common_meshtype_metadata_all(block_type, nd, sof)
        IF (sof .NE. num) THEN
          IF (rank .EQ. 0) &
              PRINT *, "Precision does not match, recompile code so &
                  &that sizeof(REAL) = ", sof
          CALL MPI_ABORT(comm, errcode, ierr)
        ENDIF

        IF (nd .NE. c_dimension_2d .AND. nd .NE. c_dimension_irrelevant ) THEN
          IF (rank .EQ. 0) &
              PRINT *, "Dimensionality does not match, file is ", nd, "D"
          CALL MPI_ABORT(comm, errcode, ierr)
        ENDIF

        SELECT CASE(block_type)
        CASE(c_var_cartesian)
          ! Grid variables
          CALL cfd_get_nd_cartesian_variable_metadata_all(nd, dims, extents, &
              stagger, mesh_name, mesh_class)

          IF (dims(1) .NE. nx_global) THEN
            IF (rank .EQ. 0) &
                PRINT *, "Number of gridpoints does not match, gridpoints &
                    &in file is", dims(1)
            CALL MPI_ABORT(comm, errcode, ierr)
          ENDIF

          IF (str_cmp(name, "Ex")) &
              CALL cfd_get_2d_cartesian_variable_parallel(ex, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "Ey")) &
              CALL cfd_get_2d_cartesian_variable_parallel(ey, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "Ez")) &
              CALL cfd_get_2d_cartesian_variable_parallel(ez, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "Bx")) &
              CALL cfd_get_2d_cartesian_variable_parallel(bx, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "By")) &
              CALL cfd_get_2d_cartesian_variable_parallel(by, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "Bz")) &
              CALL cfd_get_2d_cartesian_variable_parallel(bz, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "Jx")) &
              CALL cfd_get_2d_cartesian_variable_parallel(jx, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "Jy")) &
              CALL cfd_get_2d_cartesian_variable_parallel(jy, subtype_field, &
                  subarray_field)

          IF (str_cmp(name, "Jz")) &
              CALL cfd_get_2d_cartesian_variable_parallel(jz, subtype_field, &
                  subarray_field)

        CASE(c_var_particle)
          CALL cfd_get_nd_particle_variable_metadata_all(npart_l, extents, &
              mesh_name, mesh_class)

          IF (npart_l .NE. npart_global) THEN
            IF (rank .EQ. 0) THEN
              PRINT *, "*** ERROR ***"
              PRINT *, "Malformed restart dump. Number of particle variables", &
                  " does not match grid."
            ENDIF
            CALL MPI_ABORT(comm, errcode, ierr)
            STOP
          ENDIF

          ! particle variables
          IF (str_cmp(name, "Px")) &
              CALL cfd_get_nd_particle_variable_parallel_with_iterator(npart, &
                  npart_per_it, subtype_particle_var, it_px)

          IF (str_cmp(name, "Py")) &
              CALL cfd_get_nd_particle_variable_parallel_with_iterator(npart, &
                  npart_per_it, subtype_particle_var, it_py)

          IF (str_cmp(name, "Pz")) &
              CALL cfd_get_nd_particle_variable_parallel_with_iterator(npart, &
                  npart_per_it, subtype_particle_var, it_pz)

#ifdef PER_PARTICLE_WEIGHT
          IF (str_cmp(name, "Weight")) &
              CALL cfd_get_nd_particle_variable_parallel_with_iterator(npart, &
                  npart_per_it, subtype_particle_var, it_weight)
#else
          IF (str_cmp(name, "Weight")) THEN
            IF (rank .EQ. 0) &
                PRINT *, "Cannot load dump file with per particle weight &
                    &if the code is compiled without per particle weights. &
                    &Code terminates"
            CALL MPI_ABORT(comm, errcode, ierr)
          ENDIF
#endif
          IF (str_cmp(name, "Species")) &
              CALL cfd_get_nd_particle_variable_parallel_with_iterator(npart, &
                  npart_per_it, subtype_particle_var, it_species)
        END SELECT
      CASE(c_type_mesh)
        CALL cfd_get_common_meshtype_metadata_all(block_type, nd, sof)
        IF (block_type .EQ. c_mesh_particle) THEN
          CALL cfd_get_nd_particle_grid_metadata_all(nd, coord_type, npart_l, &
              extents)
          IF (npart_l .NE. npart_global) THEN
            npart = npart_l/nproc
            IF (npart * nproc .NE. npart_l) THEN
              IF (rank .EQ. 0) &
                  PRINT *, "Cannot evenly subdivide particles over", nproc, &
                      "processors. Trying to fix"
              IF (rank .LT. npart_l-(npart*nproc)) npart = npart+1
            ENDIF
            CALL create_subtypes_for_load(npart)
            CALL create_allocated_partlist(main_root, npart)
            ALLOCATE(species_id(npart))
            current=>main_root%head
            npart_global = npart_l
          ENDIF
          CALL cfd_get_nd_particle_grid_parallel_with_iterator(nd, &
              main_root%count, npart_l, npart_per_it, sof, &
              subtype_particle_var, it_part)
        ENDIF
      CASE(c_type_constant)
        CALL cfd_get_real_constant(weight)
        constant_weight = .TRUE.
      END SELECT
      CALL cfd_skip_block()
    ENDDO
    CALL cfd_close()

    current=>main_root%head
    ipart = 1
    DO WHILE(ASSOCIATED(current))
      next=>current%next
      CALL remove_particle_from_partlist(main_root, current)
      CALL add_particle_to_partlist(&
          particle_species(species_id(ipart))%attached_list, current)
      current=>next
      ipart = ipart+1
    ENDDO

    DEALLOCATE(species_id)

#ifdef PER_PARTICLE_WEIGHT
    IF (constant_weight) THEN
      current=>main_root%head
      DO WHILE(ASSOCIATED(current))
        current%weight = weight
        current=>current%next
      ENDDO
    ENDIF
#endif

  END SUBROUTINE restart_data



  SUBROUTINE it_part(data, npart_this_it, start, direction)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(KIND=8), INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER, INTENT(IN) :: direction

    INTEGER(KIND=8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) THEN
      cur=>main_root%head
    ENDIF
    DO ipart = 1, npart_this_it
      cur%part_pos(direction) = data(ipart)
      cur=>cur%next
    ENDDO

  END SUBROUTINE it_part



  SUBROUTINE it_px(data, npart_this_it, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(KIND=8), INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER(KIND=8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur=>main_root%head
    DO ipart = 1, npart_this_it
      cur%part_p(1) = data(ipart)
      cur=>cur%next
    ENDDO

  END SUBROUTINE it_px



  SUBROUTINE it_py(data, npart_this_it, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(KIND=8), INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER(KIND=8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur=>main_root%head
    DO ipart = 1, npart_this_it
      cur%part_p(2) = data(ipart)
      cur=>cur%next
    ENDDO

  END SUBROUTINE it_py



  SUBROUTINE it_pz(data, npart_this_it, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(KIND=8), INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER(KIND=8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur=>main_root%head
    DO ipart = 1, npart_this_it
      cur%part_p(3) = data(ipart)
      cur=>cur%next
    ENDDO

  END SUBROUTINE it_pz



#ifdef PER_PARTICLE_WEIGHT
  SUBROUTINE it_weight(data, npart_this_it, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(KIND=8), INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER(KIND=8) :: ipart
    TYPE(particle), POINTER, SAVE :: cur

    IF (start) cur=>main_root%head
    DO ipart = 1, npart_this_it
      cur%weight = data(ipart)
      cur=>cur%next
    ENDDO

  END SUBROUTINE it_weight
#endif



  SUBROUTINE it_species(data, npart_this_it, start)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(KIND=8), INTENT(INOUT) :: npart_this_it
    LOGICAL, INTENT(IN) :: start
    INTEGER(KIND=8) :: ipart
    INTEGER(KIND=8), SAVE :: ipart_total

    IF (start) ipart_total = 1
    DO ipart = 1, npart_this_it
      species_id(ipart_total) = NINT(data(ipart))
      ipart_total = ipart_total+1
    ENDDO

  END SUBROUTINE it_species

END MODULE setup
