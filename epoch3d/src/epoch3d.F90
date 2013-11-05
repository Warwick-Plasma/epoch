PROGRAM pic

  ! EPOCH3D is a Birdsall and Langdon type PIC code derived from the PSC
  ! written by Hartmut Ruhl.

  ! The particle pusher (particles.F90) and the field solver (fields.f90) are
  ! almost exact copies of the equivalent routines from PSC, modified slightly
  ! to allow interaction with the changed portions of the code and for
  ! readability. The MPI routines are exactly equivalent to those in PSC, but
  ! are completely rewritten in a form which is easier to extend with arbitrary
  ! fields and particle properties. The support code is entirely new and is not
  ! equivalent to PSC.

  ! EPOCH3D written by C.S.Brady, Centre for Fusion, Space and Astrophysics,
  ! University of Warwick, UK
  ! PSC written by Hartmut Ruhl

  USE balance
  USE deck
  USE diagnostics
  USE fields
  USE helper
  USE ic_module
  USE mpi_routines
  USE particles
  USE setup
  USE welcome
  USE window
  USE split_particle
  USE collisions
  USE particle_migration
  USE ionise
#ifdef PHOTONS
  USE photons
#endif

  IMPLICIT NONE

  INTEGER :: ispecies, ierr
  LOGICAL :: halt = .FALSE., push = .TRUE.
  LOGICAL :: force_dump
  CHARACTER(LEN=64) :: deck_file = 'input.deck'
  CHARACTER(LEN=*), PARAMETER :: data_dir_file = 'USE_DATA_DIRECTORY'
  CHARACTER(LEN=*), PARAMETER :: stop_file = 'STOP'
  CHARACTER(LEN=*), PARAMETER :: stop_file_nodump = 'STOP_NODUMP'
  CHARACTER(LEN=64) :: timestring
  REAL(num) :: runtime, real_walltime_start

  step = 0
#ifdef COLLISIONS_TEST
  ! used for testing
  CALL test_collisions
  STOP
#endif

  CALL mpi_minimal_init ! mpi_routines.f90
  real_walltime_start = MPI_WTIME()
  CALL minimal_init     ! setup.f90
  CALL setup_partlists  ! partlist.f90
  CALL get_job_id(jobid)
  CALL welcome_message  ! welcome.f90
  CALL register_objects ! custom.f90

  IF (rank .EQ. 0) THEN
    OPEN(unit=lu, status='OLD', file=TRIM(data_dir_file), iostat=ierr)
    IF (ierr .EQ. 0) THEN
      READ(lu,'(A)') data_dir
      CLOSE(lu)
      PRINT*, 'Using data directory "' // TRIM(data_dir) // '"'
    ELSE
      PRINT*, 'Specify output directory'
      READ(*,'(A)') data_dir
    ENDIF
    CALL cleanup_stop_files
  ENDIF

  CALL MPI_BCAST(data_dir, 64, MPI_CHARACTER, 0, comm, errcode)
  CALL read_deck(deck_file, .TRUE., c_ds_first)
  IF (use_exact_restart) CALL read_cpu_split
  CALL setup_particle_boundaries ! boundary.f90
  CALL mpi_initialise  ! mpi_routines.f90
  CALL after_control   ! setup.f90
  CALL open_files      ! setup.f90

  ! Re-scan the input deck for items which require allocated memory
  CALL read_deck(deck_file, .TRUE., c_ds_last)
  CALL after_deck_last

  ! restart flag is set
  IF (ic_from_restart) THEN
    CALL restart_data(step)    ! restart from data in file save.data
    IF (rank .EQ. 0) PRINT *, 'Load from restart dump OK'
  ELSE
    ! auto_load particles
    CALL auto_load
    time = 0.0_num
  ENDIF

  CALL manual_load
  CALL initialise_window ! window.f90
  CALL set_dt
  CALL deallocate_ic

  npart_global = 0
  DO ispecies = 1, n_species
    npart_global = npart_global + species_list(ispecies)%count
  ENDDO

  ! .TRUE. to over_ride balance fraction check
  IF (npart_global .GT. 0) CALL balance_workload(.TRUE.)

  CALL particle_bcs
  CALL efield_bcs
  IF (ic_from_restart) THEN
    CALL bfield_bcs(.TRUE.)
    CALL update_eb_fields_final
    IF (dt_from_restart .GT. 0) THEN
      time = time + dt_from_restart / 2.0_num
    ELSE
      time = time + dt / 2.0_num
    ENDIF
    CALL moving_window
  ELSE
    CALL bfield_final_bcs
  ENDIF

  ! Setup particle migration between species
  IF (use_particle_migration) CALL initialise_migration

  IF (rank .EQ. 0) PRINT *, 'Equilibrium set up OK, running code'
#ifdef PHOTONS
  IF (use_qed) CALL setup_qed_module()
#endif

  walltime_start = MPI_WTIME()
  IF (.NOT.ic_from_restart) CALL output_routines(step) ! diagnostics.f90
  IF (use_ionisation) CALL initialise_ionisation

  DO
    IF ((step .GE. nsteps .AND. nsteps .GE. 0) &
        .OR. (time .GE. t_end) .OR. halt) EXIT
    push = (time .GE. particle_push_start_time)
#ifdef PHOTONS
    IF (push .AND. use_qed .AND. time .GT. qed_start_time) THEN
      CALL qed_update_optical_depth()
    ENDIF
#endif
    CALL update_eb_fields_half
    IF (push) THEN
      ! .FALSE. this time to use load balancing threshold
      IF (use_balance) CALL balance_workload(.FALSE.)
      CALL push_particles
      IF (use_particle_lists) THEN
        ! After this line, the particles can be accessed on a cell by cell basis
        ! Using the particle_species%secondary_list property
        CALL reorder_particles_to_grid

        ! call collision operator
        IF (use_collisions) CALL particle_collisions

        ! Early beta version of particle splitting operator
        IF (use_split) CALL split_particles

        CALL reattach_particles_to_mainlist
      ENDIF
      IF (use_particle_migration) CALL migrate_particles(step)
      IF (use_ionisation) CALL ionise_particles
    ENDIF

    CALL check_for_stop_condition(halt, force_dump)
    IF (halt) EXIT
    step = step + 1
    time = time + dt / 2.0_num
    CALL output_routines(step)
    time = time - dt / 2.0_num

    CALL update_eb_fields_final
    time = time + dt

    CALL moving_window

    ! This section ensures that the particle count for the species_list
    ! objects is accurate. This makes some things easier, but increases
    ! communication
#ifdef PARTICLE_COUNT_UPDATE
    DO ispecies = 1, n_species
      CALL MPI_ALLREDUCE(species_list(ispecies)%attached_list%count, &
          species_list(ispecies)%count, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      species_list(ispecies)%count_update_step = step
    ENDDO
#endif
  ENDDO

  IF (rank .EQ. 0) runtime = MPI_WTIME() - walltime_start

#ifdef PHOTONS
  IF (use_qed) CALL shutdown_qed_module()
#endif

  CALL output_routines(step, force_dump)

  IF (rank .EQ. 0) THEN
    CALL create_full_timestring(runtime, timestring)
    WRITE(*,*) 'Final runtime of core = ' // TRIM(timestring)
  ENDIF

  CALL close_files
  CALL MPI_FINALIZE(errcode)

CONTAINS

  SUBROUTINE cleanup_stop_files()

    INTEGER :: ierr

    IF (rank .NE. 0) RETURN

    OPEN(unit=lu, status='OLD', &
        file=TRIM(data_dir) // '/' // TRIM(stop_file), iostat=ierr)
    IF (ierr .EQ. 0) CLOSE(lu, status='DELETE')

    OPEN(unit=lu, status='OLD', &
        file=TRIM(data_dir) // '/' // TRIM(stop_file_nodump), iostat=ierr)
    IF (ierr .EQ. 0) CLOSE(lu, status='DELETE')

  END SUBROUTINE cleanup_stop_files



  SUBROUTINE check_for_stop_condition(halt, force_dump)

    LOGICAL, INTENT(OUT) :: halt, force_dump
    INTEGER :: ierr
    INTEGER, SAVE :: check_counter = 0
    LOGICAL, SAVE :: adjust_frequency = .TRUE.
    LOGICAL, SAVE :: check_walltime_started = .FALSE.
    LOGICAL :: buffer(2), got_stop_condition, got_stop_file
    REAL(num) :: walltime

    IF (check_stop_frequency .LE. 0 &
        .AND. check_walltime_frequency .LE. 0) RETURN

    ! Once we've started checking the walltime, we may as well check for
    ! stop files as well since it is the broadcast which slows things down.
    IF (adjust_frequency) THEN
      IF (check_walltime_frequency .LE. 0) THEN
        adjust_frequency = .FALSE.
      ELSE IF (time .GE. check_walltime_start) THEN
        adjust_frequency = .FALSE.
        check_walltime_started = .TRUE.
        IF (check_walltime_frequency .LT. check_stop_frequency) &
            check_stop_frequency = check_walltime_frequency
      ENDIF
    ENDIF

    got_stop_condition = .FALSE.
    force_dump = .FALSE.
    check_counter = check_counter + 1

    IF (check_counter .LT. check_stop_frequency) RETURN
    check_counter = 0

    IF (rank .EQ. 0) THEN
      ! First check walltime, if specified
      IF (check_walltime_started) THEN
        walltime = MPI_WTIME()
        IF (walltime - real_walltime_start .GE. stop_at_walltime) THEN
          got_stop_condition = .TRUE.
          force_dump = .TRUE.
          PRINT*,'Stopping because "stop_at_walltime" has been exceeded.'
        ENDIF
      ENDIF

      ! Next check if stop file exists
      OPEN(unit=lu, status='OLD', iostat=ierr, &
          file=TRIM(data_dir) // '/' // TRIM(stop_file))
      IF (ierr .EQ. 0) THEN
        got_stop_file = .TRUE.
        got_stop_condition = .TRUE.
        force_dump = .TRUE.
        CLOSE(lu, status='DELETE')
      ELSE
        OPEN(unit=lu, status='OLD', iostat=ierr, &
            file=TRIM(data_dir) // '/' // TRIM(stop_file_nodump))
        IF (ierr .EQ. 0) THEN
          got_stop_file = .TRUE.
          got_stop_condition = .TRUE.
          force_dump = .FALSE.
          CLOSE(lu, status='DELETE')
        ELSE
          got_stop_file = .FALSE.
        ENDIF
      ENDIF

      IF (got_stop_file) PRINT*,'Stopping because "STOP" file has been found.'

      buffer(1) = got_stop_condition
      buffer(2) = force_dump
    ENDIF

    CALL MPI_BCAST(buffer, 2, MPI_LOGICAL, 0, comm, errcode)
    got_stop_condition = buffer(1)
    force_dump = buffer(2)

    IF (got_stop_condition) halt = .TRUE.

  END SUBROUTINE check_for_stop_condition

END PROGRAM pic
