PROGRAM pic

  ! EPOCH2D is a Birdsall and Langdon type PIC code derived from the PSC
  ! written by Hartmut Ruhl.

  ! The particle pusher (particles.F90) and the field solver (fields.f90) are
  ! almost exact copies of the equivalent routines from PSC, modified slightly
  ! to allow interaction with the changed portions of the code and for
  ! readability. The MPI routines are exactly equivalent to those in PSC, but
  ! are completely rewritten in a form which is easier to extend with arbitrary
  ! fields and particle properties. The support code is entirely new and is not
  ! equivalent to PSC.

  ! EPOCH2D written by C.S.Brady, Centre for Fusion, Space and Astrophysics,
  ! University of Warwick, UK
  ! PSC written by Hartmut Ruhl

  USE mpi
  USE balance
  USE boundary
  USE sdf_job_info
  USE custom_parser
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
  USE ionise
#ifdef PHOTONS
  USE photons
#endif

  IMPLICIT NONE

  INTEGER :: ispecies, step = 0
  LOGICAL :: halt = .FALSE.
  CHARACTER(LEN=64) :: deck_file = 'input.deck'

#ifdef COLLISIONS_TEST
  ! used for testing
  CALL test_collisions
  STOP
#endif

  CALL mpi_minimal_init ! mpi_routines.f90
  CALL minimal_init     ! setup.f90
  CALL setup_partlists  ! partlist.f90
  CALL get_job_id(jobid)
  CALL welcome_message  ! welcome.f90
  CALL register_objects ! custom.f90

  IF (rank .EQ. 0) THEN
    PRINT *, 'Specify output directory'
    READ(*, *) data_dir
  ENDIF

  CALL MPI_BCAST(data_dir, 64, MPI_CHARACTER, 0, MPI_COMM_WORLD, errcode)
  CALL read_deck(deck_file, .TRUE., c_ds_first)
  CALL setup_particle_boundaries ! boundary.f90
  CALL mpi_initialise  ! mpi_routines.f90
  CALL after_control   ! setup.f90
  CALL open_files      ! setup.f90

  ! restart flag is set
  IF (ic_from_restart) THEN
    ! Re-scan the input deck for items which require allocated memory
    CALL read_deck(deck_file, .TRUE., c_ds_last)
    CALL after_deck_last
    CALL restart_data(step)    ! restart from data in file save.data
    IF (rank .EQ. 0) PRINT *, 'Load from restart dump OK'
    output_file = restart_snapshot + 1
  ELSE
    ! Using autoloader
    CALL allocate_ic
    ! Re-scan the input deck for items which require allocated memory
    CALL read_deck(deck_file, .TRUE., c_ds_last)
    CALL after_deck_last
    ! auto_load particles
    CALL auto_load
    time = 0.0_num
    output_file = 0
  ENDIF

  CALL manual_load
  CALL initialise_window ! window.f90
  CALL set_dt
  IF (.NOT. ic_from_restart) CALL deallocate_ic

  npart_global = 0
  DO ispecies = 1, n_species
    npart_global = npart_global + species_list(ispecies)%count
  ENDDO

  ! .TRUE. to over_ride balance fraction check
  IF (npart_global .GT. 0) CALL balance_workload(.TRUE.)

  CALL particle_bcs
  CALL efield_bcs
  CALL bfield_final_bcs

  IF (rank .EQ. 0) PRINT *, 'Equilibrium set up OK, running code'

  walltime_start = MPI_WTIME()
  CALL output_routines(step) ! diagnostics.f90
#ifdef PHOTONS
  IF (use_qed) CALL setup_qed_module()
#endif
  IF (use_ionisation) CALL initialise_ionisation

  DO
    IF ((step .GE. nsteps .AND. nsteps .GE. 0) &
        .OR. (time .GE. t_end) .OR. halt) EXIT
    CALL update_eb_fields_half
#ifdef PHOTONS
    IF (time .GT. qed_start_time .AND. use_qed) THEN
      CALL qed_update_optical_depth()
    ENDIF
#endif
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
    IF (use_ionisation) CALL ionise_particles
    CALL update_eb_fields_final
    step = step + 1
    time = time + dt
    IF (dlb) THEN
      ! .FALSE. this time to use load balancing threshold
      CALL balance_workload(.FALSE.)
    ENDIF

    CALL moving_window

    ! This section ensures that the particle count for the species_list
    ! objects is accurate. This makes some things easier, but increases
    ! communication
#ifdef PARTICLE_COUNT_UPDATE
    DO ispecies = 1, n_species
      CALL MPI_ALLREDUCE(species_list(ispecies)%attached_list%count, &
          species_list(ispecies)%count, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
    ENDDO
#endif
    IF (halt) EXIT
    CALL output_routines(step)
  ENDDO

#ifdef PHOTONS
  IF (use_qed) CALL shutdown_qed_module()
#endif

  IF (rank .EQ. 0) &
      PRINT *, 'Final runtime of core = ', MPI_WTIME() - walltime_start

  IF (halt) CALL output_routines(step)

  CALL close_files
  CALL MPI_FINALIZE(errcode)

END PROGRAM pic
