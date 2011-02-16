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
#ifdef SPLIT_PARTICLES_AFTER_PUSH
  USE split_particle
#endif
#ifdef PARTICLE_IONISE
  USE ionise
#endif

  IMPLICIT NONE

  INTEGER :: ispecies, i = 0
  LOGICAL :: halt = .FALSE.

  CALL minimal_init ! setup.f90
  CALL setup_partlists ! partlist.f90
  CALL mpi_minimal_init ! mpi_routines.f90
  CALL get_job_id(jobid)
  CALL welcome_message ! welcome.f90
  CALL register_objects ! custom.f90

  IF (rank .EQ. 0) THEN
    PRINT *, "Specify output directory"
    READ(*, *) data_dir
  ENDIF

  CALL MPI_BCAST(data_dir, 64, MPI_CHARACTER, 0, MPI_COMM_WORLD, errcode)
  deck_state = c_ds_deck
  CALL read_deck("input.deck", .TRUE.)
  CALL setup_particle_boundaries ! boundary.f90
  CALL mpi_initialise ! mpi_routines.f90
  CALL after_control ! setup.f90
  CALL open_files    ! setup.f90

  IF (move_window) CALL allocate_window ! window.f90

  ! Read extended IO options
  deck_state = c_ds_eio
  CALL read_deck("input.deck", .TRUE.)

  ! restart flag is set
  IF (ic_from_restart) THEN
    CALL restart_data(i)    ! restart from data in file save.data
    IF (rank .EQ. 0) PRINT *, "Load from restart dump OK"
    output_file = restart_snapshot + 1
  ELSE
    ! Using autoloader
    CALL allocate_ic
    ! External initialisation
    deck_state = c_ds_ic
    CALL read_deck("input.deck", .TRUE.)
    ! auto_load particles
    CALL auto_load
    CALL deallocate_ic
    time = 0.0_num
    output_file = 0
  ENDIF

  CALL manual_load
  ! .TRUE. to over_ride balance fraction check
  CALL balance_workload(.TRUE.)

  ! npart_global isn't really used anymore, just check where it is used
  IF (npart_global .LT. 0) THEN
    npart_global = 0
    DO ispecies = 1, n_species
      npart_global = npart_global+particle_species(ispecies)%count
    ENDDO
  ENDIF

  CALL particle_bcs
  CALL efield_bcs
  CALL bfield_bcs(.FALSE.)

  IF (.NOT. ic_from_restart) THEN
    CALL set_dt
    CALL update_eb_fields_half
    CALL push_particles
    CALL update_eb_fields_final
    IF (rank .EQ. 0) PRINT *, "Equilibrium set up OK, running code"
    walltime_start = MPI_WTIME()
    CALL output_routines(i) ! diagnostics.f90
  ELSE
    walltime_start = MPI_WTIME()
  ENDIF

  DO
    IF ((i .GE. nsteps .AND. nsteps .GE. 0) &
        .OR. (time .GE. t_end) .OR. halt) EXIT
    i = i + 1
    CALL set_dt
    CALL update_eb_fields_half
    CALL push_particles
#ifdef SPLIT_PARTICLES_AFTER_PUSH
    ! After this line, the particles can be accessed on a cell by cell basis
    ! Using the particle_family%secondary_list property
    CALL reorder_particles_to_grid
    ! CALL Collisions  !An example, no collision operator yet
#ifdef PER_PARTICLE_WEIGHT
    CALL split_particles ! Early beta version of particle splitting operator
#endif
    CALL reattach_particles_to_mainlist
#endif
    CALL update_eb_fields_final
#ifdef PARTICLE_IONISE
    CALL ionise_particles
#endif
    time = time+dt
    IF (dlb) THEN
      ! .FALSE. this time to use load balancing threshold
      CALL balance_workload(.FALSE.)
    ENDIF
    IF (move_window .AND. .NOT. window_started &
        .AND. time .GE. window_start_time) THEN
      bc_field(c_bd_x_min) = bc_x_min_after_move
      bc_field(c_bd_x_max) = bc_x_max_after_move
      CALL setup_particle_boundaries
      CALL setup_communicator
      window_started = .TRUE.
    ENDIF

    ! If we have a moving window then update the window position
    IF (move_window .AND. window_started) THEN
      window_shift_fraction = window_shift_fraction + dt*window_v_x/dx
      ! Allow for posibility of having jumped two cells at once
      IF (FLOOR(window_shift_fraction) .GE. 1.0_num) THEN
        IF (use_offset_grid) THEN
          window_shift(1) = &
              window_shift(1) + REAL(FLOOR(window_shift_fraction), num)*dx
        ENDIF
        CALL shift_window
        window_shift_fraction = &
            window_shift_fraction - REAL(FLOOR(window_shift_fraction), num)
        CALL particle_bcs
      ENDIF
    ENDIF

    ! This section ensures that the particle count for the particle_species
    ! objects is accurate. This makes some things easier, but increases
    ! communication
#ifdef PARTICLE_COUNT_UPDATE
    DO ispecies = 1, n_species
      CALL MPI_ALLREDUCE(particle_species(ispecies)%attached_list%count, &
          particle_species(ispecies)%count, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
    ENDDO
#endif
    IF (halt) EXIT
    CALL output_routines(i)
  ENDDO

  IF (rank .EQ. 0) &
      PRINT *, "Final runtime of core = ", MPI_WTIME() - walltime_start

  IF (halt) CALL output_routines(i)

  CALL close_files
  CALL MPI_FINALIZE(errcode)

END PROGRAM pic
