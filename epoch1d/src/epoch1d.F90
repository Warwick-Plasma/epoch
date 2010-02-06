PROGRAM pic

  ! EPOCH1D is a Birdsall and Langdon type PIC code derived from the PSC
  ! written by Hartmut Ruhl.

  ! The particle pusher (particles.F90) and the field solver (fields.f90) are
  ! almost exact copies of the equivalent routines from PSC, modified slightly
  ! to allow interaction with the changed portions of the code and for
  ! readability. The MPI routines are exactly equivalent to those in PSC, but
  ! are completely rewritten in a form which is easier to extend with arbitrary
  ! fields and particle properties. The support code is entirely new and is not
  ! equivalent to PSC.

  ! EPOCH1D written by C.S.Brady, Centre for Fusion, Space and Astrophysics,
  ! University of Warwick, UK
  ! PSC written by Hartmut Ruhl

  USE setup
  USE ic_module
#ifdef NO_DECK
  USE control
#endif
  USE deck
  USE welcome
  USE diagnostics
  USE field
  USE particles
  USE mpi_routines
  !USE boundary
  USE balance
  USE solve_gauss
#ifdef SPLIT_PARTICLES_AFTER_PUSH
  USE split_particle
#endif
#ifdef PART_IONISE
  USE ionise
#endif
  USE window

  IMPLICIT NONE

  INTEGER :: ispecies, i = 0
  REAL(num) :: walltime_current
  LOGICAL :: halt = .FALSE.

  CALL minimal_init ! setup.f90
  CALL setup_partlists ! partlist.f90
  CALL mpi_minimal_init ! mpi_routines.f90
  CALL welcome_message ! welcome.f90
  CALL register_objects ! custom.f90

  deck_state = c_ds_deck
  IF (rank .EQ. 0) THEN
    PRINT *, "Specify output directory"
    READ(*, *) data_dir
  ENDIF

  CALL MPI_BCAST(data_dir, 64, MPI_CHARACTER, 0, MPI_COMM_WORLD, errcode)
#ifdef NO_DECK
  CALL setup_control_block
  CALL setup_boundaries_block
  CALL setup_species_block
  CALL setup_output_block
  IF (rank .EQ. 0) &
      PRINT *, "Control variables setup OK. Setting initial conditions"
#else
  CALL read_deck("input.deck", .TRUE.)
#endif
  CALL setup_particle_boundaries ! boundary.f90
  CALL mpi_initialise ! mpi_routines.f90
  CALL after_control ! setup.f90
  CALL open_files    ! setup.f90

  ! If the user has specified extended IO options then read the file
  IF (use_extended_io) THEN
    deck_state = c_ds_eio
    CALL read_deck(TRIM(extended_io_file), .TRUE.)
  ENDIF

  ! restart flag is set
  IF (IAND(ictype, c_ic_restart) .NE. 0) THEN
    CALL restart_data    ! restart from data in file SAVE.data
    IF (rank .EQ. 0) PRINT *, "Load from restart dump OK"
    output_file = restart_snapshot
  ELSE
    ! Using autoloader
    IF (IOR(ictype, c_ic_autoload) .NE. 0) THEN
      CALL allocate_ic
    ENDIF
    ! Early internal initialisation
    IF (IAND(ictype, c_ic_early_internal) .NE. 0) THEN
      CALL ic_early  ! define initial profiles
    ENDIF
    ! External initialisation
    IF (IAND(ictype, c_ic_external) .NE. 0) THEN
      deck_state = c_ds_ic
      CALL read_deck(icfile%value, .TRUE.)
    ENDIF
    ! Late internal initialisation
    IF (IAND(ictype, c_ic_late_internal) .NE. 0) THEN
      CALL ic_late    ! define initial profiles
    ENDIF
    ! auto_load particles
    IF (IOR(ictype, c_ic_autoload) .NE. 0) THEN
      CALL auto_load
      CALL deallocate_ic
    ENDIF
    time = 0.0_num
    output_file = 0
  ENDIF

  CALL distribute_particles
  IF (.NOT. neutral_background) CALL do_gauss
  CALL balance_workload(.TRUE.)

  IF (IAND(ictype, c_ic_manual) .NE. 0) THEN
    CALL manual_load
    CALL distribute_particles
    ! .TRUE. to over_ride balance fraction check
    CALL balance_workload(.TRUE.)
  ENDIF

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

  IF (.NOT. restart) THEN
    CALL set_dt
    CALL update_eb_fields_half
    CALL push_particles
    CALL update_eb_fields_final
    IF (rank .EQ. 0) PRINT *, "Equilibrium set up OK, running code"
    CALL output_routines(i) ! diagnostics.f90
  ENDIF
  walltime_current = MPI_WTIME(errcode)

  DO
    IF ((i .GE. nsteps .AND. nsteps .GE. 0) .OR. &
        (time .GE. t_end) .OR. halt) EXIT
    i = i + 1
    CALL set_dt
    CALL update_eb_fields_half
    CALL push_particles
#ifdef SPLIT_PARTICLES_AFTER_PUSH
    ! After this line, the particles can be accessed on a cell by cell basis
    ! Using the particle_family%secondary_list property
    CALL reorder_particles_to_grid
    ! CALL Collisions  !An example, no collision operator yet
    CALL split_particles ! Early beta version of particle splitting operator
    CALL reattach_particles_to_mainlist
#endif
    CALL update_eb_fields_final
#ifdef PART_IONISE
    CALL ionise_particles
#endif
    CALL output_routines(i)
    time = time+dt
    IF (dlb) THEN
      ! .FALSE. this time to use load balancing threshold
      CALL balance_workload(.FALSE.)
    ENDIF
    IF (move_window .AND. .NOT. window_started .AND. &
        time .GE. window_start_time) THEN
      xbc_left = xbc_left_after_move
      xbc_right = xbc_right_after_move
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
          window_shift = &
              window_shift + REAL(FLOOR(window_shift_fraction), num)*dx
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
  ENDDO

  IF (rank .EQ. 0) &
      PRINT *, "Final runtime of core = ", MPI_WTIME(errcode)-walltime_current

  CALL output_routines(i)

  CALL close_files
  CALL MPI_FINALIZE(errcode)

END PROGRAM pic
