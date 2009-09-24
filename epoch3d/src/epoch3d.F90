PROGRAM pic

  !EPOCH3D is a Birdsall and Langdon type PIC code derived from the PSC written by Hartmut Ruhl. 

  !The particle pusher (particles.F90) and the field solver (fields.f90) are almost exact copies of the equivalent routines from PSC,
  !modified slightly to allow interaction with the changed portions of the code and for readability. The MPI routines are exactly 
  !equivalent to those in PSC, but are completely rewritten in a form which is easier to extend with arbitrary fields and particle properties.
  !The support code is entirely new and is not equivalent to PSC.

  !EPOCH3D written by C.S.Brady, Centre for Fusion, Space and Astrophysics, University of Warwick, UK
  !PSC written by Hartmut Ruhl

  USE shared_data
  USE setup
  USE initial_conditions
  USE deck
  USE welcome
  USE diagnostics
  USE field
  USE particles
  USE mpi_routines
  USE boundary
  USE balance
  USE helper
  USE solve_gauss
#ifdef SPLIT_PARTICLES_AFTER_PUSH
  USE split_particle
#endif
  USE custom_deck
  USE window

  IMPLICIT NONE

  INTEGER :: i = 0
  REAL(num) :: walltime_current,dwalltime
  LOGICAL :: halt = .FALSE.

  CALL minimal_init !setup.f90
  CALL setup_partlists !partlist.f90
  CALL mpi_minimal_init !mpi_routines.f90
  CALL welcome_message !welcome.f90
  CALL RegisterObjects !custom.f90
  Deck_State = DS_DECK
	!Ask for output directory name on rank 0
  IF (rank .EQ. 0) THEN
     PRINT *,"Specify output directory"
     READ(*,*) data_dir
  ENDIF
  CALL MPI_BCAST(data_dir,64,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
  CALL Read_Deck("input.deck",.TRUE.)
  CALL Setup_Particle_Boundaries !boundary.f90
  CALL mpi_initialise !mpi_routines.f90
  CALL After_Control !setup.f90
  CALL open_files    !setup.f90

  !If the user has specified extended IO options then read the file
  IF (use_extended_io) THEN
     Deck_State = DS_EIO
     CALL Read_Deck(TRIM(extended_io_file),.TRUE.)
  ENDIF

  !Restart flag is set
  IF (IAND(ictype,IC_RESTART) .NE. 0) THEN
     CALL restart_data    !restart from data in file save.data          
     IF (rank .EQ. 0) PRINT *,"Load from restart dump OK"
     output_file=restart_snapshot
  ELSE
     !Using autoloader
     IF (IOR(ictype,IC_AUTOLOAD) .NE. 0) THEN
        CALL AllocateIC
     ENDIF
     !Early internal initialisation
     IF (IAND(ictype,IC_EARLY_INTERNAL) .NE. 0) THEN
        CALL IC_Early  !define initial profiles
     ENDIF
     !External initialisation
     IF (IAND(ictype,IC_EXTERNAL) .NE. 0) THEN
        Deck_State=DS_IC
        CALL Read_Deck(icfile%value,.TRUE.)
     ENDIF
     !Late internal initialisation
     IF (IAND(ictype,IC_LATE_INTERNAL) .NE. 0) THEN
        CALL IC_Late    !define initial profiles
     ENDIF
     !Autoload particles
     IF (IOR(ictype,IC_AUTOLOAD) .NE. 0) THEN
        CALL Autoload
        CALL DeallocateIC
     ENDIF
     time = 0.0_num
     output_file=0
  END IF

  CALL Distribute_Particles
  IF (.NOT. Neutral_Background) CALL Do_Gauss
  CALL Balance_Workload(.TRUE.)

  IF (IAND(ictype,IC_MANUAL) .NE. 0) THEN
     CALL ManualLoad
     CALL Distribute_Particles
     !.TRUE. to override balance fraction check
     CALL Balance_Workload(.TRUE.)
  ENDIF

  CALL Particle_Bcs
  CALL EField_BCS
  CALL BField_BCS(.FALSE.)

  CALL MPI_BARRIER(comm,errcode)

  IF (.NOT. restart) THEN
     CALL set_dt
     CALL update_eb_fields_half
     CALL push_particles
     CALL update_eb_fields_final
     IF (rank .EQ. 0) PRINT *,"Equilibrium set up OK, running code"
     CALL output_routines(i) !diagnostics.f90
  ENDIF
  walltime_current=MPI_WTIME(errcode)
  DO
     IF ((i >= nsteps .AND. nsteps >=0) .OR. (time >= t_end) .OR. halt) EXIT
     i = i + 1
     CALL set_dt
     CALL update_eb_fields_half
     CALL push_particles
#ifdef SPLIT_PARTICLES_AFTER_PUSH
     !After this line, the particles can be accessed on a cell by cell basis
     !Using the ParticleFamily%SecondaryList property
     CALL Reorder_Particles_to_grid
     !CALL Collisions  !An example, no collision operator yet
     CALL Split_Particles !Early beta version of particle splitting operator
     CALL Reattach_Particles_to_mainlist
#endif
     CALL update_eb_fields_final
     CALL output_routines(i)
     time=time+dt
     IF (DLB) THEN
        !.FALSE. this time to use load balancing threshold
        CALL Balance_Workload(.FALSE.)
     ENDIF
     IF (move_window .AND. .NOT. window_started .AND. time .GE. window_start_time) THEN
        xbc_left=xbc_left_after_move
        xbc_right=xbc_right_after_move
        CALL Setup_Particle_Boundaries
        CALL Setup_Communicator
        window_started=.TRUE.
     ENDIF
     !If we have a moving window then update the window position
     IF (move_window .AND. window_started) THEN
        window_shift_fraction=window_shift_fraction + dt*window_v_x/dx
        !Allow for posibility of having jumped two cells at once
        IF (FLOOR(window_shift_fraction) .GE. 1.0_num) THEN
           IF (use_offset_grid) THEN
              Window_shift(1)=Window_Shift(1)+REAL(FLOOR(window_shift_fraction),num)*dx
           ENDIF
           CALL Shift_Window
           window_shift_fraction=window_shift_fraction-REAL(FLOOR(window_shift_fraction),num)
           CALL Particle_BCS
        ENDIF
     ENDIF
     !This section ensures that the particle count for the ParticleSpecies objects is accurate
     !This makes some things easier, but increases communication
#ifdef PARTICLE_COUNT_UPDATE
     DO iSpecies=1,nSpecies
        CALL MPI_ALLREDUCE(ParticleSpecies(iSpecies)%AttachedList%Count,ParticleSpecies(iSpecies)%Count&
             ,1,MPI_INTEGER8,MPI_SUM,comm,errcode)
     ENDDO
#endif
     IF (halt) EXIT
  END DO
  IF (rank .EQ. 0) PRINT *,"Final runtime of core =",MPI_WTIME(errcode)-walltime_current

  CALL output_routines(i)

  CALL close_files
  CALL MPI_FINALIZE(errcode)

END PROGRAM pic









