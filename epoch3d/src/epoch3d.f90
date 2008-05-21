PROGRAM pic

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

  IMPLICIT NONE

  INTEGER :: i = 0
  REAL(num) :: walltime_current,dwalltime
  LOGICAL :: halt = .FALSE.
  TYPE(particle),POINTER :: Current
  INTEGER :: fail


  CALL minimal_init !setup.f90
  CALL mpi_minimal_init !mpi_routines.f90
  CALL welcome_message !welcome.f90
  CALL Read_Deck !deck.f90
  CALL mpi_initialise !mpi_routines.f90
  CALL After_Control !setup.f90
  CALL open_files    !setup.f90

  IF (restart) THEN            
     CALL restart_data    ! restart from data in file save.data
     IF (rank .EQ. 0) PRINT *,"Restart dump loaded OK, running code"
     output_file=restart_snapshot+1
  ELSE
     !As strange as it may seem, this is all that is done in the VLI
     !Program of PSC
     CALL equilibrium    ! define initial profiles
     CALL particle_bcs
     time = 0.0_num
     restart_snapshot=0
     output_file=0
     CALL set_dt
     CALL update_eb_fields_half
     CALL push_particles
     CALL update_eb_fields_final
     IF (rank .EQ. 0) PRINT *,"Equilibrium set up OK, running code"
  END IF




  CALL output_routines(i) !diagnostics.f90
  DO
     IF ((i >= nsteps .AND. nsteps >=0) .OR. (time >= t_end) .OR. halt) EXIT
     i = i + 1
     CALL set_dt
     CALL update_eb_fields_half
     CALL push_particles
     CALL update_eb_fields_final
     time=time+dt
     CALL output_routines(i)
  END DO

  PRINT *,"Finished with npart=",npart,"on processor",rank
  CALL mpi_close
  CALL close_files
  CALL MPI_FINALIZE(errcode)

END PROGRAM pic









