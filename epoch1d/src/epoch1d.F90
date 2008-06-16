PROGRAM pic

!EPOCH1D is a Birdsall and Langdon type PIC code derived from the PSC written by Hartmut Ruhl. 

!The particle pusher (particles.F90) and the field solver (fields.f90) are almost exact copies of the equivalent routines from PSC,
!modified slightly to allow interaction with the changed portions of the code and for readability. The MPI routines are exactly 
!equivalent to those in PSC, but are completely rewritten in a form which is easier to extend with arbitrary fields and particle properties.
!The support code is entirely new and is not equivalent to PSC.

!EPOCH1D written by C.S.Brady, Centre for Fusion, Space and Astrophysics, University of Warwick, UK
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
  USE balance

  IMPLICIT NONE

  INTEGER :: i = 0
  REAL(num) :: walltime_current,dwalltime
  CHARACTER(LEN=64) :: DeckName
  TYPE(PARTICLE),POINTER :: Current,Next
  REAL(num) :: cell_x_r
  INTEGER :: cell_x1
  INTEGER,DIMENSION(:),ALLOCATABLE :: Counted



  CALL minimal_init !setup.f90
  CALL setup_partlists !partlist.f90
  CALL mpi_minimal_init !mpi_routines.f90
  CALL welcome_message !welcome.f90
#ifdef INPUTFILE
  !Input deck file is specified at runtime
  IF (rank .EQ. 0) THEN
     PRINT *,"Specify output directory"
     READ(*,*) data_dir
  ENDIF
  CALL MPI_BCAST(data_dir,64,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
  CALL Read_Deck("input.deck",.TRUE.)
#else
  !Input deck file is prespecified to be "input.deck"
  CALL Read_Deck("input.deck",.TRUE.) !deck.f90
#endif
  CALL Setup_Particle_Boundaries
  CALL mpi_initialise !mpi_routines.f90
  CALL After_Control !setup.f90
  CALL open_files    !setup.f90

  IF (restart) THEN            
     CALL restart_data    ! restart from data in file save.data          
  ELSE
     CALL equilibrium    ! define initial profiles
     time = 0.0_num
     restart_snapshot=0
  END IF
  IF (domain == DO_DECOMPOSED) THEN 
     CALL Distribute_Particles
     CALL Balance_Workload(.TRUE.)
  ENDIF

  CALL Particle_Bcs()
  CALL EField_Bcs()
  CALL BField_Bcs()

  IF (rank .EQ. 0) PRINT *,"Equilibrium set up OK, running code"
  CALL set_dt
  CALL update_eb_fields_half
  !CALL push_particles
  CALL update_eb_fields_final

  CALL output_routines(i) !diagnostics.f90


  DO
     IF ((i >= nsteps .AND. nsteps >=0) .OR. (time >= t_end)) EXIT
     i = i + 1
     CALL set_dt
     CALL update_eb_fields_half
     CALL push_particles
     CALL update_eb_fields_final
     IF (DLB .AND. domain==DO_DECOMPOSED) CALL Balance_Workload(.FALSE.)
     time=time+dt
     CALL output_routines(i)
     IF (halt) EXIT

  END DO

  CALL mpi_close
  CALL close_files
  CALL MPI_FINALIZE(errcode)

END PROGRAM pic
