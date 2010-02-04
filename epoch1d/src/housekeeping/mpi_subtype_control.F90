MODULE mpi_subtype_control

  ! ---------------------------------------------------------------------------------
  ! This module contains the subroutines which create the subtypes used in
  ! IO
  ! ---------------------------------------------------------------------------------

  USE shared_data
  IMPLICIT NONE

CONTAINS

  ! ---------------------------------------------------------------------------------
  ! get_total_local_particles - Returns the number of particles on this
  ! processor.
  ! ---------------------------------------------------------------------------------
  FUNCTION get_total_local_particles()

    ! This subroutine describes the total number of particles on the current processor
    ! It simply sums over every particle species

    INTEGER(KIND=8) :: get_total_local_particles
    INTEGER :: ispecies

    get_total_local_particles = 0
    DO ispecies = 1, n_species
      get_total_local_particles = get_total_local_particles+particle_species(ispecies)%attached_list%count
    ENDDO

  END FUNCTION get_total_local_particles



  ! ---------------------------------------------------------------------------------
  ! get_total_local_dumped_particles - Returns the number of particles on this
  ! processor which should be written to disk. Parameter is whether this number should
  ! be calculated for a normal dump or a restart dump (all species are written at restart)
  ! ---------------------------------------------------------------------------------
  FUNCTION get_total_local_dumped_particles(force_restart)

    ! This subroutine describes the total number of particles on the current processor
    ! which are members of species with the dump = T attribute in the input deck
    ! If force_restart = .TRUE. then the subroutine simply counts all the particles

    LOGICAL, INTENT(IN) :: force_restart
    INTEGER(KIND=8) :: get_total_local_dumped_particles
    INTEGER :: ispecies

    get_total_local_dumped_particles = 0
    DO ispecies = 1, n_species
      IF (particle_species(ispecies)%dump .OR. force_restart) THEN
        get_total_local_dumped_particles = get_total_local_dumped_particles+particle_species(ispecies)%attached_list%count
      ENDIF
    ENDDO

  END FUNCTION get_total_local_dumped_particles



  ! ---------------------------------------------------------------------------------
  ! CreateSubtypes - Creates the subtypes used by the main output routines
  ! Run just before output takes place
  ! ---------------------------------------------------------------------------------
  SUBROUTINE create_subtypes(force_restart)

    ! This subroutines creates the MPI types which represent the data for the field and
    ! particles data. It is used when writing data
    LOGICAL, INTENT(IN) :: force_restart
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: npart_local
    INTEGER :: n_dump_species, ispecies, index

    ! count the number of dumped particles of each species
    n_dump_species = 0
    DO ispecies = 1, n_species
      IF (particle_species(ispecies)%dump) n_dump_species = n_dump_species+1
    ENDDO

    ALLOCATE(npart_local(1:n_dump_species))
    index = 1
    DO ispecies = 1, n_species
      IF (particle_species(ispecies)%dump) THEN
        npart_local(index) = particle_species(ispecies)%attached_list%count
        index = index+1
      ENDIF
    ENDDO

    ! Actually create the subtypes
    subtype_field = create_current_field_subtype()
    subtype_particle_var = create_ordered_particle_subtype(n_dump_species, npart_local)

  END SUBROUTINE create_subtypes



  ! ---------------------------------------------------------------------------------
  ! create_current_field_subtype - Creates the subtype corresponding to the current
  ! load balanced geometry
  ! ---------------------------------------------------------------------------------
  FUNCTION create_current_field_subtype()

    INTEGER :: create_current_field_subtype

    create_current_field_subtype = create_field_subtype(nx, cell_x_start(coordinates(1)+1))

  END FUNCTION create_current_field_subtype



  ! ---------------------------------------------------------------------------------
  ! create_subtypes_for_load - Creates subtypes when the code loads initial conditions
  ! From a file
  ! ---------------------------------------------------------------------------------
  SUBROUTINE create_subtypes_for_load(npart_local)

    ! This subroutines creates the MPI types which represent the data for the field and
    ! particles data. It is used when reading data. To this end, it takes npart_local
    ! rather than determining it from the data structures

    INTEGER(KIND=8), INTENT(IN) :: npart_local

    subtype_field = create_current_field_subtype()
    subtype_particle_var = create_particle_subtype(npart_local)

  END SUBROUTINE create_subtypes_for_load



  ! ---------------------------------------------------------------------------------
  ! create_particle_subtype - Creates a subtype representing the local particles
  ! ---------------------------------------------------------------------------------
  FUNCTION create_particle_subtype(npart_local)

    INTEGER(KIND=8), INTENT(IN) :: npart_local
    INTEGER :: create_particle_subtype

    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths, starts

    ! Create the subarray for the particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local, 1, MPI_INTEGER8, npart_each_rank, 1, MPI_INTEGER8, comm, errcode)

    ALLOCATE(lengths(1), starts(1))
    lengths = npart_local
    starts = 0
    DO ix = 1, rank
      starts = starts+npart_each_rank(ix)
    ENDDO

    CALL MPI_TYPE_INDEXED(1, lengths, starts, mpireal, create_particle_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_particle_subtype, errcode)

    DEALLOCATE(lengths, starts)

  END FUNCTION create_particle_subtype



  ! ---------------------------------------------------------------------------------
  ! create_ordered_particle_subtype - Creates a subtype representing the local particles
  ! ---------------------------------------------------------------------------------
  FUNCTION create_ordered_particle_subtype(n_dump_species, npart_local)

    INTEGER :: create_ordered_particle_subtype
    INTEGER, INTENT(IN) :: n_dump_species
    INTEGER(KIND=8), DIMENSION(n_dump_species), INTENT(IN) :: npart_local
    INTEGER :: ispecies
    INTEGER(KIND=8), DIMENSION(:, :), ALLOCATABLE :: npart_each_rank
    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths, starts

    ALLOCATE(npart_each_rank(1:n_dump_species, 1:nproc))

    ! Create the subarray for the particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local, n_dump_species, MPI_INTEGER8, npart_each_rank, n_dump_species, MPI_INTEGER8, comm, errcode)

    ALLOCATE(lengths(n_dump_species), starts(n_dump_species))
    lengths = npart_local
    DO ispecies = 1, n_dump_species
      starts(ispecies) = 0
      DO ix = 1, ispecies-1
        starts(ispecies) = starts(ispecies)+SUM(npart_each_rank(ix, :), 1)
      ENDDO
      DO ix = 1, rank
        starts(ispecies) = starts(ispecies)+npart_each_rank(ispecies, ix)
      ENDDO
    ENDDO

    CALL MPI_TYPE_INDEXED(n_dump_species, lengths, starts, mpireal, create_ordered_particle_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_ordered_particle_subtype, errcode)

    DEALLOCATE(lengths, starts)

  END FUNCTION create_ordered_particle_subtype



  ! ---------------------------------------------------------------------------------
  ! create_field_subtype - Creates a subtype representing the local processor
  ! For any arbitrary arrangement of the array. Only used directly during load balancing
  ! ---------------------------------------------------------------------------------
  FUNCTION create_field_subtype(nx_local, cell_start_x_local)

    INTEGER, INTENT(IN) :: nx_local
    INTEGER, INTENT(IN) :: cell_start_x_local
    INTEGER, DIMENSION(3) :: length, disp, types
    INTEGER :: create_field_subtype

    ! lengths = nx_local
    ! starts = cell_start_x_local-1

    ! CALL MPI_TYPE_INDEXED(1, lengths, starts, mpireal, create_field_subtype, errcode)
    ! CALL MPI_TYPE_COMMIT(create_field_subtype, errcode)

    length(1) = 1
    length(2) = nx
    length(3) = 1
    disp(1) = 0
    disp(2) = (cell_start_x_local-1)*num
    disp(3) = nx_global * num
    types(1) = MPI_LB
    types(2) = mpireal
    types(3) = MPI_UB
    CALL MPI_TYPE_STRUCT(3, length, disp, types, create_field_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_field_subtype, errcode)

  END FUNCTION create_field_subtype

END MODULE mpi_subtype_control
