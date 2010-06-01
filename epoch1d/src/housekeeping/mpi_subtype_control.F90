MODULE mpi_subtype_control

  !----------------------------------------------------------------------------
  ! This module contains the subroutines which create the subtypes used in
  ! IO
  !----------------------------------------------------------------------------

  USE shared_data

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! get_total_local_particles - Returns the number of particles on this
  ! processor.
  !----------------------------------------------------------------------------

  FUNCTION get_total_local_particles()

    ! This subroutine describes the total number of particles on the current
    ! processor. It simply sums over every particle species

    INTEGER(KIND=8) :: get_total_local_particles
    INTEGER :: ispecies

    get_total_local_particles = 0
    DO ispecies = 1, n_species
      get_total_local_particles = get_total_local_particles &
          + particle_species(ispecies)%attached_list%count
    ENDDO

  END FUNCTION get_total_local_particles



  !----------------------------------------------------------------------------
  ! get_total_local_dumped_particles - Returns the number of particles on this
  ! processor which should be written to disk. Parameter is whether this number
  ! should be calculated for a normal dump or a restart dump (all species are
  ! written at restart)
  !----------------------------------------------------------------------------

  FUNCTION get_total_local_dumped_particles(force_restart)

    ! This subroutine describes the total number of particles on the current
    ! processor which are members of species with the dump = T attribute in
    ! the input deck. If force_restart = .TRUE. then the subroutine simply
    ! counts all the particles

    LOGICAL, INTENT(IN) :: force_restart
    INTEGER(KIND=8) :: get_total_local_dumped_particles
    INTEGER :: ispecies

    get_total_local_dumped_particles = 0
    DO ispecies = 1, n_species
      IF (particle_species(ispecies)%dump .OR. force_restart) THEN
        get_total_local_dumped_particles = get_total_local_dumped_particles &
            + particle_species(ispecies)%attached_list%count
      ENDIF
    ENDDO

  END FUNCTION get_total_local_dumped_particles



  !----------------------------------------------------------------------------
  ! CreateSubtypes - Creates the subtypes used by the main output routines
  ! Run just before output takes place
  !----------------------------------------------------------------------------

  SUBROUTINE create_subtypes(force_restart)

    ! This subroutines creates the MPI types which represent the data for the
    ! field and particles data. It is used when writing data
    LOGICAL, INTENT(IN) :: force_restart
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: npart_local
    INTEGER :: n_dump_species, ispecies, index

    ! count the number of dumped particles of each species
    n_dump_species = 0
    DO ispecies = 1, n_species
      IF (particle_species(ispecies)%dump .OR. force_restart) THEN
        n_dump_species = n_dump_species + 1
      ENDIF
    ENDDO

    ALLOCATE(npart_local(n_dump_species))
    index = 1
    DO ispecies = 1, n_species
      IF (particle_species(ispecies)%dump .OR. force_restart) THEN
        npart_local(index) = particle_species(ispecies)%attached_list%count
        index = index + 1
      ENDIF
    ENDDO

    ! Actually create the subtypes
    subtype_field = create_current_field_subtype()
    subtype_particle_var = &
        create_ordered_particle_subtype(n_dump_species, npart_local)
    subarray_field = create_current_field_subarray()

  END SUBROUTINE create_subtypes



  !----------------------------------------------------------------------------
  ! create_current_field_subtype - Creates the subtype corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subtype()

    INTEGER :: create_current_field_subtype

    create_current_field_subtype = create_1d_array_subtype( &
        (/nx/), (/nx_global/), (/cell_x_min(coordinates(1)+1)/))

  END FUNCTION create_current_field_subtype



  !----------------------------------------------------------------------------
  ! create_current_field_subarray - Creates the subarray corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subarray()

    INTEGER :: create_current_field_subarray

    create_current_field_subarray = create_field_subarray(nx)

  END FUNCTION create_current_field_subarray



  !----------------------------------------------------------------------------
  ! create_subtypes_for_load - Creates subtypes when the code loads initial
  ! conditions from a file
  !----------------------------------------------------------------------------

  SUBROUTINE create_subtypes_for_load(npart_local)

    ! This subroutines creates the MPI types which represent the data for the
    ! field and particles data. It is used when reading data. To this end, it
    ! takes npart_local rather than determining it from the data structures

    INTEGER(KIND=8), INTENT(IN) :: npart_local

    subtype_field = create_current_field_subtype()
    subtype_particle_var = create_particle_subtype(npart_local)
    subarray_field = create_current_field_subarray()

  END SUBROUTINE create_subtypes_for_load



  !----------------------------------------------------------------------------
  ! create_particle_subtype - Creates a subtype representing the local
  ! particles
  !----------------------------------------------------------------------------

  FUNCTION create_particle_subtype(npart_local)

    INTEGER :: create_particle_subtype
    INTEGER(KIND=8), INTENT(IN) :: npart_local

    create_particle_subtype = &
        create_ordered_particle_subtype(1, (/npart_local/) )

  END FUNCTION create_particle_subtype



  !----------------------------------------------------------------------------
  ! create_ordered_particle_subtype - Creates a subtype representing the local
  ! particles
  !----------------------------------------------------------------------------

  FUNCTION create_ordered_particle_subtype(n_dump_species, npart_local)

    INTEGER :: create_ordered_particle_subtype
    INTEGER, INTENT(IN) :: n_dump_species
    INTEGER(KIND=8), DIMENSION(n_dump_species), INTENT(IN) :: npart_local
    INTEGER :: ispecies
    INTEGER(KIND=8), DIMENSION(:,:), ALLOCATABLE :: npart_each_rank
    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(:), ALLOCATABLE :: disp
    INTEGER(KIND=MPI_ADDRESS_KIND) :: particles_to_skip, sz

    ALLOCATE(npart_each_rank(n_dump_species, nproc))

    ! Create the subarray for the particles in this problem: subtype decribes
    ! where this process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local, n_dump_species, MPI_INTEGER8, &
        npart_each_rank, n_dump_species, MPI_INTEGER8, comm, errcode)

    ALLOCATE(lengths(n_dump_species), disp(n_dump_species))

    ! If npart_local is bigger than an integer then the data will not
    ! get written properly. This would require about 48GB per processor
    ! so it is unlikely to be a problem any time soon.
    sz = realsize
    lengths = INT(npart_local)
    particles_to_skip = 0

    DO ispecies = 1, n_dump_species
      DO ix = 1, rank
        particles_to_skip = particles_to_skip + npart_each_rank(ispecies,ix)
      ENDDO
      disp(ispecies) = particles_to_skip * sz
      DO ix = rank+1, nproc
        particles_to_skip = particles_to_skip + npart_each_rank(ispecies,ix)
      ENDDO
    ENDDO

    CALL MPI_TYPE_CREATE_HINDEXED(n_dump_species, lengths, disp, mpireal, &
        create_ordered_particle_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_ordered_particle_subtype, errcode)

    DEALLOCATE(lengths, disp)

  END FUNCTION create_ordered_particle_subtype



  !----------------------------------------------------------------------------
  ! create_field_subtype - Creates a subtype representing the local processor
  ! for any arbitrary arrangement of an array covering the entire spatial
  ! domain. Only used directly during load balancing
  !----------------------------------------------------------------------------

  FUNCTION create_field_subtype(nx_local, cell_start_x_local)

    INTEGER, INTENT(IN) :: nx_local
    INTEGER, INTENT(IN) :: cell_start_x_local
    INTEGER :: create_field_subtype

    create_field_subtype = create_1d_array_subtype( &
        (/nx_local/), (/nx_global/), (/cell_start_x_local/))

  END FUNCTION create_field_subtype



  !----------------------------------------------------------------------------
  ! create_1d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 1D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_1d_array_subtype(n_local, n_global, start)

    INTEGER, DIMENSION(1), INTENT(IN) :: n_local
    INTEGER, DIMENSION(1), INTENT(IN) :: n_global
    INTEGER, DIMENSION(1), INTENT(IN) :: start
    INTEGER, DIMENSION(1) :: lengths
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(1) :: disp
    INTEGER(KIND=MPI_ADDRESS_KIND) :: sz
    INTEGER :: create_1d_array_subtype

    sz = realsize
    lengths = n_local(1)
    disp(1) = (start(1) - 1) * sz

    CALL MPI_TYPE_CREATE_HINDEXED(1, lengths, disp, mpireal, &
        create_1d_array_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_1d_array_subtype, errcode)

  END FUNCTION create_1d_array_subtype



  !----------------------------------------------------------------------------
  ! create_2d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 2D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_2d_array_subtype(n_local, n_global, start)

    INTEGER, DIMENSION(2), INTENT(IN) :: n_local
    INTEGER, DIMENSION(2), INTENT(IN) :: n_global
    INTEGER, DIMENSION(2), INTENT(IN) :: start
    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(:), ALLOCATABLE :: disp
    INTEGER(KIND=MPI_ADDRESS_KIND) :: iy, sz
    INTEGER :: create_2d_array_subtype
    INTEGER :: ipoint

    ALLOCATE(lengths(n_local(2)), disp(n_local(2)))

    sz = realsize
    lengths = n_local(1)
    ipoint = 0

    DO iy = -1, n_local(2)-2
      ipoint = ipoint + 1
      disp(ipoint) = ((start(2) + iy) * n_global(1) + start(1) - 1) * sz
    ENDDO

    CALL MPI_TYPE_CREATE_HINDEXED(ipoint, lengths, disp, mpireal, &
        create_2d_array_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_2d_array_subtype, errcode)

    DEALLOCATE(lengths, disp)

  END FUNCTION create_2d_array_subtype



  !----------------------------------------------------------------------------
  ! create_3d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 3D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_3d_array_subtype(n_local, n_global, start)

    INTEGER, DIMENSION(3), INTENT(IN) :: n_local
    INTEGER, DIMENSION(3), INTENT(IN) :: n_global
    INTEGER, DIMENSION(3), INTENT(IN) :: start
    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(:), ALLOCATABLE :: disp
    INTEGER(KIND=MPI_ADDRESS_KIND) :: iy, iz, sz
    INTEGER :: create_3d_array_subtype
    INTEGER :: ipoint

    ALLOCATE(lengths(1:n_local(2) * n_local(3)))
    ALLOCATE(disp(1:n_local(2) * n_local(3)))

    sz = realsize
    lengths = n_local(1)
    ipoint = 0

    DO iz = -1, n_local(3)-2
      DO iy = -1, n_local(2)-2
        ipoint = ipoint + 1
        disp(ipoint) = (((start(3) + iz) * n_global(2) &
            + start(2) + iy) * n_global(1) + start(1) - 1) * sz
      ENDDO
    ENDDO

    CALL MPI_TYPE_CREATE_HINDEXED(ipoint, lengths, disp, mpireal, &
        create_3d_array_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_3d_array_subtype, errcode)

    DEALLOCATE(lengths, disp)

  END FUNCTION create_3d_array_subtype



  FUNCTION create_field_subarray(n1, n2, n3)

    INTEGER, INTENT(IN) :: n1
    INTEGER, INTENT(IN), OPTIONAL :: n2, n3
    INTEGER, DIMENSION(3) :: subsizes, sizes, starts
    INTEGER :: i, ndim, create_field_subarray

    subsizes(1) = n1
    ndim = 1
    IF (PRESENT(n2)) THEN
      subsizes(2) = n2
      ndim = 2
    ENDIF
    IF (PRESENT(n3)) THEN
      subsizes(3) = n3
      ndim = 3
    ENDIF

    DO i = 1, ndim
      starts(i) = 3
      sizes(i) = subsizes(i) + 6
    ENDDO

    CALL MPI_TYPE_CREATE_SUBARRAY(ndim, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, create_field_subarray, errcode)

    CALL MPI_TYPE_COMMIT(create_field_subarray, errcode)

  END FUNCTION create_field_subarray

END MODULE mpi_subtype_control
