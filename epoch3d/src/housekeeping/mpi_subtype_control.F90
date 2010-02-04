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
      get_total_local_particles = get_total_local_particles + &
          particle_species(ispecies)%attached_list%count
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
        get_total_local_dumped_particles = get_total_local_dumped_particles + &
            particle_species(ispecies)%attached_list%count
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
    INTEGER(KIND=8) :: npart_local

    npart_local = get_total_local_dumped_particles(force_restart)

    subtype_field = create_current_field_subtype()
    subtype_particle_var = create_particle_subtype(npart_local)

  END SUBROUTINE create_subtypes



  !----------------------------------------------------------------------------
  ! create_current_field_subtype - Creates the subtype corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subtype()

    INTEGER :: create_current_field_subtype

    create_current_field_subtype = &
        create_field_subtype(nx, ny, nz, cell_x_start(coordinates(3)+1), &
            cell_y_start(coordinates(2)+1), cell_z_start(coordinates(1)+1))

  END FUNCTION create_current_field_subtype



  !---------------------------------------------------------------------------
  ! create_subtypes_for_load - Creates subtypes when the code loads initial
  ! conditions from a file
  !---------------------------------------------------------------------------

  SUBROUTINE create_subtypes_for_load(npart_local)

    ! This subroutines creates the MPI types which represent the data for the
    ! field and particles data. It is used when reading data. To this end, it
    ! takes npart_local rather than determining it from the data structures

    INTEGER(KIND=8), INTENT(IN) :: npart_local

    subtype_field = create_current_field_subtype()
    subtype_particle_var = create_particle_subtype(npart_local)

  END SUBROUTINE create_subtypes_for_load



  !----------------------------------------------------------------------------
  ! create_particle_subtype - Creates a subtype representing the local particles
  !----------------------------------------------------------------------------

  FUNCTION create_particle_subtype(npart_local)

    INTEGER(KIND=8), INTENT(IN) :: npart_local
    INTEGER :: create_particle_subtype

    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths, starts

    ! Create the subarray for the particles in this problem: subtype decribes
    ! where this process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local, 1, MPI_INTEGER8, npart_each_rank, 1, &
        MPI_INTEGER8, comm, errcode)

    ALLOCATE(lengths(1), starts(1))
    lengths = npart_local
    starts = 0
    DO ix = 1, rank
      starts = starts+npart_each_rank(ix)
    ENDDO

    CALL MPI_TYPE_INDEXED(1, lengths, starts, mpireal, &
        create_particle_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_particle_subtype, errcode)

    DEALLOCATE(lengths, starts)

  END FUNCTION create_particle_subtype



  !----------------------------------------------------------------------------
  ! create_field_subtype - Creates a subtype representing the local processor
  ! for any arbitrary arrangement of an array covering the entire spatial
  ! domain. Only used directly during load balancing
  !----------------------------------------------------------------------------

  FUNCTION create_field_subtype(nx_local, ny_local, nz_local, &
      cell_start_x_local, cell_start_y_local, cell_start_z_local)

    INTEGER, INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER, INTENT(IN) :: cell_start_x_local
    INTEGER, INTENT(IN) :: cell_start_y_local
    INTEGER, INTENT(IN) :: cell_start_z_local
    INTEGER :: create_field_subtype

    create_field_subtype = create_3d_array_subtype((/nx_local, &
        ny_local, nz_local/), (/nx_global, ny_global, nz_global/), &
        (/cell_start_x_local, cell_start_y_local, cell_start_z_local/))

  END FUNCTION create_field_subtype



  !----------------------------------------------------------------------------
  ! create_3d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 3D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_3d_array_subtype(n_local, n_global, start)

    INTEGER, DIMENSION(3), INTENT(IN) :: n_local
    INTEGER, DIMENSION(3), INTENT(IN) :: n_global
    INTEGER, DIMENSION(3), INTENT(IN) :: start
    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths, starts
    INTEGER :: ipoint
    INTEGER :: create_3d_array_subtype

    ALLOCATE(lengths(1:n_local(2) * n_local(3)))
    ALLOCATE(starts(1:n_local(2) * n_local(3)))
    lengths = n_local(1)
    ipoint = 0
    DO iz = 0, n_local(3)-1
      DO iy = 0, n_local(2)-1
        ipoint = ipoint+1
        starts(ipoint) = (start(3)+iz-1) * n_global(1) * n_global(2) + &
            (start(2)+iy-1) * n_global(1) + start(1) -1
      ENDDO
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2)*n_local(3), lengths, starts, mpireal, &
        create_3d_array_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_3d_array_subtype, errcode)
    DEALLOCATE(lengths, starts)

  END FUNCTION create_3d_array_subtype



  !----------------------------------------------------------------------------
  ! create_2d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 2D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_2d_array_subtype(n_local, n_global, start)

    INTEGER, DIMENSION(2), INTENT(IN) :: n_local
    INTEGER, DIMENSION(2), INTENT(IN) :: n_global
    INTEGER, DIMENSION(2), INTENT(IN) :: start
    INTEGER, DIMENSION(:), ALLOCATABLE :: lengths, starts
    INTEGER :: ipoint
    INTEGER :: create_2d_array_subtype

    ALLOCATE(lengths(1:n_local(2)), starts(1:n_local(2)))
    lengths = n_local(1)
    ipoint = 0
    DO iy = 0, n_local(2)-1
      ipoint = ipoint+1
      starts(ipoint) = (start(2)+iy-1) * n_global(1) + start(1) -1
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2), lengths, starts, mpireal, &
        create_2d_array_subtype, errcode)
    CALL MPI_TYPE_COMMIT(create_2d_array_subtype, errcode)
    DEALLOCATE(lengths, starts)

  END FUNCTION create_2d_array_subtype

END MODULE mpi_subtype_control
