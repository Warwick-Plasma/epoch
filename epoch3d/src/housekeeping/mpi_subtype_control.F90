MODULE mpi_subtype_control

  !---------------------------------------------------------------------------------
  ! This module contains the subroutines which create the subtypes used in
  ! IO
  !---------------------------------------------------------------------------------

  USE shared_data
  IMPLICIT NONE

CONTAINS
  !---------------------------------------------------------------------------------
  ! Get_Total_Local_Particles - Returns the number of particles on this
  ! processor.
  !---------------------------------------------------------------------------------
  FUNCTION Get_Total_Local_Particles()

    !This subroutine describes the total number of particles on the current processor
    !It simply sums over every particle species

    INTEGER(KIND=8) :: Get_Total_Local_Particles
    INTEGER :: iSpecies
    Get_Total_Local_Particles=0
    DO iSpecies = 1,nspecies
       Get_Total_Local_Particles=Get_Total_Local_Particles+ParticleSpecies(iSpecies)%AttachedList%Count
    ENDDO

  END FUNCTION Get_Total_Local_Particles

  !---------------------------------------------------------------------------------
  ! Get_Total_Local_Dumped_Particles - Returns the number of particles on this
  ! processor which should be written to disk. Parameter is whether this number should
  ! be calculated for a normal dump or a restart dump (all species are written at restart)
  !---------------------------------------------------------------------------------
  FUNCTION Get_Total_Local_Dumped_Particles(Force_Restart)

    !This subroutine describes the total number of particles on the current processor
    !which are members of species with the Dump=T attribute in the input deck
    !If FORCE_RESTART=.TRUE. then the subroutine simply counts all the particles

    LOGICAL,INTENT(IN) :: Force_Restart
    INTEGER(KIND=8) :: Get_Total_Local_Dumped_Particles
    INTEGER :: iSpecies

    Get_Total_Local_Dumped_Particles=0
    DO iSpecies = 1,nspecies
       IF (ParticleSpecies(iSpecies)%Dump .OR. Force_Restart) THEN
          Get_Total_Local_Dumped_Particles=Get_Total_Local_Dumped_Particles+ParticleSpecies(iSpecies)%AttachedList%Count
       ENDIF
    ENDDO

  END FUNCTION Get_Total_Local_Dumped_Particles

  !---------------------------------------------------------------------------------
  ! CreateSubtypes - Creates the subtypes used by the main output routines
  ! Run just before output takes place
  !---------------------------------------------------------------------------------
  SUBROUTINE Create_Subtypes(Force_Restart)

    !This subroutines creates the MPI types which represent the data for the field and
    !particles data. It is used when writing data
    LOGICAL,INTENT(IN) :: Force_Restart
    INTEGER(KIND=8) :: npart_local

    npart_local = Get_Total_Local_Dumped_Particles(Force_Restart)

    subtype_field=Create_Current_Field_Subtype()
    subtype_particle_var=Create_Particle_Subtype(npart_local)

  END SUBROUTINE Create_Subtypes
  !---------------------------------------------------------------------------------
  ! Create_Current_Field_Subtype - Creates the subtype corresponding to the current
  ! load balanced geometry
  !---------------------------------------------------------------------------------
  FUNCTION Create_Current_Field_Subtype
    INTEGER :: Create_Current_Field_Subtype
    Create_Current_Field_Subtype=Create_Field_Subtype(nx,ny,nz,&
         cell_x_start(coordinates(3)+1),cell_y_start(coordinates(2)+1),cell_z_start(coordinates(1)+1))
  END FUNCTION Create_Current_Field_Subtype

  !---------------------------------------------------------------------------------
  ! Create_Subtypes_For_Load - Creates subtypes when the code loads initial conditions
  ! From a file
  !---------------------------------------------------------------------------------
  SUBROUTINE Create_Subtypes_For_Load(npart_local)

    !This subroutines creates the MPI types which represent the data for the field and
    !particles data. It is used when reading data. To this end, it takes npart_local
    !rather than determining it from the data structures

    INTEGER(KIND=8),INTENT(IN) :: npart_local

    subtype_field=Create_Current_Field_Subtype()
    subtype_particle_var=Create_Particle_Subtype(npart_local)

  END SUBROUTINE Create_Subtypes_For_Load

  !---------------------------------------------------------------------------------
  ! Create_Particle_Subtype - Creates a subtype representing the local particles
  !---------------------------------------------------------------------------------
  FUNCTION Create_Particle_Subtype(nPart_Local)

    INTEGER(KIND=8),INTENT(IN) :: nPart_Local
    INTEGER :: Create_Particle_Subtype

    INTEGER,DIMENSION(:),ALLOCATABLE :: lengths, starts


    ! Create the subarray for the particles in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local,1,MPI_INTEGER8,npart_each_rank,1,MPI_INTEGER8,comm,errcode)

    ALLOCATE(lengths(1),starts(1))
    lengths=npart_local
    starts=0
    DO ix=1,rank
       starts=starts+npart_each_rank(ix)
    ENDDO

    CALL MPI_TYPE_INDEXED(1,lengths,starts,mpireal,Create_Particle_Subtype,errcode)
    CALL MPI_TYPE_COMMIT(Create_Particle_Subtype,errcode)

    DEALLOCATE(lengths,starts)

  END FUNCTION Create_Particle_Subtype

  !---------------------------------------------------------------------------------
  ! Create_Field_Subtype - Creates a subtype representing the local processor
  ! For any arbitrary arrangement of an array covering the entire spatial domain. 
  ! Only used directly during load balancing
  !---------------------------------------------------------------------------------
  FUNCTION Create_Field_Subtype(nx_local,ny_local,nz_local,cell_start_x_local,cell_start_y_local,cell_start_z_local)
    INTEGER,INTENT(IN) :: nx_local, ny_local, nz_local
    INTEGER,INTENT(IN) :: cell_start_x_local, cell_start_y_local, cell_start_z_local
    INTEGER :: Create_Field_Subtype

    Create_Field_Subtype = Create_3D_Array_Subtype((/nx_local,ny_local,nz_local/),(/nx_global,ny_global,nz_global/),&
         (/cell_start_x_local,cell_start_y_local,cell_start_z_local/))
  END FUNCTION Create_Field_Subtype

  !---------------------------------------------------------------------------------
  ! Create_3D_Array_Subtype - Creates a subtype representing the local fraction of a
  ! completely arbitrary 3D array. Does not assume anything about the domain at all.
  !---------------------------------------------------------------------------------
  FUNCTION Create_3D_Array_Subtype(n_local,n_global,start)
    INTEGER,DIMENSION(3),INTENT(IN) :: n_local
    INTEGER,DIMENSION(3),INTENT(IN) :: n_global
    INTEGER,DIMENSION(3),INTENT(IN) :: start
    INTEGER, DIMENSION(:),ALLOCATABLE :: lengths,starts
    INTEGER :: iPoint
    INTEGER :: Create_3D_Array_Subtype
    ALLOCATE(lengths(1:n_local(2) * n_local(3)),starts(1:n_local(2) * n_local(3)))
    lengths=n_local(1)
    ipoint=0
    DO iz=0,n_local(3)-1
       DO iy=0,n_local(2)-1
          ipoint=ipoint+1
          starts(ipoint)=(start(3)+iz-1) * n_global(1) * n_global(2) + (start(2)+iy-1) * n_global(1) + start(1) -1
       ENDDO
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2)*n_local(3),lengths,starts,mpireal,Create_3D_Array_Subtype,errcode)
    CALL MPI_TYPE_COMMIT(Create_3D_Array_Subtype,errcode)
    DEALLOCATE(lengths,starts)
  END FUNCTION Create_3D_Array_Subtype
  !---------------------------------------------------------------------------------
  ! Create_2D_Array_Subtype - Creates a subtype representing the local fraction of a
  ! completely arbitrary 2D array. Does not assume anything about the domain at all.
  !---------------------------------------------------------------------------------
  FUNCTION Create_2D_Array_Subtype(n_local,n_global,start)
    INTEGER,DIMENSION(2),INTENT(IN) :: n_local
    INTEGER,DIMENSION(2),INTENT(IN) :: n_global
    INTEGER,DIMENSION(2),INTENT(IN) :: start
    INTEGER, DIMENSION(:),ALLOCATABLE :: lengths,starts
    INTEGER :: iPoint
    INTEGER :: Create_2D_Array_Subtype
    ALLOCATE(lengths(1:n_local(2)),starts(1:n_local(2)))
    lengths=n_local(1)
    ipoint=0
    DO iy=0,n_local(2)-1
       ipoint=ipoint+1
       starts(ipoint)=(start(2)+iy-1) * n_global(1) + start(1) -1
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2),lengths,starts,mpireal,Create_2D_Array_Subtype,errcode)
    CALL MPI_TYPE_COMMIT(Create_2D_Array_Subtype,errcode)
    DEALLOCATE(lengths,starts)
  END FUNCTION Create_2D_Array_Subtype

END MODULE mpi_subtype_control
