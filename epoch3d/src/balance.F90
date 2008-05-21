MODULE balance

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE CreateSubtypes

    INTEGER, DIMENSION(5) :: length,disp,type
    INTEGER(8), DIMENSION(:),ALLOCATABLE :: npart_each_rank

    ALLOCATE(npart_each_rank(1:nproc))

    CALL MPI_ALLGATHER(npart,1,MPI_INTEGER8,npart_each_rank,1,MPI_INTEGER8,comm,errcode)

    ! Create the subarray for the particle float properties in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    length(1)=1
    length(2)=npart
    length(3)=1
    disp(1)=0
    disp(2)=0
    DO ix=1,rank
       disp(2)=disp(2)+npart_each_rank(ix)*num
    ENDDO
    disp(3)= npart_global * num
    type(1)=MPI_LB
    type(2)=mpireal
    type(3)=MPI_UB
    IF (subtype_particle_var /= 0) CALL MPI_TYPE_FREE(subtype_particle_var,errcode)
    CALL MPI_TYPE_STRUCT(3,length,disp,type,subtype_particle_var,errcode)
    CALL MPI_TYPE_COMMIT(subtype_particle_var,errcode)

    ! Create the subarray for the particle integer properties in this problem: subtype decribes where this
    ! process's data fits into the global picture.
    length(1)=1
    length(2)=npart
    length(3)=1
    disp(1)=0
    disp(2)=0
    DO ix=1,rank
       disp(2)=disp(2)+npart_each_rank(ix)*4
    ENDDO
    disp(3)= npart_global * 4
    type(1)=MPI_LB
    type(2)=MPI_INTEGER
    type(3)=MPI_UB
    IF (subtype_particle_int /= 0) CALL MPI_TYPE_FREE(subtype_particle_int,errcode)
    CALL MPI_TYPE_STRUCT(3,length,disp,type,subtype_particle_int,errcode)
    CALL MPI_TYPE_COMMIT(subtype_particle_int,errcode)

    DEALLOCATE(npart_each_rank)

  END SUBROUTINE CreateSubtypes


END MODULE balance
