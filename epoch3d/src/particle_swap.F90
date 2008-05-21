MODULE particle_swap

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE swap_particles(upper_head,lower_head,upper_send,lower_send,rank_upper,rank_lower)

    TYPE(particle),INTENT(INOUT),POINTER :: upper_head,lower_head
    INTEGER(8), INTENT(INOUT) :: upper_send,lower_send
    INTEGER,INTENT(IN) :: rank_upper,rank_lower
    INTEGER(8),PARAMETER :: npart_per_it = 10000
    REAL(num),DIMENSION(:),ALLOCATABLE :: Data_Out,Data_in
    INTEGER :: nvar=7,current_position,count
    TYPE(particle),POINTER :: Current,Next
    INTEGER(8) :: upper_recv,lower_recv
    INTEGER :: Coords(3)


#ifdef PER_PARTICLE_WEIGHT
    nvar=nvar+1
#endif
#ifdef PER_PARTICLE_CHARGEMASS
    nvar=nvar+2
#endif

    upper_recv=0
    lower_recv=0




    !Send and receive the number of particles to go in both directions
    CALL MPI_SENDRECV(upper_send,1,MPI_INTEGER8,rank_upper,tag,lower_recv,1,MPI_INTEGER8,rank_lower,tag,comm,status,errcode)
    CALL MPI_SENDRECV(lower_send,1,MPI_INTEGER8,rank_lower,tag,upper_recv,1,MPI_INTEGER8,rank_upper,tag,comm,status,errcode)

!!$
!!$    IF (upper_send .GT. 0) PRINT *,"upper",rank,upper_send
!!$    IF (lower_send .GT. 0) PRINT *,"lower",rank,lower_send


    !Send to "upper", receive from "lower"
    Current=>upper_head
    ALLOCATE(Data_Out(1:upper_send*nvar), Data_In(1:lower_recv*nvar))
    count=0
    DO ipart=1,upper_send
       Next=>Current%Next

       !Copy data into array for sending
       current_position=count*nvar+1
       Data_Out(current_position:current_position+2)=Current%Part_Pos
       current_position=current_position+3
       Data_Out(current_position:current_position+2)=Current%part_p
       current_position=current_position+3
       Data_Out(current_position)=REAL(Current%part_species,num)+0.1
       current_position=current_position+1
#ifdef PER_PARTICLE_WEIGHT
       Data_Out(current_position)=Current%weight
       current_position=current_position+1
#endif
#ifdef PER_PARTICLE_CHARGEMASS
       Data_Out(current_position)=Current%charge
       Data_Out(current_position+1)=Current%mass
       current_position=current_position+2
#endif
       !Destroy the particle reference
       DEALLOCATE(Current)
       Current=>Next
       count=count+1
    ENDDO
    !Send/recv the data 
!!$    IF (upper_send .NE. 0) THEN
!!$       CALL MPI_CART_COORDS(comm,rank_upper,3,coords,errcode)
!!$       WRITE(10+rank,*) "Sending",upper_send*nvar,"to",coords,rank
!!$    ENDIF
!!$    IF (lower_recv .NE. 0) THEN 
!!$       CALL MPI_CART_COORDS(comm,rank_lower,3,coords,errcode)
!!$       WRITE(10+rank,*) "Receving",lower_recv*nvar,"from",coords,rank
!!$    ENDIF
    CALL MPI_SENDRECV(Data_Out,upper_send*nvar,mpireal,rank_upper,tag,Data_In,&
         lower_recv*nvar,mpireal,rank_lower,tag,comm,status,errcode)
    !Finished with Data_Out, so destroy it
    DEALLOCATE(Data_Out)


    !Now copy the new data into particles
    !Add new particles at end of array
    Current=>tail
    DO ipart=1,lower_recv
       ALLOCATE(Current%Next)
       Current=>Current%Next
       !This keeps the tail record updated
       tail=>Current
       !Copy data from array
       current_position=(ipart-1)*nvar+1
       Current%Part_pos=Data_In(current_position:current_position+2)
       current_position=current_position+3
       Current%Part_P = Data_In(current_position:current_position+2)
       current_position=current_position+3
       Current%Part_species=INT(Data_In(current_position))
       current_position=current_position+1
#ifdef PER_PARTICLE_WEIGHT
       Current%weight=Data_In(current_position)
       current_position=current_position+1
#endif
#ifdef PER_PARTICLE_CHARGEMASS
       Current%Charge=Data_In(current_position)
       Current%Mass=Data_In(current_position+1)
       current_position=current_position+2
#endif
       NULLIFY(Current%Next)
    ENDDO
    DEALLOCATE(Data_In)

    !Now completed one half of the exchange sweep, now do the other half

    !Send to "lower", receive from "upper"
    Current=>lower_head
    ALLOCATE(Data_Out(1:lower_send*nvar), Data_In(1:upper_recv*nvar))
    IF (.NOT. ALLOCATED(Data_Out) .OR. .NOT. ALLOCATED(Data_In)) THEN
       PRINT *,"Unable to malloc memory on node",rank
       CALL MPI_ABORT(comm,errcode)
    ENDIF
    count=0
    DO ipart=1,lower_send
       Next=>Current%Next

       !Copy data into array for sending
       current_position=count*nvar+1
       Data_Out(current_position:current_position+2)=Current%Part_Pos
       current_position=current_position+3
       Data_Out(current_position:current_position+2)=Current%part_p
       current_position=current_position+3
       Data_Out(current_position)=REAL(Current%part_species,num)+0.1
       current_position=current_position+1
#ifdef PER_PARTICLE_WEIGHT
       Data_Out(current_position)=Current%weight
       current_position=current_position+1
#endif
#ifdef PER_PARTICLE_CHARGEMASS
       Data_Out(current_position)=Current%charge
       Data_Out(current_position+1)=Current%mass
       current_position=current_position+2
#endif
       !Destroy the particle reference
       DEALLOCATE(Current)
       Current=>Next
       count=count+1
    ENDDO
    !Send/recv the data 
!!$    IF (lower_send .NE. 0) THEN 
!!$       CALL MPI_CART_COORDS(comm,rank_lower,3,coords,errcode)
!!$       WRITE(10+rank,*) "Sending",lower_send*nvar,"to",coords,lower_send
!!$    ENDIF
!!$    IF (upper_recv .NE. 0) THEN
!!$       CALL MPI_CART_COORDS(comm,rank_upper,3,coords,errcode)
!!$       WRITE(10+rank,*) "Receving",upper_recv*nvar,"from",coords,upper_send
!!$    ENDIF

    CALL MPI_SENDRECV(Data_Out,lower_send*nvar,mpireal,rank_lower,tag,Data_In,&
         upper_recv*nvar,mpireal,rank_upper,tag,comm,status,errcode)
    !Finished with Data_Out, so destroy it
    DEALLOCATE(Data_Out)
    !Now copy the new data into particles
    !Add new particles at end of array
    Current=>Tail
    DO ipart=1,upper_recv
       ALLOCATE(Current%Next)
       Current=>Current%Next
       !This keeps the tail record updated
       tail=>Current
       !Copy data from array
       current_position=(ipart-1)*nvar+1
       Current%Part_pos=Data_In(current_position:current_position+2)

       !IF (Current%part_pos(3) .GT. z_end_local) WRITE(rank+10,*),"Z high",current%part_pos(3),z_end_local
       !       IF (Current%part_pos(3) .LT. z_start_local) PRINT *,"Z low",Current%part_pos(3)

       current_position=current_position+3
       Current%Part_P = Data_In(current_position:current_position+2)
       current_position=current_position+3
       Current%Part_species=INT(Data_In(current_position))
       current_position=current_position+1
#ifdef PER_PARTICLE_WEIGHT
       Current%weight=Data_In(current_position)
       current_position=current_position+1
#endif
#ifdef PER_PARTICLE_CHARGEMASS
       Current%Charge=Data_In(current_position)
       Current%Mass=Data_In(current_position+1)
       current_position=current_position+2
#endif
       NULLIFY(Current%Next)
    ENDDO
    DEALLOCATE(Data_In)

    npart=npart-upper_send-lower_send
    npart=npart+upper_recv+lower_recv

    NULLIFY(upper_head,lower_head)
    upper_send=0
    lower_send=0


  END SUBROUTINE swap_particles

END MODULE particle_swap
