MODULE partlist

  USE shared_data
  IMPLICIT NONE

  SAVE

  INTEGER :: nvar

CONTAINS

  SUBROUTINE Setup_PartLists

    nvar=4

#ifdef PER_PARTICLE_WEIGHT
    nvar=nvar+1
#endif

#ifdef PER_PARTICLE_CHARGEMASS
    nvar=nvar+2
#endif

#ifdef PART_DEBUG
    nvar=nvar+2
#endif

  END SUBROUTINE Setup_PartLists

  SUBROUTINE Create_Empty_PartList(PartList)

    TYPE(ParticleList),INTENT(INOUT) :: PartList

    NULLIFY(PartList%Head)
    NULLIFY(PartList%Tail)
    PartList%Count=0
    PartList%Safe=.TRUE.

  END SUBROUTINE Create_Empty_PartList

  SUBROUTINE Create_Unsafe_PartList(PartList,aParticle,n_elements)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER :: aParticle
    INTEGER(KIND=8), INTENT(IN) :: n_elements
    TYPE(Particle),POINTER :: Current
    INTEGER(KIND=8) :: ipart

    CALL Create_Empty_PartList(PartList)

    PartList%Safe=.FALSE.
    Current=>aParticle
    ipart=1
    DO WHILE (ASSOCIATED(Current) .AND. ipart < n_elements)
       ipart=ipart+1
       Current=>Current%Next
    ENDDO
    PartList%Head=>aParticle
    PartList%Tail=>Current
    PartList%Count=ipart

  END SUBROUTINE Create_Unsafe_PartList

  SUBROUTINE Create_Unsafe_PartList_By_Tail(PartList,Head,Tail)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER :: Head,Tail
    TYPE(Particle), POINTER :: Current
    INTEGER(KIND=8) :: ipart

    CALL Create_Empty_PartList(PartList)

    PartList%Safe=.FALSE.
    PartList%Head=>Head
    PartList%Tail=>Tail

    Current=>Head
    ipart=0
    DO WHILE (ASSOCIATED(Current))
       ipart=ipart+1
       Current=>Current%Next
       IF (ASSOCIATED(Current)) THEN
          IF (ASSOCIATED(Current%Prev,TARGET=Tail)) EXIT
       ENDIF
    ENDDO

    PartList%Count=ipart

  END SUBROUTINE Create_Unsafe_PartList_By_Tail

  SUBROUTINE Create_Allocated_PartList(PartList, n_elements)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    INTEGER(8), INTENT(IN) :: n_elements
    TYPE(Particle),POINTER :: NewParticle

    CALL Create_Empty_PartList(PartList)

    DO ipart=0,n_elements-1
       ALLOCATE(NewParticle)
       CALL Add_Particle_To_PartList(PartList,NewParticle)
       NULLIFY(NewParticle)
    ENDDO

  END SUBROUTINE Create_Allocated_PartList

  SUBROUTINE Create_Filled_PartList(PartList,DataIn,n_elements)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    REAL(num),DIMENSION(:), INTENT(IN) :: DataIn
    INTEGER(8), INTENT(IN) :: n_elements
    TYPE(Particle), POINTER :: NewParticle

    INTEGER(KIND=8) :: cpos=0

    CALL Create_Empty_PartList(PartList)


    DO ipart=0,n_elements-1
       ALLOCATE(NewParticle)
       NULLIFY(NewParticle%Prev,NewParticle%Next)
       cpos=ipart*nvar+1
       CALL UnPack_Particle(DataIn(cpos:cpos+nvar),NewParticle)
#ifdef PART_DEBUG
       NewParticle%Processor = Rank
#endif
       CALL Add_Particle_To_PartList(PartList,NewParticle)
       NULLIFY(NewParticle)
    ENDDO


  END SUBROUTINE Create_Filled_PartList

  FUNCTION Test_PartList(PartList)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER :: Current
    INTEGER :: Test_PartList
    INTEGER(KIND=8) :: TestCt

    Test_PartList=0
    TestCt=0

    !Empty list is OK
    IF (.NOT. ASSOCIATED(PartList%Head) .AND. .NOT. ASSOCIATED(PartList%Tail)) THEN
       Test_PartList=0
       RETURN
    ENDIF

    !List with head or tail but not both is broken
    IF (.NOT. ASSOCIATED(PartList%Head) .OR. .NOT. ASSOCIATED(PartList%Tail)) THEN
       Test_PartList=-1
       RETURN
    ENDIF

    !Having head and tail elements which are not the end of a list are OK for
    !unsafe partlists
    IF (ASSOCIATED(PartList%Head%Prev) .AND. PartList%Safe) Test_PartList=IOR(Test_PartList,1)
    IF (ASSOCIATED(PartList%Tail%Next) .AND. PartList%Safe) Test_PartList=IOR(Test_PartList,2)

    !Since we don't KNOW that count is OK (that's what we're checking)
    !Have to check both for end of list and for having reached the tail item
    Current=>PartList%Head
    DO WHILE (ASSOCIATED(Current))
       testct=testct+1
       Current=>Current%Next
       IF (ASSOCIATED(Current)) THEN
          !This tests if we've just jumped to the tail element
          !Allows testing of unsafe partlists
          IF (ASSOCIATED(Current%Prev,TARGET=PartList%Tail)) EXIT
       ENDIF
    ENDDO

    IF (testct .NE. PartList%Count) Test_PartList=IOR(Test_PartList,4)

  END FUNCTION Test_PartList

  SUBROUTINE Destroy_PartList(PartList)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER :: NewParticle,Next
    INTEGER(KIND=8) :: ipart

    !Go through list and delete all the particles in the list
    NewParticle=>PartList%Head
    ipart=0
    DO WHILE (ipart < PartList%Count)
       Next=>NewParticle%Next
       DEALLOCATE(NewParticle)
       NewParticle=>Next
       ipart=ipart+1
    ENDDO

    CALL Create_Empty_PartList(PartList)

  END SUBROUTINE Destroy_PartList

  SUBROUTINE Copy_PartList(PartList1,PartList2)

    TYPE(ParticleList), INTENT(INOUT) :: PartList1, PartList2

    PartList2%Head=>PartList1%Head
    PartList2%Tail=>PartList2%Tail

  END SUBROUTINE Copy_PartList

  SUBROUTINE Append_PartList(Head,Tail)

    TYPE(ParticleList), INTENT(INOUT) :: Head, Tail

    IF (.NOT. Head%Safe .OR. .NOT. Tail%Safe) THEN
       IF (rank .EQ. 0) PRINT *,"Unable to append partlists because one is not safe"
       RETURN
    ENDIF

    IF (ASSOCIATED(Head%Tail)) THEN 
       Head%Tail%Next=>Tail%Head
    ELSE
       Head%Head=>Tail%Head
    ENDIF
    IF (ASSOCIATED(Tail%Head)) Tail%Head%Prev=>Head%Tail
    IF (ASSOCIATED(Tail%Tail)) Head%Tail=>Tail%Tail
    Head%Count=Head%Count+Tail%Count

    CALL Create_Empty_PartList(Tail)

  END SUBROUTINE Append_PartList

  SUBROUTINE Add_Particle_To_PartList(PartList,NewParticle)

    TYPE(ParticleList), INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER :: NewParticle

    !Note that this will work even if you are using an unsafe particle list
    !BE CAREFUL if doing so, it can cause unexpected behaviour

    !if (!Particle) return;
    IF (.NOT. ASSOCIATED(NewParticle)) RETURN
    NULLIFY(NewParticle%Next,NewParticle%Prev)

    !Add particle count
    PartList%Count=PartList%Count+1
    IF(.NOT. ASSOCIATED(PartList%Tail)) THEN
       !Partlist is empty
       PartList%Head=>NewParticle
       PartList%Tail=>NewParticle
       RETURN
    ENDIF

    PartList%Tail%Next=>NewParticle
    NewParticle%Prev=>PartList%Tail
    NULLIFY(NewParticle%Next)
    PartList%Tail=>NewParticle

  END SUBROUTINE Add_Particle_To_PartList

  SUBROUTINE Remove_Particle_From_PartList(PartList,aParticle)

    TYPE(ParticleList), INTENT(INOUT) :: PartList
    TYPE(Particle), POINTER :: aParticle

    !Note that this will work even if you are using an unsafe particle list
    !BE CAREFUL if doing so, it can cause unexpected behaviour

    !Check whether particle is head or tail of list and unlink
    IF (ASSOCIATED(PartList%Head, TARGET=aParticle)) PartList%Head=>aParticle%Next
    IF (ASSOCIATED(PartList%Tail, TARGET=aParticle)) PartList%Tail=>aParticle%Prev

    !Link particles on either side together
    IF (ASSOCIATED(aParticle%Next)) aParticle%Next%Prev=>aParticle%Prev
    IF (ASSOCIATED(aParticle%Prev)) aParticle%Prev%Next=>aParticle%Next

    NULLIFY(aParticle%Next,aParticle%Prev)

    !Decrement counter
    PartList%Count=PartList%Count-1

  END SUBROUTINE Remove_Particle_From_PartList

  SUBROUTINE Pack_Particle(Data,aParticle)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: Data
    TYPE(particle), POINTER :: aParticle
    INTEGER(kind=8) :: cpos, ct, npart_this_it, npart_left
    TYPE(particle), POINTER :: Current

    cpos=1
    Data(cpos)=aParticle%Part_Pos
    cpos=cpos+1
    Data(cpos:cpos+2)=aParticle%part_p
    cpos=cpos+3
#ifdef PER_PARTICLE_WEIGHT
    Data(cpos)=aParticle%weight
    cpos=cpos+1
#endif
#ifdef PER_PARTICLE_CHARGEMASS
    Data(cpos)=aParticle%charge
    Data(cpos+1)=aParticle%mass
    cpos=cpos+2
#endif
#ifdef PART_DEBUG
    Data(cpos)=REAL(aParticle%Processor,num)
    Data(cpos+1)=REAL(aParticle%Processor_at_t0,num)
    cpos=cpos+2
#endif

!!$    PRINT *,"In Pack",rank,aParticle%Part_Pos

  END SUBROUTINE Pack_Particle

  SUBROUTINE Unpack_Particle(Data,aParticle)

    REAL(num), DIMENSION(:), INTENT(IN) :: Data
    TYPE(particle), POINTER :: aParticle
    INTEGER(kind=8) :: cpos, ct, npart_this_it, npart_left
    TYPE(particle), POINTER :: Current

    cpos=1
    aParticle%Part_Pos=Data(cpos)
    cpos=cpos+1
    aParticle%Part_P=Data(cpos:cpos+2)
    cpos=cpos+3
#ifdef PER_PARTICLE_WEIGHT
    aParticle%Weight=Data(cpos)
    cpos=cpos+1
#endif
#ifdef PER_PARTICLE_CHARGEMASS
    aParticle%Charge=Data(cpos)
    aParticle%Mass=Data(cpos+1)
    cpos=cpos+2
#endif
#ifdef PART_DEBUG
    aParticle%Processor=rank
    aParticle%Processor_at_t0=NINT(Data(cpos+1))
    cpos=cpos+2
#endif


  END SUBROUTINE Unpack_Particle

  SUBROUTINE Display_Particle(aParticle)

    TYPE(particle),POINTER :: aParticle

    PRINT *,"Position",aParticle%Part_Pos
    PRINT *,"Momentum",aParticle%Part_P

  END SUBROUTINE Display_Particle

  FUNCTION Compare_Particles(Part1,Part2)

    TYPE(particle),POINTER :: Part1,Part2
    LOGICAL :: Compare_Particles

    Compare_Particles=.TRUE.
    IF ((ABS(Part1%Part_Pos-Part2%Part_Pos)) .NE. 0.0_num  ) Compare_Particles=.FALSE.
    IF (MAXVAL(ABS(Part1%Part_P - Part2%Part_P)) .NE. 0.0_num    ) Compare_Particles=.FALSE.

#ifdef PER_PARTICLE_WEIGHT
    IF (Part1%Weight        .NE. Part2%Weight       ) Compare_Particles=.FALSE.
#endif

#ifdef PER_PARTICLE_CHARGEMASS
    IF (Part1%Charge        .NE. Part2%Charge       ) Compare_Particles=.FALSE.
    IF (Part1%Mass          .NE. Part2%Mass         ) Compare_Particles=.FALSE.
#endif

    IF (.NOT. Compare_Particles) THEN
       CALL Display_Particle(Part1)
       CALL Display_Particle(Part2)
    ENDIF
  END FUNCTION Compare_Particles

  FUNCTION Test_Packed_Particles(PartList,Data,npart_in_data)

    TYPE(ParticleList),INTENT(IN) :: PartList
    REAL(num),DIMENSION(:), INTENT(IN) :: Data
    INTEGER(KIND=8), INTENT(IN) :: npart_in_data
    TYPE(particle),POINTER :: Current
    TYPE(Particle),POINTER :: aParticle
    LOGICAL :: Test_Packed_Particles

    Test_Packed_Particles=.FALSE.


    IF (npart_in_data * nvar .NE. SIZE(Data)) THEN
       PRINT *,"Size of data array does not match specified on",rank,npart_in_data,SIZE(Data)
       RETURN
    ENDIF
    IF (PartList%Count .NE. npart_in_data) THEN
       PRINT *,"Size of data array does not match partlist on",rank
       RETURN
    ENDIF

    ALLOCATE(aParticle)

    Current=>PartList%Head
    DO ipart=0,npart_in_data-1
       CALL Unpack_Particle(Data(ipart*nvar+1:(ipart+1)*nvar),aParticle)
       IF (.NOT. Compare_Particles(aParticle,Current)) THEN
          PRINT *,"BAD PARTICLE ",ipart,"on",rank
          RETURN
       ENDIF
       Current=>Current%Next
    ENDDO

    DEALLOCATE(aParticle)

    Test_Packed_Particles=.TRUE.

  END FUNCTION Test_Packed_Particles


  SUBROUTINE PartList_Send(PartList,Dest)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    INTEGER, INTENT(IN) :: Dest
    REAL(num), DIMENSION(:), ALLOCATABLE :: Data
    INTEGER(kind=8) :: cpos=0, npart_this_it, npart_left, ipart
    TYPE(particle), POINTER :: Current

    npart_left=PartList%Count
    npart_this_it=MIN(npart_left, npart_per_it)
    CALL MPI_SEND(PartList%Count, 1, MPI_INTEGER, Dest, tag, comm, errcode)


    !This is a reduced memory footprint algorithm
    !Try to fix later
!!$    !Copy the data for the particles into a buffer
!!$    ALLOCATE(Data(1:npart_this_it*nvar))
!!$    Current=>PartList%Head
!!$    ct=0
!!$
!!$    npart_left=PartList%Count
!!$    DO WHILE (npart_left .GT. 0)
!!$
!!$       !Copy the data into the buffer
!!$       ct=0
!!$       DO ipart=1,npart_this_it
!!$          cpos=ct*nvar+1
!!$          CALL Pack_Particle(Data(cpos:cpos+nvar),Current)
!!$          Current=>Current%Next
!!$          ct=ct+1
!!$       ENDDO
!!$
!!$       !Send the data to the receiver
!!$       CALL MPI_SEND(Data, npart_this_it*nvar, mpireal, Dest, tag, comm, errcode)
!!$       npart_left=npart_left-npart_this_it
!!$       npart_this_it=MIN(npart_left, npart_per_it)
!!$    ENDDO

    ALLOCATE(Data(1:PartList%Count*nvar))
    Data=0.0_num
    Current=>PartList%Head
    ipart=0
    DO WHILE (ipart < PartList%Count)
       cpos=ipart*nvar+1
       CALL Pack_Particle(Data(cpos:cpos+nvar),Current)
       ipart=ipart+1
       Current=>Current%Next
    ENDDO
    CALL MPI_SEND(Data, npart_left*nvar, mpireal, Dest, tag, comm, errcode)

    DEALLOCATE(Data)

  END SUBROUTINE PartList_Send

  SUBROUTINE PartList_Recv(PartList,Src)

    TYPE(ParticleList),INTENT(INOUT) :: PartList
    INTEGER, INTENT(IN) :: Src
    REAL(num), DIMENSION(:), ALLOCATABLE :: Data
    INTEGER(kind=8) :: cpos=0, ipart, npart_this_it, npart_left,count
    TYPE(ParticleList) :: PartListTemp

    CALL Create_Empty_PartList(PartList)

    Count=0
    CALL MPI_RECV(Count, 1, MPI_INTEGER, Src, tag, comm, status, errcode)

    npart_left=Count
    npart_this_it=MIN(npart_left, npart_per_it)

    !This is a reduced memory footprint algorithm
    !Try to fix later
!!$    !Copy the data for the particles into a buffer
!!$    ALLOCATE(Data(1:npart_this_it*nvar))
!!$    DO WHILE (npart_left .GT. 0)
!!$       !Receive the actual data
!!$       CALL MPI_RECV(Data, npart_this_it*nvar, mpireal, Src, tag, comm, status, errcode)
!!$
!!$       !Copy to temporary partlist and then attach that partlist to the end of the main partlist
!!$       CALL Create_Filled_PartList(PartListTemp,Data, npart_this_it*nvar)
!!$       CALL Append_PartList(PartList,PartListTemp)
!!$
!!$       !Reduce count for next iteration
!!$       npart_left=npart_left-npart_this_it
!!$       npart_this_it=MIN(npart_left, npart_per_it)
!!$    ENDDO

    ALLOCATE(Data(1:Count*nvar))
    Data=0.0_num
    CALL MPI_RECV(Data, Count*nvar, mpireal, Src, tag, comm, status, errcode)
    CALL Create_Filled_PartList(PartList, Data, Count)

    DEALLOCATE(Data)

  END SUBROUTINE PartList_Recv

  SUBROUTINE PartList_SendRecv(PartList_Send,PartList_Recv,Dest,Src)

    TYPE(ParticleList),INTENT(INOUT) :: PartList_Send, PartList_Recv
    INTEGER, INTENT(IN) :: Dest, Src
    REAL(num), DIMENSION(:), ALLOCATABLE :: Data_Send, Data_Recv,Data_Temp
    INTEGER(kind=8) :: cpos=0, ipart =0, npart_this_it
    INTEGER(kind=8) :: npart_send, npart_recv
    TYPE(particle), POINTER :: Current
    LOGICAL :: Test

    !This subroutine doesn't try to use memory efficient buffering, it sends all the particles at once
    !This should work for boundary calls, but don't try it for any other reason

    npart_send=PartList_Send%Count
    npart_recv=0
    CALL MPI_SENDRECV(npart_send, 1, MPI_INTEGER, Dest, tag, npart_recv, 1, MPI_INTEGER,&
         Src, tag, comm, status, errcode)

    !Copy the data for the particles into a buffer
    ALLOCATE(Data_Send(1:npart_send*nvar))
    ALLOCATE(Data_Recv(1:npart_recv*nvar))
    ALLOCATE(Data_Temp(1:nvar))

    !Pack particles to send into buffer
    Current=>PartList_Send%Head
    ipart=0
    DO WHILE (ipart < PartList_Send%Count)
       cpos=ipart*nvar+1
       CALL Pack_Particle(Data_Temp,Current)
       Data_Send(cpos:cpos+nvar-1)=Data_Temp
       ipart=ipart+1
       Current=>Current%Next
    ENDDO

    !No longer need the sending partlist, so destroy it to save some memory
    CALL Destroy_PartList(PartList_Send)

    !Actual MPI commands
    CALL MPI_SENDRECV(Data_Send, npart_send*nvar, mpireal, Dest, tag, Data_Recv, npart_recv*nvar, mpireal,&
         Src, tag, comm, status, errcode)

    DEALLOCATE(Data_Send)
    CALL Create_Filled_PartList(PartList_Recv, Data_recv, npart_recv)
    DEALLOCATE(Data_Recv)

  END SUBROUTINE PartList_SendRecv

END MODULE partlist
