MODULE partlist

  USE mpi
  USE shared_data

  IMPLICIT NONE

  SAVE

  INTEGER :: nvar

  TYPE pointer_item
    TYPE(particle), POINTER :: part
    TYPE(pointer_item), POINTER :: next
  END TYPE pointer_item

  TYPE pointer_list
    TYPE(pointer_item), POINTER :: head, tail
  END TYPE pointer_list

CONTAINS

  SUBROUTINE setup_partlists

    nvar = 3 + c_ndims
#ifdef PER_PARTICLE_WEIGHT
    nvar = nvar+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    nvar = nvar+2
#endif
#ifdef PARTICLE_DEBUG
    nvar = nvar+2
#endif
#if PARTICLE_ID || PARTICLE_ID4
    nvar = nvar+1
#endif

  END SUBROUTINE setup_partlists



  SUBROUTINE create_empty_partlist(partlist)

    TYPE(particle_list), INTENT(INOUT) :: partlist

    NULLIFY(partlist%head)
    NULLIFY(partlist%tail)
    partlist%count = 0
    partlist%id_update = 0
    partlist%safe = .TRUE.

  END SUBROUTINE create_empty_partlist



  SUBROUTINE create_unsafe_partlist(partlist, a_particle, n_elements)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle
    INTEGER(KIND=8), INTENT(IN) :: n_elements
    TYPE(particle), POINTER :: current
    INTEGER(KIND=8) :: ipart

    CALL create_empty_partlist(partlist)

    partlist%safe = .FALSE.
    current=>a_particle
    ipart = 1
    DO WHILE (ASSOCIATED(current) .AND. ipart .LT. n_elements)
      ipart = ipart+1
      current=>current%next
    ENDDO
    partlist%head=>a_particle
    partlist%tail=>current
    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist



  SUBROUTINE create_unsafe_partlist_by_tail(partlist, head, tail)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: head, tail
    TYPE(particle), POINTER :: current
    INTEGER(KIND=8) :: ipart

    CALL create_empty_partlist(partlist)

    partlist%safe = .FALSE.
    partlist%head=>head
    partlist%tail=>tail

    current=>head
    ipart = 0
    DO WHILE (ASSOCIATED(current))
      ipart = ipart+1
      current=>current%next
      IF (ASSOCIATED(current)) THEN
        IF (ASSOCIATED(current%prev, TARGET=tail)) EXIT
      ENDIF
    ENDDO

    partlist%count = ipart

  END SUBROUTINE create_unsafe_partlist_by_tail



  SUBROUTINE create_allocated_partlist(partlist, n_elements)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER(8), INTENT(IN) :: n_elements
    TYPE(particle), POINTER :: new_particle
    INTEGER(KIND=8) :: ipart

    CALL create_empty_partlist(partlist)

    DO ipart = 0, n_elements-1
      ALLOCATE(new_particle)
      CALL init_particle(new_particle)
      CALL add_particle_to_partlist(partlist, new_particle)
      NULLIFY(new_particle)
    ENDDO

  END SUBROUTINE create_allocated_partlist



  SUBROUTINE create_filled_partlist(partlist, data_in, n_elements)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: data_in
    INTEGER(8), INTENT(IN) :: n_elements
    TYPE(particle), POINTER :: new_particle
    INTEGER(KIND=8) :: ipart, cpos = 0

    CALL create_empty_partlist(partlist)

    DO ipart = 0, n_elements-1
      ALLOCATE(new_particle)
      cpos = ipart*nvar+1
      CALL unpack_particle(data_in(cpos:cpos+nvar-1), new_particle)
#ifdef PARTICLE_DEBUG
      new_particle%processor = rank
#endif
      CALL add_particle_to_partlist(partlist, new_particle)
      NULLIFY(new_particle)
    ENDDO

  END SUBROUTINE create_filled_partlist



  FUNCTION test_partlist(partlist)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: current
    INTEGER :: test_partlist
    INTEGER(KIND=8) :: test_ct

    test_partlist = 0
    test_ct = 0

    ! Empty list is OK
    IF (.NOT. ASSOCIATED(partlist%head) &
        .AND. .NOT. ASSOCIATED(partlist%tail)) THEN
      test_partlist = 0
      RETURN
    ENDIF

    ! List with head or tail but not both is broken
    IF (.NOT. ASSOCIATED(partlist%head) &
        .OR. .NOT. ASSOCIATED(partlist%tail)) THEN
      test_partlist = -1
      RETURN
    ENDIF

    ! Having head and tail elements which are not the end of a list are OK for
    ! unsafe partlists
    IF (ASSOCIATED(partlist%head%prev) .AND. partlist%safe) &
        test_partlist = IOR(test_partlist, 1)
    IF (ASSOCIATED(partlist%tail%next) .AND. partlist%safe) &
        test_partlist = IOR(test_partlist, 2)

    ! Since we don't KNOW that count is OK (that's what we're checking)
    ! Have to check both for end of list and for having reached the tail item
    current=>partlist%head
    DO WHILE (ASSOCIATED(current))
      test_ct = test_ct+1
      current=>current%next
      IF (ASSOCIATED(current)) THEN
        ! This tests if we've just jumped to the tail element
        ! Allows testing of unsafe partlists
        IF (ASSOCIATED(current%prev, TARGET=partlist%tail)) EXIT
      ENDIF
    ENDDO

    IF (test_ct .NE. partlist%count) test_partlist = IOR(test_partlist, 4)

  END FUNCTION test_partlist



  SUBROUTINE destroy_partlist(partlist)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: new_particle, next
    INTEGER(KIND=8) :: ipart

    ! Go through list and delete all the particles in the list
    new_particle=>partlist%head
    ipart = 0
    DO WHILE (ipart .LT. partlist%count)
      next=>new_particle%next
      DEALLOCATE(new_particle)
      new_particle=>next
      ipart = ipart+1
    ENDDO

    CALL create_empty_partlist(partlist)

  END SUBROUTINE destroy_partlist



  SUBROUTINE copy_partlist(partlist1, partlist2)

    TYPE(particle_list), INTENT(INOUT) :: partlist1, partlist2

    partlist2%head=>partlist1%head
    partlist2%tail=>partlist1%tail
    partlist2%count = partlist1%count
    partlist2%id_update = partlist1%id_update

  END SUBROUTINE copy_partlist



  SUBROUTINE append_partlist(head, tail)

    TYPE(particle_list), INTENT(INOUT) :: head, tail

    IF (.NOT. head%safe .OR. .NOT. tail%safe) THEN
      IF (rank .EQ. 0) &
          PRINT *, "Unable to append partlists because one is not safe"
      RETURN
    ENDIF

    IF (ASSOCIATED(head%tail)) THEN
      head%tail%next=>tail%head
    ELSE
      head%head=>tail%head
    ENDIF
    IF (ASSOCIATED(tail%head)) tail%head%prev=>head%tail
    IF (ASSOCIATED(tail%tail)) head%tail=>tail%tail
    head%count = head%count + tail%count
    head%id_update = head%id_update + tail%id_update

    CALL create_empty_partlist(tail)

  END SUBROUTINE append_partlist



  SUBROUTINE add_particle_to_partlist(partlist, new_particle)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: new_particle

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour

    ! if (!particle) return;
    IF (.NOT. ASSOCIATED(new_particle)) RETURN
    NULLIFY(new_particle%next, new_particle%prev)

    ! Add particle count
    partlist%count = partlist%count + 1
    partlist%id_update = 1
    IF (.NOT. ASSOCIATED(partlist%tail)) THEN
      ! partlist is empty
      partlist%head=>new_particle
      partlist%tail=>new_particle
      RETURN
    ENDIF

    partlist%tail%next=>new_particle
    new_particle%prev=>partlist%tail
    NULLIFY(new_particle%next)
    partlist%tail=>new_particle

  END SUBROUTINE add_particle_to_partlist



  SUBROUTINE remove_particle_from_partlist(partlist, a_particle)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    TYPE(particle), POINTER :: a_particle

    ! Note that this will work even if you are using an unsafe particle list
    ! BE CAREFUL if doing so, it can cause unexpected behaviour

    ! Check whether particle is head or tail of list and unlink
    IF (ASSOCIATED(partlist%head, TARGET=a_particle)) &
        partlist%head=>a_particle%next
    IF (ASSOCIATED(partlist%tail, TARGET=a_particle)) &
        partlist%tail=>a_particle%prev

    ! Link particles on either side together
    IF (ASSOCIATED(a_particle%next)) a_particle%next%prev=>a_particle%prev
    IF (ASSOCIATED(a_particle%prev)) a_particle%prev%next=>a_particle%next

    NULLIFY(a_particle%next, a_particle%prev)

    ! Decrement counter
    partlist%count = partlist%count-1

  END SUBROUTINE remove_particle_from_partlist



  SUBROUTINE pack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(KIND=8) :: cpos

    cpos = 1
    array(cpos:cpos+c_ndims-1) = a_particle%part_pos
    cpos = cpos+c_ndims
    array(cpos:cpos+2) = a_particle%part_p
    cpos = cpos+3
#ifdef PER_PARTICLE_WEIGHT
    array(cpos) = a_particle%weight
    cpos = cpos+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    array(cpos) = a_particle%charge
    array(cpos+1) = a_particle%mass
    cpos = cpos+2
#endif
#ifdef PARTICLE_DEBUG
    array(cpos) = REAL(a_particle%processor, num)
    array(cpos+1) = REAL(a_particle%processor_at_t0, num)
    cpos = cpos+2
#endif
#if PARTICLE_ID || PARTICLE_ID4
    array(cpos) = REAL(a_particle%id, num)
    cpos = cpos+1
#endif

!!$    PRINT *, "In Pack", rank, a_particle%part_pos

  END SUBROUTINE pack_particle



  SUBROUTINE unpack_particle(array, a_particle)

    REAL(num), DIMENSION(:), INTENT(IN) :: array
    TYPE(particle), POINTER :: a_particle
    INTEGER(KIND=8) :: cpos

    cpos = 1
    a_particle%part_pos = array(cpos)
    cpos = cpos+c_ndims
    a_particle%part_p = array(cpos:cpos+2)
    cpos = cpos+3
#ifdef PER_PARTICLE_WEIGHT
    a_particle%weight = array(cpos)
    cpos = cpos+1
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    a_particle%charge = array(cpos)
    a_particle%mass = array(cpos+1)
    cpos = cpos+2
#endif
#ifdef PARTICLE_DEBUG
    a_particle%processor = rank
    a_particle%processor_at_t0 = NINT(array(cpos+1))
    cpos = cpos+2
#endif
#ifdef PARTICLE_ID4
    a_particle%id = NINT(array(cpos))
    cpos = cpos+1
#elif PARTICLE_ID
    a_particle%id = NINT(array(cpos),8)
    cpos = cpos+1
#endif

  END SUBROUTINE unpack_particle



  SUBROUTINE init_particle(new_particle)

    TYPE(particle), POINTER :: new_particle

    new_particle%part_p = 0.0_num
    new_particle%part_pos = 0.0_num
#ifdef PER_PARTICLE_WEIGHT
    new_particle%weight = 0.0_num
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
    new_particle%charge = 0.0_num
    new_particle%mass = 0.0_num
#endif
#ifdef PARTICLE_DEBUG
    new_particle%processor = 0
    new_particle%processor_at_t0 = 0
#endif
#if PARTICLE_ID || PARTICLE_ID4
    new_particle%id = 0
#endif

  END SUBROUTINE init_particle



  SUBROUTINE display_particle(a_particle)

    TYPE(particle), POINTER :: a_particle

    PRINT *, "Position", a_particle%part_pos
    PRINT *, "Momentum", a_particle%part_p

  END SUBROUTINE display_particle



  FUNCTION compare_particles(part1, part2)

    TYPE(particle), POINTER :: part1, part2
    LOGICAL :: compare_particles

    compare_particles = .TRUE.
    IF ((ABS(part1%part_pos-part2%part_pos)) .NE. 0.0_num) &
        compare_particles = .FALSE.
    IF (MAXVAL(ABS(part1%part_p - part2%part_p)) .NE. 0.0_num) &
        compare_particles = .FALSE.

#ifdef PER_PARTICLE_WEIGHT
    IF (part1%weight .NE. part2%weight) compare_particles = .FALSE.
#endif

#ifdef PER_PARTICLE_CHARGE_MASS
    IF (part1%charge .NE. part2%charge) compare_particles = .FALSE.
    IF (part1%mass   .NE. part2%mass  ) compare_particles = .FALSE.
#endif

    IF (.NOT. compare_particles) THEN
      CALL display_particle(part1)
      CALL display_particle(part2)
    ENDIF

  END FUNCTION compare_particles



  FUNCTION test_packed_particles(partlist, array, npart_in_data)

    TYPE(particle_list), INTENT(IN) :: partlist
    REAL(num), DIMENSION(:), INTENT(IN) :: array
    INTEGER(KIND=8), INTENT(IN) :: npart_in_data
    TYPE(particle), POINTER :: current
    TYPE(particle), POINTER :: a_particle
    LOGICAL :: test_packed_particles
    INTEGER(KIND=8) :: ipart

    test_packed_particles = .FALSE.

    IF (npart_in_data * nvar .NE. SIZE(array)) THEN
      PRINT *, "Size of data array does not match specified on", rank, &
          npart_in_data, SIZE(array)
      RETURN
    ENDIF
    IF (partlist%count .NE. npart_in_data) THEN
      PRINT *, "Size of data array does not match partlist on", rank
      RETURN
    ENDIF

    ALLOCATE(a_particle)

    current=>partlist%head
    DO ipart = 0, npart_in_data-1
      CALL unpack_particle(array(ipart*nvar+1:(ipart+1)*nvar), a_particle)
      IF (.NOT. compare_particles(a_particle, current)) THEN
        PRINT *, "BAD PARTICLE ", ipart, "on", rank
        RETURN
      ENDIF
      current=>current%next
    ENDDO

    DEALLOCATE(a_particle)

    test_packed_particles = .TRUE.

  END FUNCTION test_packed_particles



  SUBROUTINE partlist_send(partlist, dest)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: dest
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    INTEGER(KIND=8) :: cpos = 0, npart_left, ipart, send_buf(2)
    TYPE(particle), POINTER :: current

    npart_left = partlist%count
    send_buf(1) = partlist%count
    send_buf(2) = partlist%id_update

    CALL MPI_SEND(send_buf, 2, MPI_INTEGER8, dest, tag, comm, errcode)

    ALLOCATE(array(1:partlist%count*nvar))
    array = 0.0_num
    current=>partlist%head
    ipart = 0
    DO WHILE (ipart .LT. partlist%count)
      cpos = ipart*nvar+1
      CALL pack_particle(array(cpos:cpos+nvar-1), current)
      ipart = ipart+1
      current=>current%next
    ENDDO
    CALL MPI_SEND(array, npart_left*nvar, mpireal, dest, tag, comm, errcode)

    DEALLOCATE(array)

  END SUBROUTINE partlist_send



  SUBROUTINE partlist_recv(partlist, src)

    TYPE(particle_list), INTENT(INOUT) :: partlist
    INTEGER, INTENT(IN) :: src
    REAL(num), DIMENSION(:), ALLOCATABLE :: array
    INTEGER(KIND=8) :: count, recv_buf(2)

    CALL create_empty_partlist(partlist)

    recv_buf = 0
    CALL MPI_RECV(recv_buf, 2, MPI_INTEGER8, src, tag, comm, status, errcode)
    count = recv_buf(1)
    partlist%id_update = partlist%id_update + recv_buf(2)

    ALLOCATE(array(1:count*nvar))
    array = 0.0_num
    CALL MPI_RECV(array, count*nvar, mpireal, src, tag, comm, status, errcode)
    CALL create_filled_partlist(partlist, array, count)

    DEALLOCATE(array)

  END SUBROUTINE partlist_recv



  SUBROUTINE partlist_sendrecv(partlist_send, partlist_recv, dest, src)

    TYPE(particle_list), INTENT(INOUT) :: partlist_send, partlist_recv
    INTEGER, INTENT(IN) :: dest, src
    REAL(num), DIMENSION(:), ALLOCATABLE :: data_send, data_recv, data_temp
    INTEGER(KIND=8) :: cpos = 0, ipart = 0
    INTEGER(KIND=8) :: npart_send, npart_recv, send_buf(2), recv_buf(2)
    TYPE(particle), POINTER :: current

    ! This subroutine doesn't try to use memory efficient buffering, it sends
    ! all the particles at once. This should work for boundary calls, but
    ! don't try it for any other reason

    recv_buf = 0
    send_buf(1) = partlist_send%count
    send_buf(2) = partlist_send%id_update

    CALL MPI_SENDRECV(send_buf, 2, MPI_INTEGER8, dest, tag, recv_buf, 2, &
        MPI_INTEGER8, src, tag, comm, status, errcode)

    npart_send = send_buf(1)
    npart_recv = recv_buf(1)
    partlist_recv%id_update = partlist_recv%id_update + recv_buf(2)

    ! Copy the data for the particles into a buffer
    ALLOCATE(data_send(1:npart_send*nvar))
    ALLOCATE(data_recv(1:npart_recv*nvar))
    ALLOCATE(data_temp(1:nvar))

    ! Pack particles to send into buffer
    current=>partlist_send%head
    ipart = 0
    DO WHILE (ipart .LT. partlist_send%count)
      cpos = ipart*nvar+1
      CALL pack_particle(data_temp, current)
      data_send(cpos:cpos+nvar-1) = data_temp
      ipart = ipart+1
      current=>current%next
    ENDDO

    ! No longer need the sending partlist, so destroy it to save some memory
    CALL destroy_partlist(partlist_send)

    ! Actual MPI commands
    CALL MPI_SENDRECV(data_send, npart_send*nvar, mpireal, dest, tag, &
        data_recv, npart_recv*nvar, mpireal, src, tag, comm, status, errcode)

    DEALLOCATE(data_send)
    CALL create_filled_partlist(partlist_recv, data_recv, npart_recv)
    DEALLOCATE(data_recv)

  END SUBROUTINE partlist_sendrecv



  SUBROUTINE add_particle_to_list(part, list)

    TYPE(particle), POINTER :: part
    TYPE(pointer_list) :: list
    TYPE(pointer_item), POINTER :: item

    ALLOCATE(item)
    item%part => part
    NULLIFY(item%next)

    list%tail%next => item
    list%tail => item

  END SUBROUTINE add_particle_to_list



  SUBROUTINE generate_particle_ids(partlist)

    TYPE(particle_list) :: partlist
#if PARTICLE_ID || PARTICLE_ID4
    INTEGER(KIND=8), ALLOCATABLE :: nid_all(:)
    INTEGER(KIND=8) :: nid, part_id
    INTEGER :: i, id_update
    TYPE(particle), POINTER :: current
    TYPE(pointer_list) :: idlist
    TYPE(pointer_item), POINTER :: idcurrent, idnext

    id_update = partlist%id_update

    CALL MPI_ALLREDUCE(id_update, partlist%id_update, 1, MPI_INTEGER, &
        MPI_MAX, comm, errcode)

    IF (partlist%id_update .EQ. 0) RETURN

    ALLOCATE(idlist%head)
    idlist%tail => idlist%head
    NULLIFY(idlist%head%next)
    NULLIFY(idlist%head%part)

    ! Scan through particle list and identify particles which need
    ! an ID to be assigned.
    nid = 0
    current => partlist%head
    DO WHILE(ASSOCIATED(current))
      IF (current%id .EQ. 0) THEN
        nid = nid + 1
        CALL add_particle_to_list(current, idlist)
      ENDIF
      current => current%next
    ENDDO

    ALLOCATE(nid_all(nproc))

    CALL MPI_ALLGATHER(nid, 1, MPI_INTEGER8, nid_all, 1, MPI_INTEGER8, &
        comm, errcode)

    ! Count number of particles on ranks zero to rank-1
    nid = 0
    DO i = 1, rank
      nid = nid + nid_all(i)
    ENDDO
    part_id = particles_max_id + nid

    ! Count remaining particles
    DO i = rank+1, nproc
      nid = nid + nid_all(i)
    ENDDO

    particles_max_id = particles_max_id + nid

    DEALLOCATE(nid_all)

    ! Number each particle with a unique id
    idcurrent => idlist%head%next
    DO WHILE(ASSOCIATED(idcurrent))
      part_id = part_id + 1
      idcurrent%part%id = part_id
      idnext => idcurrent%next
      DEALLOCATE(idcurrent)
      idcurrent => idnext
    ENDDO

    DEALLOCATE(idlist%head)

    partlist%id_update = 0
#endif

  END SUBROUTINE generate_particle_ids

END MODULE partlist
