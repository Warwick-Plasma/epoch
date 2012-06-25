MODULE redblack_module

  USE shared_data
  USE partlist

  IMPLICIT NONE

  REAL(num), DIMENSION(:), POINTER :: field_in1d, field_out1d
  REAL(num), DIMENSION(:,:), POINTER :: field_in2d, field_out2d
  REAL(num), DIMENSION(:,:,:), POINTER :: field_in3d, field_out3d
  INTEGER, DIMENSION(:), POINTER :: sendtypes, recvtypes
  TYPE(particle_list), DIMENSION(:), POINTER :: pointers_send, pointers_recv
  INTEGER(KIND=8), DIMENSION(:), POINTER :: sendcounts, recvcounts

  INTERFACE redblack
    MODULE PROCEDURE &
        redblackpart, &
        redblack1d, &
        redblack2d, &
        redblack3d
  END INTERFACE redblack

CONTAINS

  SUBROUTINE do_sendpart(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendcounts(iproc) .GT. 0) THEN
      CALL partlist_send_nocount(pointers_send(iproc), iproc)
      CALL destroy_partlist(pointers_send(iproc))
    ENDIF

  END SUBROUTINE do_sendpart



  SUBROUTINE do_recvpart(iproc)

    USE mpi
    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvcounts(iproc) .GT. 0) THEN
      CALL partlist_recv_nocount(pointers_recv(iproc), iproc, &
          recvcounts(iproc))
    ENDIF

  END SUBROUTINE do_recvpart



  SUBROUTINE do_send1d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) .NE. 0) THEN
      CALL MPI_SEND(field_in1d, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    ENDIF

  END SUBROUTINE do_send1d



  SUBROUTINE do_recv1d(iproc)

    USE mpi
    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) .NE. 0) THEN
      CALL MPI_RECV(field_out1d, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    ENDIF

  END SUBROUTINE do_recv1d



  SUBROUTINE do_send2d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) .NE. 0) THEN
      CALL MPI_SEND(field_in2d, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    ENDIF

  END SUBROUTINE do_send2d



  SUBROUTINE do_recv2d(iproc)

    USE mpi
    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) .NE. 0) THEN
      CALL MPI_RECV(field_out2d, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    ENDIF

  END SUBROUTINE do_recv2d



  SUBROUTINE do_send3d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) .NE. 0) THEN
      CALL MPI_SEND(field_in3d, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    ENDIF

  END SUBROUTINE do_send3d



  SUBROUTINE do_recv3d(iproc)

    USE mpi
    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) .NE. 0) THEN
      CALL MPI_RECV(field_out3d, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    ENDIF

  END SUBROUTINE do_recv3d



  SUBROUTINE redblackpart(psend, precv, sendcounts_in, recvcounts_in)

    TYPE(particle_list), DIMENSION(:), TARGET :: psend, precv
    INTEGER(KIND=8), DIMENSION(:), TARGET :: sendcounts_in, recvcounts_in

    pointers_send => psend
    pointers_recv => precv
    sendcounts => sendcounts_in
    recvcounts => recvcounts_in
    CALL redblack_main(do_sendpart, do_recvpart)

  END SUBROUTINE redblackpart



  SUBROUTINE redblack1d(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(num), DIMENSION(-2:), TARGET, INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in1d => field_in
    field_out1d => field_out
    CALL redblack_main(do_send1d, do_recv1d)

  END SUBROUTINE redblack1d



  SUBROUTINE redblack2d(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(num), DIMENSION(-2:,-2:), TARGET, INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:,-2:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in2d => field_in
    field_out2d => field_out
    CALL redblack_main(do_send2d, do_recv2d)

  END SUBROUTINE redblack2d



  SUBROUTINE redblack3d(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(num), DIMENSION(-2:,-2:,-2:), TARGET, INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:,-2:,-2:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in3d => field_in
    field_out3d => field_out
    CALL redblack_main(do_send3d, do_recv3d)

  END SUBROUTINE redblack3d



  SUBROUTINE redblack_main(do_send, do_recv)

    ! This subroutine redistributes the fields over the new processor layout
    INTEGER, PARAMETER :: ng = 3

    INTERFACE
      SUBROUTINE do_send(iproc)
        INTEGER, INTENT(IN) :: iproc
      END SUBROUTINE do_send

      SUBROUTINE do_recv(iproc)
        INTEGER, INTENT(IN) :: iproc
      END SUBROUTINE do_recv
    END INTERFACE

    INTEGER, DIMENSION(:), ALLOCATABLE :: oddlist, evnlist
    INTEGER :: nodd, nevn, iproc, i, ierr
    LOGICAL :: is_evn

    ! Split the processors into even and odd lists.
    ! Even processors send then receive, Odd processors vice-versa.
    ! Next, on even processors split the even list in two and on
    ! odd processors split the odd list in two.
    ! Repeat until the list size is equal to one.

    ! If the number of processors is not divisible by two then the
    ! even list has one extra entry.

    nevn = (nproc + 1) / 2
    nodd = nproc / 2
    ALLOCATE(evnlist(0:nevn-1), oddlist(0:nodd-1))

    DO i = 0, nodd - 1
      evnlist(i) = 2 * i
      oddlist(i) = 2 * i + 1
    ENDDO
    IF (nevn .NE. nodd) evnlist(nevn-1) = nproc - 1

    DO WHILE(nodd .GT. 0)
      is_evn = .TRUE.
      DO i = 0, nodd - 1
        IF (oddlist(i) .EQ. rank) THEN
          is_evn = .FALSE.
          EXIT
        ENDIF
      ENDDO

      IF (is_evn) THEN
        DO i = 0, nodd - 1
          iproc = oddlist(i)
          CALL do_send(iproc)
        ENDDO

        DO i = 0, nodd - 1
          iproc = oddlist(i)
          CALL do_recv(iproc)
        ENDDO

        nodd = nevn / 2
        nevn = (nevn + 1) / 2
        DO i = 0, nodd - 1
          evnlist(i) = evnlist(2*i)
          oddlist(i) = evnlist(2*i+1)
        ENDDO
        IF (nevn .NE. nodd) evnlist(nevn-1) = evnlist(2*nevn-2)
      ELSE
        DO i = 0, nevn - 1
          iproc = evnlist(i)
          CALL do_recv(iproc)
        ENDDO

        DO i = 0, nevn - 1
          iproc = evnlist(i)
          CALL do_send(iproc)
        ENDDO

        nevn = (nodd + 1) / 2
        nodd = nodd / 2
        DO i = 0, nodd - 1
          evnlist(i) = oddlist(2*i)
          oddlist(i) = oddlist(2*i+1)
        ENDDO
        IF (nevn .NE. nodd) evnlist(nevn-1) = oddlist(2*nevn-2)
      ENDIF
    ENDDO

    DEALLOCATE(evnlist)
    DEALLOCATE(oddlist)

  END SUBROUTINE redblack_main

END MODULE redblack_module
