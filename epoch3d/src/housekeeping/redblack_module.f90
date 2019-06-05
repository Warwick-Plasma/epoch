! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE redblack_module

  USE partlist

  IMPLICIT NONE

  REAL(num), DIMENSION(:), POINTER :: field_in1d, field_out1d
  REAL(num), DIMENSION(:,:), POINTER :: field_in2d, field_out2d
  REAL(num), DIMENSION(:,:,:), POINTER :: field_in3d, field_out3d
  REAL(r4), DIMENSION(:), POINTER :: field_in1dr4, field_out1dr4
  REAL(r4), DIMENSION(:,:), POINTER :: field_in2dr4, field_out2dr4
  REAL(r4), DIMENSION(:,:,:), POINTER :: field_in3dr4, field_out3dr4
  INTEGER, DIMENSION(:), POINTER :: sendtypes, recvtypes
  TYPE(particle_list), DIMENSION(:), POINTER :: pointers_send, pointers_recv
  INTEGER(i8), DIMENSION(:), POINTER :: sendcounts, recvcounts

  INTERFACE redblack
    MODULE PROCEDURE &
        redblackpart, &
        redblack1d, &
        redblack2d, &
        redblack3d, &
        redblack1d_r4, &
        redblack2d_r4, &
        redblack3d_r4
  END INTERFACE redblack

CONTAINS

  SUBROUTINE do_sendpart(iproc)

    INTEGER, INTENT(IN) :: iproc

    IF (sendcounts(iproc) > 0) THEN
      CALL partlist_send_nocount(pointers_send(iproc), iproc)
      CALL destroy_partlist(pointers_send(iproc))
    END IF

  END SUBROUTINE do_sendpart



  SUBROUTINE do_recvpart(iproc)

    INTEGER, INTENT(IN) :: iproc

    IF (recvcounts(iproc) > 0) THEN
      CALL partlist_recv_nocount(pointers_recv(iproc), iproc, &
          recvcounts(iproc))
    END IF

  END SUBROUTINE do_recvpart



  SUBROUTINE do_send1d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) /= 0) THEN
      CALL MPI_SEND(field_in1d, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    END IF

  END SUBROUTINE do_send1d



  SUBROUTINE do_recv1d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) /= 0) THEN
      CALL MPI_RECV(field_out1d, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    END IF

  END SUBROUTINE do_recv1d



  SUBROUTINE do_send2d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) /= 0) THEN
      CALL MPI_SEND(field_in2d, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    END IF

  END SUBROUTINE do_send2d



  SUBROUTINE do_recv2d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) /= 0) THEN
      CALL MPI_RECV(field_out2d, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    END IF

  END SUBROUTINE do_recv2d



  SUBROUTINE do_send3d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) /= 0) THEN
      CALL MPI_SEND(field_in3d, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    END IF

  END SUBROUTINE do_send3d



  SUBROUTINE do_recv3d(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) /= 0) THEN
      CALL MPI_RECV(field_out3d, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    END IF

  END SUBROUTINE do_recv3d



  SUBROUTINE do_send1dr4(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) /= 0) THEN
      CALL MPI_SEND(field_in1dr4, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    END IF

  END SUBROUTINE do_send1dr4



  SUBROUTINE do_recv1dr4(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) /= 0) THEN
      CALL MPI_RECV(field_out1dr4, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    END IF

  END SUBROUTINE do_recv1dr4



  SUBROUTINE do_send2dr4(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) /= 0) THEN
      CALL MPI_SEND(field_in2dr4, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    END IF

  END SUBROUTINE do_send2dr4



  SUBROUTINE do_recv2dr4(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) /= 0) THEN
      CALL MPI_RECV(field_out2dr4, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    END IF

  END SUBROUTINE do_recv2dr4



  SUBROUTINE do_send3dr4(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (sendtypes(iproc) /= 0) THEN
      CALL MPI_SEND(field_in3dr4, 1, sendtypes(iproc), iproc, tag, comm, ierr)
    END IF

  END SUBROUTINE do_send3dr4



  SUBROUTINE do_recv3dr4(iproc)

    INTEGER, INTENT(IN) :: iproc
    INTEGER :: ierr

    IF (recvtypes(iproc) /= 0) THEN
      CALL MPI_RECV(field_out3dr4, 1, recvtypes(iproc), iproc, tag, comm, &
          MPI_STATUS_IGNORE, ierr)
    END IF

  END SUBROUTINE do_recv3dr4



  SUBROUTINE redblackpart(psend, precv, sendcounts_in, recvcounts_in)

    TYPE(particle_list), DIMENSION(0:), TARGET :: psend, precv
    INTEGER(i8), DIMENSION(0:), TARGET :: sendcounts_in, recvcounts_in

    pointers_send => psend
    pointers_recv => precv
    sendcounts => sendcounts_in
    recvcounts => recvcounts_in
    CALL redblack_main(do_sendpart, do_recvpart)

  END SUBROUTINE redblackpart



  SUBROUTINE redblack1d(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(num), DIMENSION(1-ng:), TARGET, INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in1d => field_in
    field_out1d => field_out
    CALL redblack_main(do_send1d, do_recv1d)

  END SUBROUTINE redblack1d



  SUBROUTINE redblack2d(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(num), DIMENSION(1-ng:,1-ng:), TARGET, INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:,1-ng:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in2d => field_in
    field_out2d => field_out
    CALL redblack_main(do_send2d, do_recv2d)

  END SUBROUTINE redblack2d



  SUBROUTINE redblack3d(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), TARGET, INTENT(IN) :: field_in
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in3d => field_in
    field_out3d => field_out
    CALL redblack_main(do_send3d, do_recv3d)

  END SUBROUTINE redblack3d



  SUBROUTINE redblack1d_r4(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(r4), DIMENSION(1-ng:), TARGET, INTENT(IN) :: field_in
    REAL(r4), DIMENSION(1-ng:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in1dr4 => field_in
    field_out1dr4 => field_out
    CALL redblack_main(do_send1dr4, do_recv1dr4)

  END SUBROUTINE redblack1d_r4



  SUBROUTINE redblack2d_r4(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(r4), DIMENSION(1-ng:,1-ng:), TARGET, INTENT(IN) :: field_in
    REAL(r4), DIMENSION(1-ng:,1-ng:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in2dr4 => field_in
    field_out2dr4 => field_out
    CALL redblack_main(do_send2dr4, do_recv2dr4)

  END SUBROUTINE redblack2d_r4



  SUBROUTINE redblack3d_r4(field_in, field_out, sendtypes_in, recvtypes_in)

    REAL(r4), DIMENSION(1-ng:,1-ng:,1-ng:), TARGET, INTENT(IN) :: field_in
    REAL(r4), DIMENSION(1-ng:,1-ng:,1-ng:), TARGET, INTENT(OUT) :: field_out
    INTEGER, DIMENSION(0:), TARGET, INTENT(INOUT) :: sendtypes_in, recvtypes_in

    sendtypes => sendtypes_in
    recvtypes => recvtypes_in
    field_in3dr4 => field_in
    field_out3dr4 => field_out
    CALL redblack_main(do_send3dr4, do_recv3dr4)

  END SUBROUTINE redblack3d_r4



  SUBROUTINE redblack_main(do_send, do_recv)

    INTERFACE
      SUBROUTINE do_send(iproc)
        INTEGER, INTENT(IN) :: iproc
      END SUBROUTINE do_send

      SUBROUTINE do_recv(iproc)
        INTEGER, INTENT(IN) :: iproc
      END SUBROUTINE do_recv
    END INTERFACE

    INTEGER, DIMENSION(:), ALLOCATABLE :: oddlist, evnlist
    INTEGER :: nodd, nevn, iproc, i
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
    END DO
    IF (nevn /= nodd) evnlist(nevn-1) = nproc - 1

    DO WHILE(nodd > 0)
      is_evn = .TRUE.
      DO i = 0, nodd - 1
        IF (oddlist(i) == rank) THEN
          is_evn = .FALSE.
          EXIT
        END IF
      END DO

      IF (is_evn) THEN
        DO i = 0, nodd - 1
          iproc = oddlist(i)
          CALL do_send(iproc)
        END DO

        DO i = 0, nodd - 1
          iproc = oddlist(i)
          CALL do_recv(iproc)
        END DO

        nodd = nevn / 2
        nevn = (nevn + 1) / 2
        DO i = 0, nodd - 1
          evnlist(i) = evnlist(2*i)
          oddlist(i) = evnlist(2*i+1)
        END DO
        IF (nevn /= nodd) evnlist(nevn-1) = evnlist(2*nevn-2)
      ELSE
        DO i = 0, nevn - 1
          iproc = evnlist(i)
          CALL do_recv(iproc)
        END DO

        DO i = 0, nevn - 1
          iproc = evnlist(i)
          CALL do_send(iproc)
        END DO

        nevn = (nodd + 1) / 2
        nodd = nodd / 2
        DO i = 0, nodd - 1
          evnlist(i) = oddlist(2*i)
          oddlist(i) = oddlist(2*i+1)
        END DO
        IF (nevn /= nodd) evnlist(nevn-1) = oddlist(2*nevn-2)
      END IF
    END DO

    DEALLOCATE(evnlist)
    DEALLOCATE(oddlist)

  END SUBROUTINE redblack_main

END MODULE redblack_module
