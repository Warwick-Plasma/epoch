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

MODULE deck_particle_probe_block

#ifdef NO_PARTICLE_PROBES
  IMPLICIT NONE

CONTAINS

  SUBROUTINE probe_deck_dummy

  END SUBROUTINE probe_deck_dummy

#else
  USE strings_advanced
  USE probes
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: probe_deck_initialise, probe_deck_finalise
  PUBLIC :: probe_block_start, probe_block_end
  PUBLIC :: probe_block_handle_element, probe_block_check

  TYPE(particle_probe), POINTER :: working_probe
  REAL(num) :: point2(c_ndims)
  LOGICAL :: got_name, got_point, got_normal
  INTEGER :: got_x
  INTEGER, PARAMETER :: ndim = 4
  CHARACTER(LEN=*), PARAMETER :: xs(ndim) = (/'x1', 'y1', 'x2', 'y2'/)

CONTAINS

  SUBROUTINE probe_deck_initialise

  END SUBROUTINE probe_deck_initialise



  SUBROUTINE probe_deck_finalise

  END SUBROUTINE probe_deck_finalise



  SUBROUTINE probe_block_start

    IF (deck_state == c_ds_first) RETURN

    ALLOCATE(working_probe)
    CALL init_probe(working_probe)
    got_name = .FALSE.
    got_point = .FALSE.
    got_normal = .FALSE.
    got_x = 0

  END SUBROUTINE probe_block_start



  SUBROUTINE probe_block_end

    LOGICAL :: discard
    REAL(num), DIMENSION(c_ndims) :: r1
    INTEGER :: io, iu, i, scount, sarr(ndim)

    IF (deck_state == c_ds_first) RETURN

    IF (.NOT.got_name) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) '"probe" block does not have a "name" entry.'
        END DO
      END IF
      CALL abort_code(c_err_required_element_not_set)
    END IF

    discard = .FALSE.
    IF (got_point) THEN
      IF (rank == 0) THEN
        IF (got_x /= 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Both "x1", etc. and "point" were used in probe ', &
                'block, "' // TRIM(working_probe%name) // '".'
            WRITE(io,*) 'Only "point" and "normal" will be used.'
          END DO
        END IF
      END IF
      IF (.NOT. got_normal) discard = .TRUE.
    ELSE
      IF (got_x /= 2**ndim-1) THEN
        discard = .TRUE.
      ELSE
        ! Old style configuration supplied. Need to calculate the normal.
        ! The probe calculates the signed distance from a point to a plane
        ! using Hessian normal form
        r1 = point2 - working_probe%point
        ! r1 (cross) z
        working_probe%normal = (/ r1(2), -r1(1) /)

        IF (SUM(ABS(working_probe%normal)) <= c_tiny) discard = .TRUE.
      END IF
    END IF

    IF (discard) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Position of probe "' // TRIM(working_probe%name) &
              // '" ', 'not fully specified. ', 'It will be discarded.'
          IF (got_point .OR. got_normal .OR. got_x == 0) THEN
            WRITE(io,*) 'Both "point" and "normal" are required.'
          ELSE
            scount = 0
            sarr = 0
            DO i = 0,ndim-1
              IF (IAND(got_x,2**i) == 0) THEN
                scount = scount + 1
                sarr(scount) = i + 1
              END IF
            END DO
            IF (scount > 1) THEN
              DO i = 1, scount-2
                WRITE(io,'(A)',ADVANCE='NO') ' "' // xs(sarr(i)) // '",'
              END DO
              WRITE(io,*) '"' // xs(sarr(scount-1)) // '" and "' &
                  // xs(sarr(scount)) // '" not specified.'
            ELSE
              WRITE(io,*) '"' // xs(sarr(scount)) // '" not specified.'
            END IF
          END IF
        END DO
      END IF

      DEALLOCATE(working_probe)
      NULLIFY(working_probe)
    ELSE
      ! Normalise the normal. Not really necessary but doesn't hurt.
      working_probe%normal = &
          working_probe%normal / SQRT(SUM(working_probe%normal**2))

      CALL attach_probe(working_probe)
    END IF

  END SUBROUTINE probe_block_end



  FUNCTION probe_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, ispecies, io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    ! get particle probe diagnostics (rolling total of all particles which
    ! pass through a given region of real space (defined by a point on a plane
    ! and the normal to that plane.
    IF (str_cmp(element, 'dumpmask') .OR. str_cmp(element, 'dump')) THEN
      working_probe%dumpmask = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'point') .OR. str_cmp(element, 'probe_point')) THEN
      got_point = .TRUE.
      CALL get_vector(value, working_probe%point, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'normal')) THEN
      got_normal = .TRUE.
      CALL get_vector(value, working_probe%normal, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'x1')) THEN
      got_x = IOR(got_x,2**0)
      working_probe%point(1) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'y1')) THEN
      got_x = IOR(got_x,2**1)
      working_probe%point(2) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'x2')) THEN
      got_x = IOR(got_x,2**2)
      point2(1) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'y2')) THEN
      got_x = IOR(got_x,2**3)
      point2(2) = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'include_species') &
        .OR. str_cmp(element, 'probe_species')) THEN
      ispecies = as_integer_print(value, element, errcode)
      IF (errcode == c_err_none) THEN
        IF (ispecies > 0 .AND. ispecies <= n_species) THEN
          working_probe%use_species(ispecies) = .TRUE.
        ELSE
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
              WRITE(io,*) 'Unable to attach probe to non existant species ', &
                  ispecies
            END DO
          END IF
          errcode = c_err_bad_value
        END IF
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'ek_min')) THEN
      working_probe%ek_min = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'ek_max')) THEN
      working_probe%ek_max = as_real_print(value, element, errcode)
      IF (working_probe%ek_max < 0) working_probe%ek_max = HUGE(1.0_num)
      RETURN
    END IF

    IF (str_cmp(element, 'name')) THEN
      got_name = .TRUE.
      working_probe%name = TRIM(value)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION probe_block_handle_element



  FUNCTION probe_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION probe_block_check
#endif

END MODULE deck_particle_probe_block
