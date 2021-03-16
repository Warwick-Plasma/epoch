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

MODULE deck_dist_fn_block

  USE strings_advanced
  USE dist_fn
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: dist_fn_deck_initialise, dist_fn_deck_finalise
  PUBLIC :: dist_fn_block_start, dist_fn_block_end
  PUBLIC :: dist_fn_block_handle_element, dist_fn_block_check

  TYPE(distribution_function_block), POINTER :: working_block
  LOGICAL :: got_name
  INTEGER :: ndims

CONTAINS

  SUBROUTINE dist_fn_deck_initialise

  END SUBROUTINE dist_fn_deck_initialise



  SUBROUTINE dist_fn_deck_finalise

  END SUBROUTINE dist_fn_deck_finalise



  SUBROUTINE dist_fn_block_start

    IF (deck_state == c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_block)
    CALL init_dist_fn(working_block)
    ndims = 0
    got_name = .FALSE.

  END SUBROUTINE dist_fn_block_start



  SUBROUTINE dist_fn_block_end

    INTEGER :: i, dir, n, iu, io
    REAL(num) :: r1, r2
    REAL(num), PARAMETER :: pi2 = 2.0_num * pi

    IF (deck_state == c_ds_first) RETURN

    IF (.NOT.got_name) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'name not set for "dist_fn" block.'
        END DO
      END IF
      CALL abort_code(c_err_missing_elements)
      RETURN
    END IF

    IF (working_block%ndims == -1) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'ndims not set for "dist_fn" block "' &
              // TRIM(working_block%name) // '"'
        END DO
      END IF
      CALL abort_code(c_err_missing_elements)
      RETURN
    END IF

    IF (ndims > working_block%ndims) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The parameters in "dist_fn" block "' &
              // TRIM(working_block%name) // '"'
          WRITE(io,*) 'exceed the number of dimensions for the ', &
              'distribution function.'
        END DO
      END IF
      CALL abort_code(c_err_bad_value)
      RETURN
    END IF

    DO i = 1, working_block%ndims
      dir = working_block%directions(i)
      IF (dir /= c_dir_xy_angle &
          .AND. dir /= c_dir_yz_angle .AND. dir /= c_dir_zx_angle) CYCLE

      r1 = working_block%ranges(1,i)
      r2 = working_block%ranges(2,i)
      IF (ABS(r1 - r2) <= c_tiny) CYCLE

      ! If direction is an angle, set start angle to lie in the range [-pi,pi)
      n = INT(r1 / pi2)
      r1 = r1 - pi2 * n
      IF (r1 >=  pi) r1 = r1 - pi2
      IF (r1 < -pi) r1 = r1 + pi2
      working_block%ranges(1,i) = r1

      ! Set end angle to be less than 2*pi greater than start angle
      n = INT(r2 / pi2)
      r2 = r2 - pi2 * n
      IF (r2 >=  pi) r2 = r2 - pi2
      IF (r2 < -pi) r2 = r2 + pi2
      IF (r2 <=  r1) r2 = r2 + pi2
      working_block%ranges(2,i) = r2
    END DO

    CALL attach_dist_fn(working_block)

  END SUBROUTINE dist_fn_block_end



  FUNCTION dist_fn_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    CHARACTER(LEN=string_length) :: part1
    INTEGER :: part2, ispecies
    INTEGER :: work, io, iu
    REAL(num) :: work1, work2

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'name')) THEN
      working_block%name = value
      got_name = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'ndims')) THEN
      work = as_integer_print(value, element, errcode)
      IF (work >= 1 .AND. work <= 3) THEN
        working_block%ndims = work
      ELSE
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
            WRITE(io,*) 'Distribution functions can only be 1D, 2D or 3D'
          END DO
        END IF
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'dumpmask')) THEN
      working_block%dumpmask = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'output_deltaf')) THEN
      working_block%output_deltaf = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_x')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_x) = .TRUE.
      working_block%restrictions(:,c_dir_x) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_y')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_y) = .TRUE.
      working_block%restrictions(:,c_dir_y) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_z')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_z) = .TRUE.
      working_block%restrictions(:,c_dir_z) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_px')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_px) = .TRUE.
      working_block%restrictions(:,c_dir_px) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_py')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_py) = .TRUE.
      working_block%restrictions(:,c_dir_py) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_pz')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_pz) = .TRUE.
      working_block%restrictions(:,c_dir_pz) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_en') &
        .OR. str_cmp(element, 'restrict_energy')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_en) = .TRUE.
      working_block%restrictions(:,c_dir_en) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_gamma_m1') &
        .OR. str_cmp(element, 'restrict_gamma_minus_one')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_gamma_m1) = .TRUE.
      working_block%restrictions(:,c_dir_gamma_m1) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_xy_angle')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_xy_angle) = .TRUE.
      working_block%restrictions(:,c_dir_xy_angle) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_yz_angle')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_yz_angle) = .TRUE.
      working_block%restrictions(:,c_dir_yz_angle) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_zx_angle')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_zx_angle) = .TRUE.
      working_block%restrictions(:,c_dir_zx_angle) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'restrict_mod_p')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode /= c_err_none) RETURN
      working_block%use_restrictions(c_dir_mod_p) = .TRUE.
      working_block%restrictions(:,c_dir_mod_p) = (/work1, work2/)
      RETURN
    END IF

    IF (str_cmp(element, 'include_species')) THEN
      ispecies = as_integer_print(value, element, errcode)
      IF (errcode == c_err_none) THEN
        IF (ispecies > 0 .AND. ispecies <= n_species) THEN
          working_block%use_species(ispecies) = .TRUE.
        ELSE
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
              WRITE(io,*) 'Unable to apply dist_fn to non existant species ', &
                  ispecies
            END DO
          END IF
          errcode = c_err_bad_value
        END IF
      END IF
      RETURN
    END IF

    CALL split_off_int(element, part1, part2, errcode)
    IF (part2 > ndims) ndims = part2

    IF (errcode /= c_err_none) THEN
      errcode = c_err_unknown_element
      RETURN
    END IF

    IF (str_cmp(part1, 'direction')) THEN
      working_block%directions(part2) = &
          as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(part1, 'range')) THEN
      CALL split_range(TRIM(value), work1, work2, errcode)
      IF (IAND(errcode, c_err_bad_value) /= 0) THEN
        errcode = IAND(errcode, NOT(c_err_bad_value))
        errcode = IOR(errcode, c_err_warn_bad_value)
        RETURN
      END IF
      working_block%ranges(1,part2) = work1
      working_block%ranges(2,part2) = work2
      RETURN
    END IF

    IF (str_cmp(part1, 'resolution')) THEN
      working_block%resolution(part2) = &
          as_integer_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION dist_fn_block_handle_element



  FUNCTION dist_fn_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION dist_fn_block_check

END MODULE deck_dist_fn_block
