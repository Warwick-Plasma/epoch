MODULE deck_dist_fn_block

  USE strings_advanced
  USE dist_fn

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

    IF (deck_state .EQ. c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_block)
    CALL init_dist_fn(working_block)
    ndims = 0
    got_name = .FALSE.

  END SUBROUTINE dist_fn_block_start



  SUBROUTINE dist_fn_block_end

    INTEGER :: i, dir, n, iu, io, ierr
    REAL(num) :: r1, r2
    REAL(num), PARAMETER :: pi2 = 2.0_num * pi

    IF (deck_state .EQ. c_ds_first) RETURN

    IF (.NOT.got_name) THEN
      IF (rank .EQ. 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'name not set for "dist_fn" block.'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      RETURN
    ENDIF

    IF (working_block%ndims .EQ. -1) THEN
      IF (rank .EQ. 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'ndims not set for "dist_fn" block "' &
              // TRIM(working_block%name) // '"'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      RETURN
    ENDIF

    IF (ndims .GT. working_block%ndims) THEN
      IF (rank .EQ. 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The parameters in "dist_fn" block "' &
              // TRIM(working_block%name) // '"'
          WRITE(io,*) 'exceed the number of dimensions for the ', &
              'distribution function.'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      RETURN
    ENDIF

    DO i = 1, working_block%ndims
      dir = working_block%directions(i)
      IF (dir .NE. c_dir_xy_angle &
          .AND. dir .NE. c_dir_yz_angle .AND. dir .NE. c_dir_zx_angle) CYCLE

      r1 = working_block%ranges(1,i)
      r2 = working_block%ranges(2,i)
      IF (r1 .EQ. r2) CYCLE

      ! If direction is an angle, set start angle to lie in the range [-pi,pi)
      n = INT(r1 / pi2)
      r1 = r1 - pi2 * n
      IF (r1 .GE.  pi) r1 = r1 - pi2
      IF (r1 .LT. -pi) r1 = r1 + pi2
      working_block%ranges(1,i) = r1

      ! Set end angle to be less than 2*pi greater than start angle
      n = INT(r2 / pi2)
      r2 = r2 - pi2 * n
      IF (r2 .GE.  pi) r2 = r2 - pi2
      IF (r2 .LT. -pi) r2 = r2 + pi2
      IF (r2 .LE.  r1) r2 = r2 + pi2
      working_block%ranges(2,i) = r2
    ENDDO

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
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, 'name')) THEN
      working_block%name = value
      got_name = .TRUE.
      RETURN
    ENDIF

    IF (str_cmp(element, 'ndims')) THEN
      work = as_integer(value, errcode)
      IF (work .GE. 1 .AND. work .LE. 3) THEN
        working_block%ndims = work
      ELSE
        IF (rank .EQ. 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Distribution functions can only be 1D, 2D or 3D'
          ENDDO
        ENDIF
        errcode = c_err_bad_value
      ENDIF
      RETURN
    ENDIF

    IF (str_cmp(element, 'dumpmask')) THEN
      working_block%dumpmask = as_integer(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'restrict_x')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(c_dir_x) = .TRUE.
      working_block%restrictions(:,c_dir_x) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_y')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(c_dir_y) = .TRUE.
      working_block%restrictions(:,c_dir_y) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_z')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(c_dir_z) = .TRUE.
      working_block%restrictions(:,c_dir_z) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_px')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(c_dir_px) = .TRUE.
      working_block%restrictions(:,c_dir_px) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_py')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(c_dir_py) = .TRUE.
      working_block%restrictions(:,c_dir_py) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_pz')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(c_dir_pz) = .TRUE.
      working_block%restrictions(:,c_dir_pz) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'include_species')) THEN
      ispecies = as_integer(value, errcode)
      IF (errcode .EQ. c_err_none) THEN
        IF (ispecies .GT. 0 .AND. ispecies .LE. n_species) THEN
          working_block%use_species(ispecies) = .TRUE.
        ELSE
          IF (rank .EQ. 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Unable to apply dist_fn to non existant species ', &
                  ispecies
            ENDDO
          ENDIF
          errcode = c_err_bad_value
        ENDIF
      ENDIF
      RETURN
    ENDIF

    CALL split_off_int(element, part1, part2, errcode)
    IF (part2 .GT. ndims) ndims = part2

    IF (errcode .NE. c_err_none) THEN
      errcode = c_err_unknown_element
      RETURN
    ENDIF

    IF (str_cmp(part1, 'direction')) THEN
      working_block%directions(part2) = as_integer(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(part1, 'range')) THEN
      CALL split_range(TRIM(value), work1, work2, errcode)
      IF (IAND(errcode, c_err_bad_value) .NE. 0) THEN
        errcode = IAND(errcode, NOT(c_err_bad_value))
        errcode = IOR(errcode, c_err_warn_bad_value)
        RETURN
      ENDIF
      working_block%ranges(1,part2) = work1
      working_block%ranges(2,part2) = work2
      RETURN
    ENDIF

    IF (str_cmp(part1, 'resolution')) THEN
      working_block%resolution(part2) = as_integer(value, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION dist_fn_block_handle_element



  FUNCTION dist_fn_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION dist_fn_block_check

END MODULE deck_dist_fn_block
