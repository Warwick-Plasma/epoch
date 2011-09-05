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

  END SUBROUTINE dist_fn_block_start



  SUBROUTINE dist_fn_block_end

    IF (deck_state .EQ. c_ds_first) RETURN

    CALL attach_dist_fn(working_block)

  END SUBROUTINE dist_fn_block_end



  FUNCTION dist_fn_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    CHARACTER(LEN=string_length) :: part1
    INTEGER :: part2, ispecies
    INTEGER :: work, io
    REAL(num) :: work1, work2

    errcode = c_err_none
    IF (deck_state .EQ. c_ds_first) RETURN
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element, 'name')) THEN
      working_block%name = value
      RETURN
    ENDIF

    IF (str_cmp(element, 'ndims')) THEN
      work = as_integer(value, errcode)
      IF (work .GE. 1 .AND. work .LE. 3) THEN
        working_block%ndims = work
      ELSE
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Distribution functions can only be 1D, 2D or 3D'
          ENDDO
        ENDIF
        errcode = c_err_bad_value
      ENDIF
      RETURN
    ENDIF

    IF (working_block%ndims .EQ. -1) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must set number of dimensions before setting other', &
              ' distribution'
          WRITE(io,*) 'function properties.'
        ENDDO
      ENDIF
      extended_error_string = 'ndims'
      errcode = c_err_required_element_not_set
      RETURN
    ENDIF

    IF (str_cmp(element, 'dumpmask')) THEN
      working_block%dumpmask = as_integer(value, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'restrict_x')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(1) = .TRUE.
      working_block%restrictions(:,1) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_y')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(2) = .TRUE.
      working_block%restrictions(:,2) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_px')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(3) = .TRUE.
      working_block%restrictions(:,3) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_py')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(4) = .TRUE.
      working_block%restrictions(:,4) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'restrict_pz')) THEN
      CALL split_range(value, work1, work2, errcode)
      IF (errcode .NE. c_err_none) RETURN
      working_block%use_restrictions(5) = .TRUE.
      working_block%restrictions(:,5) = (/work1, work2/)
    ENDIF

    IF (str_cmp(element, 'include_species')) THEN
      ispecies = as_integer(value, errcode)
      IF (errcode .EQ. c_err_none) THEN
        IF (ispecies .GT. 0 .AND. ispecies .LE. n_species) THEN
          working_block%use_species(ispecies) = .TRUE.
        ELSE
          IF (rank .EQ. 0) THEN
            DO io = stdout, du, du - stdout ! Print to stdout and to file
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
