MODULE deck_control_block

  USE mpi
  USE strings_advanced
  USE fields

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: control_deck_initialise, control_deck_finalise
  PUBLIC :: control_block_start, control_block_end
  PUBLIC :: control_block_handle_element, control_block_check

  INTEGER, PARAMETER :: control_block_elements = 14 + 4 * c_ndims
  LOGICAL, DIMENSION(control_block_elements) :: control_block_done
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      control_block_name = (/ &
          "nx                ", &
          "ny                ", &
          "x_min             ", &
          "x_max             ", &
          "y_min             ", &
          "y_max             ", &
          "nprocx            ", &
          "nprocy            ", &
          "npart             ", &
          "nsteps            ", &
          "t_end             ", &
          "dt_multiplier     ", &
          "dlb_threshold     ", &
          "icfile            ", &
          "restart_snapshot  ", &
          "neutral_background", &
          "field_order       ", &
          "stdout_frequency  ", &
          "use_random_seed   ", &
          "smooth_currents   ", &
          "coulomb_log       ", &
          "collide           " /)
  CHARACTER(LEN=string_length), DIMENSION(control_block_elements) :: &
      alternate_name = (/ &
          "nx                ", &
          "ny                ", &
          "x_start           ", &
          "x_end             ", &
          "y_start           ", &
          "y_end             ", &
          "nprocx            ", &
          "nprocy            ", &
          "npart             ", &
          "nsteps            ", &
          "t_end             ", &
          "dt_multiplier     ", &
          "dlb_threshold     ", &
          "icfile            ", &
          "restart_snapshot  ", &
          "neutral_background", &
          "field_order       ", &
          "stdout_frequency  ", &
          "use_random_seed   ", &
          "smooth_currents   ", &
          "coulomb_log       ", &
          "collide           " /)

CONTAINS

  SUBROUTINE control_deck_initialise

    IF (deck_state .NE. c_ds_first) RETURN
    control_block_done = .FALSE.

  END SUBROUTINE control_deck_initialise



  SUBROUTINE control_deck_finalise

  END SUBROUTINE control_deck_finalise



  SUBROUTINE control_block_start

  END SUBROUTINE control_block_start



  SUBROUTINE control_block_end

  END SUBROUTINE control_block_end



  FUNCTION control_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: loop, elementselected, field_order, ierr, io

    errcode = c_err_none


    errcode = c_err_unknown_element

    elementselected = 0

    DO loop = 1, control_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(control_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF ((deck_state .NE. c_ds_first) .AND. &
        (elementselected .LE. (4*c_ndims+13))) THEN
      errcode = c_err_none
      RETURN
    ENDIF

    IF (elementselected .EQ. 0) RETURN
    IF (control_block_done(elementselected).AND. &
        (elementselected .LE. (4*c_ndims+13))) THEN
      errcode = c_err_preset_element
      RETURN
    ENDIF
    control_block_done(elementselected) = .TRUE.
    errcode = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      nx_global = as_integer(value, errcode)
    CASE(2)
      ny_global = as_integer(value, errcode)
    CASE(c_ndims+1)
      x_min = as_real(value, errcode)
    CASE(c_ndims+2)
      x_max = as_real(value, errcode)
    CASE(c_ndims+3)
      y_min = as_real(value, errcode)
    CASE(c_ndims+4)
      y_max = as_real(value, errcode)
    CASE(3*c_ndims+1)
      nprocx = as_integer(value, errcode)
    CASE(3*c_ndims+2)
      nprocy = as_integer(value, errcode)
    CASE(4*c_ndims+1)
      npart_global = as_long_integer(value, errcode)
    CASE(4*c_ndims+2)
      nsteps = as_integer(value, errcode)
    CASE(4*c_ndims+3)
      t_end = as_real(value, errcode)
    CASE(4*c_ndims+4)
      dt_multiplier = as_real(value, errcode)
    CASE(4*c_ndims+5)
      dlb_threshold = as_real(value, errcode)
      dlb = .TRUE.
    CASE(4*c_ndims+6)
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The "icfile" option is no longer supported.'
          WRITE(io,*) 'Please use the "import" directive instead'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    CASE(4*c_ndims+7)
      restart_snapshot = as_integer(value, errcode)
      ic_from_restart = .TRUE.
    CASE(4*c_ndims+8)
      neutral_background = as_logical(value, errcode)
    CASE(4*c_ndims+9)
      field_order = as_integer(value, errcode)
      IF (field_order .NE. 2 .AND. field_order .NE. 4 &
          .AND. field_order .NE. 6) THEN
        errcode = c_err_bad_value
      ELSE
        CALL set_field_order(field_order)
      ENDIF
    CASE(4*c_ndims+10)
      stdout_frequency = as_integer(value, errcode)
    CASE(4*c_ndims+11)
      use_random_seed = as_logical(value, errcode)
    CASE(4*c_ndims+12)
      smooth_currents = as_logical(value, errcode)
    CASE(4*c_ndims+13)
      CALL set_coulomb_log(value, errcode)
    CASE(4*c_ndims+14)
      CALL set_collision_matrix(TRIM(ADJUSTL(value)), errcode)
    END SELECT

  END FUNCTION control_block_handle_element



  FUNCTION control_block_check() RESULT(errcode)

    INTEGER :: errcode, index, io

    errcode = c_err_none

    ! nprocx/y and npart are optional
    control_block_done(3*c_ndims+1:4*c_ndims+1) = .TRUE.

    ! Only one of nsteps or t_end need be specified
    IF (control_block_done(4*c_ndims+2)) &
        control_block_done(4*c_ndims+3) = .TRUE.
    IF (control_block_done(4*c_ndims+3)) &
        control_block_done(4*c_ndims+2) = .TRUE.

    ! All entries after t_end are optional
    control_block_done(4*c_ndims+4:) = .TRUE.

    DO index = 1, control_block_elements
      IF (.NOT. control_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required control block element ', &
                TRIM(ADJUSTL(control_block_name(index))), &
                ' absent. Please create this entry in the input deck'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
    ENDDO

    IF (.NOT. neutral_background) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'The option "neutral_background=F" is not supported', &
              ' in this version of EPOCH.'
        ENDDO
      ENDIF
      errcode = c_err_terminate
    ENDIF

  END FUNCTION control_block_check


! #ifdef COLLISIONS
! The following code is all about reading the coll_pairs from the input deck

  FUNCTION species_index_from_name(str_in, err)

    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: species_index_from_name
    INTEGER :: ispecies

    DO ispecies=1,n_species
      IF (str_cmp(TRIM(str_in),TRIM(species_list(ispecies)%name))) THEN
        species_index_from_name = ispecies
        RETURN
      ENDIF
    ENDDO
    species_index_from_name = -1
    err = c_err_bad_value
  END FUNCTION species_index_from_name


  SUBROUTINE get_token(str_in, str_out, token_out, err)

    CHARACTER(*), INTENT(IN) :: str_in
    CHARACTER(*), INTENT(OUT) :: str_out
    CHARACTER(*), INTENT(OUT) :: token_out

    INTEGER, INTENT(INOUT) :: err
    INTEGER :: str_len, char, pos
    CHARACTER(1) :: c

    str_len = LEN(str_in)
    pos = str_len

    DO char = 1, str_len
      c = str_in(char:char)
      IF (c .EQ. ' ' .OR. ICHAR(c) .EQ. 32)  THEN
        pos = char
        EXIT
      ENDIF
    ENDDO

    IF (pos .LT. str_len) THEN
      str_out = TRIM(ADJUSTL(str_in(pos+1:str_len)))
    ELSE
      str_out = "";
    ENDIF

    token_out = TRIM(str_in(1:pos))

  END SUBROUTINE get_token



  SUBROUTINE duplicate_collmatrix_warn(species1, species2)
    CHARACTER(*), INTENT(IN) :: species1
    CHARACTER(*), INTENT(IN) :: species2
    INTEGER :: io
    DO io = stdout, du, du - stdout ! Print to stdout and to file
      WRITE(io,*)
      WRITE(io,*) '*** WARNING ***'
      WRITE(io,*) 'The collide parameter for ' &
        // TRIM(species1) // " <-> " // TRIM(species2)
      WRITE(io,*) 'has been set multiple times!'
      WRITE(io,*) 'Collisions will only be carried out once per species pair.'
      WRITE(io,*) 'Later specifications will always override earlier ones.'
      WRITE(io,*)
    ENDDO
  END SUBROUTINE duplicate_collmatrix_warn


  SUBROUTINE set_collision_matrix(str_in, errcode)
    CHARACTER(*), INTENT(IN) :: str_in
    INTEGER, INTENT(INOUT) :: errcode
    CHARACTER(LEN=string_length) :: token, tstr1, tstr2
    CHARACTER(LEN=string_length) :: species1, species2
    REAL(num) :: collstate

    INTEGER sp1, sp2
    IF (deck_state .NE. c_ds_last) RETURN

    CALL get_token(str_in, tstr1, token, errcode)

    IF (errcode .NE. 0) RETURN

    IF (str_cmp(TRIM(token), "all")) THEN
      coll_pairs = 1.0_num
      RETURN
    ENDIF

    IF (str_cmp(TRIM(token), "none")) THEN
      coll_pairs = -1.0_num
      RETURN
    ENDIF

    sp1 = species_index_from_name(token,errcode)
    IF (errcode .NE. 0) RETURN
    species1 = token

    CALL get_token(tstr1, tstr2, species2, errcode)
    IF (errcode .NE. 0) RETURN

    sp2 = species_index_from_name(species2,errcode)
    IF (errcode .NE. 0) RETURN

    collstate = 1.0_num
    IF (str_cmp(TRIM(tstr2), "on") .OR. str_cmp(TRIM(tstr2), "")) THEN
      collstate = 1.0_num
    ELSEIF (str_cmp(TRIM(tstr2), "off")) THEN
      collstate = -1.0_num
    ELSE
      collstate = as_real(tstr2, errcode)
      IF (errcode .NE. 0) RETURN
    ENDIF

    coll_pairs(sp1, sp2) = collstate
    coll_pairs(sp2, sp1) = collstate
    IF (coll_pairs_touched(sp1, sp2)) THEN
      CALL duplicate_collmatrix_warn(species1, species2)
    ENDIF
    coll_pairs_touched(sp1, sp2) = .TRUE.
    coll_pairs_touched(sp2, sp1) = .TRUE.
  END SUBROUTINE set_collision_matrix

  SUBROUTINE set_coulomb_log(value, errcode)
    CHARACTER(*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: errcode

    IF (str_cmp(value, "auto")) THEN
      coulomb_log_auto = .TRUE.
      RETURN
    ENDIF

    coulomb_log_auto = .FALSE.
    coulomb_log = as_real(value, errcode)

  END SUBROUTINE set_coulomb_log
! #endif

END MODULE deck_control_block
