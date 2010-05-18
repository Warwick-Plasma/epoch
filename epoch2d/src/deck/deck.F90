MODULE deck

  ! Basic operations
  USE strings
  ! Deck internals
  USE deck_constant_block
  USE deck_deo_block
  ! Deck Blocks
  USE deck_control_block
  USE deck_boundaries_block
  USE deck_species_block
  USE deck_io_block
  USE deck_window_block
  ! Initial Condition Blocks
  USE deck_ic_laser_block
  USE deck_ic_fields_block
  ! Extended IO Blocks
  USE deck_eio_dist_fn_block
#ifdef PARTICLE_PROBES
  USE deck_eio_particle_probe_block
#endif
  ! Custom blocks
  USE custom_deck
  USE version_data

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_deck, write_input_decks

  SAVE
  CHARACTER(LEN=string_length) :: current_block_name
  LOGICAL :: invalid_block

  INTEGER, PARAMETER :: buffer_size = 1024
  TYPE :: file_buffer
    CHARACTER(LEN=45+data_dir_max_length) :: filename
    CHARACTER(LEN=buffer_size), DIMENSION(:), POINTER :: buffer
    INTEGER :: pos, idx, length
    TYPE(file_buffer), POINTER :: next
  END TYPE file_buffer

  TYPE(file_buffer), POINTER :: file_buffer_head
  INTEGER :: nbuffers = 0

CONTAINS

  !----------------------------------------------------------------------------
  ! These subroutines actually call the routines which read the deck blocks
  !----------------------------------------------------------------------------

  ! This subroutine is called when a new block is started
  ! If a block NEEDS to do something when it starts, then
  ! The revelevant subroutine should be called here
  SUBROUTINE start_block(block_name)

    CHARACTER(LEN=*), INTENT(IN) :: block_name

    IF (str_cmp(block_name, "laser")) THEN
      IF (deck_state .EQ. c_ds_ic) CALL laser_start
    ELSE IF (str_cmp(block_name, "window")) THEN
      IF (deck_state .EQ. c_ds_deck) CALL window_start
    ELSE IF (str_cmp(block_name, "dist_fn")) THEN
      IF (deck_state .EQ. c_ds_eio) CALL dist_fn_start
#ifdef PARTICLE_PROBES
    ELSE IF (str_cmp(block_name, "probe")) THEN
      IF (deck_state .EQ. c_ds_eio) CALL probe_block_start
#endif
    ELSE IF (str_cmp(block_name, "fields")) THEN
      IF (deck_state .EQ. c_ds_ic) CALL fields_start
    ELSE IF (str_cmp(block_name, "species")) THEN
      CALL species_start
    ENDIF

  END SUBROUTINE start_block



  ! This subroutine is called when a new block is ended
  ! If a block NEEDS to do something when it ends, then
  ! The revelevant subroutine should be called here
  SUBROUTINE end_block(block_name)

    CHARACTER(LEN=*), INTENT(IN) :: block_name

    IF (str_cmp(block_name, "laser")) THEN
      IF (deck_state .EQ. c_ds_ic) CALL laser_end
    ELSE IF (str_cmp(block_name, "dist_fn")) THEN
      IF (deck_state .EQ. c_ds_eio) CALL dist_fn_end
    ELSE IF (str_cmp(block_name, "species")) THEN
      IF (deck_state .EQ. c_ds_deck) CALL species_end
#ifdef PARTICLE_PROBES
    ELSE IF (str_cmp(block_name, "probe")) THEN
      IF (deck_state .EQ. c_ds_eio) CALL probe_block_end
#endif
    ENDIF

  END SUBROUTINE end_block



  FUNCTION handle_block(block_name, block_element, block_value)

    CHARACTER(LEN=*), INTENT(IN) :: block_name, block_element, block_value
    CHARACTER(LEN=string_length) :: part1
    INTEGER :: handle_block
    INTEGER :: part2, val

    handle_block = c_err_none
    ! Constants can be defined in any deck state, so put them here
    IF (str_cmp(block_name, "constant")) THEN
      handle_block = handle_constant_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "deo")) THEN
      handle_block = handle_deo_deck(block_element, block_value)
      RETURN
    ENDIF

    ! Test for known blocks
    IF (str_cmp(block_name, "control"))  THEN
      IF (deck_state .EQ. c_ds_deck) &
          handle_block = handle_control_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "boundaries")) THEN
      IF (deck_state .EQ. c_ds_deck) &
          handle_block = handle_boundary_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "species")) THEN
      handle_block = handle_species_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "output")) THEN
      IF (deck_state .EQ. c_ds_deck) &
          handle_block = handle_io_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "window")) THEN
      IF (deck_state .EQ. c_ds_deck) &
          handle_block = handle_window_deck(block_element, block_value)
      RETURN
    ENDIF

    ! Initial conditions blocks go here
    IF (str_cmp(block_name, "fields")) THEN
      IF (deck_state .EQ. c_ds_ic) &
          handle_block = handle_ic_fields_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "laser")) THEN
      IF (deck_state .EQ. c_ds_ic) &
          handle_block = handle_ic_laser_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "dist_fn")) THEN
      IF (deck_state .EQ. c_ds_eio) &
          handle_block = handle_eio_dist_fn_deck(block_element, block_value)
      RETURN
    ENDIF

    IF (str_cmp(block_name, "probe")) THEN
#ifdef PARTICLE_PROBES
      IF (deck_state .EQ. c_ds_eio) &
          handle_block = handle_probe_deck(block_element, block_value)
      RETURN
#else
      IF (deck_state .EQ. c_ds_eio) THEN
        handle_block = IOR(handle_block, c_err_pp_options_wrong)
        extended_error_string = "-DPARTICLE_PROBES"
      ENDIF
      RETURN
#endif
    ENDIF

    handle_block = c_err_unknown_block
    ! Pass through to the custom block
    handle_block = handle_custom_block(block_name, block_element, block_value)

  END FUNCTION handle_block



  ! These subroutines are there to check for the basic minimal compulsory
  ! blocks are present. They're a bit ugly, but they seem to be the easiest
  ! way to do it without adding complexity to the code
  SUBROUTINE check_compulsory_blocks(errcode_deck)

    LOGICAL :: problem_found
    INTEGER, INTENT(INOUT) :: errcode_deck

    problem_found = .FALSE.

    errcode_deck = c_err_none

    IF (deck_state .EQ. c_ds_deck) THEN
      errcode_deck = IOR(errcode_deck, check_control_block())
      errcode_deck = IOR(errcode_deck, check_boundary_block())
      errcode_deck = IOR(errcode_deck, check_species_block())
      errcode_deck = IOR(errcode_deck, check_io_block())
      errcode_deck = IOR(errcode_deck, check_window_block())
      errcode_deck = IOR(errcode_deck, check_custom_blocks())
    ELSE IF (deck_state .EQ. c_ds_ic) THEN
      errcode_deck = IOR(errcode_deck, check_ic_fields_block())
      errcode_deck = IOR(errcode_deck, check_species_block())
    ENDIF
    errcode_deck = IOR(errcode_deck, check_custom_blocks())

    problem_found = (IAND(errcode_deck, c_err_missing_elements) .NE. 0)

    IF (problem_found) THEN
      errcode_deck = IOR(errcode_deck, c_err_terminate)
      IF (rank .EQ. 0) THEN
        WRITE(*, *)
        WRITE(*, *) '***ERROR***'
        WRITE(*, *) 'Not all required elements of input deck specified.'
        WRITE(*, *) 'Please fix input deck and rerun code'
        WRITE(40,*)
        WRITE(40,*) '***ERROR***'
        WRITE(40,*) 'Not all required elements of input deck specified.'
        WRITE(40,*) 'Please fix input deck and rerun code'
      ENDIF
    ELSE
      IF (rank .EQ. 0) THEN
        IF (deck_state .EQ. c_ds_deck) THEN
          WRITE(*, *) 'Input deck complete and valid. Attempting to set up ' &
              // 'equilibrium'
          WRITE(*, *)
          WRITE(40,*) 'Input deck complete and valid.'
        ELSE IF (deck_state .EQ. c_ds_ic) THEN
          WRITE(*, *) 'Initial conditions complete and valid. Attempting to ' &
              // 'load particles'
          WRITE(*, *)
          WRITE(40,*) 'Initial conditions complete and valid.'
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE check_compulsory_blocks



  !----------------------------------------------------------------------------
  ! These subroutines are the in depth detail of how the parser works
  !----------------------------------------------------------------------------

  FUNCTION get_free_lun()

    ! This subroutine simply cycles round until it finds a free lun between
    ! min_lun and max_lun
    INTEGER :: get_free_lun
    INTEGER :: lun, ierr
    INTEGER, PARAMETER :: min_lun = 10, max_lun = 20
    LOGICAL :: is_open

    is_open = .TRUE.

    lun = min_lun
    DO
      INQUIRE(unit=lun, opened=is_open)
      IF (.NOT. is_open) EXIT
      lun = lun+1
      IF (lun .GT. max_lun) THEN
        IF (rank .EQ. 0) THEN
          WRITE(*,*) '***FATAL ERROR*** unable to open lun for input deck read'
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF
    ENDDO

    get_free_lun = lun

  END FUNCTION get_free_lun



  RECURSIVE SUBROUTINE read_deck(filename, first_call)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: first_call
    CHARACTER :: u1
    INTEGER :: pos = 1, flip = 1, s, f, elements = 0, lun
    LOGICAL :: is_comment
    TYPE(string_type), DIMENSION(2) :: deck_values
    CHARACTER(LEN=45+data_dir_max_length) :: deck_filename, status_filename
    LOGICAL :: terminate = .FALSE., exists
    INTEGER :: errcode_deck, ierr, i
    LOGICAL :: white_space_over
    CHARACTER(LEN=buffer_size), DIMENSION(:), ALLOCATABLE :: tmp_buffer
    TYPE(file_buffer), POINTER :: fbuf
    LOGICAL :: already_parsed, got_eor

    ! No error yet
    errcode_deck = c_err_none
    ! Characteristic string which represents a "blank" string
    blank = "BLANKBLANK"

    lun = 5
    white_space_over = .FALSE.
    already_parsed = .FALSE.

    ! Make the whole filename by adding the data_dir to the filename
    deck_filename = TRIM(ADJUSTL(data_dir))// '/' // TRIM(ADJUSTL(filename))

    ! deck_state tells the code whether it's parsing the normal input deck
    ! Or the initial conditions. You can add more states if you want.
    ! Just search for deck_state
    IF (deck_state .EQ. c_ds_deck) THEN
      status_filename = TRIM(ADJUSTL(data_dir))// '/' // "deck.status"
    ELSE IF (deck_state .EQ. c_ds_ic) THEN
      status_filename = TRIM(ADJUSTL(data_dir))// '/' // "ic.status"
    ELSE IF (deck_state .EQ. c_ds_eio) THEN
      status_filename = TRIM(ADJUSTL(data_dir))// '/' // "eio.status"
    ENDIF

    ! If this is the first time that this deck has been called then do some
    ! housekeeping. Put any initialisation code that is needed in here
    IF (first_call) THEN
      control_block_done = .FALSE.
      boundary_block_done = .FALSE.
      CALL species_initialise
    ENDIF

    ! Is comment is a flag which tells the code when a # character has been
    ! found and everything beyond it is a comment
    is_comment = .FALSE.

    ! rank 0 reads the file and then passes it out to the other nodes using
    ! MPI_BCAST
    IF (rank .EQ. 0) THEN
      IF (.NOT. ASSOCIATED(file_buffer_head)) THEN
        ALLOCATE(file_buffer_head)
        fbuf=>file_buffer_head
        fbuf%filename = ""
        NULLIFY(fbuf%next)
      ELSE
        fbuf=>file_buffer_head
      ENDIF

      DO WHILE (ASSOCIATED(fbuf%next))
        fbuf=>fbuf%next
        IF (fbuf%filename .EQ. deck_filename) THEN
          already_parsed = .TRUE.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT. already_parsed) THEN
        ALLOCATE(fbuf%next)
        nbuffers = nbuffers + 1
        fbuf=>fbuf%next
        fbuf%filename = deck_filename
        fbuf%pos = 1
        fbuf%idx = 1
        fbuf%length = 1
        NULLIFY(fbuf%next)
        ALLOCATE(fbuf%buffer(fbuf%length))
      ENDIF

      ! Check whether or not the input deck file requested exists
      INQUIRE(file=deck_filename, exist=exists)
      IF (.NOT. exists) THEN
        PRINT *, '***ERROR***'
        PRINT *, 'Input deck file "' // deck_filename // '" does not exist.'
        PRINT *, 'Create the file and rerun the code.'
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF

      ! Get a free lun. Don't use a constant lun to allow for recursion
      lun = get_free_lun()
      OPEN(unit=lun, file=TRIM(ADJUSTL(deck_filename)))
      IF (first_call) THEN
        OPEN(unit=40, file=status_filename)
        WRITE(40,*) ascii_header
        WRITE(40,*)
      ENDIF
      deck_values(1)%value = ""
      deck_values(2)%value = ""

      ! Use non-advancing IO to pop characters off the deck file one at a time
      ! Use basic token parsing to split into two substrings across
      ! an "=" or ":" symbol
      DO
        errcode_deck = c_err_none
        ! Read a character
        ! For ordinary characters, f is zero
        ! When an end of line character is read then got_eor is .TRUE.
        ! When end of file is reached, f is negative and got_eor is .FALSE.
        got_eor = .TRUE.
        READ(lun, "(A1)", advance='no', size=s, iostat=f, eor=10) u1
        got_eor = .FALSE.

10      IF (.NOT. already_parsed) THEN
          ! Store character in a buffer so that we can write the input deck
          ! contents to a restart dump
          IF (f .EQ. 0) THEN
            fbuf%buffer(fbuf%idx)(fbuf%pos:fbuf%pos) = u1
          ELSE IF (got_eor) THEN
            fbuf%buffer(fbuf%idx)(fbuf%pos:fbuf%pos) = ACHAR(10)
          ELSE
            fbuf%buffer(fbuf%idx)(fbuf%pos:fbuf%pos) = ACHAR(0)
            fbuf%pos = fbuf%pos - 1
          ENDIF

          ! If we reached the end of the character string then move to the next
          ! element of the array
          IF (fbuf%pos .EQ. buffer_size) THEN
            ! If we reached the end of the array then allocate some more
            IF (fbuf%idx .EQ. fbuf%length) THEN
              ALLOCATE(tmp_buffer(fbuf%length))
              DO i = 1,fbuf%length
                tmp_buffer(i) = fbuf%buffer(i)
              ENDDO
              DEALLOCATE(fbuf%buffer)
              ALLOCATE(fbuf%buffer(2*fbuf%length))
              DO i = 1,fbuf%length
                fbuf%buffer(i) = tmp_buffer(i)
              ENDDO
              DEALLOCATE(tmp_buffer)
              fbuf%length = 2*fbuf%length
            ENDIF
            fbuf%pos = 1
            fbuf%idx = fbuf%idx + 1
          ELSE
            fbuf%pos = fbuf%pos + 1
          ENDIF
        ENDIF

        ! If the character is a # then switch to comment mode
        IF (u1 .EQ. '#') is_comment = .TRUE.

        ! If not in comment mode then use the character
        IF (.NOT. is_comment) THEN
          ! If the current character isn't a special character then just stick
          ! it in the buffer
          IF (u1 .NE. '=' .AND. u1 .NE. CHAR(9) .AND. u1 .NE. ':' &
              .AND. f .EQ. 0) THEN
            IF ((u1 .NE. ' ' .AND. u1 .NE. CHAR(32)) .OR. white_space_over) THEN
              deck_values(flip)%value(pos:pos) = u1
              pos = pos+1
              white_space_over = .TRUE.
            ENDIF
          ENDIF

          ! If it's equals or : then you're parsing the other part of the
          ! expression
          IF (u1 .EQ. '=' .OR. u1 .EQ. ':') THEN
            flip = 2
            pos = 1
          ENDIF
        ENDIF

        ! If got_eor is .TRUE. then you've reached the end of the line, so
        ! comment state is definitely false
        IF (got_eor) is_comment = .FALSE.

        ! If you've not read a blank line then
        IF (got_eor .AND. pos .GT. 1) THEN
          elements = elements+1
          flip = 1
          pos = 1
          deck_values(1)%value = TRIM(ADJUSTL(deck_values(1)%value))
          deck_values(2)%value = TRIM(ADJUSTL(deck_values(2)%value))
          CALL MPI_BCAST(1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)
          CALL MPI_BCAST(nbuffers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)
          CALL MPI_BCAST(deck_values(1)%value, string_length, MPI_CHARACTER, &
              0, MPI_COMM_WORLD, errcode)
          CALL MPI_BCAST(deck_values(2)%value, string_length, MPI_CHARACTER, &
              0, MPI_COMM_WORLD, errcode)
          CALL handle_deck_element(deck_values(1)%value, deck_values(2)%value, &
              errcode_deck)
          deck_values(1)%value = ""
          deck_values(2)%value = ""
          is_comment = .FALSE.
          white_space_over = .FALSE.
        ENDIF
        IF (f .LT. 0 .AND. .NOT.got_eor) THEN
          CALL MPI_BCAST(0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)
          CLOSE(lun)
          EXIT
        ENDIF
        terminate = terminate .OR. IAND(errcode_deck, c_err_terminate) .NE. 0
        IF (terminate) EXIT
      ENDDO
    ELSE
      DO
        errcode_deck = c_err_none
        CALL MPI_BCAST(f, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)
        IF (f .EQ. 0) EXIT
          CALL MPI_BCAST(nbuffers, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)
        CALL MPI_BCAST(deck_values(1)%value, string_length, MPI_CHARACTER, &
              0, MPI_COMM_WORLD, errcode)
        CALL MPI_BCAST(deck_values(2)%value, string_length, MPI_CHARACTER, &
              0, MPI_COMM_WORLD, errcode)
        CALL handle_deck_element(deck_values(1)%value, deck_values(2)%value, &
              errcode_deck)
        deck_values(1)%value = ""
        deck_values(2)%value = ""
        terminate = terminate .OR. IAND(errcode_deck, c_err_terminate) .NE. 0
        IF (terminate) EXIT
      ENDDO
    ENDIF

    CALL MPI_BARRIER(MPI_COMM_WORLD, errcode)
!!$    IF (.NOT. first_call)  THEN
!!$       RETURN
!!$    ENDIF

    ! Don't check compulsory blocks if going to bomb anyway, just stinks up
    ! the output file
    IF (.NOT. terminate .AND. first_call) THEN
      CALL check_compulsory_blocks(errcode_deck)
      CALL species_finalise
    ENDIF
    terminate = terminate .OR. IAND(errcode_deck, c_err_terminate) .NE. 0
    ! Fatal error, cause code to bomb
    IF (terminate .AND. rank .EQ. 0) THEN
      WRITE(*, *)
      WRITE(*, *) '***FATAL ERROR***'
      WRITE(*, *) 'The code cannot parse the input deck sufficiently to run.'
      WRITE(*, *) 'Please read the output file "deck.status" to check for ' &
          // 'errors.'
      WRITE(40,*)
      WRITE(40,*) '***FATAL ERROR***'
      WRITE(40,*) 'The code cannot parse the input deck sufficiently to run.'
      WRITE(40,*) 'Please read this file and correct any errors mentioned.'
      PRINT *
      PRINT *
      PRINT *
      PRINT *
      PRINT *
      PRINT *
    ENDIF

    IF (first_call) CLOSE(40)

    IF (terminate) CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)

    CALL MPI_BARRIER(MPI_COMM_WORLD, errcode)

  END SUBROUTINE read_deck



  SUBROUTINE handle_deck_element(element, value, errcode_deck)

    CHARACTER(*), INTENT(IN) :: element
    CHARACTER(*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: errcode_deck
    INTEGER :: rank_check
    INTEGER, SAVE :: err_count

    rank_check = 0

    IF (str_cmp(element, "import")) THEN
      invalid_block = .TRUE.
      IF (rank .EQ. rank_check) THEN
        WRITE(40,*)
        WRITE(40,*) 'Importing "' // TRIM(ADJUSTL(value)) // '" file'
        WRITE(40,*)
      ENDIF
      CALL read_deck(TRIM(ADJUSTL(value)), .FALSE.)
      RETURN
    ENDIF

    IF (str_cmp(element, "begin")) THEN
      errcode_deck = handle_block(value, blank, blank)
      invalid_block = IAND(errcode_deck, c_err_unknown_block) .NE. 0
      invalid_block = invalid_block .OR. IAND(errcode_deck, &
          c_err_pp_options_wrong) .NE. 0
      IF (invalid_block .AND. rank .EQ. rank_check) THEN
        IF(IAND(errcode_deck, c_err_pp_options_wrong) .NE. 0) THEN
          WRITE(*, *)
          WRITE(*, *) '***WARNING***'
          WRITE(*, *) 'The block "' // TRIM(value) // '" cannot be set because'
          WRITE(*, *) 'the code has not been compiled with the correct ' &
              // 'preprocessor options.'
          WRITE(*, *) 'Code will continue, but to use selected features, ' &
              // 'please recompile with the'
          WRITE(*, *) TRIM(extended_error_string) // ' option'
          WRITE(*, *)
          WRITE(40,*)
          WRITE(40,*) '***WARNING***'
          WRITE(40,*) 'The block "' // TRIM(value) // '" cannot be set because'
          WRITE(40,*) 'the code has not been compiled with the correct ' &
              // 'preprocessor options.'
          WRITE(40,*) 'Code will continue, but to use selected features, ' &
              // 'please recompile with the'
          WRITE(40,*) TRIM(extended_error_string) // ' option'
          WRITE(40,*)
        ELSE
          WRITE(*, *) '***WARNING***'
          WRITE(*, *) 'Unknown block "' // TRIM(value) // '" in input deck, ' &
              // 'ignoring', deck_state
        ENDIF
      ENDIF
      CALL start_block(value)
      err_count = 0
      current_block_name = value
      IF (rank .EQ. rank_check) THEN
        WRITE(40,*) 'Beginning "' // TRIM(ADJUSTL(value)) // '" block'
        WRITE(40,*)
      ENDIF
      ! Reset errcode_deck here because reporting c_err_unknown_element is OK
      errcode_deck = c_err_none
      RETURN
    ENDIF
    IF (str_cmp(element, "end")) THEN
      CALL end_block(current_block_name)
      invalid_block = .TRUE.
      IF (rank .EQ. rank_check) THEN
        WRITE(40,*)
        WRITE(40,*) 'Ending "' // TRIM(ADJUSTL(value)) // '" block'
        WRITE(40,*)
        IF (err_count .NE. 0) THEN
          WRITE(40,*) '***WARNING***'
          WRITE(40,*) 'block "' // TRIM(ADJUSTL(value)) // '" contains errors'
          WRITE(40,*)
        ENDIF
      ENDIF
      RETURN
    ENDIF

    ! Check invalid block to limit amount of rubbish that appears
    ! If the input deck is invalid
    IF (.NOT. invalid_block) THEN
      errcode_deck = handle_block(current_block_name, element, value)
    ELSE
      RETURN
    ENDIF

    IF (errcode_deck .EQ. c_err_none) THEN
      IF (rank .EQ. rank_check) &
          WRITE(40, *) CHAR(9), "Element ", TRIM(ADJUSTL(element)), "=", &
              TRIM(ADJUSTL(value)), " handled OK"
      RETURN
    ENDIF
    ! Test for error conditions
    ! If an error is fatal then set terminate to .TRUE.
    IF (IAND(errcode_deck, c_err_unknown_element) .NE. 0) THEN
      IF (rank .EQ. rank_check) THEN
        WRITE(*, *)
        WRITE(*, *) '***WARNING***'
        WRITE(*, *) 'Unrecognised element "' // TRIM(element) &
            // '" in input deck.'
        WRITE(*, *) 'Code will continue to run, but behaviour is undefined'
        WRITE(*, *)
        WRITE(40,*)
        WRITE(40,*) '***WARNING***'
        WRITE(40,*) 'Unrecognised element "' // TRIM(element) &
            // '" in input deck.'
        WRITE(40,*) 'Code will continue to run, but behaviour is undefined'
        WRITE(40,*)
      ENDIF
    ENDIF
    IF (IAND(errcode_deck, c_err_preset_element) .NE. 0) THEN
      IF (rank .EQ. rank_check) THEN
        WRITE(*, *)
        WRITE(*, *) '***WARNING***'
        WRITE(*, *) 'element "' // TRIM(element) &
            // '" is set multiple times in this deck.'
        WRITE(*, *) 'Code will continue using first value in deck'
        WRITE(*, *)
        WRITE(40,*)
        WRITE(40,*) '***WARNING***'
        WRITE(40,*) 'element "' // TRIM(element) &
            // '" is set multiple times in this deck.'
        WRITE(40,*) 'Code will continue using first value in deck'
        WRITE(40,*)
      ENDIF
    ENDIF
    IF (IAND(errcode_deck, c_err_preset_element_use_later) .NE. 0) THEN
      IF (rank .EQ. rank_check) THEN
        WRITE(*, *)
        WRITE(*, *) '***WARNING***'
        WRITE(*, *) 'element "' // TRIM(element) &
            // '" is set multiple times in this deck.'
        WRITE(*, *) 'Code will continue using last value in deck'
        WRITE(*, *)
        WRITE(40,*)
        WRITE(40,*) '***WARNING***'
        WRITE(40,*) 'element "' // TRIM(element) &
            // '" is set multiple times in this deck.'
        WRITE(40,*) 'Code will continue using last value in deck'
        WRITE(40,*)
      ENDIF
    ENDIF
    IF (IAND(errcode_deck, c_err_bad_value) .NE. 0) THEN
      IF (rank .EQ. rank_check) THEN
        WRITE(*, *)
        WRITE(*, *) '***ERROR***'
        WRITE(*, *) 'value "' // TRIM(value) // '" in element "' &
            // TRIM(element) // '" is'
        WRITE(*, *) 'invalid or could not be parsed. Code will terminate.'
        WRITE(*, *)
        WRITE(40,*)
        WRITE(40,*) '***ERROR***'
        WRITE(40,*) 'value "' // TRIM(value) // '" in element "' &
            // TRIM(element) // '" is'
        WRITE(40,*) 'invalid or could not be parsed. Code will terminate.'
        WRITE(40,*)
        errcode_deck = IOR(errcode_deck, c_err_terminate)
      ENDIF
    ENDIF

    IF (IAND(errcode_deck, c_err_required_element_not_set) .NE. 0) THEN
      IF (rank .EQ. rank_check) THEN
        WRITE(*, *)
        WRITE(*, *) '***ERROR***'
        WRITE(*, *) 'value "' // TRIM(value) // '" in element "' &
            // TRIM(element) // '" cannot be'
        WRITE(*, *) 'set because a prerequisite element "' &
            // TRIM(extended_error_string) // '" has'
        WRITE(*, *) 'not been set. Code will terminate'
        WRITE(*, *)
        WRITE(40,*)
        WRITE(40,*) '***ERROR***'
        WRITE(40,*) 'value "' // TRIM(value) // '" in element "' &
            // TRIM(element) // '" cannot be'
        WRITE(40,*) 'set because a prerequisite element "' &
            // TRIM(extended_error_string) // '" has'
        WRITE(40,*) 'not been set. Code will terminate'
        WRITE(40,*)
        errcode_deck = IOR(errcode_deck, c_err_terminate)
      ENDIF
    ENDIF
    IF (IAND(errcode_deck, c_err_pp_options_wrong) .NE. 0) THEN
      IF (rank .EQ. rank_check) THEN
        WRITE(*, *)
        WRITE(*, *) '***WARNING***'
        WRITE(*, *) 'The element "' // TRIM(element) // '" of block "' &
            // TRIM(current_block_name) // '" cannot be set'
        WRITE(*, *) 'because the code has not been compiled with the ' &
            // 'correct preprocessor options.'
        WRITE(*, *) 'Code will continue, but to use selected features, ' &
            // 'please recompile with the'
        WRITE(*, *) TRIM(extended_error_string) // ' option'
        WRITE(*, *)
        WRITE(40,*)
        WRITE(40,*) '***WARNING***'
        WRITE(40,*) 'The element "' // TRIM(element) // '" of block "' &
            // TRIM(current_block_name) // '" cannot be set'
        WRITE(40,*) 'because the code has not been compiled with the ' &
            // 'correct preprocessor options.'
        WRITE(40,*) 'Code will continue, but to use selected features, ' &
            // 'please recompile with the'
        WRITE(40,*) TRIM(extended_error_string) // ' option'
        WRITE(40,*)
      ENDIF
    ENDIF
    IF (IAND(errcode_deck, c_err_other) .NE. 0) THEN
      IF (rank .EQ. rank_check) THEN
        WRITE(*, *)
        WRITE(*, *) '***ERROR***'
        WRITE(*, *) 'You have managed to find an impossible situation in ' &
            // 'this code.'
        WRITE(*, *) 'Good for you. Just because of that, the code will ' &
            // 'terminate.'
        WRITE(*, *)
        WRITE(40,*)
        WRITE(40,*) '***ERROR***'
        WRITE(40,*) 'You have managed to find an impossible situation in ' &
            // 'this code.'
        WRITE(40,*) 'Good for you. Just because of that, the code will ' &
            // 'terminate.'
        WRITE(40,*)
        errcode_deck = IOR(errcode_deck, c_err_terminate)
      ENDIF
    ENDIF

    err_count = err_count+1

  END SUBROUTINE handle_deck_element



  SUBROUTINE write_input_decks

    TYPE(file_buffer), POINTER :: fbuf
    CHARACTER(LEN=1) :: dummy1(1), dummy2
    INTEGER :: i

    IF (rank .EQ. 0) THEN
      fbuf=>file_buffer_head
      DO i = 1,nbuffers
        fbuf=>fbuf%next

        CALL cfd_write_source_code(TRIM(fbuf%filename), "Embedded_input_deck", &
            fbuf%buffer(1:fbuf%idx-1), fbuf%buffer(fbuf%idx)(1:fbuf%pos-1), 0)
      ENDDO
    ELSE
      DO i = 1,nbuffers
        CALL cfd_write_source_code("", "Embedded_input_deck", dummy1, dummy2, 0)
      ENDDO
    ENDIF

  END SUBROUTINE write_input_decks

END MODULE deck
