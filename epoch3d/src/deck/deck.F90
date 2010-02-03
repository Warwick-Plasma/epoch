MODULE deck

  !Basic operations
  USE shared_data
  USE strings
  USE strings_advanced
  !Deck internals
  USE deck_constant_block
  USE deck_deo_block
  !Deck Blocks
  USE deck_control_block
  USE deck_boundaries_block
  USE deck_species_block
  USE deck_io_block
  USE deck_window_block
  !Initial Condition Blocks
  USE deck_ic_laser_block
  USE deck_ic_species_block
  USE deck_ic_fields_block
  USE deck_ic_external_block
  !Extended IO Blocks
  USE deck_eio_dist_fn_block
#ifdef PARTICLE_PROBES
  USE deck_eio_particle_probe_block
#endif
  !Custom blocks
  USE custom_deck

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_deck

  SAVE
  CHARACTER(len=string_length) :: current_block_name
  LOGICAL :: invalid_block
CONTAINS

  !---------------------------------------------------------------------------------
  !These subroutines actually call the routines which read the deck blocks
  !---------------------------------------------------------------------------------

  !This subroutine is called when a new block is started
  !If a block NEEDS to do something when it starts, then
  !The revelevant subroutine should be called here
  SUBROUTINE start_block(block_name)
    CHARACTER(len=*),INTENT(IN) :: block_name

    IF (str_cmp(block_name,"laser")) CALL laser_start
    IF (str_cmp(block_name,"window")) CALL window_start
    IF (str_cmp(block_name,"dist_fn")) CALL dist_fn_start
#ifdef PARTICLE_PROBES
    IF (str_cmp(block_name,"probe")) CALL probe_block_start
#endif

  END SUBROUTINE start_block

  !This subroutine is called when a new block is ended
  !If a block NEEDS to do something when it ends, then
  !The revelevant subroutine should be called here
  SUBROUTINE end_block(block_name)
    CHARACTER(len=*),INTENT(IN) :: block_name

    IF (str_cmp(block_name,"laser")) CALL laser_end
    IF (str_cmp(block_name,"dist_fn")) CALL dist_fn_end
#ifdef PARTICLE_PROBES
    IF (str_cmp(block_name,"probe")) CALL probe_block_end
#endif


  END SUBROUTINE end_block


  FUNCTION handle_block(block_name,block_element,block_value)

    CHARACTER(len=*),INTENT(IN) :: block_name, block_element,block_value
    CHARACTER(len=string_length) :: part1
    INTEGER :: handle_block
    INTEGER :: part2,val

    handle_block=c_err_unknown_block
    !Constants can be defined in any deck state, so put them here
    IF (str_cmp(block_name,"constant")) THEN
       handle_block=handle_constant_deck(block_element,block_value)
       RETURN
    ENDIF
    IF (str_cmp(block_name,"deo")) THEN
       handle_block=handle_deo_deck(block_element,block_value)
       RETURN
    ENDIF
    IF (deck_state .EQ. c_ds_deck) THEN
       !Test for known blocks
       IF (str_cmp(block_name,"control"))  THEN
          handle_block=handle_control_deck(block_element,block_value)
          RETURN
       ENDIF
       IF (str_cmp(block_name,"boundaries")) THEN
          handle_block=handle_boundary_deck(block_element,block_value)
          RETURN
       ENDIF
       IF (str_cmp(block_name,"species")) THEN
          handle_block=handle_species_deck(block_element,block_value)
          RETURN
       ENDIF
       IF (str_cmp(block_name,"output")) THEN
          handle_block=handle_io_deck(block_element,block_value)
          RETURN
       ENDIF
       IF (str_cmp(block_name,"window")) THEN
          handle_block=handle_window_deck(block_element,block_value)
          RETURN
       ENDIF
    ELSE IF (deck_state .EQ. c_ds_ic) THEN
       !Initial conditions blocks go here
       IF (str_cmp(block_name,"fields")) THEN
          handle_block=handle_ic_fields_deck(block_element,block_value)
          RETURN
       ENDIF
       IF (str_cmp(block_name,"fields_external")) THEN
          handle_block=handle_ic_fields_deck(block_element,block_value)
          RETURN
       ENDIF
       IF (str_cmp(block_name,"laser")) THEN
          handle_block=handle_ic_laser_deck(block_element,block_value)
          RETURN
       ENDIF
       val=c_err_none
       CALL split_off_int(block_name,part1,part2,val)
       IF (val .EQ. c_err_none) THEN
          IF (str_cmp(part1,"species")) THEN
             handle_block=handle_ic_species_deck(part2,block_element,block_value)
             RETURN
          ENDIF
          IF (str_cmp(part1,"species_external")) THEN
             handle_block=handle_ic_external_species_deck(part2,block_element,block_value)
             RETURN
          ENDIF
       ENDIF
    ELSE IF (deck_state .EQ. c_ds_eio) THEN
       IF (str_cmp(block_name,"dist_fn")) THEN
          handle_block=handle_eio_dist_fn_deck(block_element,block_value)
          RETURN
       ENDIF
       IF (str_cmp(block_name,"probe")) THEN
#ifdef PARTICLE_PROBES
          handle_block=handle_probe_deck(block_element,block_value)
          RETURN
#else
          handle_block=c_err_pp_options_wrong
          extended_error_string="-DPARTICLE_PROBES"
          RETURN
#endif
       ENDIF
    ENDIF

    !Pass through to the custom block
    handle_block=handle_custom_block(block_name,block_element,block_value)

  END FUNCTION handle_block



  !These subroutines are there to check for the basic minimal compulsory blocks are present
  !They're a bit ugly, but they seem to be the easiest way to do it without adding complexity to the code
  SUBROUTINE check_compulsory_blocks(errcode_deck)

    LOGICAL :: problem_found
    INTEGER,INTENT(INOUT) :: errcode_deck

    problem_found=.FALSE.

    errcode_deck=c_err_none

    IF (deck_state .EQ. c_ds_deck) THEN
       errcode_deck=IOR(errcode_deck,check_control_block())
       errcode_deck=IOR(errcode_deck,check_boundary_block())
       errcode_deck=IOR(errcode_deck,check_species_block())
       errcode_deck=IOR(errcode_deck,check_io_block())
       errcode_deck=IOR(errcode_deck,check_window_block())
       errcode_deck=IOR(errcode_deck,check_custom_blocks())
    ELSE IF (deck_state .EQ. c_ds_ic) THEN
       errcode_deck=IOR(errcode_deck,check_ic_fields_block())
       errcode_deck=IOR(errcode_deck,check_ic_species_block())
    ENDIF
    errcode_deck=IOR(errcode_deck,check_custom_blocks())

    problem_found =(IAND(errcode_deck,c_err_missing_elements) .NE. 0)

    IF (problem_found) THEN
       errcode_deck=IOR(errcode_deck,c_err_terminate)
       IF (rank .EQ. 0) THEN
          PRINT *,""
          PRINT *,"Not all required elements of input deck specified. Please fix input deck and rerun code"
          WRITE(40,*) ""
          WRITE(40,*) "Not all required elements of input deck specified. Please fix input deck and rerun code"
       ENDIF
    ELSE
       IF (rank .EQ. 0) THEN
          IF (deck_state .EQ. c_ds_deck) THEN
             PRINT *,"Input deck complete and valid. Attempting to set up equilibrium"
             PRINT *,""
             WRITE(40,*) "Input deck complete and valid."
          ELSE IF (deck_state .EQ. c_ds_ic) THEN
             PRINT *,"Initial conditions complete and valid. Attempting to load particles"
             PRINT *,""
             WRITE(40,*) "Initial conditions complete and valid."
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE check_compulsory_blocks

  !---------------------------------------------------------------------------------
  !These subroutines are the in depth detail of how the parser works
  !---------------------------------------------------------------------------------
  FUNCTION get_free_lun()

    !This subroutine simply cycles round until it finds a free lun between min_lun and max_lun
    INTEGER :: get_free_lun
    INTEGER :: lun
    INTEGER, PARAMETER :: min_lun=10, max_lun=20
    LOGICAL :: is_open

    is_open=.TRUE.

    lun=min_lun
    DO
       INQUIRE(unit=lun,opened=is_open)
       IF (.NOT. is_open) EXIT
       lun=lun+1
       IF (lun .GT. max_lun) THEN
          IF (rank .EQ. 0) THEN
             PRINT *,"***FATAL ERROR*** unable to open lun for input deck read"
          ENDIF
          CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
       ENDIF
    ENDDO

    get_free_lun=lun

  END FUNCTION get_free_lun

  RECURSIVE SUBROUTINE read_deck(filename,first_call)

    CHARACTER(len=*),INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: first_call
    CHARACTER :: u1
    INTEGER :: pos=1,flip=1,s,f,elements=0,lun
    LOGICAL :: is_comment
    TYPE(string_type), DIMENSION(2) :: deck_values
    CHARACTER(len=45+data_dir_max_length):: deck_filename, status_filename
    LOGICAL :: terminate =.FALSE., exists
    INTEGER :: errcode_deck
    LOGICAL :: white_space_over

    !No error yet
    errcode_deck=c_err_none
    !Characteristic string which represents a "blank" string
    blank="BLANKBLANK"

    lun=5

    !Make the whole filename by adding the data_dir to the filename
    deck_filename=TRIM(ADJUSTL(data_dir))// '/' // TRIM(ADJUSTL(filename))

    !deck_state tells the code whether it's parsing the normal input deck
    !Or the initial conditions. You can add more states if you want.
    !Just search for deck_state
    IF (deck_state .EQ. c_ds_deck) THEN
       status_filename=TRIM(ADJUSTL(data_dir))// '/' // "deck.status"
    ELSE IF (deck_state .EQ. c_ds_ic) THEN
       status_filename=TRIM(ADJUSTL(data_dir))// '/' // "ic.status"
    ELSE IF(deck_state .EQ. c_ds_eio) THEN
       status_filename=TRIM(ADJUSTL(data_dir))// '/' // "eio.status"
    ENDIF

    !If this is the first time that this deck has been called then do some housekeeping
    !Put any initialisation code that is needed in here
    IF (first_call) THEN
       control_block_done=.FALSE.
       boundary_block_done=.FALSE.
    ENDIF

    !Is comment is a flag which tells the code when a # character has been found and everything beyond it
    !Is a comment
    is_comment=.FALSE.

    !rank 0 reads the file and then passes it out to the other nodes using MPI_BCAST
    IF (rank .EQ. 0) THEN
       !Check whether or not the input deck file requested exists
       INQUIRE(file=deck_filename,exist=exists)
       IF (.NOT. exists) THEN
          PRINT *,"***ERROR*** Input deck file ",deck_filename," does not exist. Create the file and rerun the code."
          CALL MPI_ABORT(MPI_COMM_WORLD,errcode)
       ENDIF
       !Get a free lun. Don't use a constant lun to allow for recursion
       lun=get_free_lun()
       OPEN(unit=lun,file=TRIM(ADJUSTL(deck_filename)))
       IF (first_call) OPEN(unit=40,file=status_filename)
       deck_values(1)%value=""
       deck_values(2)%value=""
       !Use non-advancing IO to pop characters off the deck file one at a time
       !Use basic token parsing to split into two substrings across an "=" or ":" symbol
       DO
          errcode_deck=c_err_none
          !Read a character
          !When you reach an EOL character iostat returns -2
          !When you reach an EOF iostat returns -1
          READ(lun,"(A1)",advance='no',size=s,iostat=f),u1
          !If the character is a # then switch to comment mode
          IF (u1 .EQ. '#') is_comment=.TRUE.
          !If not in comment mode then use the character
          IF (.NOT. is_comment) THEN
             !If the current character isn't a special character then just stick it in the buffer
             !             IF (u1 .NE. '=' .AND. u1 .NE. char(32) .AND. u1 .NE. char(9) .AND. u1 .NE. ':' .AND. f .EQ. 0) THEN
             IF (u1 .NE. '=' .AND. u1 .NE. char(9) .AND. u1 .NE. ':' .AND. f .EQ. 0) THEN
                IF ((u1 .NE. ' ' .AND. u1 .NE. char(32)).OR. white_space_over) THEN
                   deck_values(flip)%value(pos:pos)=u1
                   pos=pos+1
                   white_space_over=.TRUE.
                ENDIF
             ENDIF
             !If it's equals or : then you're parsing the other part of the expression
             IF (u1 .EQ. '=' .OR. u1 .EQ. ':') THEN
                flip=2
                pos=1
             ENDIF
          ENDIF
          !If f=-2 then you've reached the end of the line, so comment state is definitely false
          IF (f .EQ. -2) is_comment=.FALSE.
          !If you've not read a blank line then
          IF (f .EQ. -2 .AND. pos .GT. 1) THEN
             elements=elements+1
             flip=1
             pos=1
             deck_values(1)%value=TRIM(ADJUSTL(deck_values(1)%value))
             deck_values(2)%value=TRIM(ADJUSTL(deck_values(2)%value))
             CALL MPI_BCAST(1,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
             CALL MPI_BCAST(deck_values(1)%value,string_length,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
             CALL MPI_BCAST(deck_values(2)%value,string_length,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
             CALL handle_deck_element(deck_values(1)%value,deck_values(2)%value,errcode_deck)
             deck_values(1)%value=""
             deck_values(2)%value=""
             is_comment=.FALSE.
             white_space_over=.FALSE.
          ENDIF
          IF (f .EQ. -1) THEN
             CALL MPI_BCAST(0,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
             CLOSE(lun)
             EXIT
          ENDIF
          terminate=terminate .OR. IAND(errcode_deck,c_err_terminate) .NE. 0
          IF (terminate) EXIT
       ENDDO
    ELSE
       DO
          errcode_deck=c_err_none
          CALL MPI_BCAST(f,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
          IF (f .EQ. 0) EXIT
          CALL MPI_BCAST(deck_values(1)%value,string_length,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
          CALL MPI_BCAST(deck_values(2)%value,string_length,MPI_CHARACTER,0,MPI_COMM_WORLD,errcode)
          CALL handle_deck_element(deck_values(1)%value,deck_values(2)%value,errcode_deck)
          deck_values(1)%value=""
          deck_values(2)%value=""
          terminate=terminate .OR. IAND(errcode_deck,c_err_terminate) .NE. 0
          IF (terminate) EXIT
       ENDDO
    ENDIF


    CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)
!!$    IF (.NOT. first_call)  THEN
!!$       RETURN
!!$    ENDIF

    !Don't check compulsory blocks if going to bomb anyway, just stinks up the output file
    IF (.NOT. terminate .AND. first_call) CALL check_compulsory_blocks(errcode_deck)
    terminate=terminate .OR. IAND(errcode_deck,c_err_terminate) .NE. 0
    !Fatal error, cause code to bomb
    IF (terminate .AND. rank .EQ. 0) THEN
       PRINT *,""
       WRITE(40,*) ""
       PRINT *,'***FATAL ERROR*** The code cannot parse the input deck sufficiently to run. Please read the output file "deck.status" in the current directory to check for errors.'
       WRITE(40,*) "***FATAL ERROR*** The code cannot parse the input deck sufficiently to run. Please read this file and correct any errors mentioned"
       PRINT *,""
       PRINT *,""
       PRINT *,""
       PRINT *,""
       PRINT *,""
       PRINT *,""
    ENDIF

    IF (first_call) CLOSE(40)

    IF (terminate) CALL MPI_ABORT(MPI_COMM_WORLD,errcode)

    CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)


  END SUBROUTINE read_deck

  SUBROUTINE handle_deck_element(element,value,errcode_deck)

    CHARACTER(*),INTENT(IN) :: element
    CHARACTER(*),INTENT(IN) :: value
    INTEGER,INTENT(INOUT) :: errcode_deck
    INTEGER :: state,rank_check
    INTEGER, SAVE :: err_count

    rank_check=0
    state=0

    IF (str_cmp(element,"import")) THEN
       invalid_block=.TRUE.
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*),""
          WRITE(40,*),"Importing ",TRIM(ADJUSTL(value)), " file"
          WRITE(40,*),""
       ENDIF
       CALL read_deck(TRIM(ADJUSTL(value)),.FALSE.)
       RETURN
    ENDIF

    IF (str_cmp(element,"begin")) THEN
       errcode_deck=handle_block(value,blank,blank)
       invalid_block=IAND(errcode_deck,c_err_unknown_block) .NE. 0
       IF(invalid_block) THEN
          IF (rank .EQ. rank_check) THEN
             PRINT *,char(9),"Unknown block ",TRIM(value)," in input deck, ignoring"
          ENDIF
       ENDIF
       CALL start_block(value)
       err_count=0
       current_block_name=value
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*),"Beginning ", TRIM(ADJUSTL(value)), " block"
          WRITE(40,*),""
       ENDIF
       !Reset errcode_deck here because reporting c_err_unknown_element is OK
       errcode_deck=c_err_none
       RETURN
    ENDIF
    IF (str_cmp(element,"end")) THEN
       CALL end_block(current_block_name)
       invalid_block=.TRUE.
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*),""
          WRITE(40,*),"Ending ",TRIM(ADJUSTL(value)), " block"
          WRITE(40,*),""
          IF (err_count .NE. 0) THEN
             WRITE(40,*) "***WARNING*** block ",TRIM(ADJUSTL(value))," contains errors"
             WRITE(40,*) ""
          ENDIF
       ENDIF
       RETURN
    ENDIF

    !Check invalid block to limit amount of rubbish that appears
    !If the input deck is invalid
    IF (.NOT. invalid_block) THEN
       errcode_deck=handle_block(current_block_name,element,value)
    ELSE
       RETURN
    ENDIF

    IF (errcode_deck==c_err_none) THEN
       IF (rank .EQ. rank_check) WRITE(40,*),char(9),"Element ", TRIM(ADJUSTL(element))," = ",TRIM(ADJUSTL(value)), " handled OK"
       RETURN
    ENDIF
    !Test for error conditions
    !If an error is fatal then set terminate to .TRUE.
    IF (IAND(errcode_deck,c_err_unknown_element) /= 0) THEN
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***WARNING*** Unrecognised element ",TRIM(element), " in input deck. Code will continue to run, but behaviour is undefined"
          WRITE(40,*) "***WARNING*** Unrecognised element ",TRIM(element), " in input deck. Code will continue to run, but behaviour is undefined"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(errcode_deck,c_err_preset_element) /= 0) THEN
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*) ""
          PRINT *,"***WARNING*** element ",TRIM(element), " is set multiple times in this deck. Code will continue using first value in deck"
          WRITE(40,*) "***WARNING*** element ",TRIM(element), " is set multiple times in this deck. Code will continue using first value in deck"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(errcode_deck, c_err_preset_element_use_later) /= 0) THEN
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***WARNING*** element ",TRIM(element), " is set multiple times in this deck. Code will continue using last value in deck"
          WRITE(40,*) "***WARNING*** element ",TRIM(element), " is set multiple times in this deck. Code will continue using last value in deck"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(errcode_deck, c_err_bad_value) /= 0) THEN
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** value ",TRIM(value)," in element ",TRIM(element)," is invalid or could not be parsed. Code will terminate."
          WRITE(40,*) "***ERROR*** value ",TRIM(value)," in element ",TRIM(element)," is invalid or could not be parsed. Code will terminate."
          PRINT *,""
          WRITE(40,*) ""
          errcode_deck=IOR(errcode_deck,c_err_terminate)
       ENDIF
    ENDIF
!!$    IF (IAND(errcode_deck, ERR_BAD_VALUE_NO_TERMINATE) /= 0) THEN
!!$       IF (rank .EQ. rank_check) THEN
!!$          WRITE(40,*) ""
!!$          PRINT *,""
!!$          PRINT *,"***ERROR*** value ",TRIM(value)," in non essential element ",TRIM(element)," is invalid or could not be parsed. Code will continue but behaviour is undefined."
!!$          WRITE(40,*) "***ERROR*** value ",TRIM(value)," in non essential element ",TRIM(element)," is invalid or could not be parsed. Code will continue but behaviour is undefined."
!!$          PRINT *,""
!!$          WRITE(40,*) ""
!!$       ENDIF
!!$    ENDIF
    IF (IAND(errcode_deck, c_err_required_element_not_set) /= 0) THEN
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** value ",TRIM(value)," in element ",TRIM(element)," cannot be set because a prerequisite element, ", TRIM(extended_error_string),",has not been set. Code will terminate"
          WRITE(40,*) "***ERROR*** value ",TRIM(value)," in element ",TRIM(element)," cannot be set because a prerequisite element, ", TRIM(extended_error_string),",has not been set. Code will terminate"
          PRINT *,""
          WRITE(40,*) ""
          errcode_deck=IOR(errcode_deck,c_err_terminate)
       ENDIF
    ENDIF
    IF (IAND(errcode_deck , c_err_pp_options_wrong) /= 0) THEN
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** The element ",TRIM(element)," of block ",TRIM(current_block_name)," cannot be set because the code has not been compiled with the correct preprocessor options.",&
               "Code will continue, but to use selected features, please recompile with ",TRIM(extended_error_string)," option"
          WRITE(40,*) "***ERROR*** The element ",TRIM(element)," of block ",TRIM(current_block_name)," cannot be set because the code has not been compiled with the correct preprocessor options.",&
               "Code will continue, but to use selected features, please recompile with ",TRIM(extended_error_string)," option"
          PRINT *,""
          WRITE(40,*) ""
       ENDIF
    ENDIF
    IF (IAND(errcode_deck , c_err_other) /= 0) THEN
       IF (rank .EQ. rank_check) THEN
          WRITE(40,*) ""
          PRINT *,""
          PRINT *,"***ERROR*** You have managed to find an impossible situation in this code. Good for you. Just because of that, code will terminate."
          WRITE(40,*) "***ERROR*** You have managed to find an impossible situation in this code. Good for you. Just because of that, code will terminate."
          PRINT *,""
          WRITE(40,*) ""
          errcode_deck=IOR(errcode_deck,c_err_terminate)
       ENDIF
    ENDIF

    err_count=err_count+1

  END SUBROUTINE handle_deck_element

END MODULE deck
