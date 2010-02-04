MODULE deck_eio_dist_fn_block

  USE shared_data
  USE strings_advanced
  USE dist_fn

  IMPLICIT NONE

  SAVE

  TYPE(distribution_function_block),POINTER :: working_block
CONTAINS

  FUNCTION handle_eio_dist_fn_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_eio_dist_fn_deck

    CHARACTER(len=string_length) :: part1
    INTEGER :: part2
    INTEGER :: work
    REAL(num) :: work1,work2

    handle_eio_dist_fn_deck=ERR_NONE
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element,"name")) THEN
       working_block%name=value
       RETURN
    ENDIF

    IF (str_cmp(element,"ndims")) THEN
       work=as_integer(value,handle_eio_dist_fn_deck)
       IF (work .EQ. 2 .OR. work .EQ. 3) THEN
          working_block%ndims=work
       ELSE
          IF (rank .EQ. 0) PRINT *,"Distribution functions can only be 2D or 3D"
          handle_eio_dist_fn_deck=ERR_BAD_VALUE
       ENDIF
       RETURN
    ENDIF
    IF (working_block%ndims .EQ. -1) THEN
       IF (rank .EQ. 0) PRINT *,"Must set number of dimensions before setting other distribution function properties."
       extended_error_string="ndims"
       handle_eio_dist_fn_deck=ERR_REQUIRED_ELEMENT_NOT_SET
       RETURN
    ENDIF

    IF (str_cmp(element,"dumpmask")) THEN
       working_block%dumpmask=as_integer(value,handle_eio_dist_fn_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"restrict_x")) THEN
       CALL split_range(value,work1,work2,handle_eio_dist_fn_deck)
       IF (handle_eio_dist_fn_deck .NE. ERR_NONE) RETURN
       working_block%use_restrictions(1)=.TRUE.
       working_block%restrictions(1,:)=(/work1,work2/)
    ENDIF
    IF (str_cmp(element,"restrict_px")) THEN
       CALL split_range(value,work1,work2,handle_eio_dist_fn_deck)
       IF (handle_eio_dist_fn_deck .NE. ERR_NONE) RETURN
       working_block%use_restrictions(2)=.TRUE.
       working_block%restrictions(2,:)=(/work1,work2/)
    ENDIF
    IF (str_cmp(element,"restrict_py")) THEN
       CALL split_range(value,work1,work2,handle_eio_dist_fn_deck)
       IF (handle_eio_dist_fn_deck .NE. ERR_NONE) RETURN
       working_block%use_restrictions(3)=.TRUE.
       working_block%restrictions(3,:)=(/work1,work2/)
    ENDIF
    IF (str_cmp(element,"restrict_pz")) THEN
       CALL split_range(value,work1,work2,handle_eio_dist_fn_deck)
       IF (handle_eio_dist_fn_deck .NE. ERR_NONE) RETURN
       working_block%use_restrictions(4)=.TRUE.
       working_block%restrictions(4,:)=(/work1,work2/)
    ENDIF


    CALL split_off_int(element,part1,part2,handle_eio_dist_fn_deck)
    IF (handle_eio_dist_fn_deck .NE. ERR_NONE) THEN
       handle_eio_dist_fn_deck = ERR_UNKNOWN_ELEMENT
       RETURN
    ENDIF
    IF (str_cmp(part1,"direction")) THEN
       working_block%directions(part2)=as_real(value,handle_eio_dist_fn_deck)
       RETURN
    ENDIF
    IF (str_cmp(part1,"range")) THEN
       CALL split_range(TRIM(value),work1,work2,handle_eio_dist_fn_deck)
       IF (handle_eio_dist_fn_deck .NE. ERR_NONE) RETURN
       working_block%ranges(part2,1)=work1
       working_block%ranges(part2,2)=work2
       RETURN
    ENDIF
    IF (str_cmp(part1,"resolution")) THEN
       working_block%resolution(part2)=as_integer(value,handle_eio_dist_fn_deck)
       RETURN
    ENDIF
    IF (str_cmp(part1,"include_species_")) THEN
       IF (part2 .LT. 1 .OR. part2 .GT. n_species) THEN
          IF (rank .EQ. 0) PRINT *,"Species ",part2," does not exist, ignoring attempt to set output state."
          handle_eio_dist_fn_deck=ERR_NONE
          RETURN
       ENDIF
       working_block%use_species(part2)=as_logical(value,handle_eio_dist_fn_deck)
       RETURN
    ENDIF


    handle_eio_dist_fn_deck=ERR_UNKNOWN_ELEMENT

  END FUNCTION handle_eio_dist_fn_deck

  FUNCTION check_eio_dist_fn_block()

    INTEGER :: check_eio_dist_fn_block

    !Should do error checking but can't be bothered at the moment
    check_eio_dist_fn_block=ERR_NONE

  END FUNCTION check_eio_dist_fn_block

  SUBROUTINE dist_fn_start
    !Every new laser uses the internal time function
    ALLOCATE(working_block)
    CALL setup_dist_fn(working_block)

  END SUBROUTINE dist_fn_start

  SUBROUTINE dist_fn_end

    CALL attach_dist_fn(working_block)

  END SUBROUTINE dist_fn_end

END MODULE deck_eio_dist_fn_block
