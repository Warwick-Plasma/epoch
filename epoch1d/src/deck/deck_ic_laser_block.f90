MODULE deck_ic_laser_block

  USE shared_data
  USE strings_advanced
  USE shared_parser_data
  USE strings
  USE shunt
  USE evaluator
  USE laser

  IMPLICIT NONE

  SAVE

  TYPE(laser_block),POINTER :: working_laser
  LOGICAL :: direction_set=.FALSE.
  INTEGER :: direction

CONTAINS

  FUNCTION handle_ic_laser_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_ic_laser_deck
    REAL(num) :: dummy

    handle_ic_laser_deck=ERR_NONE
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF(str_cmp(element,"direction")) THEN
       !If the direction has already been set, simply ignore further calls to it
       IF (direction_set) RETURN
       direction=as_direction(value,handle_ic_laser_deck)
       direction_set=.TRUE.
       CALL init_laser(direction,working_laser)
       RETURN
    ENDIF

    IF (.NOT. direction_set) THEN
       IF (rank .EQ. 0) THEN
          PRINT *,"***ERROR*** Cannot set laser properties before direction is set"
          WRITE(40,*) "***ERROR*** Cannot set laser properties before direction is set"
       ENDIF
       extended_error_string="direction"
       handle_ic_laser_deck=ERR_REQUIRED_ELEMENT_NOT_SET
       RETURN
    ENDIF

    IF (str_cmp(element,"amp")) THEN
       working_laser%amp=as_real(value,handle_ic_laser_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"freq")) THEN
       working_laser%freq=as_real(value,handle_ic_laser_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"phase")) THEN
       working_laser%phase=as_real(value,handle_ic_laser_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"t_start")) THEN
       working_laser%t_start=as_time(value,handle_ic_laser_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"t_end")) THEN
       working_laser%t_end=as_time(value,handle_ic_laser_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"t_profile")) THEN
       working_laser%use_time_function=.TRUE.
       working_laser%time_function%stack_point=0
       CALL tokenize(value,working_laser%time_function,handle_ic_laser_deck)
       !evaluate it once to check that it's a valid block
       dummy=evaluate(working_laser%time_function,handle_ic_laser_deck)
       RETURN
    ENDIF
    IF (str_cmp(element,"pol")) THEN
       !Convert from degrees to radians
       working_laser%pol=pi * as_real(value,handle_ic_laser_deck)/180.0_num
       RETURN
    ENDIF
    IF (str_cmp(element,"id")) THEN
       working_laser%id=as_integer(value,handle_ic_laser_deck)
       RETURN
    ENDIF
    handle_ic_laser_deck=ERR_UNKNOWN_ELEMENT

  END FUNCTION handle_ic_laser_deck

  FUNCTION check_ic_laser_block()

    INTEGER :: check_ic_laser_block

    !Should do error checking but can't be bothered at the moment
    check_ic_laser_block=ERR_NONE

  END FUNCTION check_ic_laser_block

  SUBROUTINE laser_start
    !Every new laser uses the internal time function
    ALLOCATE(working_laser)
    working_laser%use_time_function=.FALSE.
  END SUBROUTINE laser_start

  SUBROUTINE laser_end
    CALL attach_laser(working_laser)
    direction_set=.FALSE.
  END SUBROUTINE laser_end

  FUNCTION as_time(value,err)

    CHARACTER(len=*),INTENT(IN) :: value
    INTEGER,INTENT(INOUT) :: err
    REAL(num) :: as_time

    IF (str_cmp(value,"start")) THEN
       as_time=0.0_num
       RETURN
    ENDIF

    IF (str_cmp(value,"end")) THEN
       as_time=t_end
       RETURN
    ENDIF

    as_time=as_real(value,err)

  END FUNCTION as_time

END MODULE deck_ic_laser_block
