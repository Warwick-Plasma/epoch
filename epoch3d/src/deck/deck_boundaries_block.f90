MODULE deck_boundaries_block

  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: boundary_block_elements = 6
  LOGICAL, DIMENSION(boundary_block_elements) :: boundary_block_done
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      boundary_block_name = (/ &
          "bc_x_min", "bc_x_max", "bc_y_min", "bc_y_max", &
          "bc_z_min", "bc_z_max" /)
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      alternate_name = (/ &
          "xbc_left ", "xbc_right", "ybc_down ", "ybc_up   ", &
          "zbc_back ", "zbc_front" /)

CONTAINS

  FUNCTION handle_boundary_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_boundary_deck
    INTEGER :: loop, elementselected

    handle_boundary_deck = c_err_unknown_element

    elementselected = 0

    DO loop = 1, boundary_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(boundary_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (boundary_block_done(elementselected)) THEN
      handle_boundary_deck = c_err_preset_element
      RETURN
    ENDIF
    boundary_block_done(elementselected) = .TRUE.
    handle_boundary_deck = c_err_none

    SELECT CASE (elementselected)
    CASE(1)
      bc_x_min = as_bc(value, handle_boundary_deck)
    CASE(2)
      bc_x_max = as_bc(value, handle_boundary_deck)
    CASE(3)
      bc_y_min = as_bc(value, handle_boundary_deck)
    CASE(4)
      bc_y_max = as_bc(value, handle_boundary_deck)
    CASE(5)
      bc_z_min = as_bc(value, handle_boundary_deck)
    CASE(6)
      bc_z_max = as_bc(value, handle_boundary_deck)
    END SELECT

  END FUNCTION handle_boundary_deck



  FUNCTION check_boundary_block()

    INTEGER :: check_boundary_block
    INTEGER :: index

    check_boundary_block = c_err_none

    DO index = 1, boundary_block_elements
      IF (.NOT. boundary_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          WRITE(*, *)
          WRITE(*, *) '***ERROR***'
          WRITE(*, *) 'Required boundary block element "' &
              // TRIM(ADJUSTL(boundary_block_name(index))) // '" absent.'
          WRITE(*, *) 'Please create this entry in the input deck'
          WRITE(40,*)
          WRITE(40,*) '***ERROR***'
          WRITE(40,*) 'Required boundary block element "' &
              // TRIM(ADJUSTL(boundary_block_name(index))) // '" absent.'
          WRITE(40,*) 'Please create this entry in the input deck'
        ENDIF
        check_boundary_block = c_err_missing_elements
      ENDIF
    ENDDO

  END FUNCTION check_boundary_block

END MODULE deck_boundaries_block
