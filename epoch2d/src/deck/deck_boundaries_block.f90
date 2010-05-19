MODULE deck_boundaries_block

  USE strings

  IMPLICIT NONE

  SAVE

  INTEGER, PARAMETER :: boundary_block_nbase = 2 * c_ndims
  INTEGER, PARAMETER :: boundary_block_elements = 3 * boundary_block_nbase
  LOGICAL, DIMENSION(boundary_block_elements) :: boundary_block_done
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      boundary_block_name = (/ &
          "bc_x_min         ", &
          "bc_x_max         ", &
          "bc_y_min         ", &
          "bc_y_max         ", &
          "bc_x_min_field   ", &
          "bc_x_max_field   ", &
          "bc_y_min_field   ", &
          "bc_y_max_field   ", &
          "bc_x_min_particle", &
          "bc_x_max_particle", &
          "bc_y_min_particle", &
          "bc_y_max_particle" /)
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      alternate_name = (/ &
          "xbc_left          ", &
          "xbc_right         ", &
          "ybc_down          ", &
          "ybc_up            ", &
          "xbc_left_field    ", &
          "xbc_right_field   ", &
          "ybc_down_field    ", &
          "ybc_up_field      ", &
          "xbc_left_particle ", &
          "xbc_right_particle", &
          "ybc_down_particle ", &
          "ybc_up_particle   " /)

CONTAINS

  FUNCTION handle_boundary_deck(element, value)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: handle_boundary_deck
    INTEGER :: loop, elementselected, itmp
    INTEGER, PARAMETER :: nbase = boundary_block_nbase

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

    IF (elementselected .LE. nbase) THEN
      boundary_block_done(elementselected+  nbase) = .TRUE.
      boundary_block_done(elementselected+2*nbase) = .TRUE.
    ENDIF

    SELECT CASE (elementselected)
    CASE(1)
      itmp = as_bc(value, handle_boundary_deck)
      bc_x_min_field = itmp
      bc_x_min_particle = itmp
    CASE(2)
      itmp = as_bc(value, handle_boundary_deck)
      bc_x_max_field = itmp
      bc_x_max_particle = itmp
    CASE(3)
      itmp = as_bc(value, handle_boundary_deck)
      bc_y_min_field = itmp
      bc_y_min_particle = itmp
    CASE(4)
      itmp = as_bc(value, handle_boundary_deck)
      bc_y_max_field = itmp
      bc_y_max_particle = itmp
    CASE(nbase+1)
      bc_x_min_field = as_bc(value, handle_boundary_deck)
      boundary_block_done(1)  = .TRUE.
    CASE(nbase+2)
      bc_x_max_field = as_bc(value, handle_boundary_deck)
      boundary_block_done(2)  = .TRUE.
    CASE(nbase+3)
      bc_y_min_field = as_bc(value, handle_boundary_deck)
      boundary_block_done(3)  = .TRUE.
    CASE(nbase+4)
      bc_y_max_field = as_bc(value, handle_boundary_deck)
      boundary_block_done(4)  = .TRUE.
    CASE(2*nbase+1)
      bc_x_min_particle = as_bc(value, handle_boundary_deck)
      boundary_block_done(1)  = .TRUE.
    CASE(2*nbase+2)
      bc_x_max_particle = as_bc(value, handle_boundary_deck)
      boundary_block_done(2)  = .TRUE.
    CASE(2*nbase+3)
      bc_y_min_particle = as_bc(value, handle_boundary_deck)
      boundary_block_done(3)  = .TRUE.
    CASE(2*nbase+4)
      bc_y_max_particle = as_bc(value, handle_boundary_deck)
      boundary_block_done(4)  = .TRUE.
    END SELECT

  END FUNCTION handle_boundary_deck



  FUNCTION check_boundary_block()

    INTEGER :: check_boundary_block
    INTEGER :: index
    INTEGER, PARAMETER :: nbase = boundary_block_nbase

    check_boundary_block = c_err_none

    DO index = 1, nbase
      IF (.NOT.boundary_block_done(index+nbase) &
          .AND. .NOT.boundary_block_done(index+2*nbase)) THEN
        boundary_block_done(index+  nbase) = .TRUE.
        boundary_block_done(index+2*nbase) = .TRUE.
      ENDIF
    ENDDO

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
