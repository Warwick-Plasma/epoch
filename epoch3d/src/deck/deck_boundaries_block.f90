MODULE deck_boundaries_block

  USE strings_advanced

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: boundary_deck_initialise, boundary_deck_finalise
  PUBLIC :: boundary_block_start, boundary_block_end
  PUBLIC :: boundary_block_handle_element, boundary_block_check

  INTEGER, PARAMETER :: boundary_block_nbase = 2 * c_ndims
  INTEGER, PARAMETER :: boundary_block_elements = 3 * boundary_block_nbase + 4
  LOGICAL, DIMENSION(boundary_block_elements) :: boundary_block_done
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      boundary_block_name = (/ &
          'bc_x_min         ', &
          'bc_x_max         ', &
          'bc_y_min         ', &
          'bc_y_max         ', &
          'bc_z_min         ', &
          'bc_z_max         ', &
          'bc_x_min_field   ', &
          'bc_x_max_field   ', &
          'bc_y_min_field   ', &
          'bc_y_max_field   ', &
          'bc_z_min_field   ', &
          'bc_z_max_field   ', &
          'bc_x_min_particle', &
          'bc_x_max_particle', &
          'bc_y_min_particle', &
          'bc_y_max_particle', &
          'bc_z_min_particle', &
          'bc_z_max_particle', &
          'cpml_thickness   ', &
          'cpml_kappa_max   ', &
          'cpml_a_max       ', &
          'cpml_sigma_max   ' /)
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      alternate_name = (/ &
          'xbc_left          ', &
          'xbc_right         ', &
          'ybc_down          ', &
          'ybc_up            ', &
          'zbc_back          ', &
          'zbc_front         ', &
          'xbc_left_field    ', &
          'xbc_right_field   ', &
          'ybc_down_field    ', &
          'ybc_up_field      ', &
          'zbc_back_field    ', &
          'zbc_front_field   ', &
          'xbc_left_particle ', &
          'xbc_right_particle', &
          'ybc_down_particle ', &
          'ybc_up_particle   ', &
          'zbc_back_particle ', &
          'zbc_front_particle', &
          'cpml_thickness    ', &
          'cpml_kappa_max    ', &
          'cpml_a_max        ', &
          'cpml_sigma_max    ' /)

CONTAINS

  SUBROUTINE boundary_deck_initialise

    IF (deck_state .NE. c_ds_first) RETURN
    boundary_block_done = .FALSE.

  END SUBROUTINE boundary_deck_initialise



  SUBROUTINE boundary_deck_finalise

  END SUBROUTINE boundary_deck_finalise



  SUBROUTINE boundary_block_start

  END SUBROUTINE boundary_block_start



  SUBROUTINE boundary_block_end

  END SUBROUTINE boundary_block_end



  FUNCTION boundary_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: loop, elementselected, itmp
    INTEGER, PARAMETER :: nbase = boundary_block_nbase

    errcode = c_err_none
    IF (deck_state .NE. c_ds_first) RETURN

    errcode = c_err_unknown_element

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
      errcode = c_err_preset_element
      RETURN
    ENDIF
    boundary_block_done(elementselected) = .TRUE.
    errcode = c_err_none

    IF (elementselected .LE. nbase) THEN
      boundary_block_done(elementselected+  nbase) = .TRUE.
      boundary_block_done(elementselected+2*nbase) = .TRUE.
    ENDIF

    SELECT CASE (elementselected)
    CASE(1)
      itmp = as_bc(value, errcode)
      bc_field(c_bd_x_min) = itmp
      bc_particle(c_bd_x_min) = itmp
    CASE(2)
      itmp = as_bc(value, errcode)
      bc_field(c_bd_x_max) = itmp
      bc_particle(c_bd_x_max) = itmp
    CASE(3)
      itmp = as_bc(value, errcode)
      bc_field(c_bd_y_min) = itmp
      bc_particle(c_bd_y_min) = itmp
    CASE(4)
      itmp = as_bc(value, errcode)
      bc_field(c_bd_y_max) = itmp
      bc_particle(c_bd_y_max) = itmp
    CASE(5)
      itmp = as_bc(value, errcode)
      bc_field(c_bd_z_min) = itmp
      bc_particle(c_bd_z_min) = itmp
    CASE(6)
      itmp = as_bc(value, errcode)
      bc_field(c_bd_z_max) = itmp
      bc_particle(c_bd_z_max) = itmp
    CASE(nbase+1)
      bc_field(c_bd_x_min) = as_bc(value, errcode)
      boundary_block_done(1)  = .TRUE.
    CASE(nbase+2)
      bc_field(c_bd_x_max) = as_bc(value, errcode)
      boundary_block_done(2)  = .TRUE.
    CASE(nbase+3)
      bc_field(c_bd_y_min) = as_bc(value, errcode)
      boundary_block_done(3)  = .TRUE.
    CASE(nbase+4)
      bc_field(c_bd_y_max) = as_bc(value, errcode)
      boundary_block_done(4)  = .TRUE.
    CASE(nbase+5)
      bc_field(c_bd_z_min) = as_bc(value, errcode)
      boundary_block_done(5)  = .TRUE.
    CASE(nbase+6)
      bc_field(c_bd_z_max) = as_bc(value, errcode)
      boundary_block_done(6)  = .TRUE.
    CASE(2*nbase+1)
      bc_particle(c_bd_x_min) = as_bc(value, errcode)
      boundary_block_done(1)  = .TRUE.
    CASE(2*nbase+2)
      bc_particle(c_bd_x_max) = as_bc(value, errcode)
      boundary_block_done(2)  = .TRUE.
    CASE(2*nbase+3)
      bc_particle(c_bd_y_min) = as_bc(value, errcode)
      boundary_block_done(3)  = .TRUE.
    CASE(2*nbase+4)
      bc_particle(c_bd_y_max) = as_bc(value, errcode)
      boundary_block_done(4)  = .TRUE.
    CASE(2*nbase+5)
      bc_particle(c_bd_z_min) = as_bc(value, errcode)
      boundary_block_done(5)  = .TRUE.
    CASE(2*nbase+6)
      bc_particle(c_bd_z_max) = as_bc(value, errcode)
      boundary_block_done(6)  = .TRUE.
    CASE(3*nbase+1)
      cpml_thickness = as_integer(value, errcode)
    CASE(3*nbase+2)
      cpml_kappa_max = as_real(value, errcode)
    CASE(3*nbase+3)
      cpml_a_max = as_real(value, errcode)
    CASE(3*nbase+4)
      cpml_sigma_max = as_real(value, errcode)
    END SELECT

  END FUNCTION boundary_block_handle_element



  FUNCTION boundary_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: index, io, ierr
    INTEGER, PARAMETER :: nbase = boundary_block_nbase
    LOGICAL :: error

    errcode = c_err_none

    DO index = 1, nbase
      IF (.NOT.boundary_block_done(index+nbase) &
          .AND. .NOT.boundary_block_done(index+2*nbase)) THEN
        boundary_block_done(index+  nbase) = .TRUE.
        boundary_block_done(index+2*nbase) = .TRUE.
      ENDIF
    ENDDO

    DO index = 1, boundary_block_elements - 4
      IF (.NOT. boundary_block_done(index)) THEN
        IF (rank .EQ. 0) THEN
          DO io = stdout, du, du - stdout ! Print to stdout and to file
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required boundary block element "' &
                // TRIM(ADJUSTL(boundary_block_name(index))) // '" absent.'
            WRITE(io,*) 'Please create this entry in the input deck'
          ENDDO
        ENDIF
        errcode = c_err_missing_elements
      ENDIF
    ENDDO

    ! Sanity check on periodic boundaries
    error = .FALSE.
    DO index = 1, c_ndims
      IF (bc_field(2*index-1) .EQ. c_bc_periodic &
          .AND. bc_field(2*index) .NE. c_bc_periodic) error = .TRUE.
      IF (bc_field(2*index-1) .NE. c_bc_periodic &
          .AND. bc_field(2*index) .EQ. c_bc_periodic) error = .TRUE.
      IF (bc_particle(2*index-1) .EQ. c_bc_periodic &
          .AND. bc_particle(2*index) .NE. c_bc_periodic) error = .TRUE.
      IF (bc_particle(2*index-1) .NE. c_bc_periodic &
          .AND. bc_particle(2*index) .EQ. c_bc_periodic) error = .TRUE.
    ENDDO

    IF (error) THEN
      IF (rank .EQ. 0) THEN
        DO io = stdout, du, du - stdout ! Print to stdout and to file
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Periodic boundaries must be specified on both sides', &
              ' of the domain.'
        ENDDO
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

  END FUNCTION boundary_block_check

END MODULE deck_boundaries_block
