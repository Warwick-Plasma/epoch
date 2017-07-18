! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2017 Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE deck_injector_block

  USE strings_advanced
  USE injectors
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: injector_deck_initialise, injector_deck_finalise
  PUBLIC :: injector_block_start, injector_block_end
  PUBLIC :: injector_block_handle_element, injector_block_check

  TYPE(injector_block), POINTER :: working_injector
  LOGICAL :: boundary_set = .FALSE.
  INTEGER :: boundary

CONTAINS

  SUBROUTINE injector_deck_initialise

  END SUBROUTINE injector_deck_initialise



  SUBROUTINE injector_deck_finalise

  END SUBROUTINE injector_deck_finalise



  SUBROUTINE injector_block_start

    IF (deck_state == c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_injector)

  END SUBROUTINE injector_block_start



  SUBROUTINE injector_block_end

    IF (deck_state == c_ds_first) RETURN

    CALL attach_injector(working_injector)
    boundary_set = .FALSE.

  END SUBROUTINE injector_block_end



  FUNCTION injector_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dummy
    INTEGER :: io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'boundary')) THEN
      ! If the boundary has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary = as_boundary_print(value, element, errcode)
      boundary_set = .TRUE.
      CALL init_injector(boundary, working_injector)
      RETURN
    ENDIF

    IF (str_cmp(element, 'species')) THEN
      working_injector%species = as_integer_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 't_start')) THEN
      working_injector%t_start = as_time_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 't_end')) THEN
      working_injector%t_end = as_time_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'npart_per_cell')) THEN
      working_injector%npart_per_cell = as_integer_print(value, element, &
          errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'density')) THEN
      working_injector%density = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_x')) THEN
      working_injector%temperature(1) = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_y')) THEN
      working_injector%temperature(2) = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp_z')) THEN
      working_injector%temperature(3) = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'temp')) THEN
      working_injector%temperature = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_x')) THEN
      working_injector%drift(1) = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_y')) THEN
      working_injector%drift(2) = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift_z')) THEN
      working_injector%drift(3) = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    IF (str_cmp(element, 'drift')) THEN
      working_injector%drift = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element

  END FUNCTION injector_block_handle_element



  FUNCTION injector_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(injector_block), POINTER :: current
    INTEGER :: error, io, iu

    errcode = c_err_none
    current => injector_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%species ==-1) error = IOR(error, 1)
      current => current%next
    ENDDO

    current => injector_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%species ==-1) error = IOR(error, 1)
      current => current%next
    ENDDO

    IF (IAND(error, 1) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a species for every injector.'
        ENDDO
      ENDIF
      errcode = c_err_missing_elements
    ENDIF

  END FUNCTION injector_block_check

END MODULE deck_injector_block
