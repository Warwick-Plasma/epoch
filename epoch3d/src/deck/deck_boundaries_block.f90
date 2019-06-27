! Copyright (C) 2009-2019 University of Warwick
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

MODULE deck_boundaries_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: boundary_deck_initialise, boundary_deck_finalise
  PUBLIC :: boundary_block_start, boundary_block_end
  PUBLIC :: boundary_block_handle_element, boundary_block_check

CONTAINS

  SUBROUTINE boundary_deck_initialise

    IF (deck_state /= c_ds_first) RETURN

    bc_field(:) = -1
    bc_particle(:) = -1

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
    INTEGER :: itmp

    errcode = c_err_none
    IF (deck_state /= c_ds_first) RETURN

    IF (str_cmp(element, 'bc_x_min') &
        .OR. str_cmp(element, 'xbc_left')) THEN
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_x_min) = itmp
      bc_particle(c_bd_x_min) = itmp

    ELSE IF (str_cmp(element, 'bc_x_max') &
        .OR. str_cmp(element, 'xbc_right')) THEN
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_x_max) = itmp
      bc_particle(c_bd_x_max) = itmp

    ELSE IF (str_cmp(element, 'bc_y_min') &
        .OR. str_cmp(element, 'ybc_down')) THEN
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_y_min) = itmp
      bc_particle(c_bd_y_min) = itmp

    ELSE IF (str_cmp(element, 'bc_y_max') &
        .OR. str_cmp(element, 'ybc_up')) THEN
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_y_max) = itmp
      bc_particle(c_bd_y_max) = itmp

    ELSE IF (str_cmp(element, 'bc_z_min') &
        .OR. str_cmp(element, 'zbc_back')) THEN
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_z_min) = itmp
      bc_particle(c_bd_z_min) = itmp

    ELSE IF (str_cmp(element, 'bc_z_max') &
        .OR. str_cmp(element, 'zbc_front')) THEN
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_z_max) = itmp
      bc_particle(c_bd_z_max) = itmp

    ELSE IF (str_cmp(element, 'bc_x_min_field') &
        .OR. str_cmp(element, 'xbc_left_field')) THEN
      bc_field(c_bd_x_min) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_x_max_field') &
        .OR. str_cmp(element, 'xbc_right_field')) THEN
      bc_field(c_bd_x_max) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_y_min_field') &
        .OR. str_cmp(element, 'ybc_down_field')) THEN
      bc_field(c_bd_y_min) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_y_max_field') &
        .OR. str_cmp(element, 'ybc_up_field')) THEN
      bc_field(c_bd_y_max) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_z_min_field') &
        .OR. str_cmp(element, 'zbc_back_field')) THEN
      bc_field(c_bd_z_min) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_z_max_field') &
        .OR. str_cmp(element, 'zbc_front_field')) THEN
      bc_field(c_bd_z_max) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_x_min_particle') &
        .OR. str_cmp(element, 'xbc_left_particle')) THEN
      bc_particle(c_bd_x_min) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_x_max_particle') &
        .OR. str_cmp(element, 'xbc_right_particle')) THEN
      bc_particle(c_bd_x_max) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_y_min_particle') &
        .OR. str_cmp(element, 'ybc_down_particle')) THEN
      bc_particle(c_bd_y_min) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_y_max_particle') &
        .OR. str_cmp(element, 'ybc_up_particle')) THEN
      bc_particle(c_bd_y_max) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_z_min_particle') &
        .OR. str_cmp(element, 'zbc_back_particle')) THEN
      bc_particle(c_bd_z_min) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'bc_z_max_particle') &
        .OR. str_cmp(element, 'zbc_front_particle')) THEN
      bc_particle(c_bd_z_max) = as_bc_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'cpml_thickness')) THEN
      cpml_thickness = as_integer_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'cpml_kappa_max')) THEN
      cpml_kappa_max = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'cpml_a_max')) THEN
      cpml_a_max = as_real_print(value, element, errcode)

    ELSE IF (str_cmp(element, 'cpml_sigma_max')) THEN
      cpml_sigma_max = as_real_print(value, element, errcode)

    ELSE
      errcode = c_err_unknown_element

    END IF

  END FUNCTION boundary_block_handle_element



  FUNCTION boundary_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: idx, io, iu
    LOGICAL :: error
    CHARACTER(LEN=*), DIMENSION(6), PARAMETER :: &
        part_name =  (/ 'bc_x_min_particle', 'bc_x_max_particle', &
                        'bc_y_min_particle', 'bc_y_max_particle', &
                        'bc_z_min_particle', 'bc_z_max_particle' /)
    CHARACTER(LEN=*), DIMENSION(6), PARAMETER :: &
        field_name = (/ 'bc_x_min_field', 'bc_x_max_field', &
                        'bc_y_min_field', 'bc_y_max_field', &
                        'bc_z_min_field', 'bc_z_max_field' /)

    errcode = c_err_none

    DO idx = 1, 2 * c_ndims
      IF (bc_particle(idx) == -1) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required boundary block element "' &
                // part_name(idx) // '" absent.'
            WRITE(io,*) 'Please create this entry in the input deck'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
      IF (bc_field(idx) == -1) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required boundary block element "' &
                // field_name(idx) // '" absent.'
            WRITE(io,*) 'Please create this entry in the input deck'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
    END DO

    ! Sanity check on periodic boundaries
    error = .FALSE.
    DO idx = 1, c_ndims
      IF (bc_field(2*idx-1) == c_bc_periodic &
          .AND. bc_field(2*idx) /= c_bc_periodic) error = .TRUE.
      IF (bc_field(2*idx-1) /= c_bc_periodic &
          .AND. bc_field(2*idx) == c_bc_periodic) error = .TRUE.
      IF (bc_particle(2*idx-1) == c_bc_periodic &
          .AND. bc_particle(2*idx) /= c_bc_periodic) error = .TRUE.
      IF (bc_particle(2*idx-1) /= c_bc_periodic &
          .AND. bc_particle(2*idx) == c_bc_periodic) error = .TRUE.
    END DO

    IF (error) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Periodic boundaries must be specified on both sides', &
              ' of the domain.'
        END DO
      END IF
      CALL abort_code(c_err_bad_value)
    END IF

    ! Sanity check on heat_bath boundaries
    error = .FALSE.
    DO idx = 1, 2 * c_ndims
      IF (bc_field(idx) == c_bc_heat_bath) error = .TRUE.
    END DO

    IF (error) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'heat_bath boundaries apply to particles only'
        END DO
      END IF
      CALL abort_code(c_err_bad_value)
    END IF

  END FUNCTION boundary_block_check

END MODULE deck_boundaries_block
