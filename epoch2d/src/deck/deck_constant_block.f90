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

MODULE deck_constant_block

  USE shunt
  USE evaluator

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: constant_deck_initialise, constant_deck_finalise
  PUBLIC :: constant_block_start, constant_block_end
  PUBLIC :: constant_block_handle_element, constant_block_check

CONTAINS

  SUBROUTINE constant_deck_initialise

    INTEGER :: i

    IF (ALLOCATED(deck_constant_list)) THEN
      DO i = 1, n_deck_constants
        CALL deallocate_stack(deck_constant_list(i)%execution_stream)
      END DO
      DEALLOCATE(deck_constant_list)
    END IF
    n_deck_constants = 0

  END SUBROUTINE constant_deck_initialise



  SUBROUTINE constant_deck_finalise

    INTEGER :: i, errcode
    REAL(num) :: dc
    LOGICAL :: const_is_open

    IF (.NOT.print_deck_constants) RETURN
    IF (rank /= 0) RETURN

    IF (deck_state == c_ds_first) THEN
      WRITE(du,*) 'Constant block values after first pass:'
    ELSE
      WRITE(du,*) 'Constant block values after second pass:'
    END IF
    WRITE(du,*)

    INQUIRE(unit=duc, opened=const_is_open)

    DO i = 1, n_deck_constants
      errcode = 0
      dc = evaluate(deck_constant_list(i)%execution_stream, errcode)
      WRITE(du,'("  ", A, " = ", G18.11)') TRIM(deck_constant_list(i)%name), dc
      IF (const_is_open) THEN
        WRITE(duc,'(A, " = ", G18.11)') TRIM(deck_constant_list(i)%name), dc
      END IF
    END DO

  END SUBROUTINE constant_deck_finalise



  SUBROUTINE constant_block_start

  END SUBROUTINE constant_block_start



  SUBROUTINE constant_block_end

  END SUBROUTINE constant_block_end



  FUNCTION constant_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: ix, io, iu
    TYPE(deck_constant), DIMENSION(:), ALLOCATABLE :: buffer
    TYPE(primitive_stack) :: temp
    TYPE(stack_element) :: iblock

    errcode = c_err_none

    IF (value == blank) RETURN

    ! First check whether constant already exists
    DO ix = 1, n_deck_constants
      IF (str_cmp(TRIM(element), TRIM(deck_constant_list(ix)%name))) RETURN
    END DO

    ! If we're here then then named constant doesn't yet exist, so create it

    ! First issue a warning message if the name overrides a built-in one
    CALL load_block(element, iblock)
    IF (iblock%ptype /= c_pt_bad .AND. iblock%ptype /= c_pt_null &
        .AND. iblock%ptype /= c_pt_default_constant) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'The constant block variable "' // TRIM(element) &
              // '" conflicts with a'
          WRITE(io,*) 'built-in constant name. It will be ignored.'
          WRITE(io,*)
        END DO
      END IF
      RETURN
    END IF

    CALL initialise_stack(temp)
    CALL tokenize(value, temp, errcode)
    IF (errcode /= c_err_none) THEN
      CALL deallocate_stack(temp)
      RETURN
    END IF

    ! Take a copy of the old list
    IF (n_deck_constants > 0) THEN
      ALLOCATE(buffer(1:n_deck_constants))
      buffer = deck_constant_list

      ! Allocate the new list
      DEALLOCATE(deck_constant_list)
      ALLOCATE(deck_constant_list(1:n_deck_constants+1))

      deck_constant_list(1:n_deck_constants) = buffer
      DEALLOCATE(buffer)
    ELSE
      ! Allocate the new list
      ALLOCATE(deck_constant_list(1:n_deck_constants+1))
    END IF

    ! Add the new value
    deck_constant_list(n_deck_constants+1)%execution_stream = temp
    deck_constant_list(n_deck_constants+1)%name = element

    n_deck_constants = n_deck_constants + 1

  END FUNCTION constant_block_handle_element



  FUNCTION constant_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION constant_block_check

END MODULE deck_constant_block
