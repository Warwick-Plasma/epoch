! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE deck_hybrid_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: hybrid_deck_initialise, hybrid_deck_finalise
  PUBLIC :: hybrid_block_start, hybrid_block_end
  PUBLIC :: hybrid_block_handle_element, hybrid_block_check

CONTAINS

  SUBROUTINE hybrid_deck_initialise

  END SUBROUTINE hybrid_deck_initialise



  SUBROUTINE hybrid_deck_finalise

    INTEGER :: io, iu
#ifdef HYBRID
    IF (deck_state == c_ds_first) RETURN
#else
    IF (use_hybrid) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_hybrid=T" in the', &
              ' "hybrid" block.'
          WRITE(io,*) 'Please recompile with the -DHYBRID preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE hybrid_deck_finalise



  SUBROUTINE hybrid_block_start

  END SUBROUTINE hybrid_block_start



  SUBROUTINE hybrid_block_end

  END SUBROUTINE hybrid_block_end



  FUNCTION hybrid_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_hybrid') &
        .OR. str_cmp(element, 'hybrid')) THEN
      use_hybrid = as_logical_print(value, element, errcode)
      RETURN
    END IF

    errcode = c_err_unknown_element

  END FUNCTION hybrid_block_handle_element



  FUNCTION hybrid_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef HYBRID
    INTEGER :: io, iu
#endif

    errcode = c_err_none

  END FUNCTION hybrid_block_check

END MODULE deck_hybrid_block
