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

MODULE deck_solid_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: solid_deck_initialise, solid_deck_finalise
  PUBLIC :: solid_block_start, solid_block_end
  PUBLIC :: solid_block_handle_element, solid_block_check

CONTAINS

  SUBROUTINE solid_deck_initialise

  END SUBROUTINE solid_deck_initialise



  SUBROUTINE solid_deck_finalise

  END SUBROUTINE solid_deck_finalise



  SUBROUTINE solid_block_start

    INTEGER :: io, iu
#ifdef HYBRID
    INTEGER :: isolid

    IF (deck_state == c_ds_first) THEN
      ! Count number of solids in first pass
      solid_count = solid_count + 1
    ELSE
      IF (.NOT. made_solid_array) THEN
        ! Create array to hold solids
        ALLOCATE(solid_array(solid_count))
        made_solid_array = .TRUE.

        ! Allocate arrays for each solid
        DO isolid = 1, solid_count
          ALLOCATE(solid_array(isolid)%ion_density(1-ng:nx+ng,1-ng:ny+ng, &
              1-ng:nz+ng))
          ALLOCATE(solid_array(isolid)%el_density(1-ng:nx+ng,1-ng:ny+ng, &
              1-ng:nz+ng))
        END DO
      ELSE
        solid_index = solid_index + 1
      END IF
    END IF
#else
    IF (rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'Unable to add a solid to the simulation.'
        WRITE(io,*) 'Please recompile with the -DHYBRID preprocessor flag.'
      END DO
    END IF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE solid_block_start



  SUBROUTINE solid_block_end

    INTEGER :: io, iu

    IF (deck_state == c_ds_first) RETURN

#ifdef HYBRID
    IF (.NOT. use_hybrid) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Solids can only be modelled when running in hybrid mode.'
          WRITE(io,*) 'To use the solid block, set use_hybrid = T in the', &
              ' hybrid block.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE solid_block_end



  FUNCTION solid_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

#ifdef HYBRID

    IF (str_cmp(element, 'atomic_no') &
        .OR. str_cmp(element, 'background_Z') &
        .OR. str_cmp(element, 'hybrid_Z')) THEN
      solid_array(solid_index)%z = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mass_no') &
        .OR. str_cmp(element, 'A')) THEN
      solid_array(solid_index)%mass_no = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'I') &
        .OR. str_cmp(element, 'I_ex') &
        .OR. str_cmp(element, 'excitation_energy')) THEN
      solid_array(solid_index)%iex = &
          as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'ion_density') &
        .OR. str_cmp(element, 'rho' ) &
        .OR. str_cmp(element, 'ni')) THEN
      CALL fill_array(solid_array(solid_index)%ion_density, value)
      RETURN
    END IF

    IF (str_cmp(element, 'resistivity') &
        .OR. str_cmp(element, 'resistivity_model' ) &
        .OR. str_cmp(element, 'res_model')) THEN
      IF (str_cmp(value, 'vacuum') .OR. str_cmp(value, 'Vacuum')) THEN
        solid_array(solid_index)%res_model = c_resist_vacuum
      ELSE IF (str_cmp(value, 'milchberg') &
        .OR. str_cmp(value, 'Milchberg')) THEN
        solid_array(solid_index)%res_model = c_resist_milchberg
      ELSE IF (str_cmp(value, 'plastic') .OR. str_cmp(value, 'Plastic')) THEN
        solid_array(solid_index)%res_model = c_resist_plastic
      ELSE IF (str_cmp(value, 'rlm')) THEN
        solid_array(solid_index)%res_model = c_resist_rlm
      ELSE
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) ''
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Unrecognised resistivity model: ', TRIM(value)
            WRITE(io,*) 'Code will ignore return currents from this solid'
            WRITE(io,*) ''
          END DO
        END IF
      END IF
      RETURN
    END IF

    errcode = c_err_unknown_element
#endif

  END FUNCTION solid_block_handle_element



  FUNCTION solid_block_check() RESULT(errcode)

    INTEGER :: errcode

    errcode = c_err_none

  END FUNCTION solid_block_check



  SUBROUTINE fill_array(array, value)

    ! A simplified version of the script in deck_species_block. It evaluates the
    ! equation string stored in 'value' (for the maths parser to interpret), and
    ! writes the elements to 'array'.

    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: array
    CHARACTER(LEN=*), INTENT(IN) :: value
    TYPE(stack_element) :: iblock
    TYPE(primitive_stack) :: stack
    INTEGER :: io, iu, ix, iy, iz
    TYPE(parameter_pack) :: parameters

    CALL initialise_stack(stack)
    CALL tokenize(value, stack, errcode)

    ! Sanity check
    array(1,1,1) = evaluate(stack, errcode)
    IF (errcode /= c_err_none) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to parse input deck.'
        END DO
      END IF
      CALL abort_code(errcode)
    END IF

    DO iz = 1-ng, nz+ng
      parameters%pack_iz = iz
      DO iy = 1-ng, ny+ng
        parameters%pack_iy = iy
        DO ix = 1-ng, nx+ng
          parameters%pack_ix = ix
          array(ix,iy,iz) = evaluate_with_parameters(stack, parameters, errcode)
        END DO
      END DO
    END DO

  END SUBROUTINE fill_array

END MODULE deck_solid_block
