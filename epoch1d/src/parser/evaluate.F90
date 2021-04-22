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

MODULE evaluator

  USE shunt
  USE stack

  IMPLICIT NONE

  TYPE stack_list
    TYPE(primitive_stack) :: stack
    TYPE(stack_list), POINTER :: prev, next
  END TYPE stack_list

  TYPE(stack_list), POINTER :: sl_head, sl_tail
  INTEGER :: sl_size

CONTAINS

#ifdef SIMPLIFY_DEBUG
  SUBROUTINE basic_evaluate(input_stack, parameters, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    TYPE(parameter_pack), INTENT(IN) :: parameters
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i, n_elements
    REAL(num), ALLOCATABLE :: array(:)

    IF (input_stack%should_simplify) THEN
      CALL basic_evaluate_standard(input_stack, parameters)

      n_elements = eval_stack_stack_point
      ALLOCATE(array(1:n_elements))

      ! Pop off the final answers
      DO i = n_elements,1,-1
        array(i) = pop_off_eval()
      END DO

      CALL simplify_stack(input_stack, err)

      CALL basic_evaluate_standard(input_stack, parameters)

      ! Check the final answers
      DO i = 1, n_elements
        IF (ABS(eval_stack_entries(i) - array(i)) > c_tiny) THEN
          PRINT*,i,eval_stack_entries(i),array(i),eval_stack_entries(i)-array(i)
        END IF
      END DO

      DEALLOCATE(array)
    ELSE
      CALL basic_evaluate_standard(input_stack, parameters)
    END IF

  END SUBROUTINE basic_evaluate



  SUBROUTINE basic_evaluate_standard(input_stack, parameters)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    TYPE(parameter_pack), INTENT(IN) :: parameters
    INTEGER :: i, err, ispec
    TYPE(stack_element) :: iblock

    CALL eval_reset()

    err = c_err_none
    DO i = 1, input_stack%stack_point
      iblock = input_stack%entries(i)
      IF (iblock%ptype == c_pt_variable) THEN
        CALL push_on_eval(iblock%numerical_data)
      ELSE IF (iblock%ptype == c_pt_species) THEN
        CALL do_species(iblock%value, err)
      ELSE IF (iblock%ptype == c_pt_operator) THEN
        CALL do_operator(iblock%value, err)
      ELSE IF (iblock%ptype == c_pt_constant &
          .OR. iblock%ptype == c_pt_default_constant) THEN
        CALL do_constant(iblock%value, .FALSE., parameters, err)
      ELSE IF (iblock%ptype == c_pt_function) THEN
        CALL do_functions(iblock%value, .FALSE., parameters, err)
        IF (err < 0) THEN
          ispec = -err
          CALL do_evaluate(iblock%value, .FALSE., parameters, ispec, err)
        END IF
      END IF

      IF (err /= c_err_none) CALL abort_block(iblock, i, err)
    END DO

  END SUBROUTINE basic_evaluate_standard



#else
  SUBROUTINE basic_evaluate(input_stack, parameters, error)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    TYPE(parameter_pack), INTENT(IN) :: parameters
    INTEGER, INTENT(INOUT) :: error
    INTEGER :: i, err, ispec
    TYPE(stack_element) :: iblock

    IF (input_stack%should_simplify) CALL simplify_stack(input_stack, error)

    CALL eval_reset()

    err = c_err_none
    DO i = 1, input_stack%stack_point
      iblock = input_stack%entries(i)
      IF (iblock%ptype == c_pt_variable) THEN
        CALL push_on_eval(iblock%numerical_data)
      ELSE IF (iblock%ptype == c_pt_species) THEN
        CALL do_species(iblock%value, err)
      ELSE IF (iblock%ptype == c_pt_operator) THEN
        CALL do_operator(iblock%value, err)
      ELSE IF (iblock%ptype == c_pt_constant &
          .OR. iblock%ptype == c_pt_default_constant) THEN
        CALL do_constant(iblock%value, .FALSE., parameters, err)
      ELSE IF (iblock%ptype == c_pt_function) THEN
        CALL do_functions(iblock%value, .FALSE., parameters, err)
        IF (err < 0) THEN
          ispec = -err
          CALL do_evaluate(iblock%value, .FALSE., parameters, ispec, err)
        END IF
      END IF

      IF (err /= c_err_none) CALL abort_block(iblock, i, err)
    END DO

  END SUBROUTINE basic_evaluate



#endif
  SUBROUTINE sl_append

    TYPE(stack_list), POINTER :: sl_tmp

    IF (sl_size == 0) THEN
      ALLOCATE(sl_head)
      sl_tail => sl_head
    ELSE
      ALLOCATE(sl_tmp)
      sl_tmp%prev => sl_tail
      sl_tail%next => sl_tmp
      sl_tail => sl_tmp
    END IF
    sl_size = sl_size + 1
    CALL initialise_stack(sl_tail%stack)

  END SUBROUTINE sl_append



  SUBROUTINE simplify_stack(input_stack, error)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(INOUT) :: error
    INTEGER :: i, err, ispec
    TYPE(stack_element) :: iblock
    TYPE(primitive_stack) :: output_stack
    TYPE(parameter_pack) :: parameters

    parameters%pack_ix = 1

    ! Evaluating expressions and push the results onto eval_stack.
    ! When we reach a block which requests a time or space variable,
    ! push a special flag to the eval_stack and create a new stack containing
    ! the block (sl_tail%stack).
    ! If do_* encounters an expression whose arguments contain a space or
    ! time varying expression (eg. gauss(x+y,0,1)), it pushes garbage to the
    ! eval_stack and a flag is set. It is then fixed up by
    ! update_stack_for_block()

    ! Eg. the following expression:
    ! 2*3*x + func1(4+5,x*func2(6*7,8*9*x^2),x+y) + func3(1*2,3,func4(4+5))+1+2
    ! simplifies to:
    ! 6*x + func1(9,x*func2(42,72*x^2),x+y) + f3 + 1 + 2
    !   where f3=func3(2,3,func4(9))

    CALL eval_reset()

    sl_size = 0

    err = c_err_none
    DO i = 1, input_stack%stack_point
      iblock = input_stack%entries(i)
      IF (iblock%ptype == c_pt_variable) THEN
        CALL push_on_eval(iblock%numerical_data)
      ELSE IF (iblock%ptype == c_pt_species) THEN
        CALL do_species(iblock%value, err)
      ELSE IF (iblock%ptype == c_pt_operator) THEN
        CALL do_operator(iblock%value, err)
        CALL update_stack_for_block(iblock, err)
      ELSE IF (iblock%ptype == c_pt_constant &
          .OR. iblock%ptype == c_pt_default_constant) THEN
        CALL do_constant(iblock%value, .TRUE., parameters, err)
        CALL update_stack_for_block(iblock, err)
      ELSE IF (iblock%ptype == c_pt_function) THEN
        CALL do_functions(iblock%value, .TRUE., parameters, err)
        IF (err < 0) THEN
          ispec = -err
          CALL do_evaluate(iblock%value, .TRUE., parameters, ispec, err)
        END IF
        CALL update_stack_for_block(iblock, err)
      END IF

      IF (err /= c_err_none) THEN
        PRINT *, 'BAD block', err, iblock%ptype, i, iblock%value
        CALL abort_code(err)
        STOP
      END IF
    END DO

    ! We may now just be left with a list of values on the eval_stack
    ! If so, push them onto sl_tail%stack
    i = eval_stack_stack_point
    IF (i > 0) CALL update_stack(i)

    ! Now populate output_stack with the simplified expression
    CALL initialise_stack(output_stack)
    output_stack%should_simplify = .FALSE.

    IF (sl_size > 0) THEN
      CALL append_stack(output_stack, sl_tail%stack)
      DEALLOCATE(sl_tail)
      sl_size = 0
      eval_stack_stack_point = 0
    END IF

    CALL deallocate_stack(input_stack)
    input_stack = output_stack

  END SUBROUTINE simplify_stack



  SUBROUTINE update_stack(nvalues)

    INTEGER, INTENT(IN) :: nvalues
    INTEGER :: sp, i, n
    REAL(num), ALLOCATABLE :: entries(:)
    INTEGER, ALLOCATABLE :: flags(:)
    TYPE(stack_list), POINTER :: sl_tmp, sl_part
    TYPE(stack_element) :: new_block

    new_block%ptype = c_pt_variable
    new_block%value = 0

    sp = eval_stack_stack_point

    ALLOCATE(entries(nvalues))
    ALLOCATE(flags(nvalues))

    n = nvalues
    DO i = 1, nvalues
      entries(n) = eval_stack_entries(sp)
      flags(n) = eval_stack_flags(sp)
      IF (flags(n) /= 0) THEN
        sl_size = sl_size - 1
        sl_part => sl_tail
        sl_tail => sl_tail%prev
      END IF
      n = n - 1
      sp = sp - 1
    END DO
    eval_stack_stack_point = sp

    CALL sl_append()
    DO i = 1, nvalues
      IF (flags(i) == 0) THEN
        new_block%numerical_data = entries(i)
        CALL push_to_stack(sl_tail%stack, new_block)
      ELSE
        CALL append_stack(sl_tail%stack, sl_part%stack)
        sl_tmp => sl_part%next
        DEALLOCATE(sl_part)
        sl_part => sl_tmp
      END IF
    END DO

    DEALLOCATE(entries)
    DEALLOCATE(flags)

  END SUBROUTINE update_stack



  SUBROUTINE update_stack_for_block(iblock, err)

    TYPE(stack_element), INTENT(INOUT) :: iblock
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: nvalues

    IF (err == c_err_other .OR. err == c_err_window) THEN
      ! Operator just pushed a bogus value to stack, so we'll ignore it
      eval_stack_stack_point = eval_stack_stack_point - 1
      CALL push_eval_flag()
      CALL sl_append()
      CALL push_to_stack(sl_tail%stack, iblock)
      IF (iblock%value == c_const_time .OR. err == c_err_window) THEN
        sl_tail%stack%is_time_varying = .TRUE.
      END IF
      err = c_err_none
      RETURN
    END IF

    ! Number of eval_stack entries consumed by operator
    nvalues = eval_stack_nvalues
    IF (nvalues == 0) RETURN

    eval_stack_nvalues = 0
    ! Operator just pushed a bogus value to stack, so we'll ignore it
    eval_stack_stack_point = eval_stack_stack_point - 1

    CALL update_stack(nvalues)

    CALL push_eval_flag()
    CALL push_to_stack(sl_tail%stack, iblock)

  END SUBROUTINE update_stack_for_block



  SUBROUTINE evaluate_with_parameters_to_array(input_stack, parameters, &
      n_elements, array, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    TYPE(parameter_pack), INTENT(IN) :: parameters
    INTEGER, INTENT(IN) :: n_elements
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i

    CALL basic_evaluate(input_stack, parameters, err)

    IF (eval_stack_stack_point /= n_elements) THEN
      PRINT*, '*** ERROR ***'
      PRINT*, 'Input deck line number ', TRIM(deck_line_number)
      PRINT*, 'Incorrect number of arguments.'
      err = IOR(err, c_err_bad_value)
      CALL abort_code(err)
    END IF

    ! Pop off the final answers
    DO i = MIN(eval_stack_stack_point,n_elements),1,-1
      array(i) = pop_off_eval()
    END DO

  END SUBROUTINE evaluate_with_parameters_to_array



  SUBROUTINE evaluate_to_array(input_stack, n_elements, array, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(IN) :: n_elements
    REAL(num), DIMENSION(:), INTENT(INOUT) :: array
    INTEGER, INTENT(INOUT) :: err
    TYPE(parameter_pack) :: parameters

    parameters%pack_ix = 1

    CALL evaluate_with_parameters_to_array(input_stack, parameters, &
        n_elements, array, err)

  END SUBROUTINE evaluate_to_array



  SUBROUTINE evaluate_and_return_all_with_parameters(input_stack, parameters, &
      n_elements, array, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    TYPE(parameter_pack), INTENT(IN) :: parameters
    INTEGER, INTENT(OUT) :: n_elements
    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(INOUT) :: err
    INTEGER :: i

    IF (ASSOCIATED(array)) DEALLOCATE(array)

    CALL basic_evaluate(input_stack, parameters, err)

    n_elements = eval_stack_stack_point
    ALLOCATE(array(1:n_elements))

    ! Pop off the final answers
    DO i = n_elements,1,-1
      array(i) = pop_off_eval()
    END DO

  END SUBROUTINE evaluate_and_return_all_with_parameters



  SUBROUTINE evaluate_and_return_all(input_stack, n_elements, array, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(OUT) :: n_elements
    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(INOUT) :: err
    TYPE(parameter_pack) :: parameters

    CALL evaluate_and_return_all_with_parameters(input_stack, parameters, &
        n_elements, array, err)

  END SUBROUTINE evaluate_and_return_all



  SUBROUTINE evaluate_as_list(input_stack, array, n_elements)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, DIMENSION(:), INTENT(OUT) :: array
    INTEGER, INTENT(OUT) :: n_elements
    INTEGER :: i, err
    TYPE(stack_element) :: iblock
    TYPE(parameter_pack) :: parameters

    array(1) = 0
    n_elements = 1

    err = c_err_none
    DO i = 1, input_stack%stack_point
      iblock = input_stack%entries(i)

      IF (iblock%ptype == c_pt_subset) THEN
        n_elements = n_elements + 1
        array(n_elements) = iblock%value
      ELSE IF (iblock%ptype == c_pt_constant &
          .OR. iblock%ptype == c_pt_default_constant) THEN
        CALL do_constant(iblock%value, .FALSE., parameters, err)
        array(1) = array(1) + INT(pop_off_eval())
      ELSE IF (iblock%ptype /= c_pt_operator) THEN
        err = c_err_bad_value
      END IF

      IF (err /= c_err_none) CALL abort_block(iblock, i, err)
    END DO

  END SUBROUTINE evaluate_as_list



  FUNCTION evaluate_with_parameters(input_stack, parameters, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    TYPE(parameter_pack), INTENT(IN) :: parameters
    INTEGER, INTENT(INOUT) :: err
    REAL(num), DIMENSION(1) :: array
    REAL(num) :: evaluate_with_parameters

    CALL evaluate_with_parameters_to_array(input_stack, parameters, 1, &
        array, err)
    evaluate_with_parameters = array(1)

  END FUNCTION evaluate_with_parameters



  FUNCTION evaluate(input_stack, err)

    TYPE(primitive_stack), INTENT(INOUT) :: input_stack
    INTEGER, INTENT(INOUT) :: err
    REAL(num) :: evaluate
    TYPE(parameter_pack) :: parameters

    evaluate = evaluate_with_parameters(input_stack, parameters, err)

  END FUNCTION evaluate



  SUBROUTINE do_evaluate(opcode, simplify, parameters, ispec, err)

    ! Initial conditions arrays might not be allocated when we try to access
    ! them. In that situation, we need to evaluate the functional form instead.
    ! Routine must be defined here to remove circular dependencies.
    INTEGER, INTENT(IN) :: opcode
    TYPE(parameter_pack), INTENT(IN) :: parameters
    LOGICAL, INTENT(IN) :: simplify
    INTEGER, INTENT(INOUT) :: ispec, err
    REAL(num) :: val, val_local, convert
    INTEGER :: n, err_simplify
    INTEGER :: ix
    TYPE(parameter_pack) :: param
    TYPE(particle_species), POINTER :: species
#include "particle_head.inc"

    err = c_err_none
    err_simplify = c_err_none

    IF (simplify) err_simplify = c_err_other

    IF (opcode == c_func_rho) THEN
      species => species_list(ispec)

      IF (parameters%use_grid_position) THEN
        val_local = evaluate_with_parameters(species%density_function, &
                  parameters, errcode)
      ELSE
#include "pack_to_grid.inc"
        val_local = 0.0_num
        DO ix = sf_min, sf_max
          param%pack_ix = cell_x + ix
          val = evaluate_with_parameters(species%density_function, &
                  param, errcode)
          val_local = val_local + gx(ix) * val
        END DO
      END IF

      CALL push_on_eval(val_local)
      err = err_simplify
      RETURN
    END IF

    n = 0
    IF (opcode == c_func_tempx) THEN
      n = 1
      convert = 1.0_num
    ELSE IF (opcode == c_func_tempy) THEN
      n = 2
      convert = 1.0_num
    ELSE IF (opcode == c_func_tempz) THEN
      n = 3
      convert = 1.0_num
    ELSE IF (opcode == c_func_tempx_ev) THEN
      n = 1
      convert = kb / ev
    ELSE IF (opcode == c_func_tempy_ev) THEN
      n = 2
      convert = kb / ev
    ELSE IF (opcode == c_func_tempz_ev) THEN
      n = 3
      convert = kb / ev
    END IF

    IF (n > 0) THEN
      species => species_list(ispec)

      IF (parameters%use_grid_position) THEN
        val_local = evaluate_with_parameters(species%temperature_function(n), &
                  parameters, errcode)
      ELSE
#include "pack_to_grid.inc"
        val_local = 0.0_num
        DO ix = sf_min, sf_max
          param%pack_ix = cell_x + ix
          val = evaluate_with_parameters(species%temperature_function(n), &
                  param, errcode)
          val_local = val_local + gx(ix) * val
        END DO
      END IF

      CALL push_on_eval(convert * val_local)
      err = err_simplify
      RETURN
    END IF

    err = c_err_unknown_element

  END SUBROUTINE do_evaluate



  SUBROUTINE abort_block(block, sp, err)

    TYPE(stack_element), INTENT(IN) :: block
    INTEGER, INTENT(IN) :: sp, err
    INTEGER :: i

    PRINT*, '*** ERROR ***'
    PRINT*, 'Unable to parse block'
    PRINT*, 'Input deck line number ', TRIM(deck_line_number)
    PRINT*, 'Error code:'
    DO i = 0, c_err_max
      IF (IAND(err, 2**i) /= 0) PRINT*, '    ', TRIM(c_err_char(i))
    END DO
    PRINT*, 'Block type: ', TRIM(c_pt_char(block%ptype))
    PRINT*, 'Block value: ', block%value
    PRINT*, 'Stack point: ', sp
    CALL abort_code(err)
    STOP

  END SUBROUTINE abort_block

END MODULE evaluator
