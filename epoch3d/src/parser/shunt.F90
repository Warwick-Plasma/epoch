MODULE shunt

  USE tokenizer_blocks

  IMPLICIT NONE

  SAVE

CONTAINS

  FUNCTION char_type(char)

    CHARACTER, INTENT(IN) :: char
    INTEGER :: char_type

    char_type = c_char_unknown

    IF (char .EQ. ' ') THEN
      char_type = c_char_space
    ELSEIF (char .GE. '0' .AND. char .LE. '9' .OR. char .EQ. '.') THEN
      char_type = c_char_numeric
    ELSEIF ((char .GE. 'A' .AND. char .LE. 'Z') &
        .OR. (char .GE. 'a' .AND. char .LE. 'z') .OR. char .EQ. '_') THEN
      char_type = c_char_alpha
    ELSEIF (char .EQ. '(' .OR. char .EQ. ')' .OR. char .EQ. ',') THEN
      char_type = c_char_delimiter
    ! 92 is the ASCII code for backslash
    ELSEIF (char .EQ. '+' .OR. char .EQ. '-' .OR. ICHAR(char) .EQ. 92 &
        .OR. char .EQ. '/' .OR. char .EQ. '*' .OR. char .EQ. '^') THEN
      char_type = c_char_opcode
    ENDIF

  END FUNCTION char_type



  SUBROUTINE load_block(name, block)

    CHARACTER(LEN=*), INTENT(IN) :: name
    TYPE(stack_element), INTENT(OUT) :: block
    INTEGER :: work
    REAL(num) :: value

    block%ptype = c_pt_bad
    block%value = 0
    block%numerical_data = 0.0_num
    work = 0

    IF (LEN(TRIM(name)) .EQ. 0) THEN
      block%ptype = c_pt_null
      block%value = 0
      block%numerical_data = 0.0_num
      RETURN
    ENDIF

    work = as_constant(name)
    IF (work .NE. 0) THEN
      ! block is a named constant
      block%ptype = c_pt_constant
      block%value = work
      RETURN
    ENDIF

    work = as_deck_constant(name)
    IF (work .NE. 0) THEN
      ! block is a deck constant
      block%ptype = c_pt_deck_constant
      block%value = work
      RETURN
    ENDIF

    work = as_operator(name)
    IF (work .NE. 0) THEN
      ! block is an operator
      block%ptype = c_pt_operator
      block%value = work
      RETURN
    ENDIF

    work = as_function(name)
    IF (work .NE. 0) THEN
      ! block is a function
      block%ptype = c_pt_function
      block%value = work
      RETURN
    ENDIF

    work = as_parenthesis(name)
    IF (work .NE. 0) THEN
      ! block is a parenthesis
      block%ptype = c_pt_parenthesis
      block%value = work
      RETURN
    ENDIF

    work = as_species(name)
    IF (work .NE. 0) THEN
      ! block is a species name
      block%ptype = c_pt_species
      block%value = work
      RETURN
    ENDIF

    work = as_subset(name)
    IF (work .NE. 0) THEN
      ! block is a subset name
      block%ptype = c_pt_subset
      block%value = work
      RETURN
    ENDIF

    IF (str_cmp(name, ', ')) THEN
      block%ptype = c_pt_separator
      block%value = 0
      block%numerical_data = 0.0_num
      RETURN
    ENDIF

    value = as_real_simple(name, work)
    IF (IAND(work, c_err_bad_value) .EQ. 0) THEN
      ! block is a simple variable
      block%ptype = c_pt_variable
      block%value = 0
      block%numerical_data = value
    ENDIF

  END SUBROUTINE load_block



  FUNCTION as_parenthesis(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_parenthesis

    as_parenthesis = 0

    IF (str_cmp(name, '(')) THEN
      as_parenthesis = c_paren_left_bracket
    ENDIF

    IF (str_cmp(name, ')')) THEN
      as_parenthesis = c_paren_right_bracket
    ENDIF

  END FUNCTION as_parenthesis



  SUBROUTINE initialise_stack(stack)

    TYPE(primitive_stack), INTENT(INOUT) :: stack

    stack%stack_point = 0

    stack%stack_size = 1
    ALLOCATE(stack%entries(stack%stack_size))
    stack%init = .TRUE.

  END SUBROUTINE initialise_stack



  SUBROUTINE deallocate_stack(stack)

    TYPE(primitive_stack), INTENT(INOUT) :: stack

    IF (.NOT.stack%init .AND. rank .EQ. 0) THEN
      PRINT*,'*** WARNING ***'
      PRINT*,'deallocate_stack not initialised'
    ENDIF

    stack%stack_point = 0
    stack%stack_size = 0
    IF (stack%init) DEALLOCATE(stack%entries)
    stack%init = .FALSE.

  END SUBROUTINE deallocate_stack



  SUBROUTINE copy_stack(stack, copy)

    TYPE(primitive_stack), INTENT(IN) :: stack
    TYPE(primitive_stack), INTENT(OUT) :: copy

    copy = stack
    ALLOCATE(copy%entries(copy%stack_size))
    copy%entries(1:copy%stack_point) = stack%entries(1:copy%stack_point)

  END SUBROUTINE copy_stack



  SUBROUTINE push_to_stack(stack, value)

    TYPE(stack_element), INTENT(IN) :: value
    TYPE(primitive_stack), INTENT(INOUT) :: stack
    TYPE(stack_element), POINTER :: old_buffer(:)
    INTEGER :: old_size

    stack%stack_point = stack%stack_point + 1

    IF (.NOT.stack%init .AND. rank .EQ. 0) THEN
      PRINT*,'*** WARNING ***'
      PRINT*,'push_to_stack not initialised'
    ENDIF

    IF (stack%stack_point .GT. stack%stack_size) THEN
      old_size = stack%stack_size
      stack%stack_size = 2 * stack%stack_size
      old_buffer => stack%entries
      ALLOCATE(stack%entries(stack%stack_size))
      stack%entries(1:old_size) = old_buffer(1:old_size)
      DEALLOCATE(old_buffer)
    ENDIF

    stack%entries(stack%stack_point) = value

  END SUBROUTINE push_to_stack



  SUBROUTINE pop_to_stack(stack1, stack2)

    TYPE(primitive_stack), INTENT(INOUT) :: stack1, stack2

    CALL push_to_stack(stack2, stack1%entries(stack1%stack_point))
    stack1%stack_point = stack1%stack_point - 1

  END SUBROUTINE pop_to_stack



  SUBROUTINE pop_to_null(stack)

    TYPE(primitive_stack), INTENT(INOUT) :: stack

    IF (.NOT.stack%init .AND. rank .EQ. 0) THEN
      PRINT*,'*** WARNING ***'
      PRINT*,'pop_to_null not initialised'
    ENDIF

    stack%stack_point = stack%stack_point - 1

  END SUBROUTINE pop_to_null



  SUBROUTINE pop_from_stack(stack, value)

    TYPE(stack_element), INTENT(OUT) :: value
    TYPE(primitive_stack), INTENT(INOUT) :: stack

    IF (.NOT.stack%init .AND. rank .EQ. 0) THEN
      PRINT*,'*** WARNING ***'
      PRINT*,'pop_from_stack not initialised'
    ENDIF

    value = stack%entries(stack%stack_point)
    stack%stack_point = stack%stack_point - 1

  END SUBROUTINE pop_from_stack



  SUBROUTINE stack_snoop(stack, value, offset)

    TYPE(primitive_stack), INTENT(INOUT) :: stack
    TYPE(stack_element), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: offset

    IF (stack%stack_point-offset .LE. 0) THEN
      PRINT *, 'Unable to snoop stack', stack%stack_point
      STOP
    ENDIF

    value = stack%entries(stack%stack_point-offset)

  END SUBROUTINE stack_snoop



  SUBROUTINE tokenize(expression, output, err)

    CHARACTER(LEN=*), INTENT(IN) :: expression
    TYPE(primitive_stack), INTENT(INOUT) :: output
    INTEGER, INTENT(INOUT) :: err

    CHARACTER(LEN=500) :: current
    INTEGER :: current_type, current_pointer, i, ptype

    TYPE(primitive_stack) :: stack
    TYPE(stack_element) :: block

    CALL initialise_stack(stack)

    current(:) = ' '
    current(1:1) = expression(1:1)
    current_pointer = 2
    current_type = char_type(expression(1:1))

    err = c_err_none

    last_block_type = c_pt_null

    DO i = 2, LEN(TRIM(expression))
      ptype = char_type(expression(i:i))
      ! This is a bit of a hack.
      ! Allow numbers to follow letters in an expression *except* in the
      ! special case of a single 'e' character, to allow 10.0e5, etc.
      IF (ptype .EQ. current_type .AND. .NOT. (ptype .EQ. c_char_delimiter) &
        .OR. (ptype .EQ. c_char_numeric .AND. current_type .EQ. c_char_alpha &
        .AND. .NOT. str_cmp(current, 'e'))) THEN
        current(current_pointer:current_pointer) = expression(i:i)
        current_pointer = current_pointer+1
      ELSE
#ifndef RPN_DECK
        CALL tokenize_subexpression_infix(current, block, stack, output, err)
#else
        CALL tokenize_subexpression_rpn(current, block, stack, output, err)
#endif
        IF (err .NE. c_err_none) RETURN
        current(:) = ' '
        current_pointer = 2
        current(1:1) = expression(i:i)
        current_type = ptype
      ENDIF
    ENDDO

#ifndef RPN_DECK
    CALL tokenize_subexpression_infix(current, block, stack, output, err)
#else
    CALL tokenize_subexpression_rpn(current, block, stack, output, err)
#endif
    IF (err .NE. c_err_none) RETURN

    DO i = 1, stack%stack_point
      CALL pop_to_stack(stack, output)
    ENDDO
    CALL deallocate_stack(stack)

  END SUBROUTINE tokenize



  SUBROUTINE tokenize_subexpression_infix(current, block, stack, output, err)

    ! This subroutine tokenizes input in normal infix maths notation
    ! It uses Dijkstra's shunting yard algorithm to convert to RPN

    CHARACTER(LEN=*), INTENT(IN) :: current
    TYPE(stack_element), INTENT(INOUT) :: block
    TYPE(primitive_stack), INTENT(INOUT) :: stack, output
    INTEGER, INTENT(INOUT) :: err
    TYPE(deck_constant) :: const
    TYPE(stack_element) :: block2
    INTEGER :: ipoint

    IF (ICHAR(current(1:1)) .EQ. 0) RETURN

    ! Populate the block
    CALL load_block(current, block)
#ifdef PARSER_DEBUG
    block%text = TRIM(current)
#endif
    IF (block%ptype .EQ. c_pt_bad) THEN
      IF (rank .EQ. 0) THEN
        PRINT *, 'Unable to parse block with text ', TRIM(current)
        CALL check_deprecated(current)
      ENDIF
      err = c_err_bad_value
      CALL deallocate_stack(stack)
      RETURN
    ENDIF

    IF (block%ptype .EQ. c_pt_deck_constant) THEN
      const = deck_constant_list(block%value)
      DO ipoint = 1, const%execution_stream%stack_point
        CALL push_to_stack(output, const%execution_stream%entries(ipoint))
      ENDDO
    ENDIF

    IF (block%ptype .NE. c_pt_parenthesis &
        .AND. block%ptype .NE. c_pt_null) THEN
      last_block_type = block%ptype
    ENDIF

    IF (block%ptype .EQ. c_pt_variable &
        .OR. block%ptype .EQ. c_pt_constant &
        .OR. block%ptype .EQ. c_pt_species &
        .OR. block%ptype .EQ. c_pt_subset) THEN
      CALL push_to_stack(output, block)
    ENDIF

    IF (block%ptype .EQ. c_pt_parenthesis) THEN
      IF (block%value .EQ. c_paren_left_bracket) THEN
        CALL push_to_stack(stack, block)
      ELSE
        DO
          CALL stack_snoop(stack, block2, 0)
          IF (block2%ptype .EQ. c_pt_parenthesis &
              .AND. block2%value .EQ. c_paren_left_bracket) THEN
            CALL pop_to_null(stack)
            ! If stack isn't empty then check for function
            IF (stack%stack_point .NE. 0) THEN
              CALL stack_snoop(stack, block2, 0)
              IF (block2%ptype .EQ. c_pt_function) THEN
                CALL add_function_to_stack(block2%value, stack, output)
              ENDIF
            ENDIF
            EXIT
          ELSE
            CALL pop_to_stack(stack, output)
          ENDIF
        ENDDO
      ENDIF
    ENDIF

    IF (block%ptype .EQ. c_pt_function) THEN
      ! Just push functions straight onto the stack
      CALL push_to_stack(stack, block)
    ENDIF

    IF (block%ptype .EQ. c_pt_separator) THEN
      DO
        CALL stack_snoop(stack, block2, 0)
        IF (block2%ptype .NE. c_pt_parenthesis) THEN
          CALL pop_to_stack(stack, output)
        ELSE
          IF (block2%value .NE. c_paren_left_bracket) THEN
            PRINT *, 'Bad function expression'
            STOP
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

    IF (block%ptype .EQ. c_pt_operator) THEN
      DO
        IF (stack%stack_point .EQ. 0) THEN
          ! stack is empty, so just push operator onto stack and
          ! leave loop
          CALL push_to_stack(stack, block)
          EXIT
        ENDIF
        ! stack is not empty so check precedence etc.
        CALL stack_snoop(stack, block2, 0)
        IF (block2%ptype .NE. c_pt_operator) THEN
          ! Previous block is not an operator so push current operator
          ! to stack and leave loop
          CALL push_to_stack(stack, block)
          EXIT
        ELSE
          IF (opcode_assoc(block%value) .EQ. c_assoc_la &
              .OR. opcode_assoc(block%value) .EQ. c_assoc_a) THEN
            ! Operator is full associative or left associative
            IF (opcode_precedence(block%value) &
                .LE. opcode_precedence(block2%value)) THEN
              CALL pop_to_stack(stack, output)
              CYCLE
            ELSE
              CALL push_to_stack(stack, block)
              EXIT
            ENDIF
          ELSE
            IF (opcode_precedence(block%value) &
                .LT. opcode_precedence(block2%value)) THEN
              CALL pop_to_stack(stack, output)
              CYCLE
            ELSE
              CALL push_to_stack(stack, block)
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE tokenize_subexpression_infix



  SUBROUTINE tokenize_subexpression_rpn(current, block, stack, output, err)

    ! This routine tokenizes input which is already in Reverse Polish Notiation

    CHARACTER(LEN=*), INTENT(IN) :: current
    TYPE(stack_element), INTENT(INOUT) :: block
    TYPE(primitive_stack), INTENT(INOUT) :: stack, output
    INTEGER, INTENT(INOUT) :: err
    TYPE(deck_constant) :: const
    TYPE(stack_element) :: block2
    INTEGER :: ipoint

    IF (ICHAR(current(1:1)) .EQ. 0) RETURN

    ! Populate the block
    CALL load_block(current, block)
#ifdef PARSER_DEBUG
    block%text = TRIM(current)
    PRINT *, block%ptype, TRIM(current)
#endif
    IF (block%ptype .EQ. c_pt_bad) THEN
      IF (rank .EQ. 0) THEN
        PRINT *, 'Unable to parse block with text ', TRIM(current)
      ENDIF
      err = c_err_bad_value
      CALL deallocate_stack(stack)
      RETURN
    ENDIF

    IF (block%ptype .NE. c_pt_parenthesis &
        .AND. block%ptype .NE. c_pt_null) THEN
      last_block_type = block%ptype
      IF (debug_mode) PRINT *, 'Setting', block%ptype, TRIM(current)
    ENDIF

    IF (block%ptype .EQ. c_pt_deck_constant) THEN
      const = deck_constant_list(block%value)
      DO ipoint = 1, const%execution_stream%stack_point
        CALL push_to_stack(output, const%execution_stream%entries(ipoint))
      ENDDO
    ELSE IF (block%ptype .NE. c_pt_null) THEN
      CALL push_to_stack(output, block)
    ENDIF

  END SUBROUTINE tokenize_subexpression_rpn



  SUBROUTINE add_function_to_stack(opcode, stack, output)

    INTEGER, INTENT(IN) :: opcode
    TYPE(primitive_stack), INTENT(INOUT) :: stack, output
    TYPE(primitive_stack) :: func_stack
    TYPE(stack_element) :: block
    INTEGER :: id, i, err

    CALL stack_snoop(output, block, 0)
    id = block%value
    func_stack%stack_point = 0

    IF (opcode .EQ. c_func_rho) THEN
      CALL copy_stack(species_list(id)%density_function, func_stack)
    ELSE IF (opcode .EQ. c_func_tempx) THEN
      CALL copy_stack(species_list(id)%temperature_function(1), func_stack)
    ELSE IF (opcode .EQ. c_func_tempy) THEN
      CALL copy_stack(species_list(id)%temperature_function(2), func_stack)
    ELSE IF (opcode .EQ. c_func_tempz) THEN
      CALL copy_stack(species_list(id)%temperature_function(3), func_stack)
    ELSE IF (opcode .EQ. c_func_tempx_ev) THEN
      CALL copy_stack(species_list(id)%temperature_function(1), func_stack)
      CALL tokenize('* kb / ev', func_stack, err)
    ELSE IF (opcode .EQ. c_func_tempy_ev) THEN
      CALL copy_stack(species_list(id)%temperature_function(2), func_stack)
      CALL tokenize('* kb / ev', func_stack, err)
    ELSE IF (opcode .EQ. c_func_tempz_ev) THEN
      CALL copy_stack(species_list(id)%temperature_function(3), func_stack)
      CALL tokenize('* kb / ev', func_stack, err)
    ENDIF

    IF (func_stack%stack_point .GT. 0) THEN
      CALL pop_to_null(output)
      CALL pop_to_null(stack)
      DO i = 1, func_stack%stack_point
        CALL push_to_stack(output, func_stack%entries(i))
      ENDDO
      CALL deallocate_stack(func_stack)
    ELSE
      CALL pop_to_stack(stack, output)
    ENDIF

  END SUBROUTINE add_function_to_stack



  SUBROUTINE display_tokens(token_list)

    TYPE(primitive_stack), INTENT(IN) :: token_list
    INTEGER :: i

    IF (rank .EQ. 0) THEN
      DO i = 1, token_list%stack_point
        PRINT *, 'Type', token_list%entries(i)%ptype
        PRINT *, 'Data', token_list%entries(i)%value
        PRINT *, 'NumData', token_list%entries(i)%numerical_data
#ifdef PARSER_DEBUG
        PRINT *, 'Text :', TRIM(token_list%entries(i)%text)
#endif
        PRINT *, '---------------'
      ENDDO
    ENDIF

  END SUBROUTINE display_tokens



  SUBROUTINE set_stack_zero(stack)

    TYPE(primitive_stack), INTENT(INOUT) :: stack

    CALL deallocate_stack(stack)
    CALL initialise_stack(stack)
    CALL tokenize('0', stack, errcode)

  END SUBROUTINE set_stack_zero

END MODULE shunt
