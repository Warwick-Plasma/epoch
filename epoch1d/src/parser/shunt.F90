MODULE shunt

  USE tokenizer_blocks

  IMPLICIT NONE

  SAVE

  TYPE(primitive_stack) :: dumper

CONTAINS

  FUNCTION char_type(char)

    CHARACTER, INTENT(IN) :: char
    INTEGER, PARAMETER :: c_nops = 7
    CHARACTER(LEN=c_nops) :: operators = "+-\/*^"
    INTEGER :: char_type, i
    LOGICAL :: unary

    char_type = c_char_unknown

    unary = (last_block_type .NE. c_pt_variable .AND. &
        last_block_type .NE. c_pt_constant)

    IF (char .EQ. " " .OR. ICHAR(char) .EQ. 32) char_type = c_char_space
    IF (char .GE. "A" .AND. char .LE. "z" .OR. char .EQ. "_") &
        char_type = c_char_alpha
    IF (char .EQ. "(" .OR. char .EQ. ")" .OR. char .EQ. ",") &
        char_type = c_char_delimiter
    DO i = 1, c_nops
      IF (char .EQ. operators(i:i)) char_type = c_char_opcode
    ENDDO

  END FUNCTION char_type



  SUBROUTINE load_block(name, BLOCK)

    CHARACTER(LEN=*), INTENT(IN) :: name
    TYPE(stack_element), INTENT(OUT) :: BLOCK
    INTEGER :: work
    REAL(num) :: value

    block%ptype = c_pt_bad
    block%data = 0
    block%numerical_data = 0.0_num
    work = 0

    IF (LEN(TRIM(name)) .EQ. 0) THEN
      block%ptype = c_pt_null
      block%data = 0
      block%numerical_data = 0.0_num
      RETURN
    ENDIF

    work = as_constant(name)
    IF (work .NE. 0) THEN
      ! block is a named constant
      block%ptype = c_pt_constant
      block%data = work
      RETURN
    ENDIF

    work = as_deferred_execution_object(name)
    IF (work .NE. 0) THEN
      ! block is a deferred execution object
      block%ptype = c_pt_deferred_execution_object
      block%data = work
      RETURN
    ENDIF

    work = as_operator(name)
    IF (work .NE. 0) THEN
      ! block is an operator
      block%ptype = c_pt_operator
      block%data = work
      RETURN
    ENDIF

    work = as_function(name)
    IF (work .NE. 0) THEN
      ! block is a function
      block%ptype = c_pt_function
      block%data = work
      RETURN
    ENDIF

    work = as_parenthesis(name)
    IF (work .NE. 0) THEN
      ! block is a parenthesis
      block%ptype = c_pt_parenthesis
      block%data = work
      RETURN
    ENDIF

    IF (str_cmp(name, ", ")) THEN
      block%ptype = c_pt_separator
      block%data = 0
      block%numerical_data = 0.0_num
      RETURN
    ENDIF

    value = as_real_simple(name, work)
    IF (IAND(work, c_err_bad_value) .EQ. 0) THEN
      ! block is a simple variable
      block%ptype = c_pt_variable
      block%data = 0
      block%numerical_data = value
    ENDIF

  END SUBROUTINE load_block



  FUNCTION as_parenthesis(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: as_parenthesis

    as_parenthesis = 0

    IF (str_cmp(name, "(")) THEN
      as_parenthesis = c_paren_left_bracket
    ENDIF

    IF (str_cmp(name, ")")) THEN
      as_parenthesis = c_paren_right_bracket
    ENDIF

  END FUNCTION as_parenthesis



  SUBROUTINE push_to_stack(stack, value)

    TYPE(stack_element), INTENT(IN) :: value
    TYPE(primitive_stack), INTENT(INOUT) :: stack

    stack%stack_point = stack%stack_point+1
    stack%data(stack%stack_point) = value

  END SUBROUTINE push_to_stack



  SUBROUTINE pop_to_stack(stack1, stack2)

    TYPE(primitive_stack), INTENT(INOUT) :: stack1, stack2

    stack2%stack_point = stack2%stack_point+1
    stack2%data(stack2%stack_point) = stack1%data(stack1%stack_point)
    stack1%stack_point = stack1%stack_point-1

  END SUBROUTINE pop_to_stack



  SUBROUTINE pop_to_null(stack)

    TYPE(primitive_stack), INTENT(INOUT) :: stack

    stack%stack_point = stack%stack_point-1

  END SUBROUTINE pop_to_null



  SUBROUTINE pop_from_stack(stack, value)

    TYPE(stack_element), INTENT(OUT) :: value
    TYPE(primitive_stack), INTENT(INOUT) :: stack

    value = stack%data(stack%stack_point)
    stack%stack_point = stack%stack_point-1

  END SUBROUTINE pop_from_stack



  SUBROUTINE stack_snoop(stack, value, offset)

    TYPE(primitive_stack), INTENT(INOUT) :: stack
    TYPE(stack_element), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: offset

    IF (stack%stack_point-offset .LE. 0) THEN
      PRINT *, "Unable to snoop stack", stack%stack_point
      STOP
    ENDIF

    value = stack%data(stack%stack_point-offset)

  END SUBROUTINE stack_snoop



  SUBROUTINE tokenize(expression, output, err)

    CHARACTER(LEN=*), INTENT(IN) :: expression
    TYPE(primitive_stack), INTENT(INOUT) :: output
    INTEGER, INTENT(INOUT) :: err

#ifndef RPN_DECK
    CALL tokenize_infix(expression, output, err)
#else
    CALL tokenize_rpn(expression, output, err)
#endif

  END SUBROUTINE tokenize



  SUBROUTINE tokenize_infix(expression, output, err)

    ! This subroutine tokenizes input in normal infix maths notation
    ! It uses Dijkstra's shunting yard algorithm to convert to RPN

    CHARACTER(LEN=*), INTENT(IN) :: expression
    TYPE(primitive_stack), INTENT(INOUT) :: output
    INTEGER, INTENT(INOUT) :: err
    TYPE(deferred_execution_object) :: deo

    CHARACTER(LEN=500) :: current
    INTEGER :: current_type, current_pointer, i, ptype, ipoint

    TYPE(primitive_stack) :: stack
    TYPE(stack_element) :: BLOCK, block2

    stack%stack_point = 0

    current(:) = " "
    current(1:1) = expression(1:1)
    current_pointer = 2
    current_type = char_type(expression(1:1))

    err = c_err_none

    last_block_type = c_pt_null

    DO i = 2, LEN(expression)
      ptype = char_type(expression(i:i))
      IF (ptype .EQ. current_type .AND. &
          .NOT. (ptype .EQ. c_char_delimiter)) THEN
        current(current_pointer:current_pointer) = expression(i:i)
        current_pointer = current_pointer+1
      ELSE
        IF (ICHAR(current(1:1)) .NE. 0) THEN
          ! Populate the block
          CALL load_block(current, BLOCK)
!!$#ifdef PARSER_DEBUG
!!$          block%text = TRIM(current)
!!$#endif
          IF (block%ptype .EQ. c_pt_bad) THEN
            IF (rank .EQ. 0) THEN
              PRINT *, "Unable to parse block with text ", TRIM(current)
            ENDIF
            err = c_err_bad_value
            RETURN
          ENDIF
          IF (block%ptype .EQ. c_pt_deferred_execution_object) THEN
            deo = deferred_objects(block%data)
            DO ipoint = 1, deo%execution_stream%stack_point
              CALL push_to_stack(output, deo%execution_stream%data(ipoint))
            ENDDO
          ENDIF

          IF (block%ptype .NE. c_pt_parenthesis .AND. &
              block%ptype .NE. c_pt_null) THEN
            last_block_type = block%ptype
          ENDIF

          IF (block%ptype .EQ. c_pt_variable .OR. &
              block%ptype .EQ. c_pt_constant) THEN
            CALL push_to_stack(output, BLOCK)
          ENDIF

          IF (block%ptype .EQ. c_pt_parenthesis) THEN
            IF (block%data .EQ. c_paren_left_bracket) THEN
              CALL push_to_stack(stack, BLOCK)
            ELSE
              DO
                CALL stack_snoop(stack, block2, 0)
                IF (block2%ptype .EQ. c_pt_parenthesis .AND. &
                    block2%data .EQ. c_paren_left_bracket) THEN
                  CALL pop_to_null(stack)
                  ! If stack isn't empty then check for function
                  IF (stack%stack_point .NE. 0) THEN
                    CALL stack_snoop(stack, block2, 0)
                    IF (block2%ptype .EQ. c_pt_function) &
                        CALL pop_to_stack(stack, output)
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
            CALL push_to_stack(stack, BLOCK)
          ENDIF

          IF (block%ptype .EQ. c_pt_separator) THEN
            DO
              CALL stack_snoop(stack, block2, 0)
              IF (block2%ptype .NE. c_pt_parenthesis) THEN
                CALL pop_to_stack(stack, output)
              ELSE
                IF (block2%data .NE. c_paren_left_bracket) THEN
                  PRINT *, "Bad function expression"
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
                CALL push_to_stack(stack, BLOCK)
                EXIT
              ENDIF
              ! stack is not empty so check precedence etc.
              CALL stack_snoop(stack, block2, 0)
              IF (block2%ptype .NE. c_pt_operator) THEN
                ! Previous block is not an operator so push current operator
                ! to stack and leave loop
                CALL push_to_stack(stack, BLOCK)
                EXIT
              ELSE
                IF (opcode_assoc(block%data) .EQ. c_assoc_la .OR. &
                    opcode_assoc(block%data) .EQ. c_assoc_a) THEN
                  ! Operator is full associative or left associative
                  IF (opcode_precedence(block%data) .LE. &
                      opcode_precedence(block2%data)) THEN
                    CALL pop_to_stack(stack, output)
                    CYCLE
                  ELSE
                    CALL push_to_stack(stack, BLOCK)
                    EXIT
                  ENDIF
                ELSE
                  IF (opcode_precedence(block%data) .LT. &
                      opcode_precedence(block2%data)) THEN
                    CALL pop_to_stack(stack, output)
                    CYCLE
                  ELSE
                    CALL push_to_stack(stack, BLOCK)
                    EXIT
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF

        ENDIF
        current(:) = " "
        current_pointer = 2
        current(1:1) = expression(i:i)
        current_type = ptype
      ENDIF
    ENDDO

    DO i = 1, stack%stack_point
      CALL pop_to_stack(stack, output)
    ENDDO

  END SUBROUTINE tokenize_infix



  SUBROUTINE tokenize_rpn(expression, output, err)

    ! This routine tokenizes input which is already in Reverse Polish Notiation

    CHARACTER(LEN=*), INTENT(IN) :: expression
    TYPE(primitive_stack), INTENT(INOUT) :: output
    TYPE(deferred_execution_object) :: deo
    INTEGER, INTENT(INOUT) :: err

    CHARACTER(LEN=500) :: current
    INTEGER :: current_type, current_pointer, i, ptype, ipoint

    TYPE(primitive_stack) :: stack
    TYPE(stack_element) :: BLOCK

    stack%stack_point = 0
    last_block_type = c_pt_null

    current(:) = " "
    current(1:1) = expression(1:1)
    current_pointer = 2
    current_type = char_type(expression(1:1))

    err = c_err_none

    DO i = 2, LEN(expression)
      ptype = char_type(expression(i:i))
      IF (ptype .EQ. current_type .AND. &
          .NOT. (ptype .EQ. c_char_delimiter)) THEN
        current(current_pointer:current_pointer) = expression(i:i)
        current_pointer = current_pointer+1
      ELSE
        IF (ICHAR(current(1:1)) .NE. 0) THEN
          ! Populate the block
          CALL load_block(current, BLOCK)
!!$#ifdef PARSER_DEBUG
!!$          block%text = TRIM(current)
!!$#endif
          PRINT *, block%ptype, TRIM(current)
          IF (block%ptype .EQ. c_pt_bad) THEN
            IF (rank .EQ. 0) THEN
              PRINT *, "Unable to parse block with text ", TRIM(current)
            ENDIF
            err = c_err_bad_value
            RETURN
          ENDIF
          IF (block%ptype .NE. c_pt_parenthesis .AND. &
              block%ptype .NE. c_pt_null) THEN
            last_block_type = block%ptype
            IF (debug_mode) PRINT *, "Setting", block%ptype, TRIM(current)
          ENDIF
          IF (block%ptype .EQ. c_pt_deferred_execution_object) THEN
            deo = deferred_objects(block%data)
            DO ipoint = 1, deo%execution_stream%stack_point
              CALL push_to_stack(output, deo%execution_stream%data(ipoint))
            ENDDO
            CYCLE
          ELSE IF (block%ptype .NE. c_pt_null) THEN
            CALL push_to_stack(output, BLOCK)
          ENDIF
        ENDIF
        current(:) = " "
        current_pointer = 2
        current(1:1) = expression(i:i)
        current_type = ptype
      ENDIF
    ENDDO

    DO i = 1, stack%stack_point
      CALL pop_to_stack(stack, output)
    ENDDO

  END SUBROUTINE tokenize_rpn



  SUBROUTINE display_tokens(token_list)

    TYPE(primitive_stack), INTENT(IN) :: token_list
    INTEGER :: i

    IF (rank .EQ. 0) THEN
      DO i = 1, token_list%stack_point
        PRINT *, "Type", token_list%data(i)%ptype
        PRINT *, "Data", token_list%data(i)%data
        PRINT *, "NumData", token_list%data(i)%numerical_data
!!$#ifdef PARSER_DEBUG
!!$        PRINT *, "Text :", TRIM(token_list%data(i)%text)
!!$#endif
        PRINT *, "---------------"
      ENDDO
    ENDIF

  END SUBROUTINE display_tokens

END MODULE shunt
