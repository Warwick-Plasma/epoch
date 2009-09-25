MODULE Shunt

  USE shared_parser_data
  USE shared_data
  USE strings
  USE tokenizer_blocks
  IMPLICIT NONE

  SAVE

  TYPE(primitivestack):: dumper

CONTAINS


  FUNCTION CharKind(Char)

    CHARACTER,INTENT(IN) :: Char
    INTEGER,PARAMETER :: NOPS=7
    CHARACTER(len=NOPS) :: Operators="+-\/*^"
    INTEGER :: CharKind,i
    LOGICAL :: unary

    CharKind=CHAR_UNKNOWN

    unary=(last_block_type .NE. PT_VARIABLE .AND. last_block_type .NE. PT_CONSTANT)

    IF (Char .EQ. " " .OR. ICHAR(Char) .EQ. 32) CharKind=CHAR_SPACE
    IF (Char >= "A" .AND. Char <= "z" .OR. Char .EQ. "_") CharKind=CHAR_ALPHA
    IF (Char .EQ. "(" .OR. CHAR .EQ. ")" .OR. CHAR .EQ. ",") CharKind=CHAR_DELIMITER
    DO i=1,NOPS
       IF (Char .EQ. Operators(i:i)) CharKind=CHAR_OPCODE
    ENDDO
  END FUNCTION CharKind

  SUBROUTINE LoadBlock(name,block)

    CHARACTER(len=*), INTENT(IN) :: name
    TYPE(StackElement), INTENT(OUT) :: block
    INTEGER :: work
    REAL(num) :: Value

    Block%type=PT_BAD
    Block%Data=0
    Block%NumericalData=0.0_num
    work=0

    IF (LEN(TRIM(name)) .EQ. 0) THEN
       Block%Type=PT_NULL
       Block%Data=0
       Block%NumericalData=0.0_num
       RETURN
    ENDIF

    work=AsConstant(name)
    IF (work .NE. 0) THEN
       !Block is a named constant
       Block%Type=PT_CONSTANT
       Block%Data=work
       RETURN
    ENDIF

    work=AsDeferredExecutionObject(name)
    IF (work .NE. 0) THEN
       !Block is a deferred execution object
       Block%Type=PT_DEFERRED_EXECUTION_OBJECT
       Block%Data=work
       RETURN
    ENDIF

    work=AsOperator(name)
    IF (work .NE. 0) THEN
       !Block is an operator
       Block%Type=PT_OPERATOR
       Block%Data=work
       RETURN
    ENDIF

    work=AsFunction(name)
    IF (work .NE. 0) THEN
       !Block is a function
       Block%Type=PT_FUNCTION
       Block%Data=work 
       RETURN
    ENDIF

    work=AsParenthesis(name)
    IF (work .NE. 0) THEN
       !Block is a parenthesis
       Block%Type=PT_PARENTHESIS
       Block%Data=work
       RETURN
    ENDIF

    IF (StrCmp(name,",")) THEN
       Block%Type=PT_SEPARATOR
       Block%Data=0
       Block%NumericalData=0.0_num
       RETURN
    ENDIF

    Value=AsRealSimple(name,work)
    IF (IAND(work,ERR_BAD_VALUE) .EQ. 0) THEN
       !Block is a simple variable
       Block%Type=PT_VARIABLE
       Block%Data=0
       Block%NumericalData=Value
    ENDIF

  END SUBROUTINE LoadBlock

  FUNCTION AsParenthesis(name)
    CHARACTER(len=*), INTENT(IN) :: name
    INTEGER :: AsParenthesis

    AsParenthesis=0

    IF (StrCmp(name,"(")) THEN
       AsParenthesis=PAREN_LEFT_BRACKET
    ENDIF

    IF (StrCmp(name,")")) THEN
       AsParenthesis=PAREN_RIGHT_BRACKET
    ENDIF
  END FUNCTION AsParenthesis

  SUBROUTINE PushToStack(Stack,Value)

    TYPE(StackElement),INTENT(IN) :: Value
    TYPE(PrimitiveStack),INTENT(INOUT) :: Stack

    Stack%StackPoint=Stack%StackPoint+1
    Stack%Data(Stack%StackPoint)=Value   

  END SUBROUTINE PushToStack

  SUBROUTINE PopToStack(Stack1,Stack2)

    TYPE(PrimitiveStack),INTENT(INOUT) :: Stack1,Stack2

    Stack2%StackPoint=Stack2%StackPoint+1
    Stack2%Data(Stack2%StackPoint)=Stack1%Data(Stack1%StackPoint)
    Stack1%StackPoint=Stack1%StackPoint-1

  END SUBROUTINE PopToStack

  SUBROUTINE PopToNull(Stack)

    TYPE(PrimitiveStack),INTENT(INOUT) :: Stack

    Stack%StackPoint=Stack%StackPoint-1

  END SUBROUTINE PopToNull

  SUBROUTINE PopFromStack(Stack,Value)

    TYPE(StackElement),INTENT(OUT) :: Value
    TYPE(PrimitiveStack),INTENT(INOUT) :: Stack

    Value=Stack%Data(Stack%StackPoint)
    Stack%StackPoint=Stack%StackPoint-1

  END SUBROUTINE PopFromStack


  SUBROUTINE Stack_Snoop(stack,value,offset)

    TYPE(PrimitiveStack),INTENT(INOUT) :: Stack
    TYPE(StackElement),INTENT(OUT) :: Value
    INTEGER,INTENT(IN) :: offset

    IF (Stack%StackPoint-offset .LE. 0) THEN
       PRINT *,"Unable to snoop stack",Stack%stackpoint
       STOP
    ENDIF

    Value=Stack%Data(Stack%StackPoint-offset)

  END SUBROUTINE Stack_Snoop

  SUBROUTINE Tokenize(Expression,output,ERR)

    CHARACTER(len=*),INTENT(IN) :: Expression
    TYPE(primitivestack), INTENT(INOUT) :: output
    INTEGER,INTENT(INOUT) :: ERR

#ifndef RPN_DECK
    CALL Tokenize_Infix(Expression,output,ERR)
#else
    CALL Tokenize_RPN(Expression,output,ERR)
#endif

  END SUBROUTINE Tokenize

  SUBROUTINE Tokenize_Infix(Expression,output,ERR)

    !This subroutine tokenizes input in normal infix maths notation
    !It uses Dijkstra's shunting yard algorithm to convert to RPN

    CHARACTER(len=*),INTENT(IN) :: Expression
    TYPE(primitivestack), INTENT(INOUT) :: output
    INTEGER,INTENT(INOUT) :: ERR
    TYPE(Deferred_Execution_Object) :: DEO

    CHARACTER(len=500) :: current
    INTEGER :: Current_Type,CurrentPointer,i,Type,iPoint
    INTEGER :: OpCode,OpCode_Last
    REAL(num) :: val
    CHARACTER(len=6) :: Operators="+-*\^%"

    TYPE(PrimitiveStack) :: stack
    TYPE(StackElement) :: Block,Block2

    stack%StackPoint=0

    Current(:)=" "
    Current(1:1)=Expression(1:1)
    CurrentPointer=2
    Current_Type=CharKind(Expression(1:1))

    ERR=ERR_NONE

    last_block_type=PT_NULL

    DO i=2,LEN(Expression)
       Type=CharKind(Expression(i:i))
       IF (Type .EQ. Current_Type .AND. .NOT. (Type .EQ. CHAR_DELIMITER)) THEN
          Current(CurrentPointer:CurrentPointer) = Expression(i:i)
          CurrentPointer=CurrentPointer+1
       ELSE
          IF (ICHAR(Current(1:1)) .NE. 0) THEN
             !Populate the block
             CALL LoadBlock(Current,Block)
#ifdef PARSER_DEBUG
             Block%Text=TRIM(Current)
#endif
             IF (Block%Type .EQ. PT_BAD) THEN
                IF (rank .EQ. 0) THEN
                   PRINT *,"Unable to parse block with text ",TRIM(Current)
                ENDIF
                Err=ERR_BAD_VALUE
                RETURN
             ENDIF
             IF (Block%Type .EQ. PT_DEFERRED_EXECUTION_OBJECT) THEN
                DEO=Deferred_Objects(Block%Data)
                DO iPoint=1,DEO%Execution_Stream%stackpoint
                   CALL PushToStack(output,DEO%Execution_Stream%Data(iPoint))
                ENDDO
             ENDIF

             IF (Block%Type .NE. PT_PARENTHESIS .AND. Block%Type .NE. PT_NULL) THEN
                last_block_type=Block%Type
             ENDIF

             IF (Block%Type .EQ. PT_VARIABLE .OR. Block%Type .EQ. PT_CONSTANT) THEN
                CALL PushToStack(output,Block)
             ENDIF

             IF (Block%Type .EQ. PT_PARENTHESIS) THEN
                IF (Block%Data .EQ. PAREN_LEFT_BRACKET) THEN
                   CALL PushToStack(stack,block)
                ELSE
                   DO
                      CALL Stack_Snoop(stack,block2,0)
                      IF (block2%type .EQ. PT_PARENTHESIS .AND. block2%Data .EQ. PAREN_LEFT_BRACKET) THEN
                         CALL PopToNull(stack)
                         !If stack isn't empty then check for function
                         IF (stack%stackpoint .NE. 0) THEN
                            CALL Stack_Snoop(stack,block2,0)
                            IF (block2%Type .EQ. PT_FUNCTION) CALL PopToStack(stack,output)
                         ENDIF
                         EXIT
                      ELSE
                         CALL PopToStack(stack,output)
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF

             IF (Block%Type .EQ. PT_FUNCTION) THEN
                !Just push functions straight onto the stack
                CALL PushToStack(stack,block)
             ENDIF

             IF (Block%Type .EQ. PT_SEPARATOR) THEN
                DO
                   CALL Stack_Snoop(stack,block2,0)
                   IF (block2%Type .NE. PT_PARENTHESIS) THEN
                      CALL PopToStack(stack,output)
                   ELSE
                      IF (block2%Data .NE. PAREN_LEFT_BRACKET) THEN
                         PRINT *,"Bad function expression"
                         STOP
                      ENDIF
                      EXIT
                   ENDIF
                ENDDO
             ENDIF

             IF (Block%Type .EQ. PT_OPERATOR) THEN
                DO
                   IF(Stack%StackPoint .EQ. 0) THEN
                      !Stack is empty, so just push operator onto stack and leave loop
                      CALL PushToStack(stack,block)
                      EXIT
                   ENDIF
                   !Stack is not empty so check precedence etc.
                   CALL Stack_Snoop(stack,block2,0)
                   IF (block2%type .NE. PT_OPERATOR) THEN
                      !Previous block is not an operator so push current operator to stack and leave loop
                      CALL PushToStack(stack,block)
                      EXIT
                   ELSE
                      IF (OpCode_Assoc(Block%Data) .EQ. ASSOC_LA .OR. OpCode_Assoc(Block%Data) .EQ. ASSOC_A) THEN
                         !Operator is full associative or left associative
                         IF (OpCode_Precedence(Block%Data) .LE. OpCode_Precedence(Block2%Data)) THEN
                            CALL PopToStack(stack,output)
                            CYCLE
                         ELSE
                            CALL PushToStack(stack,block)
                            EXIT
                         ENDIF
                      ELSE
                         IF (OpCode_Precedence(Block%Data) .LT. OpCode_Precedence(Block2%Data)) THEN
                            CALL PopToStack(stack,output)
                            CYCLE
                         ELSE
                            CALL PushToStack(stack,block)
                            EXIT
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF

          ENDIF
          Current(:) =" "
          CurrentPointer=2
          Current(1:1) = Expression(i:i)
          Current_Type=Type
       ENDIF
    ENDDO

    DO i=1,Stack%StackPoint
       CALL PopToStack(Stack,output)
    ENDDO

  END SUBROUTINE Tokenize_Infix


  SUBROUTINE Tokenize_RPN(Expression,output,ERR)

    !This routine tokenizes input which is already in Reverse Polish Notiation

    CHARACTER(len=*),INTENT(IN) :: Expression
    TYPE(primitivestack), INTENT(INOUT) :: output
    TYPE(Deferred_Execution_Object) :: DEO
    INTEGER,INTENT(INOUT) :: ERR

    CHARACTER(len=500) :: current
    INTEGER :: Current_Type,CurrentPointer,i,Type,iPoint
    INTEGER :: OpCode,OpCode_Last
    REAL(num) :: val
    CHARACTER(len=6) :: Operators="+-*\^%"

    TYPE(PrimitiveStack) :: stack
    TYPE(StackElement) :: Block,Block2

    stack%StackPoint=0
    last_block_type=PT_NULL

    Current(:)=" "
    Current(1:1)=Expression(1:1)
    CurrentPointer=2
    Current_Type=CharKind(Expression(1:1))

    ERR=ERR_NONE


    DO i=2,LEN(Expression)
       Type=CharKind(Expression(i:i))
       IF (Type .EQ. Current_Type .AND. .NOT. (Type .EQ. CHAR_DELIMITER)) THEN
          Current(CurrentPointer:CurrentPointer) = Expression(i:i)
          CurrentPointer=CurrentPointer+1
       ELSE
          IF (ICHAR(Current(1:1)) .NE. 0) THEN
             !Populate the block
             CALL LoadBlock(Current,Block)
#ifdef PARSER_DEBUG
             Block%Text=TRIM(Current)
#endif
             PRINT *,Block%Type,TRIM(Current)
             IF (Block%Type .EQ. PT_BAD) THEN
                IF (rank .EQ. 0) THEN
                   PRINT *,"Unable to parse block with text ",TRIM(Current)
                ENDIF
                Err=ERR_BAD_VALUE
                RETURN
             ENDIF
             IF (Block%Type .NE. PT_PARENTHESIS .AND. Block%Type .NE. PT_NULL) THEN
                last_block_type=Block%Type
                IF (DEBUG_MODE) PRINT *,"Setting",Block%Type,TRIM(Current)
             ENDIF
             IF (Block%Type .EQ. PT_DEFERRED_EXECUTION_OBJECT) THEN
                DEO=Deferred_Objects(Block%Data)
                DO iPoint=1,DEO%Execution_Stream%stackpoint
                   CALL PushToStack(output,DEO%Execution_Stream%Data(iPoint))
                ENDDO
                CYCLE
             ELSE IF (Block%Type .NE. PT_NULL) THEN
                CALL PushToStack(output,Block)
             ENDIF
          ENDIF
          Current(:) =" "
          CurrentPointer=2
          Current(1:1) = Expression(i:i)
          Current_Type=Type
       ENDIF
    ENDDO

    DO i=1,Stack%StackPoint
       CALL PopToStack(Stack,output)
    ENDDO

  END SUBROUTINE Tokenize_RPN

  SUBROUTINE DisplayTokens(TokenList)
    TYPE(primitivestack), INTENT(IN) :: TokenList
    INTEGER :: i
    IF (rank .EQ. 0) THEN
       DO i=1,TokenList%StackPoint
          PRINT *,"Type",TokenList%Data(i)%Type
          PRINT *,"Data",TokenList%Data(i)%Data
          PRINT *,"NumData",TokenList%Data(i)%NumericalData
#ifdef PARSER_DEBUG
          PRINT *,"Text :",TRIM(TokenList%Data(i)%Text)
#endif
          PRINT *,"---------------"
       ENDDO
    ENDIF

  END SUBROUTINE DisplayTokens

END MODULE Shunt
