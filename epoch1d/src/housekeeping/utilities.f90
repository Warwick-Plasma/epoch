MODULE utilities

  USE constants

  IMPLICIT NONE

  INTERFACE grow_array
    MODULE PROCEDURE grow_real_array, grow_integer_array, grow_string_array
  END INTERFACE grow_array

  PRIVATE :: grow_real_array, grow_integer_array, grow_string_array

CONTAINS

  SUBROUTINE grow_real_array(array, idx)

    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    REAL(num), DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx .LE. old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    ENDDO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    ENDDO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_real_array



  SUBROUTINE grow_integer_array(array, idx)

    INTEGER, DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx .LE. old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    ENDDO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    ENDDO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_integer_array



  SUBROUTINE grow_string_array(array, idx)

    CHARACTER(LEN=string_length), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: idx
    CHARACTER(LEN=string_length), DIMENSION(:), ALLOCATABLE :: tmp_array
    INTEGER :: old_size, new_size, i

    old_size = SIZE(array)
    IF (idx .LE. old_size) RETURN

    ALLOCATE(tmp_array(old_size))
    DO i = 1, old_size
      tmp_array(i) = array(i)
    ENDDO

    new_size = 2 * old_size
    DEALLOCATE(array)
    ALLOCATE(array(new_size))

    DO i = 1, old_size
      array(i) = tmp_array(i)
    ENDDO

    DEALLOCATE(tmp_array)

  END SUBROUTINE grow_string_array

END MODULE utilities
