MODULE random_generator

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: random, random_init, get_random_state, set_random_state
  INTEGER :: x = 123456789, y = 362436069, z = 521288629, w = 916191069

CONTAINS

  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators,
  !     period 597273182964842497>2^59
  ! Overall period>2^123;  Default seeds x,y,z,w.
  ! Set your own seeds with random_init(ix,iy,iz,iw).
  FUNCTION random()

    DOUBLE PRECISION :: random
    INTEGER :: kiss, a1, b1, a2, b2, a3, b3

    x = 69069 * x + 1327217885

    a1 = y
    b1 = 13
    a2 = IEOR(a1, ISHFT(a1, b1))
    b2 = -17
    a3 = IEOR(a2, ISHFT(a2, b2))
    b3 = 5
    y  = IEOR(a3, ISHFT(a3, b3))

    z = 18000 * IAND(z, 65535) + ISHFT(z, - 16)
    w = 30903 * IAND(w, 65535) + ISHFT(w, - 16)

    kiss = x + y + ISHFT(z, 16) + w

    random = (DBLE(kiss) + 2147483648.d0) / 4294967296.d0

  END FUNCTION random



  SUBROUTINE random_init(seed)

    INTEGER, INTENT(IN) :: seed
    INTEGER :: i
    DOUBLE PRECISION :: dummy

    x = x + seed
    y = y + seed
    z = z + seed
    w = w + seed

    ! 'Warm-up' the generator by cycling through a few times
    DO i = 1, 1000
      dummy = random()
    ENDDO

  END SUBROUTINE random_init



  SUBROUTINE get_random_state(state)

    INTEGER, INTENT(OUT) :: state(:)

    state(1) = x
    state(2) = y
    state(3) = z
    state(4) = w

  END SUBROUTINE get_random_state



  SUBROUTINE set_random_state(state)

    INTEGER, INTENT(IN) :: state(:)

    x = state(1)
    y = state(2)
    z = state(3)
    w = state(4)

  END SUBROUTINE set_random_state

END MODULE random_generator
