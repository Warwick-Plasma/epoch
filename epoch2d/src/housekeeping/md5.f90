MODULE md5

   IMPLICIT NONE
   INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9
   INTEGER(i4) :: h0, h1, h2, h3
   INTEGER :: total_length

CONTAINS

  SUBROUTINE lowercase(string)

    CHARACTER(LEN=*), INTENT(INOUT) :: string
    INTEGER :: i, c

    DO i = 1, LEN(string)
      c = ICHAR(string(i:i))
      IF (c > 64 .AND. c < 90) string(i:i) = CHAR(c+32)
    ENDDO

  END SUBROUTINE lowercase



  SUBROUTINE md5_init()

    h0 = 1732584193  ! 0x67452301
    h1 = -271733879  ! 0xefcdab89
    h2 = -1732584194 ! 0x98badcfe
    h3 = 271733878   ! 0x10325476
    total_length = 0

  END SUBROUTINE md5_init



  FUNCTION md5_generate(string) RESULT(md5)

    CHARACTER(LEN=32) :: md5
    CHARACTER(LEN=*), INTENT(IN) :: string

    CALL md5_init

    CALL md5_calculate(string, md5, .TRUE.)

  END FUNCTION md5_generate



  FUNCTION md5_append(string) RESULT(md5)

    CHARACTER(LEN=32) :: md5
    CHARACTER(LEN=*), INTENT(IN) :: string

    CALL md5_calculate(string, md5, .FALSE.)

  END FUNCTION md5_append



  SUBROUTINE md5_calculate(string, md5, generate)

    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=32), INTENT(INOUT) :: md5
    LOGICAL, INTENT(IN) :: generate
    CHARACTER(LEN=LEN(string)+64) :: new_string
    CHARACTER(LEN=8) :: wtmp
    LOGICAL :: pad_string

    INTEGER(i4) :: j, n1, n2, n3, n4, pos, new_len
    INTEGER(i4) :: a, b, c, d, f, g, temp
    INTEGER(i4) :: w(16), i, int_len
    INTEGER(i4), PARAMETER :: r(64) = (/ &
          7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,  5, &
          9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,  4, 11, &
         16, 23,  4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,  6, 10, 15, &
         21,  6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21 /)
    INTEGER(i4), PARAMETER :: k(64) = (/ &
           -680876936,  -389564586,   606105819, -1044525330,  -176418897, &
           1200080426, -1473231341,   -45705983,  1770035416, -1958414417, &
               -42063, -1990404162,  1804603682,   -40341101, -1502002290, &
           1236535329,  -165796510, -1069501632,   643717713,  -373897302, &
           -701558691,    38016083,  -660478335,  -405537848,   568446438, &
          -1019803690,  -187363961,  1163531501, -1444681467,   -51403784, &
           1735328473, -1926607734,     -378558, -2022574463,  1839030562, &
            -35309556, -1530992060,  1272893353,  -155497632, -1094730640, &
            681279174,  -358537222,  -722521979,    76029189,  -640364487, &
           -421815835,   530742520,  -995338651,  -198630844,  1126891415, &
          -1416354905,   -57434055,  1700485571, -1894986606,    -1051523, &
          -2054922799,  1873313359,   -30611744, -1560198380,  1309151649, &
           -145523070, -1120210379,   718787259,  -343485551 /)

    new_len = LEN(string)
    new_string(:new_len) = string

    IF (generate) THEN
      pad_string = .TRUE.
    ELSE
      pad_string = (new_len == 0 .OR. MOD(new_len, 64) /= 0)
    ENDIF

    total_length = total_length + new_len

    IF (pad_string) THEN
      j = new_len + 1
      new_string(j:j) = CHAR(128) ! Append "1" bit to message
      i = MOD(j, 64)
      DO WHILE(i /= 56)
        j = j + 1
        ! Pad with "0" until message length mod 512 bits equals 448 bits
        new_string(j:j) = CHAR(0)
        i = MOD(j, 64)
      ENDDO

      int_len = total_length * 8

      DO i = 0, 3
        temp = IAND(int_len, 255)
        j = j + 1
        new_string(j:j) = CHAR(temp)
        int_len = ISHFT(int_len, -8)
      ENDDO

      DO i = 1, 4
        j = j + 1
        new_string(j:j) = CHAR(0)
      ENDDO
      new_len = j
    ENDIF

    DO i = 1, INT(new_len/64)
      DO j = 1, 16
        pos = (j-1)*4+(i-1)*64
        n1 = ICHAR(new_string(4+pos:4+pos))
        n2 = ICHAR(new_string(3+pos:3+pos))
        n3 = ICHAR(new_string(2+pos:2+pos))
        n4 = ICHAR(new_string(1+pos:1+pos))

        WRITE(wtmp,'(4(z2.2))') n1, n2, n3, n4
        READ(wtmp,'(z8)') w(j)
      ENDDO

      a = h0
      b = h1
      c = h2
      d = h3

      DO j = 1, 64
        IF (j >= 1 .AND. j <= 16) THEN
          f = IOR(IAND(b, c), IAND(NOT(b), d))
          g = j
        ELSE IF (j >= 17 .AND. j <= 32) THEN
          f = IOR(IAND(d, b), IAND(NOT(d), c))
          g = MOD(5*(j-1) + 1, 16) + 1
        ELSE IF (j >= 33 .AND. j <= 48) THEN
          f = IEOR(b, IEOR(c, d))
          g = MOD(3*(j-1) + 5, 16) + 1
        ELSE IF (j >= 49 .AND. j <= 64) THEN
          f = IEOR(c, IOR(b, NOT(d)))
          g = MOD(7*(j-1), 16) + 1
        ENDIF

        temp = d
        d = c
        c = b
        b = b + leftrotate((a + f + k(j) + w(g)), r(j))
        a = temp
      ENDDO

      h0 = h0 + a
      h1 = h1 + b
      h2 = h2 + c
      h3 = h3 + d
    ENDDO

    a = to_bytes(h0)
    b = to_bytes(h1)
    c = to_bytes(h2)
    d = to_bytes(h3)

    WRITE(md5,'(4(z8.8))') a, b, c, d
    CALL lowercase(md5)

    IF (generate) THEN
      h0 = a
      h1 = b
      h2 = c
      h3 = d
    ENDIF

  END SUBROUTINE md5_calculate



  FUNCTION leftrotate(x, c)

    INTEGER(i4) :: leftrotate
    INTEGER(i4), INTENT(IN) :: x, c

    leftrotate = IOR(ISHFT(x, c), ISHFT(x, (c-32)))

  END FUNCTION leftrotate



  FUNCTION to_bytes(bytes)

    INTEGER(i4) :: to_bytes
    INTEGER(i4), INTENT(IN) :: bytes
    INTEGER(i4) :: i, tmp, shiftbytes

    to_bytes = 0
    shiftbytes = bytes
    DO i = 1, 4
      to_bytes = ISHFT(to_bytes, 8)
      tmp = IAND(shiftbytes, 255)
      to_bytes = to_bytes + tmp
      shiftbytes = ISHFT(shiftbytes, -8)
    ENDDO

  END FUNCTION to_bytes

END MODULE md5
