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

MODULE numerics

  USE constants

  IMPLICIT NONE

CONTAINS

  REAL(num) FUNCTION factorial(x)

    REAL(num), INTENT(IN) :: x
    INTEGER :: i

    factorial = 1.0_num
    DO i = 2, INT(x)
      IF (factorial < SQRT(HUGE(factorial))) THEN
        factorial = factorial * REAL(i, num)
      ELSE
        factorial = HUGE(0.0_num)
        EXIT
      END IF
    END DO
    RETURN

  END FUNCTION factorial



  REAL(num) FUNCTION gamma_fn(X)
!----------------------------------------------------------------------
!
! This routine calculates the GAMMA function for a real argument X.
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X >= 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.
!
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!          approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS > 1.0
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable
!
!     Approximate values for some important machines are:
!
!                            beta       maxexp        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.
!
!
!  Intrinsic functions required are:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.
!
!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.
!
!
!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!  Latest modification: December 1, 2011
!
!  Modified by: A. A. Lawrence-Douglas,
!               Department of Physics,
!               University of Warwick,
!               Coventry,
!               Warwickshire,
!               England
!               CV47 0HD
!
!----------------------------------------------------------------------

    REAL(num), INTENT(IN) :: x
    INTEGER :: i, n, iy
    LOGICAL :: parity, inv
    REAL(num) :: c(7), fact, p(8), q(8), res, sum, xden, xnum, y, y1, ysq, z

!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------

    REAL(num), PARAMETER :: sqrtpi = SQRT(pi)

!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------

    REAL(num), PARAMETER :: xbig = 171.624_num, xminin = 2.23e-308_num, &
        eps = EPSILON(1.0_num)

!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1, 2).
!----------------------------------------------------------------------

    P = (/ -1.71618513886549492533811_num,  24.7656508055759199108314_num, &
           -379.804256470945635097577_num,  629.331155312818442661052_num, &
            866.966202790413211295064_num, -31451.2729688483675254357_num, &
           -36144.4134186911729807069_num,  66456.1438202405440627855_num /)
    Q = (/ -30.8402300119738975254353_num,  315.350626979604161529144_num, &
           -1015.15636749021914166146_num, -3107.77167157231109440444_num, &
            22538.1184209801510330112_num,  4755.84627752788110767815_num, &
           -134659.959864969306392456_num, -115132.259675553483497211_num /)

!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------

    C = (/ -1.910444077728e-3_num,          8.4171387781295e-4_num, &
           -5.952379913043012e-4_num,       7.93650793500350248e-4_num, &
           -2.777777777777681622553e-3_num, 8.333333333333333331554247e-2_num, &
            5.7083835261e-3_num /)

!----------------------------------------------------------------------
!  Statement functions for conversion between integer and float
!----------------------------------------------------------------------

    parity = .FALSE.
    inv = .FALSE.
    fact = 1.0_num
    n = 0
    y = x
    IF (y <= 0.0_num) THEN
!----------------------------------------------------------------------
!  Argument is negative
!----------------------------------------------------------------------
      y = -x
      y1 = AINT(y)
      iy = FLOOR(y1)
      res = y - y1
      IF (res > c_tiny) THEN
        IF (iy /= (iy / 2) * 2) parity = .TRUE.
        inv = .TRUE.
        fact = -pi / SIN(pi * res)
        y = y + 1.0_num
      ELSE
        gamma_fn = c_largest_number
        RETURN
      END IF
    END IF
!----------------------------------------------------------------------
!  Argument is positive
!----------------------------------------------------------------------
    IF (y < eps) THEN
!----------------------------------------------------------------------
!  Argument < EPS
!----------------------------------------------------------------------
      IF (y >= xminin) THEN
        res = 1.0_num / y
      ELSE
        gamma_fn = c_largest_number
        RETURN
      END IF
    ELSE IF (y < 12.0_num) THEN
      y1 = y
      IF (y < 1.0_num) THEN
!----------------------------------------------------------------------
!  0.0 < argument < 1.0
!----------------------------------------------------------------------
        z = y
        y = y + 1.0_num
      ELSE
!----------------------------------------------------------------------
!  1.0 < argument < 12.0, reduce argument if necessary
!----------------------------------------------------------------------
        n = INT(y) - 1
        y = y - REAL(n, num)
        z = y - 1.0_num
      END IF
!----------------------------------------------------------------------
!  Evaluate approximation for 1.0 < argument < 2.0
!----------------------------------------------------------------------
      xnum = 0.0_num
      xden = 1.0_num
      DO i = 1, 8
         xnum = (xnum + p(i)) * z
         xden = xden * z + q(i)
      END DO
      res = xnum / xden + 1.0_num
      IF (y1 < y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  0.0 < argument < 1.0
!----------------------------------------------------------------------
        res = res / y1
      ELSE IF (y1 > y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  2.0 < argument < 12.0
!----------------------------------------------------------------------
        DO i = 1, n
          res = res * y
          y = y + 1.0_num
        END DO
      END IF
    ELSE
!----------------------------------------------------------------------
!  Evaluate for argument >= 12.0,
!----------------------------------------------------------------------
      IF (y <= xbig) THEN
        ysq = y * y
        sum = c(7)
        DO i = 1, 6
          sum = sum / ysq + c(i)
        END DO
        sum = sum / y - y + sqrtpi
        sum = sum + (y - 0.5_num) * LOG(y)
        res = EXP(sum)
      ELSE
        gamma_fn = c_largest_number
        RETURN
      END IF
    END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
    IF (parity) res = -res
    IF (inv) res = fact / res
    gamma_fn = res
! ---------- Last line of GAMMA ----------
  END FUNCTION gamma_fn



  REAL(num) FUNCTION rkbesl(x, alpha, nb, ize, ncalc)
!-------------------------------------------------------------------
!
!  This FORTRAN 77 routine calculates modified Bessel functions
!  of the second kind, K SUB(N+ALPHA) (X), for non-negative
!  argument X, and non-negative order N+ALPHA, with or without
!  exponential scaling.
!
!  Explanation of variables in the calling sequence
!
!  Description of output values ..
!
! X     - Working precision non-negative real argument for which
!         K's or exponentially scaled K's (K*EXP(X))
!         are to be calculated.  If K's are to be calculated,
!         X must not be greater than XMAX (see below).
! ALPHA - Working precision fractional part of order for which
!         K's or exponentially scaled K's (K*EXP(X)) are
!         to be calculated.  0 <= ALPHA < 1.0.
! NB    - Integer number of functions to be calculated, NB > 0.
!         The first function calculated is of order ALPHA, and the
!         last is of order (NB - 1 + ALPHA).
! IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated,
!         and 2 if exponentially scaled K's are to be calculated.
! BK    - Working precision output vector of length NB.  If the
!         routine terminates normally (NCALC=NB), the vector BK
!         contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
!         or the corresponding exponentially scaled functions.
!         If (0 < NCALC < NB), BK(I) contains correct function
!         values for I <= NCALC, and contains the ratios
!         K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
! NCALC - Integer output variable indicating possible errors.
!         Before using the vector BK, the user should check that
!         NCALC=NB, i.e., all orders have been calculated to
!         the desired accuracy.  See error returns below.
!
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!   EPS    = The smallest positive floating-point number such that
!            1.0+EPS > 1.0
!   XMAX   = Upper limit on the magnitude of X when IZE=1;  Solution
!            to equation:
!               W(X) * (1-1/8X+9/128X**2) = beta**minexp
!            where  W(X) = EXP(-X)*SQRT(PI/2X)
!   SQXMIN = Square root of beta**minexp
!   XINF   = Largest positive machine number; approximately
!            beta**maxexp
!   XMIN   = Smallest positive machine number; approximately
!            beta**minexp
!
!
!     Approximate values for some important machines are:
!
!                          beta       minexp      maxexp      EPS
!
!  CRAY-1        (S.P.)      2        -8193        8191    7.11E-15
!  Cyber 180/185
!    under NOS   (S.P.)      2         -975        1070    3.55E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)      2         -126         128    1.19E-7
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)      2        -1022        1024    2.22D-16
!  IBM 3033      (D.P.)     16          -65          63    2.22D-16
!  VAX           (S.P.)      2         -128         127    5.96E-8
!  VAX D-Format  (D.P.)      2         -128         127    1.39D-17
!  VAX G-Format  (D.P.)      2        -1024        1023    1.11D-16
!
!
!                         SQXMIN       XINF        XMIN      XMAX
!
! CRAY-1        (S.P.)  6.77E-1234  5.45E+2465  4.59E-2467 5674.858
! Cyber 180/855
!   under NOS   (S.P.)  1.77E-147   1.26E+322   3.14E-294   672.788
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)  1.08E-19    3.40E+38    1.18E-38     85.337
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)  1.49D-154   1.79D+308   2.23D-308   705.342
! IBM 3033      (D.P.)  7.35D-40    7.23D+75    5.40D-79    177.852
! VAX           (S.P.)  5.42E-20    1.70E+38    2.94E-39     86.715
! VAX D-Format  (D.P.)  5.42D-20    1.70D+38    2.94D-39     86.715
! VAX G-Format  (D.P.)  7.46D-155   8.98D+307   5.57D-309   706.728
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  In case of an error, NCALC /= NB, and not all K's are
!  calculated to the desired accuracy.
!
!  NCALC < -1:  An argument is out of range. For example,
!       NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >=
!       XMAX.  In this case, the B-vector is not calculated,
!       and NCALC is set to MIN0(NB,0)-2  so that NCALC /= NB.
!  NCALC = -1:  Either  K(ALPHA,X) >= XINF  or
!       K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) >= XINF.  In this case,
!       the B-vector is not calculated.  Note that again
!       NCALC /= NB.
!
!  0 < NCALC < NB: Not all requested function values could
!       be calculated accurately.  BK(I) contains correct function
!       values for I <= NCALC, and contains the ratios
!       K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT
!
!
! Acknowledgement
!
!  This program is based on a program written by J. B. Campbell
!  (2) that computes values of the Bessel functions K of real
!  argument and real order.  Modifications include the addition
!  of non-scaled functions, parameterization of machine
!  dependencies, and the use of more accurate approximations
!  for SINH and SIN.
!
! References: "On Temme's Algorithm for the Modified Bessel
!              Functions of the Third Kind," Campbell, J. B.,
!              TOMS 6(4), Dec. 1980, pp. 581-586.
!
!             "A FORTRAN IV Subroutine for the Modified Bessel
!              Functions of the Third Kind of Real Order and Real
!              Argument," Campbell, J. B., Report NRC/ERB-925,
!              National Research Council, Canada.
!
!  Latest modification: May 30, 1989
!
!  Modified by: W. J. Cody and L. Stoltz
!               Applied Mathematics Division
!               Argonne National Laboratory
!               Argonne, IL  60439
!
!-------------------------------------------------------------------

    LOGICAL :: run
    INTEGER :: i, iend, itemp, ize, j, k, m, mplus1, nb, ncalc
    REAL(num) :: alpha, blpha, bk1, bk2, c, dm, d1, d2, d3, enu, &
        ex, f0, f1, f2, p0, q0, ratio, twonu, twox, t1, t2, wminf, &
        x, x2by4, bk(nb), p(8), q(7), r(5), s(4), t(6), estm(6), estf(7)
    REAL(num), PARAMETER :: tinyx = 1.0e-10_num, &
        a = 0.11593151565841244881_num, d = 0.797884560802865364_num

!---------------------------------------------------------------------
!  Machine dependent parameters
!---------------------------------------------------------------------

    REAL(num), PARAMETER :: eps = EPSILON(1.0_num), sqxmin = SQRT(c_tiny)

!---------------------------------------------------------------------
!  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
!                                         + Euler's constant
!         Coefficients converted from hex to decimal and modified
!         by W. J. Cody, 2/26/82
!  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
!  T    - Approximation for SINH(Y)/Y
!---------------------------------------------------------------------

    p = (/  0.805629875690432845_num,      0.204045500205365151e2_num, &
            0.157705605106676174e3_num,    0.536671116469207504e3_num, &
            0.900382759291288778e3_num,    0.730923886650660393e3_num, &
            0.229299301509425145e3_num,    0.822467033424113231_num /)
    q = (/  0.294601986247850434e2_num,    0.277577868510221208e3_num, &
            0.120670325591027438e4_num,    0.276291444159791519e4_num, &
            0.344374050506564618e4_num,    0.221063190113378647e4_num, &
            0.572267338359892221e3_num /)
    r = (/ -0.48672575865218401848_num,    0.13079485869097804016e2_num, &
           -0.10196490580880537526e3_num,  0.34765409106507813131e3_num, &
            0.34958981245219347820e-3_num /)
    s = (/ -0.25579105509976461286e2_num,  0.21257260432226544008e3_num, &
           -0.61069018684944109624e3_num,  0.42269668805777760407e3_num /)
    t = (/  0.16125990452916363814e-9_num, 0.25051878502858255354e-7_num, &
            0.27557319615147964774e-5_num, 0.19841269840928373686e-3_num, &
            0.83333333333334751799e-2_num, 0.16666666666666666446_num /)
    estm = (/ 52.0583_num, 5.7607_num, 2.7782_num, 14.4303_num, 185.3004_num, &
              9.3715_num /)
    estf = (/ 41.8341_num, 7.1075_num, 6.4306_num, 42.5110_num, 1.35633_num, &
              84.5096_num, 20.0_num /)
!---------------------------------------------------------------------
    ex = x
    enu = alpha
    ncalc = MIN(nb, 0) - 2
    IF ((nb > 0) .AND. ((enu >= 0.0_num) .AND. (enu < 1.0_num)) &
        .AND. ((ize >= 1) .AND. (ize <= 2)) &
        .AND. ((ize /= 1) .OR. (ex <= c_largest_exp)) &
        .AND. (ex > 0.0_num)) THEN
      k = 0
      IF (enu < sqxmin) enu = 0.0_num
      IF (enu > 0.5_num) THEN
        k = 1
        enu = enu - 1.0_num
      END IF
      twonu = enu + enu
      iend = nb + k - 1
      c = enu * enu
      d3 = -c
      IF (ex <= 1.0_num) THEN
!---------------------------------------------------------------------
!  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
!                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
!---------------------------------------------------------------------
        d1 = 0.0_num
        d2 = p(1)
        t1 = 1.0_num
        t2 = q(1)
        DO i = 2, 7, 2
          d1 = c * d1 + p(i)
          d2 = c * d2 + p(i+1)
          t1 = c * t1 + q(i)
          t2 = c * t2 + q(i+1)
        END DO
        d1 = enu * d1
        t1 = enu * t1
        f1 = LOG(ex)
        f0 = a + enu * (p(8) - enu * (d1 + d2) / (t1 + t2)) - f1
        q0 = EXP(-enu * (a - enu * (p(8) + enu * (d1 - d2) / (t1 - t2)) - f1))
        f1 = enu * f0
        p0 = EXP(f1)
!---------------------------------------------------------------------
!  Calculation of F0 =
!---------------------------------------------------------------------
        d1 = r(5)
        t1 = 1.0_num
        DO i = 1, 4
          d1 = c * d1 + r(i)
          t1 = c * t1 + s(i)
        END DO
        IF (ABS(f1) <= 0.5_num) THEN
          f1 = f1 * f1
          d2 = 0.0_num
          DO i = 1, 6
            d2 = f1 * d2 + t(i)
          END DO
          d2 = f0 + f0 * f1 * d2
        ELSE
          d2 = SINH(f1) / enu
        END IF
        f0 = d2 - enu * d1 / (t1 * p0)
        IF (ex <= tinyx) THEN
!--------------------------------------------------------------------
!  X<=1.0E-10
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
!--------------------------------------------------------------------
          bk(1) = f0 + ex * f0
          IF (ize == 1) bk(1) = bk(1) - ex * bk(1)
          ratio = p0 / f0
          c = ex * c_largest_number
          IF (i /= 0) THEN
!--------------------------------------------------------------------
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
!  ALPHA >= 1/2
!--------------------------------------------------------------------
            ncalc = -1
            IF (bk(1) >= c / ratio) THEN
              rkbesl = bk(nb)
              RETURN
            END IF
            bk(1) = ratio * bk(1) / ex
            twonu = twonu + 2.0_num
            ratio = twonu
          END IF
          ncalc = 1
          IF (nb == 1) THEN
            rkbesl = bk(nb)
            RETURN
          END IF
!--------------------------------------------------------------------
!  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1
!--------------------------------------------------------------------
          ncalc = -1
          DO i = 2, nb
            IF (ratio >= c)  THEN
              rkbesl = bk(nb)
              RETURN
            END IF
            bk(i) = ratio / ex
            twonu = twonu + 2.0_num
            ratio = twonu
          END DO
          ncalc = 1
          DO i = 2, nb
            IF (bk(ncalc) >= c_largest_number / bk(i)) THEN
              rkbesl = bk(nb)
              RETURN
            END IF
            bk(i) = bk(ncalc) * bk(i)
            ncalc = i
          END DO
          rkbesl = bk(nb)
          RETURN
        ELSE
!--------------------------------------------------------------------
!  1.0E-10 < X <= 1.0
!--------------------------------------------------------------------
          c = 1.0_num
          x2by4 = ex * ex / 4.0_num
          p0 = 0.5_num * p0
          q0 = 0.5_num * q0
          d1 = -1.0_num
          d2 = 0.0_num
          bk1 = 0.0_num
          bk2 = 0.0_num
          f1 = f0
          f2 = p0
          run = .TRUE.
          DO WHILE(run)
            d1 = d1 + 2.0_num
            d2 = d2 + 1.0_num
            d3 = d1 + d3
            c = x2by4 * c / d2
            f0 = (d2 * f0 + p0 + q0) / d3
            p0 = p0 / (d2 - enu)
            q0 = q0 / (d2 + enu)
            t1 = c * f0
            t2 = c * (p0 - d2 * f0)
            bk1 = bk1 + t1
            bk2 = bk2 + t2
            run = (ABS(t1 / (f1 + bk1)) > eps) .OR. (ABS(t2 / (f2 + bk2)) > eps)
          END DO
          bk1 = f1 + bk1
          bk2 = 2.0_num * (f2 + bk2) / ex
          IF (ize == 2) THEN
            d1 = EXP(ex)
            bk1 = bk1 * d1
            bk2 = bk2 * d1
          END IF
          wminf = estf(1) * ex + estf(2)
        END IF
      ELSE IF (eps * ex > 1.0_num) THEN
!--------------------------------------------------------------------
!  X > 1/EPS
!--------------------------------------------------------------------
        ncalc = nb
        bk1 = 1.0_num / (d * SQRT(ex))
        DO i = 1, nb
          bk(i) = bk1
        END DO
        rkbesl = bk(nb)
        RETURN
      ELSE
!--------------------------------------------------------------------
!  X > 1.0
!--------------------------------------------------------------------
        twox = ex + ex
        blpha = 0.0_num
        ratio = 0.0_num
        IF (ex <= 4.0_num) THEN
!--------------------------------------------------------------------
!  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0
!--------------------------------------------------------------------
          d2 = AINT(estm(1) / ex + estm(2))
          m = INT(d2)
          d1 = d2 + d2
          d2 = d2 - 0.5_num
          d2 = d2 * d2
          DO i = 2, m
            d1 = d1 - 2.0_num
            d2 = d2 - d1
            ratio = (d3 + d2) / (twox + d1 - ratio)
          END DO
!--------------------------------------------------------------------
!  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
!    recurrence and K(ALPHA,X) from the wronskian
!--------------------------------------------------------------------
          d2 = AINT(estm(3) * ex + estm(4))
          m = INT(d2)
          c = ABS(enu)
          d3 = c + c
          d1 = d3 - 1.0_num
          f1 = c_tiny
          f0 = (2.0_num * (c + d2) / ex + 0.5_num * ex / (c + d2 + 1.0_num)) &
              * c_tiny
          DO i = 3, m
            d2 = d2 - 1.0_num
            f2 = (d3 + d2 + d2) * f0
            blpha = (1.0_num + d1 / d2) * (f2 + blpha)
            f2 = f2 / ex + f1
            f1 = f0
            f0 = f2
          END DO
          f1 = (d3 + 2.0_num) * f0 / ex + f1
          d1 = 0.0_num
          t1 = 1.0_num
          DO i = 1, 7
            d1 = c * d1 + p(i)
            t1 = c * t1 + q(i)
          END DO
          p0 = EXP(c * (a + c * (p(8) - c * d1 / t1) - LOG(ex))) / ex
          f2 = (c + 0.5_num - ratio) * f1 / ex
          bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0
          IF (ize == 1) bk1 = bk1 * EXP(-ex)
          wminf = estf(3) * ex + estf(4)
        ELSE
!--------------------------------------------------------------------
!  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
!  recurrence, for  X > 4.0
!--------------------------------------------------------------------
          dm = AINT(estm(5) / ex + estm(6))
          m = INT(dm)
          d2 = dm - 0.5_num
          d2 = d2 * d2
          d1 = dm + dm
          DO i = 2, m
            dm = dm - 1.0_num
            d1 = d1 - 2.0_num
            d2 = d2 - d1
            ratio = (d3 + d2) / (twox + d1 - ratio)
            blpha = (ratio + ratio * blpha) / dm
          END DO
          bk1 = 1.0_num / ((d + d * blpha) * SQRT(ex))
          IF (ize == 1) bk1 = bk1 * EXP(-ex)
          wminf = estf(5) * (ex - ABS(ex - estf(7))) + estf(6)
        END IF
!--------------------------------------------------------------------
!  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
!    K(ALPHA+1,X)/K(ALPHA,X)
!--------------------------------------------------------------------
        bk2 = bk1 + bk1 * (enu + 0.5_num - ratio) / ex
      END IF
!--------------------------------------------------------------------
!  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
!  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1
!--------------------------------------------------------------------
      ncalc = nb
      bk(1) = bk1
      IF (iend == 0) THEN
        rkbesl = bk(nb)
        RETURN
      END IF
      j = 2 - k
      IF (j > 0) bk(j) = bk2
      IF (iend == 1) THEN
        rkbesl = bk(nb)
        RETURN
      END IF
      m = MIN(INT(wminf - enu), iend)
      itemp = 1
      DO i = 2, m
        t1 = bk1
        bk1 = bk2
        twonu = twonu + 2.0_num
        IF (bk1 / MIN(1.0_num, ex) &
            >= (c_largest_number / twonu) * MIN(1.0_num, ex)) EXIT
        bk2 = twonu / ex * bk1 + t1
        itemp = i
        j = j + 1
        IF (j > 0) bk(j) = bk2
      END DO
      m = itemp
      IF (m == iend) THEN
        rkbesl = bk(nb)
        RETURN
      END IF
      ratio = bk2 / bk1
      mplus1 = m + 1
      ncalc = -1
      DO i = mplus1, iend
        twonu = twonu + 2.0_num
        ratio = twonu / ex + 1.0_num / ratio
        j = j + 1
        IF (j > 1) THEN
          bk(j) = ratio
        ELSE
          IF (bk2 >= c_largest_number / ratio) THEN
            rkbesl = bk(nb)
            RETURN
          END IF
          bk2 = ratio * bk2
        END IF
      END DO
      ncalc = MAX(mplus1 - k, 1)
      IF (ncalc == 1) bk(1) = bk2
      IF (nb == 1) THEN
        rkbesl = bk(nb)
        RETURN
      END IF
      j = ncalc + 1
      DO i = j, nb
        IF (bk(ncalc) >= c_largest_number / bk(i)) THEN
          rkbesl = bk(nb)
          RETURN
        END IF
        bk(i) = bk(ncalc) * bk(i)
        ncalc = i
      END DO
    END IF
    rkbesl = bk(nb)
    RETURN
!---------- Last line of rkbesl ----------
  END FUNCTION rkbesl

END MODULE numerics
