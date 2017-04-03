! Copyright (C) 2017      Chris Brady <C.S.Brady@warwick.ac.uk>
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

#if defined(__GFORTRAN__) || defined(__INTEL_COMPILER)
#define USE_ISATTY
#endif

MODULE terminal_controls

  USE, INTRINSIC :: iso_fortran_env

  CHARACTER(LEN=5), DIMENSION(12) :: vt100_control = (/'[39m','[30m','[31m',&
      '[32m','[33m','[34m','[35m','[36m','[1m ','[2m ','[4m ','[0m '/)

  INTEGER, PARAMETER :: c_term_default_colour = 1
  INTEGER, PARAMETER :: c_term_black = 2
  INTEGER, PARAMETER :: c_term_red = 3
  INTEGER, PARAMETER :: c_term_green = 4
  INTEGER, PARAMETER :: c_term_yellow = 5
  INTEGER, PARAMETER :: c_term_blue = 6
  INTEGER, PARAMETER :: c_term_magenta = 7
  INTEGER, PARAMETER :: c_term_cyan = 8
  INTEGER, PARAMETER :: c_term_bold = 9
  INTEGER, PARAMETER :: c_term_dim = 10
  INTEGER, PARAMETER :: c_term_underline = 11
  INTEGER, PARAMETER :: c_term_reset_attributes = 12

  INTEGER, PARAMETER :: c_term_max = 12

CONTAINS

  SUBROUTINE set_term_attr(controlcode)

    INTEGER, INTENT(IN) :: controlcode
#ifdef USE_ISATTY
    LOGICAL, SAVE :: first = .TRUE.
    LOGICAL, SAVE :: tty = .FALSE.

    IF (first) THEN
      first = .FALSE.
      tty = ISATTY(output_unit)
    END IF

    IF (.NOT. tty) RETURN

    IF (controlcode < 1 .OR. controlcode > c_term_max) RETURN

    WRITE(*,'(A)',ADVANCE='NO') ACHAR(27) // TRIM(vt100_control(controlcode))
#endif

  END SUBROUTINE set_term_attr

END MODULE terminal_controls
