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

MODULE version_data

  IMPLICIT NONE
  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9

  CHARACTER(LEN=*), PARAMETER :: c_code_name = 'EPOCH2D'
  INTEGER(i4) :: c_version, c_revision, c_minor_rev
  INTEGER(i4), PARAMETER :: c_code_io_version = 1
  CHARACTER(LEN=*), PARAMETER :: c_commit_id = &
_COMMIT
  CHARACTER(LEN=*), PARAMETER :: c_compile_machine = &
_MACHINE
  CHARACTER(LEN=*), PARAMETER :: c_compile_flags = 'unknown'
  INTEGER(i4), PARAMETER :: c_compile_date = _DATE
  CHARACTER(LEN=16) :: version_string
  CHARACTER(LEN=70) :: ascii_header

END MODULE version_data
