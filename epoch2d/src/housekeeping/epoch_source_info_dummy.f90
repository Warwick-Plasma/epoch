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

MODULE epoch_source_info

IMPLICIT NONE

CHARACTER(LEN=1) :: epoch_bytes_git_version = ''
CHARACTER(LEN=1) :: epoch_bytes_compile_date_string = ''
CHARACTER(LEN=1) :: epoch_bytes_compile_machine_info = ''
CHARACTER(LEN=1) :: epoch_bytes_compiler_info = ''
CHARACTER(LEN=1) :: epoch_bytes_compiler_flags = ''
INTEGER, PARAMETER :: epoch_bytes_compile_date = 1
CHARACTER(LEN=1) :: epoch_bytes_checksum_type = ''
CHARACTER(LEN=1) :: epoch_bytes_checksum = ''
CHARACTER(LEN=1) :: epoch_bytes_mimetype = ''
INTEGER, PARAMETER :: epoch_bytes_padding = 0
INTEGER, PARAMETER :: epoch_bytes_len = 0
INTEGER(8) :: epoch_bytes(1)
CHARACTER(LEN=1) :: epoch_diff_bytes_checksum_type = ''
CHARACTER(LEN=1) :: epoch_diff_bytes_checksum = ''
CHARACTER(LEN=1) :: epoch_diff_bytes_mimetype = ''
INTEGER, PARAMETER :: epoch_diff_bytes_padding = 0
INTEGER, PARAMETER :: epoch_diff_bytes_len = 0
INTEGER(8) :: epoch_diff_bytes(1)

END MODULE epoch_source_info
