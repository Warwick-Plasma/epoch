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

MODULE custom_laser

  USE shared_data

  IMPLICIT NONE

CONTAINS

  FUNCTION custom_laser_time_profile(laser)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num) :: custom_laser_time_profile

    custom_laser_time_profile = 1.0_num

  END FUNCTION custom_laser_time_profile

END MODULE custom_laser
