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

MODULE timer

  USE constants
  USE mpi

  IMPLICIT NONE

  LOGICAL :: timer_collect
  REAL(num), PARAMETER :: avg_weight1 = 0.4_num
  REAL(num), PARAMETER :: avg_weight2 = 1.0_num - avg_weight1
  INTEGER, PARAMETER :: c_timer_step = 1
  INTEGER, PARAMETER :: c_timer_dt = 2
  INTEGER, PARAMETER :: c_timer_io = 3
  INTEGER, PARAMETER :: c_timer_balance = 4
  INTEGER, PARAMETER :: c_timer_max = 4

  REAL(num) :: timer_walltime
  REAL(num) :: timer_first(c_timer_max)
  REAL(num) :: timer_time(c_timer_max)
  REAL(num) :: timer_average(c_timer_max)
  LOGICAL :: timer_avg_first(c_timer_max)

CONTAINS

  SUBROUTINE timer_init

    timer_first = 0.0_num
    timer_time = 0.0_num
    timer_average = 0.0_num
    timer_avg_first = .TRUE.
    timer_collect = .FALSE.

  END SUBROUTINE timer_init



  SUBROUTINE timer_start(id, use_old)

    INTEGER, INTENT(IN) :: id
    LOGICAL, INTENT(IN), OPTIONAL :: use_old

    IF (PRESENT(use_old)) THEN
      IF (.NOT.use_old) timer_walltime = MPI_WTIME()
    ELSE
      timer_walltime = MPI_WTIME()
    END IF

    timer_first(id) = timer_walltime

  END SUBROUTINE timer_start



  SUBROUTINE timer_stop(id)

    INTEGER, INTENT(IN) :: id

    timer_walltime = MPI_WTIME()

    timer_time(id) = timer_walltime - timer_first(id)
    IF (timer_avg_first(id)) THEN
      timer_avg_first(id) = .FALSE.
      timer_average(id) = timer_time(id)
    ELSE
      timer_average(id) = avg_weight1 * timer_time(id) &
          + avg_weight2 * timer_average(id)
    END IF

  END SUBROUTINE timer_stop



  SUBROUTINE timer_reset

    INTEGER, PARAMETER :: id = c_timer_dt
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      first = .FALSE.
      timer_avg_first(c_timer_step) = .TRUE.
      RETURN
    END IF

    timer_time(id) = timer_time(c_timer_step) - timer_time(c_timer_io) &
        - timer_time(c_timer_balance)

    IF (timer_avg_first(id)) THEN
      timer_avg_first(id) = .FALSE.
      timer_average(id) = timer_time(id)
    ELSE
      timer_average(id) = avg_weight1 * timer_time(id) &
          + avg_weight2 * timer_average(id)
    END IF

    timer_first = 0.0_num
    timer_time = 0.0_num

  END SUBROUTINE timer_reset

END MODULE timer
