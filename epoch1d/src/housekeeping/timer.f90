MODULE timer

  USE shared_data
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
    ENDIF

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
    ENDIF

  END SUBROUTINE timer_stop



  SUBROUTINE timer_reset

    INTEGER, PARAMETER :: id = c_timer_dt
    LOGICAL, SAVE :: first = .TRUE.

    IF (first) THEN
      first = .FALSE.
      timer_avg_first(c_timer_step) = .TRUE.
      RETURN
    ENDIF

    timer_time(id) = timer_time(c_timer_step) - timer_time(c_timer_io) &
        - timer_time(c_timer_balance)

    IF (timer_avg_first(id)) THEN
      timer_avg_first(id) = .FALSE.
      timer_average(id) = timer_time(id)
    ELSE
      timer_average(id) = avg_weight1 * timer_time(id) &
          + avg_weight2 * timer_average(id)
    ENDIF

    timer_first = 0.0_num
    timer_time = 0.0_num

  END SUBROUTINE timer_reset

END MODULE timer
