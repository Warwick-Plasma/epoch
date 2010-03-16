MODULE job_info

  IMPLICIT NONE

  TYPE :: jobid_type
    INTEGER(4) :: start_seconds
    INTEGER(4) :: start_milliseconds
  END TYPE jobid_type

  CONTAINS

  FUNCTION unix_seconds(values)

    INTEGER, DIMENSION(8), INTENT(IN) :: values
    INTEGER :: unix_seconds
    INTEGER :: days, year, month, day, h, m, s
    INTEGER, PARAMETER, DIMENSION(12) :: days_since_new_year = (/ &
        0, &
        31, &
        (31+28), &
        (31+28+31), &
        (31+28+31+30), &
        (31+28+31+30+31), &
        (31+28+31+30+31+30), &
        (31+28+31+30+31+30+31), &
        (31+28+31+30+31+30+31+31), &
        (31+28+31+30+31+30+31+31+30), &
        (31+28+31+30+31+30+31+31+30+31), &
        (31+28+31+30+31+30+31+31+30+31+30) /)

    year = values(1)
    month = values(2)
    day = values(3)
    h = values(5)
    m = values(6) - values(4)
    s = values(7)

    days = (year - 1970)*365 + (year - 1969)/4 + days_since_new_year(month)

    IF (MOD(year,400) .EQ. 0 .OR. &
        (MOD(year,4) .EQ. 0 .AND. MOD(year,100) .NE. 0)) days = days + 1

    unix_seconds = (((days + day - 1)*24 + h)*60 + m)*60 + s

  END FUNCTION unix_seconds



  FUNCTION get_unix_time()

    INTEGER :: get_unix_time
    INTEGER, DIMENSION(8) :: val

    CALL DATE_AND_TIME(values = val)

    get_unix_time = unix_seconds(val)

  END FUNCTION get_unix_time



  SUBROUTINE get_job_id(jobid)

    TYPE(jobid_type), INTENT(OUT) :: jobid
    INTEGER, DIMENSION(8) :: val

    CALL DATE_AND_TIME(values = val)

    jobid%start_seconds = unix_seconds(val)
    jobid%start_milliseconds = val(8)

  END SUBROUTINE get_job_id

END MODULE job_info
