MODULE Custom_Laser
  USE shared_data
  IMPLICIT NONE

CONTAINS


  FUNCTION Custom_Laser_Time_Profile(laser)

    TYPE(Laser_Block),INTENT(IN) :: laser
    REAL(num) :: Custom_Laser_Time_Profile

    Custom_Laser_Time_Profile=1.0_num

  END FUNCTION Custom_Laser_Time_Profile

END MODULE Custom_Laser
