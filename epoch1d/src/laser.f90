MODULE laser

  !This module contributed by
  !Dr. N. J. Sircombe

  USE shared_data
  USE partlist
  USE shared_parser_data
  USE evaluator
  USE shunt
  USE custom_laser
  USE boundary
  IMPLICIT NONE

CONTAINS

  SUBROUTINE Init_Laser(direction,laser)
    INTEGER,INTENT(IN) :: direction
    TYPE(Laser_Block),INTENT(INOUT) :: laser

    Laser%Direction=Direction
    NULLIFY(Laser%Next)

  END SUBROUTINE Init_Laser


  !Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE Attach_laser(laser)

    INTEGER :: direction
    TYPE(Laser_Block),POINTER,INTENT(INOUT) :: laser

    direction=laser%direction

    IF (direction == BD_LEFT) THEN
       CALL Attach_Laser_To_List(Laser_Left,Laser,direction) 
    ENDIF
    IF (direction == BD_RIGHT) THEN
       CALL Attach_Laser_To_List(Laser_Right,Laser,direction) 
    ENDIF

  END SUBROUTINE Attach_laser

  FUNCTION Laser_Time_Profile(laser)

    TYPE(Laser_Block),POINTER,INTENT(IN) :: laser
    REAL(num) :: Laser_Time_Profile
    INTEGER :: err

    err=0
    IF (Laser%UseTimeFunction) THEN
       Laser_Time_Profile=evaluate(laser%TimeFunction,err)
       RETURN
    ENDIF

    Laser_Time_Profile=Custom_Laser_Time_Profile(laser)

  END FUNCTION Laser_Time_Profile

  !Actually does the attaching of the laser to the correct list
  SUBROUTINE Attach_Laser_To_List(list,laser,Direction)

    TYPE(Laser_Block),INTENT(INOUT),POINTER :: list
    TYPE(Laser_Block),INTENT(IN),POINTER :: laser
    INTEGER,INTENT(IN) :: Direction
    TYPE(Laser_Block),DIMENSION(:),ALLOCATABLE :: temp
    INTEGER :: iLaser
    TYPE(Laser_Block),POINTER :: Current

    IF (ASSOCIATED(list)) THEN
       Current=>list
       DO WHILE(ASSOCIATED(Current%Next))
          Current=>Current%Next
       ENDDO
       Current%Next=>laser
    ELSE
       list=>laser
    ENDIF

  END SUBROUTINE Attach_Laser_To_List

  SUBROUTINE Set_Laser_dt

    INTEGER :: iLaser
    REAL(num) :: dt_local
    TYPE(Laser_Block),POINTER :: Current

    dt_laser=1000000.0_num
    Current=>Laser_Left
    DO WHILE(ASSOCIATED(Current))
       dt_local=2.0_num*pi/Current%Freq
       dt_laser=MIN(dt_laser,dt_local)
       Current=>Current%Next
    ENDDO
    Current=>Laser_Right
    DO WHILE(ASSOCIATED(Current))
       dt_local=2.0_num*pi/Current%Freq
       dt_laser=MIN(dt_laser,dt_local)
       Current=>Current%Next
    ENDDO

    !Need at least two iterations per laser period
    !(Nyquist)
    dt_laser=dt_laser/2.0_num

  END SUBROUTINE Set_Laser_dt

  !Laser boundary for the left boundary
  SUBROUTINE laser_bcs_left
    REAL(num):: t_env
    REAL(num):: lx
    REAL(num) :: FPlus
    INTEGER :: iLaser
    INTEGER :: err

    TYPE(Laser_Block),POINTER :: Current

    lx=dt/dx
    FPlus=0.0_num
    err=0
    Bx(0) =  0.0_num


    Current=>Laser_Left
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus = Fplus + t_env * Current%Amp &
               * SIN(Current%Freq*time + Current%Phase) * SIN(Current%Pol)&
               *COS(Current%Angle)
       ENDIF
       Current=>Current%Next
    ENDDO

    !Set the y magnetic field
    By(-1:0)=(1.0_num / (C + lx*C**2)) &
         * (-4.0_num * Fplus &
         + 2.0_num * Ez(1) - (C - lx*C**2)*By(1) &
         - (dt / epsilon0) * Jz(1))


    FPlus=0.0_num
    Current=>Laser_Left
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus = Fplus + t_env * Current%Amp &
               * SIN(Current%Freq*time + Current%Phase) * COS(Current%Pol)
       ENDIF
       Current=>Current%Next
    ENDDO

    Bz(-1:0)=(1.0_num / (C + lx*C**2)) &
         * (4.0_num * Fplus &
         - 2.0_num * Ey(1) - (C - lx*C**2)*Bz(1) &
         + (dt / epsilon0) * Jy(1))
    
  END SUBROUTINE laser_bcs_left


  !Laser boundary for the right boundary
  SUBROUTINE laser_bcs_right
    REAL(num):: t_env
    REAL(num):: lx
    REAL(num):: FMinus
    INTEGER :: iLaser
    INTEGER :: err

    TYPE(Laser_Block),POINTER :: Current

    lx=dt/dx

    FMinus=0.0_num
    err=0
    Bx(nx+1) =  0.0_num

    Current=>Laser_Right
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fminus = Fminus + t_env * Current%Amp &
               * SIN(Current%Freq*time + Current%Phase) * SIN(Current%Pol)&
               *COS(Current%Angle)
       ENDIF
       Current=>Current%Next
    ENDDO

    By(nx:nx+1)=(1.0_num / (C + lx*C**2)) &
         * (4.0_num * Fminus &
         - 2.0_num * Ez(nx) - (C - lx*C**2)*By(nx) &
         + (dt / epsilon0) * Jz(nx))


    FMinus=0.0_num
    Current=>Laser_Right
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fminus = Fminus + t_env * Current%Amp * Current%Profile &
               * SIN(Current%Freq*time + Current%Phase) * COS(Current%Pol)
       ENDIF
       Current=>Current%Next
    ENDDO

    Bz(nx:nx+1)=(1.0_num / (C + lx*C**2)) &
         * (-4.0_num * Fminus &
         + 2.0_num * Ey(nx) - (C - lx*C**2)*Bz(nx) &
         - (dt / epsilon0) * Jy(nx))

  END SUBROUTINE laser_bcs_right

  SUBROUTINE outflow_bcs_left
    REAL(num):: lx

    lx=dt/dx
    Bx(-1:0) =  0.0_num
    !Set the y magnetic field
    By(-1:0)=(1.0_num / (C + lx*C**2)) &
         * (2.0_num * Ez(1) - (C - lx*C**2)*By(1) &
         - (dt / epsilon0) * Jz(1))
    Bz(-1:0)=(1.0_num / (C + lx*C**2)) &
         * (-2.0_num * Ey(1) - (C - lx*C**2)*Bz(1) &
         + (dt / epsilon0) * Jy(1))
    
  END SUBROUTINE outflow_bcs_left

  SUBROUTINE outflow_bcs_right
    REAL(num):: lx

    lx=dt/dx
    Bx(nx) =  0.0_num
    !Set the y magnetic field
    By(nx:nx+1)=(1.0_num / (C + lx*C**2)) &
         * (- 2.0_num * Ez(nx) - (C - lx*C**2)*By(nx) &
         + (dt / epsilon0) * Jz(nx))
    Bz(nx:nx+1)=(1.0_num / (C + lx*C**2)) &
         * (2.0_num * Ey(nx) - (C - lx*C**2)*Bz(nx) &
         - (dt / epsilon0) * Jy(nx))
  END SUBROUTINE outflow_bcs_right


END MODULE laser
