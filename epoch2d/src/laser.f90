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

    CALL Allocate_With_Direction(laser%profile,direction)
    CALL Allocate_With_Direction(laser%phase,direction)

    Laser%Direction=Direction
    NULLIFY(Laser%Next)

  END SUBROUTINE Init_Laser

  SUBROUTINE Delete_Laser(laser)
    TYPE(Laser_Block),INTENT(INOUT) :: laser

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)

  END SUBROUTINE Delete_Laser


  !Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE Attach_laser(laser)

    INTEGER :: direction
    TYPE(Laser_Block),POINTER,INTENT(INOUT) :: laser

    direction=laser%direction

    IF (Laser%k == 0) Laser%k=Laser%freq

    IF (direction == BD_LEFT .OR. direction == BD_RIGHT) THEN
       Laser%Phase(1:ny) = Laser%phase(1:ny) &
            - Laser%k * (y(1:ny) * TAN(laser%angle))
    ELSE IF (direction == BD_UP .OR. direction == BD_DOWN) THEN
       Laser%Phase(1:nx) = Laser%phase(1:nx) &
            - Laser%k * (x(1:nx) * TAN(laser%angle))
    ENDIF

    IF (direction == BD_LEFT) THEN
       CALL Attach_Laser_To_List(Laser_Left,Laser,direction) 
    ENDIF
    IF (direction == BD_RIGHT) THEN
       CALL Attach_Laser_To_List(Laser_Right,Laser,direction) 
    ENDIF
    IF (direction == BD_UP) THEN
       CALL Attach_Laser_To_List(Laser_Up,Laser,direction)
    ENDIF
    IF (direction == BD_DOWN) THEN
       CALL Attach_Laser_To_List(Laser_Down,Laser,direction)
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

  SUBROUTINE Allocate_With_Direction(Array,Direction)

    REAL(num),DIMENSION(:),POINTER :: Array
    INTEGER, INTENT(IN) :: Direction

    IF (direction == BD_LEFT .OR. direction == BD_RIGHT) THEN
       ALLOCATE(Array(-2:ny+3))
    ELSE IF (direction == BD_UP .OR. direction == BD_DOWN) THEN
       ALLOCATE(Array(-2:nx+3))
    ENDIF

  END SUBROUTINE Allocate_With_Direction

  SUBROUTINE Set_Laser_dt

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
    Current=>Laser_Up
    DO WHILE(ASSOCIATED(Current))
       dt_local=2.0_num*pi/Current%Freq
       dt_laser=MIN(dt_laser,dt_local)
       Current=>Current%Next
    ENDDO
    Current=>Laser_Down
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
    REAL(num),DIMENSION(:),ALLOCATABLE :: FPlus
    INTEGER :: err

    TYPE(Laser_Block),POINTER :: Current

    lx=dt/dx
    ALLOCATE(FPlus(1:ny))
    FPlus=0.0_num
    err=0
    Bx(0,1:ny) =  0.0_num


    Current=>Laser_Left
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus(1:ny) = Fplus(1:ny) + t_env * Current%Amp * Current%Profile(1:ny) &
               * SIN(Current%Freq*time + Current%Phase(1:ny)) * SIN(Current%Pol)&
               *COS(Current%Angle)
       ENDIF
       Current=>Current%Next
    ENDDO


    !Set the y magnetic field
    By(0,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (-4.0_num * Fplus &
         + 2.0_num * Ez(1,1:ny) - (C - lx*C**2)*By(1,1:ny) &
         - (dt / epsilon0) * Jz(1,1:ny))


    FPlus=0.0_num
    Current=>Laser_Left
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus(1:ny) = Fplus(1:ny) + t_env * Current%Amp * Current%Profile(1:ny) &
               * SIN(Current%Freq*time + Current%Phase(1:ny)) * COS(Current%Pol)
       ENDIF
       Current=>Current%Next
    ENDDO

    Bz(0,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (4.0_num * Fplus &
         - 2.0_num * Ey(1,1:ny) - (C - lx*C**2)*Bz(1,1:ny) &
         + (dt / epsilon0) * Jy(1,1:ny))
    DEALLOCATE(FPlus)
    
  END SUBROUTINE laser_bcs_left


  !Laser boundary for the right boundary
  SUBROUTINE laser_bcs_right
    REAL(num):: t_env
    REAL(num):: lx
    REAL(num),DIMENSION(:),ALLOCATABLE :: FMinus
    INTEGER :: err

    TYPE(Laser_Block),POINTER :: Current

    lx=dt/dx
    ALLOCATE(FMinus(1:ny))
    FMinus=0.0_num
    err=0
    Bx(nx+1,1:ny) =  0.0_num

    Current=>Laser_Right
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fminus(1:ny) = Fminus(1:ny) + t_env * Current%Amp * Current%Profile(1:ny) &
               * SIN(Current%Freq*time + Current%Phase(1:ny)) * SIN(Current%Pol)&
               *COS(Current%Angle)
       ENDIF
       Current=>Current%Next
    ENDDO

    By(nx,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (4.0_num * Fminus &
         - 2.0_num * Ez(nx,1:ny) - (C - lx*C**2)*By(nx,1:ny) &
         + (dt / epsilon0) * Jz(nx,1:ny))


    FMinus=0.0_num
    Current=>Laser_Right
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fminus(1:ny) = Fminus(1:ny) + t_env * Current%Amp * Current%Profile(1:ny) &
               * SIN(Current%Freq*time + Current%Phase(1:ny)) * COS(Current%Pol)
       ENDIF
       Current=>Current%Next
    ENDDO

    Bz(nx,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (-4.0_num * Fminus &
         + 2.0_num * Ey(nx,1:ny) - (C - lx*C**2)*Bz(nx,1:ny) &
         - (dt / epsilon0) * Jy(nx,1:ny))

    DEALLOCATE(FMinus)
 
  END SUBROUTINE laser_bcs_right


  !Laser boundary for the bottom boundary
  SUBROUTINE laser_bcs_down
    REAL(num):: t_env
    REAL(num):: ly
    REAL(num),DIMENSION(:),ALLOCATABLE :: FPlus
    INTEGER :: err

    TYPE(Laser_Block),POINTER :: Current

    ly=dt/dy
    ALLOCATE(FPlus(1:nx))
    FPlus=0.0_num
    err=0
    By(1:nx,-1) =  0.0_num


    Current=>Laser_Down
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus(1:nx) = Fplus(1:nx) + t_env * Current%Amp * Current%Profile(1:nx) &
               * SIN(Current%Freq*time + Current%Phase(1:nx)) * SIN(Current%Pol)&
               *COS(Current%Angle)
       ENDIF
       Current=>Current%Next
    ENDDO


    !Set the y magnetic field
    Bx(1:nx,0)=(1.0_num / (C + ly*C**2)) &
         * (-4.0_num * Fplus &
         - 2.0_num * Ez(1:nx,1) - (C - ly*C**2)*Bx(1:nx,1) &
         - (dt / epsilon0) * Jz(1:nx,1))


    FPlus=0.0_num
    Current=>Laser_Down
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus(1:nx) = Fplus(1:nx) + t_env * Current%Amp * Current%Profile(1:nx) &
               * SIN(Current%Freq*time + Current%Phase(1:nx)) * COS(Current%Pol)
       ENDIF
       Current=>Current%Next
    ENDDO

    Bz(1:nx,0)=(1.0_num / (C + ly*C**2)) &
         * (4.0_num * Fplus &
         + 2.0_num * Ex(1:nx,1) - (C - ly*C**2)*Bz(1:nx,1) &
         + (dt / epsilon0) * Jx(1:nx,1))
    DEALLOCATE(FPlus)

  END SUBROUTINE laser_bcs_down

  !Laser boundary for the bottom boundary
  SUBROUTINE laser_bcs_up
    REAL(num):: t_env
    REAL(num):: ly
    REAL(num),DIMENSION(:),ALLOCATABLE :: FPlus
    INTEGER :: err

    TYPE(Laser_Block),POINTER :: Current

    ly=dt/dy
    ALLOCATE(FPlus(1:nx))
    FPlus=0.0_num
    err=0
    By(1:nx,ny+1) =  0.0_num


    Current=>Laser_Down
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus(1:nx) = Fplus(1:nx) + t_env * Current%Amp * Current%Profile(1:nx) &
               * SIN(Current%Freq*time + Current%Phase(1:nx)) * SIN(Current%Pol)&
               *COS(Current%Angle)
       ENDIF
       Current=>Current%Next
    ENDDO


    !Set the x magnetic field
    Bx(1:nx,ny)=(1.0_num / (C + ly*C**2)) &
         * (-4.0_num * Fplus &
         + 2.0_num * Ez(1:nx,ny-1) - (C - ly*C**2)*Bx(1:nx,ny-1) &
         - (dt / epsilon0) * Jz(1:nx,ny-1))


    FPlus=0.0_num
    Current=>Laser_Down
    DO WHILE(ASSOCIATED(Current))
       !Evaluate the temporal evolution of the laser
       IF (time .GE. Current%t_start .AND. time .LE. Current%t_end) THEN
          t_env=Laser_Time_Profile(Current)
          Fplus(1:nx) = Fplus(1:nx) + t_env * Current%Amp * Current%Profile(1:nx) &
               * SIN(Current%Freq*time + Current%Phase(1:nx)) * COS(Current%Pol)
       ENDIF
       Current=>Current%Next
    ENDDO

    Bz(1:nx,ny)=(1.0_num / (C + ly*C**2)) &
         * (4.0_num * Fplus &
         - 2.0_num * Ex(1:nx,ny-1) + (C - ly*C**2)*Bz(1:nx,ny-1) &
         + (dt / epsilon0) * Jx(1:nx,ny-1))
    DEALLOCATE(FPlus)

  END SUBROUTINE laser_bcs_up

  SUBROUTINE outflow_bcs_left
    REAL(num):: lx

    lx=dt/dx
    Bx(0,1:ny) =  0.0_num
    !Set the y magnetic field
    By(0,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (2.0_num * Ez(1,1:ny) - (C - lx*C**2)*By(1,1:ny) &
         - (dt / epsilon0) * Jz(1,1:ny))
    Bz(0,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (-2.0_num * Ey(1,1:ny) - (C - lx*C**2)*Bz(1,1:ny) &
         + (dt / epsilon0) * Jy(1,1:ny))
    
  END SUBROUTINE outflow_bcs_left

  SUBROUTINE outflow_bcs_right
    REAL(num):: lx

    lx=dt/dx
    Bx(nx,1:ny) =  0.0_num
    !Set the y magnetic field
    By(nx,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (- 2.0_num * Ez(nx,1:ny) - (C - lx*C**2)*By(nx,1:ny) &
         + (dt / epsilon0) * Jz(nx,1:ny))
    Bz(nx,1:ny)=(1.0_num / (C + lx*C**2)) &
         * (2.0_num * Ey(nx,1:ny) - (C - lx*C**2)*Bz(nx,1:ny) &
         - (dt / epsilon0) * Jy(nx,1:ny))
  END SUBROUTINE outflow_bcs_right

  SUBROUTINE outflow_bcs_down
    REAL(num):: ly

    ly=dt/dy
    !Set the x magnetic field
    Bx(1:nx,0)=(1.0_num / (C + ly*C**2)) &
         * (-2.0_num * Ez(1:nx,1) - (C - ly*C**2)*Bx(1:nx,1) &
         + (dt / epsilon0) * Jz(1:nx,0))
    By(1:nx,0) =  0.0_num
    Bz(1:nx,0)=(1.0_num / (C + ly*C**2)) &
         * (+2.0_num * Ex(1:nx,1) - (C - ly*C**2)*Bz(1:nx,1) &
         - (dt / epsilon0) * Jx(1:nx,1))
    !CALL BField_bcs
  END SUBROUTINE outflow_bcs_down

  SUBROUTINE outflow_bcs_up
    REAL(num):: ly

    ly=dt/dy

    !Set the x magnetic field
    Bx(1:nx,ny)=(1.0_num / (C + ly*C**2)) &
         * (2.0_num * Ez(1:nx,ny) - (C - ly*C**2)*Bx(1:nx,ny) &
         - (dt / epsilon0) * Jz(1:nx,ny))
    By(1:nx,ny) =  0.0_num
    Bz(1:nx,ny)=(1.0_num / (C + ly*C**2)) &
         * (-2.0_num * Ex(1:nx,ny) + (C - ly*C**2)*Bz(1:nx,ny)&
         + (dt / epsilon0) * Jx(1:nx,ny))

    !CALL BField_bcs
  END SUBROUTINE outflow_bcs_up


END MODULE laser
