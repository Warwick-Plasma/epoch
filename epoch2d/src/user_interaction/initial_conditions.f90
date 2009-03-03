MODULE initial_conditions

  USE shared_data
  USE strings
  USE shunt
  USE dist_fn
  USE probes
  USE iocontrol
  USE laser
  USE simple_io

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: IC_Early,IC_Late,ManualLoad


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE IC_Early

    INTEGER :: iSpecies
    REAL(num) :: target_den=5.0e24_num
    REAL(num) :: spotcentre=0.0_num
    REAL(num) :: spotsize=20.0e-6_num
    TYPE(Laser_Block),POINTER :: Laser

    !Set the initial conditions
    DO iSpecies=1,nSpecies
       InitialConditions(iSpecies)%Temp=0.0_num
       DO iy=-2,ny+2
          InitialConditions(iSpecies)%Rho(:,iy)=target_den*(1.0_num-(1.0_num+&
               COS(2.0_num*pi/length_y * y(iy)))/2.0_num)
       ENDDO
    ENDDO

    !Now set the laser
    ALLOCATE(Laser)
    CALL Init_Laser(BD_LEFT,Laser)
    Laser%ID=1
    Laser%Amp=3.88191e12_num
    Laser%Freq=24.2831853e13_num
    Laser%UseTimeFunction=.FALSE.
    Laser%Phase=0.0_num
    Laser%Pol=0.0_num
    Laser%t_start=0.0_num
    Laser%t_end=4.0e-13_num
    DO iy=0,ny+1
       Laser%Profile(iy)=EXP(-((y(iy)-spotcentre)**2) / ((0.5_num*spotsize)**2))
    ENDDO

    CALL Attach_Laser(Laser)


  END SUBROUTINE IC_Early

  SUBROUTINE IC_Late


  END SUBROUTINE IC_Late

  SUBROUTINE ManualLoad
    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies



  END SUBROUTINE ManualLoad

END MODULE initial_conditions
