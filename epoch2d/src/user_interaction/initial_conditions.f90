MODULE initial_conditions

  USE shared_data
  USE strings
  USE shunt
  USE dist_fn
  USE iocontrol

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: IC_Early,IC_Late,ManualLoad


CONTAINS

  !-----------------------------------------------------------------------------
  !This function contains the equilibrium
  !-----------------------------------------------------------------------------

  SUBROUTINE IC_Early


  END SUBROUTINE IC_Early

  SUBROUTINE IC_Late

    INTEGER :: iSpecies

    move_window=.FALSE.
    window_v_x=3.0e8_num
    window_start_time=7.0e-13_num
    xbc_left_after_move=BC_SIMPLE_OUTFLOW
    xbc_right_after_move=xbc_right

    DO iSpecies=1,nSpecies
       DO ix=1,nx
          InitialConditions(iSpecies)%Rho(ix,:)=5.0e26_num*EXP(-((x(ix)/20.0e-6_num)**2))
       ENDDO
       DO iy=-1,ny+2
          DO ix=-1,nx+2
             InitialConditions(iSpecies)%Temp(ix,iy,2)=0.0e6_num!*EXP(-((y(iy)/2.0e-6_num)**2))
             InitialConditions(iSpecies)%Temp(ix,iy,1)=0.0e6_num!*(1.0_num + 99.0 * EXP(-((x(ix)/20.0e-6_num)**2)))
          ENDDO
       ENDDO

       ParticleSpecies(iSpecies)%Density=1.0_num
       ParticleSpecies(iSpecies)%Temperature=1.0e6_num
    ENDDO

  END SUBROUTINE IC_Late

  SUBROUTINE ManualLoad
    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies

    REAL(num),DIMENSION(3,2) :: ranges
    LOGICAL,DIMENSION(5) :: use_restrictions
    REAL(num),DIMENSION(5,2) :: restrictions
    INTEGER,DIMENSION(3) :: resolution=(/400,400,400/)

    use_restrictions=.FALSE.

    CALL CFD_OPEN("out.cfd",rank,comm,MPI_MODE_CREATE + MPI_MODE_WRONLY)
!!$    CALL general_3d_dist_fn("testx_px_py",(/DIR_X,DIR_PX,DIR_PY/),ranges,resolution,1,restrictions,use_restrictions)
!!$    ranges=0.0_num
!!$    resolution=100
!!$    CALL general_3d_dist_fn("testy_px_py",(/DIR_Y,DIR_PX,DIR_PY/),ranges,resolution,1,restrictions,use_restrictions)
!!$    ranges=0.0_num
!!$    resolution=100
    CALL general_3d_dist_fn("testx_y_px",(/DIR_X,DIR_Y,DIR_PX/),ranges,resolution,1,restrictions,use_restrictions)
    ranges=0.0_num
    resolution=400
    CALL general_2d_dist_fn("testx_px",(/DIR_X,DIR_PX/),ranges,resolution,1,restrictions,use_restrictions)
!!$    ranges=0.0_num
!!$    resolution=100
!!$    CALL general_2d_dist_fn("number_density",(/DIR_X,DIR_Y/),ranges,resolution,1,restrictions,use_restrictions)
    CALL cfd_close
   ! STOP


  END SUBROUTINE ManualLoad

END MODULE initial_conditions
