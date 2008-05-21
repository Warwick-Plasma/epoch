MODULE boundary

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Efield_bcs

    IF (xbc_left == periodic .OR. xbc_right == periodic) THEN
       Ex(-1:0,:)=Ex(nx-1:nx,:)
       Ey(-1:0,:)=Ey(nx-1:nx,:)
       Ez(-1:0,:)=Ez(nx-1:nx,:)
       Ex(nx+1:nx+2,:)=Ex(1:2,:)
       Ey(nx+1:nx+2,:)=Ey(1:2,:)
       Ez(nx+1:nx+2,:)=Ez(1:2,:)
       RETURN
    ELSE
       IF (xbc_left == other) THEN
          Ex(-1:0,:)=0.0_num
          Ey(-1:0,:)=0.0_num
          Ez(-1:0,:)=0.0_num
       ENDIF

       IF (xbc_right == other) THEN
          Ex(nx+1:nx+2,:)=0.0_num
          Ey(nx+1:nx+2,:)=0.0_num
          Ez(nx+1:nx+2,:)=0.0_num
       ENDIF
    ENDIF

    IF (ybc_up == periodic .OR. ybc_down == periodic) THEN
       Ex(:,-1:0)=Ex(:,ny-1:ny)
       Ey(:,-1:0)=Ey(:,ny-1:ny)
       Ez(:,-1:0)=Ez(:,ny-1:ny)
       Ex(:,ny+1:ny+2)=Ex(:,1:2)
       Ey(:,ny+1:ny+2)=Ey(:,1:2)
       Ez(:,ny+1:ny+2)=Ez(:,1:2)
       RETURN
    ELSE
       IF (ybc_down == other) THEN
          Ex(:,-1:0)=0.0_num
          Ey(:,-1:0)=0.0_num
          Ez(:,-1:0)=0.0_num
       ENDIF

       IF (ybc_up == other) THEN
          Ex(:,ny+1:ny+2)=0.0_num
          Ey(:,ny+1:ny+2)=0.0_num
          Ez(:,ny+1:ny+2)=0.0_num
       ENDIF
    ENDIF

  END SUBROUTINE Efield_bcs

  SUBROUTINE Bfield_bcs

    IF (xbc_left == periodic .OR. xbc_right == periodic) THEN
       Bx(-1:0,:)=Bx(nx-1:nx,:)
       By(-1:0,:)=By(nx-1:nx,:)
       Bz(-1:0,:)=Bz(nx-1:nx,:)

       Bx(nx+1:nx+2,:)=Bx(1:2,:)
       By(nx+1:nx+2,:)=By(1:2,:)
       Bz(nx+1:nx+2,:)=By(1:2,:)
    ELSE
       IF (xbc_left == other) THEN
          Bx(-1:0,:)=0.0_num
          By(-1,:)=0.0_num
          Bz(-1,:)=0.0_num
       ENDIF

       IF (xbc_right == other) THEN
          Bx(nx+1:nx+2,:)=0.0_num
          By(nx+1:nx+2,:)=0.0_num
          Bz(nx+1:nx+2,:)=0.0_num
       ENDIF
    ENDIF


    IF (ybc_up == periodic .OR. ybc_down == periodic) THEN
       Bx(:,-1:0)=Bx(:,ny-1:ny)
       By(:,-1:0)=By(:,ny-1:ny)
       Bz(:,-1:0)=Bz(:,ny-1:ny)

       Bx(:,ny+1:ny+2)=Bx(:,1:2)
       By(:,ny+1:ny+2)=By(:,1:2)
       Bz(:,ny+1:ny+2)=By(:,1:2)
    ELSE
       IF (ybc_down == other) THEN
          Bx(:,-1:0)=0.0_num
          By(:,-1:0)=0.0_num
          Bz(:,-1:0)=0.0_num
       ENDIF

       IF (ybc_up == other) THEN
          Bx(:,ny+1:ny+2)=0.0_num
          By(:,ny+1:ny+2)=0.0_num
          Bz(:,ny+1:ny+2)=0.0_num
       ENDIF
    ENDIF
  END SUBROUTINE Bfield_bcs

  SUBROUTINE Particle_bcs

    DO ipart=1,npart
       IF (xbc_left == periodic .OR. xbc_right == periodic) THEN
          IF (part_pos(ipart,1) .GT. x_end) THEN
             part_pos(ipart,1) = x_start + (part_pos(ipart,1)-x_end)
          ENDIF

          IF (part_pos(ipart,1) .LT. x_start) THEN
             part_pos(ipart,1) = x_end + (part_pos(ipart,1)-x_start)
          ENDIF
       ELSE
          IF (part_pos(ipart,1) .GT. x_end .AND. xbc_right == other) THEN
             part_pos(ipart,1) =  x_end - (part_pos(ipart,1)-x_end)
             part_p(ipart,1) = - part_p(ipart,1)
          ENDIF

          IF (part_pos(ipart,1) .LT. x_start .AND. xbc_left == other) THEN
             part_pos(ipart,1) =  x_start - (part_pos(ipart,1)-x_start)
             part_p(ipart,1) = - part_p(ipart,1)
          ENDIF
       ENDIF

       IF (ybc_up == periodic .OR. ybc_down == periodic) THEN
          IF (part_pos(ipart,2) .GT. y_end) THEN
             part_pos(ipart,2) = y_start + (part_pos(ipart,2)-y_end)
          ENDIF

          IF (part_pos(ipart,2) .LT. y_start) THEN
             part_pos(ipart,2) = y_end + (part_pos(ipart,2)-y_start)
          ENDIF
       ELSE
          IF (part_pos(ipart,2) .GT. y_end .AND. ybc_up == other) THEN
             part_pos(ipart,2) =  y_end - (part_pos(ipart,2)-y_end)
             part_p(ipart,2) = - part_p(ipart,2)
          ENDIF

          IF (part_pos(ipart,2) .LT. y_start .AND. ybc_down == other) THEN
             part_pos(ipart,2) =  y_start - (part_pos(ipart,2)-y_start)
             part_p(ipart,2) = - part_p(ipart,2)
          ENDIF
       ENDIF

    ENDDO

  END SUBROUTINE Particle_bcs

  SUBROUTINE Periodic_Summation_Bcs(Val)

    REAL(num),DIMENSION(-2:,-2:), INTENT(INOUT) :: Val

    IF (xbc_left .EQ. periodic .OR. xbc_right == periodic) THEN
       Val(1:3,:)=Val(1:3,:)+Val(nx+1:nx+3,:)
       Val(nx-2:nx,:)=Val(nx-2:nx,:)+Val(-2:0,:)
    ENDIF

    IF (ybc_up .EQ. periodic .OR. ybc_down == periodic) THEN
       Val(:,1:3)=Val(:,1:3)+Val(:,ny+1:ny+3)
       Val(:,ny-2:ny)=Val(:,ny-2:ny)+Val(:,-2:0)
    ENDIF

  END SUBROUTINE Periodic_Summation_Bcs

END MODULE boundary
