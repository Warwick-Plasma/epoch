MODULE window

  USE shared_data
  USE boundary
  USE partlist
  USE helper

  IMPLICIT NONE
CONTAINS

  SUBROUTINE Shift_Window

    INTEGER :: iWindow

    !Shift the window round one cell at a time.
    !Inefficient, but it works
    DO iWindow=1,FLOOR(window_shift_fraction)
       CALL Insert_Particles
       !Shift the box around
       x_starts=x_starts+dx
       x_start_local=x_start_local+dx
       x_start=x_start+dx
       x_global=x_global+dx
       x=x+dx

       x_ends=x_ends+dx
       x_end_local=x_end_local+dx
       x_end=x_end+dx
       CALL Remove_Particles

       !Shift fields around
       CALL Shift_Field(Ex)
       CALL Shift_Field(Ey)
       CALL Shift_Field(Ez)

       CALL Shift_Field(Bx)
       CALL Shift_Field(By)
       CALL Shift_Field(Bz)
    ENDDO

  END SUBROUTINE Shift_Window

  SUBROUTINE Shift_Field(Field)

    REAL(num),DIMENSION(-2:nx+3,-2:ny+3,-2:nz+3),INTENT(INOUT) :: Field
    Field(-2:nx+2,:,:)=Field(-1:nx+3,:,:)
    CALL Field_BC(Field)

  END SUBROUTINE Shift_Field

  SUBROUTINE Insert_Particles

    TYPE(Particle),POINTER :: Current
    INTEGER :: iSpecies, iPart, i, iSuby, iSubz
    REAL(num) :: rand
    INTEGER :: clock,idum
    REAL(num) :: cell_y_r,cell_frac_y
    REAL(num) :: cell_z_r,cell_frac_z
    INTEGER :: cell_y,cell_z
    REAL(num),DIMENSION(-1:1) :: gy
    REAL(num),DIMENSION(-1:1) :: gz
    REAL(num) :: weight_local, temp_local

    !This subroutine injects particles at the right hand edge of the box

    !Only processors on the right need do anything
    IF (coordinates(3) .EQ. nprocx-1) THEN
       CALL SYSTEM_CLOCK(clock)
       idum=-(clock+rank+1)
       DO iSpecies=1,nSpecies
          DO iz=1,nz
             DO iy=1,ny
                DO iPart=1,ParticleSpecies(iSpecies)%window_npart_per_cell
                   ALLOCATE(Current)
                   rand=RANDOM(idum)-0.5_num
                   Current%Part_Pos(1)=x_end+dx + rand*dx
                   rand=RANDOM(idum)-0.5_num
                   Current%Part_Pos(2)=y(iy)+dy*rand
                   rand=RANDOM(idum)-0.5_num
                   Current%Part_Pos(3)=z(iz)+dz*rand

                   cell_y_r = (Current%Part_Pos(2)-y_start_local) / dy -0.5_num
                   cell_y=NINT(cell_y_r)
                   cell_frac_y = REAL(cell_y,num) - cell_y_r
                   cell_y=cell_y+1

                   cell_z_r = (Current%Part_Pos(3)-z_start_local) / dz -0.5_num
                   cell_z=NINT(cell_z_r)
                   cell_frac_z = REAL(cell_z,num) - cell_z_r
                   cell_z=cell_z+1

                   gy(-1) = 0.5_num * (0.5_num + cell_frac_y)**2
                   gy( 0) = 0.75_num - cell_frac_y**2
                   gy( 1) = 0.5_num * (0.5_num - cell_frac_y)**2

                   gz(-1) = 0.5_num * (0.5_num + cell_frac_z)**2
                   gz( 0) = 0.75_num - cell_frac_z**2
                   gz( 1) = 0.5_num * (0.5_num - cell_frac_z)**2

                   DO i=1,3
                      temp_local=0.0_num
                      DO iSubz=-1,+1
                         DO iSuby=-1,+1
                            temp_local=temp_local+gy(iSuby)*gz(iSubz)*ParticleSpecies(iSpecies)%window_Temperature(cell_y+iSuby,cell_z+iSubz,i)
                         ENDDO
                      ENDDO
                      Current%Part_P(i)=MomentumFromTemperature(ParticleSpecies(iSpecies)%Mass,temp_local,idum)
                   ENDDO

                   weight_local=0.0_num
                   DO iSubz=-1,+1
                      DO iSuby=-1,+1
                         weight_local=weight_local+gy(iSuby)*gz(iSubz)*ParticleSpecies(iSpecies)%Window_Density(cell_y+iSuby,cell_z+iSubz)&
                              /(REAL(ParticleSpecies(iSpecies)%window_npart_per_cell,num)/(dx*dy*dz))
                      ENDDO
                   ENDDO
                   Current%Weight=weight_local

#ifdef PART_DEBUG
                   Current%Processor=rank
                   Current%Processor_at_t0=rank
#endif
                   CALL Add_Particle_To_PartList(ParticleSpecies(iSpecies)%AttachedList,Current)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE Insert_Particles

  SUBROUTINE Remove_Particles

    TYPE(Particle),POINTER :: Current, Next
    INTEGER :: iSpecies

    IF (coordinates(2) .EQ. 0) THEN
       DO iSpecies=1,nSpecies
          Current=>ParticleSpecies(iSpecies)%AttachedList%Head
          DO WHILE(ASSOCIATED(Current))
             Next=>Current%Next
             IF (Current%Part_Pos(1) .LT. x_start-0.5_num*dx) THEN
                CALL Remove_Particle_From_PartList(ParticleSpecies(iSpecies)%AttachedList,Current)
                DEALLOCATE(Current)
             ENDIF
             Current=>Next
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE Remove_Particles

END MODULE window
