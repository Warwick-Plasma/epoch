MODULE window

  USE shared_data
  USE boundary
  USE partlist
  USE helper

  IMPLICIT NONE
CONTAINS

  SUBROUTINE shift_window

    INTEGER :: iwindow

    !Shift the window round one cell at a time.
    !Inefficient, but it works
    DO iwindow=1,FLOOR(window_shift_fraction)
       CALL insert_particles
       !Shift the box around
       x_starts=x_starts+dx
       x_start_local=x_start_local+dx
       x_start=x_start+dx
       x_global=x_global+dx
       x=x+dx

       x_ends=x_ends+dx
       x_end_local=x_end_local+dx
       x_end=x_end+dx
       CALL remove_particles

       !Shift fields around
       CALL shift_field(ex)
       CALL shift_field(ey)
       CALL shift_field(ez)

       CALL shift_field(bx)
       CALL shift_field(by)
       CALL shift_field(bz)
    ENDDO

  END SUBROUTINE shift_window

  SUBROUTINE shift_field(field)

    REAL(num),DIMENSION(-2:nx+3,-2:ny+3,-2:nz+3),INTENT(INOUT) :: field
    field(-2:nx+2,:,:)=field(-1:nx+3,:,:)
    CALL field_bc(field)

  END SUBROUTINE shift_field

  SUBROUTINE insert_particles

    TYPE(particle),POINTER :: current
    INTEGER :: ispecies, ipart, i, isuby, isubz
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
       DO ispecies=1,n_species
          DO iz=1,nz
             DO iy=1,ny
                DO ipart=1,particle_species(ispecies)%window_npart_per_cell
                   ALLOCATE(current)
                   rand=random(idum)-0.5_num
                   current%part_pos(1)=x_end+dx + rand*dx
                   rand=random(idum)-0.5_num
                   current%part_pos(2)=y(iy)+dy*rand
                   rand=random(idum)-0.5_num
                   current%part_pos(3)=z(iz)+dz*rand

                   cell_y_r = (current%part_pos(2)-y_start_local) / dy -0.5_num
                   cell_y=NINT(cell_y_r)
                   cell_frac_y = REAL(cell_y,num) - cell_y_r
                   cell_y=cell_y+1

                   cell_z_r = (current%part_pos(3)-z_start_local) / dz -0.5_num
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
                      DO isubz=-1,+1
                         DO isuby=-1,+1
                            temp_local=temp_local+gy(isuby)*gz(isubz)*particle_species(ispecies)%window_temperature(cell_y+isuby,cell_z+isubz,i)
                         ENDDO
                      ENDDO
                      current%part_p(i)=momentum_from_temperature(particle_species(ispecies)%mass,temp_local,idum)
                   ENDDO

                   weight_local=0.0_num
                   DO isubz=-1,+1
                      DO isuby=-1,+1
                         weight_local=weight_local+gy(isuby)*gz(isubz)*particle_species(ispecies)%window_density(cell_y+isuby,cell_z+isubz)&
                              /(REAL(particle_species(ispecies)%window_npart_per_cell,num)/(dx*dy*dz))
                      ENDDO
                   ENDDO
                   current%weight=weight_local

#ifdef PART_DEBUG
                   current%processor=rank
                   current%processor_at_t0=rank
#endif
                   CALL add_particle_to_partlist(particle_species(ispecies)%attached_list,current)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE insert_particles

  SUBROUTINE remove_particles

    TYPE(particle),POINTER :: current, next
    INTEGER :: ispecies

    IF (coordinates(2) .EQ. 0) THEN
       DO ispecies=1,n_species
          current=>particle_species(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(current))
             next=>current%next
             IF (current%part_pos(1) .LT. x_start-0.5_num*dx) THEN
                CALL remove_particle_from_partlist(particle_species(ispecies)%attached_list,current)
                DEALLOCATE(current)
             ENDIF
             current=>next
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE remove_particles

END MODULE window
