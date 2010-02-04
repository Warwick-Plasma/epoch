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

       CALL shift_field(jx)
       CALL shift_field(jy)
       CALL shift_field(jz)

       CALL shift_field(bx)
       CALL shift_field(by)
       CALL shift_field(bz)
    ENDDO

  END SUBROUTINE shift_window

  SUBROUTINE shift_field(field)

    REAL(num),DIMENSION(-2:nx+3),INTENT(INOUT) :: field
    field(-2:nx+2)=field(-1:nx+3)
    CALL field_bc(field)

  END SUBROUTINE shift_field

  SUBROUTINE insert_particles

    TYPE(particle),POINTER :: current
    INTEGER :: ispecies, ipart,i
    REAL(num) :: rand
    INTEGER :: clock,idum
    REAL(num) :: cell_x_r,cell_frac_x
    INTEGER :: cell_x
    REAL(num),DIMENSION(-1:1) :: gx
    REAL(num) :: weight_local, temp_local


    !This subroutine injects particles at the right hand edge of the box

    !Only processors on the right need do anything
    IF (coordinates(1) .EQ. nprocx-1) THEN
       CALL SYSTEM_CLOCK(clock)
       idum=-(clock+rank+1)
       DO ispecies=1,n_species
             DO ipart=1,particle_species(ispecies)%npart_per_cell
                ALLOCATE(current)
                rand=random(idum)-0.5_num
                current%part_pos=x_end+dx + rand*dx

                cell_x_r = (current%part_pos-x_start_local) / dx -0.5_num
                cell_x=NINT(cell_x_r)
                cell_frac_x = REAL(cell_x,num) - cell_x_r
                cell_x=cell_x+1

!!$                IF (cell_y .NE. iy) PRINT *,"BAD CELL"

                gx(-1) = 0.5_num * (0.5_num + cell_frac_x)**2
                gx( 0) = 0.75_num - cell_frac_x**2
                gx( 1) = 0.5_num * (0.5_num - cell_frac_x)**2

                DO i=1,3
                   temp_local=particle_species(ispecies)%temperature(i)
                   current%part_p(i)=momentum_from_temperature(particle_species(ispecies)%mass,temp_local,idum)
                ENDDO

                weight_local=particle_species(ispecies)%density/&
						(REAL(particle_species(ispecies)%npart_per_cell,num)/(dx))
#ifdef PER_PARTICLE_WEIGHT
                current%weight=weight_local
#endif
#ifdef PART_DEBUG
                current%processor=rank
                current%processor_at_t0=rank
#endif
                CALL add_particle_to_partlist(particle_species(ispecies)%attached_list,current)
             ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE insert_particles

  SUBROUTINE remove_particles

    TYPE(particle),POINTER :: current, next
    INTEGER :: ispecies

    IF (coordinates(1) .EQ. 0) THEN
       DO ispecies=1,n_species
          current=>particle_species(ispecies)%attached_list%head
          DO WHILE(ASSOCIATED(current))
             next=>current%next
             IF (current%part_pos .LT. x_start-0.5_num*dx) THEN
                CALL remove_particle_from_partlist(particle_species(ispecies)%attached_list,current)
                DEALLOCATE(current)
             ENDIF
             current=>next
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE remove_particles

END MODULE window
