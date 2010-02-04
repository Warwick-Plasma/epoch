MODULE dist_fn

  USE shared_data
  USE partlist
  USE output_cartesian 
  USE output
  USE iocontrol
  IMPLICIT NONE

CONTAINS

  SUBROUTINE attach_dist_fn(block)
    TYPE(distribution_function_block),POINTER :: block
    TYPE(distribution_function_block),POINTER :: current

    current=>dist_fns
    IF (.NOT. ASSOCIATED(current)) THEN
       !This is the first distribution function to add
       dist_fns=>block
       RETURN
    ENDIF
    DO WHILE(ASSOCIATED(current%next))
       current=>current%next
    ENDDO
    current%next=>block

  END SUBROUTINE attach_dist_fn

  SUBROUTINE setup_dist_fn(block)
    TYPE(distribution_function_block),POINTER :: block

    NULLIFY(block%next)
    ALLOCATE(block%use_species(1:n_species))
    block%use_species=.FALSE.
    block%ranges=1.0_num
    block%use_restrictions=.FALSE.
    block%ndims=-1
  END SUBROUTINE setup_dist_fn

  SUBROUTINE clean_dist_fns

    TYPE(distribution_function_block),POINTER :: current,next

    current=>dist_fns
    DO WHILE(ASSOCIATED(current))
       next=>current%next
       DEALLOCATE(current)
       current=>next
    ENDDO

  END SUBROUTINE clean_dist_fns

  SUBROUTINE write_dist_fns(code)
    INTEGER,INTENT(IN) :: code

    INTEGER :: ispecies
    REAL(num),DIMENSION(3,2) :: ranges
    LOGICAL,DIMENSION(5) :: use_restrictions
    REAL(num),DIMENSION(5,2) :: restrictions
    INTEGER,DIMENSION(3) :: resolution=(/100,100,100/)
    TYPE(distribution_function_block),POINTER :: current

    !Write the distribution functions
    current=>dist_fns
    DO WHILE(ASSOCIATED(current))
       IF (IAND(current%dumpmask,code) .NE. 0) THEN
          DO ispecies=1,n_species
             IF (.NOT. current%use_species(ispecies)) CYCLE
             ranges=current%ranges
             resolution=current%resolution
             restrictions=current%restrictions
             use_restrictions=current%use_restrictions

             IF (current%ndims .EQ. 2) THEN
                CALL general_2d_dist_fn(current%name,current%directions(1:2),ranges(1:2,:),resolution(1:2),ispecies,restrictions,use_restrictions)
             ELSE
                CALL general_3d_dist_fn(current%name,current%directions,ranges,resolution,ispecies,restrictions,use_restrictions)
             ENDIF
             IF (current%store_ranges) current%ranges=ranges
          ENDDO
       ENDIF
       current=>current%next
    ENDDO
  END SUBROUTINE write_dist_fns


  SUBROUTINE general_3d_dist_fn(name,direction,range,resolution,species,restrictions,use_restrictions)

    CHARACTER(len=*),INTENT(IN) :: name
    INTEGER,DIMENSION(3),INTENT(IN) :: direction
    REAL(num),DIMENSION(3,2),INTENT(INOUT) :: range
    INTEGER,DIMENSION(3),INTENT(INOUT) :: resolution
    INTEGER,INTENT(IN) :: species
    REAL(num),DIMENSION(5,2),INTENT(IN) :: restrictions
    LOGICAL,DIMENSION(5),INTENT(IN) :: use_restrictions

    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: data, data2
    REAL(num),DIMENSION(:),ALLOCATABLE :: grid1,grid2,grid3
    LOGICAL,DIMENSION(3) :: parallel
    REAL(num),DIMENSION(3) :: dgrid
    REAL(num) :: current_data, temp_data
    INTEGER :: idim,idir
    LOGICAL,DIMENSION(3) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x,use_y,need_reduce
    INTEGER,DIMENSION(3) :: start_local,global_resolution
    INTEGER :: color
    INTEGER :: comm_new, type_new

    LOGICAL,DIMENSION(3,5) :: use_direction
    LOGICAL,DIMENSION(3) :: calc_mod
    INTEGER,DIMENSION(3) :: p_count
    INTEGER,DIMENSION(3) :: l_direction
    REAL(num),DIMENSION(3) :: conv 
    INTEGER,DIMENSION(3) :: cell
    LOGICAL :: use_this
    REAL(num) :: real_space_area, part_weight


    TYPE(particle),POINTER :: current
    CHARACTER(len=string_length) :: grid_name, norm_grid_name, var_name
    REAL(num),DIMENSION(3) :: stagger=0.0_num
    REAL(num),DIMENSION(5) :: particle_data

    REAL(num) :: max_p_conv

    INTEGER :: ind

    use_x=.FALSE.
    use_y=.FALSE.
    need_reduce=.TRUE.
    color=0
    global_resolution=resolution
    parallel=.FALSE.
    start_local=1
    calc_range=.FALSE.
    calc_ranges=.FALSE.
    p_count=0
    use_direction=.FALSE.
    l_direction=0

    real_space_area=1.0_num


    DO idim=1,3
       IF (direction(idim) .EQ. DIR_X) THEN
          use_x=.TRUE.
          resolution(idim) = nx
          range(idim,1)=x_start_local
          range(idim,2)=x_end_local
          start_local(idim)=cell_x_start(coordinates(2)+1)
          global_resolution(idim)=nx_global
          dgrid(idim)=dx
          parallel(idim)=.TRUE.
          use_direction(idim,1)=.TRUE.
          l_direction(idim)=1
          conv(idim)=MAX(length_x,length_y)
          CYCLE
       ENDIF
       IF (direction(idim) .EQ. DIR_Y)  THEN
          use_y=.TRUE.
          resolution(idim) = ny
          range(idim,1)=y_start_local
          range(idim,2)=y_end_local
          start_local(idim)=cell_y_start(coordinates(1)+1)
          global_resolution(idim)=ny_global
          dgrid(idim)=dy
          parallel(idim)=.TRUE.
          use_direction(idim,2)=.TRUE.
          l_direction(idim)=2
          conv(idim)=MAX(length_x,length_y)
          CYCLE
       ENDIF
       conv(idim)=c*m0
       !If we're here then this must be a momentum space direction
       !So determine which momentum space directions are needed
       IF (range(idim,1) .EQ. range(idim,2)) THEN
          calc_range(idim) = .TRUE.
          calc_ranges=.TRUE.
       ENDIF
       IF (IAND(direction(idim),DIR_PX) .NE. 0) THEN
          use_direction(idim,3)=.TRUE.
          p_count(idim)=p_count(idim)+1
          l_direction(idim)=3
       ENDIF
       IF (IAND(direction(idim),DIR_PY) .NE. 0) THEN
          use_direction(idim,4)=.TRUE.
          p_count(idim)=p_count(idim)+1
          l_direction(idim)=4
       ENDIF
       IF (IAND(direction(idim),DIR_PZ) .NE. 0) THEN
          use_direction(idim,5)=.TRUE.
          p_count(idim)=p_count(idim)+1
          l_direction(idim)=5
       ENDIF
    ENDDO

    IF (.NOT. use_x) real_space_area=real_space_area*dx
    IF (.NOT. use_y) real_space_area=real_space_area*dy

    DO idim=1,3
       IF (p_count(idim) .GT. 1) calc_mod(idim)=.TRUE.
    ENDDO
    !Calculate ranges for directions where needed
    IF (calc_ranges) THEN
       DO idim=1,3
          IF (calc_range(idim)) THEN
             range(idim,1)=1.0e6_num
             range(idim,2)=-1.0e6_num
          ENDIF
       ENDDO
       current=>particle_species(species)%attached_list%head
       ind =0
       DO WHILE(ASSOCIATED(current))
          particle_data(1:2)=current%part_pos
          particle_data(3:5)=current%part_p
          DO idim=1,3
             IF (calc_range(idim)) THEN
                IF (calc_mod(idim)) THEN
                   DO idir=1,5
                      IF (use_direction(idim,idir)) current_data=current_data+particle_data(idir)
                   ENDDO
                   current_data=SQRT(current_data)
                ELSE
                   current_data=particle_data(l_direction(idim))
                ENDIF
                use_this=.TRUE.
                DO idir=1,5
                   IF (use_restrictions(idir) .AND. (particle_data(idir) .LT. restrictions(idir,1) &
                        .OR. particle_data(idir) .GT. restrictions(idir,2))) use_this=.FALSE.
                ENDDO
                IF (use_this) THEN
                   IF (current_data .LT. range(idim,1)) range(idim,1)=current_data
                   IF (current_data .GT. range(idim,2)) range(idim,2)=current_data
                ENDIF
             ENDIF
          ENDDO
          ind=ind+1
          current=>current%next
       ENDDO
    ENDIF

    max_p_conv=-10.0_num
    DO idim=1,3
       IF (.NOT. parallel(idim)) THEN
          !If not parallel then this is a momentum dimension
          CALL MPI_ALLREDUCE(range(idim,1),temp_data,1,mpireal,MPI_MIN,comm,errcode)
          range(idim,1)=temp_data
          CALL MPI_ALLREDUCE(range(idim,2),temp_data,1,mpireal,MPI_MAX,comm,errcode)
          range(idim,2)=temp_data
       ENDIF
       !Fix so that if distribution function is zero then it picks an arbitrary scale in that direction
       IF (range(idim,1) .EQ. range(idim,2)) THEN
          range(idim,1)=-1.0_num
          range(idim,2)=1.0_num
       ENDIF
       !Calculate the maxmium range of a momentum direction
       IF (range(idim,2)-range(idim,1) .GT. max_p_conv .AND. .NOT. parallel(idim)) max_p_conv=range(idim,2)-range(idim,1)
    ENDDO

    DO idim=1,3
       IF (.NOT. parallel(idim)) conv(idim)=max_p_conv
    ENDDO

    !Setup MPI
    IF (use_x .AND. use_y) need_reduce=.FALSE.
    IF (.NOT. use_y) color=color+coordinates(2) !If using x direction need to reduce only across all y
    IF (.NOT. use_x) color=color+nprocx*coordinates(1) !If using y direction need to reduce only across all x

    IF (need_reduce) THEN
       CALL MPI_COMM_SPLIT(comm,color,rank,comm_new,errcode)
    ELSE
       comm_new=MPI_COMM_NULL
    ENDIF

    type_new=create_3d_field_subtype(resolution,global_resolution,start_local)
    !Create grids
    DO idim=1,3
       IF (.NOT. parallel(idim)) dgrid(idim)=(range(idim,2)-range(idim,1))/REAL(resolution(idim)-1,num)
    ENDDO
    ALLOCATE(grid1(0:global_resolution(1)),grid2(0:global_resolution(2)))
    ALLOCATE(grid3(0:global_resolution(3)))

    grid1(0)=-dgrid(1) + range(1,1)
    DO idir=1,global_resolution(1)
       grid1(idir)=grid1(idir-1)+dgrid(1)
    ENDDO

    grid2(0)=-dgrid(2) + range(2,1)
    DO idir=1,global_resolution(2)
       grid2(idir)=grid2(idir-1)+dgrid(2)
    ENDDO

    grid3(0)=-dgrid(3) + range(3,1)
    DO idir=1,global_resolution(3)
       grid3(idir)=grid3(idir-1)+dgrid(3)
    ENDDO

    ALLOCATE(data(1:resolution(1),1:resolution(2),1:resolution(3)))
    data=0.0_num

    current=>particle_species(species)%attached_list%head
    part_weight=weight
    DO WHILE(ASSOCIATED(current))
       particle_data(1:2)=current%part_pos
       particle_data(3:5)=current%part_p
       use_this=.TRUE.
       DO idir=1,5
          IF (use_restrictions(idir) .AND. (particle_data(idir) .LT. restrictions(idir,1) &
               .OR. particle_data(idir) .GT. restrictions(idir,2))) use_this=.FALSE.
       ENDDO
       IF (use_this) THEN
          DO idim=1,3
             IF (calc_mod(idim)) THEN
                DO idir=1,5
                   IF (use_direction(idim,idir)) current_data=current_data+particle_data(idir)
                ENDDO
                current_data=SQRT(current_data)
             ELSE
                current_data=particle_data(l_direction(idim))
             ENDIF
             cell(idim)=NINT((current_data-range(idim,1))/dgrid(idim))+1
             IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) use_this=.FALSE.
          ENDDO
#ifdef PER_PARTICLE_WEIGHT
          part_weight=current%weight
#endif
          IF (use_this) data(cell(1),cell(2),cell(3))=data(cell(1),cell(2),cell(3))+part_weight !* real_space_area
       ENDIF
       current=>current%next
    ENDDO

    IF (need_reduce) THEN
       ALLOCATE(data2(1:resolution(1),1:resolution(2),1:resolution(3)))
       data2=0.0_num
       CALL MPI_ALLREDUCE(data,data2,resolution(1)*resolution(2)*resolution(3),mpireal,MPI_SUM,comm_new,errcode)
       data=data2
       DEALLOCATE(data2)
    ENDIF

    grid_name="Grid_"//TRIM(name)//"_"//TRIM(particle_species(species)%name)
    norm_grid_name="Norm_Grid_"//TRIM(name)//"_"//TRIM(particle_species(species)%name)
    var_name=TRIM(name)//"_"//TRIM(particle_species(species)%name)

    CALL cfd_write_3d_cartesian_grid(TRIM(grid_name),"Grid",grid1(1:global_resolution(1)),grid2(1:global_resolution(2))&
         ,grid3(1:global_resolution(3)),0)
    IF (use_offset_grid) THEN
       IF (parallel(1)) grid1=grid1-range(1,1)
       IF (parallel(2)) grid2=grid2-range(2,1)
       IF (parallel(3)) grid3=grid3-range(3,1)
    ENDIF
    CALL cfd_write_3d_cartesian_grid(TRIM(norm_grid_name),"Grid",grid1(1:global_resolution(1))/conv(1),grid2(1:global_resolution(2))/conv(2)&
         ,grid3(1:global_resolution(3))/conv(3),0)

    CALL cfd_write_3d_cartesian_variable_parallel(TRIM(var_name),"dist_fn",global_resolution,stagger,TRIM(norm_grid_name),"Grid"&
         ,data,type_new)
    CALL mpi_type_free(type_new,errcode)
    IF (need_reduce) &
         CALL MPI_COMM_FREE(comm_new,errcode)

    DEALLOCATE(data)
    DEALLOCATE(grid1,grid2,grid3)

  END SUBROUTINE general_3d_dist_fn

  SUBROUTINE general_2d_dist_fn(name,direction,range,resolution,species,restrictions,use_restrictions)

    CHARACTER(len=*),INTENT(IN) :: name
    INTEGER,DIMENSION(2),INTENT(IN) :: direction
    REAL(num),DIMENSION(2,2),INTENT(INOUT) :: range
    INTEGER,DIMENSION(2),INTENT(INOUT) :: resolution
    INTEGER,INTENT(IN) :: species
    REAL(num),DIMENSION(5,2),INTENT(IN) :: restrictions
    LOGICAL,DIMENSION(5),INTENT(IN) :: use_restrictions

    REAL(num),DIMENSION(:,:),ALLOCATABLE :: data, data2
    REAL(num),DIMENSION(:),ALLOCATABLE :: grid1,grid2
    LOGICAL,DIMENSION(2) :: parallel
    REAL(num),DIMENSION(2) :: dgrid
    REAL(num) :: current_data,temp_data
    INTEGER :: idim,idir
    LOGICAL,DIMENSION(2) :: calc_range
    LOGICAL :: calc_ranges

    LOGICAL :: use_x,use_y,need_reduce
    INTEGER,DIMENSION(2) :: start_local,global_resolution
    INTEGER :: color
    INTEGER :: comm_new, type_new

    LOGICAL,DIMENSION(2,5) :: use_direction
    LOGICAL,DIMENSION(2) :: calc_mod
    INTEGER,DIMENSION(2) :: p_count
    INTEGER,DIMENSION(2) :: l_direction
    REAL(num),DIMENSION(2) :: conv 
    INTEGER,DIMENSION(2) :: cell
    LOGICAL :: use_this
    REAL(num) :: real_space_area, part_weight


    TYPE(particle),POINTER :: current
    CHARACTER(len=string_length) :: grid_name, norm_grid_name, var_name
    REAL(num),DIMENSION(2) :: stagger=0.0_num
    REAL(num),DIMENSION(5) :: particle_data
    REAL(num) :: max_p_conv

    INTEGER :: ind


    use_x=.FALSE.
    use_y=.FALSE.
    need_reduce=.TRUE.
    color=0
    global_resolution=resolution
    parallel=.FALSE.
    start_local=1
    calc_range=.FALSE.
    calc_ranges=.FALSE.
    p_count=0
    use_direction=.FALSE.
    l_direction=0
    real_space_area=1.0_num


    DO idim=1,2
       IF (direction(idim) .EQ. DIR_X) THEN
          use_x=.TRUE.
          resolution(idim) = nx
          range(idim,1)=x_start_local
          range(idim,2)=x_end_local
          !IF (use_offset_grid) range(idim,:)=range(idim,:)-x_start
          start_local(idim)=cell_x_start(coordinates(2)+1)
          global_resolution(idim)=nx_global
          dgrid(idim)=dx
          parallel(idim)=.TRUE.
          use_direction(idim,1)=.TRUE.
          l_direction(idim)=1
          conv(idim)=MAX(length_x,length_y)
          CYCLE
       ENDIF
       IF (direction(idim) .EQ. DIR_Y)  THEN
          use_y=.TRUE.
          resolution(idim) = ny
          range(idim,1)=y_start_local
          range(idim,2)=y_end_local
          !IF (use_offset_grid) range(idim,:)=range(idim,:)-y_start
          start_local(idim)=cell_y_start(coordinates(1)+1)
          global_resolution(idim)=ny_global
          dgrid(idim)=dy
          parallel(idim)=.TRUE.
          use_direction(idim,2)=.TRUE.
          l_direction(idim)=2
          conv(idim)=MAX(length_x,length_y)
          CYCLE
       ENDIF
       conv(idim)=c*m0
       !If we're here then this must be a momentum space direction
       !So determine which momentum space directions are needed
       IF (range(idim,1) .EQ. range(idim,2)) THEN
          calc_range(idim) = .TRUE.
          calc_ranges=.TRUE.
       ENDIF
       IF (IAND(direction(idim),DIR_PX) .NE. 0) THEN
          use_direction(idim,3)=.TRUE.
          p_count(idim)=p_count(idim)+1
          l_direction(idim)=3
       ENDIF
       IF (IAND(direction(idim),DIR_PY) .NE. 0) THEN
          use_direction(idim,4)=.TRUE.
          p_count(idim)=p_count(idim)+1
          l_direction(idim)=4
       ENDIF
       IF (IAND(direction(idim),DIR_PZ) .NE. 0) THEN
          use_direction(idim,5)=.TRUE.
          p_count(idim)=p_count(idim)+1
          l_direction(idim)=5
       ENDIF
    ENDDO

    IF (.NOT. use_x) real_space_area=real_space_area*dx
    IF (.NOT. use_y) real_space_area=real_space_area*dy

    DO idim=1,2
       IF (p_count(idim) .GT. 1) calc_mod(idim)=.TRUE.
    ENDDO
    !Calculate ranges for directions where needed
    IF (calc_ranges) THEN
       DO idim=1,2
          IF (calc_range(idim)) THEN
             range(idim,1)=1.0e6_num
             range(idim,2)=-1.0e6_num
          ENDIF
       ENDDO
       current=>particle_species(species)%attached_list%head
       ind =0
       DO WHILE(ASSOCIATED(current))
          particle_data(1:2)=current%part_pos
          particle_data(3:5)=current%part_p
          DO idim=1,2
             IF (calc_range(idim)) THEN
                IF (calc_mod(idim)) THEN
                   DO idir=1,5
                      IF (use_direction(idim,idir)) current_data=current_data+particle_data(idir)
                   ENDDO
                   current_data=SQRT(current_data)
                ELSE
                   current_data=particle_data(l_direction(idim))
                ENDIF
                use_this=.TRUE.
                DO idir=1,5
                   IF (use_restrictions(idir) .AND. (particle_data(idir) .LT. restrictions(idir,1) &
                        .OR. particle_data(idir) .GT. restrictions(idir,2))) use_this=.FALSE.
                ENDDO
                IF (use_this) THEN
                   IF (current_data .LT. range(idim,1)) range(idim,1)=current_data
                   IF (current_data .GT. range(idim,2)) range(idim,2)=current_data
                ENDIF
             ENDIF
          ENDDO
          ind=ind+1
          current=>current%next
       ENDDO
    ENDIF

    max_p_conv=-10.0_num
    DO idim=1,2
       IF (.NOT. parallel(idim)) THEN
          !If not parallel then this is a momentum dimension
          CALL MPI_ALLREDUCE(range(idim,1),temp_data,1,mpireal,MPI_MIN,comm,errcode)
          range(idim,1)=temp_data
          CALL MPI_ALLREDUCE(range(idim,2),temp_data,1,mpireal,MPI_MAX,comm,errcode)
          range(idim,2)=temp_data

          !Calculate the maxmium range of a momentum direction
          IF (range(idim,2)-range(idim,1) .GT. max_p_conv .AND. .NOT. parallel(idim)) max_p_conv=range(idim,2)-range(idim,1)
       ENDIF
       !Fix so that if distribution function is zero then it picks an arbitrary scale in that direction
       IF (range(idim,1) .EQ. range(idim,2)) THEN
          range(idim,1)=-1.0_num
          range(idim,2)=1.0_num
       ENDIF
    ENDDO

    DO idim=1,2
       IF (.NOT. parallel(idim)) conv(idim)=max_p_conv
    ENDDO


    !Setup MPI
    IF (use_x .AND. use_y) need_reduce=.FALSE.
    IF (.NOT. use_y) color=color+coordinates(2) !If using x direction need to reduce only across all y
    IF (.NOT. use_x) color=color+nprocx*coordinates(1) !If using y direction need to reduce only across all x

    IF (need_reduce) THEN
       CALL MPI_COMM_SPLIT(comm,color,rank,comm_new,errcode)
    ELSE
       comm_new=MPI_COMM_NULL
    ENDIF

    type_new=create_2d_field_subtype(resolution,global_resolution,start_local)
    !Create grids
    DO idim=1,2
       IF (.NOT. parallel(idim)) dgrid(idim)=(range(idim,2)-range(idim,1))/REAL(resolution(idim)-1,num)
    ENDDO
    ALLOCATE(grid1(0:global_resolution(1)),grid2(0:global_resolution(2)))

    grid1(0)=-dgrid(1) + range(1,1)
    DO idir=1,global_resolution(1)
       grid1(idir)=grid1(idir-1)+dgrid(1)
    ENDDO

    grid2(0)=-dgrid(2) + range(2,1)
    DO idir=1,global_resolution(2)
       grid2(idir)=grid2(idir-1)+dgrid(2)
    ENDDO

    ALLOCATE(data(1:resolution(1),1:resolution(2)))
    data=0.0_num

    current=>particle_species(species)%attached_list%head
    part_weight=weight
    DO WHILE(ASSOCIATED(current))
       particle_data(1:2)=current%part_pos
       particle_data(3:5)=current%part_p
       use_this=.TRUE.
       DO idir=1,5
          IF (use_restrictions(idir) .AND. (particle_data(idir) .LT. restrictions(idir,1) &
               .OR. particle_data(idir) .GT. restrictions(idir,2))) use_this=.FALSE.
       ENDDO
       IF (use_this) THEN
          DO idim=1,2
             IF (calc_mod(idim)) THEN
                DO idir=1,5
                   IF (use_direction(idim,idir)) current_data=current_data+particle_data(idir)
                ENDDO
                current_data=SQRT(current_data)
             ELSE
                current_data=particle_data(l_direction(idim))
             ENDIF
             cell(idim)=NINT((current_data-range(idim,1))/dgrid(idim))+1
             IF (cell(idim) .LT. 1 .OR. cell(idim) .GT. resolution(idim)) use_this=.FALSE.
          ENDDO
#ifdef PER_PARTICLE_WEIGHT
          part_weight=current%weight
#endif
          IF (use_this)data(cell(1),cell(2))=data(cell(1),cell(2))+part_weight!*real_space_area
       ENDIF
       current=>current%next
    ENDDO


    IF (need_reduce) THEN
       ALLOCATE(data2(1:resolution(1),1:resolution(2)))
       data2=0.0_num
       CALL MPI_ALLREDUCE(data,data2,resolution(1)*resolution(2),mpireal,MPI_SUM,comm_new,errcode)
       data=data2
       DEALLOCATE(data2)
    ENDIF

    grid_name="Grid_"//TRIM(name)//"_"//TRIM(particle_species(species)%name)
    norm_grid_name="Norm_Grid_"//TRIM(name)//"_"//TRIM(particle_species(species)%name)
    var_name=TRIM(name)//"_"//TRIM(particle_species(species)%name)

    CALL cfd_write_2d_cartesian_grid(TRIM(grid_name),"Grid",grid1(1:global_resolution(1)),grid2(1:global_resolution(2))&
         ,0)
    IF (use_offset_grid) THEN
       IF (parallel(1)) grid1=grid1-range(1,1)
       IF (parallel(2)) grid2=grid2-range(2,1)
    ENDIF
    CALL cfd_write_2d_cartesian_grid(TRIM(norm_grid_name),"Grid",grid1(1:global_resolution(1))/conv(1),grid2(1:global_resolution(2))/conv(2)&
         ,0)

    CALL cfd_write_2d_cartesian_variable_parallel(TRIM(var_name),"dist_fn",global_resolution,stagger,TRIM(norm_grid_name),"Grid"&
         ,data,type_new)
    CALL mpi_type_free(type_new,errcode)
    IF (need_reduce) &
         CALL MPI_COMM_FREE(comm_new,errcode)

    DEALLOCATE(data)
    DEALLOCATE(grid1,grid2)

  END SUBROUTINE general_2d_dist_fn

  FUNCTION create_3d_field_subtype(n_local,n_global,start)
    INTEGER,DIMENSION(3),INTENT(IN) :: n_local
    INTEGER,DIMENSION(3),INTENT(IN) :: n_global
    INTEGER,DIMENSION(3),INTENT(IN) :: start
    INTEGER, DIMENSION(:),ALLOCATABLE :: lengths,starts
    INTEGER :: ipoint
    INTEGER :: create_3d_field_subtype
    ALLOCATE(lengths(1:n_local(2) * n_local(3)),starts(1:n_local(2) * n_local(3)))
    lengths=n_local(1)
    ipoint=0
    DO iz=0,n_local(3)-1
       DO iy=0,n_local(2)-1
          ipoint=ipoint+1
          starts(ipoint)=(start(3)+iz-1) * n_global(1) * n_global(2) + (start(2)+iy-1) * n_global(1) + start(1) -1
       ENDDO
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2)*n_local(3),lengths,starts,mpireal,create_3d_field_subtype,errcode)
    CALL MPI_TYPE_COMMIT(create_3d_field_subtype,errcode)
    DEALLOCATE(lengths,starts)
  END FUNCTION create_3d_field_subtype

  FUNCTION create_2d_field_subtype(n_local,n_global,start)
    INTEGER,DIMENSION(2),INTENT(IN) :: n_local
    INTEGER,DIMENSION(2),INTENT(IN) :: n_global
    INTEGER,DIMENSION(2),INTENT(IN) :: start
    INTEGER, DIMENSION(:),ALLOCATABLE :: lengths,starts
    INTEGER :: ipoint
    INTEGER :: create_2d_field_subtype
    ALLOCATE(lengths(1:n_local(2)),starts(1:n_local(2)))
    lengths=n_local(1)
    ipoint=0
    DO iy=0,n_local(2)-1
       ipoint=ipoint+1
       starts(ipoint)=(start(2)+iy-1) * n_global(1) + start(1) -1
    ENDDO

    CALL MPI_TYPE_INDEXED(n_local(2),lengths,starts,mpireal,create_2d_field_subtype,errcode)
    CALL MPI_TYPE_COMMIT(create_2d_field_subtype,errcode)
    DEALLOCATE(lengths,starts)
  END FUNCTION create_2d_field_subtype

END MODULE dist_fn
