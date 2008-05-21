MODULE diagnostics

  USE shared_data 
  USE calc_df
  USE output_cartesian 
  USE output_particle
  USE output_arb
  USE output
  USE iocontrol

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines,iterate_charge

CONTAINS

  SUBROUTINE iterate_charge(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER(8),SAVE :: curpos
    INTEGER(8) :: npart_this_it

    IF (start) curpos=1

    npart_this_it=MIN(n_points,npart-curpos+1)

    DO ipart=curpos,curpos+npart_this_it-1
       data(ipart-curpos+1) = species(part_species(ipart),1)
    ENDDO
    curpos=curpos + npart_this_it
    n_points=npart_this_it

  END SUBROUTINE iterate_charge

  SUBROUTINE iterate_mass(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER(8),SAVE :: curpos
    INTEGER(8) :: npart_this_it

    IF (start) curpos=1

    npart_this_it=MIN(n_points,npart-curpos+1)

    IF (npart_this_it == 0) THEN
       n_points=0
       RETURN
    ENDIF

    DO ipart=curpos,curpos+npart_this_it-1
       data(ipart-curpos+1) = species(part_species(ipart),2)
    ENDDO

    curpos=curpos + npart_this_it
    n_points=npart_this_it

  END SUBROUTINE iterate_mass

  SUBROUTINE iterate_vx(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER(8),SAVE :: curpos
    INTEGER(8) :: npart_this_it
    REAL(num) :: root

    IF (start) curpos=1

    npart_this_it=MIN(n_points,npart-curpos+1)

    IF (npart_this_it == 0) THEN
       n_points=0
       RETURN
    ENDIF

    DO ipart=curpos,curpos+npart_this_it-1
       root=1.0_num/SQRT(species(part_species(ipart),2)**2 + (part_p(ipart,1)**2 + part_p(ipart,2)**2 + part_p(ipart,3)**2)/c**2)
       data(ipart-curpos+1) = part_p(ipart,1) * root
    ENDDO

    curpos=curpos + npart_this_it
    n_points=npart_this_it

  END SUBROUTINE iterate_vx

  SUBROUTINE iterate_vy(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER(8),SAVE :: curpos
    INTEGER(8) :: npart_this_it
    REAL(num) :: root

    IF (start) curpos=1

    npart_this_it=MIN(n_points,npart-curpos+1)

    IF (npart_this_it == 0) THEN
       n_points=0
       RETURN
    ENDIF

    DO ipart=curpos,curpos+npart_this_it-1
       root=1.0_num/SQRT(species(part_species(ipart),2)**2 + (part_p(ipart,1)**2 + part_p(ipart,2)**2 + part_p(ipart,3)**2)/c**2)
       data(ipart-curpos+1) = part_p(ipart,2) * root
    ENDDO

    curpos=curpos + npart_this_it
    n_points=npart_this_it

  END SUBROUTINE iterate_vy


  SUBROUTINE iterate_vz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    INTEGER(8),SAVE :: curpos
    INTEGER(8) :: npart_this_it
    REAL(num) :: root

    IF (start) curpos=1

    npart_this_it=MIN(n_points,npart-curpos+1)

    IF (npart_this_it == 0) THEN
       n_points=0
       RETURN
    ENDIF

    DO ipart=curpos,curpos+npart_this_it-1
       root=1.0_num/SQRT(species(part_species(ipart),2)**2 + (part_p(ipart,1)**2 + part_p(ipart,2)**2 + part_p(ipart,3)**2)/c**2)
       data(ipart-curpos+1) = part_p(ipart,3) * root
    ENDDO

    curpos=curpos + npart_this_it
    n_points=npart_this_it

  END SUBROUTINE iterate_vz

  SUBROUTINE output_routines(i)   ! i=step index

    INTEGER, INTENT(in) :: i
    LOGICAL :: print_arrays,last_call
    CHARACTER(LEN=13+Data_Dir_Max_Length) :: filename
    CHARACTER(LEN=50) :: Temp_Name
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Data
    REAL(num),DIMENSION(2) :: Stagger=0.0_num
    INTEGER(KIND=8) :: n_part_per_it=1000000
    INTEGER :: iSpecies,code


    CALL io_test(i,print_arrays,last_call)

    WRITE(filename, '("nfs:",a,"/",i4.4,".cfd")') TRIM(data_dir), output_file

    IF (print_arrays) THEN
       ALLOCATE(Data(-2:nx+3,-2:ny+3))
       !Open the file
       !(filename,rank_of_current_process,MPI_COMMUNICATOR (can be MPI_COMM_WORLD), MPI_FILE_MODE (passed straight to MPI_FILE_OPEN))
       CALL cfd_Open(filename,rank,comm,MPI_MODE_CREATE + MPI_MODE_WRONLY)
       !Write the snapshot information
       !If you prefer the VisIT cycles to display the dump number, change i for output_file
       !(code_time,n_iterations,rank to write)
       CALL cfd_Write_Snapshot_Data(time,i,0)

       !Always dump the variables with the "Every" attribute
       code=IO_ALWAYS
       !Only dump variables with the "FULL" attributre on full dump intervals
       IF (MOD(output_file,Full_Dump_Every) .EQ. 0)  code=IOR(code,IO_FULL)

       IF (IAND(DumpMask(1),code) .NE. 0) CALL cfd_Write_nD_Particle_Grid_All("Particles","Grid",part_pos,npart_global,PARTICLE_CARTESIAN,subtype_particle_mesh)
       !Write the cartesian mesh
       !(Mesh_Name,Mesh_Class,x_array,y_array,rank to write)
       If (IAND(DumpMask(2),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Grid("Grid","Grid",x(1:nx),y(1:ny),0)
       !(Variable_Name,Variable_Class,Array,global_npart,Mesh_Name,Mesh_Class,MPI_TYPE describing data distribution)
       IF (IAND(DumpMask(3),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_All("P_x","Particles",part_p(:,1)&
            ,npart_global,"Particles","Grid",subtype)
       IF (IAND(DumpMask(4),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_All("P_y","Particles",part_p(:,2)&
            ,npart_global,"Particles","Grid",subtype)
       IF (IAND(DumpMask(5),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_All("P_z","Particles",part_p(:,3)&
            ,npart_global,"Particles","Grid",subtype)
       IF (IAND(DumpMask(6),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("V_x","Particles"&
            ,iterate_vx,npart_global,n_part_per_it,"Particles","Grid",subtype)
       IF (IAND(DumpMask(7),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("V_y","Particles"&
            ,iterate_vy,npart_global,n_part_per_it,"Particles","Grid",subtype)
       IF (IAND(DumpMask(8),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("V_z","Particles"&
            ,iterate_vz,npart_global,n_part_per_it,"Particles","Grid",subtype)

       !Serial Cartesian write because shared_grid
       !Parallel equivalent is cfd_Write_1D_Cartesian_Variable_All
       !(Variable_Name,Variable_Class,Grid_Stagger,Mesh_Name,Mesh_Class,Array,rank of process to write data)
       !Note that serial writes must be called on each node to keep file pointer updated
       IF (IAND(DumpMask(9),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Ex","Electric Field",Stagger,"Grid"&
            ,"Grid",Ex(1:nx,1:ny),0)
       IF (IAND(DumpMask(10),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Ey","Electric Field",Stagger,"Grid"&
            ,"Grid",Ey(1:nx,1:ny),0)
       IF (IAND(DumpMask(11),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Ez","Electric Field",Stagger,"Grid"&
            ,"Grid",Ez(1:nx,1:ny),0)

       IF (IAND(DumpMask(12),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Bx","Magnetic Field",Stagger,"Grid"&
            ,"Grid",Bx(1:nx,1:ny),0)
       IF (IAND(DumpMask(13),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("By","Magnetic Field",Stagger,"Grid"&
            ,"Grid",By(1:nx,1:ny),0)
       IF (IAND(DumpMask(14),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Bz","Magnetic Field",Stagger,"Grid"&
            ,"Grid",Bz(1:nx,1:ny),0)

       IF (IAND(DumpMask(15),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Jx","Current",Stagger,"Grid","Grid"&
            ,Jx(1:nx,1:ny),0)
       IF (IAND(DumpMask(16),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Jy","Current",Stagger,"Grid","Grid"&
            ,Jy(1:nx,1:ny),0)
       IF (IAND(DumpMask(17),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable("Jz","Current",Stagger,"Grid","Grid"&
            ,Jz(1:nx,1:ny),0)

       !Since these use species lookup tables, have to use the iterator functions
       !(Variable_Name,Variable_Class,Iterator_Function,global_npart,npart_per_iteration,Mesh_Name,Mesh_Class,MPI_TYPE describing data distribution)
       IF (IAND(DumpMask(18),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Q","Particles",iterate_charge,npart_global,n_part_per_it,"Particles","Grid",subtype)
       IF (IAND(DumpMask(19),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("mass","Particles",iterate_mass,npart_global,n_part_per_it,"Particles","Grid",subtype)

       IF (IAND(DumpMask(20),code) .NE. 0) THEN
          DO iSpecies=1,nspecies
             WRITE(Temp_Name,'("EkBar_",a)') TRIM(Species_Name(iSpecies)%Value)
             CALL cfd_Write_2D_Cartesian_Variable(TRIM(ADJUSTL(Temp_Name)),"EkBar",Stagger,"Grid","Grid",temperature(1:nx,1:ny,iSpecies),0)
          ENDDO
       ENDIF
       !These are derived variables from the particles
       !Since you only dump after several particle updates it's actually quicker to
       IF (IAND(DumpMask(21),code) .NE. 0) THEN
          CALL calc_mass_density(Data)
          CALL cfd_Write_2D_Cartesian_Variable("Mass_Density","Derived",Stagger,"Grid","Grid",Data(1:nx,1:ny),0)
       ENDIF
       IF (IAND(DumpMask(22),code) .NE. 0) THEN
          CALL calc_charge_density(Data)
          CALL cfd_Write_2D_Cartesian_Variable("Charge_Density","Derived",Stagger,"Grid","Grid",Data(1:nx,1:ny),0)
       ENDIF

       IF (IAND(dumpmask(23),code) .NE. 0) THEN
          CALL cfd_Write_Real_Constant("Weight","Particles",weight,0)
       ENDIF

       IF (IAND(dumpmask(24),code) .NE. 0) THEN
          CALL cfd_Write_Arb_Block("Species","Species","PIC2D",npart_global*4,Write_Species)
       ENDIF


       !Close the file
       CALL cfd_Close()

       output_file = output_file + 1
       IF (rank .EQ. 0) WRITE(20,*) "Dumped data at",time,"at iteration",i,"for dump",output_file-1
       DEALLOCATE(Data)
    ENDIF


  END SUBROUTINE output_routines

  SUBROUTINE io_test(i, print_arrays, last_call)

    INTEGER, INTENT(IN) :: i
    LOGICAL, INTENT(OUT) :: print_arrays, last_call

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
       t1 = time
       restart = .FALSE.
    END IF

    print_arrays = .FALSE.
    last_call = .FALSE.

    IF (time >= t1) THEN
       print_arrays = .TRUE.
       t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. i == nsteps) THEN
       last_call = .TRUE.
       print_arrays = .TRUE.
    END IF

  END SUBROUTINE io_test

  SUBROUTINE set_dt        ! sets CFL limited step

    REAL(num) :: cons, dt1, dt2, a, scale,dt_global,dt_p,root
    INTEGER :: ix, iy

    dt=MIN(dx/(c*SQRT(2.0_num)),dy/(c * SQRT(2.0)),2.0_num*pi/las_omega)
    dt=dt_multiplier * dt 

  END SUBROUTINE set_dt

  SUBROUTINE energy_account()


  END SUBROUTINE energy_account

  SUBROUTINE Write_Species(filehandle,current_displacement)

    INTEGER,INTENT(IN) :: filehandle
    INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement

    CALL MPI_FILE_SET_VIEW(filehandle, current_displacement, MPI_INTEGER, subtype_int,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(filehandle, Part_Species,npart , MPI_INTEGER, status, errcode)

  END SUBROUTINE Write_Species

END MODULE diagnostics







