MODULE diagnostics

  USE shared_data 
  USE calc_df
  USE output_cartesian 
  USE output_particle
  USE output_arb
  USE output
  USE iocontrol
  USE balance

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines,iterate_charge

CONTAINS

  SUBROUTINE output_routines(i)   ! i=step index

    INTEGER, INTENT(in) :: i
    LOGICAL :: print_arrays,last_call
    CHARACTER(LEN=13+Data_Dir_Max_Length) :: filename
    CHARACTER(LEN=50) :: Temp_Name
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: Data
    REAL(num),DIMENSION(3) :: Stagger=0.0_num
    INTEGER(KIND=8) :: n_part_per_it=1000000
    INTEGER :: iSpecies,code
    INTEGER,DIMENSION(3) :: dims

    dims=(/nx_global,ny_global,nz_global/)


    CALL io_test(i,print_arrays,last_call)

    WRITE(filename, '("nfs:",a,"/",i4.4,".cfd")') TRIM(data_dir), output_file

    IF (print_arrays) THEN
       CALL CreateSubtypes
       ALLOCATE(Data(1:nx,1:ny,1:nz))
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

       IF (IAND(DumpMask(1),code) .NE. 0) CALL cfd_Write_nD_Particle_Grid_With_Iterator_All("Particles","Grid",Iterate_Particles,3,npart,npart_global,n_part_per_it,PARTICLE_CARTESIAN,subtype_particle_var)
       !Write the cartesian mesh
       !(Mesh_Name,Mesh_Class,x_array,y_array,rank to write)
       IF (IAND(DumpMask(2),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Grid("Grid","Grid",x_global(1:nx_global),y_global(1:ny_global),z_global(1:nz_global),0)
       !(Variable_Name,Variable_Class,Array,global_npart,Mesh_Name,Mesh_Class,MPI_TYPE describing data distribution)
       IF (IAND(dumpmask(3),code) .NE. 0)CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Px","Particles",&
            iterate_px,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)
       IF (IAND(dumpmask(4),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Py","Particles",&
            iterate_py,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)
       IF (IAND(dumpmask(5),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Pz","Particles",&
            iterate_pz,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)

       IF (IAND(dumpmask(6),code) .NE. 0)CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Vx","Particles",&
            iterate_vx,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)
       IF (IAND(dumpmask(7),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Vy","Particles",&
            iterate_vy,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)
       IF (IAND(dumpmask(8),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Vz","Particles",&
            iterate_vz,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)


       !Serial Cartesian write because shared_grid
       !Parallel equivalent is cfd_Write_1D_Cartesian_Variable_All
       !(Variable_Name,Variable_Class,Grid_Stagger,Mesh_Name,Mesh_Class,Array,rank of process to write data)
       !Note that serial writes must be called on each node to keep file pointer updated
       IF (IAND(DumpMask(9),code) .NE. 0) THEN
          Data=Ex(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Ex","Electric Field",Dims,Stagger,"Grid"&
               ,"Grid",Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(10),code) .NE. 0) THEN
          Data=Ey(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Ey","Electric Field",Dims,Stagger,"Grid"&
               ,"Grid",Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(11),code) .NE. 0) THEN
          Data=Ez(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Ez","Electric Field",Dims,Stagger,"Grid"&
               ,"Grid",Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(12),code) .NE. 0) THEN
          Data=Bx(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Bx","Magnetic Field",Dims,Stagger,"Grid"&
               ,"Grid",Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(13),code) .NE. 0) THEN
          Data=By(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("By","Magnetic Field",Dims,Stagger,"Grid"&
               ,"Grid",Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(14),code) .NE. 0) THEN
          Data=Bz(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Bz","Magnetic Field",Dims,Stagger,"Grid"&
               ,"Grid",Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(15),code) .NE. 0) THEN
          Data=Jx(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Jx","Current",Dims,Stagger,"Grid","Grid"&
               ,Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(16),code) .NE. 0) THEN
          Data=Jy(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Jy","Current",Dims,Stagger,"Grid","Grid"&
               ,Data,subtype_field)
       ENDIF
       IF (IAND(DumpMask(17),code) .NE. 0) THEN
          Data=Jz(1:nx,1:ny,1:nz)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Jz","Current",Dims,Stagger,"Grid","Grid"&
               ,Data,subtype_field)
       ENDIF
       !Since these use species lookup tables, have to use the iterator functions
       !(Variable_Name,Variable_Class,Iterator_Function,global_npart,npart_per_iteration,Mesh_Name,Mesh_Class,MPI_TYPE describing data distribution)
       IF (IAND(DumpMask(18),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Q","Particles",iterate_charge,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)
       IF (IAND(DumpMask(19),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("mass","Particles",iterate_mass,npart_global,n_part_per_it,"Particles","Grid",subtype_particle_var)

       IF (IAND(DumpMask(20),code) .NE. 0) THEN
          DO iSpecies=1,nspecies
             WRITE(Temp_Name,'("EkBar_",a)') TRIM(Species_Name(iSpecies)%Value)
             Data=temperature(1:nx,1:ny,1:nz,iSpecies)
             CALL cfd_Write_3D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"EkBar",Dims,Stagger,"Grid","Grid",Data,subtype_field)
          ENDDO
       ENDIF
       !These are derived variables from the particles
       !Since you only dump after several particle updates it's actually quicker to
!!$       IF (IAND(DumpMask(21),code) .NE. 0) THEN
!!$          CALL calc_mass_density(Data)
!!$          CALL cfd_Write_2D_Cartesian_Variable("Mass_Density","Derived",Stagger,"Grid","Grid",Data(1:nx,1:ny),0)
!!$       ENDIF
!!$       IF (IAND(DumpMask(22),code) .NE. 0) THEN
!!$          CALL calc_charge_density(Data)
!!$          CALL cfd_Write_2D_Cartesian_Variable("Charge_Density","Derived",Stagger,"Grid","Grid",Data(1:nx,1:ny),0)
!!$       ENDIF

       IF (IAND(dumpmask(23),code) .NE. 0) THEN
          CALL cfd_Write_Real_Constant("Weight","Particles",weight,0)
       ENDIF

!!$       IF (IAND(dumpmask(24),code) .NE. 0) THEN
!!$          CALL cfd_Write_Arb_Block("Species","Species","PIC2D",npart_global*4,Write_Species)
!!$       ENDIF


       !Close the file
       CALL cfd_Close()

       output_file = output_file + 1
       IF (rank .EQ. 0) WRITE(20,*) "Dumped data at",time,"at iteration",i,"for dump",output_file-1
       CALL FLUSH(20)
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

    dt=MIN(dx/(c*SQRT(3.0_num)),dy/(c * SQRT(3.0_num)),dz/(c * SQRT(3.0_num)),2.0_num*pi/las_omega)
    dt=dt_multiplier * dt 

  END SUBROUTINE set_dt

  SUBROUTINE energy_account()


  END SUBROUTINE energy_account

  SUBROUTINE Write_Species(filehandle,current_displacement)

    INTEGER,INTENT(IN) :: filehandle
    INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement

    CALL MPI_FILE_SET_VIEW(filehandle, current_displacement, MPI_INTEGER, subtype_particle_int,&
         "native", MPI_INFO_NULL,cfd_errcode)
    CALL MPI_FILE_WRITE_ALL(filehandle, Part_Species,npart , MPI_INTEGER, status, errcode)

  END SUBROUTINE Write_Species

  !Iterator for particle positions
  SUBROUTINE iterate_particles(data,n_points,direction,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(8),INTENT(INOUT) :: n_points
    INTEGER, INTENT(IN) :: direction
    LOGICAL,INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start) Cur=>Head
    partcount=0
    DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       Data(partcount)=Cur%Part_Pos(direction)
       Cur=>Cur%Next
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_particles

  !Iterator for particle charge
  SUBROUTINE iterate_charge(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
       data(partcount) = Cur%Charge
#else
       data(partcount) = species(Cur%Part_Species,1)
#endif
       Cur=>Cur%Next
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_charge

#ifdef PER_PARTICLE_WEIGHT

  !Iterator for particle weight
  !Only present if you are using the PER_PARTICLE_WEIGHT
  !Precompiler option
  SUBROUTINE iterate_weight(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Weight
       Cur=>Cur%Next
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_weight
#endif

  !Iterator for particle mass
  SUBROUTINE iterate_mass(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
       data(partcount) = Cur%Mass
#else
       data(partcount) = species(Cur%Part_Species,2)
#endif
       Cur=>Cur%Next
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_mass

  !Iterator for particle velocities
  SUBROUTINE iterate_vx(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       root=1.0_num/SQRT(species(Cur%part_species,2)**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
       data(partcount) = Cur%Part_p(1) * root
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_vx

  SUBROUTINE iterate_vy(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       root=1.0_num/SQRT(species(Cur%part_species,2)**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
       data(partcount) = Cur%Part_p(2) * root
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_vy

  SUBROUTINE iterate_vz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       root=1.0_num/SQRT(species(Cur%part_species,2)**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
       data(partcount) = Cur%Part_p(3) * root
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_vz

  !Iterator for particle momenta
  SUBROUTINE iterate_px(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(1)
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_px

  SUBROUTINE iterate_py(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(2)
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_py

  SUBROUTINE iterate_pz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    IF (start) Cur=>Head
    partcount=0
    DO WHILE(ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
       partcount=partcount+1
       data(partcount) = Cur%Part_p(3)
       Cur=>Cur%Next
    ENDDO

    n_points=partcount

  END SUBROUTINE iterate_pz

END MODULE diagnostics







