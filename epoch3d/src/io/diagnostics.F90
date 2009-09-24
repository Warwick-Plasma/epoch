MODULE diagnostics

  USE shared_data 
  USE calc_df
  USE output_cartesian 
  USE output_particle
  USE output_arb
  USE output
  USE iocontrol
  USE balance
  USE particlepointeradvance
  USE dist_fn
  USE probes
  USE mpi_subtype_control

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_dt, output_routines,iterate_charge

CONTAINS



  SUBROUTINE output_routines(i)   ! i=step index

    INTEGER, INTENT(in) :: i
    LOGICAL :: print_arrays,last_call
    CHARACTER(LEN=9+Data_Dir_Max_Length+n_zeros) :: filename,filenamedesc
    CHARACTER(LEN=50) :: Temp_Name
    REAL(num),DIMENSION(:,:,:),ALLOCATABLE :: Data
    REAL(num),DIMENSION(3) :: Stagger=0.0_num
    INTEGER(KIND=8) :: n_part_per_it=100000,npart_local,npart_dump_global
    INTEGER :: iSpecies,code
    INTEGER,DIMENSION(3) :: dims
    TYPE(Distribution_Function_Block),POINTER :: Current

    dims=(/nx_global,ny_global,nz_global/)

    CALL io_test(i,print_arrays,last_call)
    !Allows a maximum of 10^999 output dumps, should be enough for anyone (feel free to laugh when this isn't the case)
    WRITE(FileNameDesc,'("(''nfs:'',a,''/'',i",i3.3,".",i3.3,"''.cfd'')")'),n_Zeros,n_Zeros
    WRITE(filename, FileNameDesc) TRIM(data_dir), output_file
    IF (print_arrays) THEN
       !Always dump the variables with the "Every" attribute
       code=IO_ALWAYS
       !Only dump variables with the "FULL" attributre on full dump intervals
       IF (MOD(output_file,Full_Dump_Every) .EQ. 0)  code=IOR(code,IO_FULL)
       IF (MOD(output_file,Restart_Dump_Every) .EQ. 0 .AND. Restart_Dump_Every .GT. -1) code=IOR(code,IO_RESTARTABLE)
       IF (last_call .AND. force_final_to_be_restartable) code=IOR(code,IO_RESTARTABLE)

       npart_local=Get_Total_Local_Dumped_Particles(IAND(code,IO_RESTARTABLE) .NE. 0)
       CALL MPI_ALLREDUCE(npart_local,npart_dump_global,1,MPI_INTEGER8,MPI_SUM,comm,errcode)
       CALL Create_Subtypes(IAND(code,IO_RESTARTABLE) .NE. 0)
       ALLOCATE(Data(-2:nx+3,-2:ny+3,-2:nz+3))
       !Open the file
       !(filename,rank_of_current_process,MPI_COMMUNICATOR (can be MPI_COMM_WORLD), MPI_FILE_MODE (passed straight to MPI_FILE_OPEN))
       CALL cfd_Open(filename,rank,comm,MPI_MODE_CREATE + MPI_MODE_WRONLY)
       !Write the snapshot information
       !If you prefer the VisIT cycles to display the dump number, change i for output_file
       !(code_time,n_iterations,rank to write)
       CALL cfd_Write_Snapshot_Data(time,i,0)

       IF (IAND(DumpMask(1),code) .NE. 0) CALL cfd_Write_nD_Particle_Grid_With_Iterator_All("Particles","Part_Grid",&
            Iterate_Particles,3,npart_local,npart_dump_global,npart_per_it,PARTICLE_CARTESIAN,subtype_particle_var)
       !Write the cartesian mesh
       !(Mesh_Name,Mesh_Class,x_array,y_array,rank to write)
       IF (IAND(DumpMask(2),code) .NE. 0) THEN
          IF (.NOT. Use_Offset_Grid) THEN
             CALL cfd_Write_3D_Cartesian_Grid("Grid","Grid",x_global(1:nx_global),y_global(1:ny_global),z_global(1:nz_global),0)
          ELSE
             CALL cfd_Write_3D_Cartesian_Grid("Grid","Grid",x_offset_global(1:nx_global),y_offset_global(1:ny_global)&
                  ,z_offset_global(1:nz_global),0)      
             CALL cfd_Write_3D_Cartesian_Grid("Grid_Full","Grid",x_global(1:nx_global),y_global(1:ny_global)&
                  ,z_global(1:nz_global),0)          
          ENDIF
       ENDIF

       !(Variable_Name,Variable_Class,Array,global_npart,Mesh_Name,Mesh_Class,MPI_TYPE describing data distribution)
       IF (IAND(DumpMask(3),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Px","Particles"&
            ,iterate_px,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       IF (IAND(DumpMask(4),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Py","Particles"&
            ,iterate_py,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       IF (IAND(DumpMask(5),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Pz","Particles"&
            ,iterate_pz,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       IF (IAND(DumpMask(6),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Vx","Particles"&
            ,iterate_vx,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       IF (IAND(DumpMask(7),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Vy","Particles"&
            ,iterate_vy,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       IF (IAND(DumpMask(8),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Vz","Particles"&
            ,iterate_vz,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
#ifdef PART_DEBUG
       CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Processor","Particles"&
            ,iterate_processor,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Processor_at_t0","Particles"&
            ,iterate_processor0,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
#endif

       IF (IAND(DumpMask(9),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Ex","Electric Field",Dims,Stagger,"Grid"&
            ,"Grid",Ex(1:nx,1:ny,1:nz),subtype_field)
       IF (IAND(DumpMask(10),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Ey","Electric Field",Dims,Stagger,"Grid"&
            ,"Grid",Ey(1:nx,1:ny,1:nz),subtype_field)
       IF (IAND(DumpMask(11),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Ez","Electric Field",Dims,Stagger,"Grid"&
            ,"Grid",Ez(1:nx,1:ny,1:nz),subtype_field)

       IF (IAND(DumpMask(12),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Bx","Magnetic Field",Dims,Stagger,"Grid"&
            ,"Grid",Bx(1:nx,1:ny,1:nz),subtype_field)
       IF (IAND(DumpMask(13),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("By","Magnetic Field",Dims,Stagger,"Grid"&
            ,"Grid",By(1:nx,1:ny,1:nz),subtype_field)
       IF (IAND(DumpMask(14),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Bz","Magnetic Field",Dims,Stagger,"Grid"&
            ,"Grid",Bz(1:nx,1:ny,1:nz),subtype_field)

       IF (IAND(DumpMask(15),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Jx","Current",Dims,Stagger,"Grid","Grid"&
            ,Jx(1:nx,1:ny,1:nz),subtype_field)
       IF (IAND(DumpMask(16),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Jy","Current",Dims,Stagger,"Grid","Grid"&
            ,Jy(1:nx,1:ny,1:nz),subtype_field)
       IF (IAND(DumpMask(17),code) .NE. 0) CALL cfd_Write_3D_Cartesian_Variable_Parallel("Jz","Current",Dims,Stagger,"Grid","Grid"&
            ,Jz(1:nx,1:ny,1:nz),subtype_field)

       !Since these use species lookup tables, have to use the iterator functions
       !(Variable_Name,Variable_Class,Iterator_Function,global_npart,npart_per_iteration,Mesh_Name,Mesh_Class,MPI_TYPE describing data distribution)
       IF (IAND(DumpMask(18),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Q","Particles",iterate_charge,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       IF (IAND(DumpMask(19),code) .NE. 0) CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("mass","Particles",iterate_mass,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)

       IF (IAND(DumpMask(20),code) .NE. 0) THEN
          DO iSpecies=1,nspecies
             WRITE(Temp_Name,'("EkBar_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
             CALL cfd_Write_3D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"EkBar",Dims,Stagger,"Grid","Grid",ekbar(1:nx,1:ny,1:nz,iSpecies),subtype_field)
          ENDDO
       ENDIF
       !These are derived variables from the particles
       !Since you only dump after several particle updates it's actually quicker to
       IF (IAND(DumpMask(21),code) .NE. 0) THEN
          CALL calc_mass_density(Data,0)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Mass_Density","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny,1:nz),subtype_field)
          IF (IAND(DumpMask(21),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_mass_density(Data,iSpecies)
                WRITE(Temp_Name,'("Mass_Density_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_3D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny,1:nz),subtype_field)
             ENDDO
          ENDIF
       ENDIF
       IF (IAND(DumpMask(22),code) .NE. 0) THEN
          CALL calc_charge_density(Data,0)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Charge_Density","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny,1:nz),subtype_field)
          IF (IAND(DumpMask(22),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_charge_density(Data,iSpecies)
                WRITE(Temp_Name,'("Charge_Density_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_3D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny,1:nz),subtype_field)
             ENDDO
          ENDIF
       ENDIF

       IF (IAND(DumpMask(23),code) .NE. 0) THEN
          CALL calc_number_density(Data,0)
          CALL cfd_Write_3D_Cartesian_Variable_Parallel("Number_Density","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny,1:nz),subtype_field)
          IF (IAND(DumpMask(23),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_number_density(Data,iSpecies)
                WRITE(Temp_Name,'("Number_Density_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_3D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny,1:nz),subtype_field)
             ENDDO
          ENDIF
       ENDIF

       IF (IAND(dumpmask(24),code) .NE. 0) THEN
#ifdef PER_PARTICLE_WEIGHT
          CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Weight","Particles",iterate_weight,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
#else
          CALL cfd_Write_Real_Constant("Weight","Particles",weight,0)
#endif
       ENDIF

       IF (IAND(dumpmask(25),code) .NE. 0) THEN
          CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Species","Particles",iterate_species,npart_dump_global,n_part_per_it,"Particles","Part_Grid",subtype_particle_var)
       ENDIF

#ifdef FIELD_DEBUG
       Data=rank
       CALL cfd_Write_3D_Cartesian_Variable_Parallel("Rank","Processor",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny,1:nz),subtype_field)
#endif

       IF (IAND(dumpmask(26),code) .NE. 0) THEN
          CALL Write_Dist_Fns(code)
       ENDIF
#ifdef PARTICLE_PROBES
       IF (IAND(dumpmask(27),code) .NE. 0) THEN
          CALL Write_Probes(code)
       ENDIF
#endif

!!$       IF (IAND(DumpMask(28),code) .NE. 0) THEN
!!$          CALL calc_temperature(Data,0)
!!$          CALL cfd_Write_2D_Cartesian_Variable_Parallel("Temperature","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
!!$          IF (IAND(DumpMask(28),IO_SPECIES) .NE. 0) THEN
!!$             DO iSpecies=1,nspecies
!!$                CALL calc_temperature(Data,iSpecies)
!!$                WRITE(Temp_Name,'("Temperature_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
!!$                CALL cfd_Write_2D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
!!$             ENDDO
!!$          ENDIF
!!$       ENDIF

       !Close the file
       CALL cfd_Close()

       output_file = output_file + 1
       IF (rank .EQ. 0) THEN
          WRITE(20,*) "Dumped data at",time,"at iteration",i,"for dump",output_file-1
          CALL FLUSH(20)
       ENDIF

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

    REAL(num) :: dtx, dty, dtz, dt0
    INTEGER :: ix, iy, iz

    dtx=dx/c
    dty=dy/c
    dtz=dz/c
    dt=MIN(dtx**2,dty**2,dtz**2)/SQRT(dtx**2+dty**2+dtz**2)
	 IF (dt_plasma .NE. 0.0_num) dt=MIN(dt,dt_plasma)
    IF (dt_laser .NE. 0.0_num) dt=MIN(dt,dt_laser)
    dt=dt_multiplier * dt 

  END SUBROUTINE set_dt

  SUBROUTINE energy_account()


  END SUBROUTINE energy_account

  SUBROUTINE Write_Species(filehandle,current_displacement)

    INTEGER,INTENT(IN) :: filehandle
    INTEGER(KIND=MPI_OFFSET_KIND),INTENT(IN) :: current_displacement

!!$    CALL MPI_FILE_SET_VIEW(filehandle, current_displacement, MPI_INTEGER, subtype_particle_int,&
!!$         "native", MPI_INFO_NULL,cfd_errcode)
!!$    CALL MPI_FILE_WRITE_ALL(filehandle, Part_Species,npart , MPI_INTEGER, status, errcode)

  END SUBROUTINE Write_Species

  !Iterator for particle positions
  SUBROUTINE iterate_particles(data,n_points,direction,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: Data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL,INTENT(IN) :: start
    INTEGER,INTENT(IN) :: direction
    TYPE(Particle),POINTER,SAVE :: Cur
    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily
    INTEGER(8) :: partcount


    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             Data(partcount)=Cur%Part_Pos(direction)-Window_Shift(direction)
             !IF (Cur%Part_Pos(1) .EQ. Cur%Part_Pos(2)) PRINT *,"PATBAD"
             Cur=>Cur%Next
          ENDDO
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
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
    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             data(partcount) = Cur%Charge
#else
             data(partcount) = CurrentFamily%Charge
#endif
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
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

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=Cur%Weight
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
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

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             data(partcount) = Cur%Mass
#else
             data(partcount) = CurrentFamily%Mass
#endif
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_mass

#ifdef PART_DEBUG
  !Iterator for particle processor
  SUBROUTINE iterate_processor(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=REAL(cur%Processor,num)
             IF (cur%Processor .GE. nproc) PRINT *,"Bad Processor"
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_processor

  SUBROUTINE iterate_processor0(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=REAL(cur%Processor_at_t0,num)
             IF (cur%Processor .GE. nproc) PRINT *,"Bad Processor"
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_processor0
#endif

  !Iterator for particle processor
  SUBROUTINE iterate_species(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount)=REAL(CurrentFamily%ID,num)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_species

  !Iterator for particle velocities
  SUBROUTINE iterate_vx(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root,part_m

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             part_m = Cur%Mass
#else
             part_m = CurrentFamily%Mass
#endif
             root=SQRT(part_m**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
             IF (root .NE. 0.0_num) root=1.0_num/root
             data(partcount) = Cur%Part_p(1) * root
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_vx

  SUBROUTINE iterate_vy(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root,part_m

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             part_m = Cur%Mass
#else
             part_m = CurrentFamily%Mass
#endif
             root=SQRT(part_m**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
             IF (root .NE. 0.0_num) root=1.0_num/root
             data(partcount) = Cur%Part_p(2) * root
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_vy

  SUBROUTINE iterate_vz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount
    REAL(num) :: root,part_m
    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
#ifdef PER_PARTICLE_CHARGEMASS
             part_m = Cur%Mass
#else
             part_m = CurrentFamily%Mass
#endif
             root=SQRT(part_m**2 + (Cur%part_p(1)**2 + Cur%part_p(2)**2 + Cur%part_p(3)**2)/c**2)
             IF (root .NE. 0.0_num) root=1.0_num/root
             data(partcount) = Cur%Part_p(3) * root
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
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

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount) = Cur%Part_p(1)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_px

  SUBROUTINE iterate_py(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount) = Cur%Part_p(2)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_py

  SUBROUTINE iterate_pz(data,n_points,start)

    REAL(num),DIMENSION(:),INTENT(INOUT) :: data
    INTEGER(8),INTENT(INOUT) :: n_points
    LOGICAL, INTENT(IN) :: start
    TYPE(Particle),POINTER,SAVE :: Cur
    INTEGER(8) :: partcount

    TYPE(ParticleList),POINTER,SAVE :: CurrentList
    TYPE(ParticleFamily),POINTER,SAVE :: CurrentFamily

    IF (start)  THEN
       CALL Start_ParticleFamily(CurrentFamily,CurrentList,Cur)
    ENDIF
    partcount=0
    DO WHILE (ASSOCIATED(CurrentFamily) .AND. (partcount .LT. n_points))
       DO WHILE (ASSOCIATED(CurrentList) .AND. (partcount .LT. n_points))
          DO WHILE (ASSOCIATED(Cur) .AND. (partcount .LT. n_points))
             partcount=partcount+1
             data(partcount) = Cur%Part_p(3)
             Cur=>Cur%Next
          ENDDO
          !If the current partlist is exhausted, switch to the next one
          IF (.NOT. ASSOCIATED(Cur)) CALL Advance_ParticleList(CurrentList,Cur)
       ENDDO
       !If the current particlefamily is exhausted, then switch to the next one
       DO WHILE (.NOT. ASSOCIATED(Cur))
          CALL Advance_ParticleFamily(CurrentFamily,CurrentList,Cur)
          IF (.NOT. ASSOCIATED(CurrentFamily)) EXIT
          IF (.NOT. CurrentFamily%Dump) NULLIFY(Cur)
       ENDDO
    ENDDO
    n_points=partcount

  END SUBROUTINE iterate_pz


END MODULE diagnostics







