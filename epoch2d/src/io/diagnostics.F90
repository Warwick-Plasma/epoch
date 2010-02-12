MODULE diagnostics

  USE shared_data 
  USE calc_df
  USE output_cartesian 
  USE output_particle
  USE output_arb
  USE output
  USE iocontrol
  USE balance
  USE dist_fn
  USE probes
  USE iterators
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
    REAL(num),DIMENSION(:,:),ALLOCATABLE :: Data
    REAL(num),DIMENSION(2) :: Stagger=0.0_num
    INTEGER(KIND=8) :: n_part_per_it=100000,npart_local,npart_dump_global
    INTEGER :: iSpecies,code
    INTEGER,DIMENSION(2) :: dims

    dims=(/nx_global,ny_global/)

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
       !If the code is doing a restart dump then tell the iterators that this is a restart dump
       IF (IAND(code,IO_RESTARTABLE) .NE. 0) THEN 
          Iterator_Settings%Restart=.TRUE.
       ELSE
          Iterator_Settings%Restart=.FALSE.
       ENDIF

       ALLOCATE(Data(-2:nx+3,-2:ny+3))
       !Open the file
       !(filename,rank_of_current_process,MPI_COMMUNICATOR (can be MPI_COMM_WORLD), MPI_FILE_MODE (passed straight to MPI_FILE_OPEN))
       CALL cfd_Open(filename,rank,comm,MPI_MODE_CREATE + MPI_MODE_WRONLY)
       !Write the snapshot information
       !If you prefer the VisIT cycles to display the dump number, change i for output_file
       !(code_time,n_iterations,rank to write)
       CALL cfd_Write_Snapshot_Data(time,i,0)

       IF (IAND(DumpMask(1),code) .NE. 0) CALL cfd_Write_nD_Particle_Grid_With_Iterator_All("Particles","Part_Grid",&
            Iterate_Particles,2,npart_local,npart_dump_global,npart_per_it,PARTICLE_CARTESIAN,subtype_particle_var)
       !Write the cartesian mesh
       !(Mesh_Name,Mesh_Class,x_array,y_array,rank to write)
       IF (IAND(DumpMask(2),code) .NE. 0) THEN
          IF (.NOT. Use_Offset_Grid) THEN
             CALL cfd_Write_2D_Cartesian_Grid("Grid","Grid",x_global(1:nx_global),y_global(1:ny_global),0)
          ELSE
             CALL cfd_Write_2D_Cartesian_Grid("Grid","Grid",x_offset_global(1:nx_global),y_offset_global(1:ny_global),0)      
             CALL cfd_Write_2D_Cartesian_Grid("Grid_Full","Grid",x_global(1:nx_global),y_global(1:ny_global),0)          
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

       !Serial Cartesian write because shared_grid
       !Parallel equivalent is cfd_Write_1D_Cartesian_Variable_All
       !(Variable_Name,Variable_Class,Grid_Stagger,Mesh_Name,Mesh_Class,Array,rank of process to write data)
       !Note that serial writes must be called on each node to keep file pointer updated
       IF (IAND(DumpMask(9),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Ex","Electric Field",Dims,Stagger,"Grid"&
            ,"Grid",Ex(1:nx,1:ny),subtype_field)
       IF (IAND(DumpMask(10),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Ey","Electric Field",Dims,Stagger,"Grid"&
            ,"Grid",Ey(1:nx,1:ny),subtype_field)
       IF (IAND(DumpMask(11),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Ez","Electric Field",Dims,Stagger,"Grid"&
            ,"Grid",Ez(1:nx,1:ny),subtype_field)

       IF (IAND(DumpMask(12),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Bx","Magnetic Field",Dims,Stagger,"Grid"&
            ,"Grid",Bx(1:nx,1:ny),subtype_field)
       IF (IAND(DumpMask(13),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("By","Magnetic Field",Dims,Stagger,"Grid"&
            ,"Grid",By(1:nx,1:ny),subtype_field)
       IF (IAND(DumpMask(14),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Bz","Magnetic Field",Dims,Stagger,"Grid"&
            ,"Grid",Bz(1:nx,1:ny),subtype_field)

       IF (IAND(DumpMask(15),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Jx","Current",Dims,Stagger,"Grid","Grid"&
            ,Jx(1:nx,1:ny),subtype_field)
       IF (IAND(DumpMask(16),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Jy","Current",Dims,Stagger,"Grid","Grid"&
            ,Jy(1:nx,1:ny),subtype_field)
       IF (IAND(DumpMask(17),code) .NE. 0) CALL cfd_Write_2D_Cartesian_Variable_Parallel("Jz","Current",Dims,Stagger,"Grid","Grid"&
            ,Jz(1:nx,1:ny),subtype_field)

       !Since these use species lookup tables, have to use the iterator functions
       !(Variable_Name,Variable_Class,Iterator_Function,global_npart,npart_per_iteration,Mesh_Name,Mesh_Class,MPI_TYPE describing data distribution)
       IF (IAND(DumpMask(18),code) .NE. 0)&
 			CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("Q","Particles",&
			iterate_charge,npart_dump_global,n_part_per_it&
			,"Particles","Part_Grid",subtype_particle_var)
       IF (IAND(DumpMask(19),code) .NE. 0)&
 			CALL cfd_Write_nD_Particle_Variable_With_Iterator_All("mass","Particles",&
			iterate_mass,npart_dump_global,n_part_per_it&
			,"Particles","Part_Grid",subtype_particle_var)

       IF (IAND(DumpMask(20),code) .NE. 0) THEN
          IF (IAND(DumpMask(20),IO_NO_INTRINSIC) .EQ. 0)  THEN
          	CALL calc_ekbar(Data,0)
				CALL cfd_Write_2D_Cartesian_Variable_Parallel("EkBar","EkBar",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
			ENDIF
          IF (IAND(DumpMask(20),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_ekbar(Data,iSpecies)
                WRITE(Temp_Name,'("EkBar_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_2D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"EkBar",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
             ENDDO
       ENDIF
		ENDIF
       !These are derived variables from the particles
       !Since you only dump after several particle updates it's actually quicker to
       IF (IAND(DumpMask(21),code) .NE. 0) THEN
          CALL calc_mass_density(Data,0)
          IF (IAND(DumpMask(21),IO_NO_INTRINSIC) .EQ. 0) &
				CALL cfd_Write_2D_Cartesian_Variable_Parallel("Mass_Density","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
          IF (IAND(DumpMask(21),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_mass_density(Data,iSpecies)
                WRITE(Temp_Name,'("Mass_Density_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_2D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
             ENDDO
          ENDIF
       ENDIF
       IF (IAND(DumpMask(22),code) .NE. 0) THEN
          CALL calc_charge_density(Data,0)
  			 IF (IAND(DumpMask(22),IO_NO_INTRINSIC) .EQ. 0) &	
          	CALL cfd_Write_2D_Cartesian_Variable_Parallel("Charge_Density","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
          IF (IAND(DumpMask(22),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_charge_density(Data,iSpecies)
                WRITE(Temp_Name,'("Charge_Density_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_2D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
             ENDDO
          ENDIF
       ENDIF

       IF (IAND(DumpMask(23),code) .NE. 0) THEN
          CALL calc_number_density(Data,0)
			 IF (IAND(DumpMask(23),IO_NO_INTRINSIC) .EQ. 0) &
          	CALL cfd_Write_2D_Cartesian_Variable_Parallel("Number_Density","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
          IF (IAND(DumpMask(23),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_number_density(Data,iSpecies)
                WRITE(Temp_Name,'("Number_Density_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_2D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
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
       CALL cfd_Write_2D_Cartesian_Variable_Parallel("Rank","Processor",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
#endif

       IF (IAND(dumpmask(26),code) .NE. 0) THEN
          CALL Write_Dist_Fns(code)
       ENDIF
#ifdef PARTICLE_PROBES
       IF (IAND(dumpmask(27),code) .NE. 0) THEN
          CALL Write_Probes(code)
       ENDIF
#endif
       IF (IAND(DumpMask(28),code) .NE. 0) THEN
          CALL calc_temperature(Data,0)
          CALL cfd_Write_2D_Cartesian_Variable_Parallel("Temperature","Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
          IF (IAND(DumpMask(28),IO_SPECIES) .NE. 0) THEN
             DO iSpecies=1,nspecies
                CALL calc_temperature(Data,iSpecies)
                WRITE(Temp_Name,'("Temperature_",a)') TRIM(ParticleSpecies(iSpecies)%Name)
                CALL cfd_Write_2D_Cartesian_Variable_Parallel(TRIM(ADJUSTL(Temp_Name)),"Derived",Dims,Stagger,"Grid","Grid",Data(1:nx,1:ny),subtype_field)
             ENDDO
          ENDIF
       ENDIF

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
	 INTEGER :: ioutput

    REAL(num), SAVE :: t1 = 0.0_num

    IF (restart) THEN
       t1 = time
       restart = .FALSE.
    END IF

    print_arrays = .FALSE.
    last_call = .FALSE.

	 DO ioutput = 1, num_vars_to_dump
	 	IF (FLOOR((time - t1)/dt) .LE. averaged_data(ioutput)%average_over_iterations) THEN
			CALL average_field(ioutput)
		ENDIF
	ENDDO

    IF (time >= t1) THEN
       print_arrays = .TRUE.
       t1 = t1 + dt_snapshots
    END IF

    IF (time >= t_end .OR. i == nsteps) THEN
       last_call = .TRUE.
       print_arrays = .TRUE.
    END IF

  END SUBROUTINE io_test

  SUBROUTINE average_field(ioutput)
		INTEGER,INTENT(IN) :: ioutput
		INTEGER :: count
		
		count=averaged_data(ioutput)%average_over_iterations
			
		SELECT CASE(ioutput)
			CASE (9)
				averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						ex/REAL(count,num)
			CASE(10)
					averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						ey/REAL(count,num)
			CASE(11)
					averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						ez/REAL(count,num)
			CASE (12)
				averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						bx/REAL(count,num)
			CASE(13)
					averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						by/REAL(count,num)
			CASE(14)
					averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						bz/REAL(count,num)
			CASE (15)
				averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						jx/REAL(count,num)
			CASE(16)
					averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						jy/REAL(count,num)
			CASE(17)
					averaged_data(ioutput)%data=averaged_data(ioutput)%data + &
						jz/REAL(count,num)										
			
		END SELECT
		
  END SUBROUTINE average_field

  SUBROUTINE set_dt()        ! sets CFL limited step

    REAL(num) :: dtx,dty
    INTEGER :: ix, iy
    dtx=dx/c
    dty=dy/c
    dt=dtx*dty/SQRT(dtx**2+dty**2)
    IF (dt_laser .NE. 0.0_num) dt=MIN(dt,dt_laser)
    IF (dt_plasma_frequency .NE. 0.0_num) dt=MIN(dt,dt_plasma_frequency)
#ifdef NEWTONIAN
    dtx=dx/max_part_v
    dty=dy/max_part_v
    dt=MIN(dt,dtx*dty/SQRT(dtx**2+dty**2))
#endif
    dt=dt_multiplier * dt 

  END SUBROUTINE set_dt

  SUBROUTINE energy_account()


  END SUBROUTINE energy_account




END MODULE diagnostics







