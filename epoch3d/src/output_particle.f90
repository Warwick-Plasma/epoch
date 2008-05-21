MODULE output_particle

  USE shared_data
  USE iocommon
  USE output

  IMPLICIT NONE

CONTAINS

  !-------------------------------------------------------------------------------------------------------------------
  !Code to write a 2D Cartesian grid in serial from the node with rank {rank_write}
  !Serial operation, so no need to specify nx,ny
  !-------------------------------------------------------------------------------------------------------------------
  SUBROUTINE cfd_Write_nD_Particle_Grid_All(name,class,particles,npart_global,Particle_coord_type,Particle_Type)

    REAL(num),DIMENSION(:,:),INTENT(IN) :: particles
    CHARACTER(len=*), INTENT(IN) :: name, class
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(4), INTENT(IN) :: Particle_coord_type
    INTEGER, INTENT(IN) :: Particle_Type
    INTEGER(8) :: npart_local
    INTEGER(8) :: blocklen,mdlen
    INTEGER(4) :: ndim,i,disp0
    INTEGER(4) :: sizes(2)
    REAL(num) :: mn,mx


    sizes = SHAPE(particles)
    npart_local = sizes(1)
    ndim = sizes(2)


    !Metadata is
    !* ) MeshType (INTEGER(4)) All mesh blocks contain this
    !* ) nd    INTEGER(4)
    !* ) sof   INTEGER(4)
    !Specific to particle mesh
    !1 ) ct    INTEGER(4)
    !2 ) npart INTEGER(8)
    !3 ) d1min REAL(num)
    !4 ) d1max REAL(num)
    !5 ) d2min REAL(num)
    !6 ) d2max REAL(num)
    !.
    !.
    !.
    !n ) dnmin REAL(num)
    !n+1) dnmax REAL(num)


    mdlen=MeshType_Header_Offset + 1 * SoI + 1 * SoI8  + ndim * 2 * num !1 INT, 1 INT8, 2REAL per Dim
    blocklen=mdlen + num*ndim*npart_global

    !Now written header, write metadata
    CALL cfd_Write_Block_Header(name,class,TYPE_MESH,blocklen,mdlen,default_rank)
    disp0=current_displacement
    CALL cfd_Write_MeshType_Header(MESH_PARTICLE,ndim,num,default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, Particle_Coord_Type, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * SoI

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, MPI_INTEGER8,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * SoI8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    DO i=1,ndim
       CALL MPI_ALLREDUCE(MINVAL(particles(:,i)), mn, 1, mpireal,MPI_MIN,cfd_comm,cfd_errcode)
       CALL MPI_ALLREDUCE(MAXVAL(particles(:,i)), mx, 1, mpireal,MPI_MAX,cfd_comm,cfd_errcode)
       IF (cfd_rank == default_rank) THEN
          CALL MPI_FILE_WRITE(cfd_filehandle, mn , 1, mpireal, cfd_status, cfd_errcode)
          CALL MPI_FILE_WRITE(cfd_filehandle, mx , 1, mpireal, cfd_status, cfd_errcode)
       ENDIF
       current_displacement = current_displacement + 2 * num
    ENDDO

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, Particle_Type,&
         "native", MPI_INFO_NULL, cfd_errcode)
    !Write the real data
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle,particles,npart_local*ndim,mpireal,cfd_status,cfd_errcode)
    current_displacement = current_displacement + ndim * npart_global * num

  END SUBROUTINE cfd_Write_nD_Particle_Grid_All

  !-------------------------------------------------------------------------------------------------------------------
  !Code to write a 2D Cartesian grid in serial from the node with rank {rank_write}
  !Serial operation, so no need to specify nx,ny
  !-------------------------------------------------------------------------------------------------------------------
  SUBROUTINE cfd_Write_nD_Particle_Grid_With_Iterator_All(name,class,iterator,ndims,npart_local,npart_global,npart_per_iteration,Particle_coord_type,Particle_Type)

    CHARACTER(len=*), INTENT(IN) :: name, class
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8), INTENT(IN) :: npart_local
    INTEGER(8), INTENT(IN) :: npart_per_iteration
    INTEGER(4), INTENT(IN) :: ndims
    INTEGER(4), INTENT(IN) :: Particle_coord_type
    INTEGER, INTENT(IN) :: Particle_Type
    REAL(num),ALLOCATABLE,DIMENSION(:) :: Data
    INTERFACE
       SUBROUTINE iterator(data,npart_it,direction,start)
         USE shared_data
         REAL(num),DIMENSION(:),INTENT(INOUT) :: data
         INTEGER,INTENT(IN) :: direction
         INTEGER(8),INTENT(INOUT) :: npart_it
         LOGICAL,INTENT(IN) :: start
       END SUBROUTINE iterator
    END INTERFACE
    INTEGER(8) :: blocklen,mdlen,npart_this_cycle,min_npart_this_cycle,npart_sent
    INTEGER(4) :: idim
    INTEGER(4) :: sizes(2)
    INTEGER(MPI_OFFSET_KIND) :: OffsetForMinMax
    REAL(num) :: mn,mx
    REAL(num),ALLOCATABLE,DIMENSION(:,:) :: MinMax
    LOGICAL :: start



    !Metadata is
    !* ) MeshType (INTEGER(4)) All mesh blocks contain this
    !* ) nd    INTEGER(4)
    !* ) sof   INTEGER(4)
    !Specific to particle mesh
    !1 ) ct    INTEGER(4)
    !2 ) npart INTEGER(8)
    !3 ) d1min REAL(num)
    !4 ) d1max REAL(num)
    !5 ) d2min REAL(num)
    !6 ) d2max REAL(num)
    !.
    !.
    !.
    !n ) dnmin REAL(num)
    !n+1) dnmax REAL(num)


    mdlen=MeshType_Header_Offset + 1 * SoI + 1 * SoI8  + ndims * 2 * num !1 INT, 1 INT8, 2REAL per Dim
    blocklen=mdlen + num*ndims*npart_global

    ALLOCATE(MinMax(1:ndims,1:2))
    MinMax=0.0_num

    !Now written header, write metadata
    CALL cfd_Write_Block_Header(name,class,TYPE_MESH,blocklen,mdlen,default_rank)
    CALL cfd_Write_MeshType_Header(MESH_PARTICLE,ndims,num,default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER, MPI_INTEGER,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, Particle_Coord_Type, 1, MPI_INTEGER, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * SoI

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, MPI_INTEGER8,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * SoI8

    !This is to skip past the location for the min/max values(Just write zeros). They will be filled in later
    OffsetForMinMax=current_displacement
    CALL MPI_FILE_SET_VIEW(cfd_filehandle, OffsetForMinMax, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, MinMax , ndims*2, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 2 * ndims * num

    !Write the real data

    start=.TRUE.
    ALLOCATE(Data(1:npart_per_iteration))
    npart_sent=0
    DO idim=1,ndims
       CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, Particle_Type,&
            "native", MPI_INFO_NULL, cfd_errcode)
       npart_this_cycle=npart_per_iteration
       start=.TRUE.
       DO
          CALL Iterator(Data,npart_this_cycle,idim,start)
          IF (npart_this_cycle <=0) EXIT
          IF (Start) THEN
             MinMax(idim,1)=MINVAL(Data(1:npart_this_cycle))
             MinMax(idim,2)=MAXVAL(Data(1:npart_this_cycle))
          ELSE
             MinMax(idim,1)=MIN(MinMax(idim,1),MINVAL(Data(1:npart_this_cycle)))
             MinMax(idim,2)=MAX(MinMax(idim,2),MAXVAL(Data(1:npart_this_cycle)))
          ENDIF
          start=.FALSE.
          npart_sent=npart_sent+npart_this_cycle
          CALL MPI_FILE_WRITE(cfd_filehandle,Data,npart_this_cycle,mpireal,cfd_status,cfd_errcode)
       ENDDO
       current_displacement = current_displacement +  npart_global * num
    ENDDO
    DEALLOCATE(Data)


    CALL MPI_FILE_SET_VIEW(cfd_filehandle, OffsetForMinMax, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    DO idim=1,ndims
       CALL MPI_ALLREDUCE(MinMax(idim,1), mn, 1, mpireal,MPI_MIN,cfd_comm,cfd_errcode)
       CALL MPI_ALLREDUCE(MinMax(idim,2), mx, 1, mpireal,MPI_MAX,cfd_comm,cfd_errcode)
       IF (cfd_rank == default_rank) THEN
          CALL MPI_FILE_WRITE(cfd_filehandle, mn , 1, mpireal, cfd_status, cfd_errcode)
          CALL MPI_FILE_WRITE(cfd_filehandle, mx , 1, mpireal, cfd_status, cfd_errcode)
       ENDIF
    ENDDO
    DEALLOCATE(MinMax)

    CALL MPI_BARRIER(comm,errcode)

  END SUBROUTINE cfd_Write_nD_Particle_Grid_With_Iterator_All

  !-------------------------------------------------------------------------------------------------------------------
  !Code to write a 2D Cartesian grid in serial from the node with rank {rank_write}
  !Serial operation, so no need to specify nx,ny
  !-------------------------------------------------------------------------------------------------------------------
  SUBROUTINE cfd_Write_nD_Particle_Variable_All(name,class,particles,npart_global,meshname,meshclass,Particle_Type)

    REAL(num),DIMENSION(:),INTENT(IN) :: particles
    CHARACTER(len=*), INTENT(IN) :: name, class, meshname, meshclass
    INTEGER,INTENT(IN) :: Particle_Type
    INTEGER(8), INTENT(IN) :: npart_global
    INTEGER(8) :: npart_local
    INTEGER(8) :: blocklen,mdlen
    INTEGER(4) :: i
    INTEGER(4) :: sizes(2)
    REAL(num) :: mn,mx


    npart_local = SIZE(particles)



    !Metadata is
    !* ) MeshType (INTEGER(4)) All mesh blocks contain this
    !* ) nd     INTEGER(4)
    !* ) sof    INTEGER(4)
    !Specific to particle variable
    !1 ) npart  INTEGER(8)
    !2 ) vmin   REAL(num)
    !3 ) vmax   REAL(num)
    !4 ) mesh   CHARACTER
    !5 ) mclass CHARACTER


    mdlen=MeshType_Header_Offset + 1 * SoI8  + 2 * num  + 2 * MaxStringLen
    blocklen=mdlen + num*npart_global

    !Now written header, write metadata
    CALL cfd_Write_Block_Header(name,class,TYPE_MESH_VARIABLE,blocklen,mdlen,default_rank)
    CALL cfd_Write_MeshType_Header(VAR_PARTICLE,DIMENSION_IRRELEVANT,num,default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, MPI_INTEGER8,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * SoI8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    CALL MPI_ALLREDUCE(MINVAL(particles), mn, 1, mpireal,MPI_MIN,cfd_comm,cfd_errcode)
    CALL MPI_ALLREDUCE(MAXVAL(particles), mx, 1, mpireal,MPI_MAX,cfd_comm,cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, mn , 1, mpireal, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, mx , 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 2 * num

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL cfd_Safe_Write_String(meshname)
       CALL cfd_Safe_Write_String(meshclass)
    ENDIF
    current_displacement = current_displacement + 2 * MaxStringLen

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, Particle_Type,&
         "native", MPI_INFO_NULL, cfd_errcode)
    !Write the real data
    CALL MPI_FILE_WRITE_ALL(cfd_filehandle,particles,npart_local,mpireal,cfd_status,cfd_errcode)
    current_displacement = current_displacement + npart_global * num

  END SUBROUTINE cfd_Write_nD_Particle_Variable_All

  SUBROUTINE cfd_Write_nD_Particle_Variable_With_Iterator_All(name,class,iterator,npart_global,npart_per_iteration,meshname,meshclass,Particle_Type)

    CHARACTER(len=*), INTENT(IN) :: name,class,meshname,meshclass
    INTEGER,INTENT(IN) :: Particle_Type
    INTEGER(8), INTENT(IN) :: npart_global,npart_per_iteration
    INTEGER(8) :: npart_this_cycle,min_npart_this_cycle

    REAL(num),ALLOCATABLE,DIMENSION(:) :: Data
    INTERFACE
       SUBROUTINE iterator(data,npart_it,start)
         USE shared_data
         REAL(num),DIMENSION(:),INTENT(INOUT) :: data
         INTEGER(8),INTENT(INOUT) :: npart_it
         LOGICAL,INTENT(IN) :: start
       END SUBROUTINE iterator
    END INTERFACE
    INTEGER(8) :: npart_local
    INTEGER(8) :: blocklen,mdlen
    INTEGER(4) :: i
    INTEGER(4) :: sizes(2)
    REAL(num) :: mn,mx,mn_g,mx_g
    INTEGER(MPI_OFFSET_KIND) :: OffsetForMinMax
    LOGICAL :: start


    !Metadata is
    !* ) MeshType (INTEGER(4)) All mesh blocks contain this
    !* ) nd     INTEGER(4)
    !* ) sof    INTEGER(4)
    !Specific to particle variable
    !1 ) npart  INTEGER(8)
    !2 ) vmin   REAL(num)
    !3 ) vmax   REAL(num)
    !4 ) mesh   CHARACTER
    !5 ) mclass CHARACTER


    mdlen=MeshType_Header_Offset + 1 * SoI8  + 2 * num  + 2 * MaxStringLen
    blocklen=mdlen + num*npart_global

    !Now written header, write metadata
    CALL cfd_Write_Block_Header(name,class,TYPE_MESH_VARIABLE,blocklen,mdlen,default_rank)
    CALL cfd_Write_MeshType_Header(VAR_PARTICLE,DIMENSION_IRRELEVANT,num,default_rank)

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_INTEGER8, MPI_INTEGER8,&
         "native", MPI_INFO_NULL,cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, npart_global, 1, MPI_INTEGER8, cfd_status, cfd_errcode)
    ENDIF
    current_displacement = current_displacement + 1 * SoI8

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, 0.0_num , 1, mpireal, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, 0.0_num , 1, mpireal, cfd_status, cfd_errcode)
    ENDIF
    OffsetForMinMax = current_displacement
    current_displacement = current_displacement + 2 * num

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, MPI_CHARACTER, MPI_CHARACTER,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL cfd_Safe_Write_String(meshname)
       CALL cfd_Safe_Write_String(meshclass)
    ENDIF
    current_displacement = current_displacement + 2 * MaxStringLen

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, current_displacement, mpireal, Particle_Type,&
         "native", MPI_INFO_NULL, cfd_errcode)

    start=.TRUE.
    npart_this_cycle=npart_per_iteration
    ALLOCATE(Data(1:npart_per_iteration))
    DO
       CALL Iterator(Data,npart_this_cycle,start)
       !PRINT *,rank,npart_this_cycle
       IF (npart_this_cycle <=0) EXIT
       IF (Start) THEN
          mn=MINVAL(Data(1:npart_this_cycle))
          mx=MAXVAL(Data(1:npart_this_cycle))
       ELSE
          mn=MIN(mn,MINVAL(Data(1:npart_this_cycle)))
          mx=MAX(mx,MAXVAL(Data(1:npart_this_cycle)))
       ENDIF
       start=.FALSE.
       !IF (cfd_rank .NE. 0) PRINT *,mn
       CALL MPI_FILE_WRITE(cfd_filehandle,Data,npart_this_cycle,mpireal,cfd_status,cfd_errcode)
    ENDDO
    DEALLOCATE(Data)
    current_displacement = current_displacement + npart_global * num

    CALL MPI_ALLREDUCE(mn, mn_g, 1, mpireal,MPI_MIN,cfd_comm,cfd_errcode)
    CALL MPI_ALLREDUCE(mx, mx_g, 1, mpireal,MPI_MAX,cfd_comm,cfd_errcode)
    mn=mn_g
    mx=mx_g

    CALL MPI_FILE_SET_VIEW(cfd_filehandle, OffsetForMinMax, mpireal, mpireal,&
         "native", MPI_INFO_NULL, cfd_errcode)
    IF (cfd_rank == default_rank) THEN
       CALL MPI_FILE_WRITE(cfd_filehandle, mn , 1, mpireal, cfd_status, cfd_errcode)
       CALL MPI_FILE_WRITE(cfd_filehandle, mx , 1, mpireal, cfd_status, cfd_errcode)
    ENDIF

    CALL MPI_BARRIER(comm,errcode)

  END SUBROUTINE cfd_Write_nD_Particle_Variable_With_Iterator_All
END MODULE output_particle

