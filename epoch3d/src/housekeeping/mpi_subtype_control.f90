! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE mpi_subtype_control

  !----------------------------------------------------------------------------
  ! This module contains the subroutines which create the subtypes used in
  ! IO
  !----------------------------------------------------------------------------

  USE shared_data

  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! get_total_local_particles - Returns the number of particles on this
  ! processor.
  !----------------------------------------------------------------------------

  FUNCTION get_total_local_particles()

    ! This subroutine describes the total number of particles on the current
    ! processor. It simply sums over every particle species

    INTEGER(i8) :: get_total_local_particles
    INTEGER :: ispecies

    get_total_local_particles = 0
    DO ispecies = 1, n_species
      get_total_local_particles = get_total_local_particles &
          + species_list(ispecies)%attached_list%count
    END DO

  END FUNCTION get_total_local_particles



  !----------------------------------------------------------------------------
  ! CreateSubtypes - Creates the subtypes used by the main output routines
  ! Run just before output takes place
  !----------------------------------------------------------------------------

  SUBROUTINE create_subtypes

    INTEGER :: i

    ! This subroutines creates the MPI types which represent the data for the
    ! field and particles data. It is used when writing data

    ! Actually create the subtypes
    subtype_field = create_current_field_subtype()
    subarray_field = create_current_field_subarray(ng)
    subarray_field_big = create_current_field_subarray(jng)

    subtype_field_r4 = create_current_field_subtype(MPI_REAL4)
    subarray_field_r4 = create_current_field_subarray(ng, MPI_REAL4)
    subarray_field_big_r4 = create_current_field_subarray(jng, MPI_REAL4)

    DO i = 1, n_subsets
      CALL create_subset_subtypes(i)
    END DO

  END SUBROUTINE create_subtypes



  !----------------------------------------------------------------------------
  ! Creates the subtypes used by the main output routines for subsets
  ! Run every time a time-varying subset is considered for output
  !----------------------------------------------------------------------------

  SUBROUTINE create_subset_subtypes(subset_index)

    INTEGER, INTENT(IN) :: subset_index
    INTEGER :: j, rd, n_min, n_max, mpitype, npd, npdm, n0
    INTEGER, DIMENSION(c_ndims) :: nt_global, n_local, n_global, starts
    INTEGER, DIMENSION(2,c_ndims) :: ranges
    LOGICAL :: proc_outside_range
    TYPE(subset), POINTER :: sub

    sub => subset_list(subset_index)

    nt_global = [nx_global, ny_global, nz_global]
    starts = 0
    n0 = 1

    ranges = cell_global_ranges(sub)
    n_global = ranges(2,:) - ranges(1,:)

    ranges = cell_local_ranges(sub)
    n_local = ranges(2,:) - ranges(1,:)

    starts = cell_starts(sub)

    ! These calculations rely on the original domain size, so will be wrong
    ! for skipped sets as yet
    proc_outside_range = .FALSE.
    DO j = 1, c_ndims
      IF (n_local(j) <= c_tiny) proc_outside_range = .TRUE.
    END DO

    IF (sub%skip) THEN
      DO j = 1, c_ndims
        rd = sub%skip_dir(j)
        n_min = n_global_min(j)
        n_max = n_global_max(j)
        npd = (n_max - n0) / rd + 1
        npdm = (n_min - 1 - n0) / rd + 1
        IF (n_min < 2) npdm = 0
        n_global(j) = (nt_global(j) - n0) / rd + 1
        n_local(j) = npd - npdm
        sub%n_start(j) = n0 + npdm * rd - n_min
        starts(j) = npdm
      END DO
    END IF

    ! Just exit if the subset hasn't changed extents since last time
    IF (equal(sub%n_global, n_global)) THEN
      IF (equal(sub%n_local, n_local)) THEN
        IF (equal(sub%starts, starts)) THEN
          RETURN
        END IF
      END IF
    END IF

    sub%n_global(:) = n_global(:)
    sub%n_local(:) = n_local(:)
    sub%starts(:) = starts(:)

    IF (sub%subtype /= MPI_DATATYPE_NULL) THEN
      CALL MPI_TYPE_FREE(sub%subtype, errcode)
      CALL MPI_TYPE_FREE(sub%subtype_r4, errcode)
      CALL MPI_TYPE_FREE(sub%subarray, errcode)
      CALL MPI_TYPE_FREE(sub%subarray_r4, errcode)
    END IF

    mpitype = MPI_DATATYPE_NULL
    IF (proc_outside_range) THEN
      CALL MPI_TYPE_CONTIGUOUS(0, mpireal, mpitype, errcode)
    ELSE
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, n_global, n_local, &
          starts, MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    END IF
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    sub%subtype = mpitype

    mpitype = MPI_DATATYPE_NULL
    IF (proc_outside_range) THEN
      CALL MPI_TYPE_CONTIGUOUS(0, MPI_REAL4, mpitype, errcode)
    ELSE
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, n_global, n_local, &
          starts, MPI_ORDER_FORTRAN, MPI_REAL4, mpitype, errcode)
    END IF
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    sub%subtype_r4 = mpitype

    IF (sub%skip) THEN
      ! Skipped varables are created in a temporary array that contains no
      ! ghost cells
      n_global = n_local
      starts = 0
    ELSE
      n_global = [nx + 2 * ng, ny + 2 * ng, nz + 2 * ng]
      starts = array_starts(sub)
    END IF

    mpitype = MPI_DATATYPE_NULL
    IF (proc_outside_range) THEN
      CALL MPI_TYPE_CONTIGUOUS(0, mpireal, mpitype, errcode)
    ELSE
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, n_global, n_local, &
          starts, MPI_ORDER_FORTRAN, mpireal, mpitype, errcode)
    END IF
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    sub%subarray = mpitype

    mpitype = MPI_DATATYPE_NULL
    IF (proc_outside_range) THEN
      CALL MPI_TYPE_CONTIGUOUS(0, MPI_REAL4, mpitype, errcode)
    ELSE
      CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, n_global, n_local, &
          starts, MPI_ORDER_FORTRAN, MPI_REAL4, mpitype, errcode)
    END IF
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    sub%subarray_r4 = mpitype

  END SUBROUTINE create_subset_subtypes



  !----------------------------------------------------------------------------
  ! Frees the subtypes created by create_subtypes
  !----------------------------------------------------------------------------

  SUBROUTINE free_subtypes

    CALL MPI_TYPE_FREE(subtype_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big, errcode)

    CALL MPI_TYPE_FREE(subtype_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big_r4, errcode)

  END SUBROUTINE free_subtypes



  !----------------------------------------------------------------------------
  ! create_current_field_subtype - Creates the subtype corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subtype(basetype_in)

    INTEGER :: create_current_field_subtype
    INTEGER, OPTIONAL, INTENT(IN) :: basetype_in
    INTEGER :: basetype

    IF (PRESENT(basetype_in)) THEN
      basetype = basetype_in
    ELSE
      basetype = mpireal
    END IF

    create_current_field_subtype = &
        create_field_subtype(basetype, nx, ny, nz, nx_global_min, &
            ny_global_min, nz_global_min)

  END FUNCTION create_current_field_subtype



  !----------------------------------------------------------------------------
  ! create_current_field_subarray - Creates the subarray corresponding to the
  ! current load balanced geometry
  !----------------------------------------------------------------------------

  FUNCTION create_current_field_subarray(ng, basetype_in)

    INTEGER :: create_current_field_subarray
    INTEGER, INTENT(IN) :: ng
    INTEGER, OPTIONAL, INTENT(IN) :: basetype_in
    INTEGER :: basetype

    IF (PRESENT(basetype_in)) THEN
      basetype = basetype_in
    ELSE
      basetype = mpireal
    END IF

    create_current_field_subarray = &
        create_field_subarray(basetype, ng, nx, ny, nz)

  END FUNCTION create_current_field_subarray



  !----------------------------------------------------------------------------
  ! create_subtypes_for_load - Creates subtypes when the code loads initial
  ! conditions from a file
  !----------------------------------------------------------------------------

  SUBROUTINE create_subtypes_for_load(species_subtypes, species_subtypes_i4, &
      species_subtypes_i8)

    ! This subroutines creates the MPI types which represent the data for the
    ! field and particles data. It is used when reading data.

    INTEGER, POINTER :: species_subtypes(:)
    INTEGER, POINTER :: species_subtypes_i4(:), species_subtypes_i8(:)
    INTEGER :: i

    subtype_field = create_current_field_subtype()
    subarray_field = create_current_field_subarray(ng)
    subarray_field_big = create_current_field_subarray(jng)

    subtype_field_r4 = create_current_field_subtype(MPI_REAL4)
    subarray_field_r4 = create_current_field_subarray(ng, MPI_REAL4)
    subarray_field_big_r4 = create_current_field_subarray(jng, MPI_REAL4)

    ALLOCATE(species_subtypes(n_species))
    ALLOCATE(species_subtypes_i4(n_species))
    ALLOCATE(species_subtypes_i8(n_species))
    DO i = 1,n_species
      CALL create_particle_subtypes(species_list(i)%attached_list%count, &
          species_subtypes(i), species_subtypes_i4(i), species_subtypes_i8(i))
    END DO

  END SUBROUTINE create_subtypes_for_load



  !----------------------------------------------------------------------------
  ! free_subtypes_for_load - Frees subtypes created by create_subtypes_for_load
  !----------------------------------------------------------------------------

  SUBROUTINE free_subtypes_for_load(species_subtypes, species_subtypes_i4, &
      species_subtypes_i8)

    INTEGER, POINTER :: species_subtypes(:)
    INTEGER, POINTER :: species_subtypes_i4(:), species_subtypes_i8(:)
    INTEGER :: i

    CALL MPI_TYPE_FREE(subtype_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big, errcode)

    CALL MPI_TYPE_FREE(subtype_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_r4, errcode)
    CALL MPI_TYPE_FREE(subarray_field_big_r4, errcode)
    DO i = 1,n_species
      CALL MPI_TYPE_FREE(species_subtypes(i), errcode)
      CALL MPI_TYPE_FREE(species_subtypes_i4(i), errcode)
      CALL MPI_TYPE_FREE(species_subtypes_i8(i), errcode)
    END DO
    DEALLOCATE(species_subtypes)
    DEALLOCATE(species_subtypes_i4)
    DEALLOCATE(species_subtypes_i8)

  END SUBROUTINE free_subtypes_for_load



  !----------------------------------------------------------------------------
  ! create_particle_subtype - Creates a subtype representing the local
  ! particles
  !----------------------------------------------------------------------------

  SUBROUTINE create_particle_subtypes(npart_in, subtype, subtype_i4, subtype_i8)

    INTEGER(i8), INTENT(IN) :: npart_in
    INTEGER, INTENT(OUT) :: subtype, subtype_i4, subtype_i8
    INTEGER(i8), DIMENSION(1) :: npart_local
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: npart_each_rank
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(3) :: disp
    INTEGER(KIND=MPI_ADDRESS_KIND) :: particles_to_skip, total_particles
    INTEGER :: i, mpitype, basetype, typesize

    npart_local = npart_in

    ALLOCATE(npart_each_rank(nproc))

    ! Create the subarray for the particles in this problem: subtype decribes
    ! where this process's data fits into the global picture.
    CALL MPI_ALLGATHER(npart_local, 1, MPI_INTEGER8, &
        npart_each_rank, 1, MPI_INTEGER8, comm, errcode)

    particles_to_skip = 0
    DO i = 1, rank
      particles_to_skip = particles_to_skip + npart_each_rank(i)
    END DO

    total_particles = particles_to_skip
    DO i = rank+1, nproc
      total_particles = total_particles + npart_each_rank(i)
    END DO

    DEALLOCATE(npart_each_rank)

    basetype = mpireal
    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)

    ! If npart_in is bigger than an integer then the data will not
    ! get written/read properly. This would require about 48GB per processor
    ! so it is unlikely to be a problem any time soon.
    lengths(1) = 1
    lengths(2) = INT(npart_in)
    lengths(3) = 1
    disp(1) = 0
    disp(2) = particles_to_skip * typesize
    disp(3) = total_particles * typesize
    types(1) = MPI_LB
    types(2) = basetype
    types(3) = MPI_UB

    mpitype = 0
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    subtype = mpitype

    basetype = MPI_INTEGER4
    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)

    ! If npart_in is bigger than an integer then the data will not
    ! get written/read properly. This would require about 48GB per processor
    ! so it is unlikely to be a problem any time soon.
    lengths(1) = 1
    lengths(2) = INT(npart_in)
    lengths(3) = 1
    disp(1) = 0
    disp(2) = particles_to_skip * typesize
    disp(3) = total_particles * typesize
    types(1) = MPI_LB
    types(2) = basetype
    types(3) = MPI_UB

    mpitype = 0
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    subtype_i4 = mpitype

    basetype = MPI_INTEGER8
    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)

    ! If npart_in is bigger than an integer then the data will not
    ! get written/read properly. This would require about 48GB per processor
    ! so it is unlikely to be a problem any time soon.
    lengths(1) = 1
    lengths(2) = INT(npart_in)
    lengths(3) = 1
    disp(1) = 0
    disp(2) = particles_to_skip * typesize
    disp(3) = total_particles * typesize
    types(1) = MPI_LB
    types(2) = basetype
    types(3) = MPI_UB

    mpitype = 0
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, mpitype, errcode)
    CALL MPI_TYPE_COMMIT(mpitype, errcode)
    subtype_i8 = mpitype

  END SUBROUTINE create_particle_subtypes



  !----------------------------------------------------------------------------
  ! create_field_subtype - Creates a subtype representing the local processor
  ! for any arbitrary arrangement of an array covering the entire spatial
  ! domain. Only used directly during load balancing
  !----------------------------------------------------------------------------

  FUNCTION create_field_subtype(basetype, nx_local, ny_local, nz_local, &
      cell_start_x_local, cell_start_y_local, cell_start_z_local)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, INTENT(IN) :: nx_local
    INTEGER, INTENT(IN) :: ny_local
    INTEGER, INTENT(IN) :: nz_local
    INTEGER, INTENT(IN) :: cell_start_x_local
    INTEGER, INTENT(IN) :: cell_start_y_local
    INTEGER, INTENT(IN) :: cell_start_z_local
    INTEGER :: create_field_subtype
    INTEGER, DIMENSION(c_ndims) :: n_local, n_global, start

    n_local = (/nx_local, ny_local, nz_local/)
    n_global = (/nx_global, ny_global, nz_global/)
    start = (/cell_start_x_local, cell_start_y_local, cell_start_z_local/)

    create_field_subtype = &
        create_3d_array_subtype(basetype, n_local, n_global, start)

  END FUNCTION create_field_subtype



  !----------------------------------------------------------------------------
  ! create_1d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 1D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_1d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec1d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(1), INTENT(IN) :: n_local
    INTEGER, DIMENSION(1), INTENT(IN) :: n_global
    INTEGER, DIMENSION(1), INTENT(IN) :: start
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(1)
    INTEGER :: vec1d, vec1d_sub, typesize

    vec1d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CONTIGUOUS(n_local(1), basetype, vec1d, errcode)
    CALL MPI_TYPE_COMMIT(vec1d, errcode)

    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    starts = start - 1
    lengths = 1

    disp(1) = 0
    disp(2) = typesize * starts(1)
    disp(3) = typesize * n_global(1)
    types(1) = MPI_LB
    types(2) = vec1d
    types(3) = MPI_UB

    vec1d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec1d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec1d_sub, errcode)

    CALL MPI_TYPE_FREE(vec1d, errcode)

  END FUNCTION create_1d_array_subtype



  !----------------------------------------------------------------------------
  ! create_2d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 2D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_2d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec2d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(2), INTENT(IN) :: n_local
    INTEGER, DIMENSION(2), INTENT(IN) :: n_global
    INTEGER, DIMENSION(2), INTENT(IN) :: start
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(2)
    INTEGER :: vec2d, vec2d_sub, typesize

    vec2d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_VECTOR(n_local(2), n_local(1), n_global(1), basetype, &
        vec2d, errcode)
    CALL MPI_TYPE_COMMIT(vec2d, errcode)

    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    starts = start - 1
    lengths = 1

    disp(1) = 0
    disp(2) = typesize * (starts(1) + n_global(1) * starts(2))
    disp(3) = typesize * n_global(1) * n_global(2)
    types(1) = MPI_LB
    types(2) = vec2d
    types(3) = MPI_UB

    vec2d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec2d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec2d_sub, errcode)

    CALL MPI_TYPE_FREE(vec2d, errcode)

  END FUNCTION create_2d_array_subtype



  !----------------------------------------------------------------------------
  ! create_3d_array_subtype - Creates a subtype representing the local fraction
  ! of a completely arbitrary 3D array. Does not assume anything about the
  ! domain at all.
  !----------------------------------------------------------------------------

  FUNCTION create_3d_array_subtype(basetype, n_local, n_global, start) &
      RESULT(vec3d_sub)

    INTEGER, INTENT(IN) :: basetype
    INTEGER, DIMENSION(3), INTENT(IN) :: n_local
    INTEGER, DIMENSION(3), INTENT(IN) :: n_global
    INTEGER, DIMENSION(3), INTENT(IN) :: start
    INTEGER, DIMENSION(3) :: lengths, types
    INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(3), starts(3)
    INTEGER :: vec2d, vec2d_sub
    INTEGER :: vec3d, vec3d_sub, typesize

    vec2d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_VECTOR(n_local(2), n_local(1), n_global(1), basetype, &
        vec2d, errcode)
    CALL MPI_TYPE_COMMIT(vec2d, errcode)

    CALL MPI_TYPE_SIZE(basetype, typesize, errcode)
    starts = start - 1
    lengths = 1

    disp(1) = 0
    disp(2) = typesize * (starts(1) + n_global(1) * starts(2))
    disp(3) = typesize * n_global(1) * n_global(2)
    types(1) = MPI_LB
    types(2) = vec2d
    types(3) = MPI_UB

    vec2d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec2d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec2d_sub, errcode)

    vec3d = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CONTIGUOUS(n_local(3), vec2d_sub, vec3d, errcode)
    CALL MPI_TYPE_COMMIT(vec3d, errcode)

    disp(1) = 0
    disp(2) = typesize * n_global(1) * n_global(2) * starts(3)
    disp(3) = typesize * n_global(1) * n_global(2) * n_global(3)
    types(1) = MPI_LB
    types(2) = vec3d
    types(3) = MPI_UB

    vec3d_sub = MPI_DATATYPE_NULL
    CALL MPI_TYPE_CREATE_STRUCT(3, lengths, disp, types, vec3d_sub, errcode)
    CALL MPI_TYPE_COMMIT(vec3d_sub, errcode)

    CALL MPI_TYPE_FREE(vec2d, errcode)
    CALL MPI_TYPE_FREE(vec2d_sub, errcode)
    CALL MPI_TYPE_FREE(vec3d, errcode)

  END FUNCTION create_3d_array_subtype



  FUNCTION create_field_subarray(basetype, ng, n1, n2, n3)

    INTEGER, INTENT(IN) :: basetype, ng, n1
    INTEGER, INTENT(IN), OPTIONAL :: n2, n3
    INTEGER, DIMENSION(3) :: n_local, n_global, start
    INTEGER :: i, ndim, create_field_subarray

    n_local(1) = n1
    ndim = 1
    IF (PRESENT(n2)) THEN
      n_local(2) = n2
      ndim = 2
    END IF
    IF (PRESENT(n3)) THEN
      n_local(3) = n3
      ndim = 3
    END IF

    DO i = 1, ndim
      start(i) = 1 + ng
      n_global(i) = n_local(i) + 2 * ng
    END DO

    IF (PRESENT(n3)) THEN
      create_field_subarray = &
          create_3d_array_subtype(basetype, n_local, n_global, start)
    ELSE IF (PRESENT(n2)) THEN
      create_field_subarray = &
          create_2d_array_subtype(basetype, n_local, n_global, start)
    ELSE
      create_field_subarray = &
          create_1d_array_subtype(basetype, n_local, n_global, start)
    END IF

  END FUNCTION create_field_subarray



  !----------------------------------------------------------------------------
  ! Clips range of current subset to domain size and cell edge
  !----------------------------------------------------------------------------

  SUBROUTINE get_global_ranges(current_subset, global_ranges)

    TYPE(subset), POINTER, INTENT(IN) :: current_subset
    REAL(num), DIMENSION(2,c_ndims), INTENT(OUT) :: global_ranges
    REAL(num) :: dd
    ! fudge factor allows overshoot of the specified domain extent by about 5%
    REAL(num), PARAMETER :: fudge = 0.019_num
    INTEGER :: idim, n1, n2

    DO idim = 1, c_ndims
      global_ranges(1,idim) = dir_grid_min(idim)
      global_ranges(2,idim) = dir_grid_max(idim)

      dd = dir_d(idim)

      IF (idim == 1) THEN
        n1 = c_subset_x_min
        n2 = c_subset_x_max
      ELSE IF (idim == 2) THEN
        n1 = c_subset_y_min
        n2 = c_subset_y_max
      ELSE
        n1 = c_subset_z_min
        n2 = c_subset_z_max
      END IF

      IF (current_subset%use_restriction(n1)) &
          global_ranges(1,idim) = current_subset%restriction(n1)
      IF (current_subset%use_restriction(n2)) &
          global_ranges(2,idim) = current_subset%restriction(n2)

      IF (global_ranges(2,idim) < global_ranges(1,idim)) THEN
        global_ranges = 0
        RETURN
      END IF

      ! Correct to match cell edges
      global_ranges(1,idim) = dir_min(idim) &
          + FLOOR((global_ranges(1,idim) - dir_min(idim)) / dd + fudge) * dd
      global_ranges(2,idim) = dir_min(idim) &
          + CEILING((global_ranges(2,idim) - dir_min(idim)) / dd - fudge) * dd

      ! Correct to domain size
      global_ranges(1,idim) = MAX(global_ranges(1,idim), dir_min(idim))
      global_ranges(2,idim) = MIN(global_ranges(2,idim), dir_max(idim))
    END DO

  END SUBROUTINE get_global_ranges



  FUNCTION cell_global_ranges(current_subset)

    INTEGER, DIMENSION(2,c_ndims) :: cell_global_ranges
    TYPE(subset), POINTER, INTENT(IN) :: current_subset
    REAL(num), DIMENSION(2,c_ndims) :: ranges
    REAL(num) :: dd, lower_posn
    INTEGER :: idim

    CALL get_global_ranges(current_subset, ranges)

    DO idim = 1, c_ndims
      dd = dir_d(idim)
      lower_posn = dir_min(idim)

      cell_global_ranges(1,idim) = NINT((ranges(1,idim) - lower_posn) / dd) + 1
      cell_global_ranges(2,idim) = NINT((ranges(2,idim) - lower_posn) / dd) + 1
    END DO

  END FUNCTION cell_global_ranges



  !----------------------------------------------------------------------------
  ! Location of current processors section of global array
  !----------------------------------------------------------------------------

  FUNCTION cell_local_ranges(current_subset)

    INTEGER, DIMENSION(2,c_ndims) :: cell_local_ranges
    TYPE(subset), POINTER, INTENT(IN) :: current_subset
    REAL(num), DIMENSION(2,c_ndims) :: ranges
    REAL(num) :: dd, lower_posn
    INTEGER :: idim

    CALL get_global_ranges(current_subset, ranges)

    DO idim = 1, c_ndims
      dd = dir_d(idim)
      lower_posn = dir_min(idim)

      ranges(1,idim) = MAX(ranges(1,idim), dir_min_local(idim))
      ranges(2,idim) = MIN(ranges(2,idim), dir_max_local(idim))

      cell_local_ranges(1,idim) = NINT((ranges(1,idim) - lower_posn) / dd) + 1
      cell_local_ranges(2,idim) = NINT((ranges(2,idim) - lower_posn) / dd) + 1
    END DO

  END FUNCTION cell_local_ranges



  !----------------------------------------------------------------------------
  ! Convert from local ranges into global array, to section of local array
  ! to use
  !----------------------------------------------------------------------------

  FUNCTION cell_section_ranges(current_subset)

    INTEGER, DIMENSION(2,c_ndims) :: cell_section_ranges
    TYPE(subset), POINTER, INTENT(IN) :: current_subset
    INTEGER, DIMENSION(2,c_ndims) :: ranges
    INTEGER :: idim, min_val

    ranges = cell_local_ranges(current_subset)

    DO idim = 1, c_ndims
      min_val = n_global_min(idim)

      cell_section_ranges(1,idim) = (ranges(1,idim) - min_val)
      cell_section_ranges(2,idim) = (ranges(2,idim) - min_val)
    END DO

  END FUNCTION cell_section_ranges



  FUNCTION cell_starts(current_subset)

    INTEGER, DIMENSION(c_ndims) :: cell_starts
    TYPE(subset), POINTER, INTENT(IN) :: current_subset
    INTEGER, DIMENSION(2,c_ndims) :: ranges
    REAL(num), DIMENSION(2,c_ndims) :: global_ranges
    INTEGER :: range_global_min, idim

    CALL get_global_ranges(current_subset, global_ranges)
    ranges = cell_local_ranges(current_subset)

    DO idim = 1, c_ndims
      range_global_min = &
          NINT((global_ranges(1,idim) - dir_min(idim)) / dir_d(idim))

      ! -1 because ranges is cell indexed and global_ranges isn't
      cell_starts(idim) = ranges(1,idim) - range_global_min - 1
    END DO

  END FUNCTION cell_starts



  FUNCTION array_starts(current_subset)

    INTEGER, DIMENSION(c_ndims) :: array_starts
    TYPE(subset), POINTER, INTENT(IN) :: current_subset
    INTEGER, DIMENSION(2,c_ndims) :: ranges
    INTEGER :: idim
    INTEGER, PARAMETER :: &
        nr(3) = [c_subset_x_min, c_subset_y_min, c_subset_z_min]

    array_starts = ng
    IF (.NOT.current_subset%space_restrictions) RETURN

    ranges = cell_global_ranges(current_subset)

    DO idim = 1, c_ndims
      IF (current_subset%use_restriction(nr(idim))) THEN
        array_starts(idim) = ranges(1,idim) - n_global_min(idim) + ng
        array_starts(idim) = MAX(array_starts(idim), ng)
      END IF
    END DO

  END FUNCTION array_starts



  ! Check if two 1-dimensional arrays are equal
  PURE LOGICAL FUNCTION equal(array1, array2)

    INTEGER, DIMENSION(:), INTENT(IN) :: array1, array2
    INTEGER :: i

    equal = SIZE(array1) == SIZE(array2)

    IF (equal) THEN
      DO i = 1, SIZE(array1)
        equal = array1(i) == array2(i)
        IF (.NOT. equal) RETURN
      END DO
    END IF

  END FUNCTION equal

END MODULE mpi_subtype_control
