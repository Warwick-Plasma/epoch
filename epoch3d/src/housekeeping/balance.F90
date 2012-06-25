MODULE balance

  USE mpi
  USE partlist
  USE boundary
  USE mpi_subtype_control
  USE redblack_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE balance_workload(over_ride)

    ! This subroutine determines whether or not the code needs rebalancing,
    ! calculates where to split the domain and calls other subroutines to
    ! actually rearrange the fields and particles onto the new processors

    ! This is really, really hard to do properly
    ! So cheat

    LOGICAL, INTENT(IN) :: over_ride
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: load_x
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: load_y
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: load_z
    INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_x_min, new_cell_x_max
    INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_y_min, new_cell_y_max
    INTEGER, DIMENSION(:), ALLOCATABLE :: new_cell_z_min, new_cell_z_max
    REAL(num) :: balance_frac, npart_av
    REAL(num) :: balance_frac_x, balance_frac_y, balance_frac_z
    INTEGER(KIND=8) :: min_x, max_x, min_y, max_y, min_z, max_z
    INTEGER(KIND=8) :: npart_local, sum_npart, max_npart, wk
    INTEGER :: iproc
    INTEGER, DIMENSION(c_ndims,2) :: domain
#ifdef PARTICLE_DEBUG
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
#endif

    ! On one processor do nothing to save time
    IF (nproc .EQ. 1) RETURN

    ! This parameter allows selecting the mode of the autobalancing between
    ! leftsweep, rightsweep, auto(best of leftsweep and rightsweep) or both
    balance_mode = c_lb_all

    ! count particles
    npart_local = get_total_local_particles()

    ! The over_ride flag allows the code to force a load balancing sweep
    ! at t = 0
    IF (.NOT. over_ride) THEN
      CALL MPI_ALLREDUCE(npart_local, max_npart, 1, MPI_INTEGER8, MPI_MAX, &
          comm, errcode)
      IF (max_npart .LE. 0) RETURN
      CALL MPI_ALLREDUCE(npart_local, sum_npart, 1, MPI_INTEGER8, MPI_SUM, &
          comm, errcode)
      npart_av = REAL(sum_npart, num) / nproc
      balance_frac = (npart_av + SQRT(npart_av)) / REAL(max_npart, num)
      IF (balance_frac .GT. dlb_threshold) RETURN
      IF (rank .EQ. 0) PRINT *, 'Load balancing with fraction', balance_frac
    ENDIF

    ALLOCATE(new_cell_x_min(nprocx), new_cell_x_max(nprocx))
    ALLOCATE(new_cell_y_min(nprocy), new_cell_y_max(nprocy))
    ALLOCATE(new_cell_z_min(nprocz), new_cell_z_max(nprocz))

    ! Sweep in X
    IF (IAND(balance_mode, c_lb_x) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in X
      ALLOCATE(load_x(nx_global))
      CALL get_load_in_x(load_x)
      CALL calculate_breaks(load_x, nprocx, new_cell_x_min, new_cell_x_max)
    ELSE
      ! Just keep the original lengths
      new_cell_x_min = cell_x_min
      new_cell_x_max = cell_x_max
    ENDIF

    ! Sweep in Y
    IF (IAND(balance_mode, c_lb_y) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in Y
      ALLOCATE(load_y(ny_global))
      CALL get_load_in_y(load_y)
      CALL calculate_breaks(load_y, nprocy, new_cell_y_min, new_cell_y_max)
    ELSE
      ! Just keep the original lengths
      new_cell_y_min = cell_y_min
      new_cell_y_max = cell_y_max
    ENDIF

    ! Sweep in Z
    IF (IAND(balance_mode, c_lb_z) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in Z
      ALLOCATE(load_z(nz_global))
      CALL get_load_in_z(load_z)
      CALL calculate_breaks(load_z, nprocz, new_cell_z_min, new_cell_z_max)
    ELSE
      ! Just keep the original lengths
      new_cell_z_min = cell_z_min
      new_cell_z_max = cell_z_max
    ENDIF

    ! In the autobalancer then determine whether to balance in X, Y or Z
    ! Is this worth keeping?
    IF (IAND(balance_mode, c_lb_auto) .NE. 0 ) THEN

      ! Code is auto load balancing
      max_x = 0
      min_x = npart_global
      DO iproc = 1, nprocx
        wk = SUM(load_x(new_cell_x_min(iproc):new_cell_x_max(iproc)))
        IF (wk .GT. max_x) max_x = wk
        IF (wk .LT. min_x) min_x = wk
      ENDDO

      max_y = 0
      min_y = npart_global
      DO iproc = 1, nprocy
        wk = SUM(load_y(new_cell_y_min(iproc):new_cell_y_max(iproc)))
        IF (wk .GT. max_y) max_y = wk
        IF (wk .LT. min_y) min_y = wk
      ENDDO

      max_z = 0
      min_z = npart_global
      DO iproc = 1, nprocz
        wk = SUM(load_z(new_cell_z_min(iproc):new_cell_z_max(iproc)))
        IF (wk .GT. max_z) max_z = wk
        IF (wk .LT. min_z) min_z = wk
      ENDDO

      balance_frac_x = REAL(min_x, num) / REAL(max_x, num)
      balance_frac_y = REAL(min_y, num) / REAL(max_y, num)
      balance_frac_z = REAL(min_z, num) / REAL(max_z, num)

      IF (balance_frac_x .LT. balance_frac_y) THEN
        IF (balance_frac_x .LT. balance_frac_z) THEN
          new_cell_x_min = cell_x_min
          new_cell_x_max = cell_x_max
        ELSE
          new_cell_z_min = cell_z_min
          new_cell_z_max = cell_z_max
        ENDIF
      ELSE
        IF (balance_frac_y .LT. balance_frac_z) THEN
          new_cell_y_min = cell_y_min
          new_cell_y_max = cell_y_max
        ELSE
          new_cell_z_min = cell_z_min
          new_cell_z_max = cell_z_max
        ENDIF
      ENDIF

    ENDIF

    IF (ALLOCATED(load_x)) DEALLOCATE(load_x)
    IF (ALLOCATED(load_y)) DEALLOCATE(load_y)
    IF (ALLOCATED(load_z)) DEALLOCATE(load_z)

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor

    domain(1,:) = (/new_cell_x_min(x_coords+1), new_cell_x_max(x_coords+1)/)
    domain(2,:) = (/new_cell_y_min(y_coords+1), new_cell_y_max(y_coords+1)/)
    domain(3,:) = (/new_cell_z_min(z_coords+1), new_cell_z_max(z_coords+1)/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    ! Copy the new lengths into the permanent variables
    cell_x_min = new_cell_x_min
    cell_x_max = new_cell_x_max
    cell_y_min = new_cell_y_min
    cell_y_max = new_cell_y_max
    cell_z_min = new_cell_z_min
    cell_z_max = new_cell_z_max

    ! Set the new nx, ny, nz
    nx_global_min = cell_x_min(x_coords+1)
    nx_global_max = cell_x_max(x_coords+1)

    ny_global_min = cell_y_min(y_coords+1)
    ny_global_max = cell_y_max(y_coords+1)

    nz_global_min = cell_z_min(z_coords+1)
    nz_global_max = cell_z_max(z_coords+1)

    nx = nx_global_max - nx_global_min + 1
    ny = ny_global_max - ny_global_min + 1
    nz = nz_global_max - nz_global_min + 1

    ! Do X, Y, Z arrays separately because we already have global copies
    DEALLOCATE(x, y, z)
    ALLOCATE(x(-2:nx+3), y(-2:ny+3), z(-2:nz+3))
    x(-2:nx+3) = x_global(nx_global_min-3:nx_global_max+3)
    y(-2:ny+3) = y_global(ny_global_min-3:ny_global_max+3)
    z(-2:nz+3) = z_global(nz_global_min-3:nz_global_max+3)

    ! Recalculate x_mins and x_maxs so that rebalancing works next time
    DO iproc = 1, nprocx
      x_mins(iproc) = x_global(cell_x_min(iproc))
      x_maxs(iproc) = x_global(cell_x_max(iproc))
    ENDDO
    ! Same for y
    DO iproc = 1, nprocy
      y_mins(iproc) = y_global(cell_y_min(iproc))
      y_maxs(iproc) = y_global(cell_y_max(iproc))
    ENDDO
    ! Same for z
    DO iproc = 1, nprocz
      z_mins(iproc) = z_global(cell_z_min(iproc))
      z_maxs(iproc) = z_global(cell_z_max(iproc))
    ENDDO

    ! Set the lengths of the current domain so that the particle balancer
    ! works properly
    x_min_local = x_mins(x_coords)
    x_max_local = x_maxs(x_coords)
    y_min_local = y_mins(y_coords)
    y_max_local = y_maxs(y_coords)
    z_min_local = z_mins(z_coords)
    z_max_local = z_maxs(z_coords)

    ! Redistribute the particles onto their new processors
    CALL distribute_particles

    ! If running with particle debugging then set the t = 0 processor if
    ! over_ride = true
#ifdef PARTICLE_DEBUG
    IF (over_ride) THEN
      DO ispecies = 1, n_species
        current=>species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          current%processor_at_t0 = rank
          current=>current%next
        ENDDO
      ENDDO
    ENDIF
#endif

  END SUBROUTINE balance_workload



  SUBROUTINE redistribute_fields(new_domain)

    ! This subroutine redistributes the field variables over the new
    ! processor layout. If using a field of your own then set the
    ! redistribute_field subroutine to implement it.

    INTEGER :: nx_new, ny_new, nz_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: temp_sum
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp_slice
    TYPE(laser_block), POINTER :: current
    INTEGER :: ispecies, index, n_species_local

    nx_new = new_domain(1,2) - new_domain(1,1) + 1
    ny_new = new_domain(2,2) - new_domain(2,1) + 1
    nz_new = new_domain(3,2) - new_domain(3,1) + 1

    ALLOCATE(temp(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))

    CALL redistribute_field(new_domain, ex, temp)
    CALL redistribute_field(new_domain, ey, temp)
    CALL redistribute_field(new_domain, ez, temp)

    CALL redistribute_field(new_domain, bx, temp)
    CALL redistribute_field(new_domain, by, temp)
    CALL redistribute_field(new_domain, bz, temp)

    CALL redistribute_field(new_domain, jx, temp)
    CALL redistribute_field(new_domain, jy, temp)
    CALL redistribute_field(new_domain, jz, temp)

    IF (cpml_boundaries) THEN
      CALL redistribute_field(new_domain, cpml_psi_eyx, temp)
      CALL redistribute_field(new_domain, cpml_psi_byx, temp)
      CALL redistribute_field(new_domain, cpml_psi_ezx, temp)
      CALL redistribute_field(new_domain, cpml_psi_bzx, temp)

      CALL redistribute_field(new_domain, cpml_psi_exy, temp)
      CALL redistribute_field(new_domain, cpml_psi_bxy, temp)
      CALL redistribute_field(new_domain, cpml_psi_ezy, temp)
      CALL redistribute_field(new_domain, cpml_psi_bzy, temp)

      CALL redistribute_field(new_domain, cpml_psi_exz, temp)
      CALL redistribute_field(new_domain, cpml_psi_bxz, temp)
      CALL redistribute_field(new_domain, cpml_psi_eyz, temp)
      CALL redistribute_field(new_domain, cpml_psi_byz, temp)

      CALL deallocate_cpml_helpers
      CALL set_cpml_helpers(nx_new, new_domain(1,1), new_domain(1,2), &
          ny_new, new_domain(2,1), new_domain(2,2), &
          nz_new, new_domain(3,1), new_domain(3,2))
    ENDIF

    DEALLOCATE(temp)

    DO index = 1, num_vars_to_dump
      IF (.NOT. ASSOCIATED(averaged_data(index)%array)) CYCLE

      n_species_local = 0
      IF (IAND(dumpmask(index), c_io_no_sum) .EQ. 0) &
          n_species_local = 1
      IF (IAND(dumpmask(index), c_io_species) .NE. 0) &
          n_species_local = n_species_local + n_species

      IF (n_species_local .LE. 0) CYCLE

      ALLOCATE(temp_sum(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3, n_species_local))

      CALL redistribute_field_sum(new_domain, averaged_data(index)%array, &
          temp_sum)

      DEALLOCATE(temp_sum)
    ENDDO

    ALLOCATE(temp_slice(-2:ny_new+3, -2:nz_new+3))

    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      CALL redistribute_field_slice_ptr(new_domain, c_dir_x, current%profile, &
          temp_slice)
      CALL redistribute_field_slice_ptr(new_domain, c_dir_x, current%phase, &
          temp_slice)

      current=>current%next
    ENDDO

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      CALL redistribute_field_slice_ptr(new_domain, c_dir_x, current%profile, &
          temp_slice)
      CALL redistribute_field_slice_ptr(new_domain, c_dir_x, current%phase, &
          temp_slice)

      current=>current%next
    ENDDO

    CALL redistribute_field_slice(new_domain, c_dir_x, ex_x_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, ex_x_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, ey_x_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, ey_x_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, ez_x_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, ez_x_max, temp_slice)

    CALL redistribute_field_slice(new_domain, c_dir_x, bx_x_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, bx_x_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, by_x_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, by_x_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, bz_x_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_x, bz_x_max, temp_slice)

    DEALLOCATE(temp_slice)
    ALLOCATE(temp_slice(-2:nx_new+3, -2:nz_new+3))

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      CALL redistribute_field_slice_ptr(new_domain, c_dir_y, current%profile, &
          temp_slice)
      CALL redistribute_field_slice_ptr(new_domain, c_dir_y, current%phase, &
          temp_slice)

      current=>current%next
    ENDDO

    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      CALL redistribute_field_slice_ptr(new_domain, c_dir_y, current%profile, &
          temp_slice)
      CALL redistribute_field_slice_ptr(new_domain, c_dir_y, current%phase, &
          temp_slice)

      current=>current%next
    ENDDO

    CALL redistribute_field_slice(new_domain, c_dir_y, ex_y_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, ex_y_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, ey_y_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, ey_y_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, ez_y_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, ez_y_max, temp_slice)

    CALL redistribute_field_slice(new_domain, c_dir_y, bx_y_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, bx_y_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, by_y_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, by_y_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, bz_y_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_y, bz_y_max, temp_slice)

    DEALLOCATE(temp_slice)
    ALLOCATE(temp_slice(-2:nx_new+3, -2:ny_new+3))

    current=>laser_z_min
    DO WHILE(ASSOCIATED(current))
      CALL redistribute_field_slice_ptr(new_domain, c_dir_z, current%profile, &
          temp_slice)
      CALL redistribute_field_slice_ptr(new_domain, c_dir_z, current%phase, &
          temp_slice)

      current=>current%next
    ENDDO

    current=>laser_z_max
    DO WHILE(ASSOCIATED(current))
      CALL redistribute_field_slice_ptr(new_domain, c_dir_z, current%profile, &
          temp_slice)
      CALL redistribute_field_slice_ptr(new_domain, c_dir_z, current%phase, &
          temp_slice)

      current=>current%next
    ENDDO

    CALL redistribute_field_slice(new_domain, c_dir_z, ex_z_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, ex_z_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, ey_z_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, ey_z_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, ez_z_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, ez_z_max, temp_slice)

    CALL redistribute_field_slice(new_domain, c_dir_z, bx_z_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, bx_z_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, by_z_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, by_z_max, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, bz_z_min, temp_slice)
    CALL redistribute_field_slice(new_domain, c_dir_z, bz_z_max, temp_slice)

    DEALLOCATE(temp_slice)

    ! Re-distribute moving windows and temperature boundaries
    IF (move_window) THEN
      ALLOCATE(temp_slice(-2:ny_new+3, -2:nz_new+3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_ptr(new_domain, c_dir_x, &
            species_list(ispecies)%density, temp_slice)
      ENDDO

      DEALLOCATE(temp_slice)

      ALLOCATE(temp(-2:ny_new+3, -2:nz_new+3, 3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_sum(new_domain, c_dir_x, &
            species_list(ispecies)%temperature, temp)
      ENDDO

      DEALLOCATE(temp)
    ENDIF

    IF (bc_particle(c_bd_x_min) .EQ. c_bc_thermal) THEN
      ALLOCATE(temp(-2:ny_new+3, -2:nz_new+3, 3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_sum(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_min, temp)
      ENDDO

      DEALLOCATE(temp)
    ENDIF

    IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
      ALLOCATE(temp(-2:ny_new+3, -2:nz_new+3, 3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_sum(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_max, temp)
      ENDDO

      DEALLOCATE(temp)
    ENDIF

    IF (bc_particle(c_bd_y_min) .EQ. c_bc_thermal) THEN
      ALLOCATE(temp(-2:nx_new+3, -2:nz_new+3, 3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_sum(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_min, temp)
      ENDDO

      DEALLOCATE(temp)
    ENDIF

    IF (bc_particle(c_bd_y_max) .EQ. c_bc_thermal) THEN
      ALLOCATE(temp(-2:nx_new+3, -2:nz_new+3, 3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_sum(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_max, temp)
      ENDDO

      DEALLOCATE(temp)
    ENDIF

    IF (bc_particle(c_bd_z_min) .EQ. c_bc_thermal) THEN
      ALLOCATE(temp(-2:nx_new+3, -2:ny_new+3, 3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_sum(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_min, temp)
      ENDDO

      DEALLOCATE(temp)
    ENDIF

    IF (bc_particle(c_bd_z_max) .EQ. c_bc_thermal) THEN
      ALLOCATE(temp(-2:nx_new+3, -2:ny_new+3, 3))

      DO ispecies = 1, n_species
        CALL redistribute_field_slice_sum(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_max, temp)
      ENDDO

      DEALLOCATE(temp)
    ENDIF

  END SUBROUTINE redistribute_fields



  SUBROUTINE redistribute_field(new_domain, field, temp)

    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: field
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: temp
    INTEGER :: n_new(c_ndims)

    n_new = SHAPE(temp)

    temp = 0.0_num
    CALL redistribute_field_3d(new_domain, field, temp)
    DEALLOCATE(field)
    ALLOCATE(field(-2:n_new(1)+3, -2:n_new(2)+3, -2:n_new(3)+3))
    field = temp

  END SUBROUTINE redistribute_field



  SUBROUTINE redistribute_field_slice(new_domain, direction, field, temp)

    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: field
    REAL(num), DIMENSION(:,:), INTENT(OUT) :: temp
    INTEGER :: n_new(c_ndims-1)

    n_new = SHAPE(temp)

    temp = 0.0_num
    CALL redistribute_field_2d(new_domain, direction, field, temp)
    DEALLOCATE(field)
    ALLOCATE(field(-2:n_new(1)+3, -2:n_new(2)+3))
    field = temp

  END SUBROUTINE redistribute_field_slice



  SUBROUTINE redistribute_field_slice_ptr(new_domain, direction, field, temp)

    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(:,:), POINTER, INTENT(INOUT) :: field
    REAL(num), DIMENSION(:,:), INTENT(OUT) :: temp
    INTEGER :: n_new(c_ndims-1)

    n_new = SHAPE(temp)

    temp = 0.0_num
    CALL redistribute_field_2d(new_domain, direction, field, temp)
    DEALLOCATE(field)
    ALLOCATE(field(-2:n_new(1)+3, -2:n_new(2)+3))
    field = temp

  END SUBROUTINE redistribute_field_slice_ptr



  SUBROUTINE redistribute_field_sum(new_domain, field, temp)

    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:,:,:), POINTER, INTENT(INOUT) :: field
    REAL(num), DIMENSION(:,:,:,:), INTENT(OUT) :: temp
    INTEGER :: i, n_new(c_ndims+1)

    n_new = SHAPE(temp)

    temp = 0.0_num
    DO i = 1, n_new(c_ndims+1)
      CALL redistribute_field_3d(new_domain, field(:,:,:,i), temp(:,:,:,i))
    ENDDO

    DEALLOCATE(field)
    ALLOCATE(field(-2:n_new(1)+3,-2:n_new(2)+3,-2:n_new(3)+3,n_new(c_ndims+1)))

    field = temp

  END SUBROUTINE redistribute_field_sum



  SUBROUTINE redistribute_field_slice_sum(new_domain, direction, field, temp)

    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: field
    REAL(num), DIMENSION(:,:,:), INTENT(OUT) :: temp
    INTEGER :: i, n_new(c_ndims)

    n_new = SHAPE(temp)

    temp = 0.0_num
    DO i = 1, n_new(c_ndims)
      CALL redistribute_field_2d(new_domain, direction, field(:,:,i), &
          temp(:,:,i))
    ENDDO

    DEALLOCATE(field)
    ALLOCATE(field(-2:n_new(1)+3, -2:n_new(2)+3, n_new(c_ndims+1)))

    field = temp

  END SUBROUTINE redistribute_field_slice_sum



  SUBROUTINE redistribute_field_2d(domain, direction, field_in, field_out)

    ! This subroutine redistributes the fields over the new processor layout
    ! The current version works by writing the field to a file and then each
    ! processor loads back in it's own part. This is better than the previous
    ! version where each processor produced it's own copy of the global array
    ! and then took its own subsection
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: domain
    INTEGER, INTENT(IN) :: direction
    REAL(num), DIMENSION(-2:,-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:,-2:), INTENT(OUT) :: field_out
    INTEGER, DIMENSION(2) :: n, n_new, n_global, n_old_start, n_start
    INTEGER :: n1_dir, n2_dir
    INTEGER :: subarray_write, subarray_read
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename

    IF (direction .EQ. c_dir_x) THEN
      n1_dir = 2
      n2_dir = 3
      n(1) = ny
      n(2) = nz
      n_global(1) = ny_global
      n_global(2) = nz_global
      n_old_start(1) = cell_y_min(y_coords+1)
      n_old_start(2) = cell_z_min(z_coords+1)
    ELSE IF (direction .EQ. c_dir_y) THEN
      n1_dir = 1
      n2_dir = 3
      n(1) = nx
      n(2) = nz
      n_global(1) = nx_global
      n_global(2) = nz_global
      n_old_start(1) = cell_x_min(x_coords+1)
      n_old_start(2) = cell_z_min(z_coords+1)
    ELSE
      n1_dir = 1
      n2_dir = 2
      n(1) = nx
      n(2) = ny
      n_global(1) = nx_global
      n_global(2) = ny_global
      n_old_start(1) = cell_x_min(x_coords+1)
      n_old_start(2) = cell_y_min(y_coords+1)
    ENDIF

    WRITE(filename, '(a, ''/balance.dat'')') TRIM(data_dir)

    n_new(1) = domain(n1_dir,2) - domain(n1_dir,1) + 1
    n_new(2) = domain(n2_dir,2) - domain(n2_dir,1) + 1
    n_start(1) = domain(n1_dir,1)
    n_start(2) = domain(n2_dir,1)

    subarray_write = create_field_subarray(n(1), n(2))
    subtype_write = create_2d_array_subtype(n, n_global, n_old_start)

    subarray_read = create_field_subarray(n_new(1), n_new(2))
    subtype_read = create_2d_array_subtype(n_new, n_global, n_start)

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDWR+MPI_MODE_CREATE, &
        MPI_INFO_NULL, fh, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_write, 'native', &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(fh, field_in, 1, subarray_write, status, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_read, 'native', &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_READ_ALL(fh, field_out, 1, subarray_read, status, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

    CALL MPI_TYPE_FREE(subarray_write, errcode)
    CALL MPI_TYPE_FREE(subtype_write, errcode)
    CALL MPI_TYPE_FREE(subarray_read, errcode)
    CALL MPI_TYPE_FREE(subtype_read, errcode)

  END SUBROUTINE redistribute_field_2d



  SUBROUTINE redistribute_field_3d(domain, field_in, field_out)

    ! This subroutine redistributes the fields over the new processor layout
    ! The current version works by writing the field to a file and then each
    ! processor loads back in it's own part. This is better than the previous
    ! version where each processor produced it's own copy of the global array
    ! and then took its own subsection
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: domain
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: field_in
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(OUT) :: field_out
    INTEGER :: nx_new, ny_new, nz_new
    INTEGER :: subarray_write, subarray_read
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
    CHARACTER(LEN=9+data_dir_max_length+n_zeros) :: filename

    WRITE(filename, '(a, ''/balance.dat'')') TRIM(data_dir)

    nx_new = domain(1,2) - domain(1,1) + 1
    ny_new = domain(2,2) - domain(2,1) + 1
    nz_new = domain(3,2) - domain(3,1) + 1

    subarray_write = create_current_field_subarray()
    subtype_write  = create_current_field_subtype()

    subarray_read = create_field_subarray(nx_new, ny_new, nz_new)
    subtype_read  = create_field_subtype(nx_new, ny_new, nz_new, &
        domain(1,1), domain(2,1), domain(3,1))

    CALL MPI_FILE_OPEN(comm, TRIM(filename), MPI_MODE_RDWR+MPI_MODE_CREATE, &
        MPI_INFO_NULL, fh, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_write, 'native', &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(fh, field_in, 1, subarray_write, status, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, mpireal, subtype_read, 'native', &
        MPI_INFO_NULL, errcode)
    CALL MPI_FILE_READ_ALL(fh, field_out, 1, subarray_read, status, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

    CALL MPI_TYPE_FREE(subarray_write, errcode)
    CALL MPI_TYPE_FREE(subtype_write, errcode)
    CALL MPI_TYPE_FREE(subarray_read, errcode)
    CALL MPI_TYPE_FREE(subtype_read, errcode)

    CALL do_field_mpi_with_lengths(field_out, nx_new, ny_new, nz_new)

  END SUBROUTINE redistribute_field_3d



  SUBROUTINE get_load_in_x(load)

    ! Calculate total load across the X direction
    ! Summed in the Y,Z directions

    INTEGER(KIND=8), DIMENSION(:), INTENT(OUT) :: load
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, sz

    load = 0

    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so x_min, NOT x_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(1) - x_min) / dx) + 1
#else
        cell = FLOOR((current%part_pos(1) - x_min) / dx + 1.5_num)
#endif
        load(cell) = load(cell) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    sz = SIZE(load)
    ALLOCATE(temp(sz))
    CALL MPI_ALLREDUCE(load, temp, sz, MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * temp + ny_global * nz_global

    DEALLOCATE(temp)

  END SUBROUTINE get_load_in_x



  SUBROUTINE get_load_in_y(load)

    ! Calculate total load across the Y direction
    ! Summed in the X,Z directions

    INTEGER(KIND=8), DIMENSION(:), INTENT(OUT) :: load
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, sz

    load = 0

    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so y_min, NOT y_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(2) - y_min) / dy) + 1
#else
        cell = FLOOR((current%part_pos(2) - y_min) / dy + 1.5_num)
#endif
        load(cell) = load(cell) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    sz = SIZE(load)
    ALLOCATE(temp(sz))
    CALL MPI_ALLREDUCE(load, temp, sz, MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * temp + nx_global * nz_global

    DEALLOCATE(temp)

  END SUBROUTINE get_load_in_y



  SUBROUTINE get_load_in_z(load)

    ! Calculate total load across the Z direction
    ! Summed in the Y,Z directions

    INTEGER(KIND=8), DIMENSION(:), INTENT(OUT) :: load
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: temp
    TYPE(particle), POINTER :: current
    INTEGER :: cell, ispecies, sz

    load = 0

    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        ! Want global position, so z_min, NOT z_min_local
#ifdef PARTICLE_SHAPE_TOPHAT
        cell = FLOOR((current%part_pos(3) - z_min) / dz) + 1
#else
        cell = FLOOR((current%part_pos(3) - z_min) / dz + 1.5_num)
#endif
        load(cell) = load(cell) + 1
        current=>current%next
      ENDDO
    ENDDO

    ! Now have local densities, so add using MPI
    sz = SIZE(load)
    ALLOCATE(temp(sz))
    CALL MPI_ALLREDUCE(load, temp, sz, MPI_INTEGER8, MPI_SUM, comm, errcode)

    ! Adjust the load of pushing one particle relative to the load
    ! of updating one field cell, then add on the field load.
    ! The push_per_field factor will be updated automatically in future.
    load = push_per_field * temp + nx_global * ny_global

    DEALLOCATE(temp)

  END SUBROUTINE get_load_in_z



  SUBROUTINE calculate_breaks(load, nproc, mins, maxs)

    ! This subroutine calculates the places in a given load profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(KIND=8), INTENT(IN), DIMENSION(:) :: load
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: mins, maxs
    INTEGER :: sz, idim, proc, old
    INTEGER(KIND=8) :: total, total_old, load_per_proc_ideal

    sz = SIZE(load)
    maxs = sz

    load_per_proc_ideal = FLOOR((SUM(load) + 0.5d0) / nproc, 8)

    proc = 1
    total = 0
    DO idim = 1, sz
      total_old = total
      total = total + load(idim)
      IF (total .GE. load_per_proc_ideal) THEN
        ! Pick the split that most closely matches the load
        IF (load_per_proc_ideal - total_old &
            .LT. total - load_per_proc_ideal) THEN
          maxs(proc) = idim - 1
        ELSE
          maxs(proc) = idim
        ENDIF
        proc = proc + 1
        total = total - load_per_proc_ideal
        IF (proc .EQ. nproc) EXIT
      ENDIF
    ENDDO

    ! Sanity check. Must be one cell of separation between each endpoint.
    ! Forwards (unnecessary?)
    old = 0
    DO proc = 1, nproc
      IF (maxs(proc) - old .LE. 0) THEN
        maxs(proc) = old + 1
      ENDIF
      old = maxs(proc)
    ENDDO

    ! Backwards
    old = sz + 1
    DO proc = nproc, 1, -1
      IF (old - maxs(proc) .LE. 0) THEN
        maxs(proc) = old - 1
      ENDIF
      old = maxs(proc)
    ENDDO

    ! Set mins
    mins(1) = 1
    DO proc = 2, nproc
      mins(proc) = maxs(proc-1) + 1
    ENDDO

  END SUBROUTINE calculate_breaks



  FUNCTION get_particle_processor(a_particle)

    ! This subroutine calculates which processor a given particles resides on

    TYPE(particle), INTENT(IN) :: a_particle
    INTEGER :: get_particle_processor
    INTEGER :: iproc, coords(c_ndims)

    get_particle_processor = -1
    coords = -1

    ! This could be replaced by a bisection method, but for the moment I
    ! just don't care

    DO iproc = 1, nprocx
      IF (a_particle%part_pos(1) .GE. x_mins(iproc) - dx / 2.0_num &
          .AND. a_particle%part_pos(1) .LT. x_maxs(iproc) + dx / 2.0_num) THEN
        coords(c_ndims) = iproc - 1
        EXIT
      ENDIF
    ENDDO

    DO iproc = 1, nprocy
      IF (a_particle%part_pos(2) .GE. y_mins(iproc) - dy / 2.0_num &
          .AND. a_particle%part_pos(2) .LT. y_maxs(iproc) + dy / 2.0_num) THEN
        coords(c_ndims-1) = iproc - 1
        EXIT
      ENDIF
    ENDDO

    DO iproc = 1, nprocz
      IF (a_particle%part_pos(3) .GE. z_mins(iproc) - dz / 2.0_num &
          .AND. a_particle%part_pos(3) .LT. z_maxs(iproc) + dz / 2.0_num) THEN
        coords(c_ndims-2) = iproc - 1
        EXIT
      ENDIF
    ENDDO

    IF (MINVAL(coords) .LT. 0) THEN
      WRITE(*,*) 'UNLOCATABLE PARTICLE', coords
      RETURN
    ENDIF
    CALL MPI_CART_RANK(comm, coords, get_particle_processor, errcode)
    ! IF (get_particle_processor .NE. rank) PRINT *,

  END FUNCTION get_particle_processor



  ! This subroutine is used to rearrange particles over processors
  SUBROUTINE distribute_particles

    ! This subroutine moves particles which are on the wrong processor
    ! to the correct processor.

    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_send
    TYPE(particle_list), DIMENSION(:), ALLOCATABLE :: pointers_recv
    TYPE(particle), POINTER :: current, next
    INTEGER :: part_proc, iproc, ispecies, ierr
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: sendcounts, recvcounts

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    ALLOCATE(sendcounts(0:nproc-1), recvcounts(0:nproc-1))

    DO ispecies = 1, n_species
      current=>species_list(ispecies)%attached_list%head
      DO iproc = 0, nproc - 1
        CALL create_empty_partlist(pointers_send(iproc))
        CALL create_empty_partlist(pointers_recv(iproc))
      ENDDO

      DO WHILE(ASSOCIATED(current))
        next=>current%next
        part_proc = get_particle_processor(current)
        IF (part_proc .LT. 0) THEN
          PRINT *, 'Unlocatable particle on processor', rank, current%part_pos
          CALL MPI_ABORT(comm, errcode, ierr)
          STOP
        ENDIF
#ifdef PARTICLE_DEBUG
        current%processor = part_proc
#endif
        IF (part_proc .NE. rank) THEN
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, current)
          CALL add_particle_to_partlist(pointers_send(part_proc), current)
        ENDIF
        current=>next
      ENDDO

      DO iproc = 0, nproc - 1
        sendcounts(iproc) = pointers_send(iproc)%count
      ENDDO

      CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER8, recvcounts, 1, &
          MPI_INTEGER8, comm, errcode)

      CALL redblack(pointers_send, pointers_recv, sendcounts, recvcounts)

      DO iproc = 0, nproc - 1
        CALL append_partlist(species_list(ispecies)%attached_list, &
            pointers_recv(iproc))
      ENDDO
    ENDDO

    DEALLOCATE(sendcounts, recvcounts)
    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles

END MODULE balance
