MODULE balance

  USE mpi
  USE partlist
  USE boundary
  USE mpi_subtype_control

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
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_x, ends_x
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_y, ends_y
    INTEGER, DIMENSION(:), ALLOCATABLE :: starts_z, ends_z
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

    ALLOCATE(starts_x(1:nprocx), ends_x(1:nprocx))
    ALLOCATE(starts_y(1:nprocy), ends_y(1:nprocy))
    ALLOCATE(starts_z(1:nprocz), ends_z(1:nprocz))

    ! Sweep in X
    IF (IAND(balance_mode, c_lb_x) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in X
      ALLOCATE(load_x(nx_global))
      CALL get_load_in_x(load_x)
      CALL calculate_breaks(load_x, nprocx, starts_x, ends_x)
    ELSE
      ! Just keep the original lengths
      starts_x = cell_x_min
      ends_x = cell_x_max
    ENDIF

    ! Sweep in Y
    IF (IAND(balance_mode, c_lb_y) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in Y
      ALLOCATE(load_y(ny_global))
      CALL get_load_in_y(load_y)
      CALL calculate_breaks(load_y, nprocy, starts_y, ends_y)
    ELSE
      ! Just keep the original lengths
      starts_y = cell_y_min
      ends_y = cell_y_max
    ENDIF

    ! Sweep in Z
    IF (IAND(balance_mode, c_lb_z) .NE. 0 &
        .OR. IAND(balance_mode, c_lb_auto) .NE. 0) THEN
      ! Rebalancing in Z
      ALLOCATE(load_z(nz_global))
      CALL get_load_in_z(load_z)
      CALL calculate_breaks(load_z, nprocz, starts_z, ends_z)
    ELSE
      ! Just keep the original lengths
      starts_z = cell_z_min
      ends_z = cell_z_max
    ENDIF

    ! In the autobalancer then determine whether to balance in X, Y or Z
    ! Is this worth keeping?
    IF (IAND(balance_mode, c_lb_auto) .NE. 0 ) THEN

      ! Code is auto load balancing
      max_x = 0
      min_x = npart_global
      DO iproc = 1, nprocx
        wk = SUM(load_x(starts_x(iproc):ends_x(iproc)))
        IF (wk .GT. max_x) max_x = wk
        IF (wk .LT. min_x) min_x = wk
      ENDDO

      max_y = 0
      min_y = npart_global
      DO iproc = 1, nprocy
        wk = SUM(load_y(starts_y(iproc):ends_y(iproc)))
        IF (wk .GT. max_y) max_y = wk
        IF (wk .LT. min_y) min_y = wk
      ENDDO

      max_z = 0
      min_z = npart_global
      DO iproc = 1, nprocz
        wk = SUM(load_z(starts_z(iproc):ends_z(iproc)))
        IF (wk .GT. max_z) max_z = wk
        IF (wk .LT. min_z) min_z = wk
      ENDDO

      balance_frac_x = REAL(min_x, num) / REAL(max_x, num)
      balance_frac_y = REAL(min_y, num) / REAL(max_y, num)
      balance_frac_z = REAL(min_z, num) / REAL(max_z, num)

      IF (balance_frac_x .LT. balance_frac_y) THEN
        IF (balance_frac_x .LT. balance_frac_z) THEN
          starts_x = cell_x_min
          ends_x = cell_x_max
        ELSE
          starts_z = cell_z_min
          ends_z = cell_z_max
        ENDIF
      ELSE
        IF (balance_frac_y .LT. balance_frac_z) THEN
          starts_y = cell_y_min
          ends_y = cell_y_max
        ELSE
          starts_z = cell_z_min
          ends_z = cell_z_max
        ENDIF
      ENDIF

    ENDIF

    IF (ALLOCATED(load_x)) DEALLOCATE(load_x)
    IF (ALLOCATED(load_y)) DEALLOCATE(load_y)
    IF (ALLOCATED(load_z)) DEALLOCATE(load_z)

    ! Now need to calculate the start and end points for the new domain on
    ! the current processor

    domain(1,:) = (/starts_x(x_coords+1), ends_x(x_coords+1)/)
    domain(2,:) = (/starts_y(y_coords+1), ends_y(y_coords+1)/)
    domain(3,:) = (/starts_z(z_coords+1), ends_z(z_coords+1)/)

    ! Redistribute the field variables
    CALL redistribute_fields(domain)

    ! Copy the new lengths into the permanent variables
    cell_x_min = starts_x
    cell_y_min = starts_y
    cell_z_min = starts_z
    cell_x_max = ends_x
    cell_y_max = ends_y
    cell_z_max = ends_z

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
    DO iproc = 0, nprocx - 1
      x_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    ! Same for y
    DO iproc = 0, nprocy - 1
      y_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO
    ! Same for z
    DO iproc = 0, nprocz - 1
      z_mins(iproc) = z_global(cell_z_min(iproc+1))
      z_maxs(iproc) = z_global(cell_z_max(iproc+1))
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

    ! This subroutine redistributes the 2D field variables over the new
    ! processor layout. If using a 2D field of your own then set the
    ! redistribute_field subroutine to implement it. 1D fields, you're on
    ! your own (have global copies and use those to repopulate?)

    INTEGER :: nx_new, ny_new, nz_new
    INTEGER, DIMENSION(c_ndims,2), INTENT(IN) :: new_domain
    REAL(num), DIMENSION(:,:,:,:), ALLOCATABLE :: temp4d
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: temp
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp2d
    TYPE(laser_block), POINTER :: current
    INTEGER :: ispecies, index, n_species_local

    nx_new = new_domain(1,2) - new_domain(1,1) + 1
    ny_new = new_domain(2,2) - new_domain(2,1) + 1
    nz_new = new_domain(3,2) - new_domain(3,1) + 1

    ALLOCATE(temp(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))

    temp = 0.0_num
    CALL redistribute_field(new_domain, ex, temp)
    DEALLOCATE(ex)
    ALLOCATE(ex(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    ex = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ey, temp)
    DEALLOCATE(ey)
    ALLOCATE(ey(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    ey = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, ez, temp)
    DEALLOCATE(ez)
    ALLOCATE(ez(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    ez = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bx, temp)
    DEALLOCATE(bx)
    ALLOCATE(bx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    bx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, by, temp)
    DEALLOCATE(by)
    ALLOCATE(by(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    by = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, bz, temp)
    DEALLOCATE(bz)
    ALLOCATE(bz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    bz = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jx, temp)
    DEALLOCATE(jx)
    ALLOCATE(jx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    jx = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jy, temp)
    DEALLOCATE(jy)
    ALLOCATE(jy(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    jy = temp

    temp = 0.0_num
    CALL redistribute_field(new_domain, jz, temp)
    DEALLOCATE(jz)
    ALLOCATE(jz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
    jz = temp

    IF (cpml_boundaries) THEN
      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_eyx, temp)
      DEALLOCATE(cpml_psi_eyx)
      ALLOCATE(cpml_psi_eyx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_eyx = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_ezx, temp)
      DEALLOCATE(cpml_psi_ezx)
      ALLOCATE(cpml_psi_ezx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_ezx = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_byx, temp)
      DEALLOCATE(cpml_psi_byx)
      ALLOCATE(cpml_psi_byx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_byx = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_bzx, temp)
      DEALLOCATE(cpml_psi_bzx)
      ALLOCATE(cpml_psi_bzx(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_bzx = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_exy, temp)
      DEALLOCATE(cpml_psi_exy)
      ALLOCATE(cpml_psi_exy(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_exy = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_ezy, temp)
      DEALLOCATE(cpml_psi_ezy)
      ALLOCATE(cpml_psi_ezy(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_ezy = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_bxy, temp)
      DEALLOCATE(cpml_psi_bxy)
      ALLOCATE(cpml_psi_bxy(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_bxy = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_bzy, temp)
      DEALLOCATE(cpml_psi_bzy)
      ALLOCATE(cpml_psi_bzy(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_bzy = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_exz, temp)
      DEALLOCATE(cpml_psi_exz)
      ALLOCATE(cpml_psi_exz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_exz = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_eyz, temp)
      DEALLOCATE(cpml_psi_eyz)
      ALLOCATE(cpml_psi_eyz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_eyz = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_bxz, temp)
      DEALLOCATE(cpml_psi_bxz)
      ALLOCATE(cpml_psi_bxz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_bxz = temp

      temp = 0.0_num
      CALL redistribute_field(new_domain, cpml_psi_byz, temp)
      DEALLOCATE(cpml_psi_byz)
      ALLOCATE(cpml_psi_byz(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3))
      cpml_psi_byz = temp

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

      ALLOCATE(temp4d(-2:nx_new+3, -2:ny_new+3, -2:nz_new+3, 1:n_species_local))

      DO ispecies = 1, n_species_local
        CALL redistribute_field(new_domain, &
            averaged_data(index)%array(:,:,:,ispecies), temp4d(:,:,:,ispecies))
      ENDDO

      DEALLOCATE(averaged_data(index)%array)
      ALLOCATE(averaged_data(index)&
          %array(-2:nx_new+3,-2:ny_new+3,-2:nz_new+3,n_species_local))

      averaged_data(index)%array = temp4d

      DEALLOCATE(temp4d)
    ENDDO

    ALLOCATE(temp2d(-2:ny_new+3, -2:nz_new+3))

    current=>laser_x_min
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    current=>laser_x_max
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:ny_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_x, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:ny_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, ex_x_min, temp2d)
    DEALLOCATE(ex_x_min)
    ALLOCATE(ex_x_min(-2:ny_new+3, -2:nz_new+3))
    ex_x_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, ex_x_max, temp2d)
    DEALLOCATE(ex_x_max)
    ALLOCATE(ex_x_max(-2:ny_new+3, -2:nz_new+3))
    ex_x_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, ey_x_min, temp2d)
    DEALLOCATE(ey_x_min)
    ALLOCATE(ey_x_min(-2:ny_new+3, -2:nz_new+3))
    ey_x_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, ey_x_max, temp2d)
    DEALLOCATE(ey_x_max)
    ALLOCATE(ey_x_max(-2:ny_new+3, -2:nz_new+3))
    ey_x_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, ez_x_min, temp2d)
    DEALLOCATE(ez_x_min)
    ALLOCATE(ez_x_min(-2:ny_new+3, -2:nz_new+3))
    ez_x_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, ez_x_max, temp2d)
    DEALLOCATE(ez_x_max)
    ALLOCATE(ez_x_max(-2:ny_new+3, -2:nz_new+3))
    ez_x_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, bx_x_min, temp2d)
    DEALLOCATE(bx_x_min)
    ALLOCATE(bx_x_min(-2:ny_new+3, -2:nz_new+3))
    bx_x_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, bx_x_max, temp2d)
    DEALLOCATE(bx_x_max)
    ALLOCATE(bx_x_max(-2:ny_new+3, -2:nz_new+3))
    bx_x_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, by_x_min, temp2d)
    DEALLOCATE(by_x_min)
    ALLOCATE(by_x_min(-2:ny_new+3, -2:nz_new+3))
    by_x_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, by_x_max, temp2d)
    DEALLOCATE(by_x_max)
    ALLOCATE(by_x_max(-2:ny_new+3, -2:nz_new+3))
    by_x_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, bz_x_min, temp2d)
    DEALLOCATE(bz_x_min)
    ALLOCATE(bz_x_min(-2:ny_new+3, -2:nz_new+3))
    bz_x_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_x, bz_x_max, temp2d)
    DEALLOCATE(bz_x_max)
    ALLOCATE(bz_x_max(-2:ny_new+3, -2:nz_new+3))
    bz_x_max = temp2d

    DEALLOCATE(temp2d)
    ALLOCATE(temp2d(-2:nx_new+3, -2:nz_new+3))

    current=>laser_y_min
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    current=>laser_y_max
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:nz_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_y, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:nz_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, ex_y_min, temp2d)
    DEALLOCATE(ex_y_min)
    ALLOCATE(ex_y_min(-2:nx_new+3, -2:nz_new+3))
    ex_y_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, ex_y_max, temp2d)
    DEALLOCATE(ex_y_max)
    ALLOCATE(ex_y_max(-2:nx_new+3, -2:nz_new+3))
    ex_y_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, ey_y_min, temp2d)
    DEALLOCATE(ey_y_min)
    ALLOCATE(ey_y_min(-2:nx_new+3, -2:nz_new+3))
    ey_y_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, ey_y_max, temp2d)
    DEALLOCATE(ey_y_max)
    ALLOCATE(ey_y_max(-2:nx_new+3, -2:nz_new+3))
    ey_y_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, ez_y_min, temp2d)
    DEALLOCATE(ez_y_min)
    ALLOCATE(ez_y_min(-2:nx_new+3, -2:nz_new+3))
    ez_y_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, ez_y_max, temp2d)
    DEALLOCATE(ez_y_max)
    ALLOCATE(ez_y_max(-2:nx_new+3, -2:nz_new+3))
    ez_y_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, bx_y_min, temp2d)
    DEALLOCATE(bx_y_min)
    ALLOCATE(bx_y_min(-2:nx_new+3, -2:nz_new+3))
    bx_y_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, bx_y_max, temp2d)
    DEALLOCATE(bx_y_max)
    ALLOCATE(bx_y_max(-2:nx_new+3, -2:nz_new+3))
    bx_y_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, by_y_min, temp2d)
    DEALLOCATE(by_y_min)
    ALLOCATE(by_y_min(-2:nx_new+3, -2:nz_new+3))
    by_y_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, by_y_max, temp2d)
    DEALLOCATE(by_y_max)
    ALLOCATE(by_y_max(-2:nx_new+3, -2:nz_new+3))
    by_y_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, bz_y_min, temp2d)
    DEALLOCATE(bz_y_min)
    ALLOCATE(bz_y_min(-2:nx_new+3, -2:nz_new+3))
    bz_y_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_y, bz_y_max, temp2d)
    DEALLOCATE(bz_y_max)
    ALLOCATE(bz_y_max(-2:nx_new+3, -2:nz_new+3))
    bz_y_max = temp2d

    DEALLOCATE(temp2d)
    ALLOCATE(temp2d(-2:nx_new+3, -2:ny_new+3))

    current=>laser_z_min
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:ny_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:ny_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    current=>laser_z_max
    DO WHILE(ASSOCIATED(current))
      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%profile, temp2d)
      DEALLOCATE(current%profile)
      ALLOCATE(current%profile(-2:nx_new+3, -2:ny_new+3))
      current%profile = temp2d

      temp2d = 0.0_num
      CALL redistribute_field_2d(new_domain, c_dir_z, current%phase, temp2d)
      DEALLOCATE(current%phase)
      ALLOCATE(current%phase(-2:nx_new+3, -2:ny_new+3))
      current%phase = temp2d

      current=>current%next
    ENDDO

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, ex_z_min, temp2d)
    DEALLOCATE(ex_z_min)
    ALLOCATE(ex_z_min(-2:nx_new+3, -2:ny_new+3))
    ex_z_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, ex_z_max, temp2d)
    DEALLOCATE(ex_z_max)
    ALLOCATE(ex_z_max(-2:nx_new+3, -2:ny_new+3))
    ex_z_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, ey_z_min, temp2d)
    DEALLOCATE(ey_z_min)
    ALLOCATE(ey_z_min(-2:nx_new+3, -2:ny_new+3))
    ey_z_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, ey_z_max, temp2d)
    DEALLOCATE(ey_z_max)
    ALLOCATE(ey_z_max(-2:nx_new+3, -2:ny_new+3))
    ey_z_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, ez_z_min, temp2d)
    DEALLOCATE(ez_z_min)
    ALLOCATE(ez_z_min(-2:nx_new+3, -2:ny_new+3))
    ez_z_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, ez_z_max, temp2d)
    DEALLOCATE(ez_z_max)
    ALLOCATE(ez_z_max(-2:nx_new+3, -2:ny_new+3))
    ez_z_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, bx_z_min, temp2d)
    DEALLOCATE(bx_z_min)
    ALLOCATE(bx_z_min(-2:nx_new+3, -2:ny_new+3))
    bx_z_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, bx_z_max, temp2d)
    DEALLOCATE(bx_z_max)
    ALLOCATE(bx_z_max(-2:nx_new+3, -2:ny_new+3))
    bx_z_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, by_z_min, temp2d)
    DEALLOCATE(by_z_min)
    ALLOCATE(by_z_min(-2:nx_new+3, -2:ny_new+3))
    by_z_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, by_z_max, temp2d)
    DEALLOCATE(by_z_max)
    ALLOCATE(by_z_max(-2:nx_new+3, -2:ny_new+3))
    by_z_max = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, bz_z_min, temp2d)
    DEALLOCATE(bz_z_min)
    ALLOCATE(bz_z_min(-2:nx_new+3, -2:ny_new+3))
    bz_z_min = temp2d

    temp2d = 0.0_num
    CALL redistribute_field_2d(new_domain, c_dir_z, bz_z_max, temp2d)
    DEALLOCATE(bz_z_max)
    ALLOCATE(bz_z_max(-2:nx_new+3, -2:ny_new+3))
    bz_z_max = temp2d

    DEALLOCATE(temp2d)

    ! Re-distribute moving windows and temperature boundaries
    DO ispecies = 1, n_species

      IF (move_window) THEN
        ALLOCATE(temp2d(-2:ny_new+3,-2:nz_new+3))

        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%density, temp2d)

        DEALLOCATE(species_list(ispecies)%density)
        ALLOCATE(species_list(ispecies)%density(-2:ny_new+3,-2:nz_new+3))

        species_list(ispecies)%density = temp2d

        DEALLOCATE(temp2d)

        ALLOCATE(temp(-2:ny_new+3,-2:nz_new+3,3))

        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%temperature(-2:ny+3,-2:nz+3,1), &
            temp(-2:ny_new+3,-2:nz_new+3,1))
        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%temperature(-2:ny+3,-2:nz+3,2), &
            temp(-2:ny_new+3,-2:nz_new+3,2))
        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%temperature(-2:ny+3,-2:nz+3,3), &
            temp(-2:ny_new+3,-2:nz_new+3,3))

        DEALLOCATE(species_list(ispecies)%temperature)
        ALLOCATE(species_list(ispecies)&
            %temperature(-2:ny_new+3,-2:nz_new+3,1:3))

        species_list(ispecies)%temperature = temp

        DEALLOCATE(temp)
      ENDIF

      IF (bc_particle(c_bd_x_min) .EQ. c_bc_thermal) THEN
        ALLOCATE(temp(-2:ny_new+3,-2:nz_new+3,3))
        temp = 0.0_num

        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_min(-2:ny+3,-2:nz+3,1), &
            temp(-2:ny_new+3,-2:nz_new+3,1))
        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_min(-2:ny+3,-2:nz+3,2), &
            temp(-2:ny_new+3,-2:nz_new+3,2))
        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_min(-2:ny+3,-2:nz+3,3), &
            temp(-2:ny_new+3,-2:nz_new+3,3))

        DEALLOCATE(species_list(ispecies)%ext_temp_x_min)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_x_min(-2:ny_new+3,-2:nz_new+3,3))

        species_list(ispecies)%ext_temp_x_min = temp

        DEALLOCATE(temp)
      ENDIF

      IF (bc_particle(c_bd_x_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(temp(-2:ny_new+3,-2:nz_new+3,3))
        temp = 0.0_num

        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_max(-2:ny+3,-2:nz+3,1), &
            temp(-2:ny_new+3,-2:nz_new+3,1))
        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_max(-2:ny+3,-2:nz+3,2), &
            temp(-2:ny_new+3,-2:nz_new+3,2))
        CALL redistribute_field_2d(new_domain, c_dir_x, &
            species_list(ispecies)%ext_temp_x_max(-2:ny+3,-2:nz+3,3), &
            temp(-2:ny_new+3,-2:nz_new+3,3))

        DEALLOCATE(species_list(ispecies)%ext_temp_x_max)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_x_max(-2:ny_new+3,-2:nz_new+3,3))

        species_list(ispecies)%ext_temp_x_max = temp

        DEALLOCATE(temp)
      ENDIF

      IF (bc_particle(c_bd_y_min) .EQ. c_bc_thermal) THEN
        ALLOCATE(temp(-2:nx_new+3,-2:nz_new+3,3))
        temp = 0.0_num

        CALL redistribute_field_2d(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_min(-2:nx+3,-2:nz+3,1), &
            temp(-2:nx_new+3,-2:nz_new+3,1))
        CALL redistribute_field_2d(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_min(-2:nx+3,-2:nz+3,2), &
            temp(-2:nx_new+3,-2:nz_new+3,2))
        CALL redistribute_field_2d(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_min(-2:nx+3,-2:nz+3,3), &
            temp(-2:nx_new+3,-2:nz_new+3,3))

        DEALLOCATE(species_list(ispecies)%ext_temp_y_min)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_y_min(-2:nx_new+3,-2:nz_new+3,3))

        species_list(ispecies)%ext_temp_y_min = temp

        DEALLOCATE(temp)
      ENDIF

      IF (bc_particle(c_bd_y_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(temp(-2:nx_new+3,-2:nz_new+3,3))
        temp = 0.0_num

        CALL redistribute_field_2d(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_max(-2:nx+3,-2:nz+3,1), &
            temp(-2:nx_new+3,-2:nz_new+3,1))
        CALL redistribute_field_2d(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_max(-2:nx+3,-2:nz+3,2), &
            temp(-2:nx_new+3,-2:nz_new+3,2))
        CALL redistribute_field_2d(new_domain, c_dir_y, &
            species_list(ispecies)%ext_temp_y_max(-2:nx+3,-2:nz+3,3), &
            temp(-2:nx_new+3,-2:nz_new+3,3))

        DEALLOCATE(species_list(ispecies)%ext_temp_y_max)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_y_max(-2:nx_new+3,-2:nz_new+3,3))

        species_list(ispecies)%ext_temp_y_max = temp

        DEALLOCATE(temp)
      ENDIF

      IF (bc_particle(c_bd_z_min) .EQ. c_bc_thermal) THEN
        ALLOCATE(temp(-2:nx_new+3,-2:ny_new+3,3))
        temp = 0.0_num

        CALL redistribute_field_2d(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_min(-2:nx+3,-2:ny+3,1), &
            temp(-2:nx_new+3,-2:ny_new+3,1))
        CALL redistribute_field_2d(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_min(-2:nx+3,-2:ny+3,2), &
            temp(-2:nx_new+3,-2:ny_new+3,2))
        CALL redistribute_field_2d(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_min(-2:nx+3,-2:ny+3,3), &
            temp(-2:nx_new+3,-2:ny_new+3,3))

        DEALLOCATE(species_list(ispecies)%ext_temp_z_min)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_z_min(-2:nx_new+3,-2:ny_new+3,3))

        species_list(ispecies)%ext_temp_z_min = temp

        DEALLOCATE(temp)
      ENDIF

      IF (bc_particle(c_bd_z_max) .EQ. c_bc_thermal) THEN
        ALLOCATE(temp(-2:nx_new+3,-2:ny_new+3,3))
        temp = 0.0_num

        CALL redistribute_field_2d(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_max(-2:nx+3,-2:ny+3,1), &
            temp(-2:nx_new+3,-2:ny_new+3,1))
        CALL redistribute_field_2d(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_max(-2:nx+3,-2:ny+3,2), &
            temp(-2:nx_new+3,-2:ny_new+3,2))
        CALL redistribute_field_2d(new_domain, c_dir_z, &
            species_list(ispecies)%ext_temp_z_max(-2:nx+3,-2:ny+3,3), &
            temp(-2:nx_new+3,-2:ny_new+3,3))

        DEALLOCATE(species_list(ispecies)%ext_temp_z_max)
        ALLOCATE(species_list(ispecies)&
            %ext_temp_z_max(-2:nx_new+3,-2:ny_new+3,3))

        species_list(ispecies)%ext_temp_z_max = temp

        DEALLOCATE(temp)
      ENDIF

    ENDDO

  END SUBROUTINE redistribute_fields



  SUBROUTINE redistribute_field(domain, field_in, field_out)

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

    CALL MPI_FILE_SET_VIEW(fh, offset, subarray_write, subtype_write, &
        'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(fh, field_in, 1, subarray_write, status, errcode)

    CALL MPI_FILE_SET_VIEW(fh, offset, subarray_read, subtype_read, &
        'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_READ_ALL(fh, field_out, 1, subarray_read, status, errcode)

    CALL MPI_FILE_CLOSE(fh, errcode)

    CALL MPI_TYPE_FREE(subarray_write, errcode)
    CALL MPI_TYPE_FREE(subtype_write, errcode)
    CALL MPI_TYPE_FREE(subarray_read, errcode)
    CALL MPI_TYPE_FREE(subtype_read, errcode)

    CALL do_field_mpi_with_lengths(field_out, nx_new, ny_new, nz_new)

  END SUBROUTINE redistribute_field



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
    INTEGER :: subtype_write, subtype_read, fh
    INTEGER :: subarray_write, subarray_read
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

    CALL MPI_TYPE_FREE(subtype_write, errcode)
    CALL MPI_TYPE_FREE(subtype_read, errcode)
    CALL MPI_TYPE_FREE(subarray_write, errcode)
    CALL MPI_TYPE_FREE(subarray_read, errcode)

  END SUBROUTINE redistribute_field_2d



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



  SUBROUTINE calculate_breaks(load, nproc, starts, ends)

    ! This subroutine calculates the places in a given load profile to split
    ! The domain to give the most even subdivision possible

    INTEGER(KIND=8), INTENT(IN), DIMENSION(:) :: load
    INTEGER, INTENT(IN) :: nproc
    INTEGER, DIMENSION(:), INTENT(OUT) :: starts, ends
    INTEGER :: sz, idim, proc, old
    INTEGER(KIND=8) :: total, total_old, load_per_proc_ideal

    sz = SIZE(load)
    ends = sz

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
          ends(proc) = idim - 1
        ELSE
          ends(proc) = idim
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
      IF (ends(proc) - old .LE. 0) THEN
        ends(proc) = old + 1
      ENDIF
      old = ends(proc)
    ENDDO

    ! Backwards
    old = sz + 1
    DO proc = nproc, 1, -1
      IF (old - ends(proc) .LE. 0) THEN
        ends(proc) = old - 1
      ENDIF
      old = ends(proc)
    ENDDO

    ! Set starts
    starts(1) = 1
    DO proc = 2, nproc
      starts(proc) = ends(proc-1) + 1
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

    DO iproc = 0, nprocx - 1
      IF (a_particle%part_pos(1) .GE. x_mins(iproc) - dx / 2.0_num &
          .AND. a_particle%part_pos(1) .LT. x_maxs(iproc) + dx / 2.0_num) THEN
        coords(c_ndims) = iproc
        EXIT
      ENDIF
    ENDDO

    DO iproc = 0, nprocy - 1
      IF (a_particle%part_pos(2) .GE. y_mins(iproc) - dy / 2.0_num &
          .AND. a_particle%part_pos(2) .LT. y_maxs(iproc) + dy / 2.0_num) THEN
        coords(c_ndims-1) = iproc
        EXIT
      ENDIF
    ENDDO

    DO iproc = 0, nprocz - 1
      IF (a_particle%part_pos(3) .GE. z_mins(iproc) - dz / 2.0_num &
          .AND. a_particle%part_pos(3) .LT. z_maxs(iproc) + dz / 2.0_num) THEN
        coords(c_ndims-2) = iproc
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
    INTEGER :: part_proc, iproc, ispecies, i, oldrank, newrank, nodd, nevn, ierr
    INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: sendcounts, recvcounts
    INTEGER, DIMENSION(:), ALLOCATABLE :: oddlist, evnlist

    nodd = (nproc + 1) / 2
    nevn = nproc / 2

    ALLOCATE(pointers_send(0:nproc-1), pointers_recv(0:nproc-1))
    ALLOCATE(sendcounts(0:nproc-1), recvcounts(0:nproc-1))
    ALLOCATE(oddlist(0:nodd-1), evnlist(0:nevn-1))

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

      ! Split the processors into odd and even lists.
      ! Even processors send then receive, Odd processors vice-versa.
      ! Next, on even processors split the even list in two and on
      ! odd processors split the odd list in two.
      ! Repeat until the list size is equal to one.

      ! If the number of processors is not divisible by two then the
      ! odd list has one extra entry.

      nodd = (nproc + 1) / 2
      nevn = nproc / 2
      DO i = 0, nevn - 1
        oddlist(i) = 2 * i
        evnlist(i) = 2 * i + 1
        IF (oddlist(i) .EQ. rank .OR. evnlist(i) .EQ. rank) THEN
          newrank = i
        ENDIF
      ENDDO

      IF (nodd .NE. nevn) oddlist(nodd-1) = nproc - 1

      CALL MPI_ALLTOALL(sendcounts, 1, MPI_INTEGER8, recvcounts, 1, &
          MPI_INTEGER8, comm, errcode)

      oldrank = rank
      DO WHILE(nevn .GT. 0)
        IF (MOD(oldrank,2) .EQ. 0) THEN
          oldrank = newrank
          DO i = 0, nevn - 1
            iproc = evnlist(i)
            IF (sendcounts(iproc) .GT. 0) THEN
              CALL partlist_send_nocount(pointers_send(iproc), iproc)
              CALL destroy_partlist(pointers_send(iproc))
            ENDIF
          ENDDO

          DO i = 0, nevn - 1
            iproc = evnlist(i)
            IF (iproc .GE. nproc) EXIT
            IF (recvcounts(iproc) .GT. 0) THEN
              CALL partlist_recv_nocount(pointers_recv(iproc), iproc, &
                  recvcounts(iproc))
            ENDIF
          ENDDO

          nevn = nodd / 2
          nodd = (nodd + 1) / 2
          DO i = 0, nevn - 1
            oddlist(i) = oddlist(2*i)
            evnlist(i) = oddlist(2*i+1)
            IF (oddlist(i) .EQ. rank .OR. evnlist(i) .EQ. rank) THEN
              newrank = i
            ENDIF
          ENDDO
          IF (nodd .NE. nevn) oddlist(nodd-1) = oddlist(2*nodd-2)

        ELSE
          oldrank = newrank
          DO i = 0, nodd - 1
            iproc = oddlist(i)
            IF (iproc .LT. 0) EXIT
            IF (recvcounts(iproc) .GT. 0) THEN
              CALL partlist_recv_nocount(pointers_recv(iproc), iproc, &
                  recvcounts(iproc))
            ENDIF
          ENDDO

          DO i = 0, nodd - 1
            iproc = oddlist(i)
            IF (iproc .LT. 0) EXIT
            IF (sendcounts(iproc) .GT. 0) THEN
              CALL partlist_send_nocount(pointers_send(iproc), iproc)
              CALL destroy_partlist(pointers_send(iproc))
            ENDIF
          ENDDO

          nodd = (nevn + 1) / 2
          nevn = nevn / 2
          DO i = 0, nevn - 1
            oddlist(i) = evnlist(2*i)
            evnlist(i) = evnlist(2*i+1)
            IF (oddlist(i) .EQ. rank .OR. evnlist(i) .EQ. rank) THEN
              newrank = i
            ENDIF
          ENDDO
          IF (nodd .NE. nevn) oddlist(nodd-1) = evnlist(2*nodd-2)
        ENDIF
      ENDDO

      DO iproc = 0, nproc - 1
        CALL append_partlist(species_list(ispecies)%attached_list, &
            pointers_recv(iproc))
      ENDDO
    ENDDO

    DEALLOCATE(sendcounts, recvcounts, oddlist, evnlist)
    DEALLOCATE(pointers_send, pointers_recv)

  END SUBROUTINE distribute_particles

END MODULE balance
