MODULE multigrid

  USE mpi

  IMPLICIT NONE

  PRIVATE
  INTEGER, PARAMETER :: num = KIND(1.d0)

  TYPE t_grid
    REAL(num), DIMENSION(:), POINTER :: x, y
    REAL(num), DIMENSION(:,:), POINTER :: phi, rho, resid
    REAL(num) :: dx, dy
    INTEGER :: nx, ny, gridlevel
    LOGICAL :: converged
    LOGICAL :: finest
  END TYPE t_grid

  INTEGER :: nx = 1000, ny = 1000
  REAL(num) :: dx, dy
  INTEGER :: ngrids = 6
  INTEGER :: nu = 10
  INTEGER :: mu = 1000
  INTEGER :: eta
  TYPE(t_grid), DIMENSION(:), ALLOCATABLE :: grids
  REAL(num), PARAMETER :: length_x = 1.0_num, length_y = 1.0_num
  INTEGER :: left, right, up, down, comm
  INTEGER :: mpireal
  LOGICAL :: mpi_on
  REAL(num), PARAMETER :: epsilon0 = 8.85418782e-12_num

  PUBLIC :: setup_multigrid, run_multigrid, destroy_multigrid, setup_mpi

CONTAINS

  SUBROUTINE gauss_seidel(grid, its_max, boundary)

    TYPE(t_grid), INTENT(INOUT) :: grid
    INTEGER, INTENT(IN) :: its_max, boundary
    INTEGER :: nx_l, ny_l
    INTEGER :: ix, iy, sweep, cycles
    REAL(num) :: lamdax, lamday
    REAL(num) :: w = 1.95_num
    REAL(num) :: sum_local, sum_global
    ! Use this line to set simple clamped boundaries
    REAL(num), DIMENSION(4) :: b_val = (/0.0_num, 0.0_num, 0.0_num, 0.0_num/)
    INTEGER, SAVE :: sweep_no = 0
    INTEGER :: mnx, mny, mxx, mxy, ordx, ordy, val, errcode

    nx_l = grid%nx
    ny_l = grid%ny

    sweep_no = sweep_no + 1

    val = MOD(sweep_no, 4)
    val = 0

    IF (val .EQ. 0) THEN
      mnx = 1
      mny = 1
      mxx = nx_l
      mxy = ny_l
      ordx = 1
      ordy = 1
    ELSE IF (val .EQ. 1) THEN
      mnx = nx_l
      mny = 1
      mxx = 1
      mxy = ny_l
      ordx = -1
      ordy = 1
    ELSE IF (val .EQ. 2) THEN
      mnx = 1
      mny = ny_l
      mxx = nx_l
      mxy = 1
      ordx = 1
      ordy = -1
    ELSE IF (val .EQ. 3) THEN
      mnx = nx_l
      mny = ny_l
      mxx = 1
      mxy = 1
      ordx = -1
      ordy = -1
    ENDIF

    grid%converged = .FALSE.

    CALL boundaries(grid, b_val)
    lamdax = 1.0_num/(grid%dx**2)
    lamday = 1.0_num/(grid%dy**2)

    DO cycles = 1, its_max
      DO sweep = 0, 1
        DO iy = mny, mxy, ordy
          DO ix = mnx, mxx, ordx
            IF (MOD(ix+iy, 2) .EQ. sweep) THEN
              grid%phi(ix, iy) = (1.0-w)* grid%phi(ix, iy) + &
                  w * ((-grid%rho(ix, iy)+lamdax*(grid%phi(ix-1, iy)+&
                  grid%phi(ix+1, iy)) +lamday*(grid%phi(ix, iy-1)+&
                  grid%phi(ix, iy+1)))/(2.0_num*(lamdax+lamday)))
            ENDIF
          ENDDO
        ENDDO
        CALL boundaries(grid, b_val)
      ENDDO

      DO iy = 1, ny_l
        DO ix = 1, nx_l
          grid%resid(ix, iy) = grid%rho(ix, iy) - &
              (lamdax*(grid%phi(ix-1, iy)-2.0_num*grid%phi(ix, iy)+&
              grid%phi(ix+1, iy))+ lamday*(grid%phi(ix, iy-1)-&
              2.0_num*grid%phi(ix, iy)+grid%phi(ix, iy+1)))
        ENDDO
      ENDDO

      sum_local = MAXVAL(ABS(grid%resid(1:nx_l, 1:ny_l)))
      sum_global = sum_local
      IF (mpi_on) CALL MPI_ALLREDUCE(sum_local, sum_global, 1, mpireal, &
          MPI_MAX, comm, errcode)
      IF (sum_global .LT. 1.0e-6_num) THEN
        grid%converged = .TRUE.
        EXIT
      ELSE
      ENDIF
    ENDDO

  END SUBROUTINE gauss_seidel



  SUBROUTINE boundaries(grid, b_val)

    REAL(num), DIMENSION(4), INTENT(IN) :: b_val
    TYPE(t_grid), INTENT(INOUT) :: grid
    INTEGER :: nx, ny
    INTEGER :: tag = 100, errcode
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status

    nx = grid%nx
    ny = grid%ny

    ! At the moment, only support zero potential boundaries, this is easy to
    ! improve, but not done yet
    IF (grid%gridlevel .NE. 1 .AND. .NOT. grid%finest) THEN
      ! On all but finest grids, use clamped boundaries
      IF (right .EQ. MPI_PROC_NULL) &
          grid%phi(1:nx, ny+1)   = grid%phi(1:nx, ny) *0.0_num
      IF (left .EQ. MPI_PROC_NULL) &
          grid%phi(1:nx, 0)      = grid%phi(1:nx, 1)  *0.0_num
      IF (down .EQ. MPI_PROC_NULL) &
          grid%phi(0, 1:ny)      = grid%phi(1, 1:ny)  *0.0_num
      IF (up .EQ. MPI_PROC_NULL) &
          grid%phi(nx+1, 1:ny)   = grid%phi(nx, 1:ny) *0.0_num
    ELSE
      ! On finest grid, use real boundaries
      IF (right .EQ. MPI_PROC_NULL) grid%phi(1:nx, ny+1)   = b_val(1)
      IF (left  .EQ. MPI_PROC_NULL) grid%phi(1:nx, 0)      = b_val(2)
      IF (down  .EQ. MPI_PROC_NULL) grid%phi(0, 1:ny)      = b_val(3)
      IF (up    .EQ. MPI_PROC_NULL) grid%phi(nx+1, 1:ny)   = b_val(4)
      ! Now do MPI
      IF (mpi_on) THEN
        CALL MPI_SENDRECV(grid%phi(1:nx, ny), nx, mpireal, up, tag, &
            grid%phi(1:nx, 0), nx, mpireal, down, tag, comm, status, errcode)
        CALL MPI_SENDRECV(grid%phi(1:nx, 1), nx, mpireal, down, tag, &
            grid%phi(1:nx, ny+1), nx, mpireal, up, tag, comm, status, errcode)

        CALL MPI_SENDRECV(grid%phi(nx, 1:ny), ny, mpireal, right, tag, &
            grid%phi(0, 1:ny), ny, mpireal, left, tag, comm, status, errcode)
        CALL MPI_SENDRECV(grid%phi(1, 1:ny), ny, mpireal, left, tag, &
            grid%phi(nx+1, 1:ny), ny, mpireal, right, tag, comm, &
            status, errcode)
      ENDIF
    ENDIF

  END SUBROUTINE boundaries



  SUBROUTINE restrict(hires, lores, grid_hi, grid_lo)

    ! Restriction with full weighting
    REAL(num), DIMENSION(:,:), INTENT(INOUT)  :: hires
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: lores
    TYPE(t_grid), INTENT(IN) :: grid_hi, grid_lo
    INTEGER :: ix, iy, nx_h, ny_h, nx_l, ny_l

    nx_h = grid_hi%nx
    ny_h = grid_hi%ny
    nx_l = grid_hi%nx/2
    ny_l = grid_hi%ny/2

    DO ix = 1, nx_l
      DO iy = 1, ny_l
        lores(ix, iy) = (4.0_num*hires(ix*2-1, iy*2-1) + &
            2.0_num * (hires(ix*2, iy*2-1)+hires(ix*2-2, iy*2-1) + &
            hires(ix*2-1, iy*2)+hires(ix*2-1, iy*2-2)) + &
            (hires(ix*2, iy*2) + hires(ix*2, iy*2-2) + hires(ix*2-2, iy*2) + &
            hires(ix*2-2, iy*2-2)))/16.0_num
      ENDDO
    ENDDO

  END SUBROUTINE restrict



  SUBROUTINE prolong(lores, hires, grid_lo, grid_hi, add)

    ! Prolongation using inverse full weighting
    REAL(num), DIMENSION(:,:), INTENT(IN)  :: lores
    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: hires
    TYPE(t_grid), INTENT(IN) :: grid_lo, grid_hi
    LOGICAL, INTENT(IN) :: add
    INTEGER :: nx_l, ny_l
    INTEGER :: ix, iy, nx_h, ny_h
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temp

    nx_l = grid_lo%nx
    ny_l = grid_lo%ny
    nx_h = grid_lo%nx*2
    ny_h = grid_lo%ny*2

    ALLOCATE(temp(0:nx_h+1, 0:ny_h+1))
    temp = 0.0_num

    DO ix = 1, nx_l
      DO iy = 1, ny_l
        temp(2*ix-1, 2*iy-1) = lores(ix, iy)
        temp(2*ix, 2*iy-1)   = 0.5_num*(lores(ix, iy)+lores(ix+1, iy))
        temp(2*ix-1, 2*iy)   = 0.5_num*(lores(ix, iy)+lores(ix, iy+1))
        temp(2*ix, 2*iy)    = 0.25_num*(lores(ix, iy)+lores(ix+1, iy)+ &
            lores(ix, iy+1)+lores(ix+1, iy+1))
      ENDDO
    ENDDO

    IF (add) THEN
      hires = hires + temp
    ELSE
      hires = temp
    ENDIF

  END SUBROUTINE prolong



  SUBROUTINE v_cycle(start_grid)

    INTEGER, INTENT(IN) :: start_grid
    INTEGER :: count

    count = 0
    DO WHILE (.NOT. (grids(start_grid)%converged) .AND. count .LE. eta)
      CALL do_grid(grids(start_grid))
      count = count+1
    ENDDO

  END SUBROUTINE v_cycle



  RECURSIVE SUBROUTINE do_grid(grid)

    TYPE(t_grid), INTENT(INOUT) :: grid
    INTEGER :: curlevel
    INTEGER :: bound = 1

    IF (grid%gridlevel .EQ. 1) bound = 1
    curlevel = grid%gridlevel
    IF (grid%gridlevel .NE. ngrids) THEN
      ! Not on coarsest grid to relax and then copy down
      ! Relax
      CALL gauss_seidel(grid, nu, bound)
      ! If converged at this level, don't need to multigrid it further
      IF (.NOT. grid%converged) THEN
        ! Copy down
        grids(curlevel+1)%phi = 0.0_num
        CALL restrict(grid%resid, grids(curlevel+1)%rho, grid, &
            grids(curlevel+1))
        ! Then recursively call this routine on the coarser grid
        CALL do_grid(grids(curlevel+1))
        ! Now, after that do prolongation to get back errors
        CALL prolong(grids(curlevel+1)%phi, grid%phi, grids(curlevel+1), &
            grid, .TRUE.)
        ! Now relax again
        CALL gauss_seidel(grid, nu, bound)
      ENDIF
    ELSE
      ! Are on coarsest grid so should solve directly.
      ! I don't have a plugin LU decomp solver, so just do using a lot
      ! of gauss_seidel
      grid%converged = .FALSE.
      CALL gauss_seidel(grid, mu, bound)
      ! Do nothing else on coarsest grid
    ENDIF

  END SUBROUTINE do_grid



  SUBROUTINE allocate_grids

    INTEGER :: igrid

    ALLOCATE(grids(1:ngrids))

    DO igrid = 1, ngrids
      CALL allocate_single_grid(grids(igrid), igrid-1, dx, dy)
    ENDDO

  END SUBROUTINE allocate_grids



  SUBROUTINE allocate_single_grid(grid, level, dx, dy)

    INTEGER :: nx_l, ny_l
    TYPE(t_grid), INTENT(INOUT) :: grid
    INTEGER, INTENT(IN) :: level
    REAL(num), INTENT(IN) :: dx, dy

    nx_l = nx/(2**level)
    ny_l = ny/(2**level)

    grid%nx = nx_l
    grid%ny = ny_l
    grid%dx = dx*(2.0_num**level)
    grid%dy = dy*(2.0_num**level)
    grid%gridlevel = level+1
    grid%converged = .FALSE.

    ALLOCATE(grid%phi(0:nx_l+1, 0:ny_l+1), grid%rho(0:nx_l+1, 0:ny_l+1))
    ALLOCATE(grid%resid(0:nx_l+1, 0:ny_l+1))

  END SUBROUTINE allocate_single_grid



  SUBROUTINE setup_multigrid(nx_local, ny_local, ngrid_local, dx_local, &
      dy_local)

    INTEGER, INTENT(IN) :: nx_local, ny_local, ngrid_local
    REAL(num), INTENT(IN) :: dx_local, dy_local
    INTEGER :: nx_rec, ny_rec

    nx = nx_local
    ny = ny_local
    dx = dx_local
    dy = dy_local
    ngrids = ngrid_local

    nx_rec = (nx_local/2**(ngrids-1))*(2**(ngrids-1))
    ny_rec = (ny_local/2**(ngrids-1))*(2**(ngrids-1))

    IF (nx_rec .NE. nx_local .OR. ny_rec .NE. ny_local) THEN
      PRINT *, "Unable to form grid chain. Please try using fewer grids."
      STOP
    ENDIF

    IF (ALLOCATED(grids)) THEN
      PRINT *, "Must destroy existing multigrid before calling setup again."
      STOP
    ENDIF

    up = MPI_PROC_NULL
    down = MPI_PROC_NULL
    left = MPI_PROC_NULL
    right = MPI_PROC_NULL
    comm = MPI_COMM_NULL
    mpi_on = .FALSE.

    CALL allocate_grids

  END SUBROUTINE setup_multigrid



  SUBROUTINE setup_mpi(left_l, right_l, up_l, down_l, comm_l)

    INTEGER, INTENT(IN) :: left_l, right_l, up_l, down_l
    INTEGER, INTENT(IN) :: comm_l

    left = left_l
    right = right_l
    up = up_l
    down = down_l
    comm = comm_l

    IF (num .EQ. 4) mpireal = MPI_REAL
    IF (num .EQ. 8) mpireal = MPI_DOUBLE_PRECISION
    mpi_on = .TRUE.

  END SUBROUTINE setup_mpi



  SUBROUTINE run_multigrid(phi, rho, n_v_cycles, n_gs_cycles, n_c_cycles, &
      force_converged)

    REAL(num), DIMENSION(:,:), INTENT(INOUT) :: phi
    REAL(num), DIMENSION(:,:), INTENT(IN) :: rho
    INTEGER, INTENT(IN) :: n_v_cycles, n_gs_cycles, n_c_cycles
    INTEGER :: igrid
    LOGICAL, INTENT(INOUT) :: force_converged

    IF (.NOT. ALLOCATED(grids)) THEN
      PRINT *, "Must setup_multigrid before running"
      STOP
    ENDIF

    mu = n_c_cycles
    nu = n_gs_cycles
    eta = n_v_cycles

    grids(1)%phi(1:nx, 1:ny) = phi
    grids(1)%rho(1:nx, 1:ny) = rho

    ! First part of FMC simply copies the information down to the finest grids
    ! This is used to precondition the data
    grids(:)%finest = .TRUE.
    DO igrid = 1, ngrids-1
      CALL restrict(grids(igrid)%phi, grids(igrid+1)%phi, grids(igrid), &
          grids(igrid+1))
      CALL restrict(grids(igrid)%rho, grids(igrid+1)%rho, grids(igrid), &
          grids(igrid+1))
    ENDDO

    ! The use V_Cycles to move up towards the finest grid
    DO igrid = ngrids, 2, -1
      CALL v_cycle(igrid)
      grids(igrid)%finest = .FALSE.
      CALL prolong(grids(igrid)%phi, grids(igrid-1)%phi, grids(igrid), &
          grids(igrid-1), .FALSE.)
    ENDDO

    ! Finally v_cycle on finest grid
    IF (force_converged) THEN
      DO WHILE (.NOT. grids(1)%converged)
        CALL v_cycle(1)
      ENDDO
    ELSE
      CALL v_cycle(1)
    ENDIF

    ! Should now check that convergence has gone as Cauchy Sequence, but
    ! implement later

    ! Copy the data back
    phi = grids(1)%phi(1:nx, 1:ny)
    force_converged = grids(1)%converged

  END SUBROUTINE run_multigrid



  SUBROUTINE destroy_multigrid

    ! This subroutine deallocates all the memory used by the multigrid system
    ! Don't call this unless you've either finished or need the memory
    ! The setup routines are quite expensive

    INTEGER :: igrid

    DO igrid = 1, ngrids
      DEALLOCATE(grids(igrid)%rho, grids(igrid)%phi, grids(igrid)%resid)
    ENDDO
    DEALLOCATE(grids)

  END SUBROUTINE destroy_multigrid

END MODULE multigrid
