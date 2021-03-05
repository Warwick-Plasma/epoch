SUBROUTINE particle_sorting()

  dx_bin  = 5 * dx ! Number of cells in bin along x
  dy_bin  = 5 * dy ! Number of cells in bin along y
  
  idx_bin = 1./dx_bin ! Inverse of dx_bin
  idy_bin = 1./dy_bin ! Inverse of dy_bin

  nx_bin  = ceiling((x_grid_max_local - x_grid_min_local) * idx_bin) ! Number of bins along x
  ny_bin  = ceiling((y_grid_max_local - y_grid_min_local) * idy_bin) ! Number of bins along y

  n_bins  = nx_bin * ny_bin ! Total number of bins

  ! Calculate particle positions in terms of the bin co-ordinates

  DO ipart = 1, species_list(ispecies)%attached_list%count
    next => current%next
    part_x = (current%part_pos(1) - x_grid_min_local) * idx_bin
    part_y = (current%part_pos(2) - y_grid_min_local) * idy_bin
    
    ix = floor(part_x) ! x-coordinate of the particle bin
    iy = floor(part_y) ! y-coordinate of the particle bin

    tile_id(ipart) = iy * nx_bin + ix + 1 ! 1-D coordinate of the bins

    num(tile_id(ipart)) = num(tile_id(ipart)) + 1 ! Number of particles in each bin

    current => next

  END DO ! End do-loop for particle position in terms of bin co-ordinates

  k = 0

  ! Determine the stride of particle indices in bins  

  DO i = 1, n_bins
    g_indx(i) = k ! Starting particle index for a particular bin
    k = k + num(i)
  END DO ! End do-loop for the stride of particle indices in bins

  ! Particle sorting in 1-D bins 

  DO ipart = 1, species_list(ispecies)%attached_list%count
    next => current%next
    k = tile_id(ipart)
    g_índx(k) = g_indx(k) + 1 ! Rearranged particle index with respect to the bins

#ifndef PER_SPECIES_WEIGHT
    w(g_indx(k))         = current%weight
    fcx(g_indx(k))       = idty  * w(g_indx(k))
    fcy(g_indx(k))       = idtx  * w(g_indx(k))
    fcz(g_indx(k))       = idtxy * w(g_indx(k))
#endif

#ifndef NO_PARTICLE_PROBES
    init_x(g_indx(k))    = current%part_pos(1)
    init_y(g_indx(k))    = current%part_pos(2)
#endif

#ifdef PER_PARTICLE_CHARGE_MASS
    q(g_indx(k))         = current%charge
    m(g_indx(k))         = current%mass
    mc(g_indx(k))        = c * current%mass
    i_mc(g_indx(k))      = 1.0_num / mc
    cmratio(g_indx(k))   = q(g_indx(k)) * dtfac * i_mc(g_indx(k))
    ccmratio(g_indx(k))  = c * cmratio(g_indx(k))
#ifndef NO_PARTICLE_PROBES
    mc2(g_indx(k))       = c * mc(g_indx(k))
#endif
#endif

    ! Copy the particle properties out for sorting
    x(g_indx(k))         = current%part_pos(1) - x_grid_min_local
    y(g_indx(k))         = current%part_pos(2) - y_grid_min_local
    px(g_indx(k))        = current%part_p(1) * ipart_mc
    py(g_indx(k))        = current%part_p(2) * ipart_mc
    pz(g_indx(k))        = current%part_p(3) * ipart_mc
    pvol(g_indx(k))      = current%pvol
    gamma_rel(g_indx(k)) = SQRT(px(g_indx(k))**2 + py(g_indx(k))**2 + pz(g_indx(k))**2 + 1.0_num)
    root(g_indx(k))      = dtco2 / gamma_rel(g_indx(k))
    current => next
  END DO ! End do-loop for particle sorting

END SUBROUTINE particle_sorting

SUBROUTINE particle_pusher()

  DO i = 1, species_list(ispecies)%attached_list%count, LVEC
    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
      jj = i + j - 1
      x(jj) = x(jj) + px(jj) * root(jj)
      y(jj) = y(jj) + py(jj) * root(jj) 

#ifdef WORK_DONE_INTEGRATED
      ! This is the actual total work done by the fields: Results correspond
      ! with the electron's gamma factor

      root(jj) = cmratio(jj) / gamma_rel(jj)

      tmp_x(j) = px(jj) * root(jj)
      tmp_y(j) = py(jj) * root(jj)
      tmp_z(j) = pz(jj) * root(jj)
#endif

    END DO ! End do-loop for j

  ! Calculate fields at particle positions
  ! Grid cell position as a fraction

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      cell_x_r(j) = x(jj) * idx
      cell_y_r(j) = y(jj) * idy

    END DO ! End do-loop for grid cell position as fraction

  ! Round cell position to nearest cell
  
    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      cell_x1(jj) = FLOOR(cell_x_r(j) + 0.5_num)
      cell_y1(jj) = FLOOR(cell_y_r(j) + 0.5_num)

      cell_x2(j) = FLOOR(cell_x_r(j))
      cell_y2(j) = FLOOR(cell_y_r(j))

    END DO ! End do-loop for nearest cell position
 
  ! Calculate fraction of cell between nearest cell boundary and particle

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
  
      jj = i + j - 1
      cell_frac_x(j) = REAL(cell_x1(jj), num) - cell_x_r(j)
      cell_frac_y(j) = REAL(cell_y1(jj), num) - cell_y_r(j)
    
    END DO ! End do-loop for grid cell position fraction 

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      cfx2(j) = cell_frac_x(j)**2
      cfy2(j) = cell_frac_y(j)**2

    END DO
      
    DO j = 1, 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
      jj = i + j - 1
      cell_x1(jj) = cell_x1(jj) + 1
      cell_y1(jj) = cell_y1(jj) + 1

  ! Particle weight factors as described in Page 25 of the PSC manual
  ! These weight grid properties onto particles
  ! Also used to weight particle properties onto grid, used later to calculate J
  ! NOTE: These weights require an additional multiplication factor

  ! This weighing is for triangle shaped particles
      
      gxx(-1,j) = 0.25_num + cfx2(j) + cell_frac_x(j)
      gxx( 0,j) = 1.5_num - 2.0_num * cfx2(j)
      gxx( 1,j) = 0.25_num + cfx2(j) - cell_frac_x(j)

      gyy(-1,j) = 0.25_num - cfy2(j) + cell_frac_y(j)
      gyy( 0,j) = 1.5_num - 2.0_num * cfy2(j)
      gyy( 1,j) = 0.25_num + cfy2(j) - cell_frac_y(j) 

  ! Now redo shifted by half a cell due to grid stagger
  ! Use shifted version for ex in X, ey in Y, ez in Z
  ! And in Y&Z for bx, X&Z for by, X&Y for bz

    END DO ! End do-loop with gxx

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
  
      cell_frac_x(j) = REAL(cell_x2(j), num) - cell_x_r(j) + 0.5_num
      cell_frac_y(j) = REAL(cell_y2(j), num) - cell_y_r(j) + 0.5_num

      cell_x2(j) = cell_x2(j) + 1
      cell_y2(j) = cell_y2(j) + 1

    END DO ! End do-loop for re-doing cell_frac_(x,y)

    dcellx = 0
    dcelly = 0

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      cfx2(j) = cell_frac_x(j)**2
      cfy2(j) = cell_frac_y(j)**2

    END DO

  ! Calculating hxx
  ! NOTE: These weights require an additional multiplication factor

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      hxx(dcellx(j)-1,j) = 0.25_num - cfx2(j) + cell_frac_x(j)
      hxx(dcellx(j)  ,j) = 1.5_num - 2.0_num * cfx2(j)
      hxx(dcellx(j)+1,j) = 0.25_num + cfx2(j) - cell_frac_x(j)

      hyy(dcelly(j)-1,j) = 0.25_num + cfy2(j) + cell_frac_y(j)
      hyy(dcelly(j)  ,j) = 1.5_num - 2.0_num * cfy2(j)
      hyy(dcelly(j)+1,j) = 0.25_num + cfy2(j) - cell_frac_y(j)

    END DO ! End do-loop for hxx

  ! These are the electric and magnetic fields interpolated to the
  ! particle position. They have been checked and are correct.
  ! Actually checking this is messy

  ! Calculate e-fields at particle position for triangle particles

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
  
      ex_part(j) = &
               gyy(-1,j) * (hxx(-1,j) * ex(cell_x2(j)-1,cell_y1(j)-1) &
             +              hxx( 0,j) * ex(cell_x2(j)  ,cell_y1(j)-1) &
             +              hxx( 1,j) * ex(cell_x2(j)+1,cell_y1(j)-1)) &
             + gyy( 0,j) * (hxx(-1,j) * ex(cell_x2(j)-1,cell_y1(j)  ) &
             +              hxx( 0,j) * ex(cell_x2(j)  ,cell_y1(j)  ) &
             +              hxx( 1,j) * ex(cell_x2(j)+1,cell_y1(j)  )) &
             + gyy( 1,j) * (hxx(-1,j) * ex(cell_x2(j)-1,cell_y1(j)+1) &
             +              hxx( 0,j) * ex(cell_x2(j)  ,cell_y1(j)+1) &
             +              hxx( 1,j) * ex(cell_x2(j)+1,cell_y1(j)+1))

      ey_part(j) = &
               hyy(-1,j) * (gxx(-1,j) * ey(cell_x1(j)-1,cell_y2(j)-1) &
             +              gxx( 0,j) * ey(cell_x1(j)  ,cell_y2(j)-1) &
             +              gxx( 1,j) * ey(cell_x1(j)+1,cell_y2(j)-1)) &
             + hyy( 0,j) * (gxx(-1,j) * ey(cell_x1(j)-1,cell_y2(j)  ) &
             +              gxx( 0,j) * ey(cell_x1(j)  ,cell_y2(j)  ) &
             +              gxx( 1,j) * ey(cell_x1(j)+1,cell_y2(j)  )) &
             + hyy( 1,j) * (gxx(-1,j) * ey(cell_x1(j)-1,cell_y2(j)+1) &
             +              gxx( 0,j) * ey(cell_x1(j)  ,cell_y2(j)+1) &
             +              gxx( 1,j) * ey(cell_x1(j)+1,cell_y2(j)+1))

      ez_part(j) = &
               gyy(-1,j) * (gxx(-1,j) * ez(cell_x1(j)-1,cell_y1(j)-1) &
             +              gxx( 0,j) * ez(cell_x1(j)  ,cell_y1(j)-1) &
             +              gxx( 1,j) * ez(cell_x1(j)+1,cell_y1(j)-1)) &
             + gyy( 0,j) * (gxx(-1,j) * ez(cell_x1(j)-1,cell_y1(j)  ) &
             +              gxx( 0,j) * ez(cell_x1(j)  ,cell_y1(j)  ) &
             +              gxx( 1,j) * ez(cell_x1(j)+1,cell_y1(j)  )) &
             + gyy( 1,j) * (gxx(-1,j) * ez(cell_x1(j)-1,cell_y1(j)+1) &
             +              gxx( 0,j) * ez(cell_x1(j)  ,cell_y1(j)+1) &
             +              gxx( 1,j) * ez(cell_x1(j)+1,cell_y1(j)+1))

    END DO ! End do-loop for e-fields at particle position

  ! Calculate b-fields at particle position for triangle particles

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
  
      bx_part(j) = &
               hyy(-1,j) * (gxx(-1,j) * bx(cell_x1(j)-1,cell_y2(j)-1) &
             +              gxx( 0,j) * bx(cell_x1(j)  ,cell_y2(j)-1) &
             +              gxx( 1,j) * bx(cell_x1(j)+1,cell_y2(j)-1)) &
             + hyy( 0,j) * (gxx(-1,j) * bx(cell_x1(j)-1,cell_y2(j)  ) &
             +              gxx( 0,j) * bx(cell_x1(j)  ,cell_y2(j)  ) &
             +              gxx( 1,j) * bx(cell_x1(j)+1,cell_y2(j)  )) &
             + hyy( 1,j) * (gxx(-1,j) * bx(cell_x1(j)-1,cell_y2(j)+1) &
             +              gxx( 0,j) * bx(cell_x1(j)  ,cell_y2(j)+1) &
             +              gxx( 1,j) * bx(cell_x1(j)+1,cell_y2(j)+1))

      by_part(j) = &
               gyy(-1,j) * (hxx(-1,j) * by(cell_x2(j)-1,cell_y1(j)-1) &
             +              hxx( 0,j) * by(cell_x2(j)  ,cell_y1(j)-1) &
             +              hxx( 1,j) * by(cell_x2(j)+1,cell_y1(j)-1)) &
             + gyy( 0,j) * (hxx(-1,j) * by(cell_x2(j)-1,cell_y1(j)  ) &
             +              hxx( 0,j) * by(cell_x2(j)  ,cell_y1(j)  ) &
             +              hxx( 1,j) * by(cell_x2(j)+1,cell_y1(j)  )) &
             + gyy( 1,j) * (hxx(-1,j) * by(cell_x2(j)-1,cell_y1(j)+1) &
             +              hxx( 0,j) * by(cell_x2(j)  ,cell_y1(j)+1) &
             +              hxx( 1,j) * by(cell_x2(j)+1,cell_y1(j)+1))

      bz_part(j) = &
               hyy(-1,j) * (hxx(-1,j) * bz(cell_x2(j)-1,cell_y2(j)-1) &
             +              hxx( 0,j) * bz(cell_x2(j)  ,cell_y2(j)-1) &
             +              hxx( 1,j) * bz(cell_x2(j)+1,cell_y2(j)-1)) &
             + gyy( 0,j) * (hxx(-1,j) * bz(cell_x2(j)-1,cell_y2(j)  ) &
             +              hxx( 0,j) * bz(cell_x2(j)  ,cell_y2(j)  ) &
             +              hxx( 1,j) * bz(cell_x2(j)+1,cell_y2(j)  )) &
             + gyy( 1,j) * (hxx(-1,j) * bz(cell_x2(j)-1,cell_y2(j)+1) &
             +              hxx( 0,j) * bz(cell_x2(j)  ,cell_y2(j)+1) &
             +              hxx( 1,j) * bz(cell_x2(j)+1,cell_y2(j)+1))

    END DO ! End do-loop for b-fields at particle position


  ! Update particle momenta using weighted fields

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j -1

      uxm(j) = px(jj) + cmratio(jj) * ex_part(j)
      uym(j) = py(jj) + cmratio(jj) * ey_part(j)
      uzm(j) = pz(jj) + cmratio(jj) * ez_part(j)

    END DO

#ifdef HC_PUSH

  ! Half timestep, then use Higuera-Cary push
  ! See https://aip.scitation.org/doi/10.1063/1.4979989
      
    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
      
      jj = i + j - 1
      gamma_rel(jj) = uxm(j)**2 + uym(j)**2 + uzm(j)**2 + 1.0_num
      alpha(j)      = 0.5_num * q(jj) * dt / m(jj)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      beta_x(j) = alpha(j) * bx_part(j)
      beta_y(j) = alpha(j) * by_part(j)
      beta_z(j) = alpha(j) * bz_part(j)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      beta2(j)      = beta_x(j)**2 + beta_y(j)**2 + beta_z(j)**2
      beta_dot_u(j) = beta_x(j) * uxm(j) + beta_y(j) * uym(j) + beta_z(j) * uzm(j)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j -1
      gamma_rel(jj) = SQRT(0.5_num &
                     * (gamma_rel(jj) - beta2(j) &
                     +  SQRT((gamma_rel(jj) - beta2(j))**2 &
	             + 4.0_num * (beta2(j) + beta_dot_u(j)**2))))  

    END DO

#else

  ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      gamma_rel(jj) = SQRT(uxm(j)**2 + uym(j)**2 + uzm(j)**2)

    END DO

#endif

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
    
      jj = i + j - 1
      root(jj) = ccmratio(jj) / gamma_rel(jj)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      taux(j) = bx_part(j) * root(jj)
      tauy(j) = by_part(j) * root(jj)
      tauz(j) = bz_part(j) * root(jj)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      taux2(j) = taux(j)**2
      tauy2(j) = taux(j)**2
      tauz2(j) = tauz(j)**2

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      tau(j) = 1.0_num / (1.0_num + taux2(j) + tauy2(j) + tauz2(j))

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      uxp(j) = ((1.0_num + taux2(j) - tauy2(j) - tauz2(j)) * uxm(j) &
             + 2.0_num * ((taux(j) * tauy(j) + tauz(j)) * uym(j) &
             + (taux(j) * tauz(j) - tauy(j)) * uzm(j))) * tau(j) 
      uyp(j) = ((1.0_num - taux2(j) + tauy2(j) - tauz2(j)) * uym(j) &
             + 2.0_num * ((tauy(j) * tauz(j) + taux(j)) * uzm(j) &
             + (tauy(j) * taux(j) - tauz(j)) * uxm(j))) * tau(j)
      uzp(j) = ((1.0_num - taux2(j) - tauy2(j) + tauz2(j)) * uzm(j) &
             + 2.0_num * ((tauz(j) * taux(j) + tauy(j)) * uxm(j) &
             + (tauz(j) * tauy(j) - taux(j)) * uym(j))) * tau(j)

    END DO

  ! Rotation over, go to full timestep

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      px(jj) = uxp(j) + cmratio(jj) * ex_part(j)
      py(jj) = uyp(j) + cmratio(jj) * ey_part(j)
      pz(jj) = uzp(j) + cmratio(jj) * ez_part(j)

    END DO

  ! Calculate particle velocity from particle momentum

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      gamma_rel(jj) = SQRT(px(jj)**2 + py(jj)**2 + pz(jj)**2 + 1.0_num)
      igamma(jj)    = 1.0_num / SQRT(px(jj)**2 + py(jj)**2 + pz(jj)**2 + 1.0_num)
      root(jj)      = dtco2 / SQRT(px(jj)**2 + py(jj)**2 + pz(jj)**2 + 1.0_num)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      delta_x(jj) = px(jj) * root(jj)
      delta_y(jj) = py(jj) * root(jj)
      vz(jj)      = pz(jj) * c * igamma(jj) 

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      x(jj) = x(jj) + delta_x(j)
      y(jj) = y(jj) + delta_y(j)

    END DO

  ! Particle has now finished move to end of timestep, so copy back
  ! into particle array

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
 
      jj = i + j - 1
      next => current
      current%part_pos = (/ x(jj) + x_grid_min_local, &
          y(jj) + y_grid_min_local /)
      current%part_p   = mc(jj) * (/ px(jj), py(jj), pz(jj) /)

  ! Add particle to boundary candidate list
      IF (current%part_pos(1) < bnd_x_min &
          .OR. current%part_pos(1) > bnd_x_max &
          .OR. current%part_pos(2) < bnd_y_min &
          .OR. current%part_pos(2) > bnd_y_max ) THEN
        ALLOCATE(bnd_part_next)
        bnd_part_next%particle => current
        bnd_part_last%next => bnd_part_next
        bnd_part_last => bnd_part_next
      END IF

#ifdef WORK_DONE_INTEGRATED
  ! This is the actual total work done by the fields: Results correspond
  ! with the electron's gamma factor

      root(jj) = cmratio(jj) / gamma_rel(jj)

      work_x = ex_part(j) * (tmp_x(j) + px(jj) * root(jj))
      work_y = ex_part(j) * (tmp_y(j) + py(jj) * root(jj))
      work_z = ex_part(j) * (tmp_z(j) + pz(jj) * root(jj))

      current%work_x = work_x
      current%work_y = work_y
      current%work_z = work_z
     
      current%work_x_total = current%work_x_total + work_x
      current%work_y_total = current%work_y_total + work_y
      current%work_z_total = current%work_z_total + work_z
#endif

#ifndef NO_PARTICLE_PROBES
      final_x(jj) = current%part_pos(1)
      final_y(jj) = current%part_pos(2)
#endif

      current => next

    END DO  
      

  END DO ! End do-loop for i

END SUBROUTINE particle_pusher

SUBROUTINE triangle_current_deposition()

  DO i = 1, species_list(ispecies)%attached_list%count, LVEC
    
    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1

  ! Advance to t + 1.5dt to calculate current. This is detailed in
  ! the PSC manual between page 37 and 41. The version coded up looks
  ! completely different to that in the manual, but this is equivalent.
  ! Use t + 1.5dt so that can update J to t + dt at 2nd order

      x(jj) = x(jj) + delta_x(jj)
      y(jj) = y(jj) + delta_y(jj)

  ! Delta-f calculation: subtract background from
  ! calculated current.

#ifdef DELTAF_METHOD
    weight_back(j) = pvol(jj) * f0(ispecies, mc(jj) / c, px(jj), py(jj), pz(jj))
    fcx(j) = idty * (weight(jj) - weight_back(j))
    fcy(j) = idtx * (weight(jj) - weight_back(j))
    fcz(j) = idxy * (weight(jj) - weight_back(j))
#endif

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
      
      jj = i + j - 1
      cell_x_r(j) = x(jj) * idx
      cell_y_r(j) = y(jj) * idy

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      cell_x3(j) = FLOOR(cell_x_r(j) + 0.5_num)
      cell_y3(j) = FLOOR(cell_y_r(j) + 0.5_num)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      cell_frac_x(j) = REAL(cell_x3(j), num) - cell_x_r(j)
      cell_frac_y(j) = REAL(cell_y3(j), num) - cell_y_r(j)

      cell_x3(j) = cell_x3(j) + 1
      cell_y3(j) = cell_y3(j) + 1

    END DO

    hxx = 0.0_num
    hyy = 0.0_num

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      dcellx(j) = cell_x3(j) - cell_x1(j)
      dcelly(j) = cell_y3(j) - cell_y1(j)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      cfx2(j) = cell_frac_x(j)**2
      cfy2(j) = cell_frac_y(j)**2

    END DO

  ! Calculating hxx
  ! NOTE: These weights require an additional multiplication factor

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      hxx(dcellx(j)-1,j) = 0.25_num - cfx2(j) + cell_frac_x(j)
      hxx(dcellx(j)  ,j) = 1.5_num - 2.0_num * cfx2(j)
      hxx(dcellx(j)+1,j) = 0.25_num + cfx2(j) - cell_frac_x(j)

      hyy(dcelly(j)-1,j) = 0.25_num + cfy2(j) + cell_frac_y(j)
      hyy(dcelly(j)  ,j) = 1.5_num - 2.0_num * cfy2(j)
      hyy(dcelly(j)+1,j) = 0.25_num + cfy2(j) - cell_frac_y(j)

    END DO ! End do-loop for hxx

  ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of 
  ! the current update much simpler

    hxx = hxx - gxx
    hyy = hxx - gyy

  ! Remember that due to CFL condition particle can never cross more
  ! than one gridcell in one timestep

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jj = i + j - 1
      xmin(j) = sf_min + (dcellx(j) - 1) / 2
      ymin(j) = sf_min + (dcelly(j) - 1) / 2

      fjx(jj) = fcx(jj) * q(jj)
      fjy(jj) = fcy(jj) * q(jj)
      fjz(jj) = fcz(jj) * q(jj) * vz(jj)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      yfac10(j) = gyy(ymin(j),j) + 0.5_num * hyy(ymin(j),j)
      yfac11(j) = gyy(ymin(j) + 1,j) + 0.5_num * hyy(ymin(j) + 1,j)
      yfac12(j) = gyy(ymin(j) + 2,j) + 0.5_num * hyy(ymin(j) + 2,j)

      yfac20(j) = third * hyy(ymin(j),j) + 0.5_num * gyy(ymin(j),j)
      yfac21(j) = third * hyy(ymin(j) + 1,j) + 0.5_num * gyy(ymin(j) + 1,j)
      yfac22(j) = third * hyy(ymin(j) + 2,j) + 0.5_num * gyy(ymin(j) + 2,j)

      xfac10(j) = gxx(xmin(j),j) + 0.5_num * hxx(ymin(j),j)
      xfac11(j) = gxx(xmin(j) + 1,j) + 0.5_num * hxx(ymin(j) + 1,j)
      xfac12(j) = gxx(xmin(j) + 2,j) + 0.5_num * hxx(ymin(j) + 2,j)

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)
  
      jj = i + j - 1

      wx(1,j) = hxx(xmin(j),j) * yfac10(j)
      wx(2,j) = hxx(xmin(j) + 1,j) * yfac10(j)
      wx(3,j) = hxx(xmin(j) + 2,j) * yfac10(j)
      wx(4,j) = hxx(xmin(j),j) * yfac11(j)
      wx(5,j) = hxx(xmin(j) + 1,j) * yfac11(j)
      wx(6,j) = hxx(xmin(j) + 2,j) * yfac11(j)
      wx(7,j) = hxx(xmin(j),j) * yfac12(j)
      wx(8,j) = hxx(xmin(j) + 1,j) * yfac12(j)
      wx(9,j) = hxx(xmin(j) + 2,j) * yfac12(j)

      wy(1,j) = hyy(ymin(j),j) * xfac10(j)
      wy(2,j) = hyy(ymin(j),j) * xfac11(j)
      wy(3,j) = hyy(ymin(j),j) * xfac12(j)
      wy(4,j) = hyy(ymin(j) + 1,j) * xfac10(j)
      wy(5,j) = hyy(ymin(j) + 1,j) * xfac11(j)
      wy(6,j) = hyy(ymin(j) + 1,j) * xfac12(j)
      wy(7,j) = hyy(ymin(j) + 2,j) * xfac10(j)
      wy(8,j) = hyy(ymin(j) + 2,j) * xfac11(j)
      wy(9,j) = hyy(ymin(n) + 2,j) * xfac12(j)

      wz(1,j) = gxx(xmin(j)) * yfac10(j) + hxx(xmin(j),j) * yfac20(j)
      wz(2,j) = gxx(xmin(j) + 1,j) * yfac10(j) + hxx(xmin(j) + 1,j) * yfac20(j)
      wz(3,j) = gxx(xmin(j) + 2,j) * yfac10(j) + hxx(xmin(j) + 2,j) * yfac20(j)
      wz(4,j) = gxx(xmin(j),j) * yfac11(j) + hxx(xmin(j),j) * yfac21(j)
      wz(5,j) = gxx(xmin(j) + 1,j) * yfac11(j) + hxx(xmin(j) + 1,j) * yfac21(j)
      wz(6,j) = gxx(xmin(j) + 2,j) * yfac11(j) + hxx(xmin(j) + 2,j) * yfac21(j)
      wz(7,j) = gxx(xmin(j),j) * yfac12(j) + hxx(xmin(j),j) * yfac22(j)
      wz(8,j) = gxx(xmin(j) + 1,j) * yfac12(j) + hxx(xmin(j) + 1,j) * yfac22(j)
      wz(9,j) = gxx(xmin(j) + 2,j) * yfac12(j) + hxx(xmin(j) + 2,j) * yfac22(j)

      cx(j) = cell_x1(jj) + xmin(j)
      cy(j) = cell_y1(jj) + ymin(j)
      
    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      ic(j) = cx(j) + (cy(j) - 1) * nx

    END DO

    DO j = 1, MIN(LVEC, species_list(ispecies)%attached_list%count - i + 1)

      jxh(1,ic(j)) = jxh(1,ic(j)) - fjx(j) * wx(1,j)
      jxh(2,ic(j)) = jxh(2,ic(j)) - fjx(j) * (wx(1,j) + wx(2,j))
      jxh(3,ic(j)) = jxh(3,ic(j)) - fjx(j) * (wx(1,j) + wx(2,j) + wx(3,j))
      jxh(4,ic(j)) = jxh(4,ic(j)) - fjx(j) * wx(4,j)
      jxh(5,ic(j)) = jxh(5,ic(j)) - fjx(j) * (wx(4,j) + wx(5,j))
      jxh(6,ic(j)) = jxh(6,ic(j)) - fjx(j) * (wx(4,j) + wx(5,j) + wx(6,j))
      jxh(7,ic(j)) = jxh(7,ic(j)) - fjx(j) * wx(7,j)
      jxh(8,ic(j)) = jxh(8,ic(j)) - fjx(j) * (wx(7,j) + wx(8,j))
      jxh(9,ic(j)) = jxh(9,ic(j)) - fjx(j) * (wx(7,j) + wx(8,j) + wx(9,j))

      jyh(1,ic(j)) = jyh(1,ic(j)) - fjy(j) * wy(1,j)
      jyh(2,ic(j)) = jyh(2,ic(j)) - fjy(j) * wy(2,j)
      jyh(3,ic(j)) = jyh(3,ic(j)) - fjy(j) * wy(3,j)
      jyh(4,ic(j)) = jyh(4,ic(j)) - fjy(j) * (wy(1,j) + wy(4,j))
      jyh(5,ic(j)) = jyh(5,ic(j)) - fjy(j) * (wy(2,j) + wy(5,j))
      jyh(6,ic(j)) = jyh(6,ic(j)) - fjy(j) * (wy(3,j) + wy(6,j))
      jyh(7,ic(j)) = jyh(7,ic(j)) - fjy(j) * (wy(1,j) + wy(4,j) + wy(7,j))
      jyh(8,ic(j)) = jyh(8,ic(j)) - fjy(j) * (wy(2,j) + wy(5,j) + wy(8,j))
      jyh(9,ic(j)) = jyh(9,ic(j)) - fjy(j) * (wy(3,j) + wy(6,j) + wy(9,j))

      jzh(1,ic(j)) = jzh(1,ic(j)) + fjz(j) * wz(1,j)
      jzh(2,ic(j)) = jzh(2,ic(j)) + fjz(j) * wz(2,j)
      jzh(3,ic(j)) = jzh(3,ic(j)) + fjz(j) * wz(3,j)
      jzh(4,ic(j)) = jzh(4,ic(j)) + fjz(j) * wz(4,j)
      jzh(5,ic(j)) = jzh(5,ic(j)) + fjz(j) * wz(5,j)
      jzh(6,ic(j)) = jzh(6,ic(j)) + fjz(j) * wz(6,j)
      jzh(7,ic(j)) = jzh(7,ic(j)) + fjz(j) * wz(7,j)
      jzh(8,ic(j)) = jzh(8,ic(j)) + fjz(j) * wz(8,j)
      jzh(9,ic(j)) = jzh(9,ic(j)) + fjz(j) * wz(9,j)

    END DO

  END DO ! End do-loop for index i

  ! Deposit current on the cells

  DO j = 1, ny
    DO i = 1, nx
      iic = (j - 1) * nx + i
      
      jx(i,j)         = jx(i,j) + jxh(1,iic)
      jx(i + 1,j)     = jx(i + 1,j) + jxh(2,iic)
      jx(i + 2,j)     = jx(i + 2,j) + jxh(3,iic)
      jx(i,j + 1)     = jx(i, j + 1) + jxh(4,iic)
      jx(i + 1,j + 1) = jx(i + 1, j + 1) + jxh(5,iic)
      jx(i + 2,j + 1) = jx(i + 2, j + 1) + jxh(6,iic)
      jx(i,j + 2)     = jx(i,j + 2) + jxh(7,iic)
      jx(i + 1,j + 2) = jx(i + 1,j + 2) + jxh(8,iic)
      jx(i + 2,j + 2) = jx(i + 2,j + 2) + jxh(9,iic)

      jy(i,j)         = jy(i,j) + jyh(1,iic)
      jy(i + 1,j)     = jy(i + 1,j) + jyh(2,iic)
      jy(i + 2,j)     = jy(i + 2,j) + jyh(3,iic)
      jy(i,j + 1)     = jy(i, j + 1) + jyh(4,iic)
      jy(i + 1,j + 1) = jy(i + 1, j + 1) + jyh(5,iic)
      jy(i + 2,j + 1) = jy(i + 2, j + 1) + jyh(6,iic)
      jy(i,j + 2)     = jy(i,j + 2) + jyh(7,iic)
      jy(i + 1,j + 2) = jy(i + 1,j + 2) + jyh(8,iic)
      jy(i + 2,j + 2) = jy(i + 2,j + 2) + jyh(9,iic)

      jz(i,j)         = jz(i,j) + jzh(1,iic)
      jz(i + 1,j)     = jz(i + 1,j) + jzh(2,iic)
      jz(i + 2,j)     = jz(i + 2,j) + jzh(3,iic)
      jz(i,j + 1)     = jz(i, j + 1) + jzh(4,iic)
      jz(i + 1,j + 1) = jz(i + 1, j + 1) + jzh(5,iic)
      jz(i + 2,j + 1) = jz(i + 2, j + 1) + jzh(6,iic)
      jz(i,j + 2)     = jz(i,j + 2) + jzh(7,iic)
      jz(i + 1,j + 2) = jz(i + 1,j + 2) + jzh(8,iic)
      jz(i + 2,j + 2) = jz(i + 2,j + 2) + jzh(9,iic)

    END DO
  END DO

END SUBROUTINE triangle_current_deposition

  FUNCTION f0(ispecies, mass, px, py, pz)

    INTEGER, INTENT(IN) :: ispecies
    REAL(num), INTENT(IN) :: mass
    REAL(num), INTENT(IN) :: px, py, pz
    REAL(num) :: f0
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz, density
    REAL(num) :: f0_exponent, norm, two_kb_mass, two_pi_kb_mass3
    TYPE(particle_species), POINTER :: species

    species => species_list(ispecies)

    IF (ABS(species%initial_conditions%density_back) > c_tiny) THEN
       two_kb_mass = 2.0_num * kb * mass
       two_pi_kb_mass3 = (pi * two_kb_mass)**3

       Tx = species%initial_conditions%temp_back(1)
       Ty = species%initial_conditions%temp_back(2)
       Tz = species%initial_conditions%temp_back(3)
       driftx  = species%initial_conditions%drift_back(1)
       drifty  = species%initial_conditions%drift_back(2)
       driftz  = species%initial_conditions%drift_back(3)
       density = species%initial_conditions%density_back
       f0_exponent = ((px - driftx)**2 / Tx &
                    + (py - drifty)**2 / Ty &
                    + (pz - driftz)**2 / Tz) / two_kb_mass
       norm = density / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
       f0 = norm * EXP(-f0_exponent)
    ELSE
       f0 = 0.0_num
    END IF

  END FUNCTION f0














