

  SUBROUTINE get_particle(part_data, n, part)
    TYPE(particle_data), POINTER, INTENT(IN) :: part_data
    INTEGER(i8), INTENT(IN) :: n
    TYPE(particle), INTENT(INOUT) :: part

    IF (n >= 1 .AND. n <= part_data%count) THEN
      part%pos = part_data%pos(n)
      part%mom = part_data%mom(n)
      part%mass = part_data%mass(n)
      part%weight = part_data%weight(n)
      part%charge = part_data%charge(n)
    ENDIF

  END SUBROUTINE get_particle

  SUBROUTINE particle_sort

    ! This routine sorts the particles such that
    ! the memory access to the particle list
    ! is contigious

  END SUBROUTINE particle sort


  SUBROUTINE current_deposition_VB_triangle

#ifdef INTEL_VECTORISATION
!dir$ attributes align:64 :: gx
!dir$ attributes align:64 :: gy
!dir$ attributes align:64 :: gz  
#endif

#ifdef INTEL_VECTORISATION
!dir$ attributes align:64 :: hx
!dir$ attributes align:64 :: hy
!dir$ attributes align:64 :: hz
#endif

#ifdef IBM_VECTORISATION
!IBM* ALIGN(64, gx, gy, gz)
#endif

#ifdef IBM_VECTORISATION
!IBM* ALIGN(64, hx, hy, hz)
#endif

! Similarly all other arrays need to be aligned

    DO ispec = 1, n_species
      species => species_list(ispec)
      
      IF (species%immobile) CYCLE

      part_data => species%part_data
      npart = part_data%count

      CALL particle_sort ! I guess this routine should be called in the Boris pusher 
                         ! such that both gx and hx are in sync

      DO np = 1, npart, LVEC

      !$OMP SIMD
      DO n  = 1, MIN(LVEC, npart - np + 1)

       CALL get_particle(part_data, n, part)

#ifndef PER_SPECIES_WEIGHT
       part_weight(n) = part(n)%weight
       fcx(n) = idty * part_weight(n)
       fcy(n) = idtx * part_weight(n)
       fcz(n) = idxy * part_weight(n)
#endif
#ifndef NO_PARTICLE_PROBES
       init_part_x(n) = part(n)%part_pos(1)
       init_part_y(n) = part(n)%part_pos(2)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
       part_q(n)   = part(n)%charge
       part_m(n)   = part(n)%mass
       part_mc(n)  = c * part(n)%mass
       ipart_mc(n) = 1.0_num / part_mc(n)
       cmratio(n)  = part_q(n) * dtfac * ipart_mc(n)
       ccmratio(n) = c * cmratio(n)
#ifndef NO_PARTICLE_PROBES
       part_mc2 = c * part_mc
#endif
#endif

      !Copy the particle properties out for speed
      part_x(n)  = part(n)%part_pos(1) - x_grid_min_local
      part_y(n)  = part(n)%part_pos(2) - y_grid_min_local
      part_ux(n) = part(n)%part_p(1) * ipart_mc(n)
      part_uy(n) = part(n)%part_p(2) * ipart_mc(n)
      part_uz(n) = part(n)%part_p(3) * ipart_mc(n)

     !Now advance to t+1.5dt to calculate current
     !For efficient vectorization, I would prefer
     !part_x(n)  = part(n)%part_pos(1) - x_grid_min_local + delta_x
     !part_y(n)  = part(n)%part_pos(2) - x_grid_min_local + delta_x
     !This eliminates the dependency on previous step
      part_x(n)  = part_x(n) + delta_x(n)
      part_y(n)  = part_y(n) + delta_y(n)

     !Delta-f calculation: subtract background from calculated current
#ifdef DELTAF_METHOD
      weight_back(n) = part(n)%pvol * f0(ispecies, part_mc(n) / c, &
                       part(n)%part_p)
      fcx(n) = idty * (part_weight(n) - weight_back(n))
      fcy(n) = idtx * (part_weight(n) - weight_back(n))
      fcz(n) = idxy * (part_weight(n) - weight_back(n))
#endif

      cell_x_r(n) = part_x(n) * idx
      cell_y_r(n) = part_y(n) * idy

      cell_x3(n)  = FLOOR(cell_x_r(n) + 0.5_num)
      cell_y3(n)  = FLOOR(cell_y_r(n) + 0.5_num)

      cell_frac_x(n) = REAL(cell_x3(n), num) - cell_x_r(n)
      cell_frac_y(n) = REAL(cell_y3(n), num) - cell_y_r(n)

      fjx(n) = fcx(n) * part_q(n)
      fjy(n) = fcy(n) * part_q(n)
      fjz(n) = fcz(n) * part_q(n) * part_vz(n)

      hx = 0.0_num
      hy = 0.0_num

      dcellx(n) = cell_x3(n) - cell_x1(n)
      dcelly(n) = cell_y3(n) - cell_y1(n)

      xmin(n) = sf_min + (dcellx(n) - 1) / 2
      ymin(n) = sf_max + (dcelly(n) - 1) / 2


      hx(xmin(n))     = 0.25_num + cell_frac_x(n)**2 + cell_frac_x(n)
      hx(xmin(n) + 1) = 1.5_num - 2.0_num * cell_frac_x(n)**2
      hx(xmin(n) + 2) = 0.25_num + cell_frac_x(n)**2 - cell_frac_x(n)

      hy(ymin(n))     = 0.25_num + cell_frac_y(n)**2 + cell_frac_y(n)
      hy(ymin(n) + 1) = 1.5_num - 2.0_num * cell_frac_y(n)**2
      hy(ymin(n) + 2) = 0.25_num + cell_frac_y(n)**2 - cell_frac_y(n)

      yfac10(n) = gy(ymin(n)) + 0.5_num * hy(ymin(n))
      yfac11(n) = gy(ymin(n) + 1) + 0.5_num * hy(ymin(n) + 1)
      yfac12(n) = gy(ymin(n) + 1) + 0.5_num * hy(ymin(n) + 2)

      yfac20(n) = third * hy(ymin(n)) + 0.5 * gy(ymin(n))
      yfac21(n) = third * hy(ymin(n) + 1) + 0.5 * gy(ymin(n) + 1)
      yfac22(n) = third * hy(ymin(n) + 2) + 0.5 * gy(ymin(n) + 2)

      xfac10(n) = gx(xmin(n)) + 0.5_num * hx(xmin(n))
      xfac11(n) = gx(xmin(n) + 1) + 0.5_num * hx(xmin(n) + 1)
      xfac12(n) = gx(xmin(n) + 2) + 0.5_num * hx(xmin(n) + 2)

      wx(n,1) = hx(xmin(n)) * yfac10(n)
      wx(n,2) = hx(xmin(n) + 1) * yfac10(n)
      wx(n,3) = hx(xmin(n) + 2) * yfac10(n)
      wx(n,4) = hx(xmin(n)) * yfac11(n)
      wx(n,5) = hx(xmin(n) + 1) * yfac11(n)
      wx(n,6) = hx(xmin(n) + 2) * yfac11(n)
      wx(n,7) = hx(xmin(n)) * yfac12(n)
      wx(n,8) = hx(xmin(n) + 1) * yfac12(n)
      wx(n,9) = hx(xmin(n) + 2) * yfac12(n)

      wy(n,1) = hy(ymin(n)) * xfac10(n)
      wy(n,2) = hy(ymin(n)) * xfac11(n)
      wy(n,3) = hy(ymin(n)) * xfac12(n)
      wy(n,4) = hy(ymin(n) + 1) * xfac10(n)
      wy(n,5) = hy(ymin(n) + 1) * xfac11(n)
      wy(n,6) = hy(ymin(n) + 1) * xfac12(n)
      wy(n,7) = hy(ymin(n) + 2) * xfac10(n)
      wy(n,8) = hy(ymin(n) + 2) * xfac11(n)
      wy(n,9) = hy(ymin(n) + 2) * xfac12(n)

      wz(n,1) = gx(xmin(n)) * yfac10(n) + hx(xmin(n)) * yfac20(n)
      wz(n,2) = gx(xmin(n) + 1) * yfac10(n) + hx(xmin(n) + 1) * yfac20(n)
      wz(n,3) = gx(xmin(n) + 2) * yfac10(n) + hx(xmin(n) + 2) * yfac20(n)
      wz(n,4) = gx(xmin(n)) * yfac11(n) + hx(xmin(n)) * yfac21(n)
      wz(n,5) = gx(xmin(n) + 1) * yfac11(n) + hx(xmin(n) + 1) * yfac21(n)
      wz(n,6) = gx(xmin(n) + 2) * yfac11(n) + hx(xmin(n) + 2) * yfac21(n)
      wz(n,7) = gx(xmin(n)) * yfac12(n) + hx(xmin(n)) * yfac22(n)
      wz(n,8) = gx(xmin(n) + 1) * yfac12(n) + hx(xmin(n) + 1) * yfac22(n)
      wz(n,9) = gx(xmin(n) + 2) * yfac12(n) + hx(xmin(n) + 2) * yfac22(n)

      cx(n)   = cell_x1(n) + xmin(n)
      cy(n)   = cell_y1(n) + ymin(n)
      cell(n) = cx(n) + (cy(n) - 1) * nx

      jxh(n,cell(n))     = -fjx(n) * wx(n,1)
      jxh(n,cell(n) + 1) = -fjx(n) * (wx(n,1) + wx(n,2))
      jxh(n,cell(n) + 2) = -fjx(n) * (wx(n,1) + wx(n,2) + wx(n,3))
      jxh(n,cell(n) + 3) = -fjx(n) * wx(n,4)
      jxh(n,cell(n) + 4) = -fjx(n) * (wx(n,4) + wx(n,5))
      jxh(n,cell(n) + 5) = -fjx(n) * (wx(n,4) + wx(n,5) + wx(n,6))
      jxh(n,cell(n) + 6) = -fjx(n) * wx(n,7)
      jxh(n,cell(n) + 7) = -fjx(n) * (wx(n,7) + wx(n,8))
      jxh(n,cell(n) + 8) = -fjx(n) * (wx(n,7) + wx(n,8) + wx(n,9))

      jyh(n,cell(n))     = -fjy(n) * wy(n,1)
      jyh(n,cell(n) + 1) = -fjy(n) * wy(n,2)
      jyh(n,cell(n) + 2) = -fjy(n) * wy(n,3)
      jyh(n,cell(n) + 3) = -fjy(n) * (wy(n,1) + wy(n,4))
      jyh(n,cell(n) + 4) = -fjy(n) * (wy(n,2) + wy(n,5))
      jyh(n,cell(n) + 5) = -fjy(n) * (wy(n,3) + wy(n,6))
      jyh(n,cell(n) + 6) = -fjy(n) * (wy(n,1) + wy(n,4) + wy(n,7))
      jyh(n,cell(n) + 7) = -fjy(n) * (wy(n,2) + wy(n,5) + wy(n,8))
      jyh(n,cell(n) + 8) = -fjy(n) * (wy(n,3) + wy(n,6) + wy(n,9))

      jzh(n,cell(n))     = fjz(n) * wz(n,1)
      jzh(n,cell(n) + 1) = fjz(n) * wz(n,2)
      jzh(n,cell(n) + 2) = fjz(n) * wz(n,3)
      jzh(n,cell(n) + 3) = fjz(n) * wz(n,4)
      jzh(n,cell(n) + 4) = fjz(n) * wz(n,5)
      jzh(n,cell(n) + 5) = fjz(n) * wz(n,6)
      jzh(n,cell(n) + 6) = fjz(n) * wz(n,7)
      jzh(n,cell(n) + 7) = fjz(n) * wz(n,8)
      jzh(n,cell(n) + 8) = fjz(n) * wz(n,9)

    END DO !END LOOP  n = 1, MIN(LVEC, npart - np + 1)

    !$OMP END SIMD

    !$OMP SIMD

    DO n = 1, MIN(LVEC, npart - np + 1)

        jx(cx(n),cy(n))          = jx(cx(n),cy(n)) + jxh(n,1)
        jx(cx(n) + 1,cy(n))      = jx(cx(n) + 1,cy(n)) + jxh(n,2)
        jx(cx(n) + 2,cy(n))      = jx(cx(n) + 2,cy(n)) + jxh(n,3)
        jx(cx(n), cy(n) + 1)     = jx(cx(n),cy(n) + 1) + jxh(n,4)
        jx(cx(n) + 1, cy(n) + 1) = jx(cx(n) + 1,cy(n) + 1) + jxh(n,5)
        jx(cx(n) + 2, cy(n) + 1) = jx(cx(n) + 2,cy(n) + 1) + jxh(n,6)
        jx(cx(n), cy(n) + 2)     = jx(cx(n),cy(n) + 2) + jxh(n,7)
        jx(cx(n) + 1, cy(n) + 2) = jx(cx(n),cy(n) + 2) + jxh(n,8)
        jx(cx(n) + 2, cy(n) + 2) = jx(cx(n) + 2,cy(n) + 2) + jxh(n,9)

        jy(cx(n),cy(n))          = jy(cx(n),cy(n)) + jyh(n,1)
        jy(cx(n) + 1,cy(n))      = jy(cx(n) + 1,cy(n)) + jyh(n,2)
        jy(cx(n) + 2,cy(n))      = jy(cx(n) + 2,cy(n)) + jyh(n,3)
        jy(cx(n), cy(n) + 1)     = jy(cx(n),cy(n) + 1) + jyh(n,4)
        jy(cx(n) + 1, cy(n) + 1) = jy(cx(n) + 1,cy(n) + 1) + jyh(n,5)
        jy(cx(n) + 2, cy(n) + 1) = jy(cx(n) + 2,cy(n) + 1) + jyh(n,6)
        jy(cx(n), cy(n) + 2)     = jy(cx(n),cy(n) + 2) + jyh(n,7)
        jy(cx(n) + 1, cy(n) + 2) = jy(cx(n),cy(n) + 2) + jyh(n,8)
        jy(cx(n) + 2, cy(n) + 2) = jy(cx(n) + 2,cy(n) + 2) + jyh(n,9)

        jz(cx(n),cy(n))          = jz(cx(n),cy(n)) + jzh(n,1)
        jz(cx(n) + 1,cy(n))      = jz(cx(n) + 1,cy(n)) + jzh(n,2)
        jz(cx(n) + 2,cy(n))      = jz(cx(n) + 2,cy(n)) + jzh(n,3)
        jz(cx(n), cy(n) + 1)     = jz(cx(n),cy(n) + 1) + jzh(n,4)
        jz(cx(n) + 1, cy(n) + 1) = jz(cx(n) + 1,cy(n) + 1) + jzh(n,5)
        jz(cx(n) + 2, cy(n) + 1) = jz(cx(n) + 2,cy(n) + 1) + jzh(n,6)
        jz(cx(n), cy(n) + 2)     = jz(cx(n),cy(n) + 2) + jzh(n,7)
        jz(cx(n) + 1, cy(n) + 2) = jz(cx(n),cy(n) + 2) + jzh(n,8)
        jz(cx(n) + 2, cy(n) + 2) = jz(cx(n) + 2,cy(n) + 2) + jzh(n,9)

    END DO !END LOOP n = 1, MIN(LVEC, npart - np + 1)
    !$OMP END SIMD
         
    END DO !END LOOP  n = 1, npart
    END DO !END LOOP ispec = 1, nspecies

  END current_deposition_VB_triangle


















