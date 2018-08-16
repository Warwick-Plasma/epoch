! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
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

MODULE particles

  USE current_smooth
  USE boundary
  USE partlist
#ifdef PREFETCH
  USE prefetch
#endif

  IMPLICIT NONE

CONTAINS

  SUBROUTINE push_particles

    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_x3
    INTEGER :: cell_y1, cell_y2, cell_y3

    ! Xi (space factor see page 38 in manual)
    ! The code now uses gx and hx instead of xi0 and xi1

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
    REAL(num) :: jxh, jzh
    REAL(num), DIMENSION(sf0-1:sf1+1) :: jyh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, part_mc, ipart_mc, part_weight, part_m
#ifdef HC_PUSH
    REAL(num) :: beta_x, beta_y, beta_z, beta2, beta_dot_u, alpha, sigma
#endif

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy, part_mc2
#endif

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy

    ! Fields at particle location
    REAL(num) :: ex_part, ey_part, ez_part, bx_part, by_part, bz_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: uxp, uxm, uyp, uym, uzp, uzm
    REAL(num) :: tau, taux, tauy, tauz, taux2, tauy2, tauz2

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio, ccmratio

    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax
    REAL(num) :: wx, wy, wz

    ! Temporary variables
    REAL(num) :: idx, idy
    REAL(num) :: idty, idtx, idxy
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: fcx, fcy, fcz, fjx, fjy, fjz
    REAL(num) :: root, dtfac, gamma_rel, part_u2, third, igamma
    REAL(num) :: delta_x, delta_y, part_vz
    REAL(num) :: hy_iy, xfac1, yfac1, yfac2
    INTEGER :: ispecies, ix, iy, dcellx, dcelly, cx, cy
    INTEGER(i8) :: ipart
#ifdef WORK_DONE_INTEGRATED
    REAL(num) :: tmp_x, tmp_y, tmp_z
    REAL(num) :: work_x, work_y, work_z
#endif
#ifndef NO_PARTICLE_PROBES
    LOGICAL :: probes_for_species
    REAL(num) :: gamma_rel_m1
#endif
#ifndef NO_TRACER_PARTICLES
    LOGICAL :: not_tracer_species
#endif
    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif
#ifdef DELTAF_METHOD
    REAL(num) :: weight_back
#endif

    TYPE(particle), POINTER :: current, next

#ifdef PREFETCH
    CALL prefetch_particle(species_list(1)%attached_list%head)
#endif

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    gx = 0.0_num
    gy = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    dtfac = 0.5_num * dt * fac
    third = 1.0_num / 3.0_num

    idty = idt * idy * fac
    idtx = idt * idx * fac
    idxy = idx * idy * fac

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (species_list(ispecies)%immobile) CYCLE
      IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
#ifdef PHOTONS
        IF (photon_dynamics) CALL push_photons(ispecies)
#endif
        CYCLE
      END IF
#ifndef NO_PARTICLE_PROBES
      current_probe => species_list(ispecies)%attached_probes
      probes_for_species = ASSOCIATED(current_probe)
#endif
#ifndef NO_TRACER_PARTICLES
      not_tracer_species = .NOT. species_list(ispecies)%tracer
#endif

#ifdef PER_SPECIES_WEIGHT
      part_weight = species_list(ispecies)%weight
      fcx = idty * part_weight
      fcy = idtx * part_weight
      fcz = idxy * part_weight
#endif
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q   = species_list(ispecies)%charge
      part_m   = species_list(ispecies)%mass
      part_mc  = c * species_list(ispecies)%mass
      ipart_mc = 1.0_num / part_mc
      cmratio  = part_q * dtfac * ipart_mc
      ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
      part_mc2 = c * part_mc
#endif
#endif
      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next
#ifdef PREFETCH
        CALL prefetch_particle(next)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_weight = current%weight
        fcx = idty * part_weight
        fcy = idtx * part_weight
        fcz = idxy * part_weight
#endif
#ifndef NO_PARTICLE_PROBES
        init_part_x = current%part_pos(1)
        init_part_y = current%part_pos(2)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q   = current%charge
        part_m   = current%mass
        part_mc  = c * current%mass
        ipart_mc = 1.0_num / part_mc
        cmratio  = part_q * dtfac * ipart_mc
        ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
        part_mc2 = c * part_mc
#endif
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1) - x_grid_min_local
        part_y  = current%part_pos(2) - y_grid_min_local
        part_ux = current%part_p(1) * ipart_mc
        part_uy = current%part_p(2) * ipart_mc
        part_uz = current%part_p(3) * ipart_mc

        ! Calculate v(t) from p(t)
        ! See PSC manual page (25-27)
        gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        root = dtco2 / gamma_rel

        ! Move particles to half timestep position to first order
        part_x = part_x + part_ux * root
        part_y = part_y + part_uy * root

#ifdef WORK_DONE_INTEGRATED
        ! This is the actual total work done by the fields: Results correspond
        ! with the electron's gamma factor
        root = cmratio / gamma_rel

        tmp_x = part_ux * root
        tmp_y = part_uy * root
        tmp_z = part_uz * root
#endif

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x * idx - 0.5_num
        cell_y_r = part_y * idy - 0.5_num
#else
        cell_x_r = part_x * idx
        cell_y_r = part_y * idy
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        dcellx = 0
        dcelly = 0
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#include "tophat/b_part.inc"
#else
#include "triangle/e_part.inc"
#include "triangle/b_part.inc"
#endif

        ! update particle momenta using weighted fields
        uxm = part_ux + cmratio * ex_part
        uym = part_uy + cmratio * ey_part
        uzm = part_uz + cmratio * ez_part

#ifdef HC_PUSH
        ! Half timestep, then use Higuera-Cary push see
        ! https://aip.scitation.org/doi/10.1063/1.4979989
        gamma_rel = uxm**2 + uym**2 + uzm**2 + 1.0_num
        alpha = 0.5_num * part_q * dt / part_m
        beta_x = alpha * bx_part
        beta_y = alpha * by_part
        beta_z = alpha * bz_part
        beta2 = beta_x**2 + beta_y**2 + beta_z**2
        sigma = gamma_rel - beta2
        beta_dot_u = beta_x * uxm + beta_y * uym + beta_z * uzm
        gamma_rel = sigma + SQRT(sigma**2 + 4.0_num * (beta2 + beta_dot_u**2))
        gamma_rel = SQRT(0.5_num * gamma_rel)
#else
        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
        gamma_rel = SQRT(uxm**2 + uym**2 + uzm**2 + 1.0_num)
#endif
        root = ccmratio / gamma_rel

        taux = bx_part * root
        tauy = by_part * root
        tauz = bz_part * root

        taux2 = taux**2
        tauy2 = tauy**2
        tauz2 = tauz**2

        tau = 1.0_num / (1.0_num + taux2 + tauy2 + tauz2)

        uxp = ((1.0_num + taux2 - tauy2 - tauz2) * uxm &
            + 2.0_num * ((taux * tauy + tauz) * uym &
            + (taux * tauz - tauy) * uzm)) * tau
        uyp = ((1.0_num - taux2 + tauy2 - tauz2) * uym &
            + 2.0_num * ((tauy * tauz + taux) * uzm &
            + (tauy * taux - tauz) * uxm)) * tau
        uzp = ((1.0_num - taux2 - tauy2 + tauz2) * uzm &
            + 2.0_num * ((tauz * taux + tauy) * uxm &
            + (tauz * tauy - taux) * uym)) * tau

        ! Rotation over, go to full timestep
        part_ux = uxp + cmratio * ex_part
        part_uy = uyp + cmratio * ey_part
        part_uz = uzp + cmratio * ez_part

        ! Calculate particle velocity from particle momentum
        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        gamma_rel = SQRT(part_u2 + 1.0_num)
        igamma = 1.0_num / gamma_rel
        root = dtco2 * igamma

        delta_x = part_ux * root
        delta_y = part_uy * root
        part_vz = part_uz * c * igamma

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + delta_x
        part_y = part_y + delta_y

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = (/ part_x + x_grid_min_local, &
            part_y + y_grid_min_local /)
        current%part_p   = part_mc * (/ part_ux, part_uy, part_uz /)

#ifdef WORK_DONE_INTEGRATED
        ! This is the actual total work done by the fields: Results correspond
        ! with the electron's gamma factor
        root = cmratio / gamma_rel

        work_x = ex_part * (tmp_x + part_ux * root)
        work_y = ey_part * (tmp_y + part_uy * root)
        work_z = ez_part * (tmp_z + part_uz * root)

        current%work_x = work_x
        current%work_y = work_y
        current%work_z = work_z
        current%work_x_total = current%work_x_total + work_x
        current%work_y_total = current%work_y_total + work_y
        current%work_z_total = current%work_z_total + work_z
#endif

#ifndef NO_PARTICLE_PROBES
        final_part_x = current%part_pos(1)
        final_part_y = current%part_pos(2)
#endif
        ! Original code calculates densities of electrons, ions and neutrals
        ! here. This has been removed to reduce memory footprint

        ! If the code is compiled with tracer particle support then put in an
        ! IF statement so that the current is not calculated for this species
#ifndef NO_TRACER_PARTICLES
        IF (not_tracer_species) THEN
#endif
          ! Now advance to t+1.5dt to calculate current. This is detailed in
          ! the manual between pages 37 and 41. The version coded up looks
          ! completely different to that in the manual, but is equivalent.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x = part_x + delta_x
          part_y = part_y + delta_y

          ! Delta-f calcuation: subtract background from
          ! calculated current.
#ifdef DELTAF_METHOD
          weight_back = current%pvol * f0(ispecies, part_mc / c, current%part_p)
          fcx = idty * (part_weight - weight_back)
          fcy = idtx * (part_weight - weight_back)
          fcz = idxy * (part_weight - weight_back)
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
          cell_x_r = part_x * idx - 0.5_num
          cell_y_r = part_y * idy - 0.5_num
#else
          cell_x_r = part_x * idx
          cell_y_r = part_y * idy
#endif
          cell_x3 = FLOOR(cell_x_r + 0.5_num)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r
          cell_x3 = cell_x3 + 1

          cell_y3 = FLOOR(cell_y_r + 0.5_num)
          cell_frac_y = REAL(cell_y3, num) - cell_y_r
          cell_y3 = cell_y3 + 1

          hx = 0.0_num
          hy = 0.0_num

          dcellx = cell_x3 - cell_x1
          dcelly = cell_y3 - cell_y1
          ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

          ! Now change Xi1* to be Xi1*-Xi0*. This makes the representation of
          ! the current update much simpler
          hx = hx - gx
          hy = hy - gy

          ! Remember that due to CFL condition particle can never cross more
          ! than one gridcell in one timestep

          xmin = sf_min + (dcellx - 1) / 2
          xmax = sf_max + (dcellx + 1) / 2

          ymin = sf_min + (dcelly - 1) / 2
          ymax = sf_max + (dcelly + 1) / 2

          fjx = fcx * part_q
          fjy = fcy * part_q
          fjz = fcz * part_q * part_vz

          jyh = 0.0_num
          DO iy = ymin, ymax
            cy = cell_y1 + iy
            yfac1 =         gy(iy) + 0.5_num * hy(iy)
            yfac2 = third * hy(iy) + 0.5_num * gy(iy)

            hy_iy = hy(iy)

            jxh = 0.0_num
            DO ix = xmin, xmax
              cx = cell_x1 + ix
              xfac1 = gx(ix) + 0.5_num * hx(ix)

              wx = hx(ix) * yfac1
              wy = hy_iy  * xfac1
              wz = gx(ix) * yfac1 + hx(ix) * yfac2

              ! This is the bit that actually solves d(rho)/dt = -div(J)
              jxh = jxh - fjx * wx
              jyh(ix) = jyh(ix) - fjy * wy
              jzh = fjz * wz

              jx(cx, cy) = jx(cx, cy) + jxh
              jy(cx, cy) = jy(cx, cy) + jyh(ix)
              jz(cx, cy) = jz(cx, cy) + jzh
            END DO
          END DO
#ifndef NO_TRACER_PARTICLES
        END IF
#endif
#ifndef NO_PARTICLE_PROBES
        IF (probes_for_species) THEN
          ! Compare the current particle with the parameters of any probes in
          ! the system. These particles are copied into a separate part of the
          ! output file.

          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          current_probe => species_list(ispecies)%attached_probes

          ! Cycle through probes
          DO WHILE(ASSOCIATED(current_probe))
            ! Note that this is the energy of a single REAL particle in the
            ! pseudoparticle, NOT the energy of the pseudoparticle
            probe_energy = gamma_rel_m1 * part_mc2

            ! Unidirectional probe
            IF (probe_energy > current_probe%ek_min) THEN
              IF (probe_energy < current_probe%ek_max) THEN

                d_init  = SUM(current_probe%normal &
                    * (current_probe%point - (/init_part_x, init_part_y/)))
                d_final = SUM(current_probe%normal &
                    * (current_probe%point - (/final_part_x, final_part_y/)))
                IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                  ! this particle is wanted so copy it to the list associated
                  ! with this probe
                  ALLOCATE(particle_copy)
                  particle_copy = current
                  CALL add_particle_to_partlist(&
                      current_probe%sampled_particles, particle_copy)
                  NULLIFY(particle_copy)
                END IF

              END IF
            END IF
            current_probe => current_probe%next
          END DO
        END IF
#endif
        current => next
      END DO
      CALL current_bcs(species=ispecies)
    END DO

    CALL current_bcs
    CALL particle_bcs

    IF (smooth_currents) CALL smooth_current()

    IF (use_current_correction) THEN
      jx = jx - initial_jx
      jy = jy - initial_jy
      jz = jz - initial_jz
    END IF

  END SUBROUTINE push_particles



  ! Background distribution function used for delta-f calculations.
  ! Specialise to a drifting (tri)-Maxwellian to simplify and ensure
  ! zero density/current divergence.
  ! Can effectively switch off deltaf method by setting zero background density.

  FUNCTION f0(ispecies, mass, p)

    INTEGER, INTENT(IN) :: ispecies
    REAL(num), INTENT(IN) :: mass
    REAL(num), DIMENSION(:), INTENT(IN) :: p
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
       f0_exponent = ((p(1) - driftx)**2 / Tx &
                    + (p(2) - drifty)**2 / Ty &
                    + (p(3) - driftz)**2 / Tz) / two_kb_mass
       norm = density / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
       f0 = norm * EXP(-f0_exponent)
    ELSE
       f0 = 0.0_num
    END IF

  END FUNCTION f0



#ifdef PHOTONS
  SUBROUTINE push_photons(ispecies)

    ! Very simple photon pusher
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: delta_x, delta_y
    INTEGER,INTENT(IN) :: ispecies
    TYPE(particle), POINTER :: current

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy, dtfac, fac
    LOGICAL :: probes_for_species
#endif

#ifndef NO_PARTICLE_PROBES
    current_probe => species_list(ispecies)%attached_probes
    probes_for_species = ASSOCIATED(current_probe)
#endif
    dtfac = dt * c**2

    ! set current to point to head of list
    current => species_list(ispecies)%attached_list%head
    ! loop over photons
    DO WHILE(ASSOCIATED(current))
      ! Note that this is the energy of a single REAL particle in the
      ! pseudoparticle, NOT the energy of the pseudoparticle
      probe_energy = current%particle_energy

      fac = dtfac / probe_energy
      delta_x = current%part_p(1) * fac
      delta_y = current%part_p(2) * fac
#ifndef NO_PARTICLE_PROBES
      init_part_x = current%part_pos(1)
      init_part_y = current%part_pos(2)
#endif
      current%part_pos = current%part_pos + (/delta_x, delta_y/)
#ifndef NO_PARTICLE_PROBES
      final_part_x = current%part_pos(1)
      final_part_y = current%part_pos(2)
#endif

#ifndef NO_PARTICLE_PROBES
      IF (probes_for_species) THEN
        ! Compare the current particle with the parameters of any probes in
        ! the system. These particles are copied into a separate part of the
        ! output file.

        current_probe => species_list(ispecies)%attached_probes

        ! Cycle through probes
        DO WHILE(ASSOCIATED(current_probe))
          ! Unidirectional probe
          IF (probe_energy > current_probe%ek_min) THEN
            IF (probe_energy < current_probe%ek_max) THEN

              d_init  = SUM(current_probe%normal &
                  * (current_probe%point - (/init_part_x, init_part_y/)))
              d_final = SUM(current_probe%normal &
                  * (current_probe%point - (/final_part_x, final_part_y/)))
              IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                ! this particle is wanted so copy it to the list associated
                ! with this probe
                ALLOCATE(particle_copy)
                particle_copy = current
                CALL add_particle_to_partlist(&
                    current_probe%sampled_particles, particle_copy)
                NULLIFY(particle_copy)
              END IF

            END IF
          END IF
          current_probe => current_probe%next
        END DO
      END IF
#endif

      current => current%next
    END DO

  END SUBROUTINE push_photons
#endif

END MODULE particles
