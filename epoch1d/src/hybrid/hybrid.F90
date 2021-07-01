! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
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
!
!------------------------------------------------------------------------------
!
! hybrid.F90
!
! This is the main interface to the hybrid routines. When running in hybrid
! mode, we ignore the main PIC loop and control is diverted to the hybrid PIC
! loop. This is also the script which initialises the hybrid arrays, and checks
! solids have been correctly defined.

MODULE hybrid

#ifdef HYBRID
  USE balance
  USE boundary
  USE current_smooth
  USE diagnostics
  USE injectors
  USE particles
  USE particle_migration
  USE random_generator
#ifdef PHOTONS
  USE photons
#endif
#ifdef BREMSSTRAHLUNG
  USE bremsstrahlung
#endif
#endif

  IMPLICIT NONE

CONTAINS

#ifdef HYBRID

  SUBROUTINE run_hybrid_PIC(push, halt, force_dump)

    LOGICAL :: push, halt, force_dump

  END SUBROUTINE run_hybrid_PIC



  SUBROUTINE deallocate_hybrid

  END SUBROUTINE deallocate_hybrid



#else
  ! If hybrid precompiler option is off, let the routines called by other
  ! modules do nothing
  SUBROUTINE run_hybrid_PIC(push, halt, force_dump)

    LOGICAL :: push, halt, force_dump

  END SUBROUTINE run_hybrid_PIC



  SUBROUTINE deallocate_hybrid

  END SUBROUTINE deallocate_hybrid

#endif

END MODULE hybrid
