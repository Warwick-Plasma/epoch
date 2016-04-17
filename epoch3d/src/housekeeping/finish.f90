! Copyright (C) 2014-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
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

MODULE finish

  USE shared_data
  USE diagnostics
  USE setup
  USE shunt
  USE partlist
  USE deck
  USE window
  USE laser
  USE collisions
  USE dist_fn
  USE ionise

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finalise

CONTAINS

  SUBROUTINE finalise

    CALL close_files
    IF (done_mpi_initialise) CALL deallocate_memory
    CALL MPI_FINALIZE(errcode)
    STOP

  END SUBROUTINE finalise



  SUBROUTINE deallocate_memory

    INTEGER :: i, n, stat

    DEALLOCATE(x, xb, x_global, xb_global, xb_offset_global)
    DEALLOCATE(y, xb, y_global, yb_global, yb_offset_global)
    DEALLOCATE(z, xb, z_global, zb_global, zb_offset_global)
    DEALLOCATE(ex, ey, ez, bx, by, bz, jx, jy, jz)

    DEALLOCATE(npart_each_rank)
    DEALLOCATE(x_grid_mins, x_grid_maxs, cell_x_min, cell_x_max)
    DEALLOCATE(y_grid_mins, y_grid_maxs, cell_y_min, cell_y_max)
    DEALLOCATE(z_grid_mins, z_grid_maxs, cell_z_min, cell_z_max)

    DEALLOCATE(ex_x_min, ex_x_max, ey_x_min, ey_x_max, ez_x_min, ez_x_max)
    DEALLOCATE(bx_x_min, bx_x_max, by_x_min, by_x_max, bz_x_min, bz_x_max)
    DEALLOCATE(ex_y_min, ex_y_max, ey_y_min, ey_y_max, ez_y_min, ez_y_max)
    DEALLOCATE(bx_y_min, bx_y_max, by_y_min, by_y_max, bz_y_min, bz_y_max)
    DEALLOCATE(ex_z_min, ex_z_max, ey_z_min, ey_z_max, ez_z_min, ez_z_max)
    DEALLOCATE(bx_z_min, bx_z_max, by_z_min, by_z_max, bz_z_min, bz_z_max)

    DO i = 1, n_species
      CALL deallocate_stack(species_list(i)%density_function)
      DO n = 1, 3
        CALL deallocate_stack(species_list(i)%temperature_function(n))
        CALL deallocate_stack(species_list(i)%drift_function(n))
      ENDDO
      CALL destroy_partlist(species_list(i)%attached_list)
      DEALLOCATE(species_list(i)%ext_temp_x_min, STAT=stat)
      DEALLOCATE(species_list(i)%ext_temp_x_max, STAT=stat)
      DEALLOCATE(species_list(i)%ext_temp_y_min, STAT=stat)
      DEALLOCATE(species_list(i)%ext_temp_y_max, STAT=stat)
      DEALLOCATE(species_list(i)%ext_temp_z_min, STAT=stat)
      DEALLOCATE(species_list(i)%ext_temp_z_max, STAT=stat)
    ENDDO

    DEALLOCATE(species_list, STAT=stat)

    DEALLOCATE(io_block_list, STAT=stat)
    DEALLOCATE(io_list_data, STAT=stat)
    DEALLOCATE(initial_conditions, STAT=stat)
    DEALLOCATE(file_prefixes, STAT=stat)
    DEALLOCATE(file_numbers, STAT=stat)

    DO i = 1, n_subsets
      DEALLOCATE(subset_list(i)%dumpmask, STAT=stat)
      DEALLOCATE(subset_list(i)%use_species, STAT=stat)
    ENDDO
    DEALLOCATE(subset_list, STAT=stat)

    DO i = 1, n_deck_constants
      CALL deallocate_stack(deck_constant_list(i)%execution_stream)
    ENDDO
    DEALLOCATE(deck_constant_list, STAT=stat)

    CALL deallocate_input_deck_buffer
    CALL deallocate_window
    CALL deallocate_lasers
    CALL deallocate_collisions
    CALL deallocate_file_list
    CALL deallocate_dist_fns
    CALL deallocate_ionisation
    CALL deallocate_partlists
    CALL deallocate_eval_stack

    CALL MPI_COMM_FREE(comm, errcode)

  END SUBROUTINE deallocate_memory

END MODULE finish
