#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = (parameters%pack_pos(1) - x_grid_min_local) / dx - 0.5_num
        cell_y_r = (parameters%pack_pos(2) - y_grid_min_local) / dy - 0.5_num
#else
        cell_x_r = (parameters%pack_pos(1) - x_grid_min_local) / dx
        cell_y_r = (parameters%pack_pos(2) - y_grid_min_local) / dy
#endif
        cell_x = FLOOR(cell_x_r + 0.5_num)
        cell_y = FLOOR(cell_y_r + 0.5_num)
        cell_frac_x = REAL(cell_x, num) - cell_x_r
        cell_frac_y = REAL(cell_y, num) - cell_y_r
        cell_x = cell_x + 1
        cell_y = cell_y + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gxfac.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gxfac.inc"
#else
#include "triangle/gxfac.inc"
#endif
