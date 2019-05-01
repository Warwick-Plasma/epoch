        part_var = &
              gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2) &
            +           gx(-1) * grid_var(cell_x1-1,cell_y1-2) &
            +           gx( 0) * grid_var(cell_x1  ,cell_y1-2) &
            +           gx( 1) * grid_var(cell_x1+1,cell_y1-2) &
            +           gx( 2) * grid_var(cell_x1+2,cell_y1-2)) &
            + gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1) &
            +           gx(-1) * grid_var(cell_x1-1,cell_y1-1) &
            +           gx( 0) * grid_var(cell_x1  ,cell_y1-1) &
            +           gx( 1) * grid_var(cell_x1+1,cell_y1-1) &
            +           gx( 2) * grid_var(cell_x1+2,cell_y1-1)) &
            + gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ) &
            +           gx(-1) * grid_var(cell_x1-1,cell_y1  ) &
            +           gx( 0) * grid_var(cell_x1  ,cell_y1  ) &
            +           gx( 1) * grid_var(cell_x1+1,cell_y1  ) &
            +           gx( 2) * grid_var(cell_x1+2,cell_y1  )) &
            + gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1) &
            +           gx(-1) * grid_var(cell_x1-1,cell_y1+1) &
            +           gx( 0) * grid_var(cell_x1  ,cell_y1+1) &
            +           gx( 1) * grid_var(cell_x1+1,cell_y1+1) &
            +           gx( 2) * grid_var(cell_x1+2,cell_y1+1)) &
            + gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2) &
            +           gx(-1) * grid_var(cell_x1-1,cell_y1+2) &
            +           gx( 0) * grid_var(cell_x1  ,cell_y1+2) &
            +           gx( 1) * grid_var(cell_x1+1,cell_y1+2) &
            +           gx( 2) * grid_var(cell_x1+2,cell_y1+2))
