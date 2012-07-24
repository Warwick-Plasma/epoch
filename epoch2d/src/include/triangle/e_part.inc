        ex_part = &
              gy(-1) * (hx(-1) * ex(cell_x2-1,cell_y1-1) &
            +           hx( 0) * ex(cell_x2  ,cell_y1-1) &
            +           hx( 1) * ex(cell_x2+1,cell_y1-1)) &
            + gy( 0) * (hx(-1) * ex(cell_x2-1,cell_y1  ) &
            +           hx( 0) * ex(cell_x2  ,cell_y1  ) &
            +           hx( 1) * ex(cell_x2+1,cell_y1  )) &
            + gy( 1) * (hx(-1) * ex(cell_x2-1,cell_y1+1) &
            +           hx( 0) * ex(cell_x2  ,cell_y1+1) &
            +           hx( 1) * ex(cell_x2+1,cell_y1+1))

        ey_part = &
              hy(-1) * (gx(-1) * ey(cell_x1-1,cell_y2-1) &
            +           gx( 0) * ey(cell_x1  ,cell_y2-1) &
            +           gx( 1) * ey(cell_x1+1,cell_y2-1)) &
            + hy( 0) * (gx(-1) * ey(cell_x1-1,cell_y2  ) &
            +           gx( 0) * ey(cell_x1  ,cell_y2  ) &
            +           gx( 1) * ey(cell_x1+1,cell_y2  )) &
            + hy( 1) * (gx(-1) * ey(cell_x1-1,cell_y2+1) &
            +           gx( 0) * ey(cell_x1  ,cell_y2+1) &
            +           gx( 1) * ey(cell_x1+1,cell_y2+1))

        ez_part = &
              gy(-1) * (gx(-1) * ez(cell_x1-1,cell_y1-1) &
            +           gx( 0) * ez(cell_x1  ,cell_y1-1) &
            +           gx( 1) * ez(cell_x1+1,cell_y1-1)) &
            + gy( 0) * (gx(-1) * ez(cell_x1-1,cell_y1  ) &
            +           gx( 0) * ez(cell_x1  ,cell_y1  ) &
            +           gx( 1) * ez(cell_x1+1,cell_y1  )) &
            + gy( 1) * (gx(-1) * ez(cell_x1-1,cell_y1+1) &
            +           gx( 0) * ez(cell_x1  ,cell_y1+1) &
            +           gx( 1) * ez(cell_x1+1,cell_y1+1))
