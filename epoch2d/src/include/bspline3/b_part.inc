        bx_part = &
              hy(-2) * (gx(-2) * bx(cell_x1-2,cell_y2-2) &
            +           gx(-1) * bx(cell_x1-1,cell_y2-2) &
            +           gx( 0) * bx(cell_x1  ,cell_y2-2) &
            +           gx( 1) * bx(cell_x1+1,cell_y2-2) &
            +           gx( 2) * bx(cell_x1+2,cell_y2-2)) &
            + hy(-1) * (gx(-2) * bx(cell_x1-2,cell_y2-1) &
            +           gx(-1) * bx(cell_x1-1,cell_y2-1) &
            +           gx( 0) * bx(cell_x1  ,cell_y2-1) &
            +           gx( 1) * bx(cell_x1+1,cell_y2-1) &
            +           gx( 2) * bx(cell_x1+2,cell_y2-1)) &
            + hy( 0) * (gx(-2) * bx(cell_x1-2,cell_y2  ) &
            +           gx(-1) * bx(cell_x1-1,cell_y2  ) &
            +           gx( 0) * bx(cell_x1  ,cell_y2  ) &
            +           gx( 1) * bx(cell_x1+1,cell_y2  ) &
            +           gx( 2) * bx(cell_x1+2,cell_y2  )) &
            + hy( 1) * (gx(-2) * bx(cell_x1-2,cell_y2+1) &
            +           gx(-1) * bx(cell_x1-1,cell_y2+1) &
            +           gx( 0) * bx(cell_x1  ,cell_y2+1) &
            +           gx( 1) * bx(cell_x1+1,cell_y2+1) &
            +           gx( 2) * bx(cell_x1+2,cell_y2+1)) &
            + hy( 2) * (gx(-2) * bx(cell_x1-2,cell_y2+2) &
            +           gx(-1) * bx(cell_x1-1,cell_y2+2) &
            +           gx( 0) * bx(cell_x1  ,cell_y2+2) &
            +           gx( 1) * bx(cell_x1+1,cell_y2+2) &
            +           gx( 2) * bx(cell_x1+2,cell_y2+2))

        by_part = &
              gy(-2) * (hx(-2) * by(cell_x2-2,cell_y1-2) &
            +           hx(-1) * by(cell_x2-1,cell_y1-2) &
            +           hx( 0) * by(cell_x2  ,cell_y1-2) &
            +           hx( 1) * by(cell_x2+1,cell_y1-2) &
            +           hx( 2) * by(cell_x2+2,cell_y1-2)) &
            + gy(-1) * (hx(-2) * by(cell_x2-2,cell_y1-1) &
            +           hx(-1) * by(cell_x2-1,cell_y1-1) &
            +           hx( 0) * by(cell_x2  ,cell_y1-1) &
            +           hx( 1) * by(cell_x2+1,cell_y1-1) &
            +           hx( 2) * by(cell_x2+2,cell_y1-1)) &
            + gy( 0) * (hx(-2) * by(cell_x2-2,cell_y1  ) &
            +           hx(-1) * by(cell_x2-1,cell_y1  ) &
            +           hx( 0) * by(cell_x2  ,cell_y1  ) &
            +           hx( 1) * by(cell_x2+1,cell_y1  ) &
            +           hx( 2) * by(cell_x2+2,cell_y1  )) &
            + gy( 1) * (hx(-2) * by(cell_x2-2,cell_y1+1) &
            +           hx(-1) * by(cell_x2-1,cell_y1+1) &
            +           hx( 0) * by(cell_x2  ,cell_y1+1) &
            +           hx( 1) * by(cell_x2+1,cell_y1+1) &
            +           hx( 2) * by(cell_x2+2,cell_y1+1)) &
            + gy( 2) * (hx(-2) * by(cell_x2-2,cell_y1+2) &
            +           hx(-1) * by(cell_x2-1,cell_y1+2) &
            +           hx( 0) * by(cell_x2  ,cell_y1+2) &
            +           hx( 1) * by(cell_x2+1,cell_y1+2) &
            +           hx( 2) * by(cell_x2+2,cell_y1+2))

        bz_part = &
              hy(-2) * (hx(-2) * bz(cell_x2-2,cell_y2-2) &
            +           hx(-1) * bz(cell_x2-1,cell_y2-2) &
            +           hx( 0) * bz(cell_x2  ,cell_y2-2) &
            +           hx( 1) * bz(cell_x2+1,cell_y2-2) &
            +           hx( 2) * bz(cell_x2+2,cell_y2-2)) &
            + hy(-1) * (hx(-2) * bz(cell_x2-2,cell_y2-1) &
            +           hx(-1) * bz(cell_x2-1,cell_y2-1) &
            +           hx( 0) * bz(cell_x2  ,cell_y2-1) &
            +           hx( 1) * bz(cell_x2+1,cell_y2-1) &
            +           hx( 2) * bz(cell_x2+2,cell_y2-1)) &
            + hy( 0) * (hx(-2) * bz(cell_x2-2,cell_y2  ) &
            +           hx(-1) * bz(cell_x2-1,cell_y2  ) &
            +           hx( 0) * bz(cell_x2  ,cell_y2  ) &
            +           hx( 1) * bz(cell_x2+1,cell_y2  ) &
            +           hx( 2) * bz(cell_x2+2,cell_y2  )) &
            + hy( 1) * (hx(-2) * bz(cell_x2-2,cell_y2+1) &
            +           hx(-1) * bz(cell_x2-1,cell_y2+1) &
            +           hx( 0) * bz(cell_x2  ,cell_y2+1) &
            +           hx( 1) * bz(cell_x2+1,cell_y2+1) &
            +           hx( 2) * bz(cell_x2+2,cell_y2+1)) &
            + hy( 2) * (hx(-2) * bz(cell_x2-2,cell_y2+2) &
            +           hx(-1) * bz(cell_x2-1,cell_y2+2) &
            +           hx( 0) * bz(cell_x2  ,cell_y2+2) &
            +           hx( 1) * bz(cell_x2+1,cell_y2+2) &
            +           hx( 2) * bz(cell_x2+2,cell_y2+2))
