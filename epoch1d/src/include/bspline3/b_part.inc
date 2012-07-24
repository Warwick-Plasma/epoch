        bx_part = &
              gx(-2) * bx(cell_x1-2) &
            + gx(-1) * bx(cell_x1-1) &
            + gx( 0) * bx(cell_x1  ) &
            + gx( 1) * bx(cell_x1+1) &
            + gx( 2) * bx(cell_x1+2)

        by_part = &
              hx(-2) * by(cell_x2-2) &
            + hx(-1) * by(cell_x2-1) &
            + hx( 0) * by(cell_x2  ) &
            + hx( 1) * by(cell_x2+1) &
            + hx( 2) * by(cell_x2+2)

        bz_part = &
              hx(-2) * bz(cell_x2-2) &
            + hx(-1) * bz(cell_x2-1) &
            + hx( 0) * bz(cell_x2  ) &
            + hx( 1) * bz(cell_x2+1) &
            + hx( 2) * bz(cell_x2+2)
