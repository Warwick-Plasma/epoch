        ex_part = &
              hx(-2) * ex(cell_x2-2) &
            + hx(-1) * ex(cell_x2-1) &
            + hx( 0) * ex(cell_x2  ) &
            + hx( 1) * ex(cell_x2+1) &
            + hx( 2) * ex(cell_x2+2)

        ey_part = &
              gx(-2) * ey(cell_x1-2) &
            + gx(-1) * ey(cell_x1-1) &
            + gx( 0) * ey(cell_x1  ) &
            + gx( 1) * ey(cell_x1+1) &
            + gx( 2) * ey(cell_x1+2)

        ez_part = &
              gx(-2) * ez(cell_x1-2) &
            + gx(-1) * ez(cell_x1-1) &
            + gx( 0) * ez(cell_x1  ) &
            + gx( 1) * ez(cell_x1+1) &
            + gx( 2) * ez(cell_x1+2)
