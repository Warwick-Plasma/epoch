        cx2 = cell_frac_x**2
        gx(-2) = fac1 * (0.5_num + cell_frac_x)**4
        gx(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_x &
            + cx2 * (1.5_num - cell_frac_x - cx2))
        gx( 0) = 0.25_num * (fac3 + cx2 * (cx2 - 2.5_num))
        gx( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_x &
            + cx2 * (1.5_num + cell_frac_x - cx2))
        gx( 2) = fac1 * (0.5_num - cell_frac_x)**4

        cy2 = cell_frac_y**2
        gy(-2) = fac1 * (0.5_num + cell_frac_y)**4
        gy(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_y &
            + cy2 * (1.5_num - cell_frac_y - cy2))
        gy( 0) = 0.25_num * (fac3 + cy2 * (cy2 - 2.5_num))
        gy( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_y &
            + cy2 * (1.5_num + cell_frac_y - cy2))
        gy( 2) = fac1 * (0.5_num - cell_frac_y)**4

        cz2 = cell_frac_z**2
        gz(-2) = fac1 * (0.5_num + cell_frac_z)**4
        gz(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_z &
            + cz2 * (1.5_num - cell_frac_z - cz2))
        gz( 0) = 0.25_num * (fac3 + cz2 * (cz2 - 2.5_num))
        gz( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_z &
            + cz2 * (1.5_num + cell_frac_z - cz2))
        gz( 2) = fac1 * (0.5_num - cell_frac_z)**4
