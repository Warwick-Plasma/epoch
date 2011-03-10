        cx2 = cell_frac_x**2
        gx(-2) = fac1 * (0.5_num + cell_frac_x)**4
        gx(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_x &
            + cx2 * (1.5_num - cell_frac_x - cx2))
        gx( 0) = 0.25_num * (fac3 + cx2 * (cx2 - 2.5_num))
        gx( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_x &
            + cx2 * (1.5_num + cell_frac_x - cx2))
        gx( 2) = fac1 * (0.5_num - cell_frac_x)**4
