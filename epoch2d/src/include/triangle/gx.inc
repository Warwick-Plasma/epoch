        cf2 = cell_frac_x**2
        gx(-1) = 0.25_num + cf2 + cell_frac_x
        gx( 0) = 1.5_num - 2.0_num * cf2
        gx( 1) = 0.25_num + cf2 - cell_frac_x

        cf2 = cell_frac_y**2
        gy(-1) = 0.25_num + cf2 + cell_frac_y
        gy( 0) = 1.5_num - 2.0_num * cf2
        gy( 1) = 0.25_num + cf2 - cell_frac_y
