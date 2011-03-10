        cx2 = cell_frac_x**2
        gx(-1) = 0.5_num * (0.25_num + cx2 + cell_frac_x)
        gx( 0) = 0.75_num - cx2
        gx( 1) = 0.5_num * (0.25_num + cx2 - cell_frac_x)

        cy2 = cell_frac_y**2
        gy(-1) = 0.5_num * (0.25_num + cy2 + cell_frac_y)
        gy( 0) = 0.75_num - cy2
        gy( 1) = 0.5_num * (0.25_num + cy2 - cell_frac_y)

        cz2 = cell_frac_z**2
        gz(-1) = 0.5_num * (0.25_num + cz2 + cell_frac_z)
        gz( 0) = 0.75_num - cz2
        gz( 1) = 0.5_num * (0.25_num + cz2 - cell_frac_z)
