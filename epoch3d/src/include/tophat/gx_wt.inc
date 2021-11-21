        wt_var1 = wt_facx * (ABS(wt_dtx + cell_frac_x) &
            - ABS(wt_dtx - cell_frac_x))
        gx( 0) = 0.5_num + wt_var1
        gx( 1) = 0.5_num - wt_var1

        wt_var1 = wt_facy * (ABS(wt_dty + cell_frac_y) &
            - ABS(wt_dty - cell_frac_y))
        gy( 0) = 0.5_num + wt_var1
        gy( 1) = 0.5_num - wt_var1

        wt_var1 = wt_facz * (ABS(wt_dtz + cell_frac_z) &
            - ABS(wt_dtz - cell_frac_z))
        gz( 0) = 0.5_num + wt_var1
        gz( 1) = 0.5_num - wt_var1
