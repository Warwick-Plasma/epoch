        wt_var1 = wt_facx * (ABS(wt_dtx + cell_frac_x) &
            - ABS(wt_dtx - cell_frac_x))
        gx( 0) = 0.5_num + wt_var1
        gx( 1) = 0.5_num - wt_var1
