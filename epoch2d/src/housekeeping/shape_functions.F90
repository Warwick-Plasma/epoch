MODULE shape_functions

USE shared_data
IMPLICIT NONE
CONTAINS


SUBROUTINE particle_to_grid(cell_frac,output)
	REAL(num),INTENT(IN) :: cell_frac
	REAL(num),DIMENSION(-2:2),INTENT(INOUT) :: output
	REAL(num),DIMENSION(-2:2) :: cfs
	INTEGER :: ielement
	
	DO ielement=-2,2
		cfs(ielement)=(cell_frac+REAL(ielement,num))
	ENDDO

#ifdef SPLINE_FOUR
	output(-2) = 1.0_num/384.0_num * (5.0_num + 2.0_num * cfs(-2))**4
	output(-1) = 1.0_num/96.0_num * (55.0_num - 4.0_num * cfs(-1)*(&
		5.0_num + 2.0_num * cfs(-1) * (15.0_num + 2.0_num * cfs(-1) *(&
		cfs(-1) + 5.0_num))))
	output(+0) = 115.0_num/192.0_num - 5.0_num/8.0_num * cfs(0)**2+&
		cfs(0)**4/4.0_num
	output(+1) = 1.0_num/96.0_num * (55.0_num + 4.0_num * cfs(1)*(&
		5.0_num - 2.0_num * cfs(1) * (15.0_num + 2.0_num * cfs(1) *(&
		cfs(1)-5.0_num))))	
	output(+2) = 1.0_num/384.0_num * (5.0_num - 2.0_num * cfs(2))**4
#else	
	output(-2) = 0.0_num
	output(-1) = 0.5_num * (1.5_num - ABS(cfs(-1)))**2
	output(+0) = 0.75_num - ABS(cfs(0))**2
	output(+1) = 0.5_num * (1.5_num - ABS(cfs(1)))**2
	output(-2) = 0.0_num
#endif

END SUBROUTINE particle_to_grid

SUBROUTINE grid_to_particle(cell_frac,output)
	REAL(num),INTENT(IN) :: cell_frac
	REAL(num),DIMENSION(-2:2),INTENT(INOUT) :: output
	REAL(num),DIMENSION(-2:2) :: cfs
	INTEGER :: ielement
	
	DO ielement=-2,2
		cfs(ielement)=(cell_frac+REAL(ielement,num))
	ENDDO

#ifdef SPLINE_FOUR	
	output(-2) = 1.0_num/384.0_num * (5.0_num + 2.0_num * cfs(-2))**4
	output(-1) = 1.0_num/96.0_num * (55.0_num - 4.0_num * cfs(-1)*(&
		5.0_num + 2.0_num * cfs(-1) * (15.0_num + 2.0_num * cfs(-1) *(&
		cfs(-1) + 5.0_num))))
	output(+0) = 115.0_num/192.0_num - 5.0_num/8.0_num * cfs(0)**2+&
		cfs(0)**4/4.0_num
	output(+1) = 1.0_num/96.0_num * (55.0_num + 4.0_num * cfs(1)*(&
		5.0_num - 2.0_num * cfs(1) * (15.0_num + 2.0_num * cfs(1) *(&
		cfs(1)-5.0_num))))	
	output(+2) = 1.0_num/384.0_num * (5.0_num - 2.0_num * cfs(2))**4
#else	
	output(-2) = 0.0_num
	output(-1) = 0.5_num * (0.5_num + cell_frac)**2
   output(+0) = 0.75_num - cell_frac**2
   output(+1) = 0.5_num * (0.5_num - cell_frac)**2
	output(+2) = 0.0_num
#endif
	

END SUBROUTINE grid_to_particle

END MODULE shape_functions
