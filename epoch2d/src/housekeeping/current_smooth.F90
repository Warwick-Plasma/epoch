MODULE current_smooth

USE shared_data
USE shape_functions
IMPLICIT NONE

CONTAINS

SUBROUTINE smooth_current()
! A very simple current smoothing routine
	
	CALL smooth_array(jx)
	CALL smooth_array(jy)
	CALL smooth_array(jz)
	
END SUBROUTINE smooth_current

SUBROUTINE smooth_array(array)
	REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: array
	REAL(num),DIMENSION(:,:),ALLOCATABLE :: wk_array
	INTEGER :: ix, iy
		
	ALLOCATE(wk_array(1:nx,1:ny))
#ifdef HIGH_ORDER_SMOOTHING
	CALL particle_to_grid(0.0_num,weight_fn)
	
	wk_array=0.0_num

	DO ix=1,nx
		DO iy=1,ny
			DO icyclex=-sf_order,sf_order
				DO icycley=-sf_order,sf_order
					wk_array(ix,iy) = wk_array(ix,iy) + array(ix+icyclex,iy+icycley)&
					 	* weight_fn(icyclex) * weight_fn(icycley)
				ENDDO
			ENDDO
		ENDDO
	ENDDO
#else
	DO ix=1,nx
		DO iy=1,ny
			wk_array(ix,iy) = (array(ix,iy) + array(ix-1,iy)*0.25_num + &
				array(ix+1,iy)*0.25_num + array(ix,iy-1) *0.25_num + &
				array(ix,iy+1)*0.25_num)*0.5_num
		ENDDO
	ENDDO
#endif
	array(1:nx,1:ny)=wk_array
	
	DEALLOCATE(wk_array)
	
END SUBROUTINE smooth_array

END MODULE current_smooth
