MODULE Current_Smooth

USE shared_data
USE shape_functions
IMPLICIT NONE

CONTAINS

SUBROUTINE Smooth_Current()
! A very simple current smoothing routine
	
	CALL Smooth_Array(jx)
	CALL Smooth_Array(jy)
	CALL Smooth_Array(jz)
	
END SUBROUTINE Smooth_Current

SUBROUTINE Smooth_Array(Array)
	REAL(num),DIMENSION(-2:,-2:),INTENT(INOUT) :: Array
	REAL(num),DIMENSION(:,:),ALLOCATABLE :: WkArray
	INTEGER :: ix, iy
		
	ALLOCATE(WkArray(1:nx,1:ny))
#ifdef HIGH_ORDER_SMOOTHING
	CALL ParticleToGrid(0.0_num,Weight_Fn)
	
	WkArray=0.0_num

	DO ix=1,nx
		DO iy=1,ny
			DO icyclex=-sf_order,sf_order
				DO icycley=-sf_order,sf_order
					WkArray(ix,iy) = WkArray(ix,iy) + Array(ix+icyclex,iy+icycley)&
					 	* Weight_Fn(icyclex) * Weight_Fn(icycley)
				ENDDO
			ENDDO
		ENDDO
	ENDDO
#else
	DO ix=1,nx
		DO iy=1,ny
			WkArray(ix,iy) = (Array(ix,iy) + Array(ix-1,iy)*0.25_num + &
				Array(ix+1,iy)*0.25_num + Array(ix,iy-1) *0.25_num + &
				Array(ix,iy+1)*0.25_num)*0.5_num
		ENDDO
	ENDDO
#endif
	Array(1:nx,1:ny)=WkArray
	
	DEALLOCATE(WkArray)
	
END SUBROUTINE Smooth_Array

END MODULE Current_Smooth
