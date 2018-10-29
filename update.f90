!====================================================================
module m_UpDate

use m_DataStructures
implicit none
contains

subroutine update_robust(fibers, hinges, dt)
implicit none
type(fiber) , dimension(:)            :: fibers
type (rod), dimension(:)              :: hinges
real(8)                               :: dt, deltaX, deltaY, delta
integer(8)                            :: i,j
real(8), dimension(3)                 :: dir_vec


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, DIR_VEC)
!$OMP DO 
do i=1, ubound (fibers,1)
	do j=fibers(i)%first_hinge, (fibers(i)%first_hinge+fibers(i)%nbr_hinges-2)
		dir_vec=hinges(j+1)%X_i-hinges(j)%X_i
		hinges(j)%length2=sqrt(dot_product(dir_vec, dir_vec))
	end do
end do
!$OMP END DO

!$OMP DO
do i=1, ubound (fibers,1)
	hinges(fibers(i)%first_hinge)%X_i=hinges(fibers(i)%first_hinge)%X_i+hinges(fibers(i)%first_hinge)%v_i*dt
	do j=fibers(i)%first_hinge+1, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
		hinges(j)%X_i=hinges(j)%X_i+hinges(j)%v_i*dt
		dir_vec=hinges(j)%X_i-hinges(j-1)%X_i
		dir_vec=dir_vec/(sqrt(dot_product(dir_vec, dir_vec)))
		hinges(j)%X_i=hinges(j-1)%X_i+dir_vec*hinges(j-1)%length2
	end do
end do
!$OMP END DO 
!$OMP END PARALLEL


end subroutine update_robust


!*************************************************
subroutine update(hinges, dt)
type (rod), dimension(:)              :: hinges
real(8)                               :: dt
integer(8)                            :: i


do i=1, ubound(hinges,1)
	if (hinges(i)%is_stationary==0) then
		hinges(i)%X_i=hinges(i)%X_i+hinges(i)%v_i*dt
	end if
end do


end subroutine update
!*************************************************
           
subroutine update_periodic(fibers, hinges, dt, periodic_boundary, box_size, gamma_dot, t)
type(fiber) , dimension(:)            :: fibers
type (rod), dimension(:)              :: hinges
real(8)                               :: dt, cx, cy, cz, dX, t, gamma_dot 
integer(8)                            :: i,j,k,l
real(8), dimension(3)                 :: box_size, X_fiber
logical                               :: periodic_boundary

do i=1, ubound(hinges,1)
	if (hinges(i)%is_stationary==0) then
		if (hinges(i)%is_stationary==0) hinges(i)%X_i=hinges(i)%X_i+hinges(i)%v_i*dt
	end if
end do

dX = mod(gamma_dot * box_size(2)/2*t,box_size(1)) 

if (periodic_boundary.eqv. .true.)then
	do i=1, ubound (fibers,1)
		X_fiber=0
		k=0
        ! Compute center of mass for fiber i
		do l=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
			X_fiber=X_fiber+((hinges(l)%X_i+hinges(l+1)%X_i)/2d0)
			k=k+1
		end do
		X_fiber=X_fiber/real(k)
        
        ! Apply boundary conditions
        
        ! if the center of mass of the fiber is outside of the box, move it accorddingly. 
        ! based loosely on Non-Newtonian molecular Dynamics by Evans and Morriss, appendix B
        cy = int( X_fiber(2) / (box_size(2)/2d0) )
  
        X_fiber(1)=X_fiber(1)-dX*cy
        
        cx = int(X_fiber(1) / (box_size(1)/2d0))
        
        cz = int(X_fiber(3) / (box_size(3)/2d0))
        
        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
            if (hinges(j)%is_stationary==0) then     !2018/07/20   Error  hinges(i) change to hinges(j)
                 hinges(j)%X_i(1)=hinges(j)%X_i(1)-dX*cy-cx*box_size(1)
                 hinges(j)%X_i(2)=hinges(j)%X_i(2)-cy*box_size(2)
                 hinges(j)%X_i(3)=hinges(j)%X_i(3)-cz*box_size(3)
             end if
        end do
        
        
        ! 
        
  !      ! if the coordinates of the center of mass are larger or smaller than the limits of the periodic box, move the fibers accordingly. 
		!do l=1,3
  !          
		!	if ((X_fiber(l)).gt. (box_size(l)/2d0)) then 
		!		do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
		!			if (hinges(i)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)-box_size(l)
  !              end do
  !              
  !          else if (X_fiber(l).lt. (-box_size(l)/2d0)) then 
  !              
		!		do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
		!			if (hinges(i)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)+box_size(l)
  !              end do
  !              
  !          end if
  !          
  !      end do
        
        
	end do
end if


end subroutine update_periodic

end module m_UpDate
!====================================================================