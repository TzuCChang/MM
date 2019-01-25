!====================================================================
module m_UpDate

use m_DataStructures
implicit none
contains
           
subroutine update_periodic( fibers, hinges, simParameters )  !2018/10/09 change

type(simulationParameters)            :: simParameters
type(fiber) , dimension(:)            :: fibers
type (rod),   dimension(:)            :: hinges
logical                               :: periodic_boundary
integer(8)                            :: i,j,k,l,jp1, iStep, nStep, mStep, nStep_max
real(8),      dimension(3)            :: box_size, coord
real(8)                               :: dt, cx, cy, cz, dist2, dX, t, gamma_dot, dist2_max, displ_max, ratio, dt_Step

periodic_boundary = simParameters%periodic_boundary

box_size  = simParameters%box_size
gamma_dot = simParameters%gamma_dot
t         = simParameters%time
dt        = simParameters%dt


do i=1, ubound (fibers,1)   !2018/09/22  add
	do j=fibers(i)%first_hinge, (fibers(i)%first_hinge+fibers(i)%nbr_hinges-2)
		coord= hinges(j+1)%X_i-hinges(j)%X_i
		hinges(j)%length2= sqrt( dot_product(coord, coord) )
	end do
end do


nStep_max = simParameters%nStep_max
displ_max = simParameters%displ_max   !2018/10/27 segments length= 0.10 mm為半徑,走一個步的距離,圓周長的3,600分之一
nStep=      simParameters%nStep

mStep= nStep
dist2_max= 0.0d0

do i=1, ubound(hinges,1)       !2018/10/27 修正
	if( hinges(i)%is_stationary==0 ) then
		coord= hinges(i)%v_i*dt 
        dist2= sqrt( dot_product(coord, coord) )
        dist2_max= max( dist2_max, dist2 )
	end if
end do

if( dist2_max .lt. 1.0d-10 ) then
    dist2_max= 1.0d-10
end if

ratio= displ_max/dist2_max

if( ratio .gt. 1.d0 ) then
    ratio= 1.d0
end if

iStep= int(ratio*nStep_max)

if( iStep .lt. 1 ) then
    iStep= 1
else if( iStep .gt. nStep_max ) then
    iStep= nStep_max
end if


!   write(*,61 ) "1 iStep=  ", iStep, nStep, mStep, nStep_max
!61 format(A12, 4I6 )

mStep= nStep + iStep

if( mStep .gt. nStep_max ) then
    mStep= nStep_max
end if

!   write(*,62 ) "2 iStep=  ", iStep, nStep, mStep, nStep_max
!62 format(A12, 4I6 )

iStep= mStep - nStep

nStep= nStep + iStep

simParameters%nStep= nStep

dt_Step= dt*real(iStep)/real(nStep_max)

simParameters%time= simParameters%time + dt_Step      !2018/10/29  修正                     

!   write(*,63 ) "3 iStep=  ", iStep, nStep, mStep, nStep_max
!63 format(A12, 4I6 )
!   write(*,64 ) "@@@( dt = ", dt_Step, iStep, nStep, ratio, displ_max, dist2_max
!64 format(A12, D15.5, I5, I5, F9.5, 2D15.5 )
!   if( nStep .eq. nStep_max ) then               !2018/10/28 add
!       write(*,65 ) "2 time=  ", simParameters%time*1.e6
!65     format(A12, F12.2 )
!      pause
!   end if

do i=1, ubound(hinges,1)
	if( hinges(i)%is_stationary==0 ) then
		hinges(i)%X_i= hinges(i)%X_i + hinges(i)%v_i*dt_Step  !2018/10/27 修正
	end if
end do

do i=1, ubound (fibers,1)  !2018/09/22  add
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
       jp1= j+1
       coord= hinges(jp1)%X_i - hinges(j)%X_i
       coord= coord/(sqrt(dot_product(coord, coord)))
       hinges(jp1)%X_i= hinges(j)%X_i + coord*hinges(j)%length2
	end do
end do


dX = mod( gamma_dot*t*box_size(2)/2.d0, box_size(1) )    !2018/09/22 修正

if( periodic_boundary .eqv. .true. ) then

	do i=1, ubound (fibers,1)

		coord=0
		k=0
        ! Compute center of mass for fiber i
		do l=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
			coord= coord + ( (hinges(l)%X_i+hinges(l+1)%X_i)/2d0 )
			k=k+1
		end do
		coord= coord/real(k)
        
        ! Apply boundary conditions
        ! if the center of mass of the fiber is outside of the box, move it accorddingly. 
        ! based loosely on Non-Newtonian molecular Dynamics by Evans and Morriss, appendix B
        
        cy= int( coord(2) / (box_size(2)/2d0) )
  
        coord(1)= coord(1) - dX*cy
        
        cx= int( coord(1)/(box_size(1)/2d0) )
        cz= int( coord(3)/(box_size(3)/2d0) )
        
        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
            if( hinges(j)%is_stationary==0 ) then     !2018/07/20   Error  hinges(i) change to hinges(j)
                 hinges(j)%X_i(1)= hinges(j)%X_i(1) - dX*cy - cx*box_size(1)
                 hinges(j)%X_i(2)= hinges(j)%X_i(2) - cy*box_size(2)
                 hinges(j)%X_i(3)= hinges(j)%X_i(3) - cz*box_size(3)
            end if
        end do
       
!        
!      if the coordinates of the center of mass are larger or smaller than the limits of the periodic box, move the fibers accordingly.
!

		coord=0
		k=0
        ! Compute center of mass for fiber i
		do l=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
			coord= coord + ( (hinges(l)%X_i+hinges(l+1)%X_i)/2d0 )
			k=k+1
		end do
		coord= coord/real(k)
        
        do l=1,3
	       if( (coord(l)) .gt. (box_size(l)/2d0) ) then
                do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
                   if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)-box_size(l)  !2018/09/22  error hinges(i)
                end do
           else if (coord(l).lt. (-box_size(l)/2d0)) then 
		        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
			       if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)+box_size(l) !2018/09/22  error hinges(i)
                end do
           end if
        end do
        
		coord=0
		k=0
        ! Compute center of mass for fiber i
		do l=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
			coord= coord + ( (hinges(l)%X_i+hinges(l+1)%X_i)/2d0 )
			k=k+1
		end do
		coord= coord/real(k)
        
        do l=1,3
	       if( (coord(l)) .gt. (box_size(l)/2d0) ) then
                do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
                   if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)-box_size(l)  !2018/09/22  error hinges(i)
                end do
           else if (coord(l).lt. (-box_size(l)/2d0)) then 
		        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
			       if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)+box_size(l) !2018/09/22  error hinges(i)
                end do
           end if
        end do
        
    end do
    
end if

end subroutine update_periodic


subroutine update_periodic_Initial( fibers, hinges, simParameters ) !2018/10/09  add

type(simulationParameters)            :: simParameters
type(fiber) , dimension(:)            :: fibers
type (rod), dimension(:)              :: hinges
logical                               :: periodic_boundary
integer(8)                            :: i,j,k,l,jp1
real(8), dimension(3)                 :: box_size, coord
real(8)                               :: cx, cy, cz

periodic_boundary = simParameters%periodic_boundary
box_size          = simParameters%box_size

do i=1, ubound (fibers,1)   !2018/09/22  add
	do j=fibers(i)%first_hinge, (fibers(i)%first_hinge+fibers(i)%nbr_hinges-2)
		coord= hinges(j+1)%X_i-hinges(j)%X_i
		hinges(j)%length2= sqrt( dot_product(coord, coord) )
	end do
end do


do i=1, ubound(hinges,1)
	if( hinges(i)%is_stationary==0 ) then
		hinges(i)%X_i= hinges(i)%X_i + hinges(i)%v_i*0.d0  !2018/09/22 修正
	end if
end do

do i=1, ubound (fibers,1)  !2018/09/22  add
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
       jp1= j+1
       coord= hinges(jp1)%X_i - hinges(j)%X_i
       coord= coord/(sqrt(dot_product(coord, coord)))
       hinges(jp1)%X_i= hinges(j)%X_i + coord*hinges(j)%length2
	end do
end do


if( periodic_boundary .eqv. .true. ) then

	do i=1, ubound (fibers,1)

		coord=0
		k=0
        ! Compute center of mass for fiber i
		do l=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
			coord= coord + ( (hinges(l)%X_i+hinges(l+1)%X_i)/2d0 )
			k=k+1
		end do
		coord= coord/real(k)
        
        ! Apply boundary conditions
        ! if the center of mass of the fiber is outside of the box, move it accorddingly. 
        ! based loosely on Non-Newtonian molecular Dynamics by Evans and Morriss, appendix B

        cx= int( coord(1)/(box_size(1)/2d0) )
        cy= int( coord(2)/(box_size(2)/2d0) )
        cz= int( coord(3)/(box_size(3)/2d0) )
        
        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
            if( hinges(j)%is_stationary==0 ) then     !2018/07/20   Error  hinges(i) change to hinges(j)
                 hinges(j)%X_i(1)= hinges(j)%X_i(1) - cx*box_size(1)
                 hinges(j)%X_i(2)= hinges(j)%X_i(2) - cy*box_size(2)
                 hinges(j)%X_i(3)= hinges(j)%X_i(3) - cz*box_size(3)
            end if
        end do
       
!        
!      if the coordinates of the center of mass are larger or smaller than the limits of the periodic box, move the fibers accordingly.
!

		coord=0
		k=0
        ! Compute center of mass for fiber i
		do l=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
			coord= coord + ( (hinges(l)%X_i+hinges(l+1)%X_i)/2d0 )
			k=k+1
		end do
		coord= coord/real(k)
        
        do l=1,3
	       if( (coord(l)) .gt. (box_size(l)/2d0) ) then
                do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
                   if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)-box_size(l)  !2018/09/22  error hinges(i)
                end do
           else if (coord(l).lt. (-box_size(l)/2d0)) then 
		        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
			       if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)+box_size(l) !2018/09/22  error hinges(i)
                end do
           end if
        end do
        
		coord=0
		k=0
        ! Compute center of mass for fiber i
		do l=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
			coord= coord + ( (hinges(l)%X_i+hinges(l+1)%X_i)/2d0 )
			k=k+1
		end do
		coord= coord/real(k)
        
        do l=1,3
	       if( (coord(l)) .gt. (box_size(l)/2d0) ) then
                do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
                   if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)-box_size(l)  !2018/09/22  error hinges(i)
                end do
           else if (coord(l).lt. (-box_size(l)/2d0)) then 
		        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
			       if (hinges(j)%is_stationary==0) hinges(j)%X_i(l)=hinges(j)%X_i(l)+box_size(l) !2018/09/22  error hinges(i)
                end do
           end if
        end do
        
    end do
    
end if

end subroutine update_periodic_Initial

!*************************************************
subroutine update_robust(fibers, hinges, dt)
implicit none
type(fiber) , dimension(:)            :: fibers
type (rod), dimension(:)              :: hinges
integer(8)                            :: i,j
real(8), dimension(3)                 :: dir_vec
real(8)                               :: dt, deltaX, deltaY, delta

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
integer(8)                            :: i
real(8)                               :: dt

do i=1, ubound(hinges,1)
	if (hinges(i)%is_stationary==0) then
		hinges(i)%X_i=hinges(i)%X_i+hinges(i)%v_i*dt
	end if
end do

end subroutine update
!*************************************************

end module m_UpDate
!====================================================================