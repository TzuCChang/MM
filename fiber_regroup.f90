!====================================================================
module m_FiberRegroup   !2018/07/21  change name

use m_DataStructures
use m_UtilityLib

implicit none
contains

!*===================================================================
 subroutine fiber_regroup( fibers,&                   !2018/07/21 change name
                           hinges,&
                           ghost_segments,&
                           E_Young,&          !ERROR 2018/07/07
                           Inertia_Moment,&   !ERROR 2018/07/07
                           allow_breakage,&
                           min_curv,&
                           r_fiber,&
                           box_size,&
                           cells,&
                           nbr_neighbors,&
                           neighbor_list,&
                           distance_neighbors,&
                           gamma_dot,&
                           t,&                     !2018/07/21  time change to t
                           Nbr_bins,&              !2018/07/21  add
                           distanceFactor, &
                           simParameters )

implicit none
type(fiber),   dimension(:), allocatable :: fibers
type(rod),     dimension(:), allocatable :: hinges
type(segment), dimension(:), allocatable :: ghost_segments
type(cell),    dimension(:), allocatable :: cells
type(simulationParameters)               :: simParameters
logical                                  :: allow_breakage

integer(8), dimension(:,:), allocatable  :: neighbor_list
integer,    dimension(:,:), allocatable  :: hinge_index
integer,    dimension(:),   allocatable  :: ghost_segment_index
integer,    dimension(3)                 :: indx, Nbr_bins
integer(8)                               :: i, j, k, m, nbr_segments, nbr_neighbors, ind, numClones !error 2018/07/12 integer(8)

real(8), dimension(:,:), allocatable     :: distance_neighbors
real(8), dimension(3)                    :: r, box_size,min_coor, max_coor , center
real(8)                                  :: t, dX, gamma_dot, distanceFactor, finish2, finish3, start2,  start3
real(8)                                  :: E_Young, Inertia_Moment, min_curv, r_fiber, max_length, bin_length

!-----------------------------------------------------------------------------

if( simParameters%IsPeriodicY ) then
    numClones= 9
else
    numClones= 5
end if

nbr_segments= 0 

! Count the number of segments
do i=1, ubound(fibers,1)       
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        nbr_segments= nbr_segments+1
    end do
end do

if( allocated(ghost_segments) ) deallocate( ghost_segments )	
allocate( ghost_segments(nbr_segments*numClones) ) ! change 5 to 9 if you want periodicity in all directions

if( allocated(hinge_index) ) deallocate( hinge_index )
allocate( hinge_index(nbr_segments,2) )

if( allocated(ghost_segment_index) ) deallocate( ghost_segment_index )
allocate( ghost_segment_index(numClones*nbr_segments) ) ! change 5 to 9 if you want periodicity in all directions


k=1

dX= mod( gamma_dot*box_size(2)/2*t, box_size(1) )  ! This is how much the image boxes need to be translated for the lees edward boundaries

do i=1, ubound(fibers,1)   
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        !This is the original segment
        ghost_segments(k)%A=        hinges(j  )%X_i
        ghost_segments(k)%B=        hinges(j+1)%X_i
        ghost_segments(k)%orig_pos(1)=     j
        ghost_segments(k)%orig_pos(2)=     j+1
        ghost_segments(k)%axis_loc=        1
        k= k + 1  !2018/08/02  change

        !Image segment in +x
        ghost_segments(k)%A(2:3)= hinges(j  )%X_i(2:3)
        ghost_segments(k)%B(2:3)= hinges(j+1)%X_i(2:3)
        ghost_segments(k)%A(1)=   hinges(j  )%X_i(1)+  box_size(1)
        ghost_segments(k)%B(1)=   hinges(j+1)%X_i(1)+  box_size(1)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      2
        k= k + 1  !2018/08/02  change

        ! Image segment in -x
        ghost_segments(k)%A(2:3)= hinges(j  )%X_i(2:3)
        ghost_segments(k)%B(2:3)= hinges(j+1)%X_i(2:3)
        ghost_segments(k)%A(1)=   hinges(j  )%X_i(1)-  box_size(1)
        ghost_segments(k)%B(1)=   hinges(j+1)%X_i(1)-  box_size(1)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      3
        k= k + 1  !2018/08/02  change

        ! Image segment in +z
        ghost_segments(k)%A(1:2)= hinges(j  )%X_i(1:2)
        ghost_segments(k)%B(1:2)= hinges(j+1)%X_i(1:2)
        ghost_segments(k)%A(3)=   hinges(j  )%X_i(3)+  box_size(3)
        ghost_segments(k)%B(3)=   hinges(j+1)%X_i(3)+  box_size(3)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      4
        k= k + 1  !2018/08/02  change

        ! Image segment in +z
        ghost_segments(k)%A(1:2)= hinges(j  )%X_i(1:2)
        ghost_segments(k)%B(1:2)= hinges(j+1)%X_i(1:2)
        ghost_segments(k)%A(3)=   hinges(j  )%X_i(3)-  box_size(3)
        ghost_segments(k)%B(3)=   hinges(j+1)%X_i(3)-  box_size(3)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      5
        k= k + 1  !2018/08/02  change

        if( simParameters%IsPeriodicY ) then
        ! Image segment in +Ya
        ghost_segments(k)%A(1)=   hinges(j  )%X_i(1)+  dX
        ghost_segments(k)%B(1)=   hinges(j+1)%X_i(1)+  dX
        ghost_segments(k)%A(2)=   hinges(j  )%X_i(2)+  box_size(2)
        ghost_segments(k)%B(2)=   hinges(j+1)%X_i(2)+  box_size(2)
        ghost_segments(k)%A(3)=   hinges(j  )%X_i(3)
        ghost_segments(k)%B(3)=   hinges(j+1)%X_i(3)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      6
        k= k + 1  !2018/08/02  change

        ! Image segment in +Yb
        ghost_segments(k)%A(1)=   hinges(j  )%X_i(1)+  dX - box_size(1)
        ghost_segments(k)%B(1)=   hinges(j+1)%X_i(1)+  dX - box_size(1)
        ghost_segments(k)%A(2)=   hinges(j  )%X_i(2)+  box_size(2)
        ghost_segments(k)%B(2)=   hinges(j+1)%X_i(2)+  box_size(2)
        ghost_segments(k)%A(3)=   hinges(j  )%X_i(3)
        ghost_segments(k)%B(3)=   hinges(j+1)%X_i(3)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      7
        k= k + 1  !2018/08/02  change

        ! Image segment in -Ya
        ghost_segments(k)%A(1)=   hinges(j  )%X_i(1)-  dX
        ghost_segments(k)%B(1)=   hinges(j+1)%X_i(1)-  dX
        ghost_segments(k)%A(2)=   hinges(j  )%X_i(2)-  box_size(2)
        ghost_segments(k)%B(2)=   hinges(j+1)%X_i(2)-  box_size(2)
        ghost_segments(k)%A(3)=   hinges(j  )%X_i(3)
        ghost_segments(k)%B(3)=   hinges(j+1)%X_i(3)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      8
        k= k + 1  !2018/08/02  change

        ! Image segment in -Yb
        ghost_segments(k)%A(1)=   hinges(j  )%X_i(1)-  dX +   box_size(1)
        ghost_segments(k)%B(1)=   hinges(j+1)%X_i(1)-  dX +   box_size(1)
        ghost_segments(k)%A(2)=   hinges(j  )%X_i(2)-  box_size(2)
        ghost_segments(k)%B(2)=   hinges(j+1)%X_i(2)-  box_size(2)
        ghost_segments(k)%A(3)=   hinges(j  )%X_i(3)
        ghost_segments(k)%B(3)=   hinges(j+1)%X_i(3)
        ghost_segments(k)%orig_pos(1)=   j
        ghost_segments(k)%orig_pos(2)=   j+1
        ghost_segments(k)%axis_loc=      9
        k= k + 1  !2018/08/02  change
        end if
    end do
end do


max_length= 0
do i=1, ubound(fibers,1)
    do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        r= hinges(j+1)%X_i - hinges(j)%X_i
        max_length= max( max_length, sqrt(dot_product(r,r)) )  !2018/08/02 找出所有的Segments的長度中最大值
    end do
end do

min_coor=   huge(0d0)
max_coor=  -huge(0d0)
do k=1,3
    do i=1, ubound(ghost_segments,1)
        min_coor(k)= min( min_coor(k), ghost_segments(i)%A(k) )
        min_coor(k)= min( min_coor(k), ghost_segments(i)%B(k) )
	    max_coor(k)= max( max_coor(k), ghost_segments(i)%A(k) )
	    max_coor(k)= max( max_coor(k), ghost_segments(i)%B(k) )
	end do  
end do      !2018/08/02 找出可以容納所有的segments框架的座標範圍

!bin_length= 1.25*(max_length+2*r_fiber)!Mod 9/28/2014

bin_length= 1.5*( max_length + 2*r_fiber )

!min_coor= min_coor- 2*r_fiber
!max_coor= max_coor+ 2*r_fiber

min_coor= min_coor- 2*r_fiber-max_length  !2018/08/02 再擴大框架的座標範圍
max_coor= max_coor+ 2*r_fiber+max_length

Nbr_bins= ceiling( (max_coor-min_coor)/bin_length )

!bin_length=(max_coor(3)-min_coor(3))/Nbr_bins(3)
!min_coor=min_coor-bin_length
!max_coor=min_coor+bin_length      !Mod 9/28/2014
!Nbr_bins=Nbr_bins+2

min_coor= min_coor - 2*bin_length  !2018/08/02 再進一步擴大框架的座標範圍
max_coor= max_coor + 2*bin_length  !2018/08/02 再進一步擴大框架的座標範圍

Nbr_bins= Nbr_bins + 4

!print *,"Nbr bins", Nbr_bins

do j=1,ubound(hinges,1)
    do k=1,3
        hinges(j)%indx(k)= 0  !2018/08/02 先設定給起始值,歸零
    end do
end do   
                                !2018/08/02   floor( 3.412)=  3.00000  取整數
do i=1, ubound(fibers,1)        !2018/08/02   floor(-3.412)= -4.00000  取整數
    do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
       center= ( hinges(j)%X_i + hinges(j+1)%X_i )/2                !2018/08/02  change
       do k=1,3                           
            indx(k)= floor( (center(k)-min_coor(k))/bin_length ) + 1 !2018/08/02  change
            if( indx(k).le.2 ) indx(k)= 2    
            if( indx(k).ge.(Nbr_bins(k)-1) ) indx(k)= Nbr_bins(k)-1
	   end do
	   hinges(j)%ind=(indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
    end do
end do

m=1
do j=1, ubound(ghost_segments, 1)
   center= ( ghost_segments(j)%A + ghost_segments(j)%B )/2     !2018/08/02 change
   do k=1,3
      indx(k)= floor( (center(k)-min_coor(k))/bin_length )+1   !2018/08/02 change
      if(indx(k).le.2) indx(k)=2
      if(indx(k).ge.(Nbr_bins(k)-1)) indx(k)= Nbr_bins(k)-1
   end do   
   ghost_segment_index(m)=(indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
   m=m+1
end do

call QsortC( ghost_segment_index , ghost_segments )

if( allocated(cells) ) deallocate(cells)

allocate( cells(Nbr_bins(1)*Nbr_bins(2)*Nbr_bins(3)) )

do i=1, ubound(cells,1)
    cells(i)%ghost_limits(1)= 0
    cells(i)%ghost_limits(2)= 0
end do

do i=1, Nbr_bins(1)
    do j=1, Nbr_bins(2)
        do k=1, Nbr_bins(3)
            ind= (k-1)*Nbr_bins(1)*Nbr_bins(2) + (j-1)*Nbr_bins(1) + i
            cells(ind)%indx(1)= i
            cells(ind)%indx(2)= j
            cells(ind)%indx(3)= k
        end do
    end do
end do

cells(ghost_segment_index(1))%ghost_limits(1)= 1

do j=2, ubound(ghost_segment_index,1)
    if( ghost_segment_index(j) .ne. ghost_segment_index(j-1) ) then
        cells(ghost_segment_index(j  ))%ghost_limits(1)= j
        cells(ghost_segment_index(j-1))%ghost_limits(2)= j-1
    end if
end do

end subroutine fiber_regroup  
                                                                 
end module m_FiberRegroup   !2018/07/21  change name

!====================================================================