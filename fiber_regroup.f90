!====================================================================
module m_FiberRegroup   !2018/07/21  change name

use m_DataStructures
use m_UtilityLib

implicit none
contains

!*===================================================================
 subroutine fiber_regroup( fibers,&                   !2018/10/09  �ץ�
                           hinges,&
                           ghost_segments,&
                           neighbor_list,&
                           cells,&
                           simParameters )

implicit none
type(simulationParameters)               :: simParameters
type(fiber),   dimension(:), allocatable :: fibers
type(rod),     dimension(:), allocatable :: hinges
type(segment), dimension(:), allocatable :: ghost_segments
type(cell),    dimension(:), allocatable :: cells
logical                                  :: allow_breakage

integer(8), dimension(:,:), allocatable  :: neighbor_list
integer,    dimension(:,:), allocatable  :: hinge_index
integer,    dimension(:),   allocatable  :: ghost_segment_index
integer,    dimension(3)                 :: indx, Nbr_bins
integer(8)                               :: i, j, k, m, nbr_segments, nbr_neighbors, ind, nbr_GhostSegments !error 2018/07/12 integer(8)

real(8), dimension(3)                    :: r, box_size,min_coor, max_coor , center
real(8)                                  :: distanceFactor, finish2, finish3, start2,  start3
real(8)                                  :: r_fiber, max_length, bin_length

r_fiber  = simParameters%r_fiber            !2018/10/09  �ץ�
box_size = simParameters%box_size           !2018/10/09  �ץ�
nbr_neighbors= simParameters%nbr_neighbors  !2018/10/10  �ץ�

     nbr_GhostSegments= ubound(ghost_segments,1)

     nbr_segments= 0 
     do i=1, ubound(fibers,1)
        nbr_segments= nbr_segments + fibers(i)%nbr_hinges - 1 !2018/09/08  Count the number of segments
     end do

     if( allocated(hinge_index) ) deallocate( hinge_index )
         allocate( hinge_index(nbr_segments,2) )

     if( allocated(ghost_segment_index) ) deallocate( ghost_segment_index )
         allocate( ghost_segment_index(nbr_GhostSegments) ) 

     max_length= 0
     do i= 1, ubound(fibers,1)
     do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        r= hinges(j+1)%X_i - hinges(j)%X_i
        max_length= max( max_length, sqrt(dot_product(r,r)) )  !2018/08/02 ��X�Ҧ���Segments�����פ��̤j��
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
     end do      !2018/08/02 ��X�i�H�e�ǩҦ���segments�ج[���y�нd��

    !bin_length= 1.25*( max_length + 2*r_fiber )!Mod 9/28/2014
    !min_coor= min_coor- 2*r_fiber
    !max_coor= max_coor+ 2*r_fiber
    
     bin_length= 1.5*( max_length + 2*r_fiber )
     min_coor= min_coor- 2*r_fiber-max_length  !2018/08/02 �A�X�j�ج[���y�нd��
     max_coor= max_coor+ 2*r_fiber+max_length
     Nbr_bins= ceiling( (max_coor-min_coor)/bin_length )  !2018/08/05 CEILING ( 4.80) has the value  5
                                                          !2018/08/05 CEILING (-2.55) has the value -2
    !bin_length= (max_coor(3)-min_coor(3))/Nbr_bins(3)
    !min_coor= min_coor-bin_length
    !max_coor= min_coor+bin_length        !Mod 9/28/2014
    !Nbr_bins= Nbr_bins+2

     min_coor= min_coor - 2*bin_length  !2018/08/02 �A�i�@�B�X�j�ج[���y�нd��
     max_coor= max_coor + 2*bin_length  !2018/08/02 �A�i�@�B�X�j�ج[���y�нd��
     Nbr_bins= Nbr_bins + 4             !2018/08/05 �@�W�[4��bin_length���Ŷ�
                                        !2018/08/05 �w��Ҧ�ghost_segments,��X�]�e���Ŷ�,
                                        !2018/08/05 �åB��bin_length���פj�p,�N����e���T�Ӥ�V,
                                        !2018/08/05 ���O���Φ��\�h������p���l,
                                        !2018/08/05 �ƶq��Nbr_bin(1)*Nbr_bin(2)*Nbr_bin(3)
    !print *,"Nbr bins", Nbr_bins

     do j=1,ubound(hinges,1)
        hinges(j)%indx= 0               !2018/08/05  ���]�w���_�l��,�k�s
     end do   
                                        !2018/08/02   floor( 3.412)=  3.00000  �����
     do i= 1, ubound(fibers,1)          !2018/08/02   floor(-3.412)= -4.00000  �����
     do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        center= ( hinges(j)%X_i + hinges(j+1)%X_i )/2                 !2018/08/02  change
        do k=1,3                           
            indx(k)= floor( (center(k)-min_coor(k))/bin_length ) + 1  !2018/08/02  change
            if( indx(k) .le. 3 )           indx(k)= 2    
            if( indx(k) .ge. Nbr_bins(k) ) indx(k)= Nbr_bins(k)-1
	    end do
	    hinges(j)%ind= (indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
     end do       !2018/08/05 �P�_hinges(j)~(j+1)��segments����q���~�I��m���b���Ŷ��p���l��
     end do

     m= 1
     do j= 1, ubound(ghost_segments, 1)
        center= ( ghost_segments(j)%A + ghost_segments(j)%B )/2     !2018/08/02 change
     do k=1,3
        indx(k)= floor( (center(k)-min_coor(k))/bin_length )+1   !2018/08/02 change
        if( indx(k) .le. 3 )           indx(k)= 2                !2018/08/05 change
        if( indx(k) .ge. Nbr_bins(k) ) indx(k)= Nbr_bins(k)-1    !2018/08/05 change
     end do   
        ghost_segment_index(m)= (indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
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


     cells( ghost_segment_index(1) )%ghost_limits(1)= 1

     do j= 2, ubound(ghost_segment_index,1)
        if( ghost_segment_index(j) .ne. ghost_segment_index(j-1) ) then
            cells( ghost_segment_index(j)   )%ghost_limits(1)= j
            cells( ghost_segment_index(j-1) )%ghost_limits(2)= j-1
        end if
     end do
     
     simParameters%Nbr_bins= Nbr_bins            !2018/10/10  �ץ�

!    print *,Nbr_bins(1), Nbr_bins(2), Nbr_bins(3), ubound(cells,1)
!    do i=1, Nbr_bins(1)
!    do j=1, Nbr_bins(2)
!        do k=1, Nbr_bins(3)
!            ind= (k-1)*Nbr_bins(1)*Nbr_bins(2) + (j-1)*Nbr_bins(1) + i
!            if( cells(ind)%ghost_limits(1) .ne. 0 ) then
!              print *, j,k , cells(ind)%ghost_limits(1), cells(ind)%ghost_limits(2)-cells(ind)%ghost_limits(1)+1
!            end if
!        end do
!    end do
!    print *, i    
!    pause
!    end do
!    pause

end subroutine fiber_regroup 

                           
subroutine fiber_regroup_ShiftCenterToOrigion( fibers, hinges, simParameters ) !2018/08/05 add

type(fiber), allocatable, dimension(:)   :: fibers
type(rod)  , allocatable, dimension(:)   :: hinges
type(simulationParameters)               :: simParameters

real(8), dimension(3)     :: box_size, coord, min_coor, max_coor
integer(8)                :: i, j, k, m  

box_size= simParameters%box_size

!2018/08/05  Compute center of mass for all fibers coord

print *,"@@@"
print *,"based on all hinges"
m= 0
coord= 0
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
   coord= coord + hinges(j)%X_i
   m= m+1
end do
end do        
coord= coord/real(m)  !2018/08/05  coord= center of mass 
print *,"cen ", coord

!2018/08/05 shift center of mass to (0,0,0) 
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
   hinges(j)%X_i= hinges(j)%X_i - coord
end do
end do   

m= 0
coord= 0
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
   coord= coord + hinges(j)%X_i
   m= m+1
end do
end do        
coord= coord/real(m)  !2018/08/05  coord= center of mass 
print *,"cen ", coord
!pause

end subroutine fiber_regroup_ShiftCenterToOrigion 


end module m_FiberRegroup   !2018/07/21  change name

!====================================================================