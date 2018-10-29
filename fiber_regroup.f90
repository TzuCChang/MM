!====================================================================
module m_FiberRegroup   !2018/07/21  change name

! subroutine: fiber_regroup
! subroutine: fiber_regroup_ShiftCenterToOrigion
! subroutine: fiber_regroup_minmax_hinges
! subroutine: fiber_regroup_minmax_segments

use m_DataStructures    !m �N��module(�[�\�o��subroutine) !Datastructure�btipo��
use m_UtilityLib

implicit none   !implicit none�N���ǪF�褣�ŧi
contains

!*===================================================================
 subroutine fiber_regroup( fibers,&    !��i���               !2018/07/21 change name
                           hinges,&    !�Nhinge�s��
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
type(rod),     dimension(:), allocatable :: hinges  !�@��dimension
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
! �}�l�p��g��
! �p�Gy���g��:�n��27��box(�T�תŶ�) / �Yy�S���g��(���):�n��9��box(�����r)
if( simParameters%IsPeriodicY ) then    !%�N��type���U���ܼ� �ӳo���ܼƬOlogical
    numClones= 27    !2018/08/05 �ץ� change from 9 to 27
else
    numClones= 9     !2018/08/05 �ץ� change from 5 to 9
end if

nbr_segments= 0 !��l����=0

! Count the number of segments
do i=1, ubound(fibers,1)    !upperbound: fibers��1�� �ݽ֤j       
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !�Ĥ@��hinge����̫�@��hinge >>�]�����hinge=�@��segment
        nbr_segments= nbr_segments+1
    end do
end do

! ALLOCATE : creates space for allocatable arrays and variables with the POINTER attribute. 
! DEALLOCATE : frees space previously allocated for allocatable arrays and pointer targets.
! ���p���ŧi, �������ŧi, �A���s�ŧi�@��>>�T�O�������L�v�T  !�]��ghost_segments�|�@����s�ƭȥB�Q���ƧQ��(�Yfiber�_��, �ݭn�@�ӷs���y�z)
if( allocated(ghost_segments) ) deallocate( ghost_segments )	
allocate( ghost_segments(nbr_segments*numClones) ) ! change 5 to 9 if you want periodicity in all directions

if( allocated(hinge_index) ) deallocate( hinge_index )
allocate( hinge_index(nbr_segments,2) )

if( allocated(ghost_segment_index) ) deallocate( ghost_segment_index )
allocate( ghost_segment_index(numClones*nbr_segments) ) ! change 5 to 9 if you want periodicity in all directions

!dX ���Ȯɤ��z(�|�ݽT�{)
dX= mod( gamma_dot*box_size(2)/2*t, box_size(1) )  ! This is how much the image boxes need to be translated for the lees edward boundaries
!shear rate*y=�t�� �t��/2=�����t�� �����t��*t=�첾
m= 0    !m=�Ҧ�segment���ƥ� (��l)

do i=1, ubound(fibers,1)   !�q�Ĥ@��fiber�}�l
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !���Ҧ���segment
        !�C�@���ֺ���first hinge���@�w���O1 �ҥH�n�����fiber�Ĥ@��hinge �h�����᭱��
        k= j+1  !fiber(i)��hinge(���I)���ƶq
       
        !This is the original segment
        m= m + 1  !2018/08/05  changet  !m�Osegment �C�@��m+1�N�O�ƻssegment��Ŷ���cell
        ghost_segments(m)%axis_loc= 1   !1�N����, 2+x, 3-x, 4+z, 5_z, 6+x+z, 7-x+z, 8+x-z, 9-x-z
        ghost_segments(m)%orig_pos(1)= j    
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A= hinges(j)%X_i  !segment(rod)
        ghost_segments(m)%B= hinges(k)%X_i  !hinge���y���I(rod�Y���)

        !Image segment in +x
        m= m + 1  !2018/08/05  changet
        ghost_segments(m)%axis_loc= 2
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) + box_size(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) + box_size(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3)
        
        ! Image segment in -x
        m= m + 1  !2018/08/05  changet
        ghost_segments(m)%axis_loc= 3
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) - box_size(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) - box_size(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3)
        

        ! Image segment in +z
        m= m + 1  !2018/08/05  changet
        ghost_segments(m)%axis_loc= 4
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k        
        ghost_segments(m)%A(1)= hinges(j)%X_i(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3) + box_size(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3) + box_size(3)

        ! Image segment in -z
        m= m + 1  !2018/08/05  changet        
        ghost_segments(m)%axis_loc= 5
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k        
        ghost_segments(m)%A(1)= hinges(j)%X_i(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3) - box_size(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3) - box_size(3)
        
        ! �s�W�|�ӱר�(+x+z, +x-z, -x+z, -x-z)
        ! Image segment in +x+z
        m= m + 1  !2018/08/05  changet        
        ghost_segments(m)%axis_loc= 6
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) + box_size(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3) + box_size(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) + box_size(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3) + box_size(3)
        
        ! Image segment in -x+z
        m= m + 1  !2018/08/05  changet        
        ghost_segments(m)%axis_loc= 7
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) - box_size(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3) + box_size(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) - box_size(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3) + box_size(3)

        ! Image segment in -x-z
        m= m + 1  !2018/08/05  changet        
        ghost_segments(m)%axis_loc= 8
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) - box_size(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3) - box_size(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) - box_size(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3) - box_size(3)

        ! Image segment in +x-z
        m= m + 1  !2018/08/05  changet        
        ghost_segments(m)%axis_loc= 9
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) + box_size(1)
        ghost_segments(m)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3) - box_size(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) + box_size(1)
        ghost_segments(m)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3) - box_size(3)

        if( simParameters%IsPeriodicY ) then
        ! Image segment in +Ya
        m= m + 1  !2018/08/05  changet
        ghost_segments(m)%axis_loc= 10
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) + dX
        ghost_segments(m)%A(2)= hinges(j)%X_i(2) + box_size(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) + dX
        ghost_segments(m)%B(2)= hinges(k)%X_i(2) + box_size(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3)

        ! Image segment in +Yb
        m= m + 1  !2018/08/05  changet  
        ghost_segments(m)%axis_loc= 11
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) - box_size(1) + dX
        ghost_segments(m)%A(2)= hinges(j)%X_i(2) + box_size(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) - box_size(1) + dX
        ghost_segments(m)%B(2)= hinges(k)%X_i(2) + box_size(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3)

        ! Image segment in -Ya
        m= m + 1  !2018/08/05  changet  
        ghost_segments(m)%axis_loc= 12
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k           
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) - dX
        ghost_segments(m)%A(2)= hinges(j)%X_i(2) - box_size(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) - dX
        ghost_segments(m)%B(2)= hinges(k)%X_i(2) - box_size(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3)        

        ! Image segment in -Yb
        m= m + 1  !2018/08/05  changet
        ghost_segments(m)%axis_loc= 13
        ghost_segments(m)%orig_pos(1)= j
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A(1)= hinges(j)%X_i(1) + box_size(1) - dX
        ghost_segments(m)%A(2)= hinges(j)%X_i(2) - box_size(2)      
        ghost_segments(m)%A(3)= hinges(j)%X_i(3)
        ghost_segments(m)%B(1)= hinges(k)%X_i(1) + box_size(1) - dX
        ghost_segments(m)%B(2)= hinges(k)%X_i(2) - box_size(2)        
        ghost_segments(m)%B(3)= hinges(k)%X_i(3)
        end if
    end do
end do



max_length= 0
do i= 1, ubound(fibers,1)
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !��segment
        r= hinges(j+1)%X_i - hinges(j)%X_i  !r=segment���V�q
        max_length= max( max_length, sqrt(dot_product(r,r)) )  !2018/08/02 ��X�Ҧ���Segments�����פ��̤j��
    end do
end do

min_coor=   huge(0d0)   ! huge: returns the largest number that is not an infinity
max_coor=  -huge(0d0)
do k=1,3    !k=1 2 3 �N��x y z�y��  !�o��O��segment��(�i��|�y��������ɶ��L��>>�]���������segment��)
    do i=1, ubound(ghost_segments,1)    !��X�Y�����I(A B) ���̤j�ȩM�̤p��
        min_coor(k)= min( min_coor(k), ghost_segments(i)%A(k) ) !A�Oj(�Y) 
        min_coor(k)= min( min_coor(k), ghost_segments(i)%B(k) ) !B�Ok(��)
	    max_coor(k)= max( max_coor(k), ghost_segments(i)%A(k) )
	    max_coor(k)= max( max_coor(k), ghost_segments(i)%B(k) )
	end do  
end do      !2018/08/02 ��X�i�H�e�ǩҦ���segments�ج[���y�нd��

! �ݴ��� �ٻݥ[box����
!do k=1,3 
!do i= 1, ubound(fibers,1)
!    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2 
!        min_coor(k)= min( min_coor(k), hinges(j)%X_i(k) )
!        max_coor(k)= max( max_coor(k), hinges(j)%X_i(k) )       
!    end do
!end do
!min_coor= min_coor-box_size
!max_coor= max_coor+box_size

!bin_length= 1.25*( max_length + 2*r_fiber )!Mod 9/28/2014

bin_length= 1.5*( max_length + 2*r_fiber )  !����X�Ҧ�segment���d�� �A������ghost�b�Ŷ������d�򰵤���(����u�n�b�ت��a�����N�i�H�F)bin_length��15�Ӫ��|���k
! �q�`aspect ratio=5��max length �ҥHmax length+���|=���������| *1.5=9 ������10�Ӫ��|

!min_coor= min_coor- 2*r_fiber
!max_coor= max_coor+ 2*r_fiber

min_coor= min_coor- 2*r_fiber-max_length  !2018/08/02 �A�X�j�ج[���y�нd��
max_coor= max_coor+ 2*r_fiber+max_length    !���W���U�U�W�[�@��segment������

Nbr_bins= ceiling( (max_coor-min_coor)/bin_length )  !2018/08/05 CEILING ( 4.80) has the value  5
                                                     !2018/08/05 CEILING (-2.55) has the value -2
!��bin length�h�� ��X�`�@�n���X��bin ceiling�����̤j�����
!ceiling=�̱��񪺥k����

!bin_length= (max_coor(3)-min_coor(3))/Nbr_bins(3)
!min_coor= min_coor-bin_length
!max_coor= min_coor+bin_length        !Mod 9/28/2014
!Nbr_bins= Nbr_bins+2

!���s���@�� ���I�A�ݤ@��
!�����l���ت�>>���F����������٪Ŷ� ��neighbor�u�n��P��bin(���O�쥻�g�k�O������ �ҥH�L�k�`�ٮɶ�)
min_coor= min_coor - 2*bin_length  !2018/08/02 �A�i�@�B�X�j�ج[���y�нd��
max_coor= max_coor + 2*bin_length  !2018/08/02 �A�i�@�B�X�j�ج[���y�нd��
Nbr_bins= Nbr_bins + 4             !2018/08/05 �@�W�[4��bin_length���Ŷ�
                                   !2018/08/05 �w��Ҧ�ghost_segments,��X�]�e���Ŷ�,
                                   !2018/08/05 �åB��bin_length���פj�p,�N����e���T�Ӥ�V,
                                   !2018/08/05 ���O���Φ��\�h������p���l,
                                   !2018/08/05 �ƶq��Nbr_bin(1)*Nbr_bin(2)*Nbr_bin(3) (1 2 3 ���[4��)
!print *,"Nbr bins", Nbr_bins

do j=1,ubound(hinges,1)
   hinges(j)%indx= 0            !2018/08/05  ���]�w���_�l��,�k�s
end do

!�����u�ꪺfiber                 !2018/08/02   floor( 3.412)=  3.00000  �����
do i= 1, ubound(fibers,1)       !2018/08/02   floor(-3.412)= -4.00000  �����
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
       center= ( hinges(j)%X_i + hinges(j+1)%X_i )/2                 !2018/08/02  change    !����Ӹ`�I��(segment)�������I
       do k=1,3                           
            indx(k)= floor( (center(k)-min_coor(k))/bin_length ) + 1 !2018/08/02  change    !�����I/bin length�i�H�o�����b����bin�̭� !indx=???
            if( indx(k) .le. 3 )           indx(k)= 2    !�Ybin<3 �h���w�q�ĤG��bin�}�l(���q1�}�l) >>���O�u�g�k
            if( indx(k) .ge. Nbr_bins(k) ) indx(k)= Nbr_bins(k)-1
	   end do
	   hinges(j)%ind= (indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
    end do       !2018/08/05 �P�_hinges(j)~(j+1)��segments����q�����I��m���b���@�Ӫ��Ŷ��p���l�� ind>>�Τ@����ܤT��
end do

!�����ƻs�᪺ghost segment
m=1
do j= 1, ubound(ghost_segments, 1)
   center= ( ghost_segments(j)%A + ghost_segments(j)%B )/2     !2018/08/02 change
   do k=1,3
      indx(k)= floor( (center(k)-min_coor(k))/bin_length )+1   !2018/08/02 change
      !floor: returns the greatest integer less than or equal to X
      if( indx(k) .le. 3 )           indx(k)= 2                !2018/08/05 change
      if( indx(k) .ge. Nbr_bins(k) ) indx(k)= Nbr_bins(k)-1    !2018/08/05 change
   end do   
   ghost_segment_index(m)= (indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)  
   m=m+1
end do

call QsortC( ghost_segment_index , ghost_segments ) !�j�M �Ƨ�

if( allocated(cells) ) deallocate(cells)

allocate( cells(Nbr_bins(1)*Nbr_bins(2)*Nbr_bins(3)) )  !�T�Ӥ�Vbin�ۭ�=cell�`�ƥ�

!���ت��� limit�N��cell��limit �`�@���X��segment ����segment�O�ݩ�o��cell�̭�
do i=1, ubound(cells,1)
    cells(i)%ghost_limits(1)= 0 !�k�s
    cells(i)%ghost_limits(2)= 0 !�k�s
end do

do i=1, Nbr_bins(1)
    do j=1, Nbr_bins(2)
        do k=1, Nbr_bins(3)
            ind= (k-1)*Nbr_bins(1)*Nbr_bins(2) + (j-1)*Nbr_bins(1) + i  !�μƦr�N��i j k
            cells(ind)%indx(1)= i   !ind��ܲĴX��cell ���K�si j k
            cells(ind)%indx(2)= j
            cells(ind)%indx(3)= k
        end do
    end do
end do


cells( ghost_segment_index(1) )%ghost_limits(1)= 1

do j= 2, ubound(ghost_segment_index,1)
    if( ghost_segment_index(j) .ne. ghost_segment_index(j-1) ) then
        cells( ghost_segment_index(j)   )%ghost_limits(1)= j        !index��ܲĤ@�Ө�̫�@�Ӹ��b���@��cell
        cells( ghost_segment_index(j-1) )%ghost_limits(2)= j-1
    end if
end do

!print *,Nbr_bins(1), Nbr_bins(2), Nbr_bins(3), ubound(cells,1)
!do i=1, Nbr_bins(1)
!    do j=1, Nbr_bins(2)
!        do k=1, Nbr_bins(3)
!            ind= (k-1)*Nbr_bins(1)*Nbr_bins(2) + (j-1)*Nbr_bins(1) + i
!            if( cells(ind)%ghost_limits(1) .ne. 0 ) then
!print *, j,k , cells(ind)%ghost_limits(1), cells(ind)%ghost_limits(2)-cells(ind)%ghost_limits(1)+1
!            end if
!        end do
!    end do
!print *, i    
!pause
!end do
!pause

end subroutine fiber_regroup 

                           
subroutine fiber_regroup_ShiftCenterToOrigion( fibers, hinges, box_size ) !2018/08/05 add �Nfiber����q���߲��ܰ_�I

type(fiber), allocatable, dimension(:) :: fibers
type(rod)  , allocatable, dimension(:) :: hinges

real(8), dimension(3)                  :: coord, min_coor, max_coor, box_size
integer(8)                             :: i, j, k, m  

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

!======================================================================
subroutine fiber_regroup_minmax_hinges( fibers, hinges ) !2018/08/05 add

type(fiber), allocatable, dimension(:) :: fibers
type(rod)  , allocatable, dimension(:) :: hinges
real(8), dimension(3)                  :: coord, min_coor, max_coor
integer(8)                             :: i, j, k, m

!2018/08/05  Compute center of mass for all fibers coord
!2018/08/05  Compute center of mass for all fibers coord
print *,"@@@"
print *,"@@@ based on all hinges"
write(301,*), "@@@"
write(301,*), "@@@ based on all hinges"

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
write(301,*),"cen ", coord

min_coor=   huge(0d0)
max_coor=  -huge(0d0)
do k=1,3
do i= 1, ubound (fibers,1)  
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
       coord= hinges(j)%X_i
       min_coor(k)= min( min_coor(k), coord(k) )
       min_coor(k)= min( min_coor(k), coord(k) )
       max_coor(k)= max( max_coor(k), coord(k) )
       max_coor(k)= max( max_coor(k), coord(k) )
    end do
end do  
end do      !2018/08/05 ��X�i�H�e�ǩҦ���hinges�ج[���y�нd��

print *,"min ", min_coor
print *,"max ", max_coor
write(301,*),"min ", min_coor
write(301,*),"max ", max_coor

min_coor=   huge(0d0)
max_coor=  -huge(0d0)
do k=1,3
do i=1, ubound (fibers,1)
   m= 0
   coord= 0
  !Compute center of mass for fiber i
   do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
		 coord= coord + hinges(j)%X_i
		 m= m+1
   end do
   coord= coord/real(m)
   min_coor(k)= min( min_coor(k), coord(k) )
   min_coor(k)= min( min_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )   
end do
end do    !2018/08/05 ��X�i�H�e�ǩҦ��� segments �ج[���y�нd��

print *,"min ", min_coor
print *,"max ", max_coor
write(301,*),"min ", min_coor
write(301,*),"max ", max_coor
write(301,*),"@@@"
print *,"@@@"


end subroutine  fiber_regroup_minmax_hinges

!======================================================================
subroutine fiber_regroup_minmax_segments( fibers, hinges ) !2018/08/05 add

type(fiber), allocatable, dimension(:) :: fibers
type(rod)  , allocatable, dimension(:) :: hinges
real(8), dimension(3)                  :: coord, min_coor, max_coor
integer(8)                             :: i, j, k, m

!2018/08/05  Compute center of mass for all fibers coord
print *,"@@@"
print *,"@@@ based on all segments"
write(301,*),"@@@"
write(301,*),"@@@ based on all segments"

m= 0
coord= 0
!!**�Phinge���P
do i= 1, ubound (fibers,1)  
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
       coord= coord + (hinges(j)%X_i+hinges(j+1)%X_i)/2d0
       m= m+1
    end do
end do        
coord= coord/real(m)  !2018/08/05  coord= center of mass 
!!**

print *,"cen ", coord
write(301,*),"cen ", coord

min_coor=   huge(0d0)
max_coor=  -huge(0d0)
do k=1,3
do i= 1, ubound (fibers,1)  
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !!�Phinge���P ����hinge-1 segment-2
       coord= (hinges(j)%X_i+hinges(j+1)%X_i)/2d0   !!�Phinge���P
       min_coor(k)= min( min_coor(k), coord(k) )
       min_coor(k)= min( min_coor(k), coord(k) )
       max_coor(k)= max( max_coor(k), coord(k) )
       max_coor(k)= max( max_coor(k), coord(k) )
    end do
end do  
end do      !2018/08/05 ��X�i�H�e�ǩҦ���hinges�ج[���y�нd��

print *,"min ", min_coor
print *,"max ", max_coor
write(301,*),"min ", min_coor
write(301,*),"max ", max_coor

min_coor=   huge(0d0)
max_coor=  -huge(0d0)
do k=1,3
do i= 1, ubound (fibers,1)
   m= 0
   coord= 0
   do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
      coord= coord + (hinges(j)%X_i+hinges(j+1)%X_i)/2d0
      m= m+1
   end do
   coord= coord/real(m)
   min_coor(k)= min( min_coor(k), coord(k) )    ! ����min�Mmax�n���ƨ⦸??
   min_coor(k)= min( min_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
end do  
end do      !2018/08/05 ��X�i�H�e�ǩҦ���hinges�ج[���y�нd��

print *,"min ", min_coor
print *,"max ", max_coor
write(301,*),"min ", min_coor
write(301,*),"max ", max_coor
write(301,*),"@@@"
!pause

end subroutine  fiber_regroup_minmax_segments
                           
                           
end module m_FiberRegroup   !2018/07/21  change name

!====================================================================