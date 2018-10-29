!====================================================================
module m_FiberRegroup   !2018/07/21  change name

! subroutine: fiber_regroup
! subroutine: fiber_regroup_ShiftCenterToOrigion
! subroutine: fiber_regroup_minmax_hinges
! subroutine: fiber_regroup_minmax_segments

use m_DataStructures    !m 代表module(涵蓋這些subroutine) !Datastructure在tipo裡
use m_UtilityLib

implicit none   !implicit none代表有些東西不宣告
contains

!*===================================================================
 subroutine fiber_regroup( fibers,&    !用i表示               !2018/07/21 change name
                           hinges,&    !將hinge編號
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
type(rod),     dimension(:), allocatable :: hinges  !一維dimension
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
! 開始計算週期
! 如果y有週期:要算27個box(三度空間) / 若y沒有週期(牆壁):要算9個box(似井字)
if( simParameters%IsPeriodicY ) then    !%代表type底下的變數 而這個變數是logical
    numClones= 27    !2018/08/05 修正 change from 9 to 27
else
    numClones= 9     !2018/08/05 修正 change from 5 to 9
end if

nbr_segments= 0 !初始條件=0

! Count the number of segments
do i=1, ubound(fibers,1)    !upperbound: fibers跟1比 看誰大       
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !第一個hinge做到最後一個hinge >>因為兩個hinge=一個segment
        nbr_segments= nbr_segments+1
    end do
end do

! ALLOCATE : creates space for allocatable arrays and variables with the POINTER attribute. 
! DEALLOCATE : frees space previously allocated for allocatable arrays and pointer targets.
! 假如有宣告, 先取消宣告, 再重新宣告一次>>確保不受到其他影響  !因為ghost_segments會一直更新數值且被重複利用(若fiber斷掉, 需要一個新的描述)
if( allocated(ghost_segments) ) deallocate( ghost_segments )	
allocate( ghost_segments(nbr_segments*numClones) ) ! change 5 to 9 if you want periodicity in all directions

if( allocated(hinge_index) ) deallocate( hinge_index )
allocate( hinge_index(nbr_segments,2) )

if( allocated(ghost_segment_index) ) deallocate( ghost_segment_index )
allocate( ghost_segment_index(numClones*nbr_segments) ) ! change 5 to 9 if you want periodicity in all directions

!dX 先暫時不理(尚待確認)
dX= mod( gamma_dot*box_size(2)/2*t, box_size(1) )  ! This is how much the image boxes need to be translated for the lees edward boundaries
!shear rate*y=速度 速度/2=中間速度 中間速度*t=位移
m= 0    !m=所有segment的數目 (初始)

do i=1, ubound(fibers,1)   !從第一個fiber開始
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !找到所有的segment
        !每一根纖維的first hinge不一定都是1 所以要先找到fiber第一個hinge 去類推後面的
        k= j+1  !fiber(i)的hinge(端點)的數量
       
        !This is the original segment
        m= m + 1  !2018/08/05  changet  !m是segment 每一次m+1就是複製segment到空間的cell
        ghost_segments(m)%axis_loc= 1   !1代表中間, 2+x, 3-x, 4+z, 5_z, 6+x+z, 7-x+z, 8+x-z, 9-x-z
        ghost_segments(m)%orig_pos(1)= j    
        ghost_segments(m)%orig_pos(2)= k
        ghost_segments(m)%A= hinges(j)%X_i  !segment(rod)
        ghost_segments(m)%B= hinges(k)%X_i  !hinge為座標點(rod頭跟尾)

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
        
        ! 新增四個斜角(+x+z, +x-z, -x+z, -x-z)
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
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !找segment
        r= hinges(j+1)%X_i - hinges(j)%X_i  !r=segment的向量
        max_length= max( max_length, sqrt(dot_product(r,r)) )  !2018/08/02 找出所有的Segments的長度中最大值
    end do
end do

min_coor=   huge(0d0)   ! huge: returns the largest number that is not an infinity
max_coor=  -huge(0d0)
do k=1,3    !k=1 2 3 代表x y z座標  !這邊是用segment比(可能會造成比較的時間過長>>因為跟全部的segment比)
    do i=1, ubound(ghost_segments,1)    !找出頭尾端點(A B) 的最大值和最小值
        min_coor(k)= min( min_coor(k), ghost_segments(i)%A(k) ) !A是j(頭) 
        min_coor(k)= min( min_coor(k), ghost_segments(i)%B(k) ) !B是k(尾)
	    max_coor(k)= max( max_coor(k), ghost_segments(i)%A(k) )
	    max_coor(k)= max( max_coor(k), ghost_segments(i)%B(k) )
	end do  
end do      !2018/08/02 找出可以容納所有的segments框架的座標範圍

! 待測試 還需加box長度
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

bin_length= 1.5*( max_length + 2*r_fiber )  !先找出所有segment的範圍 再把把全部ghost在空間中的範圍做切割(之後只要在目的地附近找就可以了)bin_length約15個直徑左右
! 通常aspect ratio=5為max length 所以max length+直徑=約六倍直徑 *1.5=9 約等於10個直徑

!min_coor= min_coor- 2*r_fiber
!max_coor= max_coor+ 2*r_fiber

min_coor= min_coor- 2*r_fiber-max_length  !2018/08/02 再擴大框架的座標範圍
max_coor= max_coor+ 2*r_fiber+max_length    !往上往下各增加一個segment的長度

Nbr_bins= ceiling( (max_coor-min_coor)/bin_length )  !2018/08/05 CEILING ( 4.80) has the value  5
                                                     !2018/08/05 CEILING (-2.55) has the value -2
!用bin length去切 找出總共要切幾個bin ceiling為取最大的整數
!ceiling=最接近的右邊整數

!bin_length= (max_coor(3)-min_coor(3))/Nbr_bins(3)
!min_coor= min_coor-bin_length
!max_coor= min_coor+bin_length        !Mod 9/28/2014
!Nbr_bins= Nbr_bins+2

!重新取一次 晚點再看一次
!切盒子的目的>>為了之後比較能夠省空間 找neighbor只要找周圍的bin(但是原本寫法是全部找 所以無法節省時間)
min_coor= min_coor - 2*bin_length  !2018/08/02 再進一步擴大框架的座標範圍
max_coor= max_coor + 2*bin_length  !2018/08/02 再進一步擴大框架的座標範圍
Nbr_bins= Nbr_bins + 4             !2018/08/05 共增加4個bin_length的空間
                                   !2018/08/05 針對所有ghost_segments,找出包容的空間,
                                   !2018/08/05 並且依bin_length長度大小,將其長寬高三個方向,
                                   !2018/08/05 分別切割成許多正方體小盒子,
                                   !2018/08/05 數量為Nbr_bin(1)*Nbr_bin(2)*Nbr_bin(3) (1 2 3 都加4個)
!print *,"Nbr bins", Nbr_bins

do j=1,ubound(hinges,1)
   hinges(j)%indx= 0            !2018/08/05  先設定給起始值,歸零
end do

!此為真實的fiber                 !2018/08/02   floor( 3.412)=  3.00000  取整數
do i= 1, ubound(fibers,1)       !2018/08/02   floor(-3.412)= -4.00000  取整數
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
       center= ( hinges(j)%X_i + hinges(j+1)%X_i )/2                 !2018/08/02  change    !取兩個節點間(segment)的中間點
       do k=1,3                           
            indx(k)= floor( (center(k)-min_coor(k))/bin_length ) + 1 !2018/08/02  change    !中間點/bin length可以得知落在哪個bin裡面 !indx=???
            if( indx(k) .le. 3 )           indx(k)= 2    !若bin<3 則指定從第二個bin開始(不從1開始) >>為保守寫法
            if( indx(k) .ge. Nbr_bins(k) ) indx(k)= Nbr_bins(k)-1
	   end do
	   hinges(j)%ind= (indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
    end do       !2018/08/05 判斷hinges(j)~(j+1)的segments的質量中心點位置落在哪一個的空間小盒子內 ind>>用一維表示三維
end do

!此為複製後的ghost segment
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

call QsortC( ghost_segment_index , ghost_segments ) !搜尋 排序

if( allocated(cells) ) deallocate(cells)

allocate( cells(Nbr_bins(1)*Nbr_bins(2)*Nbr_bins(3)) )  !三個方向bin相乘=cell總數目

!此目的為 limit代表cell的limit 總共有幾個segment 哪些segment是屬於這個cell裡面
do i=1, ubound(cells,1)
    cells(i)%ghost_limits(1)= 0 !歸零
    cells(i)%ghost_limits(2)= 0 !歸零
end do

do i=1, Nbr_bins(1)
    do j=1, Nbr_bins(2)
        do k=1, Nbr_bins(3)
            ind= (k-1)*Nbr_bins(1)*Nbr_bins(2) + (j-1)*Nbr_bins(1) + i  !用數字代表i j k
            cells(ind)%indx(1)= i   !ind表示第幾個cell 順便存i j k
            cells(ind)%indx(2)= j
            cells(ind)%indx(3)= k
        end do
    end do
end do


cells( ghost_segment_index(1) )%ghost_limits(1)= 1

do j= 2, ubound(ghost_segment_index,1)
    if( ghost_segment_index(j) .ne. ghost_segment_index(j-1) ) then
        cells( ghost_segment_index(j)   )%ghost_limits(1)= j        !index表示第一個到最後一個落在哪一個cell
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

                           
subroutine fiber_regroup_ShiftCenterToOrigion( fibers, hinges, box_size ) !2018/08/05 add 將fiber的質量中心移至起點

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
end do      !2018/08/05 找出可以容納所有的hinges框架的座標範圍

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
end do    !2018/08/05 找出可以容納所有的 segments 框架的座標範圍

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
!!**與hinge不同
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
    do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2   !!與hinge不同 為何hinge-1 segment-2
       coord= (hinges(j)%X_i+hinges(j+1)%X_i)/2d0   !!與hinge不同
       min_coor(k)= min( min_coor(k), coord(k) )
       min_coor(k)= min( min_coor(k), coord(k) )
       max_coor(k)= max( max_coor(k), coord(k) )
       max_coor(k)= max( max_coor(k), coord(k) )
    end do
end do  
end do      !2018/08/05 找出可以容納所有的hinges框架的座標範圍

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
   min_coor(k)= min( min_coor(k), coord(k) )    ! 為何min和max要重複兩次??
   min_coor(k)= min( min_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
end do  
end do      !2018/08/05 找出可以容納所有的hinges框架的座標範圍

print *,"min ", min_coor
print *,"max ", max_coor
write(301,*),"min ", min_coor
write(301,*),"max ", max_coor
write(301,*),"@@@"
!pause

end subroutine  fiber_regroup_minmax_segments
                           
                           
end module m_FiberRegroup   !2018/07/21  change name

!====================================================================