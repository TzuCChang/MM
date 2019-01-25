!====================================================================
module m_FindNeighbors  !2018/07/21  change name
    
use m_DataStructures
use m_UtilityLib

use omp_lib
implicit none
contains

subroutine find_neighbors_new( fibers, hinges, ghost_segments, neighbor_list, cells, simParameters )  !2018/11/25 change

type(simulationParameters)              :: simParameters
type(cell),  dimension(:), allocatable  :: cells
type(fiber), dimension(:)               :: fibers
type(rod),   dimension(:)               :: hinges
type(rod)                               :: ghost_hingeA, ghost_hingeB
type(segment), dimension(:)             :: ghost_segments
logical                                 :: flag

integer(8), dimension(:,:), allocatable:: neighbor_list
integer(8), dimension(:),   allocatable:: IndexSegms
integer(8), dimension(3)               :: i_cell, nbr_bins, indx
integer(8)                             :: i, j, k, m, o, nbr_neighbors, ii , jj, kk, indd, iSegm, mSegm  !error 修正2018/07/14

real(8), dimension(:,:), allocatable   :: distance_neighbors
real(8), dimension(:),   allocatable   :: distanceSegms
real(8), dimension(3)                  :: rad, d, Gab, A1, A2, B1, B2
real(8)                                :: s, t, Gab_min, r_fiber, threshold, epsilon
real(8)                                :: distanceFactor, ex_vol_const, fric_coeff
integer(8), dimension(10)  ::nC


nbr_bins       = simParameters%nbr_bins              !2018/10/10
r_fiber        = simParameters%r_fiber               !2018/10/10
distanceFactor = simParameters%distanceFactor        !2018/10/10
nbr_neighbors  = simParameters%nbr_neighbors         !2018/11/25

epsilon= 3*tiny(1d0)

threshold= distanceFactor*r_fiber

ex_vol_const= 1  !Not really necessary, the ex_vol_forces_moments_segs
                 !is used here to find the distances without being concerned with the forces...

if( allocated(  neighbor_list ) ) then
	deallocate( neighbor_list )
end if

allocate( neighbor_list(      ubound(hinges,1), nbr_neighbors ) )
allocate( distance_neighbors( ubound(hinges,1), nbr_neighbors ) )

allocate( distanceSegms( ubound(hinges,1) ) )
allocate( IndexSegms(    ubound(hinges,1) ) )

neighbor_list=      0
distance_neighbors= 1d10

distanceSegms=   1d10
IndexSegms=      0

nC= 0

do i=1, ubound( fibers, 1 )
	do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        
	    i_cell= cells(hinges(j)%ind)%indx
        A1= hinges(j  )%X_i     !2018/08/02 修正改寫
        A2= hinges(j+1)%X_i     !2018/08/02 修正改寫
        
        distanceSegms=   1d10
        IndexSegms=      0
        iSegm= 0
        
	    do ii= -1,1
	        indx(1)= i_cell(1) + ii
	        do jj= -1,1
	            indx(2)= i_cell(2) + jj
	            do kk= -1,1
	               indx(3)= i_cell(3) + kk
	               indd= (indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1) 
	               do k= cells(indd)%ghost_limits(1), cells(indd)%ghost_limits(2)
                       
	                    if( cells(indd)%ghost_limits(1).ne.0 ) then 

                            B1= ghost_segments(k)%A !2018/08/02 修正改寫
                            B2= ghost_segments(k)%B !2018/08/02 修正改寫

	                        if( are_boxes_intersect( A1, A2, B1, B2, r_fiber, 6*r_fiber ) ) then
		                        if( ( abs(dot_product(B1-A1,B1-A1)) .ge. epsilon ) .and.&  
		                            ( abs(dot_product(B1-A2,B1-A2)) .ge. epsilon ) .and.&
		                            ( abs(dot_product(B2-A1,B2-A1)) .ge. epsilon ) .and.&
		                            ( abs(dot_product(B2-A2,B2-A2)) .ge. epsilon ) ) then
                                    
			                        ghost_hingeA%X_i= B1
			                        ghost_hingeB%X_i= B2
                                    
					                call dist_segments( B1, B2, A1, A2, s, t, Gab, Gab_min ) !2018/08/02 修正改寫

   if(      Gab_min .lt. 1.d0*r_fiber ) then
            nC(1)= nC(1) + 1
   else if( Gab_min .lt. 2.d0*r_fiber ) then
            nC(2)= nC(2) + 1
   else if( Gab_min .lt. 3.d0*r_fiber ) then
            nC(3)= nC(3) + 1
   else if( Gab_min .lt. 4.d0*r_fiber ) then
            nC(4)= nC(4) + 1
   else if( Gab_min .lt. 5.d0*r_fiber ) then
            nC(5)= nC(5) + 1
   else 
            nC(6)= nC(6) + 1
   end if
   
                                    if( Gab_min .lt. Threshold ) then
                                        iSegm= iSegm + 1
                                        distanceSegms( iSegm )= Gab_min
                                        IndexSegms(    iSegm )= k
                                        !print *, iSegm, k, Gab_min                                     
                                    end if
	                            end if
	                        end if
                        end if
                        
                    end do	  
	            end do
	        end do
        end do
        
mSegm= iSegm 
       
call SortOrder( distanceSegms, IndexSegms, mSegm )

if( mSegm > nbr_neighbors )  mSegm= nbr_neighbors

do iSegm= 1, mSegm
   distance_neighbors(j,iSegm)= distanceSegms(iSegm)
   neighbor_list(j,iSegm)= IndexSegms(iSegm)  
   !print *, iSegm, IndexSegms(iSegm), distanceSegms(iSegm)
end do

!print *,"j,mSegm= ",j, mSegm        
!pause

    end do
end do

       write(*,100),  "###(", nC(1), nC(2), nC(3), nC(4), nC(5), nC(6)  !2018/11/19 add
100    format( A4, 6I9 )
       write(301,101),"###(", nC(1), nC(2), nC(3), nC(4), nC(5), nC(6)  !2018/11/19 add
101    format( A4, 6I9 )
       write(308,101),"###(", nC(1), nC(2), nC(3), nC(4), nC(5), nC(6) 
102    format( A4, 6I9 )
pause

deallocate( distance_neighbors )                         !2018/11/25 delete

end subroutine find_neighbors_new   
                 
                                                  
subroutine find_neighbors( fibers, hinges, ghost_segments, neighbor_list, cells, simParameters )   !2018/11/25 change

type(simulationParameters)             :: simParameters
type(cell),  dimension(:), allocatable :: cells
type(fiber), dimension(:)              :: fibers
type(rod),   dimension(:)              :: hinges
type(rod)                              :: ghost_hingeA, ghost_hingeB
type(segment), dimension(:)            :: ghost_segments
logical                                :: flag

integer(8), dimension(:,:), allocatable:: neighbor_list
integer(8), dimension(3)               :: i_cell, nbr_bins, indx
integer(8)                             :: i, j, k, m, o, nbr_neighbors, ii , jj, kk, indd, iSegm, mSegm  !error 修正2018/07/14

real(8), dimension(:,:), allocatable   :: distance_neighbors
real(8), dimension(3)                  :: rad, d, Gab, A1, A2, B1, B2
real(8)                                :: s, t, Gab_min, r_fiber, threshold, epsilon, times
real(8)                                :: distanceFactor, ex_vol_const, fric_coeff
integer(8), dimension(10)  ::nC

times=         simParameters%time                                !2018/11/19 add
nbr_bins       = simParameters%nbr_bins
r_fiber        = simParameters%r_fiber                               
distanceFactor = simParameters%distanceFactor
nbr_neighbors  = simParameters%nbr_neighbors
                               
epsilon= 3*tiny(1d0)

threshold= distanceFactor*r_fiber

ex_vol_const= 1  !Not really necessary, the ex_vol_forces_moments_segs
                 !is used here to find the distances without being concerned with the forces...

if( allocated(  neighbor_list ) ) then
	deallocate( neighbor_list )
end if

allocate( neighbor_list(      ubound(hinges,1), nbr_neighbors ) )
allocate( distance_neighbors( ubound(hinges,1), nbr_neighbors ) )

neighbor_list=      0
distance_neighbors= 1e10

nC= 0


do i=1, ubound( fibers, 1 )
	do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        
	    i_cell= cells(hinges(j)%ind)%indx
        A1= hinges(j  )%X_i     !2018/08/02 修正改寫
        A2= hinges(j+1)%X_i     !2018/08/02 修正改寫
        
        iSegm= 0
        
	    do ii= -1,1
	        indx(1)= i_cell(1) + ii
	        do jj= -1,1
	            indx(2)= i_cell(2) + jj
	            do kk= -1,1
	               indx(3)= i_cell(3) + kk
	               indd= (indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1) 
	               do k= cells(indd)%ghost_limits(1), cells(indd)%ghost_limits(2)
                       
	                    if( cells(indd)%ghost_limits(1).ne.0 ) then 

                            B1= ghost_segments(k)%A !2018/08/02 修正改寫
                            B2= ghost_segments(k)%B !2018/08/02 修正改寫

	                        if( are_boxes_intersect( A1, A2, B1, B2, r_fiber, 6*r_fiber ) ) then
		                        if( ( abs(dot_product(B1-A1,B1-A1)) .ge. epsilon ) .and.&  
		                            ( abs(dot_product(B1-A2,B1-A2)) .ge. epsilon ) .and.&
		                            ( abs(dot_product(B2-A1,B2-A1)) .ge. epsilon ) .and.&
		                            ( abs(dot_product(B2-A2,B2-A2)) .ge. epsilon ) ) then
                                    
			                        ghost_hingeA%X_i= B1
			                        ghost_hingeB%X_i= B2
                                    
					                call dist_segments( B1, B2, A1, A2, s, t, Gab, Gab_min ) !2018/08/02 修正改寫

   if(      Gab_min .lt. 1.d0*r_fiber ) then
            nC(1)= nC(1) + 1
   else if( Gab_min .lt. 2.d0*r_fiber ) then
            nC(2)= nC(2) + 1
   else if( Gab_min .lt. 3.d0*r_fiber ) then
            nC(3)= nC(3) + 1
   else if( Gab_min .lt. 4.d0*r_fiber ) then
            nC(4)= nC(4) + 1
   else if( Gab_min .lt. 5.d0*r_fiber ) then
            nC(5)= nC(5) + 1
   else 
            nC(6)= nC(6) + 1
   end if

                                    iSegm= iSegm + 1
                                    flag= .false.
                                    m= 1
                                    
!print *, iSegm, k, Gab_min, Threshold
		               
				                    do while( (m.le.nbr_neighbors) .and. (flag.eqv. .false.) )
                                        
					                    if( (Gab_min .lt. Threshold) .and.&
                                            (Gab_min .lt. distance_neighbors(j,m)) )then
                                            
				                            if( m .ne. nbr_neighbors ) then
				  		                        distance_neighbors(j,m+1:nbr_neighbors)=&
							                    distance_neighbors(j,m:nbr_neighbors-1)
						                        neighbor_list(j,m+1:nbr_neighbors)=&
					                            neighbor_list(j,m:nbr_neighbors-1)
                                            end if
                                            
						                    flag= .true.
						                    distance_neighbors(j,m)= Gab_min
				                            neighbor_list(j,m)= k
                                            
                                        end if         
				                        m=m+1
                                    end do                    
	                            end if
	                        end if
                        end if
                        
                    end do	  
	            end do
	        end do
        end do

!do iSegm= 1, nbr_neighbors
!   print *, iSegm, neighbor_list(j,iSegm), distance_neighbors(j,iSegm)
!end do
!pause 

    end do
end do

!       write(*,100),  "###(", nC(1), nC(2), nC(3), nC(4), nC(5), nC(6)  !2018/11/19 add
!100    format( A4, 6I9 )
!       write(301,101),"###(", nC(1), nC(2), nC(3), nC(4), nC(5), nC(6)  !2018/11/19 add
!101    format( A4, 6I9 )
       write(311,102),"###(", 1.d6*times, nC(1), nC(2), nC(3), nC(4), nC(5), nC(6) 
102    format( A4, F12.0, 6I9 )
!pause

deallocate( distance_neighbors )                         !2018/11/25 delete

end subroutine find_neighbors

subroutine SortOrder( dinstSegm, indexSegm, mSegm )

integer(8), dimension(:),  allocatable :: IndexSegm, indexA
integer(8)                             :: i1, i2, i3, k1, k2, k3, mSegm

real(8), dimension(:),   allocatable   :: dinstSegm
real(8)                                :: c1, c2, c3

allocate( indexA( mSegm ) )

    do i1= 1, mSegm 
       indexA(i1)= i1   
    end do
        
    do i1= 1, mSegm
        
       c1= dinstSegm(i1)
       k1= indexSegm(i1)
       
       i3= i1
       k3= k1
       c3= c1
       
       do i2= i1+1, mSegm
           
          k2= indexSegm(i2)
          c2= dinstSegm(i2)
          if( c3>c2 )  then
                  i3= i2
                  k3= k2
                  c3= c2
          end if
       end do
       
       i2= indexA(i3 )
       indexA(i3)= indexA(i1)
       indexA(i1)= i2
       
       k2= indexSegm(i3)
       indexSegm(i3)= indexSegm(i1)
       indexSegm(i1)= k2
           
       dinstSegm(i1)= c3      
       dinstSegm(i3)= c1       
       !print *,"B01",i1, i3, indexA(i1), c3
    end do
    !pause
    
    !do i1=1,mSegm
    !   print *,"B02",i1, indexSegm(i1), dinstSegm(i1)
    !end do    
    !pause

end subroutine  SortOrder
                               
end module m_FindNeighbors    !2018/07/21  change name

!********************************************************************