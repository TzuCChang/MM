!====================================================================
module m_FindNeighbors  !2018/07/21  change name
    
use m_DataStructures
use m_UtilityLib

use omp_lib
implicit none
contains

subroutine find_neighbors_new( fibers,&
                               hinges,&
                               ghost_segments,&
                               nbr_neighbors,&
                               neighbor_list,&
                               distance_neighbors,&
                               r_fiber,&
                               cells,&
                               Nbr_bins,&
                               distance_factor )

type(fiber), dimension(:)              :: fibers
type(rod), dimension(:)                :: hinges
type(rod)                              :: ghost_hingeA, ghost_hingeB
type(segment), dimension(:)            :: ghost_segments
integer(8)                             :: i, j, k, m, o, nbr_neighbors, ii , jj, kk, indd  !error ­×¥¿2018/07/14 
integer(8), dimension(:,:), allocatable:: neighbor_list
real(8), dimension(:,:), allocatable   :: distance_neighbors
real(8)                                :: r_fiber, threshold, s, t, Gab_min, ex_vol_const,fric_coeff
logical                                :: flag
real(8)                                :: epsilon, distance_factor
type(cell), dimension(:), allocatable  :: cells
integer, dimension(3)                  :: i_cell, nbr_bins
integer, dimension(3)                  :: indx
real(8), dimension(3)                  :: rad, d, Gab

real(8), dimension(3) ::A1, A2, B1, B2 


epsilon= 3*tiny(1d0)

threshold= distance_factor*r_fiber

ex_vol_const=1  !Not really necessary, the ex_vol_forces_moments_segs
                !is used here to find the distances without being concerned with the forces...

if( allocated( neighbor_list) ) then
	deallocate(neighbor_list)
end if

if( allocated( distance_neighbors) ) then
	deallocate(distance_neighbors)
end if

allocate( neighbor_list(     ubound(hinges,1),nbr_neighbors) )
allocate( distance_neighbors(ubound(hinges,1),nbr_neighbors) )

neighbor_list  	  = 0
distance_neighbors= 1e10

do i=1, ubound (fibers,1)
   do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2 
	  i_cell=cells(hinges(j)%ind)%indx
	  do ii=-1,1
	     indx(1)=i_cell(1)+ii
	     do jj=-1,1
	        indx(2)=i_cell(2)+jj
	        do kk=-1,1
	           indx(3)=i_cell(3)+kk      
	           indd=(indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
                   
	           do k= cells(indd)%ghost_limits(1), cells(indd)%ghost_limits(2)
	              if (cells(indd)%ghost_limits(1).ne.0) then  
                            
                      A1= hinges(j  )%X_i
                      A2= hinges(j+1)%X_i
                      B1= ghost_segments(k)%A
                      B2= ghost_segments(k)%B
                            
	                  if( are_boxes_intersect( A1, A2, B1, B2, r_fiber, 6*r_fiber )) then
		                  if( (abs(dot_product(B1-A1,B1-A1)) .ge. epsilon) .and.&  
		                      (abs(dot_product(B1-A2,B1-A2)) .ge. epsilon) .and.&
		                      (abs(dot_product(B2-A1,B2-A1)) .ge. epsilon) .and.&
		                      (abs(dot_product(B2-A2,B2-A2)) .ge. epsilon)) then
                                    
			                   ghost_hingeA%X_i= B1
			                   ghost_hingeB%X_i= B2
                                    
					           call dist_segments( ghost_hingeA%X_i,&
                                                   ghost_hingeB%X_i,&
	                                        	   hinges(j)%X_i,&
	                                        	   hinges(j+1)%X_i,&
                                                   s,&
                                                   t,&
                                                   Gab,&
		                                           Gab_min )                                        
					           flag=.false.
                                    
					           m=1
                                    
				               do while( (m.le.nbr_neighbors) .and. (flag.eqv. .false.) )
					               if ((Gab_min .lt. Threshold) .and.&
                                       (Gab_min .lt. distance_neighbors(j,m)))then
				                       if (m.ne. nbr_neighbors) then
				  		                   distance_neighbors(j,m+1:nbr_neighbors)=&
							               distance_neighbors(j,m:nbr_neighbors-1)
						                   neighbor_list(j,m+1:nbr_neighbors)=&
					                       neighbor_list(j,m:nbr_neighbors-1)
						               end if
						               flag=.true.
						               distance_neighbors(j,m)= Gab_min
				                       neighbor_list(j,m)=k
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
   end do
end do

end subroutine find_neighbors_new   

end module m_FindNeighbors    !2018/07/21  change name

!********************************************************************