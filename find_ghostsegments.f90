!====================================================================
module m_FindGhostSegments   !2018/09/08  Add

use m_DataStructures
use m_UtilityLib

implicit none
contains


subroutine GhostSegments_Location( fibers, hinges, ghost_segments, box_size, box_dimension ) !2018/09/12 corrected

implicit none
type(fiber),   dimension(:), allocatable :: fibers
type(rod),     dimension(:), allocatable :: hinges
type(segment), dimension(:), allocatable :: ghost_segments

real(8), dimension(3)       :: box_size, coord
real(8)                     :: del_X, del_Z

integer                     :: mm1, mm2, mm3, ix, iz
integer, dimension(3)       :: box_dimension
integer(8)                  :: i, j, k, m, n


     mm1= box_dimension(1)
     mm2= box_dimension(2)
     mm3= box_dimension(3)

     m= 0

     do i=1, ubound(fibers,1)   
     do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        
        k= j+1
        n= 0
        do ix= -mm1, mm1
        do iz= -mm3, mm3
            
            n= n+1
            m= m+1
            
            ghost_segments(m)%ix= ix
            ghost_segments(m)%iz= iz
            
            del_X= ix*box_size(1)
            del_Z= iz*box_size(3)
            
            ghost_segments(m)%axis_loc= n

            ghost_segments(m)%orig_pos(1)= j
            ghost_segments(m)%orig_pos(2)= k
            
            ghost_segments(m)%A(1)= hinges(j)%X_i(1) + del_X
            ghost_segments(m)%A(2)= hinges(j)%X_i(2) 
            ghost_segments(m)%A(3)= hinges(j)%X_i(3) + del_Z
            
            ghost_segments(m)%B(1)= hinges(k)%X_i(1) + del_X
            ghost_segments(m)%B(2)= hinges(k)%X_i(2) 
            ghost_segments(m)%B(3)= hinges(k)%X_i(3) + del_Z
            
        end do
        end do
       
     end do
     end do

end subroutine  GhostSegments_Location

!====================================================================

subroutine GhostSegments_NewLocation( hinges, ghost_segments, box_size )  !2018/09/08  Add

type(rod),     dimension(:)              :: hinges
type(segment), dimension(:), allocatable :: ghost_segments

real(8)                                  :: del_X, del_Z
real(8),dimension(3)                     :: box_size

integer(8)                               :: i, j, k
integer                                  :: ix, iz

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO PRIVATE (i, j, k)

    
     do i= 1, ubound( ghost_segments, 1 )       !from segment til last segment

        j= ghost_segments(i)%orig_pos(1)        !beginning
        k= ghost_segments(i)%orig_pos(2)        !end
        
            ix= ghost_segments(i)%ix
            iz= ghost_segments(i)%iz
            
            del_X= ix*box_size(1)
            del_Z= iz*box_size(3)
            
            ghost_segments(i)%A(1)= hinges(j)%X_i(1) + del_X
            ghost_segments(i)%A(2)= hinges(j)%X_i(2) 
            ghost_segments(i)%A(3)= hinges(j)%X_i(3) + del_Z
            
            ghost_segments(i)%B(1)= hinges(k)%X_i(1) + del_X
            ghost_segments(i)%B(2)= hinges(k)%X_i(2) 
            ghost_segments(i)%B(3)= hinges(k)%X_i(3) + del_Z
     end do

!$OMP END DO !NOWAIT
!$OMP DO PRIVATE (i)

do i= 1, ubound(hinges,1)
	hinges(i)%T_Excl_Vol= 0d0
	hinges(i)%F_Excl_Vol= 0d0
end do

end subroutine GhostSegments_NewLocation 

!====================================================================

subroutine GhostSegments_Dimension( fibers, hinges, ghost_segments, box_size, box_dimension )  !2018/09/12 corrected
type(fiber),   dimension(:),  allocatable :: fibers
type(rod),     dimension(:),  allocatable :: hinges
type(segment), dimension(:), allocatable  :: ghost_segments

real(8),    dimension(3)    :: box_size, coord
real(8)                     :: FiberLength, maxLength


integer(8)                  :: i, j, nbr_segments, numClones, nbr_GhostSegments
integer,    dimension(3)    :: box_dimension
integer                     :: mm1, mm2, mm3


     maxLength=  -huge(0d0)
     do i= 1, ubound (fibers,1)
        FiberLength= 0
        do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
              coord= hinges(j+1)%X_i-hinges(j)%X_i
              FiberLength= FiberLength + sqrt( dot_product(coord,coord) )
        end do
        fibers(i)%Length= FiberLength
        maxLength= max( maxLength, FiberLength )
     end do

     mm1= 1 + int( 0.5 + maxLength/box_size(1) )
     mm2= 1 + int( 0.5 + maxLength/box_size(2) )
     mm3= 1 + int( 0.5 + maxLength/box_size(3) )
     
     box_dimension(1)= mm1
     box_dimension(2)= mm2
     box_dimension(3)= mm3

     nbr_segments= 0                                                !2018/09/09  Count the number of segments
     do i=1, ubound(fibers,1)
        nbr_segments= nbr_segments + fibers(i)%nbr_hinges - 1       !2018/09/08 
     end do
     
     numClones= (mm1+mm1+1)*(mm3+mm3+1)                             !2018/09/09
     
     nbr_GhostSegments= nbr_segments*numClones                      !2018/09/09   
     
     if( allocated(ghost_segments) ) deallocate( ghost_segments )
     
     allocate( ghost_segments(nbr_GhostSegments) ) 

     print *,"@@@ maxLength=", real(maxLength), real(box_size(1)), real(0.5 + maxLength/box_size(1)), mm1 
     print *,"@@@ maxLength=", real(maxLength), real(box_size(2)), real(0.5 + maxLength/box_size(2)), mm2
     print *,"@@@ maxLength=", real(maxLength), real(box_size(3)), real(0.5 + maxLength/box_size(3)), mm3
     print *,"@@@"
     print *,"@@@ nbr_GhostSegments=",int(nbr_segments), int(numClones), nbr_GhostSegments
     print *,"@@@"
     
     write(301,*), "@@@ maxLength=", real(maxLength), real(box_size(1)), real(0.5 + maxLength/box_size(1)), mm1 
     write(301,*), "@@@ maxLength=", real(maxLength), real(box_size(2)), real(0.5 + maxLength/box_size(2)), mm2
     write(301,*), "@@@ maxLength=", real(maxLength), real(box_size(3)), real(0.5 + maxLength/box_size(3)), mm3
     write(301,*), "@@@"
     write(301,*), "@@@ nbr_GhostSegments=",int(nbr_segments), int(numClones), nbr_GhostSegments
     write(301,*), "@@@"     
    !pause

end subroutine GhostSegments_Dimension



end module m_FindGhostSegments   !2018/07/21  change name

!====================================================================