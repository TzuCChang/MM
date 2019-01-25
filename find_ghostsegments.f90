!====================================================================
module m_FindGhostSegments   !2018/09/08  Add

use m_DataStructures
use m_UtilityLib

implicit none
contains


subroutine GhostSegments_Location( fibers, hinges, ghost_segments, simParameters ) !2018/10/09  修正

implicit none
type(simulationParameters)                :: simParameters
type(fiber),   dimension(:), allocatable  :: fibers
type(rod),     dimension(:), allocatable  :: hinges
type(segment), dimension(:), allocatable  :: ghost_segments

real(8), dimension(3)       :: box_size, coord
real(8)                     :: del_X, del_Y, del_Z, dX, gamma_dot, time     !2018/10/02

integer(8)                  :: mm1, mm2, mm3, ix, iy, iz    !2018/10/02
integer(8)                  :: i, j, k, m, n


box_size  = simParameters%box_size
gamma_dot = simParameters%gamma_dot
time      = simParameters%time

     mm1= simParameters%box_dimension(1)
     mm2= simParameters%box_dimension(2)
     mm3= simParameters%box_dimension(3)

     dX = mod( gamma_dot*box_size(2)/2*time,box_size(1) )  ! 2018/10/02 This is how much the image boxes need to be translated for the lees edward boundaries
     
if( simParameters%IsPeriodicY .eqv. .false. ) then  ! 2018/10/02  wall 
    mm2= 0
    dX=  0.d0 
end if

     m= 0

     do i=1, ubound(fibers,1)   
     do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        
        k= j+1
        n= 0
        
        do ix= -mm1, mm1
        do iy= -mm2, mm2                                       !2018/10/03
        do iz= -mm3, mm3
            
            n= n+1
            m= m+1
            
            ghost_segments(m)%ix= ix
            ghost_segments(m)%iy= iy                           !2018/10/02 
            ghost_segments(m)%iz= iz

            del_X= ix*box_size(1) + iy*dX                      !2018/10/03
            del_Y= iy*box_size(2)                              !2018/10/05
            del_Z= iz*box_size(3)
       
            ghost_segments(m)%axis_loc= n

            ghost_segments(m)%orig_pos(1)= j
            ghost_segments(m)%orig_pos(2)= k
            
            ghost_segments(m)%A(1)= hinges(j)%X_i(1) + del_X               !2018/10/02
            ghost_segments(m)%A(2)= hinges(j)%X_i(2) + del_Y               !2018/10/02
            ghost_segments(m)%A(3)= hinges(j)%X_i(3) + del_Z
            
            ghost_segments(m)%B(1)= hinges(k)%X_i(1) + del_X               !2018/10/02
            ghost_segments(m)%B(2)= hinges(k)%X_i(2) + del_Y               !2018/10/02
            ghost_segments(m)%B(3)= hinges(k)%X_i(3) + del_Z
            
        end do
        end do
        end do

     end do
     end do

end subroutine  GhostSegments_Location

subroutine GhostSegments_NewLocation( hinges, ghost_segments, simParameters ) !2018/10/05  修正
type(rod),     dimension(:)               :: hinges
type(segment), dimension(:), allocatable  :: ghost_segments
type(simulationParameters)                :: simParameters

real(8)                 :: del_X, del_Y, del_Z, dX, gamma_dot, time
real(8),dimension(3)    :: box_size
integer(8)              :: i, j, k
integer(8)              :: ix, iy, iz

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO PRIVATE (i, j, k)

box_size  = simParameters%box_size
gamma_dot = simParameters%gamma_dot
time      = simParameters%time

     dX = mod( gamma_dot*box_size(2)/2*time,box_size(1) )  ! 2018/10/05 增加
     
     do i= 1, ubound( ghost_segments, 1 )

        j= ghost_segments(i)%orig_pos(1)
        k= ghost_segments(i)%orig_pos(2)
        
            ix= ghost_segments(i)%ix
            iy= ghost_segments(i)%iy            
            iz= ghost_segments(i)%iz
            
            del_X= ix*box_size(1) + iy*dX
            del_Y= iy*box_size(2)
            del_Z= iz*box_size(3)
            
            ghost_segments(i)%A(1)= hinges(j)%X_i(1) + del_X
            ghost_segments(i)%A(2)= hinges(j)%X_i(2) + del_Y
            ghost_segments(i)%A(3)= hinges(j)%X_i(3) + del_Z
            
            ghost_segments(i)%B(1)= hinges(k)%X_i(1) + del_X
            ghost_segments(i)%B(2)= hinges(k)%X_i(2) + del_Y
            ghost_segments(i)%B(3)= hinges(k)%X_i(3) + del_Z
     end do

!$OMP END DO !NOWAIT
!$OMP DO PRIVATE (i)


end subroutine GhostSegments_NewLocation 


subroutine GhostSegments_Dimension( fibers, hinges, ghost_segments, simParameters )  !2018/10/02 修正

type(fiber),   dimension(:),  allocatable :: fibers
type(rod),     dimension(:),  allocatable :: hinges
type(segment), dimension(:), allocatable  :: ghost_segments
type(simulationParameters)                :: simParameters

real(8),    dimension(3)    :: box_size, coord
real(8)                     :: FiberLength, maxLength
integer(8)                  :: i, j, nbr_segments, numClones, nbr_GhostSegments
integer(8)                  :: mm1, mm2, mm3


box_size=      simParameters%box_size

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

     mm1= 1 + int( 0.1d0 + 0.5d0*maxLength/box_size(1) )
     mm2= 1 + int( 0.1d0 + 0.5d0*maxLength/box_size(2) )
     mm3= 1 + int( 0.1d0 + 0.5d0*maxLength/box_size(3) )
 
     if( simParameters%IsPeriodicY .eqv. .false. ) then  ! 2018/10/02  wall 
         mm2= 0   
     end if

     simParameters%box_dimension(1)= mm1
     simParameters%box_dimension(2)= mm2
     simParameters%box_dimension(3)= mm3

     nbr_segments= 0          !2018/09/09  Count the number of segments
     do i=1, ubound(fibers,1)
        nbr_segments= nbr_segments + fibers(i)%nbr_hinges - 1 !2018/09/08 
     end do
     
     numClones= (mm1+mm1+1)*(mm2+mm2+1)*(mm3+mm3+1)           !2018/09/09
     
     nbr_GhostSegments= nbr_segments*numClones                !2018/09/09   
     
     if( allocated(ghost_segments) ) deallocate( ghost_segments )
     
     allocate( ghost_segments(nbr_GhostSegments) ) 

!2018/12/09
     write(*,100),"@@@ box_size(1)=",box_size(1)*1000,maxLength*1000,(0.1d0+0.5d0*maxLength/box_size(1)),mm1,2*mm1+1  !2018/11/22
100  format( A17,3F12.6,I6,I6)
     write(*,101),"@@@ box_size(2)=",box_size(2)*1000,maxLength*1000,(0.1d0+0.5d0*maxLength/box_size(2)),mm2,2*mm2+1   !2018/11/22
101  format( A17,3F12.6,I6,I6)
     write(*,102),"@@@ box_size(3)=",box_size(3)*1000,maxLength*1000,(0.1d0+0.5d0*maxLength/box_size(3)),mm3,2*mm3+1   !2018/11/22
102  format( A17,3F12.6,I6,I6)
     write(*,*),"@@@"
     write(*,*),"@@@ nbr_GhostSegments=",int(nbr_segments), int(numClones), nbr_GhostSegments
     write(*,*),"@@@"

     write(301,200),"@@@ box_size(1)=",box_size(1)*1000,maxLength*1000,(0.1d0+0.5d0*maxLength/box_size(1)),mm1,2*mm1+1   !2018/11/22
200  format( A17,3F12.6,I6,I6)
     write(301,201),"@@@ box_size(2)=",box_size(2)*1000,maxLength*1000,(0.1d0+0.5d0*maxLength/box_size(2)),mm2,2*mm2+1   !2018/11/22
201  format( A17,3F12.6,I6,I6)
     write(301,202),"@@@ box_size(3)=",box_size(3)*1000,maxLength*1000,(0.1d0+0.5d0*maxLength/box_size(3)),mm3,2*mm3+1   !2018/11/22
202  format( A17,3F12.6,I6,I6)
     write(301,*),"@@@"
     write(301,*),"@@@ nbr_GhostSegments=",int(nbr_segments), int(numClones), nbr_GhostSegments
     write(301,*),"@@@"
!pause

end subroutine GhostSegments_Dimension



end module m_FindGhostSegments   !2018/07/21  change name

!====================================================================