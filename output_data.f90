module m_OutputData   !2018/07/21  change name
    
use m_DataStructures
implicit none
contains

subroutine output_data( t, fibers, hinges, frame, printVelocities )
type (rod), allocatable, dimension(:)  :: hinges
type (fiber), allocatable, dimension(:):: fibers
integer(8)                             :: i, j, k, mm, nn, frame, nbr_Segments_Total, nbr_Fibers
logical                                :: printVelocities

type(rod)                   :: hinge1, hinge2, hingeHead, hingeEnd
real(8)                     :: FiberLength, FiberLength_Total, t
real(8),   dimension(3)     :: r


  open(4,file='OUTPUT/nbr_frames.txt')
  	!write (4,*), ubound(hinges,1)
	write (4,*), frame
  close(4)

  !write (3,*), frame
  k=1
  write (3,*), ubound(fibers,1)
  do i=1, ubound(fibers,1)
  	write (3,*), fibers(i)%nbr_hinges
  	do j=k, k+fibers(i)%nbr_hinges-1
  		write (3,*), hinges(k)%X_i(1), hinges(k)%X_i(2), hinges(k)%X_i(3)
        if (printVelocities) then
            write (5,*), hinges(k)%v_i(1), hinges(k)%v_i(2), hinges(k)%v_i(3)
        end if
        
		k=k+1
  	end do
  end do
            
end subroutine output_data


subroutine output_Length( t, fibers, hinges, frame, printVelocities )  !2018/08/11修正
type(rod),   dimension(:),  allocatable :: hinges
type(fiber), dimension(:),  allocatable :: fibers
logical                                 :: printVelocities
integer(8)                  :: frame

real(8)                     :: FiberLength, SegmentLength, t
real(8), dimension(3)       :: coord
integer(8)                  :: i, j, k, mSegments
real                        :: length, lengthAvg
integer                     :: mm, nn, tt

  
mSegments= 0
coord= 0
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
   coord= hinges(j+1)%X_i-hinges(j)%X_i
   FiberLength= FiberLength + sqrt( dot_product(coord,coord) )
   mSegments= mSegments + 1
end do
end do        
coord= coord/real(mSegments)  !2018/08/11 

tt= t*1.e6 + 0.5              !2018/08/11   化成整數,單位=micro seconds
mm= mSegments                 !2018/08/11   total segments number
nn= ubound(fibers,1)          !2018/08/11   total Fiber number
length= FiberLength           !2018/08/11   total Fiber length
lengthAvg= length/nn          !2018/08/11   mean  Fiber length

  print *, tt, nn, mm, length, lengthAvg
  write(300,*), tt, nn, mm, length, lengthAvg
  write(301,*), tt, nn, mm, length, lengthAvg
  
!pause

end subroutine output_Length


end module m_OutputData
