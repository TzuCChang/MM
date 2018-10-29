!======================================================================
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


subroutine output_Length( t, fibers, hinges )               !2018/08/11 corrected
type(fiber), dimension(:),  allocatable :: fibers
type(rod),   dimension(:),  allocatable :: hinges
real(8)                     :: FiberLength, SegmentLength, t
real(8), dimension(3)       :: coord
integer(8)                  :: i, j, k, mSegments
real                        :: length, lengthAvg
integer                     :: mm, nn, tt

FiberLength= 0
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

tt= t*1.e6 + 0.5              !2018/08/11   make it integer,unit=micro seconds
mm= mSegments                 !2018/08/11   total segments number
nn= ubound(fibers,1)          !2018/08/11   total Fiber number
length= FiberLength           !2018/08/11   total Fiber length
lengthAvg= length/nn          !2018/08/11   mean  Fiber length

  write(300,*), tt, nn, mm, length, lengthAvg
  
  print *,      "@@@", tt, nn, mm, length, lengthAvg
  write(301,*), "@@@", tt, nn, mm, length, lengthAvg
  !pause

end subroutine output_Length

!======================================================================
subroutine output_LengthDistribution( t, fibers, indexA )  !2018/08/12 add

type(fiber), dimension(:), allocatable :: fibers
integer(8),  dimension(:), allocatable :: indexA
real(8)      :: t
integer(8)   :: i, j, maxSegments
integer      :: tt


if( allocated(indexA) .eq. .false. )  then
    
    maxSegments= 1
    do i= 1, ubound(fibers,1)  
       j= fibers(i)%nbr_hinges-1
       maxSegments= max( maxSegments,j )
    end do
    
    allocate( indexA(maxSegments) )
    
end if

indexA= 0
do i= 1, ubound (fibers,1)
   j= fibers(i)%nbr_hinges-1        !2018/08/12  segment number of fiber(i)
   indexA(j)= indexA(j) + 1
end do

tt= t*1.0e6 + 0.5                   !2018/08/12  +0.5 is for the rounding
maxSegments= ubound(indexA,1)

!print *,"&&&", tt, ubound(fibers,1), maxSegments
write(302,*),  tt, ubound(fibers,1), maxSegments
do j=1, maxSegments
!  print *, "&&&( ", j, indexA(j)
  write(302,*),     j, indexA(j) 
end do
!pause

end subroutine output_LengthDistribution

!======================================================================
subroutine output_OrientationTensor( t, fibers, hinges, AA )  !2018/08/12 add
type(fiber), dimension(:),   allocatable :: fibers
type(rod),   dimension(:),   allocatable :: hinges
real(8),     dimension(:,:)              :: AA
real(8),     dimension(3)                :: ra
real(8)      :: t, length, trace_A
integer(8)   :: ia, ja,  mm
integer      :: tt, i, j

AA= 0
mm= 0
do ia= 1, ubound (fibers,1)
do ja= fibers(ia)%first_hinge, fibers(ia)%first_hinge+fibers(ia)%nbr_hinges-2
   ra= hinges(ja+1)%X_i - hinges(ja)%X_i
       length= sqrt( dot_product(ra,ra) )
   ra= ra/length
   do i=1,3
   do j=1,3
      AA(i,j)= AA(i,j) + ra(i)*ra(j)*length     !2018/08/14  the length can represent the volume because the radius is the same
   end do
   end do
   mm= mm + 1 
end do
end do

tt= t*1.0e6 + 0.5                               !2018/08/12  +0.5 is for the rounding
trace_A= AA(1,1) + AA(2,2) + AA(3,3)

AA= AA/trace_A

print *, "### trace ",trace_A, mm
print *, "###", tt, ubound (fibers,1)
print *, AA

write(301,*), "###", tt, ubound (fibers,1)
write(301,*), AA

write(303,*), tt, ubound (fibers,1)
write(303,*), AA
!pause

end subroutine output_OrientationTensor

!======================================================================
subroutine output_PositionsForTheMomemt ( fibers, hinges, nbr_hinges)   !2018/08/31
type(fiber), dimension(:),   allocatable :: fibers
type(rod),   dimension(:),   allocatable :: hinges
integer(8)                               :: iii, jjj, kkk, nbr_hinges

        open(304,file='OUTPUT/PositionsForTheMoment.txt')               !2018/09/02
        kkk=1
        write (304,*), ubound(fibers,1)
        do iii=1, ubound(fibers,1)
  	        write (304,*), fibers(iii)%nbr_hinges
  	            do jjj=kkk, kkk+fibers(iii)%nbr_hinges-1
  		            write (304,*),0, real(hinges(kkk)%X_i(1),4), real(hinges(kkk)%X_i(2),4), real(hinges(kkk)%X_i(3),4)
		            kkk=kkk+1
  	        end do
        end do
        close(304)                                         !2018/09/02
!pause

end subroutine output_PositionsForTheMomemt

end module m_OutputData
!======================================================================