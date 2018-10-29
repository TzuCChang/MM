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
real(8)                     :: FiberLength, FiberLengthTotal, EndToEndDistance, EndToEndDistanceTotal, t
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
  
  call output_Length( t, fibers, hinges, frame, printVelocities )  
            
end subroutine output_data


subroutine output_Length( t, fibers, hinges, frame, printVelocities )  !2018/07/14  ·s¼W
type (rod), allocatable, dimension(:)  :: hinges
type (fiber), allocatable, dimension(:):: fibers
integer(8)                             :: i, j, k, mm, nn, frame, nbr_Segments_Total, nbr_Fibers
logical                                :: printVelocities

type(rod)                   :: hinge1, hinge2, hingeHead, hingeEnd
real(8)                     :: FiberLength, SegmentLength, FiberLengthTotal, EndToEndDistance, EndToEndDistanceTotal, t
real(8),   dimension(3)     :: r

  
  k=1
  
  nbr_Segments_Total= 0
  FiberLengthTotal= 0.d0
  EndToEndDistanceTotal= 0.d0
  
  nbr_Fibers= ubound(fibers,1)
  
  do i=1, nbr_Fibers
     mm= fibers(i)%nbr_hinges
     nn= k+mm-1 
     
     nbr_Segments_Total= nbr_Segments_Total + mm -1
     
            hingeHead = hinges(k)
            hingeEnd  = hinges(nn)

!print *,"X=", hingeHead%X_i(1), hingeHead%X_i(2), hingeHead%X_i(3)
!print *,"X=", hingeEND%X_i(1), hingeEND%X_i(2), hingeEND%X_i(3)

            r= hingeHead%X_i-hingeEnd%X_i
            
            EndToEndDistance= sqrt(dot_product(r,r))
            
            r= r/EndToEndDistance
            
!print *,"L ", EndToEndDistance          
!print *,"r ", r(1), r(2), r(3)            
!print *,"@ "

     FiberLength= 0.d0
     
  	 do j=k, k+mm-2
!print *,"Y=", hinges(k)%X_i(1), hinges(k)%X_i(2), hinges(k)%X_i(3)
!print *,"Y=", hinges(k+1)%X_i(1), hinges(k+1)%X_i(2), hinges(k+1)%X_i(3)
        
        hinge1 = hinges(k)
        hinge2 = hinges(k+1)
        
        r= hinge2%X_i-hinge1%X_i
        
        SegmentLength= sqrt(dot_product(r,r))
        
        r= r/SegmentLength
        
        FiberLength= FiberLength + SegmentLength
        
!print *,"L ", SegmentLength, FiberLength
!print *,"r ", r(1), r(2), r(3)     
!print *," "

		k=k+1
     end do
     
!print *,"FiberLength ", FiberLength,EndToEndDistance
!pause
     
     FiberLengthTotal= FiberLengthTotal + FiberLength
     EndToEndDistanceTotal= EndToEndDistanceTotal + EndToEndDistance
     
  end do	
  
  print *, nbr_Fibers, nbr_Segments_Total
  print *, FiberLengthTotal, FiberLengthTotal/nbr_Fibers, EndToEndDistanceTotal/nbr_Fibers
  
  write(301,*), nbr_Fibers, nbr_Segments_Total
  write(301,*), FiberLengthTotal, FiberLengthTotal/nbr_Fibers, EndToEndDistanceTotal/nbr_Fibers
  
  write(300,*), t, ",", FiberLengthTotal, ",", FiberLengthTotal/nbr_Fibers !, ",", EndToEndDistanceTotal/nbr_Fibers
  !write(301,*), r
           
end subroutine output_Length


end module m_OutputData
