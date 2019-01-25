!======================================================================
module m_OutputData   !2018/07/21  change name
    
use m_DataStructures
implicit none
contains

subroutine output_data( fibers, hinges, simParameters )  !2018/10/10 修正

type(simulationParameters)               :: simParameters
type (fiber), dimension(:), allocatable  :: fibers
type(rod),    dimension(:), allocatable  :: hinges
integer(8)                               :: i, j, k

  open(4,file='OUTPUT/nbr_frames.txt')
  	  !write (4,*), ubound(hinges,1)
	   write (4,*), simParameters%frame
  close(4)

  k=1
  
 !write(3,*), simParameters%frame  
  write(3,*), ubound(fibers,1)
  
  do i=1, ubound(fibers,1)
  	 write(3,*), fibers(i)%nbr_hinges
  	 do j=k, k+fibers(i)%nbr_hinges-1
  		write(3,*), hinges(k)%X_i(1), hinges(k)%X_i(2), hinges(k)%X_i(3)
        if( simParameters%printVelocities ) then 
            write(5,*), hinges(k)%v_i(1), hinges(k)%v_i(2), hinges(k)%v_i(3)
        end if
		k=k+1
  	 end do
  end do
            
end subroutine output_data


subroutine output_Length( t, fibers, hinges )  !2018/09/22  修正

type(fiber), dimension(:),  allocatable :: fibers
type(rod),   dimension(:),  allocatable :: hinges

real(8), dimension(3)       :: coord
real(8)                     :: FiberLength_Total, FiberLength, SegmentLength, t
real                        :: length, lengthAvg
integer(8)                  :: i, j, k, mSegments
integer                     :: mm, nn, tt


    FiberLength_Total= 0.d0
    mSegments= 0

    do i= 1, ubound (fibers,1)
    
       FiberLength= 0.d0
       do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
          coord= hinges(j+1)%X_i-hinges(j)%X_i
          SegmentLength= sqrt( dot_product(coord,coord) )
          FiberLength= FiberLength + SegmentLength
          mSegments= mSegments + 1      
       end do
   
       FiberLength_Total= FiberLength_Total+FiberLength
!      print *,i, FiberLength
       
    end do

    tt= t*1.e6 + 0.5                              !2018/08/11   化成整數,單位=micro seconds
    mm= mSegments                                 !2018/08/11   total segments number
    nn= ubound(fibers,1)                          !2018/08/11   total Fiber number
    length= 1000*FiberLength_Total                !2018/09/22   total Fiber length(mm)
    lengthAvg= 1000*FiberLength_Total/nn          !2018/09/22   mean  Fiber length(mm)

    write(300,*), tt, nn, mm, length, lengthAvg
  
    print *,      "@@@", tt, nn, mm, length, lengthAvg
    write(301,*), "@@@", tt, nn, mm, length, lengthAvg
!   pause

end subroutine output_Length






!======================================================================
subroutine output_LengthDistribution( t, fibers, indexA )  !2018/08/12 新增

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
   j= fibers(i)%nbr_hinges-1  !2018/08/12  segment number of fiber(i)
   indexA(j)= indexA(j) + 1
end do

tt= t*1.0e6 + 0.5   !2018/08/12  +0.5 的用意是4捨5入
maxSegments= ubound(indexA,1)

!print *,"&&&", tt, ubound(fibers,1), maxSegments
write(302,*),  tt, ubound(fibers,1), maxSegments
do j=1, maxSegments
!  print *, "&&&( ", j, indexA(j)
  write(302,*),     j, indexA(j) 
end do
!pause

end subroutine output_LengthDistribution

subroutine output_OrientationTensor( t, fibers, hinges, AA )  !2018/08/12 新增
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
      AA(i,j)= AA(i,j) + ra(i)*ra(j)*length   !2018/08/14  加入權重長度,因直徑一樣,長度代表體積
   end do
   end do
   mm= mm + 1 
end do
end do

tt= t*1.0e6 + 0.5                    !2018/08/12  +0.5 的用意是4捨5入
trace_A= AA(1,1) + AA(2,2) + AA(3,3)

AA= AA/trace_A

!print *, "### trace ",trace_A, mm
print *, "###", tt, ubound (fibers,1)
print *, AA

write(301,*), "###", tt, ubound (fibers,1)
write(301,*), AA

write(303,*), tt, ubound (fibers,1)
write(303,*), AA
!pause

end subroutine output_OrientationTensor

subroutine output_PositionsForTheMomemt ( fibers, hinges, nbr_hinges)   !2018/08/31
type(fiber), dimension(:),   allocatable :: fibers
type(rod),   dimension(:),   allocatable :: hinges
integer(8)                               :: iii, jjj, kkk, nbr_hinges

        open(304,file='OUTPUT/PositionsForTheMoment.txt')  !2018/09/02
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

subroutine output_FiberLengthModification( fibers, hinges )  !2018/09/22  add
type(fiber) , dimension(:)            :: fibers
type (rod), dimension(:)              :: hinges
real(8)                               :: length, FiberLength, ScaleFactor
integer(8)                            :: i,j,k,l,jp1
real(8), dimension(3)                 :: box_size, coord
logical                               :: periodic_boundary


do i=1, ubound (fibers,1)   !2018/09/22  add
    
    FiberLength= 0.d0
	do j=fibers(i)%first_hinge, (fibers(i)%first_hinge+fibers(i)%nbr_hinges-2)
        jp1= j+1        
		coord= hinges(jp1)%X_i-hinges(j)%X_i
        length= sqrt( dot_product(coord, coord) )
        
        FiberLength= FiberLength + length
		hinges(j)%length2= length
    end do
    fibers(i)%Length= FiberLength
    
end do

!2018/09/22 將每一根Fiber長度修正為5mm
ScaleFactor= 1.0d0
FiberLength= 5.0d-3  !2018/09/22 將每一根Fiber長度修正為5mm

do i=1, ubound (fibers,1)  !2018/09/22  add
    ScaleFactor= FiberLength/fibers(i)%Length
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        jp1= j+1
        coord= hinges(jp1)%X_i - hinges(j)%X_i
        length= sqrt( dot_product(coord, coord) )
        
        coord= coord/length
        hinges(jp1)%X_i= hinges(j)%X_i + coord*hinges(j)%length2*ScaleFactor
    end do
end do

end subroutine output_FiberLengthModification


subroutine output_Initial_Positions_New( fibers, hinges, simParameters )   !2018/09/22
        
implicit none

type(simulationParameters)             :: simParameters
type(fiber), dimension(:), allocatable :: fibers
type(rod),   dimension(:), allocatable :: hinges

real(8), dimension(3) :: box_size, coord
real(8)               :: del_X, del_Y, del_Z, x, y, z

integer               :: mm1, mm2, mm3, mm1_A, mm1_B, mm2_A, mm2_B, mm3_A, mm3_B, ix, iy, iz
integer               :: i, j, N_Fiber, N_hinge

box_size= simParameters%box_size

open(305,file='OUTPUT/Initial_Positions_New.txt')  !2018/09/22
        
     mm1_A= -1
     mm1_B=  1
     
     mm2_A=  0
     mm2_B=  1
     
     mm3_A= -1
     mm3_B=  1
     
mm1_A=  0
mm1_B=  0

mm2_A=  0
mm2_B=  0
     
mm3_A=  0
mm3_B=  0     
     
     mm1= 1 + mm1_B - mm1_A
     mm2= 1 + mm2_B - mm2_A
     mm3= 1 + mm3_B - mm3_A
     
     N_Fiber= ubound(fibers,1)
     
     N_Fiber= N_Fiber*mm1*mm2*mm3

     print *, mm1_A, mm1_B, mm1
     print *, mm2_A, mm2_B, mm2
     print *, mm3_A, mm3_B, mm3
     print *,"box_size(1)= ",box_size(1)*mm1
     print *,"box_size(2)= ",box_size(2)*mm2
     print *,"box_size(3)= ",box_size(3)*mm3    
     print *,N_Fiber
     
     write(301,*), mm1_A, mm1_B, mm1
     write(301,*), mm2_A, mm2_B, mm2
     write(301,*), mm3_A, mm3_B, mm3
     write(301,*), "box_size(1)= ",box_size(1)*mm1
     write(301,*), "box_size(2)= ",box_size(2)*mm2
     write(301,*), "box_size(3)= ",box_size(3)*mm3    
     write(301,*), N_Fiber
     
     
     write (305,*), N_Fiber

     do i=1, ubound(fibers,1)

        do ix= mm1_A, mm1_B
        do iy= mm2_A, mm2_B
        do iz= mm3_A, mm3_B
            
            del_X= ix*box_size(1)
            del_Y= ix*box_size(2)
            del_Z= iz*box_size(3)
              
            N_hinge= fibers(i)%nbr_hinges
!           print *,N_hinge
            write (305,*), N_hinge
 
           do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
            
              x= hinges(j)%X_i(1) + del_X
              y= hinges(j)%X_i(2) + del_Y
              z= hinges(j)%X_i(3) + del_Z
 
!             print *, " 0 ", real(x,4), real(y,4), real(z,4)
!             write (305,*)," 0 ", real(x,4), real(y,4), real(z,4)
              write (305,*)," 0 ", x, y, z
           end do
           
        end do
        end do
        end do
!pause
     end do    

close( 305 )
!pause     

end subroutine output_Initial_Positions_New

subroutine output_DynamicP_1848( flowcase_1848, simParameters )               !2018/10/11 

implicit none
type(simulationParameters)                 :: simParameters
type(DynamicP), dimension(:), allocatable  :: flowcase_1848                   !2018/10/11 add
integer(8)                                 :: h                               !2018/10/11 
       
        h= simParameters%h                                                    !2018/10/11 增加
        
        print *,     "###"               
        print *,     "### Dynamic Parameters CHANGE"                          !2018/10/11 增加
        print *,     "### StepNo., Time(Micro.Sec.), Gamma_dot, Viscosity"    !2018/10/11 增加
        print *,     "###",&
                      int(simParameters%h),&
                      int(0.5+1.0e6*flowcase_1848(h)%Duration),&
                      real(flowcase_1848(h)%Shearrate),&
                      real(flowcase_1848(h)%Viscosity)                        !2018/10/11 增加
        print *,     "###"                                                    !2018/10/11 增加
        write(301,*),"###"
        write(301,*),"### Dynamic Parameters CHANGE"                          !2018/10/11 增加
        write(301,*),"### StepNo., Time(Micro.Sec.), Gamma_dot, Viscosity"    !2018/10/11 增加
        write(301,*),"###",&
                      int(simParameters%h),&
                      int(0.5+1.0e6*flowcase_1848(h)%Duration),&
                      real(flowcase_1848(h)%Shearrate),&
                      real(flowcase_1848(h)%Viscosity)                        !2018/10/11 增加
        write(301,*),"###"                                                   !2018/10/11 增加

end subroutine output_DynamicP_1848  !2018/10/11 change name



end module m_OutputData
!======================================================================