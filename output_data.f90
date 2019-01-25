!======================================================================
module m_OutputData   !2018/07/21  change name
use m_DataStructures
use omp_lib

implicit none
contains

!======================================================================

subroutine output_OpenFiles(  simParameters )  !2018/11/25
type(simulationParameters) :: simParameters

open(300,file='OUTPUT/meanLength.txt')
open(301,file='OUTPUT/OutputMessage.txt')
open(302,file='OUTPUT/FiberLengthDistribution.txt')
open(303,file='OUTPUT/OrientationTensor.txt')
open(306,file='OUTPUT/a11.txt')
open(307,file='OUTPUT/Ln.txt')
open(308,file='OUTPUT/CurvatureDistrivution.txt')   !20181229 add
open(311,file='OUTPUT/DistanceDistrivution.txt')    !20181229 add
open(3,  file='OUTPUT/positions.out')
open(5,  file='OUTPUT/vels.out')
open(6,  file='OUTPUT/forces.out')

print *,"@@@"
print *, "Maximum number of threads        ", omp_get_max_threads()  
print *, "Number of threads being used     ", omp_get_num_threads()
write(301,*), "@@@"
write(301,*), "Maximum number of threads        ", omp_get_max_threads() 
write(301,*), "Number of threads being used     ", omp_get_num_threads()

end subroutine output_OpenFiles

!======================================================================

subroutine output_CloseFiles(  simParameters )  !2018/11/25
type(simulationParameters)  :: simParameters

write(*,*),   "Time Elapsed(1)",  OMP_get_wtime()-simParameters%start, "s"
write(301,*), "Time Elapsed(1)",  OMP_get_wtime()-simParameters%start, "s"

close(300)
close(301)
close(302)
close(303)
close(306)
close(307)
close(308)  !20181229 add
close(311)  !20181229 add
close (3) 
close (5)
close (6)

end subroutine output_CloseFiles

!======================================================================

subroutine output_data( fibers, hinges, simParameters )  !2018/10/10 修正
type(simulationParameters)               :: simParameters
type (fiber), dimension(:), allocatable  :: fibers
type(rod),    dimension(:), allocatable  :: hinges
integer(8)                               :: i, j, k, n

  n= simParameters%frame*simParameters%writ_period                      !2018/11/18  writ_period + 1
  
  open(4,file='OUTPUT/nbr_frames.txt')
  	  !write (4,*), ubound(hinges,1)
	   write (4,*), simParameters%frame, n, simParameters%nStep_Total   !2018/11/18  writ_period + 1
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

!======================================================================

subroutine output_Length( fibers, hinges, simParameters )  !2018/10/29 修正

type(simulationParameters)              :: simParameters
type(fiber), dimension(:),  allocatable :: fibers
type(rod),   dimension(:),  allocatable :: hinges

real(8), dimension(3)       :: coord
real(8)                     :: FiberLength_Total, FiberLength, SegmentLength, t, tt
real(8)                     :: length, lengthAvg
integer(8)                  :: i, j, k, mSegments
integer(8)                  :: mm, nn

    t= simParameters%time
    
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
    
    simParameters%current_time= OMP_get_wtime()-simParameters%start
    
      
    tt= t*1.e6 + 0.0                              !2018/08/11   化成整數,單位=micro seconds
    mm= mSegments                                 !2018/08/11   total segments number
    nn= ubound(fibers,1)                          !2018/08/11   total Fiber number
    length= 1000*FiberLength_Total                !2018/09/22   total Fiber length(mm)
    lengthAvg= 1000*FiberLength_Total/nn          !2018/09/22   mean  Fiber length(mm)
    
      write(*,*),"------------------------------------------------------------------------"
      write(*,101),"###", tt, nn, simParameters%nStep_Total, simParameters%AA(1,1), lengthAvg, simParameters%current_time, "s"
101   format( A4, F12.2, I8, 1X, I10, 2F12.6, F12.2, A2 )                  !2018/12/09
      write(*,*),"------------------------------------------------------------------------"
      
      write(300,201),tt, nn, simParameters%nStep_Total, simParameters%AA(1,1), lengthAvg, simParameters%current_time, "s"
201   format( F12.2, I8, 1X, I10, 2F12.6, F12.2, A2 )                  !2018/12/09  
      
      write(301,*),"------------------------------------------------------------------------"
      write(301,301),"###", tt, nn, simParameters%nStep_Total, simParameters%AA(1,1), lengthAvg, simParameters%current_time, "s"
301   format( A4, F12.2, I8, 1X, I10, 2F12.6, F12.2, A2 )                  !2018/12/09
      write(301,*),"------------------------------------------------------------------------"

!    write(300,100), tt, nn, mm, length, lengthAvg
!100 format(F12.0,2I10,2F15.6)
    
	write(307,110), tt,",", lengthAvg                   !2018/10/27   new output    
110 format(F14.2,A1,F12.6)
    
!    write(*,200),"@@@", tt, nn, mm, length, lengthAvg
!200 format(A4,F14.2,2I10,2F15.6)
!    write(301,210),"@@@", tt, nn, mm, length, lengthAvg
!210 format(A4,F14.2,2I10,2F15.6)
!pause

end subroutine output_Length

!======================================================================

subroutine output_LengthDistribution( fibers, indexA, simParameters )  !2018/10/29 修正

type(simulationParameters)              :: simParameters
type(fiber), dimension(:), allocatable  :: fibers
integer(8),  dimension(:), allocatable  :: indexA
integer(8)   :: i, j, maxSegments
real(8)      :: t, tt

t= simParameters%time
tt= t*1.d6

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

!tt= t*1.0e6 + 0.5    !2018/08/12  +0.5 的用意是4捨5入

maxSegments= ubound(indexA,1)

!write(*,100),  tt, ubound(fibers,1), maxSegments
!100 format( F18.8, I10, I10 )
    
write(302,200),  tt, ubound(fibers,1), maxSegments
200 format( F14.2, I8, I10 )

do j=1, maxSegments
    
    !print *, "&&&( ", j, indexA(j)
    
     write(302,300),  j, indexA(j)
300  format( I10, I10 ) 
     
end do
!pause

end subroutine output_LengthDistribution

!======================================================================

subroutine output_OrientationTensor( fibers, hinges, simParameters )  !2018/10/12 新增
type(fiber), dimension(:),   allocatable :: fibers
type(rod),   dimension(:),   allocatable :: hinges
type(simulationParameters)               :: simParameters
real(8),     dimension(3,3)              :: AA
real(8),     dimension(3)                :: ra
real(8)      :: t, tt, length, trace_A
integer(8)   :: ia, ja,  mm
integer(8)   :: i, j

t= simParameters%time
tt= t*1.d6

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

!tt= simParameters%time*1.0e6 + 0.5            !2018/10/12  +0.5 的用意是4捨5入
trace_A= AA(1,1) + AA(2,2) + AA(3,3)

AA= AA/trace_A

simParameters%AA= AA                             !2018/10/12 新增

   
!      write(*,101), tt, AA(1,1), AA(2,2), AA(3,3), AA(1,2), AA(1,3), AA(2,3)  !2018/12/09
!101   format( F10.0, 6F10.5 )                                             !2018/12/09
      
      write(303,301),tt ,ubound(fibers,1)
301   format( F14.2, I8 )                                                    !2018/11/18  new output

      write(306,401), tt, AA(1,1), AA(2,2), AA(3,3), AA(1,2), AA(1,3), AA(2,3)  !2018/12/07
401   format( F12.2, 2X, 6F12.6 )                                               !2018/12/07
!pause

end subroutine output_OrientationTensor

!======================================================================

subroutine output_PositionsForTheMomemt ( fibers, hinges )   !2018/10/12
type(fiber), dimension(:),   allocatable :: fibers
type(rod),   dimension(:),   allocatable :: hinges
integer(8)                               :: i, j, k

open(304,file='OUTPUT/PositionsForTheMoment.txt')            !2018/09/02

        k=1
        
        write (304,*), ubound(fibers,1)
        
        do i=1, ubound(fibers,1)
            
  	       write (304,*), fibers(i)%nbr_hinges
            
  	       do j=k, k+fibers(i)%nbr_hinges-1
               
  		      write (304,*), 0, real(hinges(k)%X_i(1),4), real(hinges(k)%X_i(2),4), real(hinges(k)%X_i(3),4)
              
		      k= k + 1
              
           end do

        end do
        
close(304)

end subroutine output_PositionsForTheMomemt

!======================================================================

subroutine output_FiberLengthModification( fibers, hinges )  !2018/09/22  add
type(fiber) , dimension(:)            :: fibers
type (rod), dimension(:)              :: hinges
logical                               :: periodic_boundary
integer(8)                            :: i,j,k,l,jp1
real(8), dimension(3)                 :: box_size, coord
real(8)                               :: length, FiberLength, ScaleFactor

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

!======================================================================

subroutine output_Initial_Positions_New( fibers, hinges, simParameters )   !2018/09/22    
implicit none
type(simulationParameters)             :: simParameters
type(fiber), dimension(:), allocatable :: fibers
type(rod),   dimension(:), allocatable :: hinges
real(8), dimension(3) :: box_size, coord
real(8)               :: del_X, del_Y, del_Z, x, y, z
integer(8)            :: mm1, mm2, mm3, mm1_A, mm1_B, mm2_A, mm2_B, mm3_A, mm3_B, ix, iy, iz
integer(8)            :: i, j, N_Fiber, N_hinge

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

!======================================================================

subroutine output_DynamicP_1848( flowcase_1848, simParameters )               !2018/10/11 

implicit none
type(simulationParameters)                 :: simParameters
type(DynamicP), dimension(:), allocatable  :: flowcase_1848                   !2018/10/11 add
integer(8)                                 :: h                               !2018/10/11 
       
        h= simParameters%h                                                    !2018/10/11 增加
        
        print *,     "###"               
        print *,     "### Dynamic Parameters CHANGE"                          !2018/10/11 增加
        print *,     "### StepNo., NextTime(Micro.Sec.), Gamma_dot, Viscosity"    !2018/10/11 增加
        print *,     "###",&
                      int(simParameters%h),&
                      int(0.5+1.0e6*flowcase_1848(h)%Duration),&
                      real(flowcase_1848(h)%Shearrate),&
                      real(flowcase_1848(h)%Viscosity)                        !2018/10/11 增加
        print *,     "###"                                                    !2018/10/11 增加
        write(301,*),"###"
        write(301,*),"### Dynamic Parameters CHANGE"                          !2018/10/11 增加
        write(301,*),"### StepNo., NextTime(Micro.Sec.), Gamma_dot, Viscosity"    !2018/10/11 增加
        write(301,*),"###",&
                      int(simParameters%h),&
                      int(0.5+1.0e6*flowcase_1848(h)%Duration),&
                      real(flowcase_1848(h)%Shearrate),&
                      real(flowcase_1848(h)%Viscosity)                        !2018/10/11 增加
        write(301,*),"###"                                                   !2018/10/11 增加

end subroutine output_DynamicP_1848  !2018/10/11 change name

!======================================================================

subroutine output_minmax_hinges( fibers, hinges ) !2018/10/27 change name
type(fiber), allocatable, dimension(:) :: fibers
type(rod)  , allocatable, dimension(:) :: hinges
real(8), dimension(3)                  :: coord, min_coor, max_coor
integer(8)                             :: i, j, k, m

!2018/08/05  Compute center of mass for all fibers coord
!2018/08/05  Compute center of mass for all fibers coord
print *,"@@@"
print *,"@@@ based on all hinges(mm)"
write(301,*), "@@@"
write(301,*), "@@@ based on all hinges(mm)"

m= 0
coord= 0
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
   coord= coord + hinges(j)%X_i
   m= m+1
end do
end do        
coord= coord/real(m)  !2018/08/05  coord= center of mass 

     write(*,100), "cen ", coord(1)*1000, coord(2)*1000, coord(3)*1000
100  format( A5, 3F15.6 )    !2018/12/09
     
     write(301,200), "cen ", coord(1)*1000, coord(2)*1000, coord(3)*1000
200  format( A5, 3F15.6 )    !2018/12/09

min_coor=   huge(0d0)
max_coor=  -huge(0d0)
do k=1,3
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-1
   coord= hinges(j)%X_i
   min_coor(k)= min( min_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
end do
end do  
end do      !2018/08/05 找出可以容納所有的hinges框架的座標範圍

!2018/12/09
     write(*,101), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
101  format( A5, 3F15.6 )
     write(*,102), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
102  format( A5, 3F15.6 )
     
     write(301,201), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
201  format( A5, 3F15.6 )
     write(301,202), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
202  format( A5, 3F15.6 )

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
   max_coor(k)= max( max_coor(k), coord(k) )   
end do
end do    !2018/08/05 找出可以容納所有的 segments 框架的座標範圍

!2018/12/09
     write(*,105), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
105  format( A5, 3F15.6 )
     write(*,106), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
106  format( A5, 3F15.6 )
     write(*,107), "@@@ "
107  format( A5 )
     write(301,205), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
205  format( A5, 3F15.6 )
     write(301,206), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
206  format( A5, 3F15.6 )
     write(301,207), "@@@ "
207  format( A5 )
!pause

end subroutine  output_minmax_hinges

!======================================================================

subroutine output_minmax_segments( fibers, hinges ) !2018/10/27 change name

type(fiber), allocatable, dimension(:) :: fibers
type(rod)  , allocatable, dimension(:) :: hinges
real(8), dimension(3)                  :: coord, min_coor, max_coor
integer(8)                             :: i, j, k, m

!2018/08/05  Compute center of mass for all fibers coord
print *,"@@@"
print *,"@@@ based on all segments(mm)"
write(301,*),"@@@"
write(301,*),"@@@ based on all segments(mm)"

m= 0
coord= 0
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
   coord= coord + (hinges(j)%X_i+hinges(j+1)%X_i)/2d0
   m= m+1
end do
end do        
coord= coord/real(m)  !2018/08/05  coord= center of mass 

     write(*,100), "cen ", coord(1)*1000, coord(2)*1000, coord(3)*1000
100  format( A5, 3F15.6 )
     
     write(301,200), "cen ", coord(1)*1000, coord(2)*1000, coord(3)*1000
200  format( A5, 3F15.6 )

min_coor=   huge(0d0)
max_coor=  -huge(0d0)
do k=1,3
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
   coord= (hinges(j)%X_i+hinges(j+1)%X_i)/2d0
   min_coor(k)= min( min_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
end do
end do  
end do      !2018/08/05 找出可以容納所有的hinges框架的座標範圍


!2018/12/09
     write(*,101), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
101  format( A5, 3F15.6 )
     write(*,102), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
102  format( A5, 3F15.6 )
     
     write(301,201), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
201  format( A5, 3F15.6 )
     write(301,202), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
202  format( A5, 3F15.6 )
     
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
   min_coor(k)= min( min_coor(k), coord(k) )
   max_coor(k)= max( max_coor(k), coord(k) )
end do  
end do      !2018/08/05 找出可以容納所有的hinges框架的座標範圍


!2018/12/09
     write(*,105), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
105  format( A5, 3F15.6 )
     write(*,106), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
106  format( A5, 3F15.6 )
     write(*,107), "@@@ "
107  format( A5 )
     write(301,205), "min ", min_coor(1)*1000, min_coor(2)*1000, min_coor(3)*1000
205  format( A5, 3F15.6 )
     write(301,206), "max ", max_coor(1)*1000, max_coor(2)*1000, max_coor(3)*1000
206  format( A5, 3F15.6 )
     write(301,207), "@@@ "
207  format( A5 )
!pause

end subroutine  output_minmax_segments

!======================================================================

subroutine output_OrientationTensor_OLD( fibers, hinges, simParameters )  !2018/10/31 修正
type(fiber), dimension(:),   allocatable :: fibers
type(rod),   dimension(:),   allocatable :: hinges
type(simulationParameters)               :: simParameters
real(8),     dimension(3,3)     :: AA, BB
real(8),     dimension(3)       :: ra
real(8)      :: t, length, trace_A, trace_B, B11_bar
integer(8)   :: ia, ja,  mm, nn, nbr_hinges
integer(8)   :: tt, i, j


      write(*,*),"@@@("
      write(*,  10 ) "@@@(","n0,","fiber,","nbr,","length,","a11,","a22,","a33,","A11_bar"
10    format( A5,A6,A8,A6,A11,A8,A10,A9,A12 )
      write(*,*),"@@@(" 
      write(301,*),"@@@("
      write(301,  12 ) "@@@(","n0,","fiber,","nbr,","length,","a11,","a22,","a33,","A11_bar"
12    format( A5,A6,A8,A6,A11,A8,A10,A9,A12 )
      write(301,*),"@@@("

AA= 0
mm= 0
nn= 0
do ia= 1, ubound (fibers,1)

   BB= 0
   do ja= fibers(ia)%first_hinge, fibers(ia)%first_hinge+fibers(ia)%nbr_hinges-2
      ra= hinges(ja+1)%X_i - hinges(ja)%X_i
          length= sqrt( dot_product(ra,ra) )
      ra= ra/length
      do i=1,3
      do j=1,3
         BB(i,j)= BB(i,j) + ra(i)*ra(j)*length   !2018/10/31 修正
         AA(i,j)= AA(i,j) + ra(i)*ra(j)*length   !2018/08/14  加入權重長度,因直徑一樣,長度代表體積
      end do
      end do
      mm= mm + 1 
   end do

   trace_B= BB(1,1) + BB(2,2) + BB(3,3)         !2018/10/31 修正
   BB=      BB/trace_B                          !2018/10/31 修正
   B11_bar= BB(1,1) + BB(3,3)                   !2018/10/31 修正
   nbr_hinges= fibers(ia)%nbr_hinges

   if( (BB(2,2) .lt. 0.10) .and. (nbr_hinges .gt. 20) ) then
      nn= nn + 1 
      write(*,  20 ) "@@@( ", nn, ia, nbr_hinges, trace_B, BB(1,1), BB(2,2), BB(3,3), B11_bar
20    format( A6, I5, I6, I6, F12.5, 4F10.5 ) 
      write(301,22 ) "@@@( ", nn, ia, nbr_hinges, trace_B, BB(1,1), BB(2,2), BB(3,3), B11_bar
22    format( A6, I5, I6, I6, F12.5, 4F10.5 )
   end if

end do

write(*,*),"@@@("
write(301,*),"@@@("      
!pause
      
tt= simParameters%time*1.0e6 + 0.5            !2018/10/12  +0.5 的用意是4捨5入
trace_A= AA(1,1) + AA(2,2) + AA(3,3)

AA= AA/trace_A

simParameters%AA= AA                          !2018/10/12 新增

!print *, "### trace ",trace_A, mm
print *, "###", tt, ubound (fibers,1)
print *, AA

write(301,*), "###", tt, ubound (fibers,1)
write(301,*), AA

write(303,*), tt, ubound (fibers,1)
write(303,*), AA
write(306,*), tt,",", AA(1,1)                 !2018/10/27   new output
!pause

end subroutine output_OrientationTensor_OLD


end module m_OutputData
!======================================================================