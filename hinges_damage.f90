!====================================================================
module m_HingesDamage

use m_DataStructures
use m_UtilityLib

implicit none
contains

subroutine hinges_curv_distribution( fibers, hinges, simParameters )  !2018/11/20 add
implicit none
type(simulationParameters)                :: simParameters
type (fiber), dimension(:), allocatable   :: fibers, fibers_temp
type (rod)  , dimension(:), allocatable   :: hinges, hinges_temp
real(8)        :: min_curv,curv, t
integer(8)     :: k, i, j, m, n
integer(8), dimension(10)  ::nC

min_curv = simParameters%min_curv
t=         simParameters%time                                !2018/11/19 add


do i=1, ubound(hinges,1)
	hinges(i)%is_broken    =.false.
	hinges(i)%is_separated =.false.
	hinges(i)%curv = huge(0d0)  !2018/09/05
end do

k=0

!print *, "number of fibers" ,  ubound(fibers,1)
do i=1, ubound(fibers,1)
	if(fibers(i)%nbr_hinges.ge.3) then 
		do j=fibers(i)%first_hinge+1, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2

            call find_curvature(hinges(j-1)%X_i, hinges(j)%X_i,hinges(j+1)%X_i, curv)
            
            hinges(j)%curv = curv      !2018/09/05         
            
            if (.not. isnan(curv) .and. curv.le.min_curv) then 
                !print *, "curv", curv
                !pause                
 				hinges(j)%is_broken = .true.		 				
				k=k+1
			end if
		end do
	end if
end do

simParameters%k= k

do i=1, ubound(fibers,1)
	hinges(fibers(i)%first_hinge)%is_separated=.true.
end do

nC= 0                                                      
do j=1, ubound(hinges,1)                                    
   if(      hinges(j)%curv .lt. 1.d0*min_curv ) then       
            nC(1)= nC(1) + 1                               
   else if( hinges(j)%curv .lt. 2.d0*min_curv ) then       
            nC(2)= nC(2) + 1                                
   else if( hinges(j)%curv .lt. 3.d0*min_curv ) then        
            nC(3)= nC(3) + 1                               
   else if( hinges(j)%curv .lt. 4.d0*min_curv ) then        
            nC(4)= nC(4) + 1                                
   else if( hinges(j)%curv .lt. 5.d0*min_curv ) then        
            nC(5)= nC(5) + 1                                
   else                                                     
            nC(6)= nC(6) + 1                                
   end if                                                   
end do

       write(*,100),  "@@@",1.d6*t, nC(1), nC(2), nC(3), nC(4), nC(5), nC(6)  !2018/11/19 add
100    format( A4, F12.2, 6I9 )
       write(301,101),"@@@",1.d6*t, nC(1), nC(2), nC(3), nC(4), nC(5), nC(6)  !2018/11/19 add
101    format( A4, F12.2, 6I9 )
       write(308,101),"@@@",1.d6*t, nC(1), nC(2), nC(3), nC(4), nC(5), nC(6) 
102    format( A4, F12.2, 6I9 )

!call hinges_SortOrder_curv( hinges )   !2018/10/05  關掉, 量大, 時間花太長


end subroutine  hinges_curv_distribution


subroutine hinges_damage( fibers, hinges, simParameters )

implicit none
type(simulationParameters)                :: simParameters
type (fiber), dimension(:), allocatable   :: fibers, fibers_temp
type (rod)  , dimension(:), allocatable   :: hinges, hinges_temp
real(8)        :: fac
integer(8)     :: k, i, j, m, n
integer(8), dimension(10)  ::nC

k= simParameters%k                                                 !2018/11/23

allocate(hinges_temp(ubound(hinges,1)+k))
allocate(fibers_temp(ubound(fibers,1)+k))

!print *, "Neo Ubound *****", (ubound(fibers,1)+k), (ubound(hinges,1)+k)

fac=0.000001*simParameters%r_fiber  !2018/11/23 修正 fac=0.0001*r_fiber,(2018/09/13 origion fac=0.07*r_fiber)

k=0

m=1
do i=1, ubound(hinges,1)
	if(hinges(i)%is_broken .eqv. .true.) then
		k=k+1
	end if
	if((hinges(i)%is_broken .eqv. .true.).or.(hinges(i)%is_separated .eqv. .true.)) then
		fibers_temp(m)%first_hinge=k+i
		m=m+1
	end if

end do

do i=1, ubound(fibers_temp,1)-1
	fibers_temp(i)%nbr_hinges=fibers_temp(i+1)%first_hinge-fibers_temp(i)%first_hinge
end do
fibers_temp(ubound(fibers_temp,1))%nbr_hinges=ubound(hinges_temp,1)-fibers_temp(ubound(fibers_temp,1))%first_hinge+1


k=1
do i=1, ubound(hinges,1)  !2018/09/22  修正錯誤

    hinges_temp(k)= hinges(i)    
  
	if (hinges(i)%is_broken .eqv. .true.)then
        
		!print *, "BREAKAGE"      
        !adding to alocate all properties of hinges(i) to hinges_temp(i)
        ! updating position to ensure not overlapping
                          
		hinges_temp(k)%X_i= hinges(i)%X_i - fac*(hinges(i)%X_i-hinges(i-1)%X_i)&
                                             /sqrt(dot_product(hinges(i)%X_i-hinges(i-1)%X_i,hinges(i)%X_i-hinges(i-1)%X_i))
        hinges_temp(k)%is_stationary = hinges(i)%is_stationary

        !adding to alocate all properties of hinges(i) to hinges_temp(i)
        ! updating position to ensure not overlapping
        
		k = k + 1

        hinges_temp(k)= hinges(i) !2018/08/04 修正錯誤 
                          
		hinges_temp(k)%X_i= hinges(i)%X_i + fac*(hinges(i+1)%X_i-hinges(i)%X_i)&
                                             /sqrt(dot_product( hinges(i+1)%X_i-hinges(i)%X_i, hinges(i+1)%X_i-hinges(i)%X_i))
        hinges_temp(k)%is_stationary = hinges(i)%is_stationary

    end if
    
    k = k + 1

end do

deallocate(hinges)
deallocate(fibers)

allocate(hinges(ubound(hinges_temp,1)))
allocate(fibers(ubound(fibers_temp,1)))

hinges=hinges_temp
fibers=fibers_temp

deallocate(hinges_temp)
deallocate(fibers_temp)


end subroutine hinges_damage



subroutine hinges_SortOrder_curv( hinges )

integer(8), dimension(:),  allocatable :: indexA, indexB
integer(8)                             :: i1, i2, i3, k1, k2, k3, mm

type (rod)  , dimension(:), allocatable :: hinges
real(8)                                :: c1, c2, c3

mm= ubound(hinges,1)
allocate( indexA( mm ) )
allocate( indexB( mm ) )

    do i1= 1, mm 
       indexA(i1)= i1
       indexB(i1)= i1
    end do
        
    do i1= 1, mm
        
       c1= hinges(i1)%curv
       k1= indexB(i1)
       
       i3= i1
       k3= k1
       c3= c1
       
       do i2= i1+1, mm
           
          k2= indexB(i2)
          c2= hinges(i2)%curv
          if( c3>c2 )  then
                  i3= i2
                  k3= k2
                  c3= c2
          end if
       end do
       
       i2= indexA(i3 )
       indexA(i3)= indexA(i1)
       indexA(i1)= i2
       
       k2= indexB(i3)
       indexB(i3)= indexB(i1)
       indexB(i1)= k2
           
       hinges(i1)%curv= c3      
       hinges(i3)%curv= c1       
       !print *,"B01",i1, i3, indexA(i1), c3
    end do
    !pause

    !print *,"@@@( minCurv", real(hinges(1)%curv,4), real(hinges(2)%curv,4)
    do i1=1,mm
       if( hinges(i1)%curv .lt. 0.1 ) then
           !print *,"B02",i1, indexB(i1), hinges(i1)%curv
           !write(301,*),"B02",i1, indexB(i1), hinges(i1)%curv
       end if
    end do    
    !pause
    
deallocate(indexA)
deallocate(indexB)

end subroutine  hinges_SortOrder_curv



end module m_HingesDamage
!====================================================================