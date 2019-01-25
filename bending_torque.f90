!====================================================================
module m_BendingTorque

use m_DataStructures
use m_UtilityLib

implicit none
contains

subroutine bending_torque_whole( fibers, hinges, simParameters )  !2018/10/08  修正

implicit none
type(simulationParameters)     :: simParameters
type (fiber), dimension(:)     :: fibers
type (rod), dimension(:)       :: hinges
real(8), dimension(3)          :: T_b
real(8)                        :: E_Young, Inertia_Moment
integer(8)                     :: i, j

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(FIBERS, HINGES, E_YOUNG, INERTIA_MOMENT )
!$OMP DO 

E_Young=        simParameters%E_Young
Inertia_Moment= simParameters%Inertia_Moment


do i=1, ubound(hinges, 1)
	hinges(i)%T=0d0
	hinges(i)%alpha=0d0
end do

!$OMP END DO
!$OMP DO !SCHEDULE(DYNAMIC)

do i=1, ubound(fibers,1)
	if(fibers(i)%nbr_hinges.GE.3) then 
        
		!print *, "TEST SECOND HINGE", fibers(i)%first_hinge+1
		!print *, "TEST PENULTIMA HINGE", fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
		!STOP
        
		do j=fibers(i)%first_hinge+1, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
		!do j=2,14
			call bending_torque_localhinge(T_b,&
                             	 hinges(j-1)%X_i,&
                             	 hinges(j)%X_i,&
                             	 hinges(j+1)%X_i,&
                             	 E_Young,&
                             	 Inertia_Moment,&
                                 hinges(j)%alpha)
        		!print *, "TB", T_b

			hinges(j-1)%T= hinges(j-1)%T + T_b
            hinges(j)%T  = hinges(j)  %T - T_b
		end do
	end if
end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!do i=1, ubound(hinges,1)
!	print *,"BENDING TORQUE", hinges(i)%T
!end do

end subroutine bending_torque_whole

!====================================================================

subroutine bending_torque_localhinge (T_b, pt1, pt2, pt3, E_Young, Inertia_Moment, alpha)
implicit none
real(8), dimension(3) :: pt1, pt2, pt3, vec1, vec2, e_theta, cross_product, T_b
real(8)               :: l1, l2, k_bending, E_Young, Inertia_Moment, alpha, beta,  eps
real(8),parameter     :: pi=3.14159265358979

eps=10*tiny(E_Young)

vec2= pt3-pt2
vec1= pt2-pt1

l1= sqrt(abs(dot_product(vec1,vec1)))
l2= sqrt(abs(dot_product(vec2,vec2)))
beta= dot_product(vec1, vec2)/(l1*l2)  !2018/07/31  修正 

k_bending= (2d0*E_Young* Inertia_Moment)/(l1+l2)

if (beta.GE.1.0) then            !2018/07/31  修正 
	alpha= 0d0
	else if (beta.le.-1.0) then  !2018/07/31  修正 
	alpha= pi                    !Changed 9/28/2014
	else
	alpha= acos( beta )          !2018/07/31  修正 
end if

!print *, "ALPHA ANGLE", alpha

 cross_product= cross(vec1, vec2)
 e_theta= cross_product/ (sqrt(dot_product(cross_product ,cross_product)))
  
!if (alpha .GT. pi/2)	then
!	print *,"Error: buckled fiber: the angle between two segments is greater than 90 degrees. Program must stop."
! 	stop
!end if
!NOTE INSERT WARNING FOR ALPHA=180
!print *,"EPS", eps

T_b=0
if (abs(cross_product(1))<=eps .and.&
    abs(cross_product(2))<=eps .and.&
    abs(cross_product(3))<=eps) then 
	T_b =0
else if (l1+l2<=eps) then
	T_b=0
else if (l1*l2<=eps) then
	T_b=0
else
	T_b= k_bending*alpha*e_theta  !2018/07/31  equation (3.21)
end if

!print *, "CROSS PRODUCT", cross(vec1, vec2), "ANGLE", alpha, "Torque", T_b 

end subroutine bending_torque_localhinge

end module m_BendingTorque
!====================================================================