!=============================================================================

module m_FiberCalc

use m_DataStructures
use m_UtilityLib
use m_SetVelocity !2018/07/22 add
implicit none
contains

!=============================================================================

subroutine fiber_calc_tensor( fibers, hinges, simParameters )  !2018/10/08  change

type(simulationParameters)  :: simParameters
type(fiber), dimension(:)   :: fibers
type(rod),   dimension(:)   :: hinges

integer(8)                  :: i, j, n1, n2

do i=1, ubound(fibers,1)
    n1= fibers(i)%first_hinge
	n2= fibers(i)%first_hinge + fibers(i)%nbr_hinges - 2  !Changed JULY 31 2012
    
	do j=n1, n2
		call fiber_calc_tensor_hinge( hinges(j), hinges(j+1), simParameters )
    end do 
end do

end subroutine fiber_calc_tensor   !2018/07/21 change name
                          
!=============================================================================

subroutine fiber_calc( fibers, hinges, simParameters ) !2018/12/01   change name

type(simulationParameters)  :: simParameters
type(fiber), dimension(:)   :: fibers
type(rod),   dimension(:)   :: hinges

integer(8)                  :: i, j, n1, n2

do i=1, ubound(hinges,1)                 !2018/12/07 add

hinges(i)%rk_sum    = 0             !2018/12/07 add
hinges(i)%uoo_sum   = 0             !2018/12/07 add
hinges(i)%omega_sum = 0             !2018/12/07 add
hinges(i)%r_sum     = 0             !2018/12/07 add
hinges(i)%r_prod_sum= 0             !2018/12/07 add

hinges(i)%rk_vel    = 0             !2018/12/07 add

end do


do i=1, ubound(fibers,1)
    n1= fibers(i)%first_hinge
	n2= fibers(i)%first_hinge + fibers(i)%nbr_hinges - 2               !Changed JULY 31 2012
	do j=n1, n2
		call fiber_calc_hinge( hinges(j), hinges(j+1), simParameters ) !2018/12/01   change name
	end do
end do

end subroutine fiber_calc   !2018/07/21 change name

!=============================================================================

subroutine fiber_calc_tensor_hinge( hinge1, hinge2, simParameters )    !2018/12/01  change

type(simulationParameters)  :: simParameters
type (rod)                  :: hinge1, hinge2

real(8), dimension(3,3,3) :: eps
real(8), dimension(3,3)   :: Id, dd, E_oo

real(8), dimension(3)     :: middle_position, X_j, vel, omega, d, part
real(8)                   :: viscosity, segment_length, pi, c1, c2, XA, YA, XC, YC, YH       !2018/12/08
integer(8)                :: i, j, k, l_


hinge1%A= 0
hinge1%C= 0
hinge1%H= 0

XA= simParameters%X_A
YA= simParameters%Y_A
XC= simParameters%X_C
YC= simParameters%Y_C
YH= simParameters%Y_H
E_oo= simParameters%E_oo

hinge1%r=        hinge2%X_i - hinge1%X_i
hinge1%length=   sqrt(dot_product(hinge1%r,hinge1%r))
d=               hinge1%r/hinge1%length
middle_position= ( hinge1%X_i + hinge2%X_i )/2d0

pi            = simParameters%pi
viscosity     = simParameters%viscosity       !2018/12/01  change
segment_length=  hinge1%length

c1= 3.0*pi*viscosity*( segment_length )
c2=     pi*viscosity*( segment_length )**3


call set_velocity( middle_position, vel, omega, simParameters )    !§ó´« 2018/12/01


!2018/12/01 new version, eps(1,2,3)= eps(2,3,1)= eps(3,1,2)= 1 ;  eps(3,2,1)= eps(2,1,3)= eps(1,3,2)= -1
part= 0
do k =1, 3
   part(1)= part(1) +  d(3)*d(k)*E_oo(2,k) - d(2)*d(k)*E_oo(3,k)
   part(2)= part(2) +  d(1)*d(k)*E_oo(3,k) - d(3)*d(k)*E_oo(1,k)
   part(3)= part(3) +  d(2)*d(k)*E_oo(1,k) - d(1)*d(k)*E_oo(2,k)
enddo

Id = 0
forall(j = 1:3) Id(j,j) = 1  !2018/12/01  unit tensor

dd= outerProd(d,d)

hinge1%A=  c1*( XA*dd + YA*(Id-dd) )       !2018/12/07 p.27 equation (3.6)
hinge1%C=  c2*( XC*dd + YC*(Id-dd) )       !2018/12/07 p.27 equation (3.6)
hinge1%H= -c2*( YH*part )                  !2018/12/07 p.27 equation (3.6)
hinge1%Auf_oo  = AdotU( hinge1%A, vel   )  !2018/12/07 p.27 equation (3.6)
hinge1%omega_oo= AdotU( hinge1%C, omega )  !2018/12/07 p.27 equation (3.6)

hinge1%ave_viscosity= viscosity

hinge1%T= 0 

end subroutine fiber_calc_tensor_hinge

!=============================================================================

subroutine fiber_calc_hinge( hinge1, hinge2, simParameters )  !2018/12/01  change

type(simulationParameters)  :: simParameters                  !2018/12/01  change
type (rod)                  :: hinge1, hinge2

integer(8)              :: i, j

real(8), dimension(3)   :: middle_position, X_start, X_j, X_local, vel, omega, unitar_vec, X_k, rv  !2018/12/04
real(8), dimension(3,1) :: X_local_mat, X_k_mat                            !2018/12/04 add
real(8), dimension(1,3) :: X_local_mat_T, u_fluid_mat_T, X_k_mat_T         !2018/12/04 add
real(8), dimension(3,3) :: XK                                                      !2018/12/07 add
real(8)                 :: r_fiber, viscosity, segment_length, pi, c1, c2, trace_XK   !2018/12/07


pi          = simParameters%pi              !2018/12/07  add
r_fiber     = simParameters%r_fiber         !2018/12/01  change
viscosity   = simParameters%viscosity       !2018/12/01  change

c1= (6.d0*pi*r_fiber)*viscosity
c2= (8.d0*pi*r_fiber**3d0)*viscosity
    

hinge1%r        = hinge2%X_i - hinge1%X_i

segment_length  = sqrt(dot_product(hinge1%r,hinge1%r))

unitar_vec      = (hinge2%X_i-hinge1%X_i)/segment_length

hinge1%nbr_beads= floor((segment_length)/(2*r_fiber)) !Changed 9/28/2014 To acommodate odd number of beads
 
middle_position= (hinge1%X_i+hinge2%X_i)/2d0

!if modulo(hinge1%nbr_beads,2).le.1e-20) then !Changed 9/28/2014 To acommodate odd number of beads
!    X_start=middle_position-unitar_vec*(2*r_fiber*real(hinge1%nbr_beads)/2-r_fiber)
!else
    X_start= middle_position-unitar_vec*(2*r_fiber*real(hinge1%nbr_beads)/2-r_fiber)
!end if

hinge1%c1_nbeads= c1*hinge1%nbr_beads
hinge1%c2_nbeads= c2*hinge1%nbr_beads

do i=1, hinge1%nbr_beads

	X_j= X_start + (i-1)*unitar_vec*(2*r_fiber)
    
    X_local= X_j - hinge1%X_i
    X_k    = X_j - middle_position                         !2018/12/04
		
	call set_velocity( X_j, vel, omega, simParameters )    !§ó´« 2018/12/01
	 
    do j=1, 3
        
         X_k_mat(  j,1)= X_k(j)
         X_k_mat_T(1,j)= X_k(j)
         
		 X_local_mat(  j,1)= X_local(j)		
		 X_local_mat_T(1,j)= X_local(j)
		 u_fluid_mat_T(1,j)= vel(j)
    end do
        
    hinge1%uoo_sum= hinge1%uoo_sum + c1*vel
    
	!print *, "Vx", Vx, "Vy", Vy, "Vort", Vorticity_z,"Viscosity"

    hinge1%rk_sum    = hinge1%rk_sum     + c1*X_k                     !2018/12/04  add

	hinge1%omega_sum = hinge1%omega_sum  + c2*omega
    hinge1%r_sum     = hinge1%r_sum      + c1*X_local
    
    XK= c1*matmul(X_local_mat, X_k_mat_T)

        trace_XK= XK(1,1) + XK(2,2) + XK(3,3)
        
        XK(1,1)= trace_XK - XK(1,1)
        XK(2,2)= trace_XK - XK(2,2)
        XK(3,3)= trace_XK - XK(3,3)

    hinge1%r_prod_sum= hinge1%r_prod_sum + XK                                             !2018/12/04  add
    
    rv= cross( X_k, vel )*0.5d0
    
    hinge1%rk_vel= hinge1%rk_vel + c1*rv                              !2018/12/07  add
    hinge2%rk_vel= hinge2%rk_vel + c1*rv                              !2018/12/07  add

end do

hinge1%ave_viscosity= viscosity

hinge1%T= 0 

end subroutine fiber_calc_hinge    !2018/07/21 change name

end module m_FiberCalc             !2018/07/21 change name
!=============================================================================