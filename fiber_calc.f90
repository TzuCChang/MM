!====================================================================
module m_FiberCalc
    
use m_DataStructures
use m_UtilityLib

use m_SetVelocity !2018/07/22 add

implicit none
contains

subroutine fiber_calc_tensor( fibers, hinges, simParameters )  !2018/10/08  change

type(simulationParameters)  :: simParameters
type (rod),   dimension(:)  :: hinges
type (fiber), dimension(:)  :: fibers

integer(8)                  :: i, j, k, first, last, flow_case
real(8)                     :: r_fiber, viscosity, gamma_dot, epsilon_dot

!!call omp_set_num_threads(6)
!$OMP PARALLEL default(private) SHARED(r_fiber, viscosity, gamma_dot, epsilon_dot, flow_case, simParameters, hinges, fibers )
!$OMP DO !SCHEDULE(DYNAMIC) ! all of tehm should be shared by default anyway

flow_case   = simParameters%flow_case       !2018/10/08  change
r_fiber     = simParameters%r_fiber         !2018/10/08  change
viscosity   = simParameters%viscosity       !2018/10/08  change
gamma_dot   = simParameters%gamma_dot       !2018/10/08  change
epsilon_dot = simParameters%epsilon_dot     !2018/10/08  change


do i=1, ubound(fibers,1)
    first=fibers(i)%first_hinge
	last =fibers(i)%first_hinge+fibers(i)%nbr_hinges-2 !Changed JULY 31 2012
    
	do j=first, last
		!print *,"indices",i, j
		call fiber_calc_tensor_hinge( hinges(j),&
                                      hinges(j+1),&
                                      r_fiber,&
                                      viscosity,&
                                      gamma_dot,&
                                      epsilon_dot,&
                                      flow_case,&
                                      simParameters )
    end do 
end do

!$OMP END DO NOWAIT
!$OMP END PARALLEL

end subroutine fiber_calc_tensor   !2018/07/21 change name
                          
!====================================================================                          
subroutine fiber_calc( fibers, hinges, simParameters ) !2018/10/08   change name

type(simulationParameters)  :: simParameters
type (rod),   dimension(:)  :: hinges
type (fiber), dimension(:)  :: fibers

integer(8)                  :: i, j, k, first, last, flow_case
real(8)                     :: r_fiber, viscosity, gamma_dot, epsilon_dot

!!call omp_set_num_threads(6)
!$OMP PARALLEL
!$OMP DO PRIVATE(first, last, i, j)

flow_case   = simParameters%flow_case       !2018/10/08  change
r_fiber     = simParameters%r_fiber         !2018/10/08  change
viscosity   = simParameters%viscosity       !2018/10/08  change
gamma_dot   = simParameters%gamma_dot       !2018/10/08  change
epsilon_dot = simParameters%epsilon_dot     !2018/10/08  change

do i=1, ubound(fibers,1)
    first=fibers(i)%first_hinge
	last =fibers(i)%first_hinge+fibers(i)%nbr_hinges-2 !Changed JULY 31 2012
	do j=first, last
		!print *,"indices",i, j
		call fiber_calc_hinge( hinges(j),&
                               hinges(j+1),&
                               r_fiber,&
                               viscosity,&
                               gamma_dot,&
                               epsilon_dot,&
                               flow_case )
	end do
end do

!$OMP END DO
!$OMP END PARALLEL

end subroutine fiber_calc   !2018/07/21 change name
                          
!====================================================================
subroutine fiber_calc_tensor_hinge( hinge1,&    !2018/07/21 change name
                                    hinge2,&
                                    r_fiber,&
                                    viscosity,&
                                    gamma_dot,&
                                    epsilon_dot,&
                                    flow_case,&
                                    simParameters )

type (rod)             :: hinge1, hinge2
real(8), dimension(3)  :: middle_position, X_start, X_j, X_local, vel, omega, unitar_vec, d, part
real(8), dimension(3,1):: X_local_mat
real(8), dimension(1,3):: X_local_mat_T, u_fluid_mat_T
real(8), dimension(3,3):: Id
real(8), dimension(3,3,3):: eps
integer(8)             :: i, j, nbr_fibers_per_segment, flow_case, aa, k, l_
real(8)                :: r_fiber, viscosity, gamma_dot, epsilon_dot, segment_length
real(8), parameter     :: pi=3.141592
type(simulationParameters)  :: simParameters

hinge1%A =0
hinge1%H =0
hinge1%C =0

hinge1%r= hinge2%X_i-hinge1%X_i

hinge1%length = sqrt(dot_product(hinge1%r,hinge1%r))
segment_length= hinge1%length
d = hinge1%r/hinge1%length

middle_position= (hinge1%X_i+hinge2%X_i)/2d0

if ( flow_case == 1 ) then   !2018/10/20

    call set_velocity( middle_position,vel,omega,gamma_dot,epsilon_dot,flow_case )    !§ó´« 2018/07/14

else if ( flow_case == 1848 ) then

    call set_velocity_NEW(middle_position,vel,omega,gamma_dot,epsilon_dot,flow_case) !·s¼W 2018/07/14

end if   !2018/10/20 

hinge1%u_oo =    vel
hinge1%omega_oo= omega

!for the force we need 2 Constants

Id = 0

forall(j = 1:3) Id(j,j) = 1

hinge1%A= 6.0*pi*viscosity*segment_length/2*( simParameters%X_A*outerProd(d,d) + simParameters%Y_A*(Id-outerProd(d,d)) )  !2018/07/31 p.27 equation (3.6)

!for the Torque we need 3 Constants
 
hinge1%C = 8.0*pi*viscosity*(segment_length/2)**3*( simParameters%X_C*outerProd(d,d) + simParameters%Y_C*(Id - outerProd(d,d)) )  !2018/07/31 p.27 equation (3.6)


!print *, ' a ', a

part = 0

do i = 1 ,3
    part(i) =0
    do j =1 ,3
        if(i .eq. j) then
            cycle
        endif
        do l_ =1 , 3
            if (l_.eq.j .or. l_.eq.i) then
                cycle
            endif
            do k =1, 3
              part(i) = part(i) +  simParameters%eps(i,j,l_)*d(l_)*d(k)*simParameters%E_oo(j,k)                
            enddo
        enddo
    enddo
enddo

! Spheroid_T =  matmul(C , (omega - hinge1%omega)) -8*pi*viscosity* segment_length**3 * Y_H *part 

hinge1%H = -8*pi*viscosity* (segment_length/2)**3 * simParameters%Y_H *part       !2018/07/31 p.27 equation (3.6)


!print *, ''
!print *, ' segment_length ', segment_length
!!print *, ' r_e ', r_e
!print *, 'hinge1%X_i' , hinge1%X_i
!print *, 'hinge2%X_i' , hinge2%X_i
!print *, 'r' , hinge1%r
!print *, 'd' , d
!print *, sqrt(dot_product(d,d))
!print *, 'd d ',outerProd(d,d)
!print *, "A ", hinge1%A
!print *, "X_A ",simParameters%X_A
!print *, "Y_A ",simParameters%Y_A
!print *,"C ", hinge1%C
!print *, ' segment_length ', segment_length
!print *, ' d ', d
!print *, ' part ' , part
!print *, "H ",hinge1%H

hinge1%ave_viscosity=viscosity

!Bending moments calc

hinge1%T=0 

end subroutine fiber_calc_tensor_hinge   !2018/07/21 change name

!====================================================================
subroutine fiber_calc_hinge( hinge1,&     !2018/07/21 change name
                             hinge2,&
                             r_fiber,&
                             viscosity,&
                             gamma_dot,&
                             epsilon_dot,&
                             flow_case )

type (rod)             :: hinge1, hinge2
real(8), dimension(3)  :: middle_position, X_start, X_j, X_local, vel, omega, unitar_vec
real(8), dimension(3,1):: X_local_mat
real(8), dimension(1,3):: X_local_mat_T, u_fluid_mat_T
integer(8)             :: i, j, nbr_fibers_per_segment, flow_case
real(8)                :: r_fiber, viscosity, gamma_dot, epsilon_dot, segment_length

hinge1%u_fluid_sum    = 0
hinge1%omega_fluid_sum= 0
hinge1%r_sum          = 0
hinge1%r_prod_sum     = 0
hinge1%r_times_u_sum  = 0

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

do i=1, hinge1%nbr_beads

	X_j= X_start + (i-1)*unitar_vec*(2*r_fiber)
    
    X_local= X_j - hinge1%X_i
		

    if (flow_case == 1 ) then     !2018/10/20
	call set_velocity( X_j,&
                       vel,&
                       omega,&
                       gamma_dot,&
                       epsilon_dot,&
                       flow_case )
    
else if ( flow_case == 1848 ) then
    call set_velocity_NEW( X_j,&
                             vel,&
                             omega,&
                             gamma_dot,&
                             epsilon_dot,&
                             flow_case )
    end if     !2018/10/20
	 
    do j=1, 3
		 X_local_mat(  j,1)= X_local(j)		
		 X_local_mat_T(1,j)= X_local(j)
		 u_fluid_mat_T(1,j)= vel(j)
    end do
        
    hinge1%u_fluid_sum= hinge1%u_fluid_sum + vel
    
	!print *, "Vx", Vx, "Vy", Vy, "Vort", Vorticity_z,"Viscosity"
    
	hinge1%omega_fluid_sum= hinge1%omega_fluid_sum + omega
    hinge1%r_sum          = hinge1%r_sum           + X_local
    hinge1%r_prod_sum     = hinge1%r_prod_sum      + matmul(X_local_mat, X_local_mat_T)
	hinge1%r_times_u_sum  = hinge1%r_times_u_sum   + matmul(X_local_mat, u_fluid_mat_T)
    
end do

hinge1%ave_viscosity= viscosity

hinge1%T= 0 

end subroutine fiber_calc_hinge    !2018/07/21 change name

end module m_FiberCalc             !2018/07/21 change name
