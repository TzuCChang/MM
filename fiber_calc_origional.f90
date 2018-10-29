!====================================================================
module m_FiberCalc
    
use m_DataStructures
use m_UtilityLib

use m_SetVelocity !2018/07/22 add

implicit none
contains

subroutine fiber_calc_tensor( fibers,&  !2018/07/21 change name
                              hinges,&
                              r_fiber,&
                              viscosity,&
                              gamma_dot,&
                              epsilon_dot,&
                              flow_case,&
                              simParameters )

type (rod),   dimension(:)  :: hinges
type (fiber), dimension(:)  :: fibers
integer(8)                  :: i, j, k, first, last, flow_case
real(8)                     :: r_fiber, viscosity, gamma_dot, epsilon_dot
type(simulationParameters)  :: simParameters

!!call omp_set_num_threads(6)
!$OMP PARALLEL default(private) SHARED(r_fiber, viscosity, gamma_dot, epsilon_dot, flow_case, simParameters, hinges, fibers )
!$OMP DO !SCHEDULE(DYNAMIC) ! all of tehm should be shared by default anyway

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
subroutine fiber_calc( fibers,&   !2018/07/21 change name
                       hinges,&
                       r_fiber,&
                       viscosity,&
                       gamma_dot,&
                       epsilon_dot,&
                       flow_case )

type (rod),   dimension(:)  :: hinges
type (fiber), dimension(:)  :: fibers
integer(8)                  :: i, j, k, first, last, flow_case
real(8)                     :: r_fiber, viscosity, gamma_dot, epsilon_dot

!!call omp_set_num_threads(6)
!$OMP PARALLEL
!$OMP DO PRIVATE(first, last, i, j)

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

call set_velocity( middle_position,vel,omega,gamma_dot,epsilon_dot,flow_case )    !§ó´« 2018/07/14

!call set_velocity_NEW(middle_position,vel,omega,gamma_dot,epsilon_dot,flow_case) !·s¼W 2018/07/14

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
		
	call set_velocity( X_j,&
                       vel,&
                       omega,&
                       gamma_dot,&
                       epsilon_dot,&
                       flow_case )
	 
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

!Bending moments calc

hinge1%T= 0 

end subroutine fiber_calc_hinge    !2018/07/21 change name
 

!====================================================================
subroutine computeTandF_bead( hinge1,&
                              hinge2,&
                              r_fiber,&
                              viscosity,&
                              gamma_dot,&
                              epsilon_dot,&
                              flow_case, &
                              bead_F,&
                              bead_T )

type (rod)             :: hinge1, hinge2
real(8), dimension(3)  :: middle_position, X_start, X_j, X_local, vel, omega, unitar_vec, bead_F, bead_T, F_h
real(8), dimension(3,1):: X_local_mat
real(8), dimension(1,3):: X_local_mat_T, u_fluid_mat_T
integer(8)             :: i, j, nbr_fibers_per_segment, flow_case
real(8)                :: r_fiber, viscosity, gamma_dot, epsilon_dot, segment_length, drag_coeff_vel, drag_coeff_omega
real(8), parameter     :: pi=3.141592

drag_coeff_vel = 6d0*pi*r_fiber*viscosity
drag_coeff_omega= 8d0*pi*r_fiber**3d0*viscosity
bead_F =0;
bead_T =0;
hinge1%r        = hinge2%X_i-hinge1%X_i

segment_length  = sqrt(dot_product(hinge1%r,hinge1%r))

unitar_vec      = (hinge2%X_i-hinge1%X_i)/segment_length

hinge1%nbr_beads= floor((segment_length)/(2*r_fiber))
 
middle_position=(hinge1%X_i+hinge2%X_i)/2d0

!if modulo(hinge1%nbr_beads,2).le.1e-20) then !Changed 9/28/2014 To acommodate odd number of beads
!    X_start=middle_position-unitar_vec*(2*r_fiber*real(hinge1%nbr_beads)/2-r_fiber)
!else

    X_start=middle_position-unitar_vec*(2*r_fiber*real(hinge1%nbr_beads)/2-r_fiber)
!end if
do i=1, hinge1%nbr_beads

	X_j=X_start+(i-1)*unitar_vec*(2*r_fiber)
    X_local=X_j-hinge1%X_i
    
    ! interpolate fluid velocity and vorticity
    
	call set_velocity( X_j,vel,omega,gamma_dot,epsilon_dot,flow_case )
    
    !print *, "vel ", vel
    F_h = drag_coeff_vel * (vel )
    ! force on current hinge 
    !F_h = drag_coeff_vel * (vel - hinge1%v_i - cross(hinge1%omega, X_local))
    
    bead_F = bead_F +   F_h   
    
    !bead_T = bead_T +  drag_coeff_omega * (omega -hinge1%omega) + cross(X_local,F_h)
    bead_T = bead_T +  drag_coeff_omega * (omega) + cross(X_local,F_h)

end do

!print *, " mu " ,  drag_coeff_omega
!print *, "omega" , omega
!print *, " drag_coeff_omega " ,  drag_coeff_omega
!print *, "omega" , omega
!print *, " drag_coeff_omega * (omega) " , drag_coeff_omega * (omega)

end subroutine computeTandF_bead

!====================================================================
subroutine computeTandF_Spheroid( hinge1,&
                                  hinge2,&
                                  r_fiber,&
                                  viscosity,&
                                  gamma_dot,&
                                  epsilon_dot,&
                                  flow_case, &
                                  Spheroid_F,&
                                  Spheroid_T )

type (rod)             :: hinge1, hinge2
real(8), dimension(3)  :: middle_position, X_start, X_j, X_local, vel, omega, unitar_vec, Spheroid_F, Spheroid_T, F_h, d, part
real(8), dimension(3,1):: X_local_mat
real(8), dimension(1,3):: X_local_mat_T, u_fluid_mat_T
real(8), dimension(3,3):: Id, A, E_oo, C
real(8), dimension(3,3,3):: eps
integer(8)             :: i, j, k, l_, nbr_fibers_per_segment, flow_case, aa
real(8)                :: r_fiber, viscosity, gamma_dot, epsilon_dot, segment_length, drag_coeff_vel, drag_coeff_omega
real(8), parameter                   :: pi=3.141592
real(8)                :: r_e, e_, L, X_A, Y_A, X_C, Y_C, Y_H 

drag_coeff_vel = 6d0*pi*r_fiber
drag_coeff_omega= 8d0*pi*r_fiber**3d0
Spheroid_F =0;
Spheroid_T =0;
hinge1%r        = hinge2%X_i-hinge1%X_i

segment_length  = sqrt(dot_product(hinge1%r,hinge1%r))

middle_position=(hinge1%X_i+hinge2%X_i)/2d0
call set_velocity( middle_position,vel,omega,gamma_dot,epsilon_dot,flow_case )

r_e = 0.5 * segment_length / r_fiber

r_e = 1.14 * r_e **0.884  ! adjust aspect ratio

e_= (1-r_e**-2)**0.5;
L = log((1+e_)/(1-e_));
d = abs(hinge1%r) / segment_length

! for the force we need 2 Constants

Id = 0
forall(j = 1:3) Id(j,j) = 1

X_A = 8.0/3.0 * e_**3*(-2*e_ + (1+e_**2)*L)**-1
Y_A = 16.0/3.0 * e_**3*(2*e_ + (3*e_**2-1)*L)**-1

! print *, outerProd(d,d)
A = 6.0*pi*viscosity * segment_length/2 * (X_A* outerProd(d,d) + Y_A * (Id - outerProd(d,d) ))

!Spheroid_F = matmul(A , (vel -  hinge1%v_i - cross(hinge1%omega, hinge1%r /2 ) ) )

Spheroid_F = matmul(A , (vel ) )


! for the Torque we need 3 Constants
X_C = 4.0/3.0 * e_**3*(1-e_**2)*(2*e_ - (1-e_**2)*L)**-1
Y_C = 4.0/3.0 * e_**3*(2-e_**2)*(-2*e_ + (1+e_**2)*L)**-1
Y_H = 4.0/3.0 * e_**5*(-2*e_ + (1+e_**2)*L)**-1
 
C = 8.0*pi*viscosity * (segment_length/2)**3 * ( X_C * outerProd(d,d) + Y_C * (Id - outerProd(d,d)))

!print *, ' segment_length ', segment_length
!print *, ' r_e ', r_e, ' e ', e_, ' l ', l
!print *, ' d ', d
!print *, ' r_e ', r_e
!print *, ' x_a ', x_a , ' y_a ',y_a, ' x_c ',x_c, ' y_c ',y_c , ' y_h ' , y_h
!print *, ' c ', c
!print *, ' a ', a
!for a shear rate the rate of strain tensor is:

E_oo = 0
E_oo(1,2)= 0.5 *gamma_dot
E_oo(2,1)= 0.5 *gamma_dot
! generate permutation tensor

do i=1,3
    do j=1,3
        do k =1,3
            if( (i .eq. j) .or. (j.eq.k) .or. (k .eq. i)) then
               eps(i,j,k) =0
            else
                aa = (k-j) + (j-i);
                if ( (aa.eq.2) .or. (aa.eq.-1)) then
                    eps(i,j,k) = 1
                else
                    eps(i,j,k) = -1
                endif              
            endif          
        enddo
    enddo
enddo

part = 0

do i = 1 ,3
    part(i) =0
    do j =1 ,3
        if(i .eq. j) then
            exit
        endif
        do l_ =1 , 3
            if (l_.eq.j .or. l_.eq.i) then
                exit
            endif
            do k =1, 3
              part(i) = part(i) +  eps(i,j,l_)*d(l_)*d(k)*E_oo(j,k)                
            enddo
        enddo
    enddo
enddo

!print *, ' part ' , part
! Spheroid_T =  matmul(C , (omega - hinge1%omega)) -8*pi*viscosity* segment_length**3 * Y_H *part  

 Spheroid_T =  matmul(C , (omega )) -8*pi*viscosity* (segment_length/2)**3 * Y_H *part     

end subroutine computeTandF_Spheroid
                          
!====================================================================                    
subroutine ComputeOutputFandT( fibers, hinges, r_fiber, viscosity, gamma_dot, epsilon_dot, flow_case )
type (rod), allocatable, dimension(:)  :: hinges
type (fiber), allocatable, dimension(:):: fibers
integer(8)                             :: i, j, k,frame, first, last, flow_case
logical                                :: printVelocities
real(8)                                :: r_fiber, viscosity, gamma_dot, epsilon_dot
real(8), dimension(3)  :: bead_F, bead_T, spheroid_F, spheroid_T

! Compute hydrodynamic forces with the beads
do i=1, ubound(fibers,1)
    first=fibers(i)%first_hinge
	last =fibers(i)%first_hinge+fibers(i)%nbr_hinges-2 !Changed JULY 31 2012
	do j=first, last
		!print *,"indices",i, j
		call computeTandF_bead( hinges(j),&
                                hinges(j+1),&
                                r_fiber,&
                                viscosity,&
                                gamma_dot,&
                                epsilon_dot,&
                                flow_case, &
                                bead_F,&
                                bead_T )

	end do
end do

!print *, 'bead_F ' , bead_F(1), bead_F(2), bead_F(3)
!print *, 'bead_T ' , bead_T(1), bead_T(2), bead_T(3)
!Compute hydrodynamic forces with Resistance tensors

do i=1, ubound(fibers,1)
    first=fibers(i)%first_hinge
	last =fibers(i)%first_hinge+fibers(i)%nbr_hinges-2 !Changed JULY 31 2012
	do j=first, last
        
        call computeTandF_Spheroid( hinges(j),&
                                    hinges(j+1),&
                                    r_fiber,&
                                    viscosity,&
                                    gamma_dot,&
                                    epsilon_dot,&
                                    flow_case, &
                                    spheroid_F,&
                                    spheroid_T )
	end do
end do

write (6,*), bead_F(1), bead_F(2), bead_F(3)
write (6,*), bead_T(1), bead_T(2), bead_T(3)
write (6,*), spheroid_F(1), spheroid_F(2), spheroid_F(3)
write (6,*), spheroid_T(1), spheroid_T(2), spheroid_T(3)
write (6,*), hinges(1)%r(1), hinges(1)%r(2), hinges(1)%r(3)
write (6,*), hinges(1)%omega(1), hinges(1)%omega(2), hinges(1)%omega(3)
write (6,*), hinges(1)%omega_fluid_sum(1)/hinges(1)%nbr_beads, hinges(1)%omega_fluid_sum(2)/hinges(1)%nbr_beads, hinges(1)%omega_fluid_sum(3)/hinges(1)%nbr_beads

!print *, 'bead_F ' , bead_F(1), bead_F(2), bead_F(3)
!print *, 'bead_T ' , bead_T(1), bead_T(2), bead_T(3)
!print *, 'spheroid_F ' , spheroid_F(1), spheroid_F(2), spheroid_F(3)
!print *, 'spheroid_T ' , spheroid_T(1), spheroid_T(2), spheroid_T(3)

end subroutine ComputeOutputFandT


end module m_FiberCalc             !2018/07/21 change name
