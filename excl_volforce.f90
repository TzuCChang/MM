!====================================================================
module m_ExclVolForces  !2018/08/01  change name
    
use m_DataStructures
use m_UtilityLib

use omp_lib
implicit none
contains

subroutine excl_VolForceMomentsTotal( fibers,&
                                      hinges,&
                                      ghost_segments,&
                                      r_fiber,&
                                      ex_vol_const,&
                                      nbr_neighbors,&
                                      neighbor_list,&
                                      distance_neighbors,&
                                      fric_coeff,&
                                      box_size,&
                                      distanceFactor,&
                                      gamma_dot,&
                                      time )

logical                                 :: is_for_neigh_list, flag
type(fiber), dimension(:)               :: fibers
type(rod),   dimension(:)               :: hinges
type(rod)                               :: ghost_hingeA, ghost_hingeB
type(segment), dimension(:), allocatable   :: ghost_segments

integer(8)                              :: i, j, k, l, m, nbr_neighbors
integer(8), dimension(:,:), allocatable :: neighbor_list
real(8), dimension(:,:), allocatable    :: distance_neighbors
real(8), dimension(3)                   :: box_size, Gab
real(8)                                 :: gamma_dot, viscosity, ex_vol_const, distanceFactor
real(8)                                 :: r_fiber, s, t, Gab_min, fric_coeff, time, dX


is_for_neigh_list= .false.

dX= mod( gamma_dot*box_size(2)/2*time, box_size(1) )

call excl_VolGhostSegments( hinges, ghost_segments, box_size, dX )

do i=1, ubound( fibers, 1 )
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        
		m= 0
		flag= .true.
        
		do while ( flag.eqv..true. )
            
			m= m+1

			l= neighbor_list(j,m)
            
			if( l.ne. 0 ) then
                
			    ghost_hingeA%X_i= ghost_segments(l)%A
			    ghost_hingeB%X_i= ghost_segments(l)%B
                
				call dist_segments( ghost_hingeA%X_i,&
	                                ghost_hingeB%X_i,&
	                                hinges(j)%X_i,&
	                                hinges(j+1)%X_i,&
                                    s,&
                                    t,&
                                    Gab,&
		                            Gab_min )                
                
				call excl_VolForceSegments( ghost_hingeA,&
                                            ghost_hingeB,&
								            hinges(j),&
								            hinges(j+1),&
							                ex_vol_const,&
                                            r_fiber,&
                                            s,&
                                            t,&
                                            Gab,&
                                            Gab_min,&
			                                is_for_neigh_list,&
			                                fric_coeff )
			end if
			
			if ((m.eq.nbr_neighbors).or.(neighbor_list(j,m).eq. 0)) then
				flag= .false.
			end if
		end do			
	end do
end do
 
end subroutine excl_VolForceMomentsTotal

                                        
subroutine excl_VolForceMomentsWalls2( fibers,&    !2018/08/01 change name
                                       hinges,&
                                       r_fiber,&
                                       ex_vol_const,&
                                       box_size,&
                                       fric_coeff,&
                                       is_fric_wall,&
				                       gamma_dot )  
                                        
integer                  :: i, j, k
type(fiber), dimension(:):: fibers
type(rod), dimension(:)  :: hinges
real(8)                  :: ex_vol_const, dist, r_fiber, vec_to_wall, threshold
real(8), dimension(3)    :: box_size, Exc_Vol_Force_Partial, r, v_rel
real(8), dimension(2)    :: wall_position
real(8)                  :: FAC, gamma_dot
logical                  :: is_fric_wall
real(8)                  :: fric_coeff

threshold= 1*r_fiber

wall_position(1)=  box_size(2)/2
wall_position(2)= -box_size(2)/2  !wall_position(1)

do k= 1,2
    !if (k.eq.1) then
    if (k.eq.2) then
        FAC= +1;
    else
        FAC= -1;
    endif

    do i=1, ubound (fibers,1)
    
        j= fibers(i)%first_hinge
        
        Exc_Vol_Force_Partial= 0;
        
        vec_to_wall=  hinges(j)%X_i(2) - wall_position(k)
        
        dist= abs(vec_to_wall)
        
        if (dist .le. threshold) then
            
            r=hinges(j+1)%X_i-hinges(j)%X_i
            Exc_Vol_Force_Partial(2) =FAC*dist*ex_vol_const*exp(-2*(dist-1))
            v_rel=0
            
            if (hinges(j)%X_i(2).gt.0) then
                v_rel(1)= gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
            else
                v_rel(1)= -gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
            end if 
            
            v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))
        
            if ((v_rel(1) .eq. v_rel(1)).and.(v_rel(2) .eq. v_rel(2)) .and.(v_rel(3) .eq. v_rel(3)).and. is_fric_wall ) then !Check for NaN. If NaN, vrel
                 
                Exc_Vol_Force_Partial  = Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
            end if
            
            hinges(j)%F_Excl_Vol   = hinges(j)%F_Excl_Vol+ Exc_Vol_Force_Partial
#ifdef TENSOR            
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(-0.5*r, Exc_Vol_Force_Partial)
#else
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(r, Exc_Vol_Force_Partial)
#endif
         end if

        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
            
            Exc_Vol_Force_Partial=0;
            vec_to_wall= hinges(j+1)%X_i(2)-wall_position(k)
            dist=abs(vec_to_wall)
	        r=hinges(j+1)%X_i-hinges(j)%X_i
            
	        if (dist .le. threshold) then
	            v_rel=0
                if (hinges(j)%X_i(2).gt.0) then
                    v_rel(1)= gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
                else
                    v_rel(1)= -gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
                end if 
                v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))
                Exc_Vol_Force_Partial(2) =FAC*dist*ex_vol_const*exp(-2*(dist-1))
                
               
                if ((v_rel(1) .eq. v_rel(1)).and.(v_rel(2) .eq. v_rel(2)) .and.(v_rel(3) .eq. v_rel(3)).and. is_fric_wall ) then !Check for NaN. If NaN, vrel
                    Exc_Vol_Force_Partial  = Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
	            end if
	     
       		    hinges(j)%F_Excl_Vol   = hinges(j)%F_Excl_Vol+ Exc_Vol_Force_Partial
#ifdef TENSOR            
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(0.5*r, Exc_Vol_Force_Partial)
#else
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(r, Exc_Vol_Force_Partial)
#endif

	        end if  
        end do
    end do
end do

end subroutine excl_VolForceMomentsWalls2  !2018/08/01 change name
                                        
                                           
subroutine excl_VolForceSegments( hingesa1,&
                                  hingesa2,&
                                  hingesb1,&
                                  hingesb2,&
                                  ex_vol_const,&
                                  r_fiber,&
                                  s,&
                                  t,&
                                  Gab,&
                                  Gab_min,&
				                  is_for_neigh_list,&
				                  fric_coeff )

type(vec_arr)       , dimension(7)  :: vec
type(simple_vec_arr), dimension(7)  :: vec_a

type(rod)               :: hingesb1, hingesb2, hingesa1, hingesa2

logical                 :: are_near, is_for_neigh_list

integer(8)              :: i,j, k,l, ind1, ind2

real(8), dimension(3)   :: void, Exc_Vol_Force_Partial, v_rel,  r1, r2, Gab
real(8)                 :: threshold, r_fiber, ex_vol_const, fric_coeff, Sab, fac,  s, t, Gab_min
real(8), parameter      :: pi=3.141592


threshold= 2*r_fiber

r1 = (hingesa2%X_i - hingesa1%X_i) 
r2 = (hingesb2%X_i - hingesb1%X_i) 

!print *, " gab_norm1  ", gab_min

if( (is_for_neigh_list .eqv. .false.) .and.&
    (Gab_min .lt. threshold) .and. &
    (Gab_min .gt. 1e-8) ) then
   
	are_near= .true.
    
	!Exc_Vol_Force_Partial =-vec(k)%Gab*ex_vol_const*exp(-2*(vec(k)%Gab_norm/R-2))!CORRECT ONE
    !fac = ex_vol_const * pi * hingesb1%ave_viscosity * hingesb1%length * r * hingesb1%gamma_dot
    
    fac = ex_vol_const*2*r_fiber
    
    Exc_Vol_Force_Partial =-(Gab/(Gab_min))*fac*exp(-2*(Gab_min/r_fiber-2))
      
	if (hingesb1%is_stationary==1) hingesb1%v_i=0
	if (hingesa1%is_stationary==1) hingesa1%v_i=0
	
    v_rel= -hingesb1%v_i+hingesa1%v_i -cross(hingesb1%omega,r2*t)&
                                      +cross(hingesa1%omega,r1*s)
    
    if ((v_rel(1) .eq. 0.0) .and. (v_rel(2) .eq. 0.0) .and. (v_rel(3) .eq. 0.0) ) then
        v_rel =0                                        !error 2018/0714 �ץ�
        !v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))   !error 2018/0714 �ץ�
    else
        !v_rel =0                                       !error 2018/0714 �ץ�
        v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))    !error 2018/0714 �ץ�
    end if
    
    if ((v_rel(1) .eq. v_rel(1)).and.(v_rel(2) .eq. v_rel(2)) .and.(v_rel(3) .eq. v_rel(3)) ) then !Check for NaN. If NaN, vrel
        Exc_Vol_Force_Partial  = Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
	end if

	hingesb1%F_Excl_Vol(1)= hingesb1%F_Excl_Vol(1)+ Exc_Vol_Force_Partial(1)
    hingesb1%F_Excl_Vol(2)= hingesb1%F_Excl_Vol(2)+ Exc_Vol_Force_Partial(2)
    hingesb1%F_Excl_Vol(3)= hingesb1%F_Excl_Vol(3)+ Exc_Vol_Force_Partial(3)
    
#ifdef TENSOR    
    void = cross((r2*t -r2/2 ), Exc_Vol_Force_Partial)
#else
    void = cross((r2*t  ), Exc_Vol_Force_Partial)
#endif

	hingesb1%T_Excl_Vol(1)= hingesb1%T_Excl_Vol(1) + void(1)
    hingesb1%T_Excl_Vol(2)= hingesb1%T_Excl_Vol(2) + void(2)
    hingesb1%T_Excl_Vol(3)= hingesb1%T_Excl_Vol(3) + void(3)
    
    !print *, "collision happened 2"
    
    end if
    
end subroutine excl_VolForceSegments                                         
                                                                                     

subroutine excl_VolGhostSegments( hinges,&
                                  ghost_segments,&
                                  box_size,&
                                  dX )

type(rod),     dimension(:)              :: hinges
type(segment), dimension(:), allocatable :: ghost_segments

integer(8)                               :: i, j,k, l, m
real(8)                                  :: dX
real(8),dimension(3)                     :: box_size

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO PRIVATE (i, j, k)

do i=1, ubound( ghost_segments, 1 )
!print *,"i", i   ," of ",   ubound(ghost_segments,1)
    j= ghost_segments(i)%orig_pos(1)
    k= ghost_segments(i)%orig_pos(2)

    select case (ghost_segments(i)%axis_loc)
        case(1) ! image +X
            ghost_segments(i)%A= hinges(j)%X_i
            ghost_segments(i)%B= hinges(k)%X_i
        case(2) ! image +X
            ghost_segments(i)%A(1)= hinges(j)%X_i(1)+box_size(1)
            ghost_segments(i)%B(1)= hinges(k)%X_i(1)+box_size(1)
        case(3) ! image -X
            ghost_segments(i)%A(1)= hinges(j)%X_i(1)-box_size(1)
            ghost_segments(i)%B(1)= hinges(k)%X_i(1)-box_size(1)
        case(4) ! image +Z
            ghost_segments(i)%A(3)= hinges(j)%X_i(3)+box_size(3)
            ghost_segments(i)%B(3)= hinges(k)%X_i(3)+box_size(3)
        case(5) ! image -Z
            ghost_segments(i)%A(3)= hinges(j)%X_i(3)-box_size(3)
            ghost_segments(i)%B(3)= hinges(k)%X_i(3)-box_size(3)
        case(6) ! image +Ya
            ghost_segments(i)%A(1)= hinges(j)%X_i(1)+dX
            ghost_segments(i)%B(1)= hinges(k)%X_i(1)+dX
            ghost_segments(i)%A(2)= hinges(j)%X_i(2)+box_size(2)
            ghost_segments(i)%B(2)= hinges(k)%X_i(2)+box_size(2)
        case(7) ! image +Yb
            ghost_segments(i)%A(1)= hinges(j)%X_i(1)+dX - box_size(1)
            ghost_segments(i)%B(1)= hinges(k)%X_i(1)+dX - box_size(1)
            ghost_segments(i)%A(2)= hinges(j)%X_i(2)+box_size(2)
            ghost_segments(i)%B(2)= hinges(k)%X_i(2)+box_size(2)
        case(8) ! image -Ya
            ghost_segments(i)%A(1)= hinges(j)%X_i(1)-dX
            ghost_segments(i)%B(1)= hinges(k)%X_i(1)-dX
            ghost_segments(i)%A(2)= hinges(j)%X_i(2)-box_size(2)
            ghost_segments(i)%B(2)= hinges(k)%X_i(2)-box_size(2)
        case(9) ! image -Yb
            ghost_segments(i)%A(1)= hinges(j)%X_i(1)- dX + box_size(1)
            ghost_segments(i)%B(1)= hinges(k)%X_i(1)- dX + box_size(1)
            ghost_segments(i)%A(2)= hinges(j)%X_i(2)-box_size(2)
            ghost_segments(i)%B(2)= hinges(k)%X_i(2)-box_size(2)    
            
        end select
end do

!$OMP END DO !NOWAIT
!$OMP DO PRIVATE (i)

do i=1, ubound(hinges,1)
	hinges(i)%T_Excl_Vol=0d0
	hinges(i)%F_Excl_Vol=0d0
end do

end subroutine excl_VolGhostSegments 

end module m_ExclVolForces    !2018/08/01  Add

!********************************************************************