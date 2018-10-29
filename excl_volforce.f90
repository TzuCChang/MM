!====================================================================
module m_ExclVolForces  !2018/08/01  change name
    
use m_DataStructures
use m_UtilityLib

use omp_lib
implicit none
    contains

!subroutine: excl_VolForceMomentsTotal  - fiber and fiber  interaction
!subroutine: excl_VolForceMomentsWalls2 - fiber and surrounding interaction
!subroutine: excl_VolForceSegments      - intra-fiber forces 
!subroutine: excl_VolGhostSegments
    
subroutine excl_VolForceMomentsTotal( fibers,&        !2018/09/08  changed 
                                      hinges,&
                                      ghost_segments,&
                                      neighbor_list,&
                                      nbr_neighbors,&    
                                      r_fiber,&
                                      fric_coeff,&
                                      ex_vol_const) 


logical                                 :: flag
type(fiber), dimension(:)               :: fibers
type(rod),   dimension(:)               :: hinges
type(rod)                               :: ghost_hingeA, ghost_hingeB, ghost_hingeC, ghost_hingeD
type(segment), dimension(:), allocatable   :: ghost_segments

integer(8)                              :: i, j, k, m, j1, j2, nbr_neighbors   !2018/08/04  add
integer(8), dimension(:,:), allocatable :: neighbor_list
real(8), dimension(3)                   :: Gab
real(8)                                 :: viscosity, ex_vol_const
real(8)                                 :: r_fiber, s, t, Gab_min, fric_coeff, threshold

threshold= 2*r_fiber    !the smallest distance between two fibers


do i=1, ubound( fibers, 1 )
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        
		m= 0
		flag= .true.
        
		do while ( flag.eqv..true. )
            
			m= m+1

			k= neighbor_list(j,m)   !2018/08/04  replace "l" with "k"
            
			if( k .ne. 0 ) then
                
                j1= ghost_segments(k)%orig_pos(1)   !2018/08/04  add ghost_segment (inherit j1)
                j2= ghost_segments(k)%orig_pos(2)   !2018/08/04  add ghost_segment (inherit j2)
                
                ghost_hingeA= hinges(j1)            !2018/08/04  add copy information from (j1): hinge(j1) inherit ghost_hingeA (ghost_segment(k) includes start and end point)
                ghost_hingeB= hinges(j2)            !2018/08/04  add copy information from (j2)
                
			    ghost_hingeA%X_i= ghost_segments(k)%A  !2018/08/04 call excl_VolGhostSegments(hinges,ghost_segments,box_size,dX)
			    ghost_hingeB%X_i= ghost_segments(k)%B  !2018/08/04 to calculate including neighbor cells and itself (total =9) (every ghost position is different)

				call dist_segments( ghost_hingeA%X_i,&
	                                ghost_hingeB%X_i,&
	                                hinges(j  )%X_i,&
	                                hinges(j+1)%X_i,&
                                    s,&
                                    t,&
                                    Gab,&
		                            Gab_min )
                
                if( (Gab_min .lt. threshold) .and. (Gab_min .gt. 1e-8) ) then   !2018/08/03 correct !if gab min < threshold, calculate force (calculate after they touch)
                                                                                !contact force: use excl_volforcesegment to calculate contact force, compare neighbor ghost
				     call excl_VolForceSegments( ghost_hingeA,&                 !2018/08/04  add
                                                 ghost_hingeB,&                 !2018/08/04  add
								                 hinges(j  ),&
                                                 hinges(j+1),&
							                     ex_vol_const,&
                                                 r_fiber,&
                                                 s,&
                                                 t,&
                                                 Gab,&
                                                 Gab_min,&
			                                     fric_coeff )
                end if
			end if
			
			if ((m.eq.nbr_neighbors).or.(neighbor_list(j,m).eq. 0)) then
				flag= .false.
			end if
		end do			
	end do
end do

end subroutine excl_VolForceMomentsTotal

!======================================================================    
!the interaction force between fibers and the wall (call in main.f90): main calls total and wall2   
                                      
subroutine excl_VolForceMomentsWalls2( fibers,&    !2018/08/01 change name
                                       hinges,&
                                       box_size,&
                                       is_fric_wall,&
                                       gamma_dot,&
                                       r_fiber,&
                                       fric_coeff,&
                                       ex_vol_const)
integer                  :: i, j, k
type(fiber), dimension(:):: fibers
type(rod), dimension(:)  :: hinges
real(8)                  :: ex_vol_const, dist, r_fiber, vec_to_wall, threshold
real(8), dimension(3)    :: box_size, Exc_Vol_Force_Partial, r, v_rel, X_fiber
real(8), dimension(2)    :: wall_position
real(8)                  :: FAC, gamma_dot
logical                  :: is_fric_wall
real(8)                  :: fric_coeff

!2018/08/05  Compute center of mass for all fibers X_fiber
k= 0
X_fiber= 0
do i= 1, ubound (fibers,1)  
do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
   X_fiber= X_fiber + ((hinges(j)%X_i+hinges(j+1)%X_i)/2d0)
   k= k+1
end do
end do        
X_fiber= X_fiber/real(k)  !2018/08/05  X_fiber= center of mass

threshold= 1*r_fiber
wall_position(1)= X_fiber(2) + box_size(2)/2  !2018/08/10  correct wall_position(1)
wall_position(2)= X_fiber(2) - box_size(2)/2  !2018/08/10  correct wall_position(1)

do k= 1,2
    !if (k.eq.1) then
    if( k .eq. 2 ) then
        FAC= +1; 
    else
        FAC= -1;
    endif

    do i=1, ubound (fibers,1)
    
        j= fibers(i)%first_hinge
        
        Exc_Vol_Force_Partial= 0;
        
        vec_to_wall= hinges(j)%X_i(2) - wall_position(k)  !error 2018/08/04 the distance between hinhes and the wall is wrong
        
if( (k .eq. 1) .and. (vec_to_wall .gt. 0) ) then  !2018/08/10 correct upper wall
     vec_to_wall= -0.1*r_fiber                    !2018/08/10
end if

if( (k .eq. 2) .and. (vec_to_wall .lt. 0) ) then  !2018/08/10 correct lower Wall
     vec_to_wall=  0.1*r_fiber                    !2018/08/10
end if
        
        dist= abs(vec_to_wall)   !2018/08/04 error 2018/08/04 the distance between hinhes and the wall is wrong in the past
                                 !2018/08/04 it only workks unless the y-position for every fibers and hinges are located within -box_size(2)/2 and +box_size(2)/2
                                 !2018/08/04 below the box, because wall_position(2)= -box_size(2)/2= -2e-4, so when k=2, dist is always greater than threshold=5e-6,2. and never touch the wall
        
        if (dist .le. threshold) then !calculter velocity first, and then torque

            r= hinges(j+1)%X_i - hinges(j)%X_i
            
            Exc_Vol_Force_Partial(2)= FAC*dist*ex_vol_const*exp(-2*(dist-1))
            
            v_rel= 0            !Error 2018/08/04  only v_rel(2)= 0 set to zero
            
            if( hinges(j)%X_i(2) .gt. 0 ) then
                v_rel(1)= gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
            else
                v_rel(1)= -gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
            end if 
            
            v_rel= v_rel/(sqrt(dot_product(v_rel,v_rel)))
        
            if( is_fric_wall ) then               !2018/08/04  Check for NaN. If NaN, vrel
                Exc_Vol_Force_Partial= Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
            end if
            
            hinges(j)%F_Excl_Vol= hinges(j)%F_Excl_Vol + Exc_Vol_Force_Partial
#ifdef TENSOR            
            hinges(j)%T_Excl_Vol= hinges(j)%T_Excl_Vol + cross(-0.5*r, Exc_Vol_Force_Partial)
#else
            hinges(j)%T_Excl_Vol= hinges(j)%T_Excl_Vol + cross(r, Exc_Vol_Force_Partial)
#endif
         end if

        do j= fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
            
            Exc_Vol_Force_Partial= 0;
            
            vec_to_wall= hinges(j+1)%X_i(2)-wall_position(k)
            
if( (k .eq. 1) .and. (vec_to_wall .gt. 0) ) then  !2018/08/10 correct exceeding upper wall
     vec_to_wall= -0.1*r_fiber                    !2018/08/10 
end if

if( (k .eq. 2) .and. (vec_to_wall .lt. 0) ) then  !2018/08/10 correct exceeding lower wall
     vec_to_wall=  0.1*r_fiber                    !2018/08/10
end if

            dist= abs(vec_to_wall)
                     
	        if( dist .le. threshold ) then
    
	            r= hinges(j+1)%X_i - hinges(j)%X_i
            
	            v_rel= 0
                
                if (hinges(j)%X_i(2).gt.0) then
                    v_rel(1)=  gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
                else
                    v_rel(1)= -gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
                end if 
                
                v_rel= v_rel/(sqrt(dot_product(v_rel,v_rel)))
                
                Exc_Vol_Force_Partial(2)= FAC*dist*ex_vol_const*exp(-2*(dist-1))
               
                if( is_fric_wall ) then !2018/08/04  Check for NaN. If NaN, vrel
                    Exc_Vol_Force_Partial= Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
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

end subroutine excl_VolForceMomentsWalls2       !2018/08/01 change name
                                        
!======================================================================                                           
subroutine excl_VolForceSegments( hingesa1,&    !2018/08/03 changed
                                  hingesa2,&
                                  hingesb1,&
                                  hingesb2,&
                                  ex_vol_const,&
                                  r_fiber,&
                                  s,&
                                  t,&
                                  Gab,&
                                  Gab_min,&
				                  fric_coeff )

type(rod)               :: hingesa1, hingesa2, hingesb1, hingesb2

integer(8)              :: i,j, k,l, ind1, ind2

real(8), dimension(3)   :: void, Exc_Vol_Force_Partial, v_rel,  r1, r2, Gab
real(8)                 :: r_fiber, ex_vol_const, fric_coeff, Sab, fac,  s, t, Gab_min
real(8), parameter      :: pi=3.141592


    r1= hingesa2%X_i - hingesa1%X_i 
    r2= hingesb2%X_i - hingesb1%X_i 
    
	!Exc_Vol_Force_Partial =-vec(k)%Gab*ex_vol_const*exp(-2*(vec(k)%Gab_norm/R-2))!CORRECT ONE
    !fac = ex_vol_const * pi * hingesb1%ave_viscosity * hingesb1%length * r * hingesb1%gamma_dot
    
    fac= ex_vol_const*2*r_fiber
    
    Exc_Vol_Force_Partial= -(Gab/(Gab_min))*fac*exp(-2*(Gab_min/r_fiber-2))

    
	if( hingesb1%is_stationary==1 ) hingesb1%v_i= 0
	if( hingesa1%is_stationary==1 ) hingesa1%v_i= 0                     !error/double check 2018/08/03  both hingesa1 and hingesa2 comes from segment: lacks the info of "is_stationary"
	
    v_rel= -hingesb1%v_i + hingesa1%v_i - cross(hingesb1%omega,r2*t)& 	!v_rel = the diffence between the relative velocity and the relative angular velocity	
                                        + cross(hingesa1%omega,r1*s)    !error/double check 2018/08/03  the information of hingesa1 and hingesa2 are from the segment(no omega data)
!print *,"vel",v_rel
!pause

    if( (v_rel(1) .eq. 0.0) .and. (v_rel(2) .eq. 0.0) .and. (v_rel(3) .eq. 0.0) ) then
        v_rel= 0                                        !error 2018/0714 corrected
        !v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))   !error 2018/0714 corrected
    else
        !v_rel =0                                       !error 2018/0714 corrected
        v_rel= v_rel/(sqrt(dot_product(v_rel,v_rel)))   !error 2018/0714 corrected
    end if

    Exc_Vol_Force_Partial= Exc_Vol_Force_Partial + v_rel*abs(fric_coeff*Exc_Vol_Force_Partial) !2018/08/04 remove if, calculate the force(normal direction))
    
    hingesb1%F_Excl_Vol= hingesb1%F_Excl_Vol + Exc_Vol_Force_Partial !2018/08/04 corrected
    
#ifdef TENSOR    
    void= cross( (r2*t -r2/2 ), Exc_Vol_Force_Partial )
#else
    void= cross( (r2*t       ), Exc_Vol_Force_Partial )
#endif

    hingesb1%T_Excl_Vol= hingesb1%T_Excl_Vol + void                  !2018/08/04 corrected

    !print *, "collision happened 2"
      
end subroutine excl_VolForceSegments                                         
                                                                                     

end module m_ExclVolForces    !2018/08/01  Add
!======================================================================