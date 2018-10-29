!====================================================================
module m_ExclVolForces  !2018/08/01  change name

!subroutine: excl_VolForceMomentsTotal  - fiber and fiber  interaction
!subroutine: excl_VolForceMomentsWalls2 - fiber and surrounding interaction
!subroutine: excl_VolForceSegments      - intra-fiber forces 
!subroutine: excl_VolGhostSegments

use m_DataStructures
use m_UtilityLib

use omp_lib
implicit none
contains

subroutine excl_VolForceMomentsTotal( fibers,&   !2018/08/03  修正  is_for_neigh_list= .false.
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

logical                                 :: flag
type(fiber), dimension(:)               :: fibers
type(rod),   dimension(:)               :: hinges
type(rod)                               :: ghost_hingeA, ghost_hingeB, ghost_hingeC, ghost_hingeD
type(segment), dimension(:), allocatable   :: ghost_segments

integer(8)                              :: i, j, k, m, j1, j2, nbr_neighbors   !2018/08/04  add
integer(8), dimension(:,:), allocatable :: neighbor_list
real(8), dimension(:,:), allocatable    :: distance_neighbors
real(8), dimension(3)                   :: box_size, Gab
real(8)                                 :: gamma_dot, viscosity, ex_vol_const, distanceFactor
real(8)                                 :: r_fiber, s, t, Gab_min, fric_coeff, time, dX, threshold

threshold= 2*r_fiber !兩個纖維能夠靠近最近的距離

dX= mod( gamma_dot*box_size(2)/2*time, box_size(1) )
! mod(a,b) = Remainder function of a/b
call excl_VolGhostSegments( hinges, ghost_segments, box_size, dX )

do i=1, ubound( fibers, 1 )
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        
		m= 0    !m=什麼??
		flag= .true.
        
		do while ( flag.eqv..true. )
            
			m= m+1

			k= neighbor_list(j,m)   !2018/08/04  將 l 改成 k
            
			if( k .ne. 0 ) then !.ne.= /=
                
                j1= ghost_segments(k)%orig_pos(1)   !2018/08/04  add ghost_segment繼承j1
                j2= ghost_segments(k)%orig_pos(2)   !2018/08/04  add ghost_segment繼承j2
                
                ghost_hingeA= hinges(j1)            !2018/08/04  add  複製 (j1) ghost_hingeA繼承給hinge(j1) 因為ghost_segment(k)為頭尾兩端
                ghost_hingeB= hinges(j2)            !2018/08/04  add  複製 (j2)
                
			    ghost_hingeA%X_i= ghost_segments(k)%A  !2018/08/04 因為前面 call excl_VolGhostSegments(hinges,ghost_segments,box_size,dX)
			    ghost_hingeB%X_i= ghost_segments(k)%B  !2018/08/04 目的是將其他周圍和本身的cell共9個納進來計算,同時每個 ghost的位置不一樣
               
!print *,"k,j1,j2 ",k,j1,j2                
!print *,"A ", ghost_hingeA%is_stationary
!print *,"A ", ghost_hingeA%X_i
!print *,"A ", ghost_hingeA%v_i
!print *,"A ", ghost_hingeA%omega

!print *,"B ", ghost_hingeB%is_stationary
!print *,"B ", ghost_hingeB%X_i
!print *,"B ", ghost_hingeB%v_i
!print *,"B ", ghost_hingeB%omega              

!print *,"C ", ghost_hingeC%is_stationary
!print *,"C ", ghost_hingeC%X_i
!print *,"C ", ghost_hingeC%v_i
!print *,"C ", ghost_hingeC%omega

!print *,"D ", ghost_hingeD%is_stationary
!print *,"D ", ghost_hingeD%X_i
!print *,"D ", ghost_hingeD%v_i
!print *,"D ", ghost_hingeD%omega
!pause


				call dist_segments( ghost_hingeA%X_i,&
	                                ghost_hingeB%X_i,&
	                                hinges(j  )%X_i,&
	                                hinges(j+1)%X_i,&
                                    s,&
                                    t,&
                                    Gab,&
		                            Gab_min )
                
                if( (Gab_min .lt. threshold) .and. (Gab_min .gt. 1e-8) ) then !2018/08/03 修正 !gab min < threshold便進來呼叫計算force (碰到才計算)
                    
!print *,"k,j1,j2 ",k,j1,j2     
!print *,"A ", ghost_hingeC%is_stationary
!print *,"A ", ghost_hingeC%X_i
!print *,"A ", ghost_hingeC%v_i
!print *,"A ", ghost_hingeC%omega
!print *,"B ", ghost_hingeD%is_stationary
!print *,"B ", ghost_hingeD%X_i
!print *,"B ", ghost_hingeD%v_i
!print *,"B ", ghost_hingeD%omega
!pause
!contact force>> 利用excl_volforcesegment計算contact force, neighbor ghost進來比
				     call excl_VolForceSegments( ghost_hingeA,&                !2018/08/04  add
                                                 ghost_hingeB,&                !2018/08/04  add
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

!纖維與牆壁之間的作用力, 有在main呼叫 (main裡呼叫total跟wall2)                                          
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
wall_position(1)= X_fiber(2) + box_size(2)/2
wall_position(2)= X_fiber(2) - box_size(2)/2  !wall_position(1)

do k= 1,2
    !if (k.eq.1) then
    if( k .eq. 2 ) then
        FAC= +1;    !此處的FAC跟hinge_damage.f90的fac有什麼關聯?
    else
        FAC= -1;
    endif

    do i=1, ubound (fibers,1)
    
        j= fibers(i)%first_hinge
        
        Exc_Vol_Force_Partial= 0;
        
        vec_to_wall= hinges(j)%X_i(2) - wall_position(k)  !error 2018/08/04 hinges到Wall的距離,這是錯誤的,
        
        dist= abs(vec_to_wall)   !2018/08/04 error 2018/08/04 hinges到Wall距離的算法是錯誤的,
                                 !2018/08/04 除非輸入的所有fibers或hinges的y座標,全部定在-box_size(2)/2和+box_size(2)/2之間才適用
                                 !2018/08/04 box下方,因為wall_position(2)= -box_size(2)/2= -2e-4, 因此當k=2時, dist永遠大於threshold=5e-6,2.都不會碰到Wall
!print *,"@@@( k= ",k
!print *,"@@@( i= ",i
!print *,"@@@( j= ",j
!print *,"@@@( wall=      ", wall_position(k)
!print *,"@@@( X(2)=      ", hinges(j)%X_i(2)
!print *,"@@@( dist=      ", dist
!print *,"@@@( threshold= ", threshold
!print *,"***("
!pause         
        
        if (dist .le. threshold) then !先算速度 在算Torque

!print *,"@@@( k= ",k
!print *,"@@@( i= ",i
!print *,"@@@( j= ",j
!print *,"@@@( wall=      ", wall_position(k)
!print *,"@@@( X(2)=      ", hinges(j)%X_i(2)
!print *,"@@@( dist=      ", dist
!print *,"@@@( threshold= ", threshold       
!print *,"***("
!pause
            r=hinges(j+1)%X_i-hinges(j)%X_i
            Exc_Vol_Force_Partial(2) =FAC*dist*ex_vol_const*exp(-2*(dist-1))
            v_rel=0            !Error 2018/08/04  只有y方向速度設定為零即可 v_rel(2)= 0
            
            if (hinges(j)%X_i(2).gt.0) then
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
            dist= abs(vec_to_wall)
	        r= hinges(j+1)%X_i - hinges(j)%X_i

!print *,"@@@( k= ",k
!print *,"@@@( i= ",i
!print *,"@@@( j= ",j
!print *,"@@@( wall=      ", wall_position(k)
!print *,"@@@( X(2)=      ", hinges(j)%X_i(2)
!print *,"@@@( dist=      ", dist
!print *,"@@@( threshold= ", threshold       
!print *,"***("
!pause            
	        if( dist .le. threshold ) then

!print *,"@@@( k= ",k
!print *,"@@@( i= ",i
!print *,"@@@( j= ",j
!print *,"@@@( wall=      ", wall_position(k)
!print *,"@@@( X(2)=      ", hinges(j)%X_i(2)
!print *,"@@@( dist=      ", dist
!print *,"@@@( threshold= ", threshold       
!print *,"***("
!pause     

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

end subroutine excl_VolForceMomentsWalls2  !2018/08/01 change name
                                        
                                           
subroutine excl_VolForceSegments( hingesa1,&  !2018/08/03 修正
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

!print *,"a1s", hingesa1%is_stationary
!print *,"a1X", hingesa1%X_i
!print *,"a1v", hingesa1%v_i
!print *,"a1w", hingesa1%omega
!print *,"a2s", hingesa2%is_stationary
!print *,"a2X", hingesa2%X_i
!print *,"a2v", hingesa2%v_i
!print *,"a2w", hingesa2%omega

!print *,"b1s", hingesb1%is_stationary
!print *,"b1X", hingesb1%X_i
!print *,"b1v", hingesb1%v_i
!print *,"b1w", hingesb1%omega
!print *,"b2s", hingesb2%is_stationary
!print *,"b2X", hingesb2%X_i
!print *,"b2v", hingesb2%v_i
!print *,"b2w", hingesb2%omega
!pause
    
	if( hingesb1%is_stationary==1 ) hingesb1%v_i= 0
	if( hingesa1%is_stationary==1 ) hingesa1%v_i= 0  !error/double check 2018/08/03  因為hingesa1和hingesa2均來自segment 沒有is_stationary資訊
	!v_rel=兩個相減為相對速度和旋轉相對速度的差	
    v_rel= -hingesb1%v_i + hingesa1%v_i - cross(hingesb1%omega,r2*t)&
                                        + cross(hingesa1%omega,r1*s)    !error/double check 2018/08/03  因為hingesa1和hingesa2均來自segment 沒有omega資訊
!print *,"vel",v_rel
!pause

    if( (v_rel(1) .eq. 0.0) .and. (v_rel(2) .eq. 0.0) .and. (v_rel(3) .eq. 0.0) ) then
        v_rel= 0                                        !error 2018/0714 修正
        !v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))   !error 2018/0714 修正
    else
        !v_rel =0                                       !error 2018/0714 修正
        v_rel= v_rel/(sqrt(dot_product(v_rel,v_rel)))   !error 2018/0714 修正
    end if

    Exc_Vol_Force_Partial= Exc_Vol_Force_Partial + v_rel*abs(fric_coeff*Exc_Vol_Force_Partial) !2018/08/04 拿掉 if 判斷 !得到force(正向力)
    
    hingesb1%F_Excl_Vol= hingesb1%F_Excl_Vol + Exc_Vol_Force_Partial !2018/08/04 修正
    
#ifdef TENSOR    
    void= cross( (r2*t -r2/2 ), Exc_Vol_Force_Partial )
#else
    void= cross( (r2*t       ), Exc_Vol_Force_Partial )
#endif

    hingesb1%T_Excl_Vol= hingesb1%T_Excl_Vol + void                  !2018/08/04 修正

    !print *, "collision happened 2"
      
end subroutine excl_VolForceSegments                                         
                                                                                     

subroutine excl_VolGhostSegments( hinges, ghost_segments, box_size, dX )

type(rod),     dimension(:)              :: hinges
type(segment), dimension(:), allocatable :: ghost_segments

integer(8)                               :: i, j, k
real(8)                                  :: dX
real(8),dimension(3)                     :: box_size

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO PRIVATE (i, j, k)

do i= 1, ubound( ghost_segments, 1 ) !從第一個segment做到最後一個segment
   j= ghost_segments(i)%orig_pos(1) !首
   k= ghost_segments(i)%orig_pos(2) !尾

   select case ( ghost_segments(i)%axis_loc )
            
        !This is the original segment
        case(1)  !2018/08/05  changet
        ghost_segments(i)%A= hinges(j)%X_i
        ghost_segments(i)%B= hinges(k)%X_i

        !Image segment in +x
        case(2)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1) + box_size(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1) + box_size(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3)
        
        ! Image segment in -x
        case(3)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1) - box_size(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1) - box_size(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3)

        ! Image segment in +z
        case(4)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3) + box_size(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3) + box_size(3)

        ! Image segment in -z
        case(5)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3) - box_size(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3) - box_size(3)
        

        ! Image segment in +x+z
        case(6)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1) + box_size(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3) + box_size(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1) + box_size(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3) + box_size(3)
        
        ! Image segment in -x+z
        case(7)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1) - box_size(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3) + box_size(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1) - box_size(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3) + box_size(3)

        ! Image segment in -x-z
        case(8)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1) - box_size(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3) - box_size(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1) - box_size(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3) - box_size(3)

        ! Image segment in +x-z
        case(9)  !2018/08/05  changet
        ghost_segments(i)%A(1)= hinges(j)%X_i(1) + box_size(1)
        ghost_segments(i)%A(2)= hinges(j)%X_i(2)      
        ghost_segments(i)%A(3)= hinges(j)%X_i(3) - box_size(3)
        ghost_segments(i)%B(1)= hinges(k)%X_i(1) + box_size(1)
        ghost_segments(i)%B(2)= hinges(k)%X_i(2)        
        ghost_segments(i)%B(3)= hinges(k)%X_i(3) - box_size(3)            
            
   end select
end do

!$OMP END DO !NOWAIT
!$OMP DO PRIVATE (i)

do i= 1, ubound(hinges,1)
	hinges(i)%T_Excl_Vol= 0d0
	hinges(i)%F_Excl_Vol= 0d0
end do

end subroutine excl_VolGhostSegments 

end module m_ExclVolForces    !2018/08/01  Add

!********************************************************************