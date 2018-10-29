!Units in SI (m, N, Pa, etc)
program cube_periodic

use m_DataStructures

use m_ReadData       !2018/07/21  change name
use m_OutputData     !2018/07/21  change name
use m_FiberCalc      !2018/07/21  change name
use m_FindNeighbors  !2018/07/21  change name
use m_FiberRegroup   !2018/07/21  change name
use m_HingesDamage   !2018/07/21  change name
use m_BendingTorque  !2018/07/21  change name
use m_Motion         !2018/07/21  change name
use m_UpDate         !2018/07/21  change name
use m_ExclVolForces  !2018/08/01  add   
use m_SimulationParameter  !2018/07/21  change name

use omp_lib

implicit none

real(8)                                :: E_Young, min_curv, r_fiber, viscosity, ex_vol_const,&
                                          gamma_dot, epsilon_dot, dt, inertia_moment, fric_coeff, distanceFactor,alpha
integer(8)                             :: nbr_neighbors, flow_case, nbr_intgr, writ_period, break_period,&
                                          i, j, k, n, frame, nbr_Fibers_OLD, nbr_Fibers_NEW, nbr_Fibers_INC
logical                                :: periodic_boundary, allow_breakage
real(8), dimension(3)                  :: box_size

integer, dimension(3)                   :: Nbr_bins

type(rod)  , allocatable, dimension(:) :: hinges
type(fiber), allocatable, dimension(:) :: fibers
                                           
real(8), parameter                     :: pi=3.141592

integer(8), dimension(:,:), allocatable:: neighbor_list

real(8), dimension(:,:), allocatable   :: distance_neighbors
logical                                :: recover_simulation
type(segment), dimension(:),allocatable:: ghost_segments
type(cell), dimension(:), allocatable  :: cells
real(8)                                :: start, finish, start2, finish2, timex, t
logical                                :: is_fric_wall, printVelocities 
type(simulationParameters)             :: simParameters
!*******************************************************************
! Default values

open(301,file='OUTPUT/OutputMessage.txt')

simParameters%IsPeriodicY =.false. 
print *, "Maximum number of threads" , omp_get_max_threads( )  
write(301,*), "Maximum number of threads" , omp_get_max_threads( ) 

#ifdef TENSOR            
print *, "Hydrodynamic representation being used is TENSOR"
write(301,*), "Hydrodynamic representation being used is TENSOR"

#else
print *, "Hydrodynamic representation being used is BEAD"
write(301,*), "Hydrodynamic representation being used is BEAD"

#endif

!$OMP PARALLEL
print *, "Number of threads being used" , omp_get_num_threads( )  
!$OMP END PARALLEL
!print *,"1"
 call  read_data( frame,&
                  recover_simulation,&
                  fric_coeff,&
                  is_fric_wall,&
                  E_Young,&
                  min_curv,&
                  r_fiber,&
                  viscosity,&
                  ex_vol_const,&
                  nbr_neighbors,&
                  gamma_dot,&
                  epsilon_dot,&
                  flow_case,&
                  periodic_boundary,&
                  box_size,&
                  dt,&
                  nbr_intgr,&
                  writ_period,&
                  allow_breakage,&
                  break_period,&
                  fibers,&
                  hinges,&
                  printVelocities,&
                  distanceFactor,&
                  simParameters )

call fiber_regroup_minmax_hinges(   fibers, hinges )                !2018/08/05 add
call fiber_regroup_minmax_segments( fibers, hinges )                !2018/08/05 add
call fiber_regroup_ShiftCenterToOrigion( fibers, hinges, box_size ) !2018/08/05 add

open(300,file='OUTPUT/end_to_end_distance.txt')


open(3,file='OUTPUT/positions.out')
open(5,file='OUTPUT/vels.out')
open(6,file='OUTPUT/forces.out')

Inertia_Moment=(pi/4.d0)*r_fiber**4d0

if(recover_simulation.eqv..true.) then 
	n=frame*writ_period+1
else 
	n=1
	frame =1
end if
start = OMP_get_wtime()

call simulation_parameter( hinges, r_fiber, gamma_dot, epsilon_dot, flow_case, simParameters ) !2018/07/15 修正字串

i= 0 
t= 0.d0                                                          !2018/07/14 增加
nbr_Fibers_INC= 1                                                !2018/07/14 修正
nbr_Fibers_NEW= ubound(fibers,1)                                 !2018/07/14 修正

nbr_Fibers_OLD= nbr_Fibers_NEW                                   !2018/07/14 修正

print *, "@@@@*( ", i, nbr_Fibers_NEW, nbr_Fibers_OLD            !2018/07/14 增加
write(301,*), "@@@@*( ", i, nbr_Fibers_NEW, nbr_Fibers_OLD    !2018/08/05 增加

call output_Length( t, fibers, hinges, frame, printVelocities )  !2018/07/14 增加

do i=n,  nbr_intgr
 
    t = dt*i                                                    !2018/07/14 修正

    !print *,"A t= ", t, i

    call cpu_time (start2)

    !print *, "Main 1"
 	if (MODULO(i,break_period)==0 .or. (i.eq.n) ) then 
 	   print *,"Integration", i, t
       write(301,*), "Integration", i, t

       call bending_torque_whole(fibers, hinges, E_Young, Inertia_Moment)

       if (allow_breakage)then
           call hinges_damage(fibers, hinges, min_curv, r_fiber)
       end if
       
 	   call fiber_regroup( fibers,&           !2018/07/21 change name
 	                       hinges,&
 	                       ghost_segments,&
                           E_Young,&          !ERROR 2018/07/07
                           Inertia_Moment,&   !ERROR 2018/07/07       
 	                       allow_breakage,&
 	                       min_curv,&
 	                       r_fiber,&
 	                       box_size,&
 	                       cells,&                                  
 	                       nbr_neighbors,&
                           neighbor_list,&
                           distance_neighbors,&
                           gamma_dot,&
                           t,&
                           Nbr_bins,&          !2018/07/21 add
                           distanceFactor,&
                           simParameters )
                                               !print *,"out config"
        !call cpu_time(finish)
        !print *, "Neighbors", finish-start
        
         call find_neighbors_new( fibers,&  !2018/07/21 change
                                  hinges,&
                                  ghost_segments,&
                                  nbr_neighbors,&
                                  neighbor_list,&
                                  distance_neighbors,&
                                  r_fiber,&
                                  cells,&
                                  Nbr_bins,&
                                  distanceFactor )
         
         call fiber_regroup_minmax_hinges(   fibers, hinges )     !2018/08/05 add
         call fiber_regroup_minmax_segments( fibers, hinges )     !2018/08/05 add

 	end if
 	
 	!do j=1, ubound (neighbor_list,1)
 	!    
 	!    print *, j, "Vecino", neighbor_list(j,:)
 	!end do
 	!stop
 	!print *,"Main 2"
	!print *,"5"
	
	!call cpu_time(start)
#ifdef TENSOR    
    call fiber_calc_tensor( fibers,&
                            hinges,&
                            r_fiber,&
                            viscosity,&
                            gamma_dot,&
                            epsilon_dot,&
                            flow_case,&
                            simparameters )
#else    
    call fiber_calc( fibers,&
                      hinges,&
                      r_fiber,&
                      viscosity,&
                      gamma_dot,&
                      epsilon_dot,&
                      flow_case )
#endif
    !call cpu_time(finish)
    !print *, "Fiber Parameters", finish-start                    
	!print *,"6"
	!do j=1,ubound(neighbor_list,1)
	!    print *,"Neighbor List", neighbor_list(j,:)
	!end do
	
    !call cpu_time(start)
 	call excl_VolForceMomentsTotal( fibers,&
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
                                    t )
    !call cpu_time(finish)
    !print *, "Interactions", finish-start
    !call cpu_time(start)
    
    if( .NOT. simParameters%IsPeriodicY ) then
        
         !call fiber_regroup_minmax_hinges( fibers, hinges )        !2018/08/05 add
         !call fiber_regroup_minmax_segments( fibers, hinges )      !2018/08/05 add
              
         call excl_VolForceMomentsWalls2( fibers,&    !2018/07/21 change name
                                          hinges,&
                                          r_fiber,&
                                          ex_vol_const,&
                                          box_size,&
                                          fric_coeff,&
                                          is_fric_wall,&
                                          gamma_dot )  
    end if
    
    !call cpu_time(finish)
    !print *, "Interactions walls", finish-start
    !call cpu_time(start)
    
 	call bending_torque_whole(fibers, hinges, E_Young, Inertia_Moment)
    
 	!call cpu_time(finish)
    !print *, "time for bending", finish-start
 	!call cpu_time(start)
    
    timex=0;
    call motion_fiber(fibers, hinges, r_fiber)
    
    !call cpu_time(finish)
    !print *, "Dealing with matrix ", finish-start
    !call cpu_time(start)
    
 	call update_periodic(fibers, hinges, dt, periodic_boundary, box_size, gamma_dot,dt* (i-1) )
    
    !call cpu_time(finish)
    !print *, "Update", finish-start
! 	if (MODULO(i,writ_period)==0 .or. (i .le. 5) ) then

 	if (MODULO(i,writ_period)==0 ) then
 		call output_data(t, fibers, hinges, frame,printVelocities)
 		frame=frame+1
    else 
        nbr_Fibers_NEW= ubound(fibers,1)                                         !2018/07/14 增加
        if ( nbr_Fibers_NEW .GT. (nbr_Fibers_OLD+nbr_Fibers_INC) ) then          !2018/07/14 增加
              print *, "@@@@@( ", i, nbr_Fibers_NEW, nbr_Fibers_OLD              !2018/07/14 增加
              call output_Length( t, fibers, hinges, frame, printVelocities )  !2018/07/14 增加   
              nbr_Fibers_OLD= nbr_Fibers_NEW                                     !2018/07/14 增加 
        end if
    end if
    
 	!call cpu_time (finish2)
 	!print *, "total incluye busqueda", finish2-start2	
    
end do 

call cpu_time(finish)
print *, "Time Elpased",  OMP_get_wtime()-start, "s"
write(301,*), "Time Elpased",  OMP_get_wtime()-start, "s"
close(300)
close(301)
close (3) 
close (5)
close (6)
write( *, * ) 'Press Enter to continue' 
read( *, * ) 

end program cube_periodic

! testing for github