!======================================================================
!Units in SI (m, N, Pa, etc)
program cube_periodic

use m_DataStructures

use m_ReadData       !2018/07/21  change name
use m_OutputData     !2018/07/21  change name
use m_FiberCalc      !2018/07/21  change name
use m_FindNeighbors  !2018/07/21  change name
use m_FindGhostSegments  !2018/09/08  add
use m_FiberRegroup   !2018/07/21  change name
use m_HingesDamage   !2018/07/21  change name
use m_BendingTorque  !2018/07/21  change name
use m_Motion         !2018/07/21  change name
use m_UpDate         !2018/07/21  change name
use m_ExclVolForces  !2018/08/01  add   
use m_SetVelocity    !2018/10/10  add 
use m_SimulationParameter  !2018/07/21  change name

use omp_lib

implicit none
type(rod),      dimension(:), allocatable  :: hinges
type(fiber),    dimension(:), allocatable  :: fibers
type(segment),  dimension(:), allocatable  :: ghost_segments
type(cell),     dimension(:), allocatable  :: cells
type(DynamicP), dimension(:), allocatable  :: flowcase_1848
type(simulationParameters)                 :: simParameters

logical                                   :: isOutputMessage

integer(8), dimension(:,:), allocatable   :: neighbor_list
integer(8), dimension(:),   allocatable   :: indexA                                                    !2018/08/12
integer(8)                                :: i, j, k, n, nbr_Fibers_OLD, nbr_Fibers_NEW, nbr_Fibers_INC
integer(8)                                :: iii, jjj, kkk, nbr_hinges                                 !2018/10/10 

real(8), dimension(:,:), allocatable      :: distance_neighbors, AA
real(8)                                   :: start, finish, start2, finish2, dX, t    !2018/10/05 by HAKAN
real(8), parameter                        :: pi=3.141592
!*******************************************************************

! Default values

open(300,file='OUTPUT/meanLength.txt')
open(301,file='OUTPUT/OutputMessage.txt')
open(302,file='OUTPUT/FiberLengthDistribution.txt')   !2018/08/12
open(303,file='OUTPUT/OrientationTensor.txt')         !2018/08/12
open(306,file='OUTPUT/a11.txt')                     !2018/10/07
open(307,file='OUTPUT/Ln.txt')                     !2018/10/07
open(3,  file='OUTPUT/positions.out')
open(5,  file='OUTPUT/vels.out')
open(6,  file='OUTPUT/forces.out')


print *,      "Maximum number of threads",    omp_get_max_threads()  
print *,      "Number of threads being used", omp_get_num_threads()
#ifdef TENSOR            
print *,      "Hydrodynamic representation being used is TENSOR"
#else
print *,      "Hydrodynamic representation being used is BEAD"
#endif
print *,      "Time(micro sec.), N_Fiber, N_Segment, T_Length(mm), Avg_Length(mm)"  !2018/09/22

write(301,*), "Maximum number of threads",    omp_get_max_threads() 
write(301,*), "Number of threads being used", omp_get_num_threads()
#ifdef TENSOR            
write(301,*), "Hydrodynamic representation being used is TENSOR"
#else
write(301,*), "Hydrodynamic representation being used is BEAD"
#endif
write(301,*), "Time(micro sec.), N_Fiber, N_Segment, T_Length(mm), Avg_Length(mm)"  !2018/09/22
write(300,*), "Time(micro sec.), N_Fiber, N_Segment, T_Length(mm), Avg_Length(mm)"  !2018/09/22

write(306,*), "Time(micro sec.), a11"                                     !2018/10/07
write(307,*), "Time(micro sec.), Ln(mm)"                                    !2018/10/07

call  read_data( fibers, hinges, simParameters )                              !2018/10/05 add by Hakan

if( simParameters%flow_case==1848 ) then                                      !2018/10/10 HAKAN SHEARRATE OVER T
    simParameters%h= 1                                                        !2018/10/10 add
    call read_data_1848(       flowcase_1848, simParameters )                 !2018/10/10 add
    call set_velocity_1848(    flowcase_1848, simParameters )                 !2018/10/10 add 
    call output_DynamicP_1848( flowcase_1848, simParameters )                 !2018/10/11 增加        
end if


start = OMP_get_wtime()
                                                                              !2018/07/14 增加
i= 0 
t= 0.d0 

allocate( AA(3,3) )

call fiber_regroup_minmax_hinges(   fibers, hinges )                          !2018/08/05 add
call fiber_regroup_minmax_segments( fibers, hinges )                          !2018/08/05 add
call fiber_regroup_ShiftCenterToOrigion( fibers, hinges, simParameters )      !2018/08/05 add
call update_periodic_Initial( fibers, hinges, simParameters )                 !2018/09/22 add
call fiber_regroup_minmax_hinges(   fibers, hinges )                          !2018/09/22 add
call fiber_regroup_minmax_segments( fibers, hinges )                          !2018/09/22 add
call GhostSegments_Dimension( fibers, hinges, ghost_segments, simParameters )  !2018/10/02 修正

!call output_FiberLengthModification( fibers, hinges )                        !2018/09/22 修正
!call output_Initial_Positions_New( fibers, hinges, box_size )                !2018/09/22 增加
!call output_PositionsForTheMomemt ( fibers, hinges, nbr_hinges)              !2018/09/22 因為剛開始多一樣,不用輸出

call output_Length(             t, fibers, hinges )                           !2018/08/12 修正
call output_LengthDistribution( t, fibers, indexA )                           !2018/08/12 增加
call output_OrientationTensor(  t, fibers, hinges, AA )                       !2018/08/12 增加
call simulation_parameter( hinges, simParameters )                            !2018/10/08    修正字串和移動位置

if( simParameters%recover_simulation.eqv..true. ) then 
	n= simParameters%frame*simParameters%writ_period + 1
else 
	n=1
	simParameters%frame= 1
end if

nbr_Fibers_INC= 0                                                             !2018/07/14 修正和移動位置
nbr_Fibers_NEW= ubound(fibers,1)                                              !2018/07/14 修正和移動位置
nbr_Fibers_OLD= nbr_Fibers_NEW                                                !2018/07/14 修正和移動位置

do i=n, simParameters%nbr_intgr
 
    t = simParameters%dt*i                                                    !2018/07/14 修正

    simParameters%time = t                                                    !2018/10/09  add

    if( simParameters%flow_case==1848 ) then 
        
        call set_velocity_1848( flowcase_1848, simParameters )                !2018/10/10 add  Hakan - DynamicParameters OVER T
    
    end if

    
    nbr_Fibers_OLD= ubound(fibers,1)                                          !2018/08/11

    isOutputMessage= .false.                                                  !2018/08/11  add

    call cpu_time (start2)

 	if (MODULO(i,simParameters%break_period)==0 .or. (i.eq.n) ) then 
        
       call bending_torque_whole( fibers, hinges, simParameters )             !2018/10/08  修正
       
       if ( simParameters%allow_breakage ) then
           call hinges_damage(fibers, hinges, simParameters )                 !2018/10/10  修正
       end if

       call GhostSegments_Location( fibers, hinges, ghost_segments, simParameters )               !2018/10/09  修正
       
       call fiber_regroup( fibers, hinges, ghost_segments, neighbor_list, cells, simParameters )  !2018/10/09 修正

       call find_neighbors_new_original( fibers,&
                                         hinges,&
                                         ghost_segments,&
                                         neighbor_list,&
                                         distance_neighbors,&
                                         cells,&
                                         simParameters )                      !2018/10/10  時間花太長,更換為原來的程式
         
!      call fiber_regroup_minmax_hinges(   fibers, hinges )                   !2018/08/05 add
!      call fiber_regroup_minmax_segments( fibers, hinges )                   !2018/08/05 add
!      call output_OrientationTensor( t, fibers, hinges, AA )                 !2018/08/12 增加
!      call output_PositionsForTheMomemt ( fibers, hinges, nbr_hinges)        !2018/09/01 因為沒有用到斷裂
       
    end if

    nbr_Fibers_NEW= ubound(fibers,1)                                          !2018/08/12 增加

    if ( nbr_Fibers_NEW .GT. (nbr_Fibers_OLD + nbr_Fibers_INC) ) then         !2018/08/12 增加
         call output_Length( t, fibers, hinges )                              !2018/08/12 修正       
!        call output_LengthDistribution( t, fibers, indexA )                  !2018/08/12 增加
         isOutputMessage= .true.                                              !2018/08/12 增加
    end if    

#ifdef TENSOR    
    call fiber_calc_tensor( fibers, hinges, simparameters )                   !2018/10/08 change
#else    
    call fiber_calc(        fibers, hinges, simparameters )                   !2018/10/08 change
#endif

    call GhostSegments_NewLocation( hinges, ghost_segments, simParameters )   !2018/10/05  修正  

 	call excl_VolForceMomentsTotal( fibers, hinges, ghost_segments, neighbor_list, simParameters )  !2018/10/09  修正
                                       
    if( .NOT. simParameters%IsPeriodicY ) then
        
        !call fiber_regroup_minmax_hinges( fibers, hinges )                   !2018/08/05 add
        !call fiber_regroup_minmax_segments(fibers, hinges )                  !2018/08/05 add
        !call output_OrientationTensor( t, fibers, hinges, AA )               !2018/08/12 增加              
         
         call excl_VolForceMomentsWalls2( fibers, hinges, simParameters )     !2018/10/09  change
    end if
    
 	call bending_torque_whole( fibers, hinges, simParameters )                !2018/10/08  修正

    call motion_fiber( fibers, hinges, simParameters )                        !2018/10/10 change
    
 	call update_periodic( fibers, hinges, simParameters )  !2018/10/09 change
    
! 	if (MODULO(i,simParameters%writ_period)==0 .or. (i .le. 5) ) then
 	if (MODULO(i,simParameters%writ_period)==0 ) then

 		call output_data( fibers, hinges, simParameters )
        call output_LengthDistribution(  t, fibers, indexA )                  !2018/08/12 增加
        call output_OrientationTensor(   t, fibers, hinges, AA )              !2018/08/12 增加 
        call output_PositionsForTheMomemt ( fibers, hinges, nbr_hinges)       !2018/09/01 跟writ_period一起輸出,可以在Fibers.in給定
        
        call fiber_regroup_minmax_segments( fibers, hinges )                  !2018/08/05 add        
        
        if( isOutputMessage .eq. .false. ) then
            call output_Length( t, fibers, hinges )                           !2018/08/12 修正
        end if
        
 		simParameters%frame = simParameters%frame + 1
        
    end if
   
end do 

call cpu_time(finish)

print *,      "Time Elpased",  OMP_get_wtime()-start, "s"
write(301,*), "Time Elpased",  OMP_get_wtime()-start, "s"

close(300)
close(301)
close(302)                     !2018/08/12
close(303)                     !2018/08/12
close (3) 
close (5)
close (6)
write( *, * ) 'Press Enter to continue' 
read( *, * ) 

end program cube_periodic

! testing for github
!======================================================================