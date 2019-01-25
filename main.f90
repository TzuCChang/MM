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
!======================================================================

implicit none
type(rod),      dimension(:), allocatable  :: hinges
type(fiber),    dimension(:), allocatable  :: fibers
type(segment),  dimension(:), allocatable  :: ghost_segments
type(cell),     dimension(:), allocatable  :: cells
type(DynamicP), dimension(:), allocatable  :: flowcase_1848
type(simulationParameters)                 :: simParameters
logical                                    :: isOutputMessage, isProgramStop, icheck_Breakage

integer(8), dimension(:,:), allocatable    :: neighbor_list
integer(8), dimension(:),   allocatable    :: indexA
integer(8)                                 :: i, j, k, n, nbr_Fibers_OLD, nbr_Fibers_NEW, nbr_Fibers_INC

!======================================================================

simParameters%start = OMP_get_wtime()                                         !2018/12/09

call output_OpenFiles( simParameters )

call read_data( fibers, hinges, simParameters )                               !2018/11/25 
                                                       
!======================================================================

if( simParameters%flow_case==1848 ) then                                      !2018/10/10 HAKAN SHEARRATE OVER T
    
    simParameters%h= 1                                                        !2018/10/10 add
    call read_data_1848(       flowcase_1848, simParameters )                 !2018/10/10 add
    call set_velocity_1848(    flowcase_1848, simParameters )                 !2018/10/10 add 
    
         call output_DynamicP_1848( flowcase_1848, simParameters )            !2018/10/11 add
    
end if

         call output_minmax_hinges(   fibers, hinges )                         !2018/10/27 change name
         call output_minmax_segments( fibers, hinges )                         !2018/10/27 change name

!======================================================================

call fiber_regroup_ShiftCenterToOrigion( fibers, hinges, simParameters )       !2018/08/05 add
call update_periodic_Initial( fibers, hinges, simParameters )                  !2018/09/22 add

         call output_minmax_hinges(    fibers, hinges )                        !2018/10/27 change name
         call output_minmax_segments(  fibers, hinges )                        !2018/10/27 change name

call GhostSegments_Dimension( fibers, hinges, ghost_segments, simParameters )  !2018/10/02 修正

         call output_LengthDistribution( fibers, indexA, simParameters )       !2018/10/29
         call output_OrientationTensor(  fibers, hinges, simParameters )       !2018/10/12 增加
         call output_Length(             fibers, hinges, simParameters )       !2018/12/09

call simulation_parameter( hinges, simParameters )                             !2018/10/08    修正字串和移動位置
!======================================================================

nbr_Fibers_INC= 0                                                             !2018/07/14 修正和移動位置
nbr_Fibers_NEW= ubound(fibers,1)                                              !2018/07/14 修正和移動位置
nbr_Fibers_OLD= nbr_Fibers_NEW                                                !2018/07/14 修正和移動位置

simParameters%nStep_max= 128                                                  !2018/10/29
!simParameters%displ_max= 1.0d-4*(3.14159d0/180.d0)/10.d0                     !2018/12/09 segments length= 0.10 mm為半徑,走一個步的距離,圓周長的3,600分之一
simParameters%nStep=     0                                                    !2018/10/29 add


isProgramStop  = .false.                                                      !2018/10/29
icheck_Breakage= .true.                                                       !2018/10/29

n= simParameters%n                                                            !2018/11/25 修正
i= simParameters%n                                                            !2018/11/25 修正

do while ( isProgramStop .eqv. .false. )

    simParameters%nStep_Total= simParameters%nStep_Total + 1                  !2018/11/18

    if( simParameters%nStep .eq. simParameters%nStep_max ) then               !2018/10/29
        
        i= i + 1                                                              !2018/10/29 
        simParameters%nStep= 0                                                !2018/10/29
        icheck_Breakage= .true.                                               !2018/10/29

    end if
    
    if( i .eq. simParameters%nbr_intgr ) then                                 !2018/10/29
        
        isProgramStop= .true.                                                 !2018/10/29
        
    end if

    
    if( simParameters%flow_case==1848 ) then 
        
        call set_velocity_1848( flowcase_1848, simParameters )                !2018/10/10 add  Hakan - DynamicParameters OVER T
        
    end if

    nbr_Fibers_OLD= ubound(fibers,1)                                          !2018/08/11

    isOutputMessage= .false.                                                  !2018/08/11  add

    if( icheck_Breakage .eqv. .true. ) then
 	if (MODULO(i,simParameters%break_period)==0 .or. (i.eq.n) ) then 

        icheck_Breakage= .false.                                              !2018/10/29
       
        call hinges_curv_distribution( fibers, hinges, simParameters )        !2019/11/23  add
       
        if ( simParameters%allow_breakage ) then
           call hinges_damage( fibers, hinges, simParameters )                !2018/10/10  修正
        end if

        call GhostSegments_Location( fibers, hinges, ghost_segments, simParameters )               !2018/10/09  修正
       
        call fiber_regroup( fibers, hinges, ghost_segments, neighbor_list, cells, simParameters )  !2018/10/09 修正

        call find_neighbors( fibers, hinges, ghost_segments, neighbor_list, cells, simParameters ) !2018/11/25  change

    end if
    end if

    nbr_Fibers_NEW= ubound(fibers,1)                                          !2018/08/12 增加

    if ( nbr_Fibers_NEW .GT. (nbr_Fibers_OLD + nbr_Fibers_INC) ) then         !2018/08/12 增加
        
             call output_Length( fibers, hinges, simParameters )              !2018/10/29      
             call output_LengthDistribution( fibers, indexA, simParameters )  !2018/10/29
             call output_OrientationTensor(  fibers, hinges, simParameters )  !2018/10/12 增加             

         isOutputMessage= .true.                                              !2018/08/12 增加
         
    end if    

if( simParameters%TensorOrBeads .eq. 1 ) then                                 !2018/12/02
    
    call fiber_calc_tensor( fibers, hinges, simparameters )                   !2018/10/08 change
else    
    call fiber_calc(        fibers, hinges, simparameters )                   !2018/10/08 change
endif

    call GhostSegments_NewLocation( hinges, ghost_segments, simParameters )   !2018/10/05  修正  

 	call excl_VolForceMomentsTotal( fibers, hinges, ghost_segments, neighbor_list, simParameters )  !2018/10/09  修正
                                       
    if( .NOT. simParameters%IsPeriodicY ) then        
         
        call excl_VolForceMomentsWalls2( fibers, hinges, simParameters )      !2018/10/09  change
        
    end if
    
 	call bending_torque_whole( fibers, hinges, simParameters )                !2018/10/08  修正

    call motion_fiber( fibers, hinges, simParameters )                        !2018/10/10 change
    
 	call update_periodic( fibers, hinges, simParameters )                     !2018/10/09 change
    
 	if (MODULO(i+1,simParameters%writ_period)==0 .and. &
       (simParameters%nStep .eq. simParameters%nStep_max) ) then              !2018/11/18  add
        
 		simParameters%frame = simParameters%frame + 1                         !2018/11/18  add
        
 		call output_data( fibers, hinges, simParameters )
        call output_LengthDistribution(  fibers, indexA, simParameters )      !2018/10/29
        call output_OrientationTensor(   fibers, hinges, simParameters )      !2018/10/12 增加 
        call output_PositionsForTheMomemt ( fibers, hinges )                  !2018/10/12 跟writ_period一起輸出,可以在Fibers.in給定
        
        call output_minmax_segments( fibers, hinges )                         !2018/10/27 change name        
        
        if( isOutputMessage .eq. .false. ) then
            call output_Length( fibers, hinges, simParameters )               !2018/10/29
        end if

    end if
   
end do 

call output_CloseFiles( simParameters )

write(*,*) 'Press Enter to continue' 
read(*,*)

end program cube_periodic

! End of Program  2018/11/24 
!======================================================================