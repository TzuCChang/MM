!*===================================================================
!*  Units in SI ( m, N, Pa, etc )
!*
program cube_periodic

use m_DataStructures

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

call simulation_parameter( hinges, r_fiber, gamma_dot, epsilon_dot, flow_case, simParameters )

call output_Length( t, fibers, hinges, frame, printVelocities )

do i=n,  nbr_intgr
 
    t = dt*i

 	if (MODULO(i,break_period)==0 .or. (i.eq.n) ) then 

       call bending_torque_whole(fibers, hinges, E_Young, Inertia_Moment)

       if (allow_breakage)then
           call hinges_damage(fibers, hinges, min_curv, r_fiber)
       end if
       
 	   call fiber_regroup( fibers,&    
 	                       hinges,&
 	                       ghost_segments,&
                           E_Young,&        
                           Inertia_Moment,&      
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
                           Nbr_bins,&        
                           distanceFactor,&
                           simParameters )
                                                      
         call find_neighbors_new( fibers,& 
                                  hinges,&
                                  ghost_segments,&
                                  nbr_neighbors,&
                                  neighbor_list,&
                                  distance_neighbors,&
                                  r_fiber,&
                                  cells,&
                                  Nbr_bins,&
                                  distanceFactor )
        
 	end if

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

 	call excl_vol_forces_moments_total( fibers,&
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

    if(.NOT. simParameters%IsPeriodicY ) then
    call excl_vol_forces_moments_walls2( fibers,&   
                                         hinges,&
                                         r_fiber,&
                                         ex_vol_const,&
                                         box_size,&
                                         fric_coeff,&
                                         is_fric_wall,&
                                         gamma_dot )  
    end if
    
 	call bending_torque_whole(fibers, hinges, E_Young, Inertia_Moment)
    
    call motion_fiber(fibers, hinges, r_fiber)
    
 	call update_periodic(fibers, hinges, dt, periodic_boundary, box_size, gamma_dot,dt* (i-1) )
    
 	if (MODULO(i,writ_period)==0 ) then
 		call output_data(t, fibers, hinges, frame,printVelocities)
    else 
        nbr_Fibers_NEW= ubound(fibers,1)                                     
        if ( nbr_Fibers_NEW .GT. (nbr_Fibers_OLD+nbr_Fibers_INC) ) then               
              call output_Length( t, fibers, hinges, frame, printVelocities ) 
              nbr_Fibers_OLD= nbr_Fibers_NEW                                   
        end if
    end if

end do 

end program cube_periodic
!*
!*===================================================================