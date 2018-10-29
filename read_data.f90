module m_ReadData   !2018/07/21  change name
    
use m_DataStructures
implicit none
contains
        
subroutine read_data(frame,&
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
                     simParameters)

real(8)                                :: E_Young, min_curv, r_fiber, viscosity, ex_vol_const,&
                                          gamma_dot, epsilon_dot, dt, void, fric_coeff, distanceFactor
integer(8)                             :: nbr_neighbors, flow_case, nbr_intgr, writ_period, break_period,&
                                          i, j, k, n, nbr_fibers, nbr_hinges, nbr_hinges_total, nbr_Segments_total
logical, intent(out)                   :: periodic_boundary
real(8), dimension(3)                  :: box_size,boxSize,boxOrigin
type(rod)  , allocatable, dimension(:) :: hinges
type(fiber), allocatable, dimension(:) :: fibers
logical                                :: recover_simulation,allow_breakage
integer(8)                             :: frame
real(8), dimension(3)                  :: coord
logical                                :: is_fric_wall,printVelocities,isPeriodicX,isPeriodicY,isPeriodicZ
type(simulationParameters)             :: simParameters

namelist /input/ recover_simulation,&
                 fric_coeff,&
                 is_fric_wall,&
                 E_Young,&
                 allow_breakage,&
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
                 break_period,&
                 printVelocities,&
                !boxSize,&
                !boxOrigin,&
                !isPeriodicX,&
                !isPeriodicY,&
                !isPeriodicZ,&
                 distanceFactor
 isPeriodicY = .false.                        
	open(1,file='INPUT/Fibers.in', status='old')

    read(1,nml = input)
	close(1)
    print *,"------------------------------------------"
	print *,"INPUT SUMMARY"
	print *,"recover_simulation", recover_simulation
	print *,"Fric Coeff", fric_coeff
	print *,"is_fric_wall", is_fric_wall
    print *,"E_Young", E_Young
    print *,"min_curv", min_curv
    print *,"r_fiber", r_fiber
    print *,"viscosity", viscosity
    print *,"ex_vol_const", ex_vol_const
    print *,"nbr_neighbors", nbr_neighbors
    print *,"gamma_dot", gamma_dot
    print *,"epsilon_dot", epsilon_dot
    print *,"flow_case", flow_case
    print *,"periodic_boundary", periodic_boundary
    print *,"box_size", box_size(1)
    print *,"box_size", box_size(2)
    print *,"box_size", box_size(3)    
    print *,"dt      ", dt 
    print *,"nbr_intgr       ", nbr_intgr
    print *,"writ_period     ", writ_period
    print *,"break_period    ", break_period
    print *,"allow breakage  ",allow_breakage
    print *,"printVelocities ", printVelocities
    print *,"distance Factor ", distanceFactor
	print *,"------------------------------------------"
   	close(1)
 simParameters%IsPeriodicY = isPeriodicY
    open(3, file='OUTPUT/positions.out')
	if (recover_simulation.eqv..true.) then
		open(4,file='OUTPUT/nbr_frames.txt')
		read(4,*), frame
                frame=frame-1
                close (4)
		n=frame*writ_period
		do i=1, frame-1
			read (3,*) nbr_fibers
			do j=1, nbr_fibers
				read (3,*) nbr_hinges
				do k=1,nbr_hinges
					read (3,*) coord(1), coord(2), coord(3)
				end do
		    end do
		end do

		open (44, file='OUTPUT/temp.txt')
		read (3,*) nbr_fibers
		write (44,*) nbr_fibers
		do j=1, nbr_fibers
			read (3,*) nbr_hinges
			write (44,*) nbr_hinges
			do k=1,nbr_hinges
				read (3,*), coord(1), coord(2), coord(3)
				write (44,*), "0", coord(1), coord(2), coord(3)
			end do
		end do

                close(44)
	else
		frame=1	
	end if
    !print *,"test 1"
		
        
	!Counting the number of segments
    if (recover_simulation.eqv..false.) then
		open (2, file="INPUT/Initial_Positions.txt")
	else
		open (2, file="OUTPUT/temp.txt")
	end if

	read (2,*), nbr_fibers

	allocate (fibers(nbr_fibers))

	nbr_hinges_total=0

	do i=1, nbr_fibers
		read(2,*), fibers(i)%nbr_hinges
		!print *,"hinge_has", fibers(i)%nbr_hinges
		do j=1, fibers(i)%nbr_hinges
			read(2,*), void, void, void
		        nbr_hinges_total=nbr_hinges_total+1
		end do
	end do
	close (2)

    nbr_Segments_total= nbr_hinges_total-nbr_fibers
    
    print *,"Total fibers   are ", nbr_fibers
    print *,"Total Segments are ", nbr_Segments_total
    
	allocate (hinges(nbr_hinges_total))

	
	if (recover_simulation.eqv..false.) then
		open (99, file="INPUT/Initial_Positions.txt")
	else
		open (99, file="OUTPUT/temp.txt")
	end if

	k=1
	read (99,*), nbr_fibers
	do i=1, nbr_fibers
		read (99,*), fibers(i)%nbr_hinges
		if (i==1) then
			fibers(i)%first_hinge=1
	        else
			fibers(i)%first_hinge=fibers(i-1)%first_hinge+fibers(i-1)%nbr_hinges
		end if
		do j=1, fibers(i)%nbr_hinges
			hinges(k)%in_fiber=i
			read(99,*), hinges(k)%is_stationary, hinges(k)%X_i(1), hinges(k)%X_i(2), hinges(k)%X_i(3)!,&
             hinges(k)%v_i =0
             hinges(k)%omega =0
			            !hinges(k)%fric
			k=k+1
		end do
	end do
	close (99)
    
      
	!do i=1, nbr_fibers
	!	print*,"Fiber Info:", fibers(i)%first_hinge,fibers(i)%nbr_hinges
	!end do

	!do i=1, ubound(hinges,1)
	!	print *, hinges(i)%X_i(1), hinges(i)%X_i(2), hinges(i)%X_i(3)
	!end do
	!print *,"FRAME", frame
    print *, "periodic_boundary    ", periodic_boundary
    
end subroutine read_data
                     
end module m_ReadData
