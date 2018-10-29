!====================================================================
module m_SetVelocity
    
use m_DataStructures
use m_UtilityLib

implicit none
contains

subroutine set_velocity( X,&    !2018/07/21 change name
                         vel,&
                         omega,&
                         gamma_dot,&
                         epsilon_dot,&
                         flow_case )

implicit none
real(8), dimension (3):: X
real(8), dimension (3):: vel
real(8), dimension (3):: omega
real(8)               :: b, gamma_dot, epsilon_dot
integer(8)            :: flow_case

vel=0
omega=0
b=0

if( flow_case==3 ) then 
	b=1
end if

if( flow_case==1 ) then
      !vel(2)=X(1)*gamma_dot
      !omega(3)=+gamma_dot
      vel(1)=X(2)*gamma_dot
      !omega(3)=-gamma_dot
      omega(3)=-0.5*gamma_dot !CORRECTED TS
else
	vel(1)=-0.5*epsilon_dot*(1+b)*X(1)
    vel(2)=-0.5*epsilon_dot*(1-b)*X(2)
    vel(3)=epsilon_dot*X(3)
end if
		
end subroutine set_velocity   !2018/07/21 change name



!2018/07/14 新增加訪問學者
subroutine set_velocity_NEW( X,&   !2018/07/21 change name
                             vel,&
                             omega,&
                             gamma_dot,&
                             epsilon_dot,&
                             flow_case )

implicit none
real(8), dimension (3):: X
real(8), dimension (3):: vel
real(8), dimension (3):: omega
real(8)               :: b, gamma_dot, epsilon_dot
integer(8)            :: flow_case

vel=0
omega=0
b=0

if( flow_case==3 ) then 
	b=1
end if

if( flow_case==1 ) then
    
      !vel(2)=X(1)*gamma_dot
      !omega(3)=+gamma_dot
      vel(1)=X(2)*gamma_dot
      !omega(3)=-gamma_dot
      omega(3)=-0.5*gamma_dot !CORRECTED TS

      !!! 06/07/2017 
      !!! This is written by Hakan Celik     
      !!! My own Velocity for translational validation

else if( flow_case==61 ) then
            
        ! Step Constant Step  0.3
        ! Ramp Constant       0.2
        ! Constant            0.1    
        ! Interval 1
        ! Constant Velocity
        
        if (X(1) <= 20.00E-03 ) then
            vel(1)= 0.1

        ! Interval 2
        ! Ramp Velocity
        
        else if (X(1) > 20.00E-03 .AND. X(1) <= 40.00E-03) then
             vel(1) = ((0.2-0.1)/(40.00E-03 - 20.00E-03))*(X(1)-20.00E-03)+0.1

        ! Interval 3
        ! Constant Velocity
        
        else if (X(1) > 40.00E-03 .AND. X(1) <= 60.00E-03) then
             vel(1) = 0.2
             
        ! Interval 4
        ! Step Constant Step Velocity
        
        else if (X(1) > 60.00E-03 .AND. X(1) <= 80.00E-03) then
             vel(1) = 0.3

        ! Interval 5
        ! Constant Velocity
        
        else if (X(1) > 80.00E-03 .AND. X(1) <= 100.00E-03) then
             vel(1)= 0.2

        ! Interval 6
        ! Ramp Velocity 
        
        else if (X(1) > 100.00E-03 .AND. X(1) <= 120.00E-03) then
             vel(1)= ((0.1-0.2)/(120.00E-03 - 100.00E-03))*(X(1)-100.00E-03)+0.2

        ! Interval 7
        
        else
            vel(1) = 0.1
        end if
        
else if( flow_case==34 ) then 
   
        if (X(1) <= 20.00E-03 ) then
            vel(1)= 0.1
        
        ! Ramp Velocity
        else if (X(1) > 20.00E-03 .AND. X(1) <= 40.00E-03) then
             vel(1) = ((0.2-0.1)/(40.00E-03 - 20.00E-03))*(X(1)-20.00E-03)+0.1

        else
            vel(1) = 0.1
        end if
        vel(1) = 1
        
else if( flow_case==1461 ) then 
        vel(1) = epsilon_dot
   
        if (X(1) <= 0 ) then
            vel(1)= -epsilon_dot
        else
            vel(1)= epsilon_dot
        end if
        
! START - HAKAN Shearrate over TIME

else if( flow_case==1848 ) then 
     
      !print *, gamma_dot
      !vel(2)=X(1)*gamma_dot
      !omega(3)=+gamma_dot
      
      vel(1)=X(2)*gamma_dot
      
      !omega(3)=-gamma_dot
      
      omega(3)=-0.5*gamma_dot !CORRECTED TS       
      
! END - HAKAN Shearrate over TIME
! This is written by Hakan Celik  
! 06/07/2017
 
else
	vel(1)=-0.5*epsilon_dot*(1+b)*X(1)
    vel(2)=-0.5*epsilon_dot*(1-b)*X(2)
    vel(3)=epsilon_dot*X(3)
end if
		
end subroutine set_velocity_NEW   !2018/07/21 change name

end module m_SetVelocity             !2018/07/21 change name
