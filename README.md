# MechanisticModelFortran

This is a particle level simulation of a fiber suspension. 
The fibers are 
simulated as rigid bodies and a flow field is imposed onto them. It is designed to be used with a 
periodic shear cell. 
It requires ifort to compile. 

This branch represents hydrodynamically the rods
as ellipsoids. 

To run simulations the project has to be set up as follows:



- Project
  
	- INPUT
    
		- Fibers.in
    
		- Initial_Positions.txt
  
	- OUTPUT




**Fibers.in** contains information of all the simulation parameters. 

```
&input
recover_simulation =          ! If you want to restart a simulation, use .true. or .false

fric_coeff         =          ! Friction coefficient [-]

is_fric_wall       =          ! .true. if the fibers have friction with wall

E_Young            =          ! fibers youngs modulus. It will be used if fibers have more than one segment [Pa]

allow_breakage     =          ! .true. to allow fibers to break

min_curv           =          ! Radius of curvature at which the fiber will break [m]

r_fiber            =          ! Fiber Radius  [m]

viscosity          =          ! Matrix viscosity [Pa-s] 

ex_vol_const       =          ! Constant used for the penalty force in the collision response [N]

nbr_neighbors      =          ! Maximum allowed neighbors per segment 

gamma_dot          =          ! Shear Rate [1/s] 

epsilon_dot        =          ! Elongational rate [1/s]

flow_case          =          ! 

periodic_boundary  =          ! .true. if periodicity in x and z is wanted

box_size(1)        =          ! Box size in the x directino [m]

box_size(2)        =          ! Box size in the y directino [m]

box_size(3)        =          ! Box size in the z directino [m]

dt                 =          ! time step [s]

nbr_intgr          =          ! total number of integrations

writ_period        =          ! frequency at which results will be printed in timesteps

break_period       =          ! frequency at which neighbor lists and breakage will be checked in timesteps

printVelocities    =          ! .true. if you want to print the velocities at each node

isPeriodicY        =          ! .true. if lees edwards conditions will be used

distanceFactor     =          ! multiplier of the fiber radius for region size where neighbors will be searched for

/

```



**Initial_Positions.txt** contains the positions of the fibers in the following format: 

```
#fibers
# 
	nodes in fiber 1
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node 1 fiber 1
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node 2 fiber 1
  .
  .
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node n fiber 1

	# nodes in fiber 2
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node 1 fiber 2
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node 2 fiber 2
  .
  .
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node n fiber 2
.
.
	
# nodes in fiber m
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node 1 fiber m
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node 2 fiber m
  .
  .
  
	is_node_stationary Xcoord Ycoord Zcoord   !Node n fiber m
```