module m_DataStructures

type pts
logical                 :: eval_hydro_force
logical                 :: eval_hydro_torque
integer(8)              :: n_pto_rod
integer(8)              :: n_fib
integer(8)              :: n_rod_abs
integer(8)              :: conn1
integer(8)              :: conn2
real(8), dimension(3)   :: coord
real(8), dimension(3)   :: fluid_velocity
real(8), dimension(3)   :: fluid_omega
real(8)                 :: diameter
real(8)                 :: equl_angle
end type pts

type segment
integer(8), dimension(3)   :: ind
integer(8), dimension(2)   :: orig_pos
integer(8)                 :: axis_loc
integer(8)                 :: ix
integer(8)                 :: iy
integer(8)                 :: iz
real(8),    dimension(3)   :: A
real(8),    dimension(3)   :: B
end type segment

type fiber
integer(8)              :: nbr_hinges
integer(8)              :: first_hinge
real(8)                 :: Length      !2019/09/09  add
end type fiber

type rod
logical                   :: is_broken
logical                   :: is_separated
logical                   :: is_segment
integer(8), dimension(3)  :: indx
integer(8)                :: is_stationary
integer(8)                :: ind
integer(8)                :: in_fiber
integer(8)                :: nbr_beads
real(8)   , dimension(3)  :: X_i
real(8)   , dimension(3)  :: T
real(8)   , dimension(3)  :: F_excl_vol
real(8)   , dimension(3)  :: T_excl_vol
real(8)   , dimension(3)  :: v_i
real(8)   , dimension(3)  :: v_i_old   !2019/09/05
real(8)   , dimension(3)  :: omega_i
real(8)   , dimension(3)  :: omega_fluid_sum
real(8)   , dimension(3)  :: u_fluid_sum
real(8)   , dimension(3)  :: u_oo
real(8)   , dimension(3)  :: omega_oo
real(8)   , dimension(3)  :: r
real(8)   , dimension(3)  :: r_unit
real(8)   , dimension(3)  :: r_sum
real(8)   , dimension(3,3):: r_prod_sum
real(8)   , dimension(3,3):: r_times_u_sum
real(8)   , dimension(3,3):: A
real(8)   , dimension(3,3):: C
real(8)   , dimension(3)  :: H
real(8),    dimension(3)  :: omega
real(8),    dimension(3)  :: omega_old   !2019/09/05
real(8)                   :: length
real(8)                   :: length2
real(8)                   :: alpha
real(8)                   :: curv
real(8)                   :: fric
real(8)			          :: ave_viscosity
real(8)			          :: gamma_dot
end type rod

type vec_arr
logical              :: coll_course
real(8), dimension(3):: pb
real(8), dimension(3):: Gab
real(8), dimension(3):: r
real(8), dimension(3):: vec
real(8)              :: Sba
real(8)              :: Gab_norm
end type vec_arr

type simple_vec_arr
real(8), dimension(3)::r
end type simple_vec_arr

type simulationParameters

logical                  :: IsPeriodicY, periodic_boundary
logical                  :: recover_simulation, allow_breakage
logical                  :: is_fric_wall, printVelocities

integer(8), dimension(3) :: Nbr_bins, box_dimension
integer(8)               :: flow_case, nbr_Dynamic, frame, h, nStep_max, nStep_Total, k, n     !2018/11/25
integer(8)               :: nbr_neighbors, nbr_intgr, writ_period, break_period, nStep         !2018/10/27

real(8), dimension(3,3,3):: eps
real(8), dimension(3,3)  :: E_oo
real(8), dimension(3,3)  :: AA                                                !2018/10/12
real(8), dimension(3)    :: box_size
real(8)    :: X_A, Y_A, X_C, Y_C, Y_H, pi                                     !2018/11/25
real(8)    :: FiberLength, FiberVolume, BoxVolume, VolumeFraction
real(8)    :: E_Young, min_curv, r_fiber, viscosity, ex_vol_const
real(8)    :: gamma_dot, epsilon_dot, fric_coeff, distanceFactor
real(8)    :: dt, time, Inertia_Moment, Controltime, displ_max                !2018/10/27
real(8)    :: start, finish                                                   !2018/11/25

end type simulationParameters


type cell
integer(8), dimension(3)::indx
integer(8), dimension(2)::ghost_limits
end type cell

type bound_box
real(8), dimension(3):: X
real(8), dimension(3):: W
end type bound_box

type DynamicP                              !2018/10/10  new add
real(8)     :: Duration
real(8)     :: Shearrate
real(8)     :: Viscosity
end type DynamicP


end module m_DataStructures

