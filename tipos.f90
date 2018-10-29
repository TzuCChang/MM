module m_DataStructures

type pts
integer(8)              :: n_pto_rod
integer(8)              :: n_fib
integer(8)              :: n_rod_abs
integer(8)              :: conn1
integer(8)              :: conn2
real   (8)              :: equl_angle
real(8), dimension(3)   :: coord
logical                 :: eval_hydro_force
logical                 :: eval_hydro_torque
real(8), dimension(3)   :: fluid_velocity
real(8), dimension(3)   :: fluid_omega
real(8)                 :: diameter
end type pts

type segment
real(8), dimension(3)   :: A
real(8), dimension(3)   :: B
integer, dimension(3)   :: ind
integer, dimension(2)   :: orig_pos
integer                 :: axis_loc
integer                 :: ix
integer                 :: iy
integer                 :: iz
end type segment

type fiber
integer(8)              :: nbr_hinges
integer(8)              :: first_hinge
real(8)                 :: Length      !2019/09/09  add
end type fiber

type rod
integer(8)                :: in_fiber
integer(8)                :: nbr_beads
real(8)                   :: length
real(8)                   :: length2
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
real(8)			          :: ave_viscosity
real(8)			          :: gamma_dot
real(8)   , dimension(3)  :: r
real(8)   , dimension(3)  :: r_unit
real(8)   , dimension(3)  :: r_sum
real(8)   , dimension(3,3):: r_prod_sum
real(8)   , dimension(3,3):: r_times_u_sum
real(8)   , dimension(3,3):: A
real(8)   , dimension(3,3):: C
real(8)   , dimension(3)  :: H
real(8)                   :: alpha
real(8)                   :: curv
logical                   :: is_broken
logical                   :: is_separated
logical                   :: is_segment
integer                   :: is_stationary
real(8)   , dimension(3)  :: omega

real(8)   , dimension(3)  :: omega_old   !2019/09/05

real(8)                   :: fric
integer   , dimension(3)  :: indx
integer                   :: ind
end type rod

type vec_arr
real(8), dimension(3)::pb
real(8)              ::Gab_norm
logical              ::coll_course
real(8), dimension(3):: Gab
real(8), dimension(3):: r
real(8)              :: Sba
real(8), dimension(3):: vec
end type vec_arr

type simple_vec_arr
real(8), dimension(3)::r
end type simple_vec_arr

type simulationParameters
real(8)                  ::X_A, Y_A, X_C, Y_C, Y_H
real(8), dimension(3,3)  ::E_oo
real(8), dimension(3,3,3)::eps
logical                  :: IsPeriodicY

end type simulationParameters


type cell
integer, dimension(3)::indx
integer, dimension(2)::ghost_limits
end type cell

type bound_box
real(8), dimension(3):: X
real(8), dimension(3):: W
end type bound_box

end module m_DataStructures

