module m_Motion
    
    use m_DataStructures
    use MATRIX_SOLVERS
    use omp_lib
    implicit none
    contains

subroutine motion_fiber(fibers, hinges, r_bead)
    type(fiber) , dimension(:):: fibers
    type(rod), dimension(:)   :: hinges
    integer(8)                :: i
    real(8)                   :: r_bead, viscosity,timex, timeA

!   timeA = OMP_get_wtime()
!   call omp_set_num_threads(6)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE (i) !SCHEDULE(DYNAMIC)

    do i=1, ubound(fibers,1)   
        
        call motion_matrix( hinges(fibers(i)%first_hinge:&
                            fibers(i)%first_hinge+fibers(i)%nbr_hinges-1),&
                            r_bead )
    end do
    
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!   print *,"V_i"
!   do i=1, ubound(hinges,1)
!	print *, hinges(i)%V_i
!	
!   end do
!   print *, " time elpased mot " , OMP_get_wtime()- timeA

end subroutine motion_fiber


subroutine copyToBanded(AB, KL, KU, nbRows, i1, i2, j1, j2, A)
    ! Copy a range of a regular matrix to a banded matrix
    real(8), dimension(:,:), allocatable :: AB
    real(8)                              ::  A(:,:)
    integer(8)                           :: i1, i2, j1, j2, KL, KU, nbRows, ii, jj, j, i
    
    jj =1
    do j=j1, j2
        ii =1
        do i= i1 , i2
            if( i >= max(1,j-KU) .and. i <= min(nbRows,j+KL)) then
                AB(KL+KU+1+i-j,j) = A(ii,jj) !for max(1,j-KU)<=i<=min(N,j+KL)
            end if
            ii = ii+1
        end do
        jj = jj+1
    end do

end subroutine copyToBanded

!*****************************************************************************************
subroutine motion_matrix(fiber_hinges, r_bead)
    implicit none
    type (rod), dimension(:)             :: fiber_hinges
    real(8), dimension(:,:), allocatable :: Amat, bvec, xvec, AB
    real(8), dimension(9,15)             :: mat
    real(8), dimension(9,1)              :: vec
    real(8)                              :: r_bead, drag_coeff_vel, drag_coeff_omega, finish2, start2, timex
    real(8), parameter                   :: pi=3.141592
    integer(8)                           :: i, k, j, l, kl, ku, minSegs 
    integer                              :: n, ii, jj

    drag_coeff_vel = 6d0*pi*r_bead
    drag_coeff_omega= 8d0*pi*r_bead**3d0
    minSegs =4

    n=9*(ubound(fiber_hinges,1)-1)   !9 eqs per every rod
    if ((ubound(fiber_hinges,1)-1) .LT. minSegs) then
        !print *, " Regular Solver "
        allocate(Amat(n,n))
        Amat = 0
    else
        !print *, " Banded Solver "
        KL = 11
        KU = 9
        allocate(AB(2*KL+KU+1,n))
        AB = 0
    end if

    !allocate(xvec(n,1))
    !Amat=0
    !xvec=0
    
    allocate(bvec(n,1))   
    bvec=0

        do i=1, ubound (fiber_hinges,1)-1                      !2018/09/05  修正
           fiber_hinges(i)%v_i_OLD= fiber_hinges(i)%v_i        !2018/09/05  修正
        end do
        
        do i=1, ubound (fiber_hinges,1)-2                      !2018/09/05  修正
           fiber_hinges(i)%omega_OLD= fiber_hinges(i)%omega    !2018/09/05  修正
        end do    

    !print *, "DIMENSION OF FIBER HINGES", ubound (fiber_hinges,1)
    if ((ubound(fiber_hinges,1)-1) .LT. minSegs) then
        if (ubound (fiber_hinges,1)==2) then
#ifdef TENSOR
            !call mini_mat(mat, vec, fiber_hinges(1), fiber_hinges(2), drag_coeff_vel, drag_coeff_omega)
            call mini_mat_tensor(mat, vec, fiber_hinges(1), drag_coeff_vel, drag_coeff_omega)
#else
            call mini_mat(mat, vec, fiber_hinges(1), fiber_hinges(2), drag_coeff_vel, drag_coeff_omega)
#endif
            !mini_mat_tensor
            Amat(1:9,1:6) =mat(1:9, 4:9 )
            Amat(1:9,7:9) =mat(1:9,13:15)

            bvec(:,1)=vec(:,1)
        else
            do i=1, ubound (fiber_hinges,1)-1
                !print*, "I", i
                !call mini_mat(mat, vec, fiber_hinges(i), fiber_hinges(i+1), drag_coeff_vel, drag_coeff_omega)
#ifdef TENSOR                
                call mini_mat_tensor(mat, vec, fiber_hinges(i), drag_coeff_vel, drag_coeff_omega)
#else
                call mini_mat(mat, vec, fiber_hinges(i), fiber_hinges(i+1), drag_coeff_vel, drag_coeff_omega)
#endif               
                !mat=i
                !vec=10*i
                if (i==1) then
                    Amat(1:9,1:12)=mat(1:9,4:15)
                else if (i==ubound(fiber_hinges,1)-1) then
                    Amat(9*(i-1)+1: 9*i, 9*(i-1)-2: 9*(i-1)+6)= mat(1:9,   1:9)
                    Amat(9*(i-1)+1: 9*i, 9*(i-1)+7: 9*(i-1)+9)= mat(1:9, 13:15)
                else
                    Amat(9*(i-1)+1: 9*i, 9*(i-1)-2: 9*(i-1)+12)= mat
                end if
                bvec(9*(i-1)+1: 9*i,1)=vec(1:9,1)

            end do
        end if

    else ! Use banded matrix! 
        if (ubound (fiber_hinges,1)==2) then

#ifdef TENSOR
            !call mini_mat(mat, vec, fiber_hinges(1), fiber_hinges(2), drag_coeff_vel, drag_coeff_omega)
            call mini_mat_tensor(mat, vec, fiber_hinges(1), drag_coeff_vel, drag_coeff_omega)
#else
            call mini_mat(mat, vec, fiber_hinges(1), fiber_hinges(2), drag_coeff_vel, drag_coeff_omega)
#endif
            call copyToBanded(AB, KL, KU, N, 1, 9, 1, 6, mat(1:9, 4:9 ))
            call copyToBanded(AB, KL, KU, N, 1, 9, 7, 9, mat(1:9,13:15))
            !Amat(1:9,1:6) =mat(1:9, 4:9 )
            !Amat(1:9,7:9) =mat(1:9,13:15)

            bvec(:,1)=vec(:,1)
        else

            do i=1, ubound (fiber_hinges,1)-1
                !print*, "I", i
#ifdef TENSOR                
                call mini_mat_tensor(mat, vec, fiber_hinges(i), drag_coeff_vel, drag_coeff_omega)
#else
                call mini_mat(mat, vec, fiber_hinges(i), fiber_hinges(i+1), drag_coeff_vel, drag_coeff_omega)
#endif      
                !mat=i
                !vec=10*i
                if (i==1) then
                    call copyToBanded(AB, KL, KU, N, 1, 9, 1, 12, mat(1:9,4:15))
                    !Amat(1:9,1:12)=mat(1:9,4:15)
                else if (i==ubound(fiber_hinges,1)-1) then
                    call copyToBanded(AB, KL, KU, n, 9*(i-1)+1, 9*i, 9*(i-1)-2, 9*(i-1)+6, mat(1:9,   1:9))
                    !call copyToBanded(AB, KL, KU, n, 9*(i-1)+1, 9*i, 9*(i-1)+7, 9*(i-1)+9, mat(1:9, 13:15))
                    ! small optimization
                    ii = KL+KU+1+(9*(i-1)+1)-(9*(i-1)+7)
                    jj = (9*(i-1)+7)
                    AB(ii,jj) = -1
                    AB(ii,jj+1) = -1
                    AB(ii,jj+2) = -1
                    !Amat(9*(i-1)+1: 9*i, 9*(i-1)-2: 9*(i-1)+6)= mat(1:9,   1:9)
                    !Amat(9*(i-1)+1: 9*i, 9*(i-1)+7: 9*(i-1)+9)= mat(1:9, 13:15)
                else
                    ! Small optimization
                    call copyToBanded(AB, KL, KU, n, 9*(i-1)+1, 9*i, 9*(i-1)-2, 9*(i-1)+12, mat)
                    !call copyToBanded(AB, KL, KU, n, 9*(i-1)+1, 9*i, 9*(i-1)-2, 9*(i-1)+12-3, mat)
                    !ii = KL+KU+1+(9*(i-1)+1)-(9*(i-1)+12-2)
                    !jj = (9*(i-1)+12-2)
                    !AB(ii,jj) = -1
                    !AB(ii+1,jj+1) = -1
                    !AB(ii+2,jj+2) = -1
                    !Amat(9*(i-1)+1: 9*i, 9*(i-1)-2: 9*(i-1)+12)= mat
                end if
                bvec(9*(i-1)+1: 9*i,1)=vec(1:9,1)

            end do
        end if
    end if
    !print *,"MATRIX A   and VECTOR B"

    !do i=1, ubound (Amat,2)
    !	print *, Amat(i,:), bvec(i,:)
    !end do
    ! Uncomment below to print the matrix A in a file. It is easier to visualize.
    !open(7,file='../OUTPUT/matrix.out')
    
    !print *, ' Amat '
    !    do j=1,n
    !        write(7,"(I8,A)",advance="no") j, ' '
    !        enddo
    !        write(7,*)  ' '
    !
    !    do i=1,n
    !        do j=1,n
    !        write(7,"(e8.2,A)",advance="no") Amat(i,j), ' '
    !        enddo
    !        write(7,*)  ' '
    !    enddo


        !print *, ' Amat '
        !do j=1,n
        !    write(*,"(I8,A)",advance="no") j, ' '
        !    enddo
        !    write(*,*)  ' '
        !
        !do i=1,n
        !    do j=1,n
        !    write(*,"(e8.2,A)",advance="no") Amat(i,j), ' '
        !    enddo
        !    write(*,*)  ' '
        !enddo
    
    !print *, "Fila 1 Matriz A"
    !call cpu_time(start2)
    ! Banded matrix only makes sense if the number of segments is 4 or more.
    if ((ubound(fiber_hinges,1)-1) .LT. minSegs) then
        !print *, " Regular Solver "
        !print *, ' bvec1 ', bvec
        call SOLVER_1(n, Amat, bvec)
        !print *, ' bvec2 ', bvec
        !call SOLVER_4(n, 11, 9, Amat, bvec)
    else
        !print *, " Banded Solver "
        call SOLVER_3(n, 11, 9, AB, bvec)
    end if

    !close (7)
    !call cpu_time(finish2)
    !timex=timex+ finish2-start2
    !call SOLVER_1(n, Amat, bvec, xvec)
    !stop
    if (ubound (fiber_hinges,1)==2) then

        do i=1,3
            fiber_hinges(1)%v_i(i)=bvec(i,1)
            fiber_hinges(1)%omega(i)=bvec(i+3,1)
        end do

        do i=1,3
            fiber_hinges(2)%v_i(i)=bvec(i+6,1)
        end do
        
        !print *, ' bvec ', bvec
    else
        
        do i=1, ubound (fiber_hinges,1)-1
            do j=1,3
                fiber_hinges(i)%v_i(j)=bvec(9*(i-1)+j,1)
                fiber_hinges(i)%omega(j)=bvec(9*(i-1)+j+3,1)
            end do
        end do

        n=ubound(bvec,1)
        l=ubound(fiber_hinges,1)
        do i=-2,0
            fiber_hinges(l)%v_i(i+3)=bvec(n+i,1)
        end do
    end if

        do i=1, ubound (fiber_hinges,1)-1                      !2018/09/05  修正
           fiber_hinges(i)%v_i= fiber_hinges(i)%v_i_old + 0.5*(fiber_hinges(i)%v_i - fiber_hinges(i)%v_i_old)    !2018/09/05  修正
        end do
        
        do i=1, ubound (fiber_hinges,1)-2                      !2018/09/05  修正
           fiber_hinges(i)%omega= fiber_hinges(i)%omega_old + 0.5*(fiber_hinges(i)%omega - fiber_hinges(i)%omega_old)    !2018/09/05  修正
        end do        
    
    
    !do j=1,ubound(bvec,1)
    !        write(*,"(e8.2)") bvec(j,1)
    !    enddo
    
    if ((ubound(fiber_hinges,1)-1) .LT. minSegs) then
        !print *, " Regular Solver "
        deallocate(Amat)
    else
        !print *, " Banded Solver "
        deallocate(AB)
    end if
    
    deallocate(bvec)
    !print *," Motion matrix 2"

end subroutine motion_matrix

!*******************************************************************
subroutine mini_mat(mat, vec, conn, conn2, drag_coeff_vel, drag_coeff_omega)
    implicit none
    real(8), dimension(9,15) :: mat
    real(8), dimension(9,1)  :: vec
    real(8)                  :: drag_coeff_vel, drag_coeff_omega
    type(rod)                :: conn, conn2

    mat      =  0
    vec      =  0

    !********************
    mat(1, 4)=  1
    mat(1, 8)=  conn%r(3)
    mat(1, 9)= -conn%r(2)
    mat(1,13)= -1
    !********************
    mat(2, 5)=  1
    mat(2, 7)= -conn%r(3)
    mat(2, 9)= +conn%r(1)
    mat(2,14)= -1
    !********************
    mat(3, 6)=  1
    mat(3, 7)= +conn%r(2)
    mat(3, 8)= -conn%r(1)
    mat(3,15)= -1
    !********************
    !mat(4, 4)=  1
    mat(4, 1)=  1
    mat(4, 4)= -drag_coeff_vel*conn%ave_viscosity*conn%nbr_beads
    mat(4, 8)= -drag_coeff_vel*conn%ave_viscosity*conn%r_sum(3)
    mat(4, 9)= +drag_coeff_vel*conn%ave_viscosity*conn%r_sum(2)
    mat(4,10)= -1
    !------------
    vec(4, 1)= -drag_coeff_vel*conn%ave_viscosity*conn%u_fluid_sum(1)-conn%F_excl_vol(1)
    !********************
    !mat(5, 5)=  1
    mat(5, 2)=  1
    mat(5, 5)= -drag_coeff_vel*conn%ave_viscosity*conn%nbr_beads
    mat(5, 7)= +drag_coeff_vel*conn%ave_viscosity*conn%r_sum(3)
    mat(5, 9)= -drag_coeff_vel*conn%ave_viscosity*conn%r_sum(1)
    mat(5,11)= -1
    !------------
    vec(5, 1)= -drag_coeff_vel*conn%ave_viscosity*conn%u_fluid_sum(2)-conn%F_excl_vol(2)
    !********************
    !mat(6,6 )=  1
    mat(6, 3)=  1
    mat(6, 6)= -drag_coeff_vel*conn%ave_viscosity*conn%nbr_beads
    mat(6, 7)= -drag_coeff_vel*conn%ave_viscosity*conn%r_sum(2)
    mat(6, 8)= +drag_coeff_vel*conn%ave_viscosity*conn%r_sum(1)
    mat(6,12)= -1
    !------------
    vec(6, 1)= -drag_coeff_vel*conn%ave_viscosity*conn%u_fluid_sum(3)-conn%F_excl_vol(3)
    !********************
    mat(7,7)=1
    mat(7, 5)=  drag_coeff_vel*conn%ave_viscosity  *   conn%r_sum(3)
    mat(7, 6)= -drag_coeff_vel*conn%ave_viscosity  *   conn%r_sum(2)
    mat(7, 7)= -drag_coeff_omega*conn%ave_viscosity*   conn%nbr_beads&
        -drag_coeff_vel*conn%ave_viscosity  *  (conn%r_prod_sum(3,3)+conn%r_prod_sum(2,2))
    mat(7, 8)=  drag_coeff_vel*conn%ave_viscosity  *   conn%r_prod_sum(1,2)
    mat(7, 9)= +drag_coeff_vel*conn%ave_viscosity  *   conn%r_prod_sum(1,3)!CORR
    mat(7,11)= +conn%r(3)
    mat(7,12)= -conn%r(2)
    !------------
    vec(7, 1)= -drag_coeff_omega*conn%ave_viscosity* conn%omega_fluid_sum(1)&
        +drag_coeff_vel*conn%ave_viscosity  *(conn%r_times_u_sum(3,2)-conn%r_times_u_sum(2,3))&
        -conn%T(1)-conn%T_excl_vol(1)
    !********************
    !mat(8,8)=1
    mat(8, 4)= -drag_coeff_vel*conn%ave_viscosity  * conn%r_sum(3)
    mat(8, 6)=  drag_coeff_vel*conn%ave_viscosity  * conn%r_sum(1)
    mat(8, 7)= +drag_coeff_vel*conn%ave_viscosity  * conn%r_prod_sum(1,2)
    mat(8, 8)= -drag_coeff_omega*conn%ave_viscosity* conn%nbr_beads&
        -drag_coeff_vel*conn%ave_viscosity  *(conn%r_prod_sum(3,3)+conn%r_prod_sum(1,1))
    mat(8, 9)=  drag_coeff_vel*conn%ave_viscosity  * (conn%r_prod_sum(2,3))
    mat(8,10)= -conn%r(3)
    mat(8,12)= +conn%r(1)
    !------------
    vec(8, 1)= -drag_coeff_omega*conn%ave_viscosity* conn%omega_fluid_sum(2)&
        +drag_coeff_vel*conn%ave_viscosity  *(conn%r_times_u_sum(1,3)-conn%r_times_u_sum(3,1))&
        -conn%T(2)-conn%T_excl_vol(2)
    !********************
    !mat(9,9)=1
    mat(9, 4)=  drag_coeff_vel*conn%ave_viscosity  * conn%r_sum(2)
    mat(9, 5)= -drag_coeff_vel*conn%ave_viscosity  * conn%r_sum(1)
    mat(9, 7)= +drag_coeff_vel*conn%ave_viscosity  * conn%r_prod_sum(1,3)
    mat(9, 8)= +drag_coeff_vel*conn%ave_viscosity  * conn%r_prod_sum(2,3)
    mat(9, 9)= -drag_coeff_omega*conn%ave_viscosity* conn%nbr_beads&
        -drag_coeff_vel*conn%ave_viscosity  *(conn%r_prod_sum(1,1)+conn%r_prod_sum(2,2))
    mat(9,10)= +conn%r(2)
    mat(9,11)= -conn%r(1)
    !------------
    vec(9, 1)= -drag_coeff_omega*conn%ave_viscosity* conn%omega_fluid_sum(3)&
        -drag_coeff_vel*conn%ave_viscosity  *(conn%r_times_u_sum(1,2)-conn%r_times_u_sum(2,1))&
        -conn%T(3)-conn%T_excl_vol(3)
    !********************
    !print *, "conn%r", conn%r
    !print *, "drag omega", drag_omega
    !print *, "drag vel " , drag_vel
    !print *, "conn_nbr beads", conn%nbr_beads
    !print *, "conn %r_sum", conn%r_sum

end subroutine mini_mat

subroutine mini_mat_tensor(mat, vec, conn,  drag_coeff_vel, drag_coeff_omega)
    implicit none
    real(8), dimension(9,15) :: mat
    real(8), dimension(9,1)  :: vec
    real(8)                  :: drag_coeff_vel, drag_coeff_omega
    type(rod)                :: conn

    mat      =  0
    vec      =  0

    !********************
    mat(1, 4)=  1
    mat(1, 8)=  conn%r(3)
    mat(1, 9)= -conn%r(2)
    mat(1,13)= -1
    !********************
    mat(2, 5)=  1
    mat(2, 7)= -conn%r(3)
    mat(2, 9)= +conn%r(1)
    mat(2,14)= -1
    !********************
    mat(3, 6)=  1
    mat(3, 7)= +conn%r(2)
    mat(3, 8)= -conn%r(1)
    mat(3,15)= -1
    !********************
    !mat(4, 4)=  1
    mat(4, 1)=  1
    mat(4, 4)= -conn%A(1,1)
    mat(4, 5)= -conn%A(1,2)
    mat(4, 6)= -conn%A(1,3)
    mat(4, 7)=  0.5*conn%A(1,2) *conn%r(3)  - 0.5*conn%A(1,3) *conn%r(2)
    mat(4, 8)=  -0.5*conn%A(1,1) *conn%r(3)  + 0.5*conn%A(1,3) *conn%r(1)
    mat(4, 9)=  0.5*conn%A(1,1) *conn%r(2)  - 0.5*conn%A(1,2) *conn%r(1)
    mat(4,10)= -1
    !------------
    vec(4, 1)= -(conn%A(1,1)*conn%u_oo(1)+ conn%A(1,2)*conn%u_oo(2)+conn%A(1,3)*conn%u_oo(3) ) -conn%F_excl_vol(1)
    !********************
    !mat(5, 5)=  1
    mat(5, 2)=  1
    mat(5, 4)= -conn%A(2,1)
    mat(5, 5)= -conn%A(2,2)
    mat(5, 6)= -conn%A(2,3)
    mat(5, 7)=  0.5*conn%A(2,2) *conn%r(3)  - 0.5*conn%A(2,3) *conn%r(2)
    mat(5, 8)= -0.5*conn%A(2,1) *conn%r(3)  + 0.5*conn%A(2,3) *conn%r(1)
    mat(5, 9)=  0.5*conn%A(2,1) *conn%r(2)  - 0.5*conn%A(2,2) *conn%r(1)
    mat(5,11)= -1
    !------------
    vec(5, 1)= -(conn%A(2,1)*conn%u_oo(1)+ conn%A(2,2)*conn%u_oo(2)+conn%A(2,3)*conn%u_oo(3) )-conn%F_excl_vol(2)
    !********************
    !mat(6,6 )=  1
    mat(6, 3)=  1
    mat(6, 4)= -conn%A(3,1)
    mat(6, 5)= -conn%A(3,2)
    mat(6, 6)= -conn%A(3,3)
    mat(6, 7)=  0.5*conn%A(3,2) *conn%r(3)  - 0.5*conn%A(3,3) *conn%r(2)
    mat(6, 8)= -0.5*conn%A(3,1) *conn%r(3)  + 0.5*conn%A(3,3) *conn%r(1)
    mat(6, 9)=  0.5*conn%A(3,1) *conn%r(2)  - 0.5*conn%A(3,2) *conn%r(1)
    mat(6,12)= -1
    !------------
    vec(6, 1)= -(conn%A(3,1)*conn%u_oo(1)+ conn%A(3,2)*conn%u_oo(2)+conn%A(3,3)*conn%u_oo(3) )-conn%F_excl_vol(3)
    !********************
    mat(7, 7)=  -conn%C(1,1)
    mat(7, 8)=  -conn%C(1,2)
    mat(7, 9)=  -conn%C(1,3)
    mat(7,2)= +conn%r(3)*0.5
    mat(7,3)= -conn%r(2)*0.5
    mat(7,11)= +conn%r(3)*0.5
    mat(7,12)= -conn%r(2)*0.5
    !------------
    vec(7, 1)= -(conn%C(1,1)*conn%omega_oo(1)+ conn%C(1,2)*conn%omega_oo(2)+conn%C(1,3)*conn%omega_oo(3)+ conn%H(1) )  -conn%T(1)-conn%T_excl_vol(1)
    !********************
    !mat(8,8)=1

    mat(8, 7)=  -conn%C(2,1)
    mat(8, 8)=  -conn%C(2,2)
    mat(8, 9)=  -conn%C(2,3)
    mat(8,1)= -conn%r(3)*0.5
    mat(8,3)= +conn%r(1)*0.5
    mat(8,10)= -conn%r(3)*0.5
    mat(8,12)= +conn%r(1)*0.5
    !------------
    vec(8, 1)= -(conn%C(2,1)*conn%omega_oo(1)+ conn%C(2,2)*conn%omega_oo(2)+conn%C(2,3)*conn%omega_oo(3)+ conn%H(2) )  -conn%T(2)-conn%T_excl_vol(2)
    
    !********************
    !mat(9,9)=1
    mat(9, 7)=  -conn%C(3,1)
    mat(9, 8)=  -conn%C(3,2)
    mat(9, 9)=  -conn%C(3,3)
    mat(9,10)= +conn%r(2)*0.5
    mat(9,11)= -conn%r(1)*0.5
    mat(9,1)= +conn%r(2)*0.5
    mat(9,2)= -conn%r(1)*0.5
    !------------
    vec(9, 1)= -(conn%C(3,1)*conn%omega_oo(1)+ conn%C(3,2)*conn%omega_oo(2)+conn%C(3,3)*conn%omega_oo(3)+ conn%H(3) ) -conn%T(3)-conn%T_excl_vol(3)
    !********************
    !print *, "conn%r", conn%r
    !print *, "drag omega", drag_omega
    !print *, "drag vel " , drag_vel
    !print *, "conn_nbr beads", conn%nbr_beads
    !print *, "conn %r_sum", conn%r_sum

end subroutine mini_mat_tensor


end module m_Motion
