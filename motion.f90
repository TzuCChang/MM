!=============================================================================
module m_Motion
    
use m_DataStructures
use MATRIX_SOLVERS
use omp_lib
implicit none
contains
!=============================================================================

subroutine motion_fiber( fibers, hinges, simParameters )
type(simulationParameters) :: simParameters  
type(fiber) , dimension(:) :: fibers
type(rod), dimension(:)    :: hinges
integer(8)                 :: i
real(8)                    :: r_bead, viscosity, timex, timeA

    r_bead= simParameters%r_fiber   !2018/10/10 add

    do i=1, ubound(fibers,1)   
        
        call motion_matrix( hinges(fibers(i)%first_hinge:&
                            fibers(i)%first_hinge+fibers(i)%nbr_hinges-1),&
                            r_bead, simParameters )                            !2018/12/02
    end do

end subroutine motion_fiber

!=============================================================================

subroutine motion_matrix( fiber_hinges, r_bead, simParameters )  !2018/12/02
implicit none
type(simulationParameters)           :: simParameters
type (rod), dimension(:)             :: fiber_hinges

real(8), dimension(:,:), allocatable :: Amat, bvec, xvec, AB
real(8), dimension(9,15)             :: mat
real(8), dimension(9,1)              :: vec
real(8)                              :: r_bead, finish2, start2, timex
integer(8)                           :: i, k, j, l, kl, ku, minSegs 
integer(8)                           :: n, ii, jj


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

        do i=1, ubound (fiber_hinges,1)-1                      !2018/09/05  タ
           fiber_hinges(i)%v_i_OLD= fiber_hinges(i)%v_i        !2018/09/05  タ
        end do
        
        do i=1, ubound (fiber_hinges,1)-2                      !2018/09/05  タ
           fiber_hinges(i)%omega_OLD= fiber_hinges(i)%omega    !2018/09/05  タ
        end do    

    !print *, "DIMENSION OF FIBER HINGES", ubound (fiber_hinges,1)
    if ((ubound(fiber_hinges,1)-1) .LT. minSegs) then
        if (ubound (fiber_hinges,1)==2) then

if( simParameters%TensorOrBeads .eq. 1 ) then       !2018/12/02
    
            call mini_mat_tensor( mat, vec, fiber_hinges(1) )
else
            call mini_mat( mat, vec, fiber_hinges(1), fiber_hinges(2) )
end if
            !mini_mat_tensor
            Amat(1:9,1:6) =mat(1:9, 4:9 )
            Amat(1:9,7:9) =mat(1:9,13:15)

            bvec(:,1)=vec(:,1)
        else
            do i=1, ubound (fiber_hinges,1)-1
                
if( simParameters%TensorOrBeads .eq. 1 ) then       !2018/12/02

                call mini_mat_tensor( mat, vec, fiber_hinges(i) )
else
                call mini_mat( mat, vec, fiber_hinges(i), fiber_hinges(i+1) )
end if               
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

if( simParameters%TensorOrBeads .eq. 1 ) then       !2018/12/02

            call mini_mat_tensor( mat, vec, fiber_hinges(1) )
else
            call mini_mat( mat, vec, fiber_hinges(1), fiber_hinges(2) )
endif
            call copyToBanded(AB, KL, KU, N, 1, 9, 1, 6, mat(1:9, 4:9 ))
            call copyToBanded(AB, KL, KU, N, 1, 9, 7, 9, mat(1:9,13:15))
            !Amat(1:9,1:6) =mat(1:9, 4:9 )
            !Amat(1:9,7:9) =mat(1:9,13:15)

            bvec(:,1)=vec(:,1)
        else

            do i=1, ubound (fiber_hinges,1)-1

if( simParameters%TensorOrBeads .eq. 1 ) then       !2018/12/02

                call mini_mat_tensor( mat, vec, fiber_hinges(i) )
else
                call mini_mat( mat, vec, fiber_hinges(i), fiber_hinges(i+1) )
endif      
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

        do i=1, ubound (fiber_hinges,1)-1                      !2018/09/05  タ
           fiber_hinges(i)%v_i= fiber_hinges(i)%v_i_old + 0.5*(fiber_hinges(i)%v_i - fiber_hinges(i)%v_i_old)    !2018/09/05  タ
        end do
        
        do i=1, ubound (fiber_hinges,1)-2                      !2018/09/05  タ
           fiber_hinges(i)%omega= fiber_hinges(i)%omega_old + 0.5*(fiber_hinges(i)%omega - fiber_hinges(i)%omega_old)    !2018/09/05  タ
        end do        
    
    !do j=1,ubound(bvec,1)
    !        write(*,"(e8.2)") bvec(j,1)
    !    enddo
    
    if ((ubound(fiber_hinges,1)-1) .LT. minSegs) then
        deallocate(Amat)
    else
        deallocate(AB)
    end if
    
    deallocate(bvec)

end subroutine motion_matrix

!=============================================================================

subroutine mini_mat( mat, vec, conn, conn2 )
implicit none
real(8), dimension(9,15) :: mat
real(8), dimension(9,1)  :: vec
type(rod)                :: conn, conn2
    
    mat=  0
    vec=  0

    !********************
    mat(1,4)=   1               !2018/12/01, 跑计[X1,u1,w1,X2,u2],  X1(force)竚(1,2,3), u1竚(4,5,6)
    mat(1,8)=   conn%r(3)       !2018/12/01, 跑计[X1,u1,w1,X2,u2],  w1(force)竚(7,8,9), X2竚(10,11,12)
    mat(1,9)=  -conn%r(2)       !2018/12/01, 跑计[X1,u1,w1,X2,u2],  u2(force)竚(13,14,15)
    mat(1,13)= -1               !2018/12/01, 跑计[X1,u1,w1,X2,u2],  u2(force)竚(13,14,15)
    !********************
    mat(2,5)=   1
    mat(2,7)=  -conn%r(3)
    mat(2,9)=   conn%r(1)
    mat(2,14)= -1
    !********************
    mat(3,6)=   1
    mat(3,7)=   conn%r(2)
    mat(3,8)=  -conn%r(1)
    mat(3,15)= -1
    !********************
    mat(4,1)=   1
    mat(4,4)=  -conn%c1_nbeads    !2018/12/07   タ
    mat(4,8)=  -conn%r_sum(3)
    mat(4,9)=   conn%r_sum(2)
    mat(4,10)= -1
    !------------
    vec(4, 1)= -( conn%uoo_sum(1) + conn%F_excl(1) )
    !********************
    mat(5, 2)=  1
    mat(5, 5)= -conn%c1_nbeads    !2018/12/07   タ
    mat(5, 7)=  conn%r_sum(3)
    mat(5, 9)= -conn%r_sum(1)
    mat(5,11)= -1
    !------------
    vec(5, 1)= -( conn%uoo_sum(2) + conn%F_excl(2) )
    !********************
    mat(6,3)=   1
    mat(6,6)=  -conn%c1_nbeads    !2018/12/07   タ
    mat(6,7)=  -conn%r_sum(2)
    mat(6,8)=   conn%r_sum(1)
    mat(6,12)= -1
    !------------
    vec(6,1)=  -( conn%uoo_sum(3) + conn%F_excl(3) )
    !********************
    mat(7,5)=   conn%rk_sum(3)
    mat(7,6)=  -conn%rk_sum(2)
    mat(7,7)=  -conn%r_prod_sum(1,1) - conn%c2_nbeads          !2018/12/07   タ
    mat(7,8)=   conn%r_prod_sum(1,2)
    mat(7,9)=   conn%r_prod_sum(1,3)
    mat(7,2)=   conn%r(3)*0.5d0                                !2018/12/02   タ
    mat(7,3)=  -conn%r(2)*0.5d0                                !2018/12/02   タ
    mat(7,11)=  conn%r(3)*0.5d0                                !2018/12/02   タ
    mat(7,12)= -conn%r(2)*0.5d0                                !2018/12/02   タ
    !------------
    vec(7,1)=  -( conn%T(1) + conn%T_excl(1) )&                !2018/12/07   タ
               -( conn%rk_vel(1) + conn%omega_sum(1) )         !2018/12/07   タ 
    !********************
    mat(8,4)=  -conn%rk_sum(3)
    mat(8,6)=   conn%rk_sum(1)
    mat(8,7)=   conn%r_prod_sum(2,1)
    mat(8,8)=  -conn%r_prod_sum(2,2) - conn%c2_nbeads          !2018/12/07   タ
    mat(8,9)=   conn%r_prod_sum(2,3)
    mat(8,1)=  -conn%r(3)*0.5d0                                !2018/12/02   タ
    mat(8,3)=   conn%r(1)*0.5d0                                !2018/12/02   タ
    mat(8,10)= -conn%r(3)*0.5d0                                !2018/12/02   タ
    mat(8,12)=  conn%r(1)*0.5d0                                !2018/12/02   タ
    !------------
    vec(8, 1)= -( conn%T(2) + conn%T_excl(2) )&                !2018/12/07   タ
               -( conn%rk_vel(2) + conn%omega_sum(2) )         !2018/12/07   タ
    !********************
    mat(9, 4)=   conn%rk_sum(2)
    mat(9, 5)=  -conn%rk_sum(1)
    mat(9, 7)=   conn%r_prod_sum(3,1)
    mat(9, 8)=   conn%r_prod_sum(3,2)
    mat(9, 9)=  -conn%r_prod_sum(3,3) - conn%c2_nbeads         !2018/12/07   タ
    mat(9,1)=    conn%r(2)*0.5d0                               !2018/12/07   タ
    mat(9,2)=   -conn%r(1)*0.5d0                               !2018/12/07   タ    
    mat(9,10)=   conn%r(2)*0.5d0                               !2018/12/07   タ
    mat(9,11)=  -conn%r(1)*0.5d0                               !2018/12/07   タ
    !------------
    vec(9, 1)=  -( conn%T(3)  + conn%T_excl(3) )&              !2018/12/07   タ 
                -( conn%rk_vel(3) + conn%omega_sum(3) )        !2018/12/07   タ
    !********************

end subroutine mini_mat

!=============================================================================

subroutine mini_mat_tensor( mat, vec, conn )  !2018/12/07  タ
    implicit none
    real(8), dimension(9,15) :: mat
    real(8), dimension(9,1)  :: vec
    type(rod)                :: conn

    mat=  0
    vec=  0

    !********************
    mat(1, 4)=  1
    mat(1, 8)=  conn%r(3)
    mat(1, 9)= -conn%r(2)
    mat(1,13)= -1
    !********************
    mat(2, 5)=  1
    mat(2, 7)= -conn%r(3)
    mat(2, 9)=  conn%r(1)
    mat(2,14)= -1
    !********************
    mat(3, 6)=  1
    mat(3, 7)=  conn%r(2)
    mat(3, 8)= -conn%r(1)
    mat(3,15)= -1
    !********************
    mat(4, 1)=  1
    mat(4, 4)= -conn%A(1,1)
    mat(4, 5)= -conn%A(1,2)
    mat(4, 6)= -conn%A(1,3)

    mat(4, 7)= -0.5d0*( -conn%A(1,2)*conn%r(3) + conn%A(1,3)*conn%r(2) )
    mat(4, 8)= -0.5d0*(  conn%A(1,1)*conn%r(3) - conn%A(1,3)*conn%r(1) )
    mat(4, 9)= -0.5d0*( -conn%A(1,1)*conn%r(2) + conn%A(1,2)*conn%r(1) )
    mat(4,10)= -1
    !------------
    vec(4, 1)= -( conn%Auf_oo(1) + conn%F_excl(1) )
    !********************
    mat(5, 2)=  1
    mat(5, 4)= -conn%A(2,1)
    mat(5, 5)= -conn%A(2,2)
    mat(5, 6)= -conn%A(2,3)
    
    mat(5, 7)= -0.5d0*( -conn%A(2,2)*conn%r(3) + conn%A(2,3)*conn%r(2) )
    mat(5, 8)= -0.5d0*(  conn%A(2,1)*conn%r(3) - conn%A(2,3)*conn%r(1) )
    mat(5, 9)= -0.5d0*( -conn%A(2,1)*conn%r(2) + conn%A(2,2)*conn%r(1) )
    mat(5,11)= -1
    !------------
    vec(5, 1)= -( conn%Auf_oo(2) + conn%F_excl(2) )
    !********************
    mat(6, 3)=  1
    mat(6, 4)= -conn%A(3,1)
    mat(6, 5)= -conn%A(3,2)
    mat(6, 6)= -conn%A(3,3)

    mat(6, 7)= -0.5d0*( -conn%A(3,2)*conn%r(3) + conn%A(3,3)*conn%r(2) )
    mat(6, 8)= -0.5d0*(  conn%A(3,1)*conn%r(3) - conn%A(3,3)*conn%r(1) )
    mat(6, 9)= -0.5d0*( -conn%A(3,1)*conn%r(2) + conn%A(3,2)*conn%r(1) )
    mat(6,12)= -1
    !------------
    vec(6, 1)= -( conn%Auf_oo(3) + conn%F_excl(3) )
    !********************
    mat(7, 7)= -conn%C(1,1)
    mat(7, 8)= -conn%C(1,2)
    mat(7, 9)= -conn%C(1,3)

    mat(7,2)=   0.5d0*conn%r(3)
    mat(7,3)=  -0.5d0*conn%r(2)
    mat(7,11)=  0.5d0*conn%r(3)
    mat(7,12)= -0.5d0*conn%r(2)
    !------------
    vec(7, 1)= -( conn%omega_oo(1) + conn%H(1) + conn%T(1) + conn%T_excl(1) )
    !********************
    mat(8, 7)=  -conn%C(2,1)
    mat(8, 8)=  -conn%C(2,2)
    mat(8, 9)=  -conn%C(2,3)

    mat(8,1)=   -0.5d0*conn%r(3)
    mat(8,3)=    0.5d0*conn%r(1)
    mat(8,10)=  -0.5d0*conn%r(3)
    mat(8,12)=   0.5d0*conn%r(1)
    !------------
    vec(8, 1)= -( conn%omega_oo(2) + conn%H(2) + conn%T(2) + conn%T_excl(2) )
    !********************
    mat(9, 7)=  -conn%C(3,1)
    mat(9, 8)=  -conn%C(3,2)
    mat(9, 9)=  -conn%C(3,3)

    mat(9,1)=    0.5d0*conn%r(2)
    mat(9,2)=   -0.5d0*conn%r(1)    
    mat(9,10)=   0.5d0*conn%r(2)
    mat(9,11)=  -0.5d0*conn%r(1)
    !------------
    vec(9, 1)= -( conn%omega_oo(3) + conn%H(3) + conn%T(3) + conn%T_excl(3) )
    !********************

end subroutine mini_mat_tensor

!=============================================================================

subroutine mini_mat_tensor_02( mat, vec, conn )  !2018/12/07  タ
    implicit none
    real(8), dimension(9,15) :: mat
    real(8), dimension(9,1)  :: vec
    type(rod)                :: conn

    mat=  0
    vec=  0

    !********************
    mat(1, 4)=  1
    mat(1, 8)=  conn%r(3)
    mat(1, 9)= -conn%r(2)
    mat(1,13)= -1
    !********************
    mat(2, 5)=  1
    mat(2, 7)= -conn%r(3)
    mat(2, 9)=  conn%r(1)
    mat(2,14)= -1
    !********************
    mat(3, 6)=  1
    mat(3, 7)=  conn%r(2)
    mat(3, 8)= -conn%r(1)
    mat(3,15)= -1
    !********************
    mat(4, 1)=  1
    mat(4, 4)= -0.5d0*conn%A(1,1)
    mat(4, 5)= -0.5d0*conn%A(1,2)
    mat(4, 6)= -0.5d0*conn%A(1,3)

    mat(4,10)= -1

    mat(4,13)= -0.5d0*conn%A(1,1)
    mat(4,14)= -0.5d0*conn%A(1,2)
    mat(4,15)= -0.5d0*conn%A(1,3)
    !------------
    vec(4, 1)= -( conn%Auf_oo(1) + conn%F_excl(1) )
    !********************
    mat(5, 2)=  1
    mat(5, 4)= -0.5d0*conn%A(2,1)
    mat(5, 5)= -0.5d0*conn%A(2,2)
    mat(5, 6)= -0.5d0*conn%A(2,3)
    
    mat(5,13)= -0.5d0*conn%A(2,1)
    mat(5,14)= -0.5d0*conn%A(2,2)
    mat(5,15)= -0.5d0*conn%A(2,3)
    mat(5,11)= -1
    !------------
    vec(5, 1)= -( conn%Auf_oo(2) + conn%F_excl(2) )
    !********************
    mat(6, 3)=  1
    mat(6, 4)= -0.5d0*conn%A(3,1)
    mat(6, 5)= -0.5d0*conn%A(3,2)
    mat(6, 6)= -0.5d0*conn%A(3,3)
    
    mat(6,13)= -0.5d0*conn%A(3,1)
    mat(6,14)= -0.5d0*conn%A(3,2)
    mat(6,15)= -0.5d0*conn%A(3,3)
    mat(6,12)= -1
    !------------
    vec(6, 1)= -( conn%Auf_oo(3) + conn%F_excl(3) )
    !********************
    mat(7, 7)= -conn%C(1,1)
    mat(7, 8)= -conn%C(1,2)
    mat(7, 9)= -conn%C(1,3)

    mat(7,2)=   0.5d0*conn%r(3)
    mat(7,3)=  -0.5d0*conn%r(2)
    mat(7,11)=  0.5d0*conn%r(3)
    mat(7,12)= -0.5d0*conn%r(2)
    !------------
    vec(7, 1)= -( conn%omega_oo(1) + conn%H(1) + conn%T(1) + conn%T_excl(1) )
    !********************
    mat(8, 7)=  -conn%C(2,1)
    mat(8, 8)=  -conn%C(2,2)
    mat(8, 9)=  -conn%C(2,3)

    mat(8,1)=   -0.5d0*conn%r(3)
    mat(8,3)=    0.5d0*conn%r(1)
    mat(8,10)=  -0.5d0*conn%r(3)
    mat(8,12)=   0.5d0*conn%r(1)
    !------------
    vec(8, 1)= -( conn%omega_oo(2) + conn%H(2) + conn%T(2) + conn%T_excl(2) )
    !********************
    mat(9, 7)=  -conn%C(3,1)
    mat(9, 8)=  -conn%C(3,2)
    mat(9, 9)=  -conn%C(3,3)

    mat(9,1)=    0.5d0*conn%r(2)
    mat(9,2)=   -0.5d0*conn%r(1)    
    mat(9,10)=   0.5d0*conn%r(2)
    mat(9,11)=  -0.5d0*conn%r(1)
    !------------
    vec(9, 1)= -( conn%omega_oo(3) + conn%H(3) + conn%T(3) + conn%T_excl(3) )
    !********************

end subroutine mini_mat_tensor_02

!=============================================================================
subroutine copyToBanded(AB, KL, KU, nbRows, i1, i2, j1, j2, A)
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

end module m_Motion
!=============================================================================