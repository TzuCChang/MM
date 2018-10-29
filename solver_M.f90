
MODULE MATRIX_SOLVERS
IMPLICIT NONE
CONTAINS
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
SUBROUTINE SOLVER_0(nbr_rows, RBF_MAT, Vector_b, Vector_x)
INTEGER, DIMENSION(nbr_rows)                           :: IPIV
REAL(8), DIMENSION(nbr_rows, nbr_rows), INTENT(IN)     :: RBF_MAT
REAL(8), DIMENSION(nbr_rows, 1)       , INTENT(OUT)    :: Vector_x
REAL(8), DIMENSION(nbr_rows, 1)       , INTENT(IN)     :: Vector_b
REAL(4), DIMENSION(nbr_rows, nbr_rows)                 :: RBF_MAT_S
REAL(4), DIMENSION(nbr_rows, 1)                        :: Vector_x_S
INTEGER                                                :: INFO
INTEGER                                                :: nbr_rows

!Vector_x=Vector_b


RBF_MAT_S=real(RBF_MAT,4)
Vector_x_S=real(VECTOR_b,4)

 CALL SGESV(nbr_rows,1, RBF_MAT_S, nbr_rows, IPIV, Vector_x_S, nbr_rows, INFO)
Vector_x=real(VECTOR_x_s,8)
 !PRINT *, "DGESV", INFO
END SUBROUTINE SOLVER_0

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
SUBROUTINE SOLVER_1(nbr_rows, RBF_MAT, Vector_b)
INTEGER, DIMENSION(nbr_rows)                           :: IPIV
REAL(8), DIMENSION(nbr_rows, nbr_rows), INTENT(IN)     :: RBF_MAT
!REAL(8), DIMENSION(nbr_rows, 1)       , INTENT(OUT)    :: Vector_x
REAL(8), DIMENSION(nbr_rows, 1)              :: Vector_b
INTEGER                                                :: INFO
INTEGER                                                :: nbr_rows

!Vector_x=Vector_b

 CALL DGESV(nbr_rows,1, RBF_MAT, nbr_rows, IPIV, Vector_b, nbr_rows, INFO)
 !PRINT *, "DGESV", INFO
END SUBROUTINE SOLVER_1

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************

SUBROUTINE SOLVER_2(nbr_rows, RBF_MAT, Vector_b, Vector_x)
INTEGER, DIMENSION(nbr_rows)                           :: IPIV
REAL(8), DIMENSION(nbr_rows, nbr_rows), INTENT(IN)     :: RBF_MAT
!REAL(8), DIMENSION(nbr_rows, nbr_rows)                 :: RBF_MAT_TR
REAL(8), DIMENSION(nbr_rows, 1)       , INTENT(OUT)    :: Vector_x
REAL(8), DIMENSION(nbr_rows, 1)       , INTENT(IN)     :: Vector_b
REAL(8), DIMENSION(nbr_rows, 1)                        :: Vector_b_TR
INTEGER                                                :: nbr_rows
REAL(8), DIMENSION(nbr_rows, nbr_rows)                 :: AF
 CHARACTER(LEN=1)                                      :: EQUED
REAL(8), DIMENSION(nbr_rows)                           :: R, C
REAL(8)                                                :: rcond
REAL(8), DIMENSION(1)                                  :: FERR, BERR
REAL(8), DIMENSION(4*nbr_rows)                         :: WORK
INTEGER, DIMENSION(nbr_rows)                           :: IWORK
INTEGER                                                :: INFO

!RBF_MAT_TR=RBF_MAT
Vector_b_TR=Vector_b

!SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
!     $                   WORK, IWORK, INFO )

 CALL DGESVX('E',&
             'N',&
              nbr_rows,&
              1,&
              RBF_MAT,&
              nbr_rows,&
              AF,&
              nbr_rows,&
              IPIV,&
              EQUED,&
              R,&
              C,&
              Vector_b_TR,&
              nbr_rows,&
              Vector_x,&
              nbr_rows,&
              rcond,&
              FERR,&
              BERR,&
              WORK,&
              IWORK,&
              INFO)
              
 PRINT *, "DGESV ", INFO
 PRINT *, "RCOND", rcond
 PRINT *, "BERR", berr
 PRINT *, "FERR", ferr 
END SUBROUTINE SOLVER_2
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
SUBROUTINE SOLVER_3(nbr_rows, KL, KU, AB, Vector_b)
IMPLICIT NONE
INTEGER, DIMENSION(nbr_rows)                           :: IPIV
!REAL(8), DIMENSION(nbr_rows, nbr_rows), INTENT(IN)    :: RBF_MAT
REAL(8)                                                :: AB(:,:)
!REAL(8), DIMENSION(nbr_rows, 1)       , INTENT(OUT)   :: Vector_x
REAL(8), DIMENSION(nbr_rows, 1)            :: Vector_b
!REAL(8), DIMENSION(2*KL+KU+1,nbr_rows)                :: AB
INTEGER                                                :: INFO
INTEGER                                                :: nbr_rows
INTEGER                                                :: KL
INTEGER                                                :: KU
INTEGER                                                :: i, j
!Vector_x=Vector_b

!AB=0
!do j=1, UBOUND(RBF_MAT,2)
!	do i=max(1,j-KU), min(nbr_rows,j+KL)
!		AB(KL+KU+1+i-j,j) = RBF_MAT(i,j) !for max(1,j-KU)<=i<=min(N,j+KL)
!	end do
!end do
! Uncomment below to print the matrix AB in a file. It is easier to visualize. 

!write(7,*)  ' '
!print *, ' AB '
!    do j=1,nbr_rows
!        write(7,"(I8,A)",advance="no") j, ' '
!        enddo
!        write(7,*)  ' '
!    
!    do i=1,2*KL+KU+1
!        write(7,"(I2,A)",advance="no") i, ' '
!        do j=1,nbr_rows
!        write(7,"(e8.2,A)",advance="no") AB(i,j), ' '
!        enddo
!        write(7,*)  ' '
!    enddo

CALL DGBSV(nbr_rows,KL, KU, 1,AB, 2*KL+KU+1, IPIV, Vector_b, nbr_rows, INFO)
! !PRINT *, "DGESV", INFO
! print *, 
END SUBROUTINE SOLVER_3

SUBROUTINE SOLVER_4(nbr_rows, KL, KU, RBF_MAT, Vector_b)
IMPLICIT NONE
INTEGER, DIMENSION(nbr_rows)                           :: IPIV
REAL(8), DIMENSION(nbr_rows, nbr_rows), INTENT(IN)     :: RBF_MAT
!REAL(8)                                               :: AB(:,:)
!REAL(8), DIMENSION(nbr_rows, 1)       , INTENT(OUT)   :: Vector_x
REAL(8), DIMENSION(nbr_rows, 1)            :: Vector_b
REAL(8), DIMENSION(2*KL+KU+1,nbr_rows)                 :: AB
INTEGER                                                :: INFO
INTEGER                                                :: nbr_rows
INTEGER                                                :: KL
INTEGER                                                :: KU
INTEGER                                                :: i, j
!Vector_x=Vector_b

AB=0
do j=1, UBOUND(RBF_MAT,2)
	do i=max(1,j-KU), min(nbr_rows,j+KL)
		AB(KL+KU+1+i-j,j) = RBF_MAT(i,j) !for max(1,j-KU)<=i<=min(N,j+KL)
	end do
end do
! Uncomment below to print the matrix AB in a file. It is easier to visualize. 

write(7,*)  ' '
print *, ' AB '
    do j=1,nbr_rows
        write(7,"(I8,A)",advance="no") j, ' '
        enddo
        write(7,*)  ' '
    
    do i=1,2*KL+KU+1
        print  *, i, ' '
        write(7,"(I2,A)",advance="no") i, ' '
        do j=1,nbr_rows
        write(7,"(e9.2,A)",advance="no") AB(i,j), ' ' !2018/07/12  error e8.2
        enddo
        write(7,*)  ' '
    enddo

CALL DGBSV(nbr_rows,KL, KU, 1,AB, 2*KL+KU+1, IPIV, Vector_b, nbr_rows, INFO)
! !PRINT *, "DGESV", INFO
! print *, 
END SUBROUTINE SOLVER_4

!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
END MODULE MATRIX_SOLVERS
