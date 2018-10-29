!  DllTest.f90 
!
!  FUNCTIONS/SUBROUTINES exported from DllTest.dll:
!  DllTest - subroutine 
!
subroutine DllTest(N) BIND(C,NAME="DllTest") 
implicit none
  ! Expose subroutine DllTest to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::DllTest
  Integer N
  print *, "Number in Fortran:", N
   print *, "this really works"
  
  ! Variables
  

 ! Body of DllTest

end subroutine DllTest


subroutine DllTest2(N) BIND(C,NAME="DllTest2") 
implicit none
  ! Expose subroutine DllTest to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::DllTest2
  Integer N
  print *, "Number in Fortran:", N*4
  
  ! Variables
  

 ! Body of DllTest

end subroutine DllTest2

subroutine Solve1(m,KL,KU,AB,B) BIND(C,NAME="Solve1") 
implicit none
  ! Expose subroutine DllTest to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::Solve1
  Integer                       m,i,INFO,NRHS,KL,KU
  INTEGER, DIMENSION(m)                           :: IPIV
  Real, dimension((KL+KU+1)*m) ::   AB
  Real, dimension(m) ::             B
  
  
  
  
  CALL DGBSV(m, KL, KU, 1 , AB , 2*KL+KU+1, IPIV,B, m, INFO)
  
  
  !print *, "Number in Fortran:", n 
  ! Variables
  
  !CALL DGESV(n,NRHS, A, n, IPIV, B, n, INFO)

end subroutine Solve1


subroutine Solve2(m,n,KL,KU,AB,B,IPIV,INFO) BIND(C,NAME="Solve2") 
use mkl_service
implicit none
  ! Expose subroutine DllTest to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::Solve2
  Integer                       m,n,i,INFO,NRHS,KL,KU,nt
  INTEGER, DIMENSION(m)                           :: IPIV
  Real, dimension((2*KL+KU+1)*n) ::   AB
  Real, dimension(m) ::             B
  
  !n = number of columns
  nt = mkl_get_max_threads()
  ! print *, "Number of threads:", nt 
  call mkl_set_num_threads( 4 )
  
  !CALL DGBTRS ('N', n, KL, KU, 1 , AB , 2*KL+KU+1, IPIV,B, m, INFO)
  
  !CALL DGBTRF (m, n, KL, KU, AB , 2*KL+KU+1, IPIV, INFO)
  call dgbsv( m, KL, KU, 1, AB, 2*KL+KU+1, IPIV, B, m, INFO )
  !print *, "Number in Fortran:", n 
  ! Variables
  
  !CALL DGESV(n,NRHS, A, n, IPIV, B, n, INFO)

end subroutine Solve2

subroutine Solve3(m,n,KL,KU,AB,B,IPIV,INFO) BIND(C,NAME="Solve3") 
implicit none
  ! Expose subroutine DllTest to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::Solve3
  Integer                       m,n,i,INFO,NRHS,KL,KU
  INTEGER, DIMENSION(m)                           :: IPIV
  Real, dimension((2*KL+KU+1)*n) ::   AB
  Real, dimension(m) ::             B
  
  !n = number of columns
  
  
  
  CALL DGBTRS ('N', n, KL, KU, 1 , AB , 2*KL+KU+1, IPIV,B, m, INFO)
  
  !CALL DGBTRF (m, n, KL, KU, AB , 2*KL+KU+1, IPIV, INFO)
  !call dgbsv( m, KL, KU, 1, AB, 2*KL+KU+1, IPIV, B, m, INFO )
  !print *, "Number in Fortran:", n 
  ! Variables
  
  !CALL DGESV(n,NRHS, A, n, IPIV, B, n, INFO)

end subroutine Solve3

subroutine Solve4(m,n,A,B,WORK,INFO) BIND(C,NAME="Solve4") 
use mkl_service
implicit none
  ! Expose subroutine DllTest to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::Solve4
  Integer                       m,n,L,INFO,nt
  Real, dimension(m*n) ::   A
  Real, dimension(m) ::             B, WORK

  !n = number of columns
  nt = mkl_get_max_threads()
  ! print *, "Number of threads:", nt 
  call mkl_set_num_threads( 4 )
  
  !CALL DGBTRS ('N', n, KL, KU, 1 , AB , 2*KL+KU+1, IPIV,B, m, INFO)
  
  !CALL DGBTRF (m, n, KL, KU, AB , 2*KL+KU+1, IPIV, INFO)
  call dgels('N', m, n, 1 , A , m, B, m,  WORK, -1, INFO )
 
  !print *, "Number in Fortran:", n 
  ! Variables
  
  !CALL DGESV(n,NRHS, A, n, IPIV, B, n, INFO)

end subroutine Solve4