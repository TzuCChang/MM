!====================================================================
module m_HingesDamage

!subroutine: hinges_damage

use m_DataStructures
use m_UtilityLib

implicit none
contains

subroutine hinges_damage(fibers, hinges, min_curv, r_fiber)
implicit none
type (fiber), dimension(:), allocatable :: fibers, fibers_temp
type (rod)  , dimension(:), allocatable :: hinges, hinges_temp
real(8)                    :: max_alpha, r_fiber, fac, alpha,min_curv,curv
integer                    :: k, i, j, m, n


fac=0.07*r_fiber ! �]������Ӿa����� �����w���]��(���g��o�X���ƭ�) �Y�O�Ӥj�|�y���t�Τ�í�w 

do i=1, ubound(hinges,1)
	hinges(i)%is_broken    =.false.
	hinges(i)%is_separated =.false.
end do

k=0

!print *, "number of fibers" ,  ubound(fibers,1)
do i=1, ubound(fibers,1)
	if(fibers(i)%nbr_hinges.ge.3) then ! .ge. = >= !��fiber(i)��hinge�j��3��
		do j=fibers(i)%first_hinge+1, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2 !�Oj=fiber(i)�Ĥ@��hinge+1��˼ƲĤG��hinge
            
            call find_curvature(hinges(j-1)%X_i, hinges(j)%X_i,hinges(j+1)%X_i, curv) !�I�sfind_curvature
            if (.not. isnan(curv) .and. curv.le.min_curv) then 
                !print *, "curv", curv
 				hinges(j)%is_broken = .true.		 				
				k=k+1
			end if
		end do
	end if
end do

do i=1, ubound(fibers,1)
	hinges(fibers(i)%first_hinge)%is_separated=.true.
end do

!PRINT *, "Monitor Hinges"
!do i=1, ubound(hinges,1)
!	print *, hinges(i)%X_i, hinges(i)%is_separated, hinges(i)%is_broken
!end do
!PRINT *, "Monitor Fibers"
!do i=1, ubound(fibers,1)
!	print*, fibers(i)%first_hinge, fibers(i)%nbr_hinges
!end do

allocate(hinges_temp(ubound(hinges,1)+k))
allocate(fibers_temp(ubound(fibers,1)+k))
!print *, "Neo Ubound *****", (ubound(fibers,1)+k), (ubound(hinges,1)+k)

k=0

m=1
do i=1, ubound(hinges,1)
	if(hinges(i)%is_broken .eqv. .true.) then
		k=k+1
	end if
	if((hinges(i)%is_broken .eqv. .true.).or.(hinges(i)%is_separated .eqv. .true.)) then
		fibers_temp(m)%first_hinge=k+i
		m=m+1
	end if

end do

do i=1, ubound(fibers_temp,1)-1
	fibers_temp(i)%nbr_hinges=fibers_temp(i+1)%first_hinge-fibers_temp(i)%first_hinge
end do
fibers_temp(ubound(fibers_temp,1))%nbr_hinges=ubound(hinges_temp,1)-fibers_temp(ubound(fibers_temp,1))%first_hinge+1

!do i=1, ubound(fibers_temp,1)
!	print *," NEON FIBERS", fibers_temp(i)%first_hinge
!end do

k=1
do i=1, ubound(hinges,1)
	if (hinges(i)%is_broken .eqv. .true.)then
        
		!print *, "BREAKAGE"      
        !adding to alocate all properties of hinges(i) to hinges_temp(i)
        !hinges_temp(k) = hinges(i)
        ! updating position to ensure not overlapping
        
hinges_temp(k)= hinges(i) !2018/08/04 �ץ����~ �]��hinges_temp(k) �O�� hinges(i) �������_���Ӫ�(k),�]���n�~��(i)���Ҧ���T,���O��m�n���ʤ@��       
                          !2018/08/04 �ץ����~ ����{���|�F�o�@��,�]��hinges_temp(k) �u�~�ӤF�y�� X_i �Mis_stationary
! i�_���ܦ�k, k+1
! fac�O0.07�ӥb�|�A�ҥHhinge-fac�ȭץ��@�I�I�A���`���רS���v�T
!�~��hinge(i)<<rod  �t�׸򨤳t�����ӭn�����ǻ���temp k �٬��~��                              
		hinges_temp(k)%X_i= hinges(i)%X_i - fac*(hinges(i)%X_i-hinges(i-1)%X_i)&
                                             /sqrt(dot_product(hinges(i)%X_i-hinges(i-1)%X_i,hinges(i)%X_i-hinges(i-1)%X_i))
        hinges_temp(k)%is_stationary = hinges(i)%is_stationary !stationary �]���~��
		k = k + 1

        !adding to alocate all properties of hinges(i) to hinges_temp(i)
        !hinges_temp(k) = hinges(i)
        ! updating position to ensure not overlapping

hinges_temp(k)= hinges(i) !2018/08/04 �ץ����~ �]��hinges_temp(k) �]�O�� hinges(i) �������_���Ӫ�(k+1),�]���n�~��(i)���Ҧ���T,���O��m�n���ʤ@��
                          !2018/08/04 �ץ����~ ����{���|�F�o�@��,�]��hinges_temp(k) �u�~�ӤF�y�� X_i �Mis_stationary
                          
		hinges_temp(k)%X_i= hinges(i)%X_i + fac*(hinges(i+1)%X_i-hinges(i)%X_i)&
                                             /sqrt(dot_product( hinges(i+1)%X_i-hinges(i)%X_i, hinges(i+1)%X_i-hinges(i)%X_i))
        hinges_temp(k)%is_stationary = hinges(i)%is_stationary
		k = k + 1
	else
		hinges_temp(k)= hinges(i)  !2018/08/04  �ƻshinges(i)����T
		k = k + 1
	end if

end do

deallocate(hinges)
deallocate(fibers)

allocate(hinges(ubound(hinges_temp,1)))
allocate(fibers(ubound(fibers_temp,1)))

hinges=hinges_temp
fibers=fibers_temp

deallocate(hinges_temp)
deallocate(fibers_temp)

!PRINT *, "Monitor Hinges"
!do i=1, ubound(hinges,1)
!	print *, hinges(i)%X_i
!end do
!PRINT *, "Monitor Fibers"
!do i=1, ubound(fibers,1)
!	print *, fibers(i)%first_hinge, fibers(i)%nbr_hinges
!end do

do i=1, ubound(hinges,1)
	hinges(i)%alpha    =0
end do

end subroutine hinges_damage

end module m_HingesDamage
!====================================================================