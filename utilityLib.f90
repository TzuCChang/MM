!====================================================================
module m_UtilityLib  !2018/07/21 new

use m_DataStructures

implicit none
contains

FUNCTION outerProd(a, b)
  real(8), DIMENSION(3,3)           :: outerProd
  real(8), DIMENSION(3), INTENT(IN) :: a, b

  outerProd =spread(a, dim=2, ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
END FUNCTION outerProd


function cross(vec1, vec2)
real(8), dimension (3):: cross, vec1, vec2

 cross(1)=  vec1(2)*vec2(3)-vec1(3)*vec2(2)
 cross(2)=-(vec1(1)*vec2(3)-vec1(3)*vec2(1))
 cross(3)=  vec1(1)*vec2(2)-vec1(2)*vec2(1)
end function cross


function clamp( a, b, c)
real(8)              ::a, b, c, clamp

if ( a < b ) then
clamp = b
return
end if

if ( a > c ) then
clamp = c
return
end if

clamp = a
end function clamp

!====================================================================

function are_boxes_intersect(hinge_pt1, hinge_pt2,segment_pt1, segment_pt2, radius, margin)

logical               ::are_boxes_intersect
real(8), dimension(3) ::hinge_pt1, hinge_pt2,segment_pt1, segment_pt2 
type(bound_box)       ::a, b
real(8)               :: radius, margin
are_boxes_intersect=.false.
a.X=(hinge_pt1+hinge_pt2)/2
a.W=abs(hinge_pt1-hinge_pt2)+2*radius+margin

b.X=(segment_pt1+segment_pt2)/2
b.W=abs(segment_pt1-segment_pt2)+2*radius+margin

are_boxes_intersect= (abs(a.X(1) - b.X(1)) * 2 < (a.W(1) + b.W(1))) .and.&
                     (abs(a.X(2) - b.X(2)) * 2 < (a.W(2) + b.W(2))) .and.&
                     (abs(a.X(3) - b.X(3)) * 2 < (a.W(3) + b.W(3)))  

end function are_boxes_intersect

!====================================================================

subroutine dist_pt_seg(ra,&
                       rb,&
                       rb_end,&
		               pb,&
                       Gab,&
                       Gab_norm,&
                       Sba,&
                       coll_course)

real(8), dimension(3):: ra, rb, rb_end, pb, Gab
real(8)              :: Sba, Gab_norm, lb
logical              :: coll_course


pb=rb_end-rb

lb=sqrt(dot_product(pb,pb))

pb=pb/lb

Sba=dot_product(ra-rb, pb)

Gab=ra-rb-Sba*pb
Gab_norm=sqrt(dot_product(Gab,Gab))


coll_course=.false.

if (Sba.le.lb .and. Sba .ge. 0.0) then
	coll_course= .true.
end if

end subroutine dist_pt_seg

!====================================================================
                       
subroutine dist_pt_pt(ra,&
                      rb,&
                      Gab,&
                      Gab_norm )

real(8), dimension(3):: ra, rb, Gab
real(8)              :: Gab_norm

Gab=ra-rb
Gab_norm=sqrt(dot_product (Gab, Gab))

end subroutine dist_pt_pt
                      

subroutine dist_segs(ra,&
                     ra_end,&
                     rb,&
                     rb_end,&
                     pa,&
                     pb,&
                     Gab,&
                     Gab_norm,&
                     Sba,&
                     coll_course,&
                     Sab)

real(8), dimension(3):: ra, ra_end, rb, rb_end, pa, pb, Gab
real(8)              :: Sab, Sba, la, lb, Gab_norm
logical              :: coll_course

pa=ra_end-ra
pb=rb_end-rb

la=sqrt(dot_product(pa,pa))
lb=sqrt(dot_product(pb,pb))

pa=pa/la
pb=pb/lb

Sab=(dot_product(ra-rb, pb)*dot_product(pa,pb)-dot_product(ra-rb,pa))&
    /(1d0-dot_product(pa, pb)**2d0)

Sba=(dot_product(rb-ra, pa)*dot_product(pb,pa)-dot_product(rb-ra,pb))&
    /(1d0-dot_product(pb, pa)**2d0)

Gab=ra+Sab*pa-rb-Sba*pb
Gab_norm=sqrt(dot_product(Gab,Gab))

coll_course=.false.
if (Sab.le.la .and. Sab .ge. 0.0) then
	if (Sba.le.lb .and. Sba .ge. 0.0) then
		coll_course= .true.
	end if
end if

end subroutine dist_segs
                     
!====================================================================                    

subroutine dist_segments( p1,&
                          q1,&
                          p2,&
                          q2,&
                          s,&
                          t,&
                          Gab,&
                          Gab_norm )

real(8), dimension(3):: p1, q1, p2, q2, d1, d2, Gab, r, c1, c2
real(8)              ::Gab_norm, epsilon, s, t
real(8)              ::a, e, f, c, b, denom
logical              :: coll_course

epsilon  =1D-40

d1=q1-p1
d2=q2-p2

r = p1 - p2

a = dot_product(d1,d1)
e = dot_product(d2,d2)
f = dot_product(d2,r)

if ((a<=epsilon) .and. (e<=epsilon)  ) then
    s = 0.0D0
    t = 0.0D0
    
    c1 = p1
    c2 = p2
    
    Gab_norm = sqrt(dot_product(c1 - c2, c1 - c2))
    Gab =  (c1 - c2)
    return
end if
    
if ( a<= epsilon ) then
    s = 0D0
    t = f/e
    t = clamp(t, 0.0D0, 1.0D0)    
else
    c = dot_product(d1,r)
    
    if(e <= epsilon)then
        t= 0.0D0
        s = clamp(-c /a, 0.0D0, 1.0D0)
    else
        b = dot_product(d1, d2)
        denom = a*e-b*b
        
        if(denom /=	0.0D0 ) then
            s = clamp((b*f -c*e) / denom  , 0.0D0, 1.0D0) 
        else
         s =0.0D0
        end if
        
        t =(b*s +f) / e
        
        if (t < 0D0) then
            t = 0.0D0
            s =clamp(-c / a, 0.0D0 , 1.0D0)
        else if ( t> 1.0) then
            t= 1.0D0
            s =clamp((b-c)/a, 0.0D0, 1.0D0)
       end if
   end if
    
end if

    c1 = p1 + d1*s
    c2 = p2 + d2*t
    
   !print *, "c1 " , c1
   !print *, "c2 " , c2
   !print *, "s " , s
   !print *, "t " , t 
   !print *, "p1 " , p1 
   !print *, "p2 " , p2 
   !print *, "d1 " , d1 
   !print *, "d2 " , d2 
   !print *, " "
   
    Gab_norm = sqrt(dot_product(c1 - c2, c1 - c2))
    Gab =  (c1 - c2)

end subroutine dist_segments
       
!====================================================================
                     
subroutine find_curvature(p1, p2, p3, curv)
implicit none 
real(8)                :: a,b,c,l1,l2,l2nb,n,rad,dotp,curv,scale1,scale2,t
real(8), dimension(3)  :: center, p1,p2,p3,v1,v2,v1n,v2n,v2nb
real(8), dimension(2)  :: p3_2d
integer                :: i

center = 0
rad    = 0
v1n    = 0
v2nb   = 0
n = size(p1,1)
v1 = p2 - p1 
v2 = p3 - p1
l1 = sqrt((v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)))
l2 = sqrt((v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3)))
v1n = v1
do i=1,3 
    v1n(i) = v1n(i)/l1
end do
v2n = v2
do i=1,3 
   v2n(i) = v2n(i)/l2
end do
dotp = v2n(1)*v1n(1) + v2n(2)*v1n(2) + v2n(3)*v1n(3)
v2nb = v2n
do i=1,3 
    v2nb(i) = v2nb(i) - dotp*v1n(i)
end do
l2nb = sqrt((v2nb(1)*v2nb(1)+v2nb(2)*v2nb(2)+v2nb(3)*v2nb(3)))
do i=1,3
    v2nb(i) = v2nb(i)/l2nb
end do
do i = 1,2
p3_2d(i) = 0
end do
do i = 1,3
    p3_2d(1) = p3_2d(1) + v2(i)*v1n(i)
    p3_2d(2) = p3_2d(2) + v2(i)*v2nb(i)
end do
a = l1
b = p3_2d(1)
c = p3_2d(2)
t = 0.5*(a-b)/c
scale1 = b/2 + c*t
scale2 = c/2 - b*t
do i = 1,3
center(i) = 0
end do 
do i=1,3
    center(i) = p1(i) + scale1*v1n(i) + scale2*v2nb(i)
end do
rad = sqrt((center(1)-p1(1))**2+(center(2)-p1(2))**2+(center(3)-p1(3))**2)
curv = rad

end subroutine find_curvature

!====================================================================

                      
end module m_UtilityLib

!====================================================================