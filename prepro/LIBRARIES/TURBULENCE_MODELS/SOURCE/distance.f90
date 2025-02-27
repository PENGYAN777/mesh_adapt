!============================================================ 
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!   Copyright: 1998-2003 Alberto Guardone
!              See COPYING file for copyright notice
!
!============================================================ 

MODULE distances

   USE lin_algebra


 CONTAINS


   FUNCTION distance_from_triangle(P, V0, V1, V2,  a_) RESULT(d)
  
      ! Note: the coefficients of equation 
      !   ax + by + cz + d = 0
      ! defining the plane are given by
      ! a = n(1), b = n(2), c = n(3), d = -SUM(n*V0)

      IMPLICIT NONE
      
      REAL(KIND=8), DIMENSION(3), INTENT(IN)  ::  P
      REAL(KIND=8), DIMENSION(3), INTENT(IN)  ::  V0, V1, V2
      REAL(KIND=8), DIMENSION(3), OPTIONAL    ::  a_
      
      REAL(KIND=8)  ::  d
      
      REAL(KIND=8), DIMENSION(2,3)    ::  dn 
      REAL(KIND=8), DIMENSION(3)      ::  I  
      REAL(KIND=8), DIMENSION(3)      ::  n,  u, v, w
      REAL(KIND=8), DIMENSION(3)      ::  dd
      REAL(KIND=8), DIMENSION(0:2)    ::  a
      REAL(KIND=8), DIMENSION(0:2,3)  ::  as
      REAL(KIND=8), DIMENSION(0:2)    ::  as_
      REAL(KIND=8)  ::  uu, uv, vv, wu, wv, den
      REAL(KIND=8)  ::  s, t
      INTEGER, DIMENSION(1)  ::  is_minv
      INTEGER  ::  is_min
            
      u = V1 - V0;  v = V2 - V0
      
      dn(1,:) = u;  dn(2,:) = v 
      n = vec_prod(dn)        ! Normal to plane Pi
      n = n / SQRT(SUM(n**2)) ! --> to unit normal
      

      ! Intersection point I(x,y,z)
      I = P - SUM( n*(P-V0) ) * n
      
      ! Parametric coordinates (s, t) of the intersection point I
      w = I - V0
      
      uu = SUM( u*u );  vv = SUM( v*v )
      uv = SUM( u*v );  wu = SUM( w*u );   wv = SUM( w*v )
      
      den = uv**2  -  uu*vv
      
      s = (uv*wv - vv*wu) / den;  t = (uv*wu - uu*wv) / den
      
      ! Areal baricentric coordinates
      a(0) = 1 - s - t;  a(1) = s;  a(2) = t
      
      is_min = 0; dd = 0
      IF (ANY(a < 0)) THEN ! Point I lies outside the triangle
         
         as = 0.d0
         
         dd(1) = distance_from_segment(P, V0, V1,  as_ )
         as(0,1) = as_(0);  as(1,1) = as_(1);  as(2,1) = 0.d0

         dd(2) = distance_from_segment(P, V0, V2,  as_ )
         as(0,2) = as_(0);  as(2,2) = as_(1);  as(1,2) = 0.d0

         dd(3) = distance_from_segment(P, V1, V2,  as_ )
         as(1,3) = as_(0);  as(2,3) = as_(1);  as(0,3) = 0.d0

         
         is_minv = MINLOC(dd); is_min = is_minv(1)
         
         d = dd(  is_min)
         a = as(:,is_min)
                                   
         ELSE ! Point I is inside the triangle
      
         d = SQRT( SUM(P - I)**2 )
      
      ENDIF
      
      IF ( PRESENT(a_) ) a_ = a
        
   END FUNCTION distance_from_triangle
  
  
   FUNCTION distance_from_segment(P, V0, V1,  a_) RESULT(d)
  
      IMPLICIT NONE
      
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  P
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  V0, V1
      REAL(KIND=8), DIMENSION(2), INTENT(OUT), OPTIONAL    ::  a_
      
      REAL(KIND=8)  ::  d
      
      REAL(KIND=8), DIMENSION(SIZE(P))  ::  I
      REAL(KIND=8), DIMENSION(SIZE(P))  ::  v, w
      REAL(KIND=8)  ::  c1, c2, s
      
 
      v = V1 - V0;  w = P - V0
      
      c1 = SUM(w*v)
      IF (c1 <= 0) THEN
         d = distance_from_point(P, V0)
         IF (PRESENT(a_)) THEN
            a_(1) = 1;  a_(2) = 0  
         ENDIF      
         RETURN
      ENDIF
      c2 = SUM(v*v)
      IF (c2 <= c1) THEN
         d = distance_from_point(P, V1)
         IF (PRESENT(a_)) THEN
            a_(1) = 0;  a_(2) = 1  
         ENDIF      
         RETURN
      ENDIF
      
      ! Parametric coordinates s of the intersection point I
      s = c1 / c2
      I = V0  +  s*v
      
      d = distance_from_point(P, I)
       
      ! Areal baricentric coordinates
      IF (PRESENT(a_)) THEN
         a_(1) = 1-s;  a_(2) = s  
      ENDIF      
             
   END FUNCTION distance_from_segment
  
    
   FUNCTION distance_from_point(P, P0) RESULT(d)
  
      IMPLICIT NONE
      
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  P
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  P0
      
      REAL(KIND=8)  ::  d
      
      d = SQRT(SUM((P-P0)**2))
             
   END FUNCTION distance_from_point
  

END MODULE distances
