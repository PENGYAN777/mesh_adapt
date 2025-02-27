!============================================================ 
!
!      Module: non_lin_algebra
!
! Description: Procedures to solve (simple) non linear 
!              equations  
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

MODULE non_lin_algebra



 !============================================================ 
 CONTAINS
 !============================================================ 


   !============================================================ 
   FUNCTION solve_eq_3(p,q,r) RESULT( xx )
   !============================================================ 


      ! Solve an algebraic equation of third order in the form
      !
      !       x^3 + p*x^2 + q*x + r = 0
      !


      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)       :: p, q, r
      COMPLEX(KIND=8), DIMENSION(3)  :: xx

      REAL(KIND=8)  ::  C1, C2, C3, CC1, CC2 , det, phi
      REAL(KIND=8), PARAMETER  ::  pi = 3.14159265358979323846d0 


      C1 = (3.d0*q - p*p)/3.d0 
      C2 = (2.d0*p*p*p - 9.d0*p*q + 27.d0*r)/27.d0

      det = (C2*C2)/4.d0 + (C1*C1*C1)/27

      IF ( det >= 0. ) THEN

         det = SQRT(det)
      
         CC1 = -0.5d0*C2 + det  
         CC2 = -0.5d0*C2 - det  

         CC1 = SIGN(1.d0,CC1) * ABS(CC1)**(1.d0/3.d0)  
         CC2 = SIGN(1.d0,CC2) * ABS(CC2)**(1.d0/3.d0)

         xx(1) = (1.d0,0.d0)* (CC1+CC2)
         xx(2) = (1.d0,0.d0)* (-0.5d0*(CC1+CC2))  &
            + (0.d0,1.d0)* (0.5d0*(CC1-CC2)*SQRT(3.d0))  
         xx(3) = (1.d0,0.d0)* (-0.5d0*(CC1+CC2))  &
            + (0.d0,1.d0)* (-0.5d0*(CC1-CC2)*SQRT(3.d0))  

      ELSE 

         C3 = 2.d0 * SQRT(-C1/3.d0)
         ! phi =  (C2*C2/4.d0) / (-C1*C1*C1/27.d0) 
         ! phi =  SQRT(phi)
         phi =  (3.d0*C2)/(C1*C3)
         phi =  ACOS(phi) 
         !xx(1) = SIGN(-2.d0,C1) * SQRT( -C1/3.d0) * COS(phi/3.d0)
         !xx(2) = SIGN(-2.d0,C1) * SQRT( -C1/3.d0) * COS(phi/3.d0 &
         !      + (2.d0/3.d0)*pi)
         !xx(3) = SIGN(-2.d0,C1) * SQRT( -C1/3.d0) * COS(phi/3.d0 &
         !      + (4.d0/3.d0)*pi)
         xx(1) = 2.d0 * SQRT( -C1/3.d0) * COS(phi/3.d0)
         xx(2) = 2.d0 * SQRT( -C1/3.d0) * COS(phi/3.d0 + (2.d0/3.d0)*pi)
         xx(3) = 2.d0 * SQRT( -C1/3.d0) * COS(phi/3.d0 + (4.d0/3.d0)*pi)

      ENDIF 

      xx(1) = xx(1) - p/3.d0
      xx(2) = xx(2) - p/3.d0
      xx(3) = xx(3) - p/3.d0


   END FUNCTION solve_eq_3
   !============================================================ 



END MODULE non_lin_algebra
