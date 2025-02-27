!============================================================ 
!
!      Module: fem_norms
!
! Description: procedure for the evaluation of the norm in 
!              the finite element method framework 
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

MODULE fem_norms



   !============================================================ 
   USE fem_ele_types
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 

!============================================================ 
!************************************************************ 
!============================================================ 
!
!  ================================
!  NORM OF AN INTERPOLATED QUANTITY 
!  ================================
!
!  
!  || U ||_L1  = SUM_i int [ |N_i U_i| dx ]
!
!  || U ||_L2  = SQRT ( SUM_i int [ |N_i U_i|^2 dx ] )
!
!  || U ||_Loo = MAX_i |N_i U_i| 
!
!
!============================================================ 
!************************************************************ 
!============================================================ 


   !============================================================ 
   FUNCTION fem_L1_norm(ele_d, uu)  RESULT (norm)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8)                                      ::  norm
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      INTEGER       ::  m, l
      !------------------------------------------------------------ 
      

      norm = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            norm  =  norm  +  ( ABS(SUM(uu(nn) * ww(:,l))) ) * rjp(l)
            
         ENDDO
      ENDDO


   END FUNCTION fem_L1_norm
   !============================================================ 



   !============================================================ 
   FUNCTION fem_L2_norm(ele_d, uu)  RESULT (norm)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8)                                      ::  norm
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      INTEGER       ::  m, l
      !------------------------------------------------------------ 
      

      norm = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            norm  =  norm  +  ( SUM(uu(nn) * ww(:,l))**2 ) * rjp(l)
            
         ENDDO
      ENDDO

      norm = SQRT(norm)


   END FUNCTION fem_L2_norm
   !============================================================ 



   !============================================================ 
   FUNCTION fem_Loo_norm(ele_d, uu)  RESULT (norm)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8)                                      ::  norm
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      !------------------------------------------------------------ 
      INTEGER       ::  m, l
      !------------------------------------------------------------ 
      

      norm = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww

         DO l = 1, ele_d(m) % l_G

            norm  =  MAX(norm, MAXVAL( ABS(uu(nn) * ww(:,l)) ) ) 
            
         ENDDO
         
      ENDDO


   END FUNCTION fem_Loo_norm
   !============================================================ 



!============================================================ 
!************************************************************ 
!============================================================ 
!
!  ===============================================
!  NORM OF A SCALAR FUNCTION KNOWN AT GAUSS POINTS 
!  ===============================================
!
!  
!  || u ||_L1  = int [ |u| dx ]
!
!  || u ||_L2  = SQRT ( int [ |u|^2 dx ] )
!
!  || u ||_Loo = MAX |u| 
!
!
!============================================================ 
!************************************************************ 
!============================================================ 



   !============================================================ 
   FUNCTION fem_L1_norm_fs(ele_d, ele_fs)  RESULT (norm)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:),   INTENT(IN)  ::  ele_d
      TYPE(scalar_func_gauss), DIMENSION(:),   INTENT(IN)  ::  ele_fs

      REAL(KIND=8)                                      ::  norm
      !------------------------------------------------------------ 
      ! Elemental quantities
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      INTEGER       ::  m, l
      !------------------------------------------------------------ 
      

      norm = 0.d0

      DO m = 1, SIZE(ele_d)

         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            norm  =  norm  +  ABS(ele_fs(m)%ss(l)) * rjp(l)
            
         ENDDO

      ENDDO


   END FUNCTION fem_L1_norm_fs
   !============================================================ 


   !============================================================ 
   FUNCTION fem_L2_norm_fs(ele_d, ele_fs)  RESULT (norm)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:),   INTENT(IN)  ::  ele_d
      TYPE(scalar_func_gauss), DIMENSION(:),   INTENT(IN)  ::  ele_fs

      REAL(KIND=8)                                      ::  norm
      !------------------------------------------------------------ 
      ! Elemental quantities
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      INTEGER       ::  m, l
      !------------------------------------------------------------ 
      

      norm = 0.d0

      DO m = 1, SIZE(ele_d)

         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            norm  =  norm  +  ( (ele_fs(m)%ss(l))**2 ) * rjp(l)
            
         ENDDO

      ENDDO

      norm = SQRT(norm)


   END FUNCTION fem_L2_norm_fs
   !============================================================ 


   !============================================================ 
   FUNCTION fem_Loo_norm_fs(ele_d, ele_fs)  RESULT (norm)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:), INTENT(IN)  ::  ele_d
      TYPE(scalar_func_gauss), DIMENSION(:), INTENT(IN)  ::  ele_fs

      REAL(KIND=8)                                       ::  norm
      !------------------------------------------------------------ 
      INTEGER       ::  m, l
      !------------------------------------------------------------ 
      

      norm = 0.d0

      DO m = 1, SIZE(ele_d)

         DO l = 1, ele_d(m) % l_G

            norm  =  MAX(norm, ABS(ele_fs(m)%ss(l)) ) 
            
         ENDDO
         
      ENDDO


   END FUNCTION fem_Loo_norm_fs
   !============================================================ 



END MODULE fem_norms
