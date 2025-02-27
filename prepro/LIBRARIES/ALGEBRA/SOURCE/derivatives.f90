!============================================================ 
!
!      Module: derivatives
!
! Description: Numerical evaluation of [partial] derivatives
!              of given [vector] functions 
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

MODULE derivatives


!   REAL(KIND=8), PARAMETER  ::  eps_M = 1.e-8
   REAL(KIND=8), PARAMETER  ::  eps_M = 1.e-14
!   REAL(KIND=8), PARAMETER  ::  eps_M = 1.e-4

   PRIVATE  ::  eps_M

 !============================================================ 
 CONTAINS
 !============================================================ 



   !============================================================ 
   FUNCTION df_dx (f, x, xref) RESULT (df)
   !============================================================ 

      IMPLICIT NONE

      REAL(KIND=8), EXTERNAL    ::  f
      REAL(KIND=8), INTENT(IN)  ::  x, xref
      REAL(KIND=8)              ::  df

      REAL(KIND=8)              ::  h
      
      
      h = SQRT(eps_M) * MAX( ABS(x), xref ) * SIGN(1.d0, x)
      
      df =  ( f(x+h) -  f(x-h) ) / (2*h)
     
      
   END FUNCTION df_dx 
   !============================================================ 



   !============================================================ 
   FUNCTION df_dx_pp (f, x, xref, param) RESULT (df)
   !============================================================ 

      IMPLICIT NONE

      REAL(KIND=8), EXTERNAL                  ::  f
      REAL(KIND=8),               INTENT(IN)  ::  x, xref
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  param
      REAL(KIND=8)                            ::  df

      REAL(KIND=8)              ::  h
      
      
      h = SQRT(eps_M) * MAX( ABS(x), xref ) * SIGN(1.d0, x)
      
      df =  ( f(x+h, param) -  f(x-h, param) ) / (2*h)
     
      
   END FUNCTION df_dx_pp 
   !============================================================ 


   !============================================================ 
   FUNCTION d2f_dx2 (f, x, xref) RESULT (d2f)
   !============================================================ 

      IMPLICIT NONE

      REAL(KIND=8), EXTERNAL    ::  f
      REAL(KIND=8), INTENT(IN)  ::  x, xref
      REAL(KIND=8)              ::  d2f

      REAL(KIND=8)              ::  h
      
      
      h = SQRT(SQRT(eps_M)) * MAX( ABS(x), xref ) * SIGN(1.d0, x)
           
      d2f =  ( f(x-h) - 2*f(x) + f(x+h) ) / h**2
     
      
   END FUNCTION d2f_dx2 
   !============================================================ 


   !============================================================ 
   FUNCTION df_dxx (f, xx, xxref) RESULT (df)
   !============================================================ 

      IMPLICIT NONE

      REAL(KIND=8), EXTERNAL                  ::  f
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  xx, xxref
      REAL(KIND=8), DIMENSION(SIZE(xx))       ::  df

      REAL(KIND=8), DIMENSION(SIZE(xx))       ::  hh
      INTEGER  :: i
      
      
      DO i = 1, SIZE(xx)
         
         hh = 0.d0
         hh(i) = SQRT(eps_M) * MAX( ABS(xx(i)), xxref(i) ) &
              * SIGN(1.d0, xx(i))
         
         df(i) =  ( f(xx+hh) -  f(xx-hh) ) / (2*hh(i))         

      ENDDO
      
      
   END FUNCTION df_dxx 
   !============================================================ 


   !============================================================ 
   FUNCTION df_dxx_pp (f, xx, xxref, param) RESULT (df)
   !============================================================ 

      IMPLICIT NONE

      REAL(KIND=8), EXTERNAL                  ::  f
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  xx, xxref
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  param
      REAL(KIND=8), DIMENSION(SIZE(xx))       ::  df

      REAL(KIND=8), DIMENSION(SIZE(xx))       ::  hh
      INTEGER  :: i
      
      
      DO i = 1, SIZE(xx)
         
         hh = 0.d0
         hh(i) = SQRT(eps_M) * MAX( ABS(xx(i)), xxref(i) ) &
              * SIGN(1.d0, xx(i))
         
         df(i) =  ( f(xx+hh, param) -  f(xx-hh, param) ) / (2*hh(i))
         
      ENDDO
      
      
   END FUNCTION df_dxx_pp 
   !============================================================ 


   !============================================================ 
   FUNCTION d2f_dxx2 (f, xx, xxref) RESULT (d2f)
   !============================================================ 

      IMPLICIT NONE

      REAL(KIND=8), EXTERNAL                  ::  f
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  xx, xxref
      REAL(KIND=8), DIMENSION(SIZE(xx),SIZE(xx))       ::  d2f

      REAL(KIND=8), DIMENSION(SIZE(xx))  ::  hh
      REAL(KIND=8), DIMENSION(SIZE(xx))  ::  hpp, hmp, hpm, hmm
      REAL(KIND=8)  ::  di, dj
      INTEGER  :: i, j
   
   
      DO i = 1, SIZE(xx)

         di = SQRT(SQRT(eps_M))* MAX( ABS(xx(i)), xxref(i) ) &
              * SIGN(1.d0, xx(i))

         hh = 0.d0; hh(i) = di
           
         d2f(i,i) =  ( f(xx-hh) - 2*f(xx) + f(xx+hh) ) / (hh(i)**2)


         DO j = i+1, SIZE(xx)

            dj = SQRT(SQRT(eps_M))* MAX( ABS(xx(j)), xxref(j) ) &
                 * SIGN(1.d0, xx(j))

            hpp = 0.d0; hpp(i) =   di; hpp(j) =   dj
            hpm = 0.d0; hpm(i) =   di; hpm(j) = - dj
            hmp = 0.d0; hmp(i) = - di; hmp(j) =   dj
            hmm = 0.d0; hmm(i) = - di; hmm(j) = - dj

            d2f(i,j) =  ((f(xx+hpp) - f(xx+hmp))/(2*di) &
                        -(f(xx+hpm) - f(xx+hmm))/(2*di) )/ (2*dj)

         ENDDO
      
         DO j = 1, i-1
         
           d2f(i,j) = d2f(j, i)
         
         ENDDO
      
      ENDDO
   
           
   END FUNCTION d2f_dxx2 
   !============================================================ 
   
  
END MODULE derivatives
