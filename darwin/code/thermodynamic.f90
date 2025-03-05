MODULE thermodynamic

  IMPLICIT NONE
  
  REAL(KIND=8), PUBLIC :: mu_ref, T_ref, T_char, T_scaling, L_scaling, P_scaling


  CONTAINS


   FUNCTION c__ww(p, ww)   RESULT(c)
   !===================================================================
   IMPLICIT NONE
   INTEGER,                     INTENT(IN) ::  p
   REAL(KIND=8), DIMENSION(:),  INTENT(IN) ::  ww

   REAL(KIND=8)                            :: c  
   REAL(KIND=8)                            :: Pi, rho, re
   REAL(KIND=8), PARAMETER                 :: gamma = 1.4  
   !===================================================================
   
   rho = ww(1)
   re  = ww(p) + ww(1)/(gamma-1) - 0.5*SUM((ww(2:p-1))**2)/ww(1)
      
   Pi  = (gamma - 1.d0)*re

   c = SQRT(gamma*Pi/rho)

   END FUNCTION c__ww   




   
   FUNCTION Mach__ww(p, ww)   RESULT(Mach)
   !===================================================================
   IMPLICIT NONE

   INTEGER,                     INTENT(IN) ::  p
   REAL(KIND=8), DIMENSION(:),  INTENT(IN) ::  ww

   REAL(KIND=8)                            ::  Mach 
   REAL(KIND=8)                            ::  mod_u, c
   !===================================================================
     
   mod_u = SQRT( SUM(ww(2:p-1)**2) ) / ww(1)

   c = c__ww(p, ww)

   Mach = mod_u / c

   END FUNCTION Mach__ww

   
   
   

   FUNCTION Pressure__ww(p, ww)   RESULT(Pressure)
   !===================================================================
   IMPLICIT NONE

   INTEGER,                     INTENT(IN) ::  p
   REAL(KIND=8), DIMENSION(:),  INTENT(IN) ::  ww

   REAL(KIND=8)                            ::  Pressure, rho, re 
   REAL(KIND=8)                            ::  gamma = 1.4
   !===================================================================

   rho = ww(1)
   re  = ww(p) + ww(1)/(gamma-1) - 0.5*SUM((ww(2:p-1))**2)/ww(1)

   Pressure  = (gamma - 1.d0)*re  
   
   END FUNCTION Pressure__ww




   FUNCTION Temperature__ww(p, ww)   RESULT(Temperature)
   !===================================================================
   IMPLICIT NONE

   INTEGER,                     INTENT(IN) ::  p
   REAL(KIND=8), DIMENSION(:),  INTENT(IN) ::  ww

   REAL(KIND=8) ::  Temperature, rho, re, Pi, M, c
   REAL(KIND=8) ::  gamma = 1.4
   !===================================================================

   rho = ww(1)
   re  = ww(p) + ww(1)/(gamma-1) - 0.5*SUM((ww(2:p-1))**2)/ww(1)

   Pi = (gamma - 1.d0)*re  
   
   c = SQRT(gamma*Pi/rho)
   
   M = SQRT(SUM(ww(2:3)**2)) / c
   
   Temperature = SUM(ww(2:3)**2) / (gamma*M**2)
   
   END FUNCTION Temperature__ww




   FUNCTION Entropy__ww(ww)   RESULT(Entropy)
   !==================================================================
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),  INTENT(IN) ::  ww
   REAL(KIND=8)                            ::  Entropy 
   !===================================================================
  
   Entropy = ww(1)
    
   END FUNCTION Entropy__ww





   FUNCTION  Knudsen__ww(grid, ww)   RESULT(Kn)
   !---------------------------------------------------------------------------
   USE derivatives
   
   IMPLICIT NONE
   
   TYPE(grid_type),     INTENT(IN) :: grid
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
   
   REAL(KIND=8), DIMENSION(SIZE(ww,2)) :: Kn
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1)-2,SIZE(ww,2)) :: gradD
   
   REAL(KIND=8), DIMENSION(SIZE(ww,2)) :: r, P, T, mu
   
   REAL(KIND=8) :: mFpath, PI=3.141592
   
   INTEGER :: j
   !---------------------------------------------------------------------------

   DO j = 1, SIZE(ww,2)
     r(j) = ww(1,j)
     P(j) = Pressure__ww(SIZE(ww,1), ww(:,j))
   ENDDO
   
   ! Sutherland law for viscosity
   T = (P / (r * 1.d0)) * T_scaling
   
   mu = mu_ref * (T/T_ref)**1.5  &
      *  (T_ref + T_char) / (T + T_char)
   
   mu = mu / (L_scaling*P_Scaling/T_scaling*SQRT(T_scaling))
   
   
   ! Computation of the gradient of the density velocity 
   ! modulus and temperature   
   CALL gradient (grid, r, gradD)

   DO j = 1, SIZE(ww,2)

     mFpath = 16.d0/(5*SQRT(PI)) * mu(j) / (r(j)*SQRT(2*P(j)/r(j)))
   
     Kn(j) = mFpath * SQRT(SUM(gradD(:,j)**2)) / r(j)

   ENDDO

   
   END FUNCTION  Knudsen__ww


END MODULE thermodynamic
