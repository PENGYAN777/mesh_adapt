!============================================================ 
!
!      Module: roe_average
!
! Description: Procedures for the Roe scheme: evaluation of
!              the Roe-averaged state and entropy fix of
!              the Roe-averaged eigenvalues 
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

MODULE roe_average


   !============================================================ 
   USE euler_equations
   USE thermodynamics
   !============================================================ 

   INTEGER  ::  ef_type
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE  ::  ef_thresold

   INTEGER, PARAMETER  ::   EF_HH1    = 1, &
                            EF_HH2    = 2, &
                            EF_HH2_VS = 3, &
                            EF_CONST_DISS = 4
   CHARACTER(*), DIMENSION(1:4), PARAMETER  ::  ef_names = &
                        (/ 'HARTEN AND HYMAN 1              ', &
                           'HARTEN AND HYMAN 2              ', & 
                           'HARTEN AND HYMAN 2, MODIFIED VS ', & 
                           'CONSTANT DISSIPATION COEFFICIENT' /)
   INTEGER  ::  average_type
   INTEGER, PARAMETER  ::   RHO_MEAN        = 1, &
                            RHO_PWG         = 2, & ! Not implemented
                            RHO_E_r         = 3, &
                            RHO_E_r__T_NUM  = 4, & ! Not implemented
                            RHO_E_r__T      = 5, & ! Not implemented
                            RHO_T_v         = 6, &
                            VINOKUR_AVERAGE = 7
                            
   CHARACTER(*), DIMENSION(1:7), PARAMETER  ::  average_names = &
                        (/ 'ARITHMETIC MEAN                            ', &
                           'EXACT SOLUTION FOR PWG                     ', & 
                           'UNKNWONS E AND RHO                         ', & 
                           'UNKNWONS E AND RHO, NUM. DERIV, T AUXILIARY', &
                           'UNKNWONS E AND RHO, T AUXILIARY            ', &
                           'UNKNWONS T AND V                           ', &
                           'VINOKUR & MONTAGNE                         ' /)

 !============================================================ 
 CONTAINS
 !============================================================ 

   !============================================================ 
   SUBROUTINE  read_param_roe_average(idf)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE

      INTEGER  ::  idf
      !------------------------------------------------------------
      INTEGER  ::  ww_size
      !------------------------------------------------------------


         READ (idf,*) ef_type
         READ (idf,*) ww_size
         ALLOCATE (ef_thresold(ww_size))
         READ (idf,*) ef_thresold
         READ (idf,*) average_type


   END SUBROUTINE  read_param_roe_average
   !============================================================ 




   !============================================================ 
   SUBROUTINE init_roe_average
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------



   END SUBROUTINE init_roe_average
   !============================================================ 




   !============================================================ 
   SUBROUTINE write_param_roe_average(idf)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE

      INTEGER  ::  idf
      !------------------------------------------------------------


      WRITE(idf,*) '   PARAMETERS FOR ROE AVERAGE'
      WRITE(idf,*) '   Entropy fix type: ', ef_names(ef_type)
      WRITE(idf,*) '   Thresolds for entropy fix: ', ef_thresold
      WRITE(idf,*) '   Solution of the supplementary equation: ', average_names(average_type)
      WRITE(idf,*)
      

   END SUBROUTINE write_param_roe_average
   !============================================================ 




   !============================================================ 
   FUNCTION entropy_fix(lambda, ww, nn) RESULT(lambda_ef)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  lambda
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  ww, nn
      REAL(KIND=8), DIMENSION(SIZE(lambda))   ::  lambda_ef
      REAL(KIND=8), DIMENSION(SIZE(lambda))   ::  delta

      REAL(KIND=8)  ::  Mach_n
      REAL(KIND=8), PARAMETER  ::  eps = 1.e-16


      SELECT CASE( ef_type )
      
         CASE( EF_HH1)
      
            delta = ef_thresold
            
            lambda_ef = MAX(ABS(lambda), delta)
      
      
         CASE( EF_HH2)

            delta = ef_thresold
           
            lambda_ef = (0.5 + SIGN(0.5d0, lambda - delta))    &
                        * lambda                               &
                      + (0.5 - SIGN(0.5d0, lambda - delta))    &
                        * (lambda**2 + delta**2) / (2* delta + eps)
      
         CASE( EF_HH2_VS)

            Mach_n = Mach_n__ww(ww, nn) 
            delta = ef_thresold * (ABS(Mach_n) + 1) 

            lambda_ef = (0.5 + SIGN(0.5d0, lambda - delta))    &
                        * lambda                               &
                      + (0.5 - SIGN(0.5d0, lambda - delta))    &
                        * (lambda**2 + delta**2) / (2* delta + eps)

         CASE( EF_CONST_DISS)

            lambda_ef = ABS(ef_thresold)

        CASE DEFAULT
         
           WRITE(*,*) ' Entropy fix of unknown type. STOP'
           STOP
           
      END SELECT


   END FUNCTION entropy_fix
   !============================================================ 



   !============================================================ 
   FUNCTION entropy_fix_eos(lambda, ww, eos, nn) RESULT(lambda_ef)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  lambda
      TYPE(eos_ext_type),         INTENT(IN)  ::  eos
      REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  ww, nn
      
      REAL(KIND=8), DIMENSION(SIZE(lambda))   ::  lambda_ef
     
      REAL(KIND=8), DIMENSION(SIZE(lambda))   ::  delta

      REAL(KIND=8)  ::  Mach_n
      REAL(KIND=8), PARAMETER  ::  eps = 1.e-16


      SELECT CASE( ef_type )
      
         CASE( EF_HH1)
      
            delta = ef_thresold
            
            lambda_ef = MAX(ABS(lambda), delta)
      
      
         CASE( EF_HH2)

            delta = ef_thresold
           
            lambda_ef = (0.5 + SIGN(0.5d0, lambda - delta))    &
                        * lambda                               &
                      + (0.5 - SIGN(0.5d0, lambda - delta))    &
                        * (lambda**2 + delta**2) / (2* delta + eps)
      
         CASE( EF_HH2_VS)

            Mach_n = Mach__ww_eos_nn(ww, eos, nn) 
            delta = ef_thresold * (ABS(Mach_n) + 1) 

            lambda_ef = (0.5 + SIGN(0.5d0, lambda - delta))    &
                        * lambda                               &
                      + (0.5 - SIGN(0.5d0, lambda - delta))    &
                        * (lambda**2 + delta**2) / (2* delta + eps)

         CASE( EF_CONST_DISS)

            lambda_ef = ABS(ef_thresold)

         CASE DEFAULT
         
           WRITE(*,*) ' Entropy fix of unknown type. STOP'
           STOP
           
      END SELECT


   END FUNCTION entropy_fix_eos
   !============================================================ 



   !============================================================ 
   SUBROUTINE intermediate_ww (ww_l, eos_l, ww_r, eos_r, &
                               ww, eos)
   !============================================================ 

      !------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN)   ::   ww_l,  ww_r
      TYPE(eos_type),             INTENT(IN)   ::  eos_l, eos_r

      REAL(KIND=8), DIMENSION(:), INTENT(OUT)  ::  ww
      TYPE(eos_ext_type),         INTENT(OUT)  ::  eos
      
      !------------------------------------------------------------
      REAL(KIND=8), DIMENSION(SIZE(ww_l))   ::  vv 
      REAL(KIND=8)  ::  s_l, s_r, P_l, P_r
      INTEGER  ::  p
      !------------------------------------------------------------

     
      !------------------------------------------------------------
      ! Roe averaging in the variable vv = (rho, u(:), ht)
      ! --------------------------------------------------
      p = SIZE(ww_l)
      
      P_l = eos_l%P;        P_r = eos_r%P
      
      IF (ww_l(1) < 0.d0  .OR.  ww_r(1) < 0.d0) THEN
        WRITE(*,*)			       
        WRITE(*,*) 'ERROR. INTERMEDIATE_WW:'   
        WRITE(*,*) 'Negative density obtained.'
        WRITE(*,*)			       
        STOP
      ENDIF				      
      
      s_l = SQRT(ww_l(1));  s_r = SQRT(ww_r(1))

      vv(2:p-1) = (ww_l(2:p-1)/s_l + ww_r(2:p-1)/s_r)            / (s_l + s_r)
      vv(p)     = ( (ww_l(p) + P_l)/s_l  + (ww_r(p) + P_r)/s_r ) / (s_l + s_r)  
      CALL rho_roe(ww_l, eos_l, ww_r, eos_r,  vv, eos)
      !------------------------------------------------------------


      !------------------------------------------------------------
      ! Intermediate state in conservative variables
      ! --------------------------------------------
      ww(1)   = vv(1)
      ww(2:p) = vv(1)*vv(2:p)
      ww(p)   = vv(1)*vv(p) - eos%P
      !------------------------------------------------------------


   END SUBROUTINE intermediate_ww
   !============================================================ 


   !============================================================ 
   SUBROUTINE  rho_roe(ww_l, eos_l, ww_r, eos_r,  vv, eos)
   !============================================================ 
  

      !------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::   ww_l,  ww_r
      TYPE(eos_type),             INTENT(IN)     ::  eos_l, eos_r
      
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  vv
      TYPE(eos_ext_type),         INTENT(OUT)    ::  eos
      
      REAL(KIND=8)                            ::  rho, T, h
      INTEGER  ::  p
      !------------------------------------------------------------
                                                               
      
      SELECT CASE(average_type)

         CASE(RHO_MEAN)
            p = SIZE(vv)
            rho = 0.5*(ww_l(1) + ww_r(1))
            vv(1) = rho
            h = vv(p) - 0.5*SUM(vv(2:p-1)**2)           
            T = T__h_v(h, 1/rho)
            eos = eos_ext__T_v(T, 1/rho)

         CASE(RHO_E_r)
            CALL rho_roe_newton_re_r(ww_l, eos_l, ww_r, eos_r,  vv, eos)
      
         CASE(RHO_T_v)
            CALL rho_roe_newton_T_v(ww_l, eos_l, ww_r, eos_r,  vv, eos)
            
         CASE(VINOKUR_AVERAGE)
            CALL vinokur_montagne_average(ww_l, eos_l, ww_r, eos_r,  vv, eos)
            
         CASE DEFAULT
            WRITE(*,*) ' Unknown kind of intermediate density'
            WRITE(*,*) ' STOP'
            STOP
            
      END SELECT


   END SUBROUTINE  rho_roe
   !============================================================ 


   
   !============================================================ 
   SUBROUTINE rho_roe_newton_re_r(ww_l, eos_l, ww_r, eos_r,  vv, eos)
   !============================================================ 


      !------------------------------------------------------------
      USE lin_algebra
      !------------------------------------------------------------


      !------------------------------------------------------------
      IMPLICIT NONE

      ! Left and right state (ww_l and ww_r) in conservative 
      ! variables, intermediate state in the vv = (rho, velocity, htot) 
      ! variables
      REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::   ww_l,  ww_r
      TYPE(eos_type),             INTENT(IN)     ::  eos_l, eos_r
      
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  vv
      TYPE(eos_ext_type),         INTENT(OUT)    ::  eos
      !------------------------------------------------------------
      REAL(KIND=8)  ::  rho_l, re_l, P_l 
      REAL(KIND=8)  ::  rho_r, re_r, P_r    
      REAL(KIND=8)  ::  Dr,    Dre,  DP  ! Jumps
      REAL(KIND=8)  ::  rho, re, h, P        ! Intermediate quantities
      REAL(KIND=8), DIMENSION(2)  ::  FF, DD
      REAL(KIND=8), DIMENSION(2, 2)  ::  JF
      REAL(KIND=8), DIMENSION(2)  ::  dP_d
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P_d2
      ! Parameters for Newton method
      REAL(KIND=8)  ::  rho_max
      REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-8
      INTEGER ,     PARAMETER  ::  NEWTON_MAX_ITER = 200
      INTEGER  ::  it, n
      !------------------------------------------------------------
      
      
      
      !------------------------------------------------------------
      ! Retrieve useful thermodynamic quantities at the left and 
      ! right state and compute jumps
      n = SIZE(ww_l)

      rho_l = ww_l(1);  rho_r = ww_r(1)               ! Density

      re_l = ww_l(n) - 0.5*SUM(ww_l(2:n-1)**2)/ww_l(1) ! Internal energy
      re_r = ww_r(n) - 0.5*SUM(ww_r(2:n-1)**2)/ww_r(1) ! (per unit volume)

      P_l = eos_l%P;  P_r = eos_r%P                ! Pressure

      Dr   = rho_r - rho_l;  Dre  =  re_r - re_l;  DP   =   P_r - P_l                           !   "  

      h = vv(n) - 0.5*SUM(vv(2:n-1)**2)           ! Specific enthalpy
      !------------------------------------------------------------

      ! Initial guess
      re  = ( re_l +  re_r) / 2
      rho = (rho_l + rho_r) / 2

      ! The supplementary equation may be an identity,
      ! thus leading to a singular problem, and this is 
      ! to be checked before starting ther Newton solver.  
      ! However, the supplementary equation is to be evaluated 
      ! using an internal energy and density satisfying the
      ! constraint on the enthalpy... since it may be not an
      ! identity but accidentally satisfyed by the initial guess.            

!      E_ = rho_*h_ - P__h_rho(h_, rho_)      
      DO it = 1, NEWTON_MAX_ITER
      
         vv(1) = rho
         eos = eos_ext__e_v(re/rho, 1/rho)
         P = eos%P
         dP_d = eos%dP
         
         FF(2) =  (re + P)/rho - h
                  
         IF ( ABS(FF(2)) < NEWTON_EPS ) EXIT

         JF(2,1) = (1 + dP_d(1)) / rho 
         
         re = re  -  FF(2) / JF(2,1)
      
      ENDDO

      rho_max = 1/minimum_specific_volume()

      DO it = 1, NEWTON_MAX_ITER
         
         vv(1) = rho
         eos = eos_ext__e_v(re/rho, 1/rho)
         P = eos%P
         dP_d = eos%dP
         
         d2P_d2 = d2P__dre2_dr2(re, rho)

         ! Evaluation of the system
         ! ------------------------
         
         ! First equation: supplementary equation
         FF(1) = dP_d(1)*Dre + dP_d(2)*Dr - DP
         ! Second equation: constraint on enthalpy
         FF(2) = (re + P)/rho - h

         IF ( MAXVAL(ABS(FF)) < NEWTON_EPS ) RETURN                  

         ! Jacobian of the system
         ! ----------------------
         
         ! First equation: supplementary equation
         JF(1,1) = d2P_d2(1,1)*Dre  +  d2P_d2(1,2)*Dr ! d2P_d2(1,2) = d2P_d2(2,1)
         JF(1,2) = d2P_d2(2,1)*Dre  +  d2P_d2(2,2)*Dr ! from thermodynamics
         ! Second equation: constraint on enthalpy
         JF(2,1) = (1 + dP_d(1)) / rho 
         JF(2,2) = (dP_d(2)*rho - re - P) / rho**2

         ! Computation of the new solution
         ! -------------------------------
         DD = - MATMUL( inverse(JF), FF )
                  
         re = re + DD(1);  rho = MIN(rho_max, rho + DD(2))
      
      ENDDO
      
      WRITE(*,*)
      WRITE(*,*) 'I am not able to find the appropriate'
      WRITE(*,*) 'intermediate density. '
      WRITE(*,*) 'rho_l:', rho_l
      WRITE(*,*) 'rho_r:', rho_r
      WRITE(*,*) 'rho found: ',rho
      WRITE(*,*) 'resid supp. eq.: ', FF(1)
      WRITE(*,*) 'resid e + P/rho = ht: ', FF(2)
      WRITE(*,*) 'h: ', h
      WRITE(*,*) 'u: ', vv(2:n-1)
      WRITE(*,*) 'What about STOPping here?'
      STOP
 

   END SUBROUTINE rho_roe_newton_re_r
   !============================================================ 
   
   
   
   !============================================================ 
   SUBROUTINE rho_roe_newton_T_v(ww_l, eos_l, ww_r, eos_r,  vv, eos)
   !============================================================ 

      !------------------------------------------------------------
      USE lin_algebra
      !------------------------------------------------------------


      !------------------------------------------------------------
      IMPLICIT NONE

      ! Left and right state (ww_l and ww_r) in conservative 
      ! variables, intermediate state in the vv = (rho, velocity, htot) 
      ! variables
      REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::   ww_l,  ww_r
      TYPE(eos_type),             INTENT(IN)     ::  eos_l, eos_r
      
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  vv
      TYPE(eos_ext_type),         INTENT(OUT)    ::  eos
      !------------------------------------------------------------
      REAL(KIND=8)  ::  rho_l, re_l, P_l, T_l
      REAL(KIND=8)  ::  rho_r, re_r, P_r, T_r    
      REAL(KIND=8)  :: Dr,    Dre,  DP  ! Jumps
      REAL(KIND=8)  :: rho, re, h, P, T , v       ! Intermediate quantities
      REAL(KIND=8), DIMENSION(2)  ::  FF, DD
      REAL(KIND=8), DIMENSION(2, 2)  ::  JF
      REAL(KIND=8), DIMENSION(2)  ::  dP_d, de_d
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P_d2
      ! Parameters for Newton method
      REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-8
      REAL(KIND=8)  :: v_min
      INTEGER ,     PARAMETER  ::  NEWTON_MAX_ITER = 1000
      INTEGER  ::  it, n
      !------------------------------------------------------------
 
      
      
      !------------------------------------------------------------
      ! Retrieve useful thermodynamic quantities at the left and 
      ! right state and compute jumps
      n = SIZE(ww_l)

      rho_l = ww_l(1);  rho_r = ww_r(1)            ! Density

      re_l = ww_l(n) - 0.5*SUM(ww_l(2:n-1)**2)/ww_l(1) ! Internal energy
      re_r = ww_r(n) - 0.5*SUM(ww_r(2:n-1)**2)/ww_r(1) ! (per unit volume)

      P_l = eos_l%P;  P_r = eos_r%P                ! Pressure
      T_l = eos_l%T;  T_r = eos_r%T                ! Temperature

      Dr   = rho_r - rho_l;  Dre  =  re_r - re_l;  DP   =   P_r - P_l                           !   "  

      h = vv(n) - 0.5*SUM(vv(2:n-1)**2)           ! Specific enthalpy
      !------------------------------------------------------------

      ! Initial guess
      re  = (  re_l +  re_r) / 2
      rho = ( rho_l + rho_r) / 2
      T   = (   T_l +   T_r) / 2
      v   = 1/rho
      ! The supplementary equation may be an identity,
      ! thus leading to a singular problem, and this is 
      ! to be checked before starting ther Newton solver.  
      ! However, the supplementary equation is to be evaluated 
      ! using an internal energy and density satisfying the
      ! constraint on the enthalpy... since it may be not an
      ! identity but accidentally satisfyed by the initial guess.            

!      re = rho*h - P(h, rho) 
      DO it = 1, NEWTON_MAX_ITER
      
         vv(1) = rho
         eos = eos_ext__T_v(T, 1/rho)
         re  = eos%e * rho
         P   = eos%P
       
         FF(2) =  (re + P)/rho - h
                  
         IF ( ABS(FF(2)) < NEWTON_EPS ) EXIT

         dP_d  =  dP__dT_dv(T, v)
         de_d  =  de__dT_dv(T, v)
         
         JF(2,1) = de_d(1) + dP_d(1)*v
         
         T = MAX(1.d-12, T  -  FF(2) / JF(2,1))
      
      ENDDO

      v_min = minimum_specific_volume()


      DO it = 1, NEWTON_MAX_ITER

         
         vv(1) = rho
         eos = eos_ext__T_v(T, v)
         re  = eos%e * rho
         P   = eos%P
         dP_d = eos%dP
               
         ! First equation: supplementary equation
         FF(1) = dP_d(1)*Dre + dP_d(2)*Dr - DP
         ! Second equation: constraint on enthalpy
         FF(2) = (re + P)/rho - h

         IF ( MAXVAL(ABS(FF)) < NEWTON_EPS ) RETURN

         ! Jacobian of the system
         ! ----------------------
         ! Evaluate first order derivatives of dP(E,r)
         ! with respect to T and v
         d2P_d2 = d2P__dre_dr__dT_dv(T, v)                 

         ! First equation: supplementary equation
         JF(1,1) = d2P_d2(1,1)*Dre  +  d2P_d2(2,1)*Dr 
         JF(1,2) = d2P_d2(1,2)*Dre  +  d2P_d2(2,2)*Dr 
         ! Second equation: constraint on enthalpy
         dP_d  =  dP__dT_dv(T, v)
         de_d  =  de__dT_dv(T, v)

         JF(2,1) = de_d(1) + dP_d(1)*v 
         JF(2,2) = de_d(2) + dP_d(2)*v  +  P 

         DD = - MATMUL( inverse(JF), FF )
                  
         T = MAX(1.d-12, T + DD(1));  v = MAX(v_min, v + DD(2))
         rho = 1/v

      ENDDO
            
      WRITE(*,*)
      WRITE(*,*) 'I am not able to find the appropriate'
      WRITE(*,*) 'intermediate density. '
      WRITE(*,*) 'rho_l:', rho_l
      WRITE(*,*) 'rho_r:', rho_r
      WRITE(*,*) 'rho found: ',rho
      WRITE(*,*) 'resid supp. eq.: ', FF(1)
      WRITE(*,*) 'resid e + P/rho = ht: ', FF(2)
      WRITE(*,*) 'h: ', h
      WRITE(*,*) 'u: ', vv(2:n-1)
      WRITE(*,*) 'What about STOPping here?'
      STOP
      
     
   END SUBROUTINE rho_roe_newton_T_v
   !============================================================ 

   
   
   !============================================================ 
   SUBROUTINE vinokur_montagne_average(ww_l, eos_l, ww_r, eos_r,  vv, eos)
   !============================================================ 


      !------------------------------------------------------------
      USE lin_algebra
      !------------------------------------------------------------


      !------------------------------------------------------------
      IMPLICIT NONE

      ! Left and right state (ww_l and ww_r) in conservative 
      ! variables, intermediate state in the vv = (rho, velocity, htot) 
      ! variables
      REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::   ww_l,  ww_r
      TYPE(eos_type),             INTENT(IN)     ::  eos_l, eos_r
      
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  vv
      TYPE(eos_ext_type),         INTENT(OUT)    ::  eos
      !------------------------------------------------------------
      REAL(KIND=8)  ::  rho_l, re_l, P_l, T_l, h_l
      REAL(KIND=8)  ::  rho_r, re_r, P_r, T_r, h_r    
      REAL(KIND=8)  :: Dr,    Dre,  DP  ! Jumps
      REAL(KIND=8)  :: rho, h, T , dP_dre, dP_dr, c2, c! Intermediate quantities
      REAL(KIND=8), DIMENSION(2)  ::  dP_d_l, dP_d_r
      INTEGER  ::  n
      REAL(KIND=8)  :: dP_dre_cap, dP_dr_cap, dP_dre_h_cap, &
                       s_cap, Di, delP, denom
      !------------------------------------------------------------
 
      
      
      !------------------------------------------------------------
      ! Retrieve useful thermodynamic quantities at the left and 
      ! right state and compute jumps
      n = SIZE(ww_l)

      rho_l = ww_l(1);  rho_r = ww_r(1)            ! Density

      re_l = ww_l(n) - 0.5*SUM(ww_l(2:n-1)**2)/ww_l(1) ! Internal energy
      re_r = ww_r(n) - 0.5*SUM(ww_r(2:n-1)**2)/ww_r(1) ! (per unit volume)

      P_l = eos_l%P;              P_r = eos_r%P                ! Pressure
      T_l = eos_l%T;              T_r = eos_r%T                ! Temperature
      h_l = (re_l + P_l) / rho_l; h_r = (re_r + P_r) / rho_r   ! Enthalpy

      Dr   = rho_r - rho_l;  Dre  =  re_r - re_l;  DP   =   P_r - P_l                           !   "  

      h = vv(n) - 0.5*SUM(vv(2:n-1)**2)           ! Specific enthalpy
      !------------------------------------------------------------

      dP_d_l = dP__dre_dr(re_l, rho_l)
      dP_d_r = dP__dre_dr(re_r, rho_r)
      
        dp_dre_cap = (    dP_d_l(1) +     dP_d_r(1)) / 2  ! Trapezoidal rule
         dp_dr_cap = (    dP_d_l(2) +     dP_d_r(2)) / 2  ! Trapezoidal rule
      dp_dre_h_cap = (h_l*dP_d_l(1) + h_r*dP_d_r(1)) / 2  ! Trapezoidal rule
             
             s_cap = dp_dr_cap + dp_dre_h_cap

      IF( ABS(Dr/rho_l) <= 1.d-11 .AND. ABS(Dre/re_l) <= 1.d-11 ) THEN
         dP_dre = dp_dre_cap
          dP_dr = dp_dr_cap
      ELSE
      
           Di  = (s_cap*Dr)**2 + DP**2
         delP  = DP - dp_dr_cap*Dr - dp_dre_cap*Dre
         denom = Di - DP*delP
      
         IF( ABS(denom) < 1.d-25 ) THEN
            WRITE(*,*) ' Null denominator (', denom,'). STOP'
            WRITE(*,*) ' In Vinokur-Montagne scheme, MODULE roe_average'
            WRITE(*,*) ' Dr/r_l = ', Dr/rho_l
            WRITE(*,*) ' Dre/re_l = ', Dre/re_l
            WRITE(*,*) ' DP/P_l = ', DP/P_l
            WRITE(*,*) ' dP = ', delP
            WRITE(*,*) ' Di = ', Di
            STOP
         ENDIF
          dP_dr = (Di*dP_dr_cap + s_cap**2*Dr*delP) / denom
         dP_dre = (Di*dP_dre_cap)                   / denom
      ENDIF

      c2 = dP_dr + h*dP_dre
      c  = SQRT(c2)
      
      rho  = (rho_l + rho_r) / 2
      vv(1) = rho
      T = T__h_v(h, 1/rho)
      eos = eos_ext__T_v(T, 1/rho)

      eos%dP(1) = dP_dre
      eos%dP(2) = dP_dr
      eos%c     = c       
      
     
   END SUBROUTINE vinokur_montagne_average
   !============================================================ 


END MODULE roe_average

