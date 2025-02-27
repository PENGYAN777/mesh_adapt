!=====================================================================
!
!      Module: thermodynamics
!
! Description: Driver for all thermodynamic functions.
!              Procedures for the determination of 
!              thermodynamic quantities using different
!              models for the pressure equation of state
!              and the specific heat in the dilute limit
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!
!    Modified: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: marco.fossati@polimi.it
!
!   Copyright: 1998-2007 Alberto Guardone, Marco Fossati
!              See COPYING file for copyright notice
!
!=====================================================================

   MODULE thermodynamics

   !------------------------------------------------------------------   
   USE ideal_gas
   USE van_der_waals_gas
   USE martin_hou_gas
   USE soave_redlich_kwong_gas
   USE peng_robinson_gas
   USE redlich_kwong_gas
   USE clausius_II_gas

   IMPLICIT NONE

   INTEGER :: pressure_EOS

   INTEGER, PARAMETER ::   IG = 1, &
                         vdWG = 2, &
                          MHG = 3, &
                         SRKG = 4, &
                          PRG = 5, &
                          RKG = 6, &
                         CIIG = 7 

   CHARACTER(*), DIMENSION(1:7), PARAMETER :: gas_names = &
                        (/ 'IDEAL GAS MODEL              ',  &
                           'VAN DER WAALS GAS MODEL      ',  &
                           'MARTIN-HOU GAS MODEL         ',  &
                           'SOAVE-REDLICH-KWONG GAS MODEL',  &
                           'PENG-ROBINSON GAS MODEL      ',  &
                           'REDLICH-KWONG GAS MODEL      ',  &
                           'CLASIUS II GAS MODEL         '  /)

   REAL(KIND=8), PARAMETER :: NEWTON_EPS      = 1.e-12
   INTEGER,      PARAMETER :: NEWTON_MAX_ITER = 100
   !------------------------------------------------------------------

   CONTAINS


   FUNCTION  eos__e_v(e, v)   RESULT(eos)
   !--------------------------------------------------------
   USE structures

   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: e, v
   TYPE(eos_type)	    :: eos

   REAL(KIND=8) :: T, P
   !--------------------------------------------------------

   T = T__e_v(e, v)
   P = P__T_v(T, v)
   
   eos%T = T;  eos%v = v
   eos%e = e;  eos%P = P
   
   END FUNCTION  eos__e_v





   FUNCTION  eos__T_v(T, v)   RESULT(eos)
   !--------------------------------------------------------
   USE structures
   
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: T, v
   TYPE(eos_type)	    :: eos

   REAL(KIND=8) :: e, P
   !--------------------------------------------------------

   e = e__T_v(T, v)	 
   P = P__T_v(T, v)
   
   eos%T = T;  eos%v = v
   eos%e = e;  eos%P = P

   END FUNCTION  eos__T_v



   

   FUNCTION  eos_ext__e_v(e, v)   RESULT(eos)
   !--------------------------------------------------------
   USE structures
   
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: e, v
   TYPE(eos_ext_type)	    :: eos

   REAL(KIND=8) :: T, P, c
   REAL(KIND=8), DIMENSION(2) :: dP
   !--------------------------------------------------------
   
   T = T__e_v(e, v)
   P = P__T_v(T, v)
   c = c__T_v(T, v)
   dP = dP__dre_dr__T_v(T, v)
   
   eos%T = T;  eos%v = v
   eos%e = e;  eos%P = P
   
   eos%c = c;  eos%dP = dP 

   END FUNCTION  eos_ext__e_v





   FUNCTION  eos_ext__T_v(T, v)   RESULT(eos)
   !--------------------------------------------------------
   USE structures
   
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: T, v
   TYPE(eos_ext_type)	    :: eos

   REAL(KIND=8) :: e, P, c
   REAL(KIND=8), DIMENSION(2) :: dP
   !--------------------------------------------------------

   e = e__T_v(T, v)	 
   P = P__T_v(T, v)
   c = c__T_v(T, v)
   
   dP = dP__dre_dr__T_v(T, v)
   
   eos%T = T;  eos%v = v
   eos%e = e;  eos%P = P
   
   eos%c = c;  eos%dP = dP 

   END FUNCTION  eos_ext__T_v





   FUNCTION  eos_ext__eos(eos)   RESULT(eos_ext)
   !--------------------------------------------------------
   USE structures
   
   IMPLICIT NONE

   TYPE(eos_type), INTENT(IN) :: eos
   TYPE(eos_ext_type)	      :: eos_ext

   REAL(KIND=8) ::  T, e, P, c, v  
   REAL(KIND=8), DIMENSION(2) :: dP
   !--------------------------------------------------------

   e = eos%e;  P = eos%P
   T = eos%T;  v = eos%v
   
   c = c__T_v(T, v)
   dP = dP__dre_dr__T_v(T, v)
   
   eos_ext%T = T;  eos_ext%v = v
   eos_ext%e = e;  eos_ext%P = P
 
   eos_ext%c = c;  eos_ext%dP = dP 

   END FUNCTION  eos_ext__eos





   FUNCTION  P__e_r(e, r)   RESULT(P)
   
   ! Pressure as a function of specific internal energy per 
   ! unit MASS and density
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: e, r
   REAL(KIND=8) :: P

   REAL(KIND=8) :: T, v
   !----------------------------------------------------------------------------

   v = 1.d0/r
   T = T__e_v(e, v)

   P = P__T_v(T, v)

   END FUNCTION  P__e_r





   FUNCTION  P__s_T(s, T)   RESULT(P)
   
   ! Pressure as a function of specific entropy per 
   ! unit MASS and temperature
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: s, T 
   REAL(KIND=8) ::  P

   REAL(KIND=8) :: v
   !----------------------------------------------------------------------------

   v = v__s_T(s,T)
   
   P = P__T_v(T, v)

   END FUNCTION  P__s_T





   FUNCTION  P__c_r(c, r)   RESULT(P)
   
   ! Pressure as a function of speed of sound and density
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: c, r
   REAL(KIND=8) :: P
   
   REAL(KIND=8) :: T, v
   !----------------------------------------------------------------------------
   
   v = 1.d0/r
   T = T__c_v(c, v)

   P = P__T_v(T, v)

   END FUNCTION  P__c_r





   FUNCTION  P__h_r(h, r)   RESULT(P)
   
   ! Pressure as a function of specific enthalpy per 
   ! unit MASS and temperature
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN) :: h, r   
   REAL(KIND=8) ::  P

   REAL(KIND=8) :: T, v
   !----------------------------------------------------------------------------
   
   v = 1.d0/r
   T = T__h_v(h, v)
   
   P = P__T_v(T, v)

   END FUNCTION  P__h_r






 
   FUNCTION  dP__de_dr(e, r)   RESULT(dP)

   ! First order partial derivatives of the pressure as a function
   ! of specific internal energy per unit MASS and density
   !--------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8),   INTENT(IN) :: e, r  	
   REAL(KIND=8), DIMENSION(2) :: dP

   REAL(KIND=8) :: T, v
   !--------------------------------------------------------------------
      
   v = 1.d0/r
   T = T__e_v(e, v)
   
   dP = dP__de_dr__T_v(T, v)

   END FUNCTION  dP__de_dr






   
   FUNCTION  dP__dre_dr(re, r)   RESULT(dP) 

   ! First order partial derivatives of the pressure as a function
   ! of specific internal energy per unit VOLUME and density
   !--------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)   :: re, r 
   REAL(KIND=8), DIMENSION(2) :: dP    
   					 
   REAL(KIND=8) :: T, v
   !--------------------------------------------------------------------
   
   v = 1.d0/r
   T = T__e_v(re/r, v)
   
   dP = dP__dre_dr__T_v(T, v)
 
   END FUNCTION  dP__dre_dr


   
   



   FUNCTION  dP__ds_dv(s, v)   RESULT(dP) 

   ! First order partial derivatives of the pressure as a function
   ! of specific entropy per unit MASS and specific volume
   !--------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)   :: s, v
   REAL(KIND=8), DIMENSION(2) :: dP
   			        
   REAL(KIND=8), DIMENSION(2) :: dP_Tv, ds_Tv
   REAL(KIND=8) :: T
   !--------------------------------------------------------------------
   
   T = T__s_v(s, v)
   
   dP_Tv = dP__dT_dv(T,v);  ds_Tv = ds__dT_dv(T,v)
   
   dP(1) = dP_Tv(1) / ds_Tv(1)
   dP(2) = dP_Tv(2) - dP_Tv(1)*(ds_Tv(2)/ds_Tv(1))
 
   END FUNCTION  dP__ds_dv

   
   
   



   FUNCTION  dP__ds_dT(s, T)   RESULT(dP) 

   ! First order partial derivatives of the pressure as a function
   ! of specific entropy per unit MASS and temperature
   !--------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)   :: s, T
   REAL(KIND=8), DIMENSION(2) :: dP
   					 
   REAL(KIND=8), DIMENSION(2) :: dP_Tv, ds_Tv
   REAL(KIND=8) :: v
   !--------------------------------------------------------------------
   
   v = v__s_T(s, T)
   
   dP_Tv = dP__dT_dv(T,v);  ds_Tv = ds__dT_dv(T,v)
   
   dP(1) = dP_Tv(2) / ds_Tv(2)
   dP(2) = dP_Tv(1) - dP_Tv(2)*(ds_Tv(1)/ds_Tv(2))
 
   END FUNCTION  dP__ds_dT


   
   
   


   FUNCTION  d2P__dre2_dr2(re, r)   RESULT(d2P)

   ! Second order partial derivatives of the pressure as a function
   ! of specific internal energy per unit VOLUME and density
   !--------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)	:: re, r 
   REAL(KIND=8), DIMENSION(2,2) :: d2P
   					   
   REAL(KIND=8) :: T, v
   !--------------------------------------------------------------------

   v = 1/r
   T = T__e_v(re/r, v)
   
   d2P = d2P__dre2_dr2__T_v(T, v)
    
   END FUNCTION  d2P__dre2_dr2







   !============================================================ 
   FUNCTION  e__P_r(P, r)   RESULT(e)
   !============================================================ 

      !----------------------------------------
      ! Specific internal energy per unit mass
      !----------------------------------------
    
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)  ::  P, r  ! Pressure and density 
      REAL(KIND=8)              ::  e     ! Specific internal energy
                                          ! per unit mass
      REAL(KIND=8) ::  T, v  ! Temperature and specific volume

      
      v = 1/r      
      
      T = T__P_v(P, v)

      e = e__T_v(T, v)

   END FUNCTION  e__P_r
   !============================================================ 



   !============================================================ 
   FUNCTION  s__P_r(P, r) RESULT(s)
   !============================================================ 

      !--------------------------------
      ! Specific entropy per unit mass
      !--------------------------------
    
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)  ::  P, r  ! Pressure and density 
      REAL(KIND=8)              ::  s     ! Specific entropy

      REAL(KIND=8) ::  T, v  ! Temperature and specific volume


      v = 1/r
      T = T__P_v(P, v)
      
      s = s__T_v(T, v)
      

   END FUNCTION  s__P_r
   !============================================================ 



   !============================================================ 
   FUNCTION  c__P_r(P, r) RESULT(c)
   !============================================================ 
    
      !----------------
      ! Speed of sound
      !----------------
    
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)  ::  P, r ! Pressure and density 
      REAL(KIND=8)              ::  c    ! Speed of sound

      REAL(KIND=8) ::  T, v  ! Temperature and specific volume


      v = 1/r
      T = T__P_v(P, v)

      c = c__T_v(T, v)

   END FUNCTION  c__P_r
   !============================================================ 



   !============================================================ 
   FUNCTION  G__P_r(P, r) RESULT(G)
   !============================================================ 
  
      !---------------------------------------
      ! Fundamental derivative of gasdynamics 
      ! (dimensionless form)
      !---------------------------------------
    
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)  ::  P, r  ! Pressure and density 
      REAL(KIND=8)              ::  G     ! Fundamental derivative

      REAL(KIND=8) ::  T, v  ! Temperature and specific volume

      
      v = 1/r
      T = T__P_v(P, v)
      
      G = G__T_v(T, v)


   END FUNCTION  G__P_r
   !============================================================ 





   ! Thermodynamic quantities as functions of T, v
   ! Functions are INDEPENDENT from the gas model

   !============================================================ 
   FUNCTION  dP__de_dr__T_v(T, v) RESULT(dP)
   !============================================================ 
    
      !-------------------------------------------------
      ! Partial derivative of the pressure with respect 
      ! to specific internal energy per unit mass and 
      ! the density 
      !-------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)    ::  T, v  ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2)  ::  dP    ! Partial derivatives of 
                                            ! the pressure P(e,r)
      REAL(KIND=8), DIMENSION(2)  ::  dP_Tv, & ! Partial derivatives of P(T,v)
                                      de_Tv    ! Partial derivatives of e(T,v)
      REAL(KIND=8)  :: P_T, P_v, e_T, e_v, r_v, P_e, P_r 
                       ! Shorthands for partial deriv.                              

      
      dP_Tv = dP__dT_dv(T, v)
      de_Tv = de__dT_dv(T, v)
      
      P_T  = dP_Tv(1);  P_v  =  dP_Tv(2)
      e_T  = de_Tv(1);  e_v  =  de_Tv(2)
      r_v  = -1 / (v**2)

      P_e  = P_T / e_T
      P_r  = (P_v - (e_v/e_T)*P_T) / r_v
      
      dP(1) = P_e;  dP(2) = P_r


   END FUNCTION  dP__de_dr__T_v
   !============================================================ 



   !============================================================ 
   FUNCTION dP__dre_dr__T_v(T, v) RESULT( dP ) 
   !============================================================ 

      !-------------------------------------------------
      ! Partial derivative of the pressure with respect 
      ! to specific internal energy per unit volume and 
      ! the density 
      !-------------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)    ::  T, v  ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2)  ::  dP    ! Partial derivatives of 
                                            ! the pressure P(re,r)
      REAL(KIND=8), DIMENSION(2)  ::  dP_Tv, & ! Partial derivatives of P(T,v)
                                      de_Tv    ! Partial derivatives of e(T,v)  

      REAL(KIND=8)  ::  e, P_T, P_v, e_T, e_v, r_v, re_T, re_v, P_re, P_r
                       ! Shorthands for partial deriv.                              

      ! Energy per unit mass
      e = e__T_v(T,v)
      
      ! Partial derivatives of P(T,v) and e(T,v) 
      dP_Tv =  dP__dT_dv(T, v);  de_Tv =  de__dT_dv(T, v)

      P_T =  dP_Tv(1);   P_v =  dP_Tv(2)
      e_T =  de_Tv(1);   e_v =  de_Tv(2)

      ! Partial derivatives of re(T,v) and r(v) 
      re_T = e_T/v;  re_v = e_v/v - e/v**2
      r_v  = - 1/(v**2)

      ! Partial derivatives of P(re,r) 
      P_re =  P_T / re_T
      P_r  = (P_v - (re_v/re_T)*P_T) / r_v
      
      dP(1) =  P_re;  dP(2) =  P_r


   END FUNCTION dP__dre_dr__T_v
   !============================================================
 
 
 
   !============================================================ 
   FUNCTION d2P__dre2_dr2__T_v(T, v) RESULT( d2P ) 
   !============================================================ 

      !-------------------------------------------------
      ! Second order partial derivative of the pressure 
      ! with respect to specific internal energy per 
      ! unit volume and the density 
      !-------------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)      ::  T, v ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P  ! Second partial derivatives of 
                                             ! the pressure P(re,r)
      REAL(KIND=8), DIMENSION(2)    ::  dP_Tv    ! First and second partial 
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P_Tv2  ! derivatives of P(T,v)
      REAL(KIND=8), DIMENSION(2)    ::  de_Tv    ! First and second partial 
      REAL(KIND=8), DIMENSION(2,2)  ::  d2e_Tv2  ! derivatives of e(T,v) 

      REAL(KIND=8)  ::  e,                               &
                         P_T,  P_v,  P_TT,  P_Tv,  P_vv, &
                         e_T,  e_v,  e_TT,  e_Tv,  e_vv, &
                        re_T, re_v, re_TT, re_Tv, re_vv, &
                        r_v,                       r_vv, &
                        P_re, P_r, P_rere, P_rer, P_rr
                        
       
      e = e__T_v(T,v)                  
                        
      ! First order partial derivatives of P(T,v) and e(T,v)       
      dP_Tv =  dP__dT_dv(T, v)
      de_Tv =  de__dT_dv(T, v)

      P_T = dP_Tv(1);   P_v = dP_Tv(2)
      e_T = de_Tv(1);   e_v = de_Tv(2)
      
      ! Second order partial derivatives of P(T,v) and e(T,v)       
      d2P_Tv2 =  d2P__dT2_dv2(T, v)
      d2e_Tv2 =  d2e__dT2_dv2(T, v)

      P_TT = d2P_Tv2(1,1);  P_Tv = d2P_Tv2(2,1);  P_vv = d2P_Tv2(2,2)
      e_TT = d2e_Tv2(1,1);  e_Tv = d2e_Tv2(2,1);  e_vv = d2e_Tv2(2,2)
      
      ! First and second order partial derivatives of re(T,v) and r(v)      
      re_T = e_T/v;  re_v = e_v/v - e/v**2
      r_v = -1 / (v**2)
      
      re_TT = e_TT/v  
      re_Tv = e_Tv/v - e_T/v**2  
      re_vv = (2*e)/v**3 - (2*e_v)/v**2 + e_vv/v 
      r_vv  =  2 / (v**3)

      ! First order partial derivatives of P(re,r)       
      P_re =  P_T / re_T
      P_r  = (P_v - (re_v/re_T)*P_T) / r_v
      
      ! Second order partial derivatives of P(re,r)       
      P_rere = (P_TT - P_re*re_TT) / (re_T**2) 
      P_rer  = (P_Tv - P_rere*re_T*re_v - P_re*re_Tv) / (re_T*r_v)
      P_rr   = (P_vv - (P_rere*re_v + 2*P_rer*r_v)*re_v &
                - P_re*re_vv - P_r*r_vv) / (r_v**2)

      d2P(1,1) = P_rere;    d2P(1,2) = P_rer
      d2P(2,1) = d2P(1,2);  d2P(2,2) = P_rr 


   END FUNCTION d2P__dre2_dr2__T_v
   !============================================================



   !============================================================ 
   FUNCTION d2P__dre_dr__dT_dv(T, v) RESULT( d2P ) 
   !============================================================ 

      !-----------------------------------------------------
      ! Partial derivatives with respect to the temperature
      ! T and the specific volume v of the partial 
      ! derivative of the pressure P with respect to the 
      ! specific internal energy per unit volume re and the 
      ! density r
      !-----------------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)      ::  T, v ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P  ! Second partial derivatives of 
                                             ! the pressure P(re,r) with
                                             ! respect to T and r
      REAL(KIND=8), DIMENSION(2)    ::  dP_Tv    ! First and second partial 
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P_Tv2  ! derivatives of P(T,v)
      REAL(KIND=8), DIMENSION(2)    ::  de_Tv    ! First and second partial 
      REAL(KIND=8), DIMENSION(2,2)  ::  d2e_Tv2  ! derivatives of e(T,v) 

      REAL(KIND=8)  ::  e,                               &
                         P_T,  P_v,  P_TT,  P_Tv,  P_vv, &
                         e_T,  e_v,  e_TT,  e_Tv,  e_vv, &
                        re_T, re_v, re_TT, re_Tv, re_vv, &
                        r_v,                       r_vv, &
                        P_reT, P_rT, P_rev, P_rv
                        
       
      e = e__T_v(T,v)                  
                        
      ! First order partial derivatives of P(T,v) and e(T,v)       
      dP_Tv =  dP__dT_dv(T, v)
      de_Tv =  de__dT_dv(T, v)

      P_T = dP_Tv(1);   P_v = dP_Tv(2)
      e_T = de_Tv(1);   e_v = de_Tv(2)
      
      ! Second order partial derivatives of P(T,v) and e(T,v)       
      d2P_Tv2 =  d2P__dT2_dv2(T, v)
      d2e_Tv2 =  d2e__dT2_dv2(T, v)

      P_TT = d2P_Tv2(1,1);  P_Tv = d2P_Tv2(2,1);  P_vv = d2P_Tv2(2,2)
      e_TT = d2e_Tv2(1,1);  e_Tv = d2e_Tv2(2,1);  e_vv = d2e_Tv2(2,2)
      
      ! First and second order partial derivatives of re(T,v) and r(v)      
      re_T = e_T/v;  re_v = e_v/v - e/v**2
      r_v = -1 / (v**2)
      
      re_TT = e_TT/v  
      re_Tv = e_Tv/v - e_T/v**2  
      re_vv = (2*e)/v**3 - (2*e_v)/v**2 + e_vv/v 
      r_vv  =  2 / (v**3)
      
      ! P_re =  P_T / re_T
      ! P_r  = (P_v - (re_v/re_T)*P_T) / r_v
      ! Partial of P_re and P_r with respect to T and v
      P_reT = (P_TT*re_T - re_TT*P_T) / (re_T**2)
      P_rev = (P_Tv*re_T - re_Tv*P_T) / (re_T**2)
      P_rT  = (P_Tv - ( (re_Tv*P_T + re_v*P_TT)*re_T  &
                         - re_TT*re_v*P_T)/(re_T**2)) / r_v
      P_rv  = (P_vv - ( (re_vv*P_T + re_v*P_Tv)*re_T &
                         - re_Tv*re_v*P_T)/(re_T**2) ) / r_v &
            - (r_vv/(r_v**2))*(P_v - (re_v*P_T)/re_T)

      d2P(1,1) = P_reT;  d2P(1,2) = P_rev
      d2P(2,1) = P_rT;   d2P(2,2) = P_rv


   END FUNCTION d2P__dre_dr__dT_dv
   !============================================================


   !============================================================ 
   FUNCTION  ds__dT_dv(T, v) RESULT(ds)
   !============================================================ 
    
      !---------------------------------
      ! First order partial derivatives 
      ! of the entropy per unit mass
      !---------------------------------

      IMPLICIT NONE

      REAL(KIND=8),   INTENT(IN)  ::  T, v  ! Temperature and specific volume       
      REAL(KIND=8), DIMENSION(2)  ::  ds    ! Partial derivatives of 
                                            ! the entropy s(T,v)

      REAL(KIND=8), DIMENSION(2)  ::  de    ! Partial derivatives of 
                                            ! the energy e(T,v)
                                            
      de = de__dt_dv(T,v)

      ds(1) = cv__T_v(T,v) / T      ! by integrating the definition of T
      ds(2) = (P__T_v(T,v) + de(2)) / T   ! from the Helmholtz potential   


   END FUNCTION  ds__dT_dv
   !============================================================ 


   !============================================================ 
   FUNCTION  dh__dT_dv(T, v) RESULT(dh)
   !============================================================ 
    
      !---------------------------------
      ! First order partial derivatives 
      ! of the enthalpy per unit mass
      !---------------------------------

      IMPLICIT NONE

      REAL(KIND=8),   INTENT(IN)  ::  T, v   ! Temperature and specific volume       
      REAL(KIND=8), DIMENSION(2)  ::  dh     ! Partial derivatives of 
                                             ! the enthalpy h(T,v)

      REAL(KIND=8), DIMENSION(2)  ::  de, dP ! Partial derivatives of 
                                             ! the energy e(T,v) 
                                             ! and pressure P(T,v)
                                            
      dP = dP__dT_dv(T,v);   de = de__dt_dv(T,v)

      dh(1) = dP(1)*v + de(1)           ! From definition h = P*v + e
      dh(2) = dP(2)*v + P__T_v(T,v) + de(2)  
 

   END FUNCTION  dh__dT_dv
   !============================================================ 


   !============================================================ 
   FUNCTION h__T_v(T, v) RESULT (h)
   !============================================================ 

      !----------------
      ! Speed of sound
      !----------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  h     ! Specific enthalpy per unit mass


      h = e__T_v(T, v) + P__T_v(T,v)*v


   END FUNCTION h__T_v
   !============================================================ 



   !============================================================ 
   FUNCTION c__T_v(T, v) RESULT (c)
   !============================================================ 

      !----------------
      ! Speed of sound
      !----------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  c     ! Speed of sound

      REAL(KIND=8)  :: cv                   ! Specific heat at const. vol.
      REAL(KIND=8), DIMENSION(2)    ::  dP  ! Partial derivatives of P(T,v)
      REAL(KIND=8)  :: P_T, P_v   ! Shorthands for partial deriv.

      cv  = cv__T_v(T,v)
      dP  = dP__dT_dv(T,v)
      
      P_T = dP(1)
      P_v = dP(2)

      c = v * SQRT( (T/cv)*P_T**2 - P_v )


   END FUNCTION c__T_v
   !============================================================ 


   !============================================================ 
   FUNCTION c2__T_v(T, v) RESULT (c2)
   !============================================================ 

      !------------------------
      ! Speed of sound squared
      !------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  c2    ! Speed of sound squared

      REAL(KIND=8)  :: cv                   ! Specific heat at const. vol.
      REAL(KIND=8), DIMENSION(2)    ::  dP  ! Partial derivatives of P(T,v)
      REAL(KIND=8)  :: P_T, P_v
                       ! Shorthands for partial deriv.                              


      cv   = cv__T_v(T,v)

      dP  = dP__dT_dv(T,v)
      P_T  =  dP(1);  P_v = dP(2)

      c2 = (v**2)*( (T/cv)*P_T**2 - P_v) 


   END FUNCTION c2__T_v
   !============================================================ 



   !============================================================ 
   FUNCTION G__T_v( T, v) RESULT (G)
   !============================================================ 

      !---------------------------------------
      ! Fundamental derivative of gasdynamics 
      ! (dimensionless form)
      !---------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  G     ! Fundamental derivatives

      REAL(KIND=8)  :: c, cv  ! speed of sound and specific heat at const. vol.
      REAL(KIND=8), DIMENSION(2)    ::  dP    ! First and second partial
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P   ! derivatives of P(T,v)
      REAL(KIND=8), DIMENSION(2)    ::  dcv  ! Partial derivatives of cv(T,v)
      REAL(KIND=8)  :: P_T, P_TT, P_Tv, P_vv, cv_T
                       ! Shorthands for partial deriv.                              

      
      c   = c__T_v(T,v)
      
      cv  = cv__T_v(T,v)
      dcv = dcv__dT_dv(T,v)
      cv_T = dcv(1)
      
      dP  = dP__dT_dv(T,v)
      d2P = d2P__dT2_dv2(T,v)
      
      P_T  = dP(1)
      P_TT = d2P(1,1);  P_Tv = d2P(2,1);  P_vv = d2P(2,2)
            
      G = (v**3/(2*c**2)) &
        * (   P_vv - (3*T*P_T*P_Tv)/cv      &
            + ((T*P_T)/cv)**2 * (3*P_TT + (P_T/T)*(1-(T*cv_T)/cv))  )


   END FUNCTION G__T_v
   !============================================================ 




   FUNCTION  cp__T_v(T, v)   RESULT (cp)
   !------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
   REAL(KIND=8) 	     ::  cp    ! Specific heat at const. pressure

   REAL(KIND=8)  ::  cv                ! Specific heat at const. vol.
   REAL(KIND=8), DIMENSION(2) :: dP    ! Partial derivatives of P(T,v)
   REAL(KIND=8) :: P_T, P_v            ! Shorthands for partial deriv.
   !------------------------------------------------------------
   
   cv  = cv__T_v(T,v)
   dP  = dP__dT_dv(T,v)
   
   P_T = dP(1)
   P_v = dP(2)
   
   cp = cv - T*((P_T**2)/P_v)
   
   END FUNCTION  cp__T_v

  
  




   !============================================================
   !============================================================
   !
   !  Temperature T as a function of the specific volume v and 
   !  of another thermodynamic variable.
   !  Functions are INDEPENDENT from the gas model.
   !  Newton iterative method used for the solution.
   !
   !============================================================
   !============================================================
 
   !============================================================ 
   FUNCTION  T__P_v(P, v)   RESULT(T)
   !============================================================ 

      !-------------------------------------
      ! Temperature T(P,v), with P pressure
      ! and v specific volume
      !-------------------------------------

      USE structures,   ONLY: REF

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN) :: P, v
      REAL(KIND=8)             :: T

      REAL(KIND=8), DIMENSION(2) :: dP
      REAL(KIND=8) :: f
      INTEGER :: it
      
            
      ! Polytropic ideal gas 
      T = P*v / Rgas

      IF (pressure_EOS /= IG) THEN
      
        T = MAX(1.1d0/REF%Z, T)

        DO it = 1, NEWTON_MAX_ITER

          f = P__T_v(T, v) - P         

          IF (ABS(f) < NEWTON_EPS  .AND.  T > 0.d0) RETURN

          dP = dP__dT_dv(T,v)
          T  = T  -  f / dP(1)

        ENDDO
     
        WRITE(*,*) ''
        WRITE(*,*) 'ERROR. T__P_V:'
        WRITE(*,*) 'Newton method failed convergence'
        WRITE(*,*) 'Pressure	    =', P
        WRITE(*,*) 'Specific Volume =', v	     
        WRITE(*,*) 'Temperature     =', T
        WRITE(*,*) ''
        STOP

      ENDIF

   END FUNCTION T__P_v
   !============================================================ 
     
   
  
   !============================================================ 
   FUNCTION  T__e_v(e, v)   RESULT(T)
   !============================================================ 

      !-----------------------------------------------------
      ! Temperature T(e,v), with e internal specific energy 
      ! per unit mass and v specific volume
      !-----------------------------------------------------
      USE structures,      ONLY: REF
      USE gas_properties,  ONLY: gas

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN) :: e, v  ! Energy per unit mass and sp. volume
      REAL(KIND=8)             :: T     ! Temperature  

      REAL(KIND=8), DIMENSION(2) :: de
      REAL(KIND=8) :: f, P
      INTEGER :: it

 
      ! Polytropic ideal gas
      P = ABS(e - gas % e_0)*(gas % gamma - 1) / v
      T = P*v / Rgas
      
      IF (pressure_EOS /= IG) THEN
      
        T = MAX(1.1d0/REF % Z, T)  

        DO it = 1, NEWTON_MAX_ITER

           f = e__T_v(T, v) - e 	

           IF (ABS(f) < NEWTON_EPS  .AND.  T > 0.d0) RETURN

           de = de__dT_dv(T,v)
        	   
           T = T  -  f / de(1)

        ENDDO

        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. T__e_v:'
        WRITE(*,*) 'Newton failed convergence'
        WRITE(*,*) 'Internal Energy =', e
        WRITE(*,*) 'Specific Volume =', v      
        WRITE(*,*) 'Temperature     =', T      
        WRITE(*,*) ' '
        STOP

      ENDIF

   END FUNCTION T__e_v
   !============================================================ 
   


   !============================================================ 
   FUNCTION T__h_v(h, v) RESULT (T)
   !============================================================ 

      !----------------------------------------------
      ! Temperature T(h,v), with h specific enthalpy
      ! per unit mass and v specific volume
      !----------------------------------------------

      USE structures,      ONLY: REF
      USE gas_properties,  ONLY: gas
      
      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  h, v  ! Enthalpy and specific volume
      REAL(KIND=8)              ::  T     ! Temperature  

      REAL(KIND=8), DIMENSION(2)  ::  dh
      REAL(KIND=8)  ::  f, P
      INTEGER  ::  it


      ! Polytropic ideal gas
      P = ABS(h - gas%e_0) * (gas % gamma - 1) / (gas % gamma * v)
      T = P*v / Rgas
      
      IF (pressure_EOS /= IG) THEN 
           
        T = MAX(1.d0/REF%Z, T)   
              
        DO it = 1, NEWTON_MAX_ITER
              
          f = h__T_v(T, v)  -  h			   

          IF (ABS(f) < NEWTON_EPS  .AND.  T > 0.d0) RETURN 

          dh = dh__dT_dv(T,v)
          T = T  -  f / dh(1)

          IF (T < 0) T = 1.d0/REF%Z
              
        ENDDO
      
        WRITE(*,*) ''
        WRITE(*,*) 'ERROR. T__h_v:'
        WRITE(*,*) 'Newton failed convergence.'
        WRITE(*,*) 'Enthalpy	    =', h
        WRITE(*,*) 'Specific Volume =', v	     
        WRITE(*,*) 'Temperature     =', T
        WRITE(*,*) ''
        STOP
      
      ENDIF  

   END FUNCTION T__h_v
   !============================================================ 
    
    
    
   !============================================================ 
   FUNCTION T__s_v(s, v) RESULT (T)
   !============================================================ 

      !----------------------------------------------
      ! Temperature T(h,v), with s specific entropy
      ! per unit mass and v specific volume
      !----------------------------------------------
      USE structures,  ONLY: REF
      USE gas_properties,  ONLY: gas      

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  s, v  ! Entropy and specific volume
      REAL(KIND=8)              ::  T     ! Temperature  

      REAL(KIND=8), DIMENSION(2)  ::  ds
      REAL(KIND=8)  ::  f
      INTEGER  ::  it
     
       
      ! Polytropic ideal gas 
      T = (EXP((s - gas % s_0)/Rgas) / v)**(gas % gamma - 1)
      
      IF (pressure_EOS /= IG) THEN      

        T = (EXP(s)/v)**(1/gas%cv_0)     
        T = MAX(1.d0/REF % Z,T)   
      
        DO it = 1, NEWTON_MAX_ITER
      
           f = s__T_v(T,v) - s         

           IF (ABS(f) < NEWTON_EPS  .AND.  T > 0.d0) RETURN

           ds = ds__dT_dv(T,v)
           
           T = T  -  f / ds(1)

           IF (T < 0) T = 1.d0/REF%Z
      
        ENDDO

        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. T__s_v:'
        WRITE(*,*) 'Newton failed convergence'
        WRITE(*,*) 'Entropy	    =', s
        WRITE(*,*) 'Specific Volume =', v      
        WRITE(*,*) 'Temperature     =', T      
        WRITE(*,*) ' '
        STOP 

      ENDIF

   END FUNCTION T__s_v
   !============================================================ 
    
    
    
   !============================================================ 
   FUNCTION T__c_v(c, v) RESULT (T)
   !============================================================ 

      !-----------------------------------------------
      ! Temperature T(c,v), with c speed of sound and
      ! v specific volume
      !------------------------------------------------
      USE structures,  ONLY: REF
      USE gas_properties,  ONLY: gas

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  c, v  ! Speed of sound and sp. volume 
      REAL(KIND=8)              ::  T     ! Temperature  

      REAL(KIND=8), DIMENSION(2)    ::  dP, dcv
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P
      REAL(KIND=8)  ::  cv, cv_T, c2, P_TT, P_Tv, P_T  
      REAL(KIND=8)  ::  f, fp
      INTEGER  ::  it

      c2 = c**2
      
      ! Polytropic ideal gas 
      T = c2 / (gas % gamma * Rgas)
      
      IF (pressure_EOS /= IG) THEN
      
        T = MAX(1.1d0/REF%Z, T)
      
        DO it = 1, NEWTON_MAX_ITER
      
           f = c2__T_v(T, v) - c2	  

           IF (ABS(f) < NEWTON_EPS  .AND.  T > 0.d0) RETURN

           cv	= cv__T_v(T,v)
           dcv  = dcv__dt_dv(T,v);  cv_T = dcv(1) 

           dP	= dP__dT_dv(T,v)
           d2P  = d2P__dT2_dv2(T,v)
           P_T  = dP(1)
           P_TT = d2P(1,1);  P_Tv = d2P(1,2)
           
           fp = (v**2)*( (P_T/cv)*( P_T*(1 - (T*cv_T)/cv) + 2*T*P_TT ) - P_Tv) ! dc2/dT
           
           T = T  -  f/fp

           IF (T < 0) T = 1.d0/REF % Z
               
        ENDDO

        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. T__c_v:'
        WRITE(*,*) 'Newton failed convergence'
        WRITE(*,*) 'Speed of Sound  =', c
        WRITE(*,*) 'Specific Volume =', v
        WRITE(*,*) 'Temperature     =', T      
        WRITE(*,*) ' '
        STOP
      
      ENDIF

   END FUNCTION T__c_v
   !============================================================ 
 
 
  
   !============================================================
   !============================================================
   !
   !  Specific volume v as a function of temperature T and of
   !  another thermodynamic variable.
   !  Functions are INDEPENDENT from the gas model.
   !  Newton iterative method used for the solution.
   !
   !============================================================
   !============================================================
 
   !============================================================ 
   FUNCTION v__P_T(P, T) RESULT (v)
   !============================================================ 

      !-------------------------------------
      ! Specific volume v(P,T), with P pressure
      ! and T temperature
      !-------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  P, T  ! Pressure and temperature
      REAL(KIND=8)              ::  v     ! Specific volume 

      REAL(KIND=8), DIMENSION(2)  ::  dP
      REAL(KIND=8)  ::  f
      INTEGER  ::  it
      
            
      ! Polytropic ideal gas 
      v = Rgas*T / P

      IF (pressure_EOS /= IG) THEN      

        DO it = 1, NEWTON_MAX_ITER
      
          f = P__T_v(T, v) - P         

          IF (ABS(f) < NEWTON_EPS  .AND.  v > 0.d0) RETURN
          
          dP = dP__dT_dv(T,v)
          v = v  -  f / dP(2)
        	
        ENDDO
     
        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. v__P_T:'
        WRITE(*,*) 'Newton failed convergence'    
        WRITE(*,*) 'Pressure	    =', P
        WRITE(*,*) 'Temperature     =', T      
        WRITE(*,*) 'Specific Volume =', v      
        WRITE(*,*) ' '
        STOP
      
      ENDIF

   END FUNCTION v__P_T
   !============================================================ 
  
     
   !============================================================ 
   FUNCTION v__s_T(s, T) RESULT (v)
   !============================================================ 

      !-------------------------------------
      ! Specific volume v(s,T), with s specific entropy
      ! per unit volume and T temperature
      !-------------------------------------

      USE gas_properties,        ONLY: gas
      USE ideal_specific_heat,   ONLY: psi__T

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  s, T  ! Entropy and temperature
      REAL(KIND=8)              ::  v     ! Specific volume 

      REAL(KIND=8), DIMENSION(2)  ::  ds
      REAL(KIND=8)  ::  f
      INTEGER  ::  it
      
            
      ! Polytropic ideal gas 
      v = EXP((s - gas%s_0)/Rgas) / T**(1/(gas % gamma - 1))

      IF (pressure_EOS /= IG) THEN      

        v = EXP(s - gas%s_0 - psi__T(T))
      
        DO it = 1, NEWTON_MAX_ITER
      
           f = s__T_v(T, v) - s 	

           IF (ABS(f) < NEWTON_EPS  .AND.  v > 0.d0) RETURN
           
           ds = ds__dT_dv(T,v)
        	      
           v = v  -  f / ds(2)
        	
        ENDDO
     
        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. v__s_T:'
        WRITE(*,*) 'Newton failed convergence'
        WRITE(*,*) 'Entropy	    =', s      
        WRITE(*,*) 'Temperature     =', T      
        WRITE(*,*) 'Specific Volume =', v	     
        WRITE(*,*) ' '
        STOP
      
      ENDIF

   END FUNCTION v__s_T
   !============================================================ 
  
     
   !============================================================ 
   FUNCTION dv__ds_dT(s, T)  RESULT(dv) 
   !============================================================ 

      !---------------------------------
      ! First order partial derivatives 
      ! of the pressure
      !---------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)    ::  s, T  ! Entropy (unit mass) and 
                                            ! temperature
      REAL(KIND=8), DIMENSION(2)  ::  dv    ! Partial derivatives of 
                                            ! the pressure v(s,T)
                                            
      REAL(KIND=8), DIMENSION(2)  ::  ds_Tv
      REAL(KIND=8)  ::  v  ! Specific volume

      
      v = v__s_T(s,T)
      
      ds_Tv = ds__dT_dv(T,v)
      
      dv(1) = 1 / ds_Tv(2)
      dv(2) = - ds_Tv(1)/ds_Tv(2)
       
       
   END FUNCTION dv__ds_dT
   !============================================================ 
   
   
   
   !============================================================ 
   !============================================================ 
   !
   !  Specific volume v as a function of (h,P) or (h,s)
   !  Functions are INDEPENDENT from the gas model.
   !  Newton iterative method used for the solution.
   !
   !============================================================ 
   !============================================================ 
    
    
   !============================================================ 
   FUNCTION v__h_P(h, P) RESULT(v)
   !============================================================ 
    
      !----------------
      ! Specific volume
      !----------------
      USE gas_properties,  ONLY: gas

      IMPLICIT NONE

      REAL(KIND=8),   INTENT(IN)  ::  h, P  ! Enthalpy and pressure       
      REAL(KIND=8)                ::  v     ! Specific volume 

      REAL(KIND=8), DIMENSION(2)  ::  dP, dh  ! Partial derivatives of 
                                              ! the pressure P(T,v) and 
                                              ! enthalpy h(T,v)
      REAL(KIND=8)  ::  f, df, T
      INTEGER  ::  it
      
            
      ! Polytropic ideal gas 
      v = (h - gas%e_0)*(gas % gamma - 1)/ (P*gas % gamma)

      IF (pressure_EOS /= IG) THEN            

        T = (h - gas%e_0)/(Rgas + gas%cv_0)
        v = T/P

        DO it = 1, NEWTON_MAX_ITER
      
           T = T__h_v(h,v)
           
           f = P__T_v(T,v) - P      
           
           IF (ABS(f) < NEWTON_EPS  .AND.  v > 0.d0) RETURN
           
           dh = dh__dT_dv(T,v);  dP = dP__dT_dv(T,v)  
           
           df = -(dh(2)/dh(1))*dP(1) + dP(2)	       
        	      
           v = v  -  f / df
     
        ENDDO
     
        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. v__h_P:'
        WRITE(*,*) 'Newton failed convergence'
        WRITE(*,*) 'Enthalpy	    =', h      
        WRITE(*,*) 'Pressure	    =', P      
        WRITE(*,*) 'Specific Volume =', v	     
        WRITE(*,*) ' '
        STOP
      
      ENDIF

   END FUNCTION  v__h_P
   !============================================================ 



   !============================================================ 
   FUNCTION v__s_P(s, P) RESULT(v)
   !============================================================ 
    
      !----------------
      ! Specific volume
      !----------------
      USE gas_properties,  ONLY: gas


      IMPLICIT NONE

      REAL(KIND=8),   INTENT(IN)  ::  s, P  ! Entropy and pressure       
      REAL(KIND=8)                ::  v     ! Specific volume 

      REAL(KIND=8), DIMENSION(2)  ::  dP, ds  ! Partial derivatives of 
                                              ! the pressure P(T,v) and 
                                              ! entropy s(T,v)
      REAL(KIND=8)  ::  f, df, T
      INTEGER  ::  it
      
            
      ! Polytropic ideal gas 
      T = (P * EXP((s - gas%s_0)/Rgas))**(gas%gamma - 1)/gas%gamma
      v = Rgas*T / P

      IF (pressure_EOS /= IG) THEN      
 
        T = (EXP(s)*P)**(1/(gas%cv_0+Rgas)) 
        v = T/P
      
        DO it = 1, NEWTON_MAX_ITER
      
           T = T__s_v(s,v)
           
           f = P__T_v(T,v) - P      
           
           IF (ABS(f) < NEWTON_EPS  .AND.  v > 0.d0) RETURN
           
           ds = ds__dT_dv(T,v);  dP = dP__dT_dv(T,v)  
           
           df = -(ds(2)/ds(1))*dP(1) + dP(2)	       
        	      
           v = v  -  f / df
     
        ENDDO
     
        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. v__s_P:'
        WRITE(*,*) 'Newton failed convergence'
        WRITE(*,*) 'Entropy	    =', s      
        WRITE(*,*) 'Pressure	    =', P
        WRITE(*,*) 'Specific Volume =', v	     
        WRITE(*,*) ' '
        STOP
      
      ENDIF

   END FUNCTION  v__s_P
   !============================================================ 



   !============================================================ 
   FUNCTION v__h_s(h, s) RESULT(v)
   !============================================================ 
    
      !----------------
      ! Specific volume
      !----------------

      USE structures,        ONLY: REF
      USE gas_properties,        ONLY: gas
      USE ideal_specific_heat,   ONLY: psi__T
      
      IMPLICIT NONE

      REAL(KIND=8),   INTENT(IN)  ::  h, s  ! Enthalpy and entropy       
      REAL(KIND=8)                ::  v     ! Specific volume 

      REAL(KIND=8), DIMENSION(2)  ::  ds, dh  ! Partial derivatives of 
                                              ! the entropy e(T,v) and 
                                              ! enthalpy h(T,v)
      REAL(KIND=8)  ::  f, df, T
      INTEGER  ::  it
      
            
      ! Initial guess: polytropic ideal gas 
      T = (h - gas%e_0)/(Rgas + gas%cv_0)
      T = MAX(T,1.d0/REF%Z)
      v = EXP(s - gas%s_0 - psi__T(T))/Rgas

      
!       IF (pressure_EOS /= IG) THEN
! 
!         T = (h - gas%e_0)/(Rgas + gas%cv_0)
!         T = MAX(T,1.d0/REF%Z)
!         v = EXP(s - gas%s_0 - psi__T(T))/Rgas

        DO it = 1, NEWTON_MAX_ITER
      
           T = T__h_v(h,v)
           
           f = s__T_v(T,v) - s     

           IF (ABS(f) < NEWTON_EPS  .AND.  v > 0.d0) RETURN
           
           dh = dh__dT_dv(T,v);  ds = ds__dT_dv(T,v)  
           
           df = -(dh(2)/dh(1))*ds(1) + ds(2)	       
        	      
           v = v  -  f / df
        	
        ENDDO
     
        WRITE(*,*) ' '     
        WRITE(*,*) 'ERROR. v__h_s:'
        WRITE(*,*) 'Newton failed convergence'
        WRITE(*,*) 'Enthalpy	    =', h
        WRITE(*,*) 'Entropy	    =', s
        WRITE(*,*) 'Specific Volume =', v	     
        WRITE(*,*) ' '
        STOP
      
!       ENDIF

   END FUNCTION  v__h_s
   !============================================================ 

  





   !============================================================
   !============================================================
   !
   !  Thermodynamic quantities as functions of T, v
   !  Functions are DEPENDENT from the gas model
   !
   !============================================================
   !============================================================

   !============================================================ 
   FUNCTION P__T_v(T, v) RESULT (P)
   !============================================================ 

      !------------------------------------
      ! Pressure P(T,v) with T temperature 
      ! and v specific volume
      !------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  P     ! Pressure 


      SELECT CASE(pressure_EOS)

         CASE(IG)
            P = P__T_v__IG( T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. P__T_V:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP
	    	    
!          CASE(vdWG)
!             P = P__T_v__vdWG( T, v)
!            
!          CASE(MHG)
!             P = P__T_v__MHG( T, v)
!            
!          CASE(SRKG)
!             P = P__T_v__SRKG( T, v)
!            
!          CASE(PRG)
!             P = P__T_v__PRG( T, v)
!            
!          CASE(RKG)
!             P = P__T_v__RKG( T, v)
!            
!          CASE(CIIG)
!             P = P__T_v__CIIG( T, v)
           
      END SELECT


   END FUNCTION P__T_v
   !============================================================ 
 
 
  
   !============================================================ 
   FUNCTION dP__dT_dv(T, v)  RESULT(dP)
   !============================================================ 

      !-------------------------------------------------
      ! Partial derivative of the pressure P(T,v) with 
      ! respect to the temperature T and the specific 
      ! volume v
      !-------------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)    ::  T, v ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2)  ::  dP   ! Partial derivatives of 
                                           ! the pressure P(T,v)

      SELECT CASE(pressure_EOS)

         CASE(IG)
            dP = dP__dT_dv__IG(T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. DP__DT_DV:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             dP = dP__dT_dv__vdWG(T, v)
!             
!          CASE(MHG)
!             dP = dP__dT_dv__MHG(T, v)
!             
!          CASE(SRKG)
!             dP = dP__dT_dv__SRKG(T, v)
!             
!          CASE(PRG)
!             dP = dP__dT_dv__PRG(T, v)
!             
!          CASE(RKG)
!             dP = dP__dT_dv__RKG(T, v)
!             
!          CASE(CIIG)
!             dP = dP__dT_dv__CIIG(T, v)
            
      END SELECT

       
   END FUNCTION dP__dT_dv
   !============================================================ 



   !============================================================ 
   FUNCTION d2P__dT2_dv2(T, v)  RESULT(d2P) 
   !============================================================ 

      !-------------------------------------------------
      ! Second order partial derivative of the pressure 
      ! P(T,v) with respect to the temperature T and the 
      ! specific volume v 
      !-------------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)      ::  T, v ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2,2)  ::  d2P  ! Second partial derivatives of 
                                             ! the pressure P(T,v)

      SELECT CASE(pressure_EOS)

         CASE(IG)
            d2P = d2P__dT2_dv2__IG(T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. D2P__DT2_DV2:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             d2P = d2P__dT2_dv2__vdWG(T, v)
!             
!          CASE(MHG)
!             d2P = d2P__dT2_dv2__MHG(T, v)
!             
!          CASE(SRKG)
!             d2P = d2P__dT2_dv2__SRKG(T, v)
!             
!          CASE(PRG)
!             d2P = d2P__dT2_dv2__PRG(T, v)
!             
!          CASE(RKG)
!             d2P = d2P__dT2_dv2__RKG(T, v)
!             
!          CASE(CIIG)
!             d2P = d2P__dT2_dv2__CIIG(T, v)
            
      END SELECT


   END FUNCTION d2P__dT2_dv2
   !============================================================ 
   


   !============================================================ 
   FUNCTION e__T_v(T, v) RESULT(e)
   !============================================================ 

      !-----------------------------------------------
      ! Specific internal energy per unit mass e(T,v)
      ! with T temperature and v specific volume
      !-----------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  e     ! Internal energy per unit mass   

      SELECT CASE(pressure_EOS)

         CASE(IG)
            e = e__T_v__IG(T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. E__T_V:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             e = e__T_v__vdWG(T, v)
!            
!          CASE(MHG)
!             e = e__T_v__MHG(T, v)
!            
!          CASE(SRKG)
!             e = e__T_v__SRKG(T, v)
!            
!          CASE(PRG)
!             e = e__T_v__PRG(T, v)
!            
!          CASE(RKG)
!             e = e__T_v__RKG(T, v)
!            
!          CASE(CIIG)
!             e = e__T_v__CIIG(T, v)
           
      END SELECT


   END FUNCTION e__T_v
   !============================================================
   
   

   !============================================================ 
   FUNCTION de__dT_dv( T, v ) RESULT( de ) 
   !============================================================ 

      !---------------------------------------------
      ! Partial derivative of the specific internal 
      ! energy per unit mass e(T,v) with respect to 
      ! the temperature T and the specific volume v
      !---------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)    ::  T, v  ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2)  ::  de    ! Partial derivatives of the
                                            ! energy per unit mass e(T,v)

      SELECT CASE(pressure_EOS)

         CASE(IG)
            de = de__dT_dv__IG(T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. DE__DT_DV:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             de = de__dT_dv__vdWG(T, v)
!             
!          CASE(MHG)
!             de = de__dT_dv__MHG(T, v)
!             
!          CASE(SRKG)
!             de = de__dT_dv__SRKG(T, v)
!             
!          CASE(PRG)
!             de = de__dT_dv__PRG(T, v)
!             
!          CASE(RKG)
!             de = de__dT_dv__RKG(T, v)
!             
!          CASE(CIIG)
!             de = de__dT_dv__CIIG(T, v)
            
      END SELECT


   END FUNCTION de__dT_dv
   !============================================================ 



   !============================================================ 
   FUNCTION d2e__dT2_dv2(T, v) RESULT( d2e ) 
   !============================================================ 

      !----------------------------------------
      ! Second order partial derivative of the 
      ! specific internal energy per unit mass 
      ! e(T,v) with respect to the temperature 
      ! T and the specific volume v
      !----------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)      ::  T, v  ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2,2)  ::  d2e   ! Second partial derivatives of the
                                              ! energy per unit mass e(T,v)

      SELECT CASE(pressure_EOS)

         CASE(IG)
            d2e = d2e__dT2_dv2__IG(T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. D2E__DT2_DV2:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             d2e = d2e__dT2_dv2__vdWG(T, v)
!             
!          CASE(MHG)
!             d2e = d2e__dT2_dv2__MHG(T, v)
!             
!          CASE(SRKG)
!             d2e = d2e__dT2_dv2__SRKG(T, v)
!             
!          CASE(PRG)
!             d2e = d2e__dT2_dv2__PRG(T, v)
!             
!          CASE(RKG)
!             d2e = d2e__dT2_dv2__RKG(T, v)
!             
!          CASE(CIIG)
!             d2e = d2e__dT2_dv2__CIIG(T, v)
            
      END SELECT


   END FUNCTION d2e__dT2_dv2
   !============================================================ 
 

    
   !============================================================ 
   FUNCTION s__T_v( T, v) RESULT (s)
   !============================================================ 

      !--------------------------------------------
      ! Specific entropy per unit mass s(T,v) with
      ! T temperature and v specific volume
      !--------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  s     ! Specific entropy per unit mass


      SELECT CASE(pressure_EOS)

         CASE(IG)
            s = s__T_v__IG( T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. S__T_V:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             s = s__T_v__vdWG( T, v)
!            
!          CASE(MHG)
!             s = s__T_v__MHG( T, v)
!            
!          CASE(SRKG)
!             s = s__T_v__SRKG( T, v)
!            
!          CASE(PRG)
!             s = s__T_v__PRG( T, v)
!            
!          CASE(RKG)
!             s = s__T_v__RKG( T, v)
!            
!          CASE(CIIG)
!             s = s__T_v__CIIG( T, v)
           
      END SELECT


   END FUNCTION s__T_v
   !============================================================ 

  
  
   !============================================================ 
   FUNCTION  cv__T_v(T, v)   RESULT (cv)
   !============================================================ 

      !------------------------------------------
      ! Specific heat at constant volume cv(T,v),
      ! with T temperature and v specific volume
      !------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
      REAL(KIND=8)              ::  cv    ! Specific heat at constant volume


      SELECT CASE(pressure_EOS)

         CASE(IG)
            cv = cv__T_v__IG(T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. CV__T_V:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             cv = cv__T_v__vdWG(T, v)
!            
!          CASE(MHG)
!             cv = cv__T_v__MHG(T, v)
!            
!          CASE(SRKG)
!             cv = cv__T_v__SRKG(T, v)
!            
!          CASE(PRG)
!             cv = cv__T_v__PRG(T, v)
!            
!          CASE(RKG)
!             cv = cv__T_v__RKG(T, v)
!            
!          CASE(CIIG)
!             cv = cv__T_v__CIIG(T, v)
           
      END SELECT


   END FUNCTION cv__T_v
   !============================================================ 

  
  
   !============================================================ 
   FUNCTION dcv__dT_dv( T, v) RESULT (dcv)
   !============================================================ 

      !------------------------------------------
      ! Partial derivatives of the specific heat 
      ! at constant volume cv(T,v) with respect
      ! to the temperature T and the specific
      ! volume v
      !------------------------------------------

      IMPLICIT NONE 

      REAL(KIND=8), INTENT(IN)    ::  T, v  ! Temperature and specific volume
      REAL(KIND=8), DIMENSION(2)  ::  dcv   ! Partial derivative of the
                                            ! specific heat at constant volume


      SELECT CASE(pressure_EOS)

         CASE(IG)
            dcv = dcv__dT_dv__IG(T, v)

         CASE DEFAULT
	    WRITE(*,*) ''
	    WRITE(*,*) 'ERROR. DCV__DT_DV:'
	    WRITE(*,*) 'Pressure EOS, not known.'
	    WRITE(*,*) ''
	    STOP


!          CASE(vdWG)
!             dcv = dcv__dT_dv__vdWG(T, v)
!            
!          CASE(MHG)
!             dcv = dcv__dT_dv__MHG(T, v)
!            
!          CASE(SRKG)
!             dcv = dcv__dT_dv__SRKG(T, v)
!            
!          CASE(PRG)
!             dcv = dcv__dT_dv__PRG(T, v)
!            
!          CASE(RKG)
!             dcv = dcv__dT_dv__RKG(T, v)
!            
!          CASE(CIIG)
!             dcv = dcv__dT_dv__CIIG(T, v)
           
      END SELECT


   END FUNCTION dcv__dT_dv
   !============================================================ 

  
  
   !============================================================ 
   !============================================================ 
   !
   !  Initialization procedures
   !
   !============================================================ 
   !============================================================
   
   !============================================================ 
   FUNCTION minimum_specific_volume() RESULT (v_min)
   !============================================================ 

      !------------------------------------------
      ! Retrieve the minimum specific volume
      ! allowed for the considered thermodynamic 
      ! model.  To be used in Newton iterative
      ! solutions.
      !------------------------------------------

      IMPLICIT NONE 

      ! Minimum specific volume allowed
      REAL(KIND=8) :: v_min   


      SELECT CASE(pressure_EOS)

         CASE(IG)
         v_min = minimum_specific_volume__IG()

         CASE DEFAULT
	 WRITE(*,*) ''
	 WRITE(*,*) 'ERROR. MINIMUM_SPECIFIC_VOLUME:'
	 WRITE(*,*) 'Pressure EOS, not known.'
	 WRITE(*,*) ''
	 STOP

!          CASE(vdWG)
!          v_min = minimum_specific_volume__vdWG()
!          CASE(MHG)
!          v_min = minimum_specific_volume__MHG()
!          CASE(SRKG)
!          v_min = minimum_specific_volume__SRKG()
!          CASE(PRG)
!          v_min = minimum_specific_volume__PRG()
!          CASE(RKG)
!          v_min = minimum_specific_volume__RKG()
!          CASE(CIIG)
!          v_min = minimum_specific_volume__CIIG()
           
      END SELECT
     

   END FUNCTION minimum_specific_volume
   !============================================================ 





   SUBROUTINE  read_param_thermodynamics(idf)
   !------------------------------------------------------------
   USE gas_properties,        ONLY: read_param_gas_properties
   USE ideal_specific_heat,   ONLY: read_param_ideal_specific_heat
   
   IMPLICIT NONE

   INTEGER :: idf
   !------------------------------------------------------------

   READ(idf,*);  READ(idf,*);  READ(idf,*)
   CALL read_param_gas_properties(idf)
   
   READ(idf,*);  READ(idf,*);  READ(idf,*)
   READ(idf,*) pressure_EOS
   
   READ(idf,*);  READ(idf,*);  READ(idf,*)
   CALL read_param_ideal_specific_heat(idf)

   END SUBROUTINE  read_param_thermodynamics





   SUBROUTINE  init_thermodynamics
   !------------------------------------------------------------
   USE gas_properties,        ONLY: gas
   USE structures,            ONLY: REF
   USE ideal_specific_heat,   ONLY: init_ideal_specific_heat,  &
                                    phip__T
   
   IMPLICIT NONE
   !------------------------------------------------------------

   SELECT CASE (pressure_EOS)

     CASE (IG)
     CALL init_IG

!      CASE (vdWG)
!      CALL init_vdWG
!      
!      CASE (MHG)
!      CALL init_MHG
!      
!      CASE (SRKG)
!      CALL init_SRKG
!      
!      CASE (PRG)
!      CALL init_PRG
!      
!      CASE (RKG)
!      CALL init_RKG
!      
!      CASE (CIIG)
!      CALL init_CIIG

     CASE DEFAULT
     WRITE(*,*) ''
     WRITE(*,*) 'ERROR. INIT_THERMODYNAMICS:'
     WRITE(*,*) 'Pressure EOS, not known.'
     WRITE(*,*) ''
     STOP
         
   END SELECT
   
   CALL init_ideal_specific_heat

   ! Reference State (not to be confused with Scaling)
   gas % e_0  = 0.d0
   gas % s_0  = 0.d0


   IF (pressure_EOS == IG) THEN
     gas % cv_0 = phip__T(1.d-20)
   ELSE  
     gas % cv_0 = phip__T(1/REF%Z)
   ENDIF

   END SUBROUTINE  init_thermodynamics





   SUBROUTINE  write_param_thermodynamics(idf)
   !--------------------------------------------------------------------
   USE gas_properties,        ONLY: write_param_gas_properties
   USE ideal_specific_heat,   ONLY: write_param_ideal_specific_heat
   
   IMPLICIT NONE
   INTEGER :: idf
   !--------------------------------------------------------------------

   CALL write_param_gas_properties(idf)

   WRITE(idf,*) '  ', pressure_EOS, ' ---> ', gas_names(pressure_EOS)

   CALL write_param_ideal_specific_heat(idf)
   
   END SUBROUTINE  write_param_thermodynamics 


END MODULE thermodynamics
