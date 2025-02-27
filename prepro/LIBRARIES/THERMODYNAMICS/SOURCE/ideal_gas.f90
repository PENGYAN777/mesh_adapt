!=====================================================================
!
!      Module: ideal_gas
!
! Description: Ideal gas pressure equation of state and 
!              (compatible) energy equation of state
!
!      Author: Alberto Guardone, Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it,
!                      fossati@aero.polimi.it
!
!   Copyright: 1998-2003 Alberto Guardone, Marco Fossati
!              See COPYING file for copyright notice
!
!=====================================================================

   MODULE  ideal_gas
   
   !------------------------------------------------------------------
   IMPLICIT NONE
   PRIVATE

   ! Module private dimensionless Gas Constant
   REAL(KIND=8) :: Rgas


   PUBLIC :: init_IG, minimum_specific_volume__IG,          &
             P__T_v__IG,  dP__dT_dv__IG, d2P__dT2_dv2__IG,  &
             e__T_v__IG,  de__dT_dv__IG, d2e__dT2_dv2__IG,  &
             s__T_v__IG,                                    &
            cv__T_v__IG, dcv__dT_dv__IG, Rgas
   !------------------------------------------------------------------

   CONTAINS


   SUBROUTINE  init_IG
   !------------------------------------------------------------------
   USE structures,       ONLY: REF
   USE gas_properties,   ONLY: gas
   
   IMPLICIT NONE
   
   REAL(KIND=8) :: Z_ref, P_ref, T_ref, v_ref
   REAL(KIND=8) :: Zc,    Pc,	 Tc,	vc
   !------------------------------------------------------------------

   IF (.NOT. gas % vapor) THEN
      
      Z_ref = 1.d0
      P_ref = REF % P
      T_ref = REF % T
      
      v_ref = Z_ref*REF % Rgas*T_ref / P_ref

      REF % Z = Z_ref
      REF % v = v_ref
      REF % u = SQRT(P_ref*v_ref)

   ELSE

      REF % Rgas = gas % Rgas
      
      Zc = 1.d0
      Pc = gas % Pc  
      Tc = gas % Tc
 
      vc = Zc*REF%Rgas*Tc / Pc

      gas % Zc = Zc

      REF % Z = Zc
      REF % v = vc
      REF % u = SQRT(Pc*vc)

      REF % P = Pc
      REF % T = Tc
      
   ENDIF

   ! Dimensionless Gas constant   
   Rgas = gas % Rgas / REF % Rgas
   
   END SUBROUTINE  init_IG





   FUNCTION  minimum_specific_volume__IG()   RESULT(v_min)

   ! Retrieve the minimum specific volume allowed for the 
   ! considered thermodynamic model.  To be used in Newton 
   ! iterative solutions.
   !------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8) :: v_min
   !------------------------------------------------------------------

   ! Minimum specific volume allowed
   v_min = 1.d-12

   END FUNCTION  minimum_specific_volume__IG




  
  

   !------------------------------------------------------------------
   !  Thermodynamic quantities as functions of T, v
   !------------------------------------------------------------------  

   FUNCTION  P__T_v__IG(T, v)   RESULT(P)

   ! Pressure P(T,v) with T temperature 
   ! and v specific volume
   !-------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T, v  ! Temperature and specific volume
   REAL(KIND=8) 	    :: P     ! Pressure
   !-------------------------------------------------------------------
      
   P = (Rgas*T) / v

   END FUNCTION  P__T_v__IG

   
   
   

   FUNCTION  dP__dT_dv__IG(T, v)   RESULT(dP)

   ! Partial derivative of the pressure P(T,v) with 
   ! respect to the temperature T and the specific volume v
   !-------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)   :: T, v   ! Temperature and specific volume
   REAL(KIND=8), DIMENSION(2) :: dP	! Partial derivatives of 
   					! the pressure P(T,v)
   REAL(KIND=8) :: P_T, P_v
   !-------------------------------------------------------------------		    

   P_T  =      Rgas    / v     
   P_v  =   - (Rgas*T) / v**2 

   dP(1) =  P_T;  dP(2) =  P_v

   END FUNCTION  dP__dT_dv__IG

   
   
   

   FUNCTION  d2P__dT2_dv2__IG(T, v)   RESULT(d2P)

   ! Second order partial derivative of the pressure 
   ! P(T,v) with respect to the temperature T and the specific volume v 
   !-------------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)	:: T, v   ! Temperature and specific volume
   REAL(KIND=8), DIMENSION(2,2) :: d2P    ! Second partial derivatives of 
   					  ! the pressure P(T,v) 
   REAL(KIND=8) :: P_TT, P_Tv, P_vv
   !-------------------------------------------------------------------    

   P_TT =   0.d0
   P_Tv =   -  Rgas    / v**2  
   P_vv =   (2*Rgas*T) / v**3  

   d2P(1,1) = P_TT;	 d2P(1,2) = P_Tv
   d2P(2,1) = d2P(1,2);  d2P(2,2) = P_vv 

   END FUNCTION  d2P__dT2_dv2__IG

   
   
   

   FUNCTION  e__T_v__IG(T, v)   RESULT(e)

   ! Specific internal energy per unit mass e(T,v)
   ! with T temperature and v specific volume
   !------------------------------------------------------------------------
   USE gas_properties,        ONLY: gas
   USE ideal_specific_heat,   ONLY: phi__T
   
   
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T, v  ! Temperature and specific volume
   REAL(KIND=8)             :: e     ! Internal energy per unit mass

   REAL(KIND=8) :: void
   !------------------------------------------------------------------------

   ! Dummy assignment to avoid compilation warnings 
   ! due to the unnecessity of v in ideal gas eos for internal energy
   void = v

   e = gas % e_0  +  phi__T(T)  

   END FUNCTION  e__T_v__IG


  
  

   FUNCTION  de__dT_dv__IG(T, v)   RESULT(de)

   ! Partial derivative of the specific internal 
   ! energy per unit mass e(T,v) with respect to 
   ! the temperature T and the specific volume v
   !------------------------------------------------------------------------
   USE ideal_specific_heat,   ONLY: phip__T
   
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)   :: T, v   ! Temperature and specific volume
   REAL(KIND=8), DIMENSION(2) :: de     ! Partial derivatives of the
                                        ! energy per unit mass e(T,v)
   REAL(KIND=8) :: e_T, e_v, phip

   REAL(KIND=8) :: void
   !------------------------------------------------------------------------
       
   ! Dummy assignment to avoid compilation warnings 
   ! due to the unnecessity of v in ideal gas eos
   ! for internal energy
   void = v


   phip  = phip__T(T)
        
   e_T  =   phip  
   e_v  =   0.d0
          
   de(1) = e_T;  de(2) = e_v

   END FUNCTION  de__dT_dv__IG


   

   FUNCTION  d2e__dT2_dv2__IG(T, v)   RESULT(d2e)

   ! Second order partial derivative of the 
   ! specific internal energy per unit mass 
   ! e(T,v) with respect to the temperature 
   ! T and the specific volume v
   !------------------------------------------------------------------------
   USE ideal_specific_heat,   ONLY: phipp__T
   
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)     :: T, v   ! Temperature and specific volume
   REAL(KIND=8), DIMENSION(2,2) :: d2e    ! Second partial derivatives of the
                                          ! energy per unit mass e(T,v)
   REAL(KIND=8) :: e_TT, e_Tv, e_vv, phipp

   REAL(KIND=8) :: void
   !------------------------------------------------------------------------

   ! Dummy assignment to avoid compilation warnings 
   ! due to the unnecessity of v in ideal gas eos for internal energy
   void = v

   phipp  = phipp__T(T)
                   
   e_TT =   phipp  
   e_Tv =   0.d0
   e_vv =   0.d0
          
   d2e(1,1) = e_TT;  d2e(1,2) = e_Tv
   d2e(2,1) = e_Tv;  d2e(2,2) = e_vv

   END FUNCTION  d2e__dT2_dv2__IG





   FUNCTION  s__T_v__IG(T, v)   RESULT(s)

   ! Specific entropy per unit mass s(T,v) with
   ! T temperature and v specific volume
   !------------------------------------------------------------------------
   USE gas_properties,        ONLY: gas
   USE ideal_specific_heat,   ONLY: psi__T

   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T,v  ! Temperature and specific volume
   REAL(KIND=8)             :: s    ! Specific entropy per unit mass   
   !------------------------------------------------------------------------

   s = gas % s_0 + Rgas*LOG(v) + psi__T(T)    

   END FUNCTION  s__T_v__IG


 
 

   FUNCTION  cv__T_v__IG(T, v)   RESULT(cv)

   ! Specific heat at constant volume cv(T,v),
   ! with T temperature and v specific volume
   !------------------------------------------------------------------------
   USE ideal_specific_heat,   ONLY: phip__T
   
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T, v  ! Temperature and specific volume
   REAL(KIND=8)             :: cv    ! Specific heat at constant volume   

   REAL(KIND=8) :: void
   !------------------------------------------------------------------------

   ! Dummy assignment to avoid compilation warnings 
   ! due to the unnecessity of v in ideal gas eos for entropy
   void = v
 
   cv = phip__T(T)   

   END FUNCTION  cv__T_v__IG





   FUNCTION  dcv__dT_dv__IG(T, v)   RESULT(dcv)

   ! Partial derivatives of the specific heat at constant volume cv(T,v) 
   ! with respect to the temperature T and the specific volume v   
   !------------------------------------------------------------------------
   USE ideal_specific_heat,   ONLY: phipp__T
   
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN)   :: T, v    ! Temperature and specific volume
   REAL(KIND=8), DIMENSION(2) :: dcv     ! Partial derivative of the 
                                         ! specific heat at constant volume   
   REAL(KIND=8) :: cv_T, cv_v
   
   REAL(KIND=8) :: void
   !------------------------------------------------------------------------

   ! Dummy assignment to avoid compilation warnings 
   ! due to the unnecessity of v in ideal gas eos for entropy
   void = v

   cv_T = phipp__T(T)  
   cv_v = 0.d0

   dcv(1) = cv_T;  dcv(2) = cv_v
   
   END FUNCTION  dcv__dT_dv__IG


   END MODULE  ideal_gas
