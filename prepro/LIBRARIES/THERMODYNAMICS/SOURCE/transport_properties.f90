!===============================================================================
!
!      Module: transport_properties
!
! Description: Procedures for the evaluation of transport
!              properties (viscosity and thermal 
!              conductivity) 
!
!      Author: Alberto Guardone, Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!                      fossati@aero.polimi.it
!
!   Copyright: 1998-2003 Alberto Guardone
!              See COPYING file for copyright notice
!
!=============================================================================== 

   MODULE  transport_properties

   !----------------------------------------------------------------------------
   IMPLICIT NONE
   PRIVATE

   ! Gas parameters (stored locally)
   REAL(KIND=8) :: Rgas, M_mole, omega, mur, kH
   
   ! Scaling reference values (stored locally)
   REAL(KIND=8) :: T_ref, P_ref, v_ref, Z_ref, cp_ref,  &
                   mu_ref, kk_ref, Tc, Pc, vc, Zc

   ! Transport models
   INTEGER :: viscosity_model
   INTEGER :: thermal_cond_model
   
   INTEGER, PARAMETER  ::  NO_VISCOSITY         = 0, &
                           CONSTANT_VISC        = 1, &
                           DILUITE_GAS_VISC     = 2, &
                           DENSE_GAS_VISC       = 3, &
                           CORRESP_STATE_VISC   = 4, &
                           SUTHERLAND_LAW_VISC  = 5, &
                           POWER_LAW_VISC       = 6
   
   INTEGER, PARAMETER  ::  NO_TCOND             = 0, &
                           CONSTANT_TCOND       = 1, &
                           DILUITE_GAS_TCOND    = 2, &
                           DENSE_GAS_TCOND      = 3, &
                           CONST_PRANDTL_TCOND  = 4

   ! Coefficients for transport models (DIMENSIONAL)
   REAL(KIND=8) :: mu_v
   REAL(KIND=8) :: mu_const, kk_const
   REAL(KIND=8) :: mu_ref_S, T_ref_S, Tchar_S	
   REAL(KIND=8) :: mu_ref_L, T_ref_L, expo
   REAL(KIND=8) :: Pr	
   REAL(KIND=8) :: Fc, beta
   
   REAL(KIND=8), DIMENSION(10) :: AA
   REAL(KIND=8), DIMENSION(7)  :: BB
   
   
   PUBLIC :: read_param_transport_properties, &
	     init_transport_properties,       &
             write_param_transport_propert,   &
             transport_coefficients,	      &	     
             viscosity_model,                 &
	     thermal_cond_model, Pr,          &
	     mu_ref_S, T_ref_S,               &
	     mu_ref_L, T_ref_L, mu_const,     &
	     mu_v
   !----------------------------------------------------------------------------

   CONTAINS


   SUBROUTINE  read_param_transport_properties(idf)   
   !----------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: idf
   !----------------------------------------------------------------------------

   READ(idf,*);  READ(idf,*);  READ(idf,*)

   ! Viscosity Models
   READ(idf,*) viscosity_model 

   SELECT CASE (viscosity_model)

     CASE (NO_VISCOSITY)
     
     CASE (CONSTANT_VISC)
     READ(idf,*)  mu_const
     	
     CASE (DILUITE_GAS_VISC)
     CASE (DENSE_GAS_VISC)
     CASE (CORRESP_STATE_VISC)
     	
     CASE (SUTHERLAND_LAW_VISC)
     READ(idf,*)  mu_ref_S, T_ref_S, Tchar_S

     CASE (POWER_LAW_VISC)
     READ(idf,*)  mu_ref_L, T_ref_L, expo

   END SELECT

   ! Bulk viscosity
   IF (viscosity_model /= 0) READ(idf,*) mu_v
   
   ! Thermal conductivity Models
   READ(idf,*) thermal_cond_model 

   SELECT CASE (thermal_cond_model)

     CASE (NO_TCOND)

     CASE (CONSTANT_TCOND)
     READ(idf,*)  kk_const
     	
     CASE (DILUITE_GAS_TCOND)
     CASE (DENSE_GAS_TCOND)

     CASE (CONST_PRANDTL_TCOND)
     READ(idf,*)  Pr
		  
   END SELECT
   
   END SUBROUTINE  read_param_transport_properties





   SUBROUTINE  init_transport_properties   
   !----------------------------------------------------------------------------
   USE structures,       ONLY: REF
   USE gas_properties,   ONLY: gas
   
   IMPLICIT NONE

   ! Constants for Chung et al. method
   ! Ind. Eng. Chem. Res. 1988, 27, 671-679

   ! Coefficients for the correction Fc for polyatomic polar moleculas
   REAL(KIND=8), PARAMETER  ::  f0 =  1.d0,	  f1 = -0.2756d0,   &
   				f2 =  0.59035d-1, f3 =  1.d0		
   
   ! Dense gas viscosity correlation coefficients 
   REAL(KIND=8), DIMENSION(10), PARAMETER  ::  a0 = &
      (/  6.32402d0,  0.12102d-2, 5.28346d0,  6.62263d0,  1.97454d1,  &
   	 -1.89992d0,  2.42745d1 , 0.79716d0, -0.23816d0,  0.68629d-1 /)

   REAL(KIND=8), DIMENSION(10), PARAMETER  ::  a1 = & 
      (/  5.04119d1, -0.11536d-2, 2.54209d2,  3.80957d1,  7.63034d0,  &
   	 -1.25367d1,  3.44945d0 , 1.11764d0,  0.67695d-1, 0.34793d0  /)

   REAL(KIND=8), DIMENSION(10), PARAMETER  ::  a2 = & 
      (/ -5.16801d1, -0.62571d-2,-1.68481d2, -8.46414d0, -1.43544d1,  &
   	  4.98529d0, -1.12913d1,  0.12348d-1,-0.81630d0,  0.59256d0  /)

   REAL(KIND=8), DIMENSION(10), PARAMETER  ::  a3 = & 
      (/  1.18902d3,  0.37283d-1, 3.89827d3,  3.14178d1,  3.15267d1,  &
   	 -1.81507d1,  6.93466d1, -4.11661d0,  4.02528d0, -0.72663d0  /)


   ! Dense gas thermal conductivity correlation coefficients 
   REAL(KIND=8), DIMENSION(7),  PARAMETER  ::  b0 = &
      (/  2.41657d0, -0.50924d0,  6.61069d0,  1.45425d1,  0.79274d0,  &
   	 -5.86340d0,  8.11710d1 				     /)
   REAL(KIND=8), DIMENSION(7),  PARAMETER  ::  b1 = &
      (/  0.74824d0, -1.50936d0,  5.62073d0, -8.91387d0,  0.82019d0,  &
   	  1.28005d1,  1.14158d2 				     /)
   REAL(KIND=8), DIMENSION(7),  PARAMETER  ::  b2 = &
      (/ -0.91858d0, -4.99912d1,  6.47599d1, -5.63794d0, -0.69369d0,  &
   	  9.58926d0, -6.08410d1 				     /)
   REAL(KIND=8), DIMENSION(7),  PARAMETER  ::  b3 = &
      (/  1.21721d2,  6.99834d1,  2.70389d1,  7.43435d1,  6.31734d0,  &
   	 -6.55292d1,  4.66775d2 				     /)
   !----------------------------------------------------------------------------


   IF ( ( (viscosity_model    .EQ.  DENSE_GAS_VISC)   .AND. &
   	  (thermal_cond_model .NE. DENSE_GAS_TCOND) ) .OR.   &
   	( (viscosity_model    .EQ.  DILUITE_GAS_VISC)	.AND. &
   	  (thermal_cond_model .NE. DILUITE_GAS_TCOND) ) ) THEN
   	  
      WRITE(*,*) 'Viscosity and thermal conductivity model'
      WRITE(*,*) 'do not match.  In init_trasport_properties. STOP'
      STOP
      
   ENDIF     

   IF ( (viscosity_model    .EQ.  DILUITE_GAS_VISC)   .AND. &
   	(.NOT. gas % vapor) ) THEN
   
      WRITE(*,*) 'The gas ', gas%name,' is not a vapor:'
      WRITE(*,*) 'Diluite gas viscosity model is not applicable. '
      WRITE(*,*) 'in init_transport_properties.  STOP. '
      STOP
      
   ENDIF     

   IF ( (viscosity_model    .EQ.  DENSE_GAS_VISC)   .AND. &
   	(.NOT. gas % vapor) ) THEN
   
      WRITE(*,*) 'The gas ', gas%name,' is not a vapor:'
      WRITE(*,*) 'Dense gas viscosity model is not applicable. '
      WRITE(*,*) 'in init_transport_properties.  STOP. '
      STOP
      
   ENDIF

   ! If Euler computations then skip initialization
   IF (viscosity_model == 0  .AND.  thermal_cond_model == 0) THEN
     
     REF % Re = 0.d0
     REF % Pr = 0.d0
     RETURN
     
   ENDIF


   ! Gas properties
   M_mole = gas % M_mole
   Rgas   = gas % Rgas

   IF (gas % vapor) THEN

     Tc = gas % Tc
     Pc = gas % Pc
     Zc = gas % Zc     
     vc = Zc*Rgas*Tc / Pc

     omega  = gas % omega
     mur    = gas % mur
     kH     = gas % kH

   ELSE

     T_ref  = REF % T
     P_ref  = REF % P
     v_ref  = REF % v
     Z_ref  = REF % Z
     cp_ref = REF % cp
     
   ENDIF


   REF % Re = 1.d0
   REF % Pr = 1.d0

   REF % mu = REF % u * REF % L / (REF % Re * REF % v)
   REF % k  = REF % mu * REF % cp / REF % Pr

   ! Local
   mu_ref = REF % mu
   kk_ref = REF % k



   ! Viscosity
   SELECT CASE (viscosity_model)
      
     CASE (CONSTANT_VISC)
   	
     CASE (DILUITE_GAS_VISC)
     Fc =  f0  +  f1*omega  +  f2*mur**4  +  f3*kH
   	
     CASE (DENSE_GAS_VISC)
     Fc =  f0  +  f1*omega  +  f2*mur**4  +  f3*kH
     AA =  a0  +  a1*omega  +  a2*mur**4  +  a3*kH

     CASE (CORRESP_STATE_VISC)
     CASE (SUTHERLAND_LAW_VISC)

   END SELECT


   ! Thermal conductivity
   SELECT CASE (thermal_cond_model)
   
     CASE (CONSTANT_TCOND)
   	
     CASE (DILUITE_GAS_TCOND)
     beta  = 0.7862d0 - 0.7109d0*omega + 1.3168d0*omega**2
   	
     CASE (DENSE_GAS_TCOND)
     beta  = 0.7862d0 - 0.7109d0*omega + 1.3168d0*omega**2
     BB    =  b0  +  b1*omega  +  b2*mur**4  +  b3*kH
   				
   END SELECT

!    IF (gas % vapor) THEN
!      CALL transport_coefficients_dim(Tc,    vc,     mu_ref, kk_ref)
!    ELSE
!      CALL transport_coefficients_dim(T_ref, v_ref,  mu_ref, kk_ref)
!    ENDIF
! 
!    gas % mu_ref = mu_ref
!    gas % kk_ref = kk_ref
   
   END SUBROUTINE  init_transport_properties





   SUBROUTINE  write_param_transport_propert(idf)
   !----------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: idf
   !----------------------------------------------------------------------------

   ! Viscosity
   WRITE(idf,*) '  ', viscosity_model 

   SELECT CASE (viscosity_model)
      
     CASE (CONSTANT_VISC)
     WRITE(idf,*) '', mu_const, '[Pa s]   -   Reference dinamical viscosity'
     	
     CASE (DILUITE_GAS_VISC)
     CASE (DENSE_GAS_VISC)
     CASE (CORRESP_STATE_VISC)
     	
     CASE (SUTHERLAND_LAW_VISC)
     WRITE(idf,*) '', mu_ref_S, T_ref_S, Tchar_S, '[Pa s] [K] [K]   -	Reference mu and T at Shuterland state'

     CASE (POWER_LAW_VISC)
     WRITE(idf,*) '', mu_ref_L, T_ref_L, expo, '[Pa s] [K]   -   Reference viscosity and Temperature'

   END SELECT
   
   IF (viscosity_model /= 0)  WRITE(idf,*) mu_v


   ! Thermal conductivity
   WRITE(idf,*) '  ', thermal_cond_model 

   SELECT CASE (thermal_cond_model)

     CASE (CONSTANT_TCOND)
     WRITE(idf,*) '', kk_const
     	
     CASE (DILUITE_GAS_TCOND)
     CASE (DENSE_GAS_TCOND)

     CASE (CONST_PRANDTL_TCOND)
     WRITE(idf,*) '', Pr
     			       
   END SELECT
   
   END SUBROUTINE  write_param_transport_propert
   
   
   
   
   
   SUBROUTINE  transport_coefficients(T_adim, v_adim, mu_adim, kk_adim)
   
   ! The viscosity and the thermal conductivity are computed 
   ! using dimensional temperature and specific volume
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN)  :: T_adim, v_adim
   REAL(KIND=8), INTENT(OUT) :: mu_adim, kk_adim
   
   REAL(KIND=8) :: T, v, mu0, kk0, mu, kk
   !----------------------------------------------------------------------------
   
   ! Dimensional temperature and specific volume
   T = T_adim * T_ref
   v = v_adim * v_ref

   ! Viscosity
   SELECT CASE (viscosity_model)
      
     CASE (CONSTANT_VISC)
     mu  = mu_const
        
     CASE (DILUITE_GAS_VISC)
     mu  = diluite_gas_viscosity(T)
        
     CASE (DENSE_GAS_VISC)
     mu0 = diluite_gas_viscosity(T)
     mu  = dense_gas_viscosity(T,v,mu0)
     
     CASE (CORRESP_STATE_VISC)
     mu  = corresp_state_viscosity(T)
        
     CASE (SUTHERLAND_LAW_VISC)
     mu  = Sutherland_law_viscosity(T)

     CASE (POWER_LAW_VISC)
     mu  = power_law_viscosity(T)

   END SELECT


   ! Thermal conductivity
   SELECT CASE (thermal_cond_model)
   
     CASE (CONSTANT_TCOND)
     kk  = kk_const
        
     CASE (DILUITE_GAS_TCOND)
     kk  = diluite_gas_thermal_cond(T, mu)
        
     CASE (DENSE_GAS_TCOND)
     kk0 = diluite_gas_thermal_cond(T, mu0)
     kk  = dense_gas_thermal_cond(T,v, kk0) 

     CASE (CONST_PRANDTL_TCOND)
     kk  = constant_Prandtl_thermal_cond(T, v, mu)

   END SELECT


   ! Dimensionless viscosity and thermal conductivity
   mu_adim = mu / mu_ref
   kk_adim = kk / kk_ref
   	      
   END SUBROUTINE  transport_coefficients





   SUBROUTINE  transport_coefficients_dim(T_dim, v_dim, mu_dim, kk_dim)
   
   ! The viscosity and the thermal conductivity are computed 
   ! using dimensional temperature and specific volume
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), INTENT(IN)  ::  T_dim,  v_dim
   REAL(KIND=8), INTENT(OUT) :: mu_dim, kk_dim
   
   REAL(KIND=8) :: T, v, mu0, kk0, mu, kk
   !----------------------------------------------------------------------------

   ! Dimensional temperature and specific volume
   T = T_dim
   v = v_dim


   ! Viscosity
   SELECT CASE (viscosity_model)
      
     CASE (CONSTANT_VISC)
     mu  = mu_const
        
     CASE (DILUITE_GAS_VISC)
     mu  = diluite_gas_viscosity(T)

     CASE (DENSE_GAS_VISC)
     mu0 = diluite_gas_viscosity(T)
     mu  = dense_gas_viscosity(T,v,mu0)
     
     CASE (CORRESP_STATE_VISC)
     mu  = corresp_state_viscosity(T)
        
     CASE (SUTHERLAND_LAW_VISC)
     mu  = Sutherland_law_viscosity(T)

     CASE (POWER_LAW_VISC)
     mu  = power_law_viscosity(T)

   END SELECT


   ! Thermal conductivity
   SELECT CASE (thermal_cond_model)

     CASE (CONSTANT_TCOND)
     kk  = kk_const
        
     CASE (DILUITE_GAS_TCOND)
     kk  = diluite_gas_thermal_cond(T, mu)
        
     CASE (DENSE_GAS_TCOND)
     kk0 = diluite_gas_thermal_cond(T, mu0)
     kk  = dense_gas_thermal_cond(T,v, kk0) 

     CASE (CONST_PRANDTL_TCOND)
     kk  = constant_Prandtl_thermal_cond(T, v, mu)
   			       
   END SELECT
   
   mu_dim = mu 
   kk_dim = kk 
   	      
   END SUBROUTINE  transport_coefficients_dim







   !----------------------------------------------------------------------------
   ! VISCOSITY LAWS (DIMENSIONAL)
   !----------------------------------------------------------------------------

   FUNCTION  diluite_gas_viscosity(T)   RESULT(mu)

   ! Uses M_mole, Tc, Fc
   !----------------------------------------------------------------------------   
   IMPLICIT NONE
   
   REAL(KIND=8), INTENT(IN) :: T
   REAL(KIND=8) 	    :: mu

   REAL(KIND=8) :: Ts, Omegas
   ! Coefficients for the computation of the reduced collision integral
   REAL(KIND=8), DIMENSION(10), PARAMETER  ::  CC_O = & 
      (/  1.16145d0,  0.14874d0,  0.52487d0,  0.77320d0,  2.16178d0,  &
   	  2.43787d0, -6.43500d-4, 7.27371d0,  1.80323d0, -0.76830d0  /)
   !----------------------------------------------------------------------------
   
   ! Reference temperature.  Ts = kT/e, k Boltzmann constant 
   ! and e potential energy paramter.  e/k = Tc/1.2593
   Ts = 1.2593d0 * (T/Tc)
   
   ! Reduced collision integral   
   Omegas =  (CC_O(1) / Ts**CC_O(2))	 &
   	  +   CC_O(3) / EXP(CC_O(4)*Ts)  &
   	  +   CC_O(5) / EXP(CC_O(6)*Ts)  &
   	  +   CC_O(7) * (Ts**CC_O(2))	 &
   	  *   SIN( CC_O(9)*Ts**CC_O(10) - CC_O(8))

   mu = 4.0785d-5 * Fc &
      * (SQRT(M_mole*T) / ((vc*1000*M_mole)**(2.d0/3.d0)*Omegas)) 

   END FUNCTION  diluite_gas_viscosity

   
   
          

   FUNCTION  dense_gas_viscosity(T, v, mu0)   RESULT(mu)
   
   ! Uses M_mole, Tc, vc, A   
   !----------------------------------------------------------------------------
   IMPLICIT NONE       
   REAL(KIND=8), INTENT(IN) :: T, v
   REAL(KIND=8) 	    :: mu
   
   REAL(KIND=8) :: Ts, Y, G1, G2 
   REAL(KIND=8) :: mu0, muk, mup
   !----------------------------------------------------------------------------  

   ! Compute the viscosity in the diluite gas limit
   mu0 = diluite_gas_viscosity(T)
   
   ! Reference temperature.  Ts = kT/e, k Boltzmann constant 
   ! and e potential energy paramter.  e/k = Tc/1.2593
   Ts = 1.2593d0 * (T/Tc)
   
   ! G = G(v)
   Y	 = (1.d0/6.d0)*(vc/v)
   G1 = (1.0d0 - 0.5d0*Y)/(1.d0 - Y)**3
   
   G2 = ( AA(1)*( 1.d0 - EXP(-AA(4)*Y))/Y    &
   	+ AA(2)*G1*exp(AA(5)*Y) 	     &
   	+ AA(3)*G1)			     &
        / (AA(1)*AA(4) + AA(2) + AA(3))
   
   ! muk = muk(v),  mup = mup(T,v)
   muk = mu0 * ( 1.d0/G2 + AA(6)*Y)
   mup = 36.344d-6*SQRT(M_mole*Tc)*AA(7)*(Y**2)*G2  &
       * EXP(AA(8) + AA(9)/Ts + AA(10)/Ts**2)	    &
       *(vc*1000*M_mole)**(-2.d0/3.d0)

   mu = 0.1d0*(muk + mup)

   END FUNCTION dense_gas_viscosity

   
   
   

   FUNCTION  corresp_state_viscosity(T)   RESULT(mu)

   ! Uses Pc, Tc
   !----------------------------------------------------------------------------
   IMPLICIT NONE       
   REAL(KIND=8), INTENT(IN)  ::  T
   REAL(KIND=8) 	     ::  mu
   
   REAL(KIND=8)  ::  Tr, Psi 
   !----------------------------------------------------------------------------
   
   ! Reduced temperature
   Tr = T/Tc  
   
   Psi = SQRT(M_mole) * ((Pc/101325.d0)**(2.d0/3.d0)) / Tc**(1.d0/6.d0)	 
   mu  = 3.4d-7 * Psi * Tr**0.94

   END FUNCTION  corresp_state_viscosity

      
      
      

   FUNCTION  Sutherland_law_viscosity(T)   RESULT(mu)

   ! Uses mu_ref_S, T_ref_S, Tchar_S
   !----------------------------------------------------------------------------
   IMPLICIT NONE       
   
   REAL(KIND=8), INTENT(IN) :: T
   REAL(KIND=8) 	    :: mu
   !----------------------------------------------------------------------------   

   mu = mu_ref_S * (T/T_ref_S)**1.5  &
      *  (T_ref_S + Tchar_S) / (T + Tchar_S)
            
   END FUNCTION  Sutherland_law_viscosity

   
   


   FUNCTION  power_law_viscosity(T)   RESULT(mu)

   ! Uses mu_ref_L, T_ref_L, expo
   !----------------------------------------------------------------------------   
   IMPLICIT NONE       
   
   REAL(KIND=8), INTENT(IN) :: T
   REAL(KIND=8) 	    :: mu
   !----------------------------------------------------------------------------      

   mu = mu_ref_L * (T/T_ref_L)**expo
         
   END FUNCTION  power_law_viscosity





   
   
   !----------------------------------------------------------------------------   
   ! THERMAL CONDUCTIVITY LAWS (DIMENSIONAL)
   !----------------------------------------------------------------------------

   FUNCTION  diluite_gas_thermal_cond(T, mu)   RESULT(kk)

   ! Uses M_mole, cv_oo_Tdim, beta
   !----------------------------------------------------------------------------
   USE ideal_specific_heat,   ONLY: cv_oo__Tdim
   
   IMPLICIT NONE
   
   REAL(KIND=8), INTENT(IN) :: T, mu 
   REAL(KIND=8) 	    :: kk
   
   REAL(KIND=8) :: cv_oo, Z, alpha, Psi
   !----------------------------------------------------------------------------   

   cv_oo = cv_oo__Tdim(T)	     
   Z	 = 2.d0 + 10.5d0*(T/Tc)**2
   alpha = cv_oo/Rgas - 1.5d0
     
   Psi   = 1.d0 + alpha &
   	 * ( (0.215d0 + 0.28288d0*alpha - 1.061d0*beta + 0.26665d0*Z)  &
   	   / (0.6366d0 + beta*Z + 1.061d0*alpha*beta))

   kk	= 7.452d0*(mu*Psi/M_mole)*(100.d0/4.184d0)

   END FUNCTION  diluite_gas_thermal_cond

   


   
   FUNCTION  dense_gas_thermal_cond(T, v, kk0)   RESULT(kk)

   ! Uses M_mole, Tc, vc, BB
   !----------------------------------------------------------------------------
   IMPLICIT NONE       
   REAL(KIND=8), INTENT(IN) :: T, v, kk0
   REAL(KIND=8) 	    :: kk
   
   REAL(KIND=8) :: Tr, Y, G1, H2
   REAL(KIND=8) :: kkk, kkp
   !----------------------------------------------------------------------------   

   ! Reduced temperature
   Tr = T/Tc
   
   Y	 = (1.d0/6.d0)*(vc/v)
   G1 = (1.0d0 - 0.5d0*Y)/(1.d0 - Y)**3
   
   H2 = ( BB(1)*(1.d0-EXP(-BB(4)*Y))/Y  &
   	+ BB(2)*G1*EXP(BB(5)*Y)         &
   	+ BB(3)*G1)                     &
        / (BB(1)*BB(4) + BB(2) + BB(3))

   kkk = kk0 * ( 1.d0/H2 + BB(6)*Y)
   kkp = 3.039d-4*SQRT(Tc/M_mole)*BB(7)*(Y**2)*H2*SQRT(Tr)&
       *(vc*1000*M_mole)**(-2.d0/3.d0)*(100.d0/4.184d0)

   kk  = kkk + kkp
   
   END FUNCTION  dense_gas_thermal_cond





   FUNCTION  constant_Prandtl_thermal_cond(T, v, mu)   RESULT(kk)
   !----------------------------------------------------------------------------
   USE thermodynamics,   ONLY: cp__T_v
   
   IMPLICIT NONE       
   REAL(KIND=8), INTENT(IN) :: T, v, mu
   REAL(KIND=8) 	    :: kk   
   !----------------------------------------------------------------------------   
  
   kk = (mu * cp__T_v(T/T_ref,v/v_ref) * cp_ref) / Pr 
  
   END FUNCTION  constant_Prandtl_thermal_cond


   END MODULE  transport_properties
