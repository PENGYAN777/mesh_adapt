! !============================================================ 
! !
! !      Module: van_der_waals_gas
! !
! ! Description: van der Waals pressure equation of state
! !              and (compatible) energy equation of state
! !
! !      Author: Alberto Guardone
! !              Dipartimento di Ingegneria Aerospaziale
! !              Politecnico di Milano
! !              Via La Masa 34, 20156 Milano, ITALY
! !              e-mail: guardone@aero.polimi.it
! !
! !   Copyright: 1998-2003 Alberto Guardone
! !              See COPYING file for copyright notice
! !
! !============================================================ 
! 
! !============================================================ 
! !
! ! Remarks: 
! !
! ! in the module, variables are made dimensionless
! ! by critical value, so that one can use the same form
! ! of the equations of state for both dimensional or 
! ! dimensionless variable, changing the coefficients
! !
! ! Dimensionless thermodynamic variables are defined by
! !
! !                                                   P_c
! ! rho* = rho rho_c ,       P* =  P P_c,   (RT)* = ------- RT ,
! !                                                  rho_c
! !                        
! !                                            T* = T Zc Tc 
! !
! !          P_c                    P_c          
! !   e* = ------- e ,     c*^2 = ------- c^2    
! !         rho_c                  rho_c
! ! 
! ! where the * means variables WITH dimension. 
! !
! ! Dimensionless length, time and momentum are defined by
! !
! !                      ( rho_c )1/2               
! !  x* = L x ,   t* = L (-------)     t , 
! !                      (  P_c  )             
! ! 
! !  m* = (rho_c P_c)^1/2 m 
! !
! ! The coefficients are to be changed as follows
! !
! !  Rgas = 1.d0
! !  a    = 3.d0
! !  b    = 1.d0/3.d0
! !
! !============================================================ 
! 
 MODULE van_der_waals_gas
! 
!    !============================================================ 
!    USE gas_properties,      ONLY: gas
!    USE ideal_specific_heat, ONLY: phi__T, phip__T, phipp__T,   &
!                                   psi__T
!    !============================================================ 
! 
! 
! 
!    !============================================================
!    IMPLICIT NONE
!    
!    PRIVATE
!    ! Gas constant (stored locally)
!    REAL(KIND=8)  ::  Rgas
!    ! van der Waals' contants
!    REAL(KIND=8)  ::  a, b
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    ! Procedures: private/public policy
!    ! ---------------------------------
!    PUBLIC   ::  read_param_vdWG, write_param_vdWG, init_vdWG,      & 
!                 minimum_specific_volume__vdWG,                     &
!                 P__T_v__vdWG, dP__dT_dv__vdWG, d2P__dT2_dv2__vdWG, &
!                 e__T_v__vdWG, de__dT_dv__vdWG, d2e__dT2_dv2__vdWG, &
!                 s__T_v__vdWG, cv__T_v__vdWG, dcv__dT_dv__vdWG
!    !============================================================ 
! 
!    
!  !============================================================   
!  CONTAINS
!  !============================================================ 
! 
! 
! 
! !============================================================ 
! !============================================================ 
! !
! !  Initialization and coefficients
! !
! !============================================================ 
! !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE read_param_vdWG(idf)
!    !============================================================ 
! 
! 
!       !------------------------------------------------------------
!       IMPLICIT NONE
! 
!       INTEGER, INTENT(IN)  ::  idf
!       INTEGER  ::  dummy 
!       !------------------------------------------------------------
! 
!       dummy = idf  ! does nothing
!       
! 
!    END SUBROUTINE read_param_vdWG
!    !============================================================ 
!    
!    
!    
!    !============================================================ 
!    SUBROUTINE write_param_vdWG(idf)
!    !============================================================ 
! 
! 
!       !------------------------------------------------------------
!       IMPLICIT NONE
! 
!       INTEGER, INTENT(IN)  ::  idf
!       INTEGER  ::  dummy 
!       !------------------------------------------------------------
! 
!       dummy = idf  ! does nothing
!       
! 
!    END SUBROUTINE write_param_vdWG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE init_vdWG
!    !============================================================ 
! 
! 
!       IMPLICIT NONE
! 
!       REAL(KIND=8)  :: Zc, Pc, Tc, vc
! 
!       IF( .NOT. gas%vapor) THEN
!          WRITE(*,*) 'The gas ', gas%name,' is not a vapor:'
!          WRITE(*,*) 'van der Waals model is not applicable. '
!          WRITE(*,*) 'in init_vdWG, module van_der_waals_gas.  STOP. '
!          STOP
!       ENDIF
! 
!       ! Retrieve information on the gas (dimensional)   
!       Zc   = 3.d0/8.d0  ! Fixed for van der Waals
!       Rgas = gas%Rgas
!       Pc   = gas%Pc  
!       Tc   = gas%Tc 
!       vc   = (Zc*Rgas*Tc) / Pc
!       
!       ! Store information on the gas structure  
!       gas%Zc    = Zc
!       gas%vc    = vc
!       gas%Z_ref = Zc
! 
!       gas%v_ref = vc
!       gas%T_ref = Zc*Tc
!       gas%uu_ref = SQRT(gas%Pc*gas%v_ref)
!       
!       ! Dimensionless coefficients      
!       Rgas = 1.d0
!          a = 3.d0
!          b = 1.d0/3.d0
! 
!    END SUBROUTINE init_vdWG
!    !============================================================ 
! 
!   
!    
!    !============================================================ 
!    FUNCTION minimum_specific_volume__vdWG() RESULT (v_min)
!    !============================================================ 
! 
!       !------------------------------------------
!       ! Retrieve the minimum specific volume
!       ! allowed for the considered thermodynamic 
!       ! model.  To be used in Newton iterative
!       ! solutions.
!       !------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8)  ::  v_min   ! Minimum specific volume allowed
! 
! 
!       v_min = (1 + 1.d-12)*b
!      
! 
!    END FUNCTION minimum_specific_volume__vdWG
!    !============================================================ 
! 
!   
!   
! !============================================================ 
! !============================================================ 
! !
! !  Thermodynamic quantities as functions of T, v
! !
! !============================================================ 
! !============================================================ 
! 
!   
!    
!    !============================================================ 
!    FUNCTION P__T_v__vdWG(T, v) RESULT (P)
!    !============================================================ 
! 
!       !------------------------------------
!       ! Pressure P(T,v) with T temperature 
!       ! and v specific volume
!       !------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
!       REAL(KIND=8)              ::  P     ! Pressure
! 
!       
!       P = (Rgas*T)/(v - b) - a/v**2
! 
! 
!    END FUNCTION P__T_v__vdWG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION dP__dT_dv__vdWG(T, v) RESULT(dP)
!    !============================================================ 
! 
!       !-------------------------------------------------
!       ! Partial derivative of the pressure P(T,v) with 
!       ! respect to the temperature T and the specific 
!       ! volume v
!       !-------------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)    ::  T, v ! Temperature and specific volume
!       REAL(KIND=8), DIMENSION(2)  ::  dP   ! Partial derivatives of 
!                                            ! the pressure P(T,v)
!       REAL(KIND=8)  ::  P_T, P_v
!                                
! 
!       P_T  =      Rgas    / (v-b)     
!       P_v  =   - (Rgas*T) / (v-b)**2  + (2*a)/v**3
! 
!       dP(1) =  P_T;  dP(2) =  P_v
! 
! 
!    END FUNCTION dP__dT_dv__vdWG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION d2P__dT2_dv2__vdWG(T, v) RESULT(d2P)
!    !============================================================ 
! 
!       !-------------------------------------------------
!       ! Second order partial derivative of the pressure 
!       ! P(T,v) with respect to the temperature T and the 
!       ! specific volume v 
!       !-------------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)      ::  T, v ! Temperature and specific volume
!       REAL(KIND=8), DIMENSION(2,2)  ::  d2P  ! Second partial derivatives of 
!                                              ! the pressure P(T,v) 
!       REAL(KIND=8)  ::  P_TT, P_Tv, P_vv
!                                
! 
!       P_TT =   0.d0
!       P_Tv =   -  Rgas    / (v-b)**2  
!       P_vv =   (2*Rgas*T) / (v-b)**3  - (6*a)/v**4
! 
!       d2P(1,1) = P_TT;      d2P(1,2) = P_Tv
!       d2P(2,1) = d2P(1,2);  d2P(2,2) = P_vv 
! 
! 
!    END FUNCTION d2P__dT2_dv2__vdWG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION e__T_v__vdWG(T, v) RESULT (e)
!    !============================================================ 
! 
!       !-----------------------------------------------
!       ! Specific internal energy per unit mass e(T,v)
!       ! with T temperature and v specific volume
!       !-----------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
!       REAL(KIND=8)              ::  e     ! Internal energy per unit mass
! 
! 
!       e = gas%e_0  +  phi__T(T)  -  a/v
!       
! 
!    END FUNCTION e__T_v__vdWG
!    !============================================================ 
!   
!   
! 
!    !============================================================ 
!    FUNCTION de__dT_dv__vdWG(T, v) RESULT(de)
!    !============================================================ 
! 
!       !---------------------------------------------
!       ! Partial derivative of the specific internal 
!       ! energy per unit mass e(T,v) with respect to 
!       ! the temperature T and the specific volume v
!       !---------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)    ::  T, v ! Temperature and specific volume
!       REAL(KIND=8), DIMENSION(2)  ::  de   ! Partial derivatives of the
!                                            ! energy per unit mass e(T,v)
!       REAL(KIND=8)  ::  e_T, e_v, phip
!                         
! 
!       phip  = phip__T(T)
!            
!       e_T  =   phip  
!       e_v  =   a/v**2
!              
!       de(1) = e_T;  de(2) = e_v
! 
! 
!    END FUNCTION de__dT_dv__vdWG
!    !============================================================
! 
! 
!    
!    !============================================================ 
!    FUNCTION d2e__dT2_dv2__vdWG(T, v) RESULT(d2e)
!    !============================================================ 
! 
!       !----------------------------------------
!       ! Second order partial derivative of the 
!       ! specific internal energy per unit mass 
!       ! e(T,v) with respect to the temperature 
!       ! T and the specific volume v
!       !----------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)      ::  T, v ! Temperature and specific volume
!       REAL(KIND=8), DIMENSION(2,2)  ::  d2e  ! Second partial derivatives of the
!                                              ! energy per unit mass e(T,v)
!       REAL(KIND=8)  ::  e_TT, e_Tv, e_vv, phipp
!                         
! 
!       phipp  = phipp__T(T)
!                       
!       e_TT =   phipp  
!       e_Tv =   0.d0
!       e_vv =   -(2*a)/v**3
!              
!       d2e(1,1) = e_TT;  d2e(1,2) = e_Tv
!       d2e(2,1) = e_Tv;  d2e(2,2) = e_vv
! 
! 
!    END FUNCTION d2e__dT2_dv2__vdWG
!    !============================================================
! 
! 
! 
!    !============================================================ 
!    FUNCTION s__T_v__vdWG(T, v) RESULT (s)
!    !============================================================ 
! 
!       !--------------------------------------------
!       ! Specific entropy per unit mass s(T,v) with
!       ! T temperature and v specific volume
!       !--------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)  ::  T,v  ! Temperature and specific volume
!       REAL(KIND=8)              ::  s    ! Specific entropy per unit mass   
! 
! 
!       s = gas%s_0 + Rgas*LOG(v-b) + psi__T(T) 
!       
! 
!    END FUNCTION s__T_v__vdWG
!    !============================================================ 
! 
!  
!  
! 
!    FUNCTION  cv__T_v__vdWG(T, v)   RESULT(cv)
! 
!    ! Specific heat at constant volume cv(T,v),
!    ! with T temperature and v specific volume
!    !------------------------------------------------------------------------
!    IMPLICIT NONE 
! 
!    REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
!    REAL(KIND=8)              ::  cv    ! Specific heat at constant volume   
! 
!    REAL(KIND=8) :: void
!    !------------------------------------------------------------------------
! 
!    ! Dummy assignment to avoid compilation warnings 
!    ! due to the unnecessity of v in ideal gas eos
!    ! for internal energy
!    void = v
! 
!    cv = phip__T(T)
! 
!    END FUNCTION cv__T_v__vdWG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION  dcv__dT_dv__vdWG(T, v)   RESULT(dcv)
! 
!    ! Partial derivatives of the specific heat 
!    ! at constant volume cv(T,v) with respect
!    ! to the temperature T and the specific
!    ! volume v
!    !------------------------------------------------------------------------
!    IMPLICIT NONE 
! 
!    REAL(KIND=8), INTENT(IN)    ::  T, v  ! Temperature and specific volume
!    REAL(KIND=8), DIMENSION(2)  ::  dcv   ! Partial derivative of the 
!                                          ! specific heat at constant volume   
!    REAL(KIND=8)  ::  cv_T, cv_v
!    
!    REAL(KIND=8) :: void
!    !------------------------------------------------------------------------
! 
!    ! Dummy assignment to avoid compilation warnings 
!    ! due to the unnecessity of v in ideal gas eos
!    ! for internal energy
!    void = v
! 
!    cv_T = phipp__T(T)  
!    cv_v = 0.d0
! 
!    dcv(1) = cv_T;  dcv(2) = cv_v
!    
!    END FUNCTION dcv__dT_dv__vdWG
!    !============================================================ 
! 
!     
! 
 END MODULE van_der_waals_gas
