! !============================================================ 
! !
! !      Module: clausius_II_gas
! !
! ! Description: Clausius II pressure equation of state and 
! !              (compatible) energy equation of state
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
! !============================================================ 
! 
 MODULE clausius_II_gas
! 
! 
!    !============================================================ 
!    USE gas_properties,      ONLY: gas
!    USE ideal_specific_heat, ONLY: phi__T, phip__T, phipp__T,      &
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
!    ! Clausius-II coefficients
!    REAL(KIND=8)  ::  a, b,c 
!    ! Clausius-II functions Q(T) and F(v)
!    REAL(KIND=8), DIMENSION(0:3)  ::  QQ, FF
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    ! Procedures: private/public policy
!    ! ---------------------------------
!    PRIVATE  ::  Q_F__T_v__CIIG
!    PUBLIC   ::  read_param_CIIG, write_param_CIIG, init_CIIG,      & 
!                 minimum_specific_volume__CIIG,                     &
!                 P__T_v__CIIG, dP__dT_dv__CIIG, d2P__dT2_dv2__CIIG, &
!                 e__T_v__CIIG, de__dT_dv__CIIG, d2e__dT2_dv2__CIIG, &
!                 s__T_v__CIIG, cv__T_v__CIIG, dcv__dT_dv__CIIG
!    !============================================================ 
!                 
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
!    !============================================================ 
!    SUBROUTINE read_param_CIIG(idf)
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
!    END SUBROUTINE read_param_CIIG
!    !============================================================ 
!    
!    
!    
!    !============================================================ 
!    SUBROUTINE write_param_CIIG(idf)
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
!    END SUBROUTINE write_param_CIIG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE init_CIIG()
!    !============================================================ 
! 
! 
!       IMPLICIT NONE
! 
!       REAL(KIND=8)  ::  Zc, Pc, Tc, vc
!       
! 
!       IF( .NOT. gas%vapor) THEN
!          WRITE(*,*) 'The gas ', gas%name,' is not a vapor:'
!          WRITE(*,*) 'Clausius_II model is not applicable. '
!          WRITE(*,*) 'in init_CIIG, module clausius_II_gas.  STOP. '
!          STOP
!       ENDIF
!             
!       ! Retrieve information on the gas (dimensional)   
!       Rgas  = gas%Rgas
!       Pc    = gas%Pc  
!       Tc    = gas%Tc 
!       Zc    = gas%Zc
!       vc   = (Zc*Rgas*Tc) / Pc
! 
!       ! Store information on the gas structure  
!       gas%vc    = vc
!       gas%Z_ref = Zc
!       gas%v_ref = vc
!       gas%T_ref = Zc*Tc
!       gas%uu_ref = SQRT(gas%Pc*gas%v_ref)
! 
!       ! Clausius-II's coefficients
!       a = (27.d0/64.d0)*((Rgas*(Tc**2)*vc)/Zc)
!       b = vc*(1.d0-1/(4*Zc))
!       c = vc*(-1.d0 + 3.d0/(Zc*8.d0))
!       
!       ! Coefficients are made dimensionless and scaled for
!       ! using Zc Tc instead of Tc as reference temperature.
!       Rgas = 1.d0
!          a = a/(Pc*Zc*Tc*vc**2)
!          b = b/vc      
!          c = c/vc      
!       !------------------------------------------------------------
! 
!    END SUBROUTINE init_CIIG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE Q_F__T_v__CIIG(T, v,  QQ, FF)
!    !============================================================ 
! 
!       !------------------------------------------------------------
!       ! Computes the functions QQ(T), T temperature, and FF(v), v
!       ! specific volume, and their derivatives up to the third 
!       ! order.
!       ! Output matrices are QQ(k) and FF(k) k is the order of the 
!       ! derivative.
!       !------------------------------------------------------------
! 
!       IMPLICIT NONE 
!                                   
!       REAL(KIND=8),                 INTENT(IN)   ::  T, v  
!       REAL(KIND=8), DIMENSION(0:3), INTENT(OUT)  ::  QQ, FF
! 
!       
!       ! Q(T): Q, Q', Q'', Q'''
!       QQ(0) =      a/T
!       QQ(1) =     -a/(T**2)
!       QQ(2) =  (2*a)/(T**3)
!       QQ(3) = -(6*a)/(T**4)
! 
!       ! F(T): F, F', F'', F'''
!       FF(0) =   1/(v+c)  
!       FF(1) =  -1/(v+c)**2
!       FF(2) =   2/(v+c)**3   
!       FF(3) =  -6/(v+c)**4  
!    
!    
!    END SUBROUTINE Q_F__T_v__CIIG
!    !============================================================ 
!  
!  
!    
!    !============================================================ 
!    FUNCTION minimum_specific_volume__CIIG() RESULT (v_min)
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
!    END FUNCTION minimum_specific_volume__CIIG
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
!    FUNCTION P__T_v__CIIG(T, v) RESULT (P)
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
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
! 
!       P = (Rgas*T)/(v - b) + QQ(0)*FF(1) 
! 
! 
!    END FUNCTION P__T_v__CIIG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION dP__dT_dv__CIIG(T, v) RESULT(dP)
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
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
! 
!       P_T  =      Rgas    / (v-b)     +  QQ(1)*FF(1) 
!       P_v  =   - (Rgas*T) / (v-b)**2  +  QQ(0)*FF(2) 
! 
!       dP(1) =  P_T;  dP(2) =  P_v
! 
! 
!    END FUNCTION dP__dT_dv__CIIG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION d2P__dT2_dv2__CIIG(T, v) RESULT(d2P)
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
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
! 
!       P_TT =                             QQ(2)*FF(1) 
!       P_Tv =   -  Rgas    / (v-b)**2  +  QQ(1)*FF(2) 
!       P_vv =   (2*Rgas*T) / (v-b)**3  +  QQ(0)*FF(3) 
! 
!       d2P(1,1) = P_TT;      d2P(1,2) = P_Tv
!       d2P(2,1) = d2P(1,2);  d2P(2,2) = P_vv 
! 
! 
!    END FUNCTION d2P__dT2_dv2__CIIG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION e__T_v__CIIG(T, v) RESULT (e)
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
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
!       
!       e = gas%e_0  +  phi__T(T)  +  (T*QQ(1) - QQ(0))*FF(0) 
!       
! 
!    END FUNCTION e__T_v__CIIG
!    !============================================================ 
!   
!   
! 
!    !============================================================ 
!    FUNCTION de__dT_dv__CIIG(T, v) RESULT(de)
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
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
!       
!       phip  = phip__T(T)
!            
!       e_T  =   phip  +   T*QQ(2)        * FF(0) 
!       e_v  =            (T*QQ(1)-QQ(0)) * FF(1) 
!              
!       de(1) = e_T;  de(2) = e_v
! 
! 
!    END FUNCTION de__dT_dv__CIIG
!    !============================================================
! 
! 
!    
!    !============================================================ 
!    FUNCTION d2e__dT2_dv2__CIIG(T, v) RESULT(d2e)
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
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
!       
!       phipp  = phipp__T(T)
!                       
!       e_TT =   phipp  +  (T*QQ(3)+QQ(2)) * FF(0) 
!       e_Tv =              T*QQ(2)        * FF(1) 
!       e_vv =             (T*QQ(1)-QQ(0)) * FF(2) 
!              
!       d2e(1,1) = e_TT;  d2e(1,2) = e_Tv
!       d2e(2,1) = e_Tv;  d2e(2,2) = e_vv
! 
! 
!    END FUNCTION d2e__dT2_dv2__CIIG
!    !============================================================
! 
! 
! 
!    !============================================================ 
!    FUNCTION s__T_v__CIIG(T, v) RESULT (s)
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
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF) 
!      
!       s = gas%s_0 + Rgas*LOG(v-b) + psi__T(T)  +  QQ(1)*FF(0) 
! 
! 
!    END FUNCTION s__T_v__CIIG
!    !============================================================ 
! 
!  
!  
!    !============================================================ 
!    FUNCTION cv__T_v__CIIG(T, v) RESULT (cv)
!    !============================================================ 
! 
!       !------------------------------------------
!       ! Specific heat at constant volume cv(T,v),
!       ! with T temperature and v specific volume
!       !------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)  ::  T, v  ! Temperature and specific volume
!       REAL(KIND=8)              ::  cv    ! Specific heat at constant volume   
! 
!       
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
! 
!       cv = phip__T(T)    +  T*QQ(2)*FF(0)
! 
! 
!    END FUNCTION cv__T_v__CIIG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION dcv__dT_dv__CIIG(T, v) RESULT (dcv)
!    !============================================================ 
! 
!       !------------------------------------------
!       ! Partial derivatives of the specific heat 
!       ! at constant volume cv(T,v) with respect
!       ! to the temperature T and the specific
!       ! volume v
!       !------------------------------------------
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), INTENT(IN)    ::  T, v  ! Temperature and specific volume
!       REAL(KIND=8), DIMENSION(2)  ::  dcv   ! Partial derivative of the 
!                                             ! specific heat at constant volume   
!       REAL(KIND=8)  ::  cv_T, cv_v
!       
!       
!       CALL Q_F__T_v__CIIG(T, v,  QQ, FF)
! 
!       cv_T = phipp__T(T)  +  (T*QQ(3) + QQ(2))*FF(0) 
!       cv_v =                  T*QQ(2)         *FF(1) 
! 
!       dcv(1) = cv_T;  dcv(2) = cv_v
!       
! 
!    END FUNCTION dcv__dT_dv__CIIG
!    !============================================================ 
! 
! 
! 
 END MODULE clausius_II_gas
