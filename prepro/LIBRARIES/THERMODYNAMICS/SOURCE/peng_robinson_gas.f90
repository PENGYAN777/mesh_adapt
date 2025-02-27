! 
! !============================================================ 
! !
! !      Module: peng_robinson_gas
! !
! ! Description: Peng-Robinson pressure equation of state and 
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
 MODULE peng_robinson_gas
! 
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
!    ! Peng-Robinson coefficients
!    REAL(KIND=8)  ::  b, fw, C_a
!    REAL(KIND=8)  ::  Zc ! Needed to adimensionalize 
!                         ! the temperature
!    ! Peng-Robinson functions Q(T) and F(v)
!    REAL(KIND=8), DIMENSION(0:3)  ::  QQ, FF
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    ! Procedures: private/public policy
!    ! ---------------------------------
!    PRIVATE  ::  Q_F__T_v__PRG
!    PUBLIC   ::  read_param_PRG, write_param_PRG, init_PRG,      & 
!                 minimum_specific_volume__PRG,                   &
!                 P__T_v__PRG, dP__dT_dv__PRG, d2P__dT2_dv2__PRG, &
!                 e__T_v__PRG, de__dT_dv__PRG, d2e__dT2_dv2__PRG, &
!                 s__T_v__PRG, cv__T_v__PRG, dcv__dT_dv__PRG
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
!    SUBROUTINE read_param_PRG(idf)
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
!    END SUBROUTINE read_param_PRG
!    !============================================================ 
!    
!    
!    
!    !============================================================ 
!    SUBROUTINE write_param_PRG(idf)
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
!    END SUBROUTINE write_param_PRG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE init_PRG
!    !============================================================ 
! 
! 
!       IMPLICIT NONE
! 
!       REAL(KIND=8)  ::  Pc, Tc, vc, omega
! 
!       IF( .NOT. gas%vapor) THEN
!          WRITE(*,*) 'The gas ', gas%name,' is not a vapor:'
!          WRITE(*,*) 'Peng-Robinson model is not applicable. '
!          WRITE(*,*) 'in init_PRG, module peng_robinson_gas.  STOP. '
!          STOP
!       ENDIF
!             
!       ! Retrieve information on the gas (dimensional)   
!       ! Zc    = 0.32137894290237984 ! Fixed for PR
!       Rgas  = gas%Rgas
!       Pc    = gas%Pc  
!       Tc    = gas%Tc 
!       omega = gas%omega
! 
!       ! Peng-Robinson's coefficients
!       C_a = 0.45724d0 * ((Rgas*Tc)**2)/Pc
!         b = 0.07780d0 *  (Rgas*Tc)    /Pc
!        fw = 0.37464d0  +  1.54226d0*omega  -  0.26992d0*omega**2
! 
!       ! Critical specific volume and compressibility factor 
!       vc = v__P_T__PRG(Pc, Tc, Tc)
!       Zc = (Pc*vc) / (Rgas*Tc) ! Fixed for PRG, is 0.321379025173
! 
!       ! Store information on the gas structure  
!       gas%Zc    = Zc
!       gas%vc    = vc
!       gas%Z_ref = Zc
!       gas%v_ref = vc
!       gas%T_ref = Zc*Tc
!       gas%uu_ref = SQRT(gas%Pc*gas%v_ref)
!       
!       ! Dimensionless coefficients are scaled for
!       ! using Zc Tc instead of Tc as reference temperature.
!       Rgas = 1.d0
!          b = b/vc      
!        C_a = C_a/(Pc*vc**2)
!       !------------------------------------------------------------
!       
!       CONTAINS
!       
!       
!       !============================================================ 
!       FUNCTION v__P_T__PRG(P, T, Tc) RESULT (v)
!       !============================================================ 
! 
!          ! Specific volume v(P,T) as a function of the
!          ! pressure P and temperature T. Dimensional 
!          ! form of the equation for computing Zc.
!          ! For PR, dimensional and dimensionless equation
!          ! have DIFFERENT form.
! 
!          IMPLICIT NONE 
! 
!          REAL(KIND=8), INTENT(IN)  ::  P, T, Tc
!          REAL(KIND=8)              ::  v
! 
!          REAL(KIND=8)  ::  a, f, fp
!          INTEGER  ::  it
!          REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-10
!          INTEGER,      PARAMETER  ::  NEWTON_MAX_ITER = 100
! 
!          !Initial guess: polytropic ideal gas 
!          v = (Rgas*T)/P  
!          
!          DO it = 1, NEWTON_MAX_ITER
!      
!             a =   C_a * (1.d0 + fw*(1.d0 - SQRT(T/Tc)))
! 
!             f =  (Rgas*T)/(v-b) - a/(v**2 + 2*v*b - b**2) - P
!             
!             IF ( ABS(f) < NEWTON_EPS ) RETURN 
! 
!             fp =  -(Rgas*T)/(v-b)**2 &
!                + (2*a*(v+b))/(v**2 + 2*v*b - b**2)**2
!             
!             v = v  -  f / fp
! 
!          ENDDO
! 
!          WRITE(*,*) ' ERROR:  v__P_T__RKG. STOP'
!          STOP
! 
! 
!       END FUNCTION v__P_T__PRG
!       !============================================================ 
!       
!    END SUBROUTINE init_PRG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE Q_F__T_v__PRG(T, v,  QQ, FF)
!    !============================================================ 
! 
!       !------------------------------------------------------------
!       ! Computes the functions QQ(T), T temperature, and FF(v), v
!       ! specific volume, and their derivatives up to the third 
!       ! order.
!       ! Output matrices are QQ(k) and FF(k) where k is the 
!       ! order of the derivative.
!       !------------------------------------------------------------
! 
!       IMPLICIT NONE 
!                                   
!       REAL(KIND=8),                 INTENT(IN)   ::  T, v  
!       REAL(KIND=8), DIMENSION(0:3), INTENT(OUT)  ::  QQ, FF
!       REAL(KIND=8)  ::  beta, sq2 
! 
! 
!       sq2 = SQRT(2.d0)
!       beta = 1.d0 + fw*(1.d0 - SQRT(Zc*T))
!       
!       ! Q(T): Q, Q', Q'', Q'''
!       QQ(0) = C_a * (beta**2.d0) 
!       QQ(1) = - C_a*fw * SQRT(Zc/T) * beta 
!       QQ(2) = ((C_a*fw)/(2*T))   *(SQRT(Zc/T)*beta + fw*Zc)
!       QQ(3) = - ((3.d0*C_a*fw)/(4.d0*T**2))*(SQRT(Zc/T)*beta + fw*Zc)
! 
!       ! F(T): F, F', F'', F'''
!       FF(0) =   LOG( (v + b*(1+sq2))/(v + b*(1-sq2)) ) / (2*b*sq2) 
!       FF(1) =       -1/(v**2 + 2*v*b - b**2)  
!       FF(2) =  2*(v+b)/(v**2 + 2*v*b - b**2)**2     
!       FF(3) =  -2*(3*v**2 + 6*v*b + 5*b**2)/(v**2 + 2*v*b - b**2)**3
!    
!    
!    END SUBROUTINE Q_F__T_v__PRG
!    !============================================================ 
!  
!  
!    
!    !============================================================ 
!    FUNCTION minimum_specific_volume__PRG() RESULT (v_min)
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
!    END FUNCTION minimum_specific_volume__PRG
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
!    FUNCTION P__T_v__PRG(T, v) RESULT (P)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
! 
!       P = (Rgas*T)/(v - b) + QQ(0)*FF(1) 
! 
! 
!    END FUNCTION P__T_v__PRG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION dP__dT_dv__PRG(T, v) RESULT(dP)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
! 
!       P_T  =      Rgas    / (v-b)     +  QQ(1)*FF(1) 
!       P_v  =   - (Rgas*T) / (v-b)**2  +  QQ(0)*FF(2) 
! 
!       dP(1) =  P_T;  dP(2) =  P_v
! 
! 
!    END FUNCTION dP__dT_dv__PRG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION d2P__dT2_dv2__PRG(T, v) RESULT(d2P)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
! 
!       P_TT =                             QQ(2)*FF(1) 
!       P_Tv =   -  Rgas    / (v-b)**2  +  QQ(1)*FF(2) 
!       P_vv =   (2*Rgas*T) / (v-b)**3  +  QQ(0)*FF(3) 
! 
!       d2P(1,1) = P_TT;      d2P(1,2) = P_Tv
!       d2P(2,1) = d2P(1,2);  d2P(2,2) = P_vv 
! 
! 
!    END FUNCTION d2P__dT2_dv2__PRG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION e__T_v__PRG(T, v) RESULT (e)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
!       
!       e = gas%e_0  +  phi__T(T)  +  (T*QQ(1) - QQ(0))*FF(0) 
!       
! 
!    END FUNCTION e__T_v__PRG
!    !============================================================ 
!   
!   
! 
!    !============================================================ 
!    FUNCTION de__dT_dv__PRG(T, v) RESULT(de)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
!       
!       phip  = phip__T(T)
!            
!       e_T  =   phip  +   T*QQ(2)        * FF(0) 
!       e_v  =            (T*QQ(1)-QQ(0)) * FF(1) 
!              
!       de(1) = e_T;  de(2) = e_v
! 
! 
!    END FUNCTION de__dT_dv__PRG
!    !============================================================
! 
! 
!    
!    !============================================================ 
!    FUNCTION d2e__dT2_dv2__PRG(T, v) RESULT(d2e)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
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
!    END FUNCTION d2e__dT2_dv2__PRG
!    !============================================================
! 
! 
! 
!    !============================================================ 
!    FUNCTION s__T_v__PRG(T, v) RESULT (s)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF) 
!      
!       s = gas%s_0 + Rgas*LOG(v-b) + psi__T(T)  +  QQ(1)*FF(0) 
! 
! 
!    END FUNCTION s__T_v__PRG
!    !============================================================ 
! 
!  
!  
!    !============================================================ 
!    FUNCTION cv__T_v__PRG(T, v) RESULT (cv)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
! 
!       cv = phip__T(T)    +  T*QQ(2)*FF(0)
! 
! 
!    END FUNCTION cv__T_v__PRG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION dcv__dT_dv__PRG(T, v) RESULT (dcv)
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
!       CALL Q_F__T_v__PRG(T, v,  QQ, FF)
! 
!       cv_T = phipp__T(T)  +  (T*QQ(3) + QQ(2))*FF(0) 
!       cv_v =                  T*QQ(2)         *FF(1) 
! 
!       dcv(1) = cv_T;  dcv(2) = cv_v
!       
! 
!    END FUNCTION dcv__dT_dv__PRG
!    !============================================================ 
! 
! 
! 
 END MODULE peng_robinson_gas
