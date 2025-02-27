! !============================================================ 
! !
! !      Module: martin_hou_gas
! !
! ! Description: Martin-Hou pressure equation of state and 
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
! ! The coefficients are to be changed as follows
! !
! !  Rgas = 1.d0
! !  kk   = kk*Zc
! !  BB   = BB*Zc
! !
! ! where on the right hand side the coefficient are evaluated
! ! using the standard reduced variables.
! !
! !============================================================ 
! 
 MODULE martin_hou_gas
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
!    ! Martin-Hou coefficients
!    REAL(KIND=8), DIMENSION(2:5)  ::  AA, BB, CC 
!    REAL(KIND=8)  ::  b, kk
!    ! Martin-Hou functions Q(T) and F(v)
!    REAL(KIND=8), DIMENSION(2:5,0:3)  ::  QQ, FF
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    ! Procedures: private/public policy
!    ! ---------------------------------
!    PRIVATE  ::  compute_coeff_adim, compute_coeff_dim,          &
!                 Q_F__T_v__MHG
!    PUBLIC   ::  read_param_MHG, write_param_MHG, init_MHG,      & 
!                 minimum_specific_volume__MHG,                   &
!                 P__T_v__MHG, dP__dT_dv__MHG, d2P__dT2_dv2__MHG, &
!                 e__T_v__MHG, de__dT_dv__MHG, d2e__dT2_dv2__MHG, &
!                 s__T_v__MHG, cv__T_v__MHG, dcv__dT_dv__MHG
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
!    !============================================================ 
!    SUBROUTINE read_param_MHG(idf)
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
!    END SUBROUTINE read_param_MHG
!    !============================================================ 
!    
!    
!    
!    !============================================================ 
!    SUBROUTINE write_param_MHG(idf)
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
!    END SUBROUTINE write_param_MHG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE init_MHG
!    !============================================================ 
! 
! 
!       IMPLICIT NONE
! 
!       REAL(KIND=8)  ::  Zc, Pc, Tc, vc, Tb
! 
!       IF( .NOT. gas%vapor) THEN
!          WRITE(*,*) 'The gas ', gas%name,' is not a vapor:'
!          WRITE(*,*) 'Martin-Hou model is not applicable. '
!          WRITE(*,*) 'in init_MHG, module martin_hou_gas.  STOP. '
!          STOP
!       ENDIF
!             
!       ! Retrieve information on the gas (dimensional)   
!       Rgas = gas%Rgas
!       Pc   = gas%Pc  
!       Tc   = gas%Tc 
!       Zc   = gas%Zc 
!       Tb   = gas%Tb
!       vc   = (Zc*Rgas*Tc) / Pc
!       
!       ! Store information on the gas structure  
!       gas%vc    = vc
!       gas%v_ref = vc
!       gas%uu_ref = SQRT(Pc*vc)
!       
!       !------------------------------------------------------------
!       ! Computes Martin-Hou's coefficients
! 	    CALL compute_coeff_adim(Zc, Pc, Tc, Tb)
!       ! CALL compute_coeff_dim(Zc, Pc, Tc, Tb)  ! <<< ALTERNATIVE
!       
!       ! Dimensionless coefficients are scaled for
!       ! using Zc Tc instead of Tc as reference temperature.
!       BB   = BB * Zc    
!       kk   = 5.475d0*Zc
!       Rgas = 1.d0
!       !------------------------------------------------------------
!       
!    END SUBROUTINE init_MHG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE compute_coeff_adim(Zc, Pc, Tc, Tb)
!    !============================================================ 
!      
!      
!       !------------------------------------------------------------
!       ! Computes the coeffiecients of the Martin-Hou equation of
!       ! state using DIMENSIONLESS variables.  Output coefficient 
!       ! are made dimensionless by critical values.
!       !------------------------------------------------------------
!      
!      
!       !------------------------------------------------------------
!       IMPLICIT NONE
!    
!       REAL(KIND=8), INTENT(IN)  ::  Zc, Pc, Tc, Tb
!       
!       REAL(KIND=8)  :: rZc, Tb_, Pb_, TBoyle_, Tp_, xi, m, beta
!       REAL(KIND=8)  :: f2, f3, f4, f5
!       REAL(KIND=8), PARAMETER  ::  n = 2.d0
!       !------------------------------------------------------------
!       
!       
!       kk = 5.475d0
!       ! Shorthand
!       rZc   = 1/Zc
!       ! Covolume constant
!       beta  = -31.883d0*Zc**2 + 20.533d0*Zc
!       b     = 1.d0  -  beta / (15.d0*Zc)
!       ! Temperature at Boyle point
!       TBoyle_  = 9.654921d0/Tc + 2.505308d0 - 6.492324d-4*Tc
!       ! T'
!       Tp_   = -0.6751d0*Zc + 0.9869d0      
!       ! Slope of reduced vapor pressure curve at critical point
!       ! Use the reduced pressure (Pb_) and reduced temperature 
!       ! (Tb_) at boiling point, and it is equal to the Reidel 
!       ! parameter, namely, alpha_c  ==  -m 
!       Tb_   =  Tb / Tc
!       Tp_   = -0.6751d0*Zc + 0.9869d0
!       Pb_   = 101325.d0/ Pc  ! Standard boiling pressure is 1 atm
!       xi    = -35.d0 + 36.d0/Tb_ + 42.d0*log(Tb_) - (Tb_)**6.d0
!       m     = (0.315d0*xi - log(Pb_)) / (0.0838d0*xi - log(Tb_))
! 
!       ! f functions
!       f2 =   - 3.8d0*rZc*(1.d0-b)        +   9.d0*(1.d0-b)**2.d0
!       f3 =     5.4d0*rZc*(1.d0-b)**2.d0  -  17.d0*(1.d0-b)**3.d0
!       f4 =   - 3.4d0*rZc*(1.d0-b)**3.d0  +  12.d0*(1.d0-b)**4.d0 
!       f5 =     0.8d0*rZc*(1.d0-b)**4.d0  -   3.d0*(1.d0-b)**5.d0
!  
!       ! Coefficients AA, BB, CC
!       AA = 0.d0; BB = 0.d0; CC = 0.d0
! 
!       CC(2) =  (   (f2 + b*rZc*Tp_ + ((rZc*Tp_)**2)*(1.d0-Zc)) * (TBoyle_ - 1.d0) &
!                  - (f2 + b*rZc*TBoyle_)                        * (Tp_     - 1.d0) ) &
!             /  (   (EXP(-kk) - EXP(-kk*Tp_ ))                  * (TBoyle_ - 1.d0) &
!                  - (EXP(-kk) - EXP(-kk*TBoyle_))               * (Tp_     - 1.d0) )
!       
!       BB(2) =  (- f2 - b*rZc*TBoyle_ + CC(2)*(EXP(-kk) - EXP(-kk*TBoyle_)) ) &
!             /  (TBoyle_ - 1.d0)
! 
!       CC(3) =  CC(2) * ((1.d0 - b)**3 - (1.d0/n - b)**3) &
!             /  ((1.d0/n - b)**2 - (1.d0 - b)**2)
! 
!       AA(2) =  f2  -  BB(2)  -  CC(2) * EXP(-kk)
! 
!       AA(4) =  f4
!  
!       CC(5) = -CC(2)*(1.d0-b)**3  -  CC(3)*(1.d0-b)**2  
!  
!       BB(5) = f5 - CC(5)*EXP(-kk)
!  
!       BB(3) =  m     * (1.d0-b)**3  -  rZc   * (1.d0-b)**2  &
!             -  BB(2) * (1.d0-b)     -  BB(5) / (1.d0-b)**2
!             
!       AA(3) =  f3  -  BB(3) -  CC(3)*EXP(-kk)
!       
! 
!    END SUBROUTINE compute_coeff_adim
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE compute_coeff_dim(Zc, Pc, Tc, Tb)
!    !============================================================ 
!      
!       !------------------------------------------------------------
!       ! Computes the coeffiecients of the Martin-Hou equation of
!       ! state using DIMENSIONAL variables.  Output coefficient 
!       ! are made dimensionless by critical values.
!       !------------------------------------------------------------
! 
!       !------------------------------------------------------------
!       IMPLICIT NONE
!       
!       REAL(KIND=8), INTENT(IN)  ::  Zc, Pc, Tc, Tb
! 
!       REAL(KIND=8)  :: vc, Tb_, Pb_, TBoyle_, Tp_, xi, m, beta, &
!                        RTc, vc_b, Tp, TBoyle, MM
!       REAL(KIND=8)  :: f2, f3, f4, f5
!       REAL(KIND=8), PARAMETER  ::  n = 2.d0
!       INTEGER  ::  i
!       !------------------------------------------------------------
!  
!       
!       kk = 5.475d0
!       ! Shorthand
!       RTc = Rgas*Tc      
!       vc   = (Zc*RTc) / Pc
!       ! Covolume constant
!       beta = -31.883d0 * Zc**2  +  20.533d0 * Zc
!       b  = vc*(1 - beta/(15.d0*Zc))
!       vc_b = vc - b  ! Shorthand for vc - b
!       ! Temperature at Boyle point      
!       TBoyle  = 9.654291d0  +  2.505308 * Tc  -  6.492324d-4 * Tc**2 
!       TBoyle_ = TBoyle / Tc
!       ! T'
!       Tp  = (-0.6751d0 * Zc  +  0.9869d0) * Tc  
!       Tp_ = Tp / Tc
!       ! Slope of reduced vapor pressure curve at critical point
!       ! Use the reduced pressure (Pb_) and reduced Temperature 
!       ! (Tb_) at boiling point, and it is equal to the Reidel 
!       ! parameter, namely, alpha_c  ==  -m 
!       Tb_ = Tb / Tc
!       Pb_   = 101325.d0/ Pc ! Standard boiling pressure is 1 atm
!       xi  = - 35.d0  +  36.d0 / Tb_  +  42.d0 * LOG(Tb_)  -  Tb_**6.d0
!       MM  = (0.315d0 * xi  -  LOG(Pb_)) / (0.0838d0 * xi  -  LOG(Tb_)) 
!       m = MM*(Pc/Tc)
! 
!       ! f functions      
!       f2 = - 3.8d0 * RTc * vc_b        +   9.d0 * Pc * vc_b**2.d0 
!       f3 =   5.4d0 * RTc * vc_b**2.d0  -  17.d0 * Pc * vc_b**3.d0 
!       f4 = - 3.4d0 * RTc * vc_b**3.d0  +  12.d0 * Pc * vc_b**4.d0 
!       f5 =   0.8d0 * RTc * vc_b**4.d0  -   3.d0 * Pc * vc_b**5.d0 
!       
!       ! Coefficients AA, BB, CC 
!       AA = 0; BB = 0; CC = 0
!       
!       CC(2) =  (   (f2 + b*Rgas*Tp + ((Rgas*Tp)**2)*(1.d0-Zc)/Pc) * (TBoyle - Tc) &
!                  - (f2 + b*Rgas*TBoyle)                           * (Tp     - Tc) ) &
!             /  (   (EXP(-kk) - EXP(-kk*Tp_ ))                     * (TBoyle - Tc) &
!                  - (EXP(-kk) - EXP(-kk*TBoyle_))                  * (Tp     - Tc) )
!       
!       BB(2) =  (- f2 - b*Rgas*TBoyle + CC(2)*(EXP(-kk) - EXP(-kk*TBoyle_)) ) &
!             /  (TBoyle - Tc)
! 
!       CC(3) =  (CC(2) * (vc_b**3 - (vc/n - b)**3)) &
!             /  ((vc/n - b)**2-vc_b**2)
! 
!       AA(2) =  f2  -  BB(2) * Tc  - CC(2) * EXP(-kk)
! 
!       AA(4) =  f4
!  
!       CC(5) = -CC(2)*(vc-b)**3  -  CC(3)*(vc-b)**2  
!  
!       BB(5) = (f5 - CC(5)*EXP(-kk)) / Tc
!  
!       BB(3) =  m     * vc_b**3  -  Rgas  * vc_b**2  &
!             -  BB(2) * vc_b     -  BB(5) / vc_b**2
!             
!       AA(3) =  f3  -  BB(3)*Tc  -  CC(3)*EXP(-kk)
! 
!       ! Coefficient are made dimensionless by critical values
!       DO i = 2, 5
!          AA(i) = AA(i) / ( vc**i * Pc      )
!          BB(i) = BB(i) / ((vc**i * Pc) / Tc)
!          CC(i) = CC(i) / ( vc**i * Pc      )
!       ENDDO
!       b = b / vc
! 
! 
!    END SUBROUTINE compute_coeff_dim
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE Q_F__T_v__MHG(T, v,  QQ, FF)
!    !============================================================ 
! 
!       !------------------------------------------------------------
!       ! Computes the functions QQ(T), T temperature, and FF(v), v
!       ! specific volume, and their derivatives up to the third 
!       ! order.
!       ! Output matrices are QQ(i,k) and FF(i,k) where i is the 
!       ! index of the function and k the order of the derivative.
!       !------------------------------------------------------------
! 
!       IMPLICIT NONE 
!                                   
!       REAL(KIND=8),                      INTENT(IN)  ::  T, v  
!       REAL(KIND=8), DIMENSION(2:5,0:3), INTENT(OUT)  ::  QQ, FF
!       INTEGER  ::  i
! 
! 
!       DO i = 2, 5
! 
!          ! Q(T): Q, Q', Q'', Q'''
!          QQ(i,0) =  AA(i)  +  BB(i)*T  +  CC(i)*        EXP(-kk*T)
!          QQ(i,1) =            BB(i)    -  CC(i)* kk    *EXP(-kk*T)
!          QQ(i,2) =                        CC(i)*(kk**2)*EXP(-kk*T)
!          QQ(i,3) =                     -  CC(i)*(kk**3)*EXP(-kk*T)
! 
!          ! F(T): F, F', F'', F'''
!          FF(i,0) =  - (1.d0/(i-1))    / (v-b)**(i-1)      
!          FF(i,1) =     1              / (v-b)**i      
!          FF(i,2) =  -  i              / (v-b)**(i+1)      
!          FF(i,3) =    (i*(i+1))       / (v-b)**(i+2)      
! 
!       ENDDO
!    
!    
!    END SUBROUTINE Q_F__T_v__MHG
!    !============================================================ 
! 
!   
!    
!    !============================================================ 
!    FUNCTION minimum_specific_volume__MHG() RESULT (v_min)
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
!    END FUNCTION minimum_specific_volume__MHG
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
!    FUNCTION P__T_v__MHG(T, v) RESULT (P)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
! 
!       P = (Rgas*T)/(v - b) + SUM( QQ(:,0)*FF(:,1) )
! 
! 
!    END FUNCTION P__T_v__MHG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION dP__dT_dv__MHG(T, v) RESULT(dP)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
! 
!       P_T  =      Rgas    / (v-b)     +  SUM( QQ(:,1)*FF(:,1) )
!       P_v  =   - (Rgas*T) / (v-b)**2  +  SUM( QQ(:,0)*FF(:,2) )
! 
!       dP(1) =  P_T;  dP(2) =  P_v
! 
! 
!    END FUNCTION dP__dT_dv__MHG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION d2P__dT2_dv2__MHG(T, v) RESULT(d2P)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
! 
!       P_TT =                             SUM( QQ(:,2)*FF(:,1) )
!       P_Tv =   -  Rgas    / (v-b)**2  +  SUM( QQ(:,1)*FF(:,2) )
!       P_vv =   (2*Rgas*T) / (v-b)**3  +  SUM( QQ(:,0)*FF(:,3) )
! 
!       d2P(1,1) = P_TT;      d2P(1,2) = P_Tv
!       d2P(2,1) = d2P(1,2);  d2P(2,2) = P_vv 
! 
! 
!    END FUNCTION d2P__dT2_dv2__MHG
!    !============================================================
!    
!    
!    
!    !============================================================ 
!    FUNCTION e__T_v__MHG(T, v) RESULT (e)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
!       
!       e = gas%e_0  +  phi__T(T)  +  SUM( ( T*QQ(:,1) - QQ(:,0))*FF(:,0) )
!       
! 
!    END FUNCTION e__T_v__MHG
!    !============================================================ 
!   
!   
! 
!    !============================================================ 
!    FUNCTION de__dT_dv__MHG(T, v) RESULT(de)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
!       
!       phip  = phip__T(T)
!            
!       e_T  =   phip  +  SUM(  T*QQ(:,2)          * FF(:,0) )
!       e_v  =            SUM( (T*QQ(:,1)-QQ(:,0)) * FF(:,1) )
!              
!       de(1) = e_T;  de(2) = e_v
! 
! 
!    END FUNCTION de__dT_dv__MHG
!    !============================================================
! 
! 
!    
!    !============================================================ 
!    FUNCTION d2e__dT2_dv2__MHG(T, v) RESULT(d2e)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
!       
!       phipp  = phipp__T(T)
!                       
!       e_TT =   phipp  +  SUM( (T*QQ(:,3)+QQ(:,2)) * FF(:,0) )
!       e_Tv =             SUM(  T*QQ(:,2)          * FF(:,1) )
!       e_vv =             SUM( (T*QQ(:,1)-QQ(:,0)) * FF(:,2) )
!              
!       d2e(1,1) = e_TT;  d2e(1,2) = e_Tv
!       d2e(2,1) = e_Tv;  d2e(2,2) = e_vv
! 
! 
!    END FUNCTION d2e__dT2_dv2__MHG
!    !============================================================
! 
! 
! 
!    !============================================================ 
!    FUNCTION s__T_v__MHG(T, v) RESULT (s)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF) 
!      
!       s = gas%s_0 + Rgas*LOG(v-b) + psi__T(T) &
!         + SUM( QQ(:,1)*FF(:,0) )
! 
! 
!    END FUNCTION s__T_v__MHG
!    !============================================================ 
! 
!  
!  
!    !============================================================ 
!    FUNCTION cv__T_v__MHG(T, v) RESULT (cv)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
! 
!       cv = phip__T(T)    + SUM( T*QQ(:,2)*FF(:,0) )
! 
! 
!    END FUNCTION cv__T_v__MHG
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION dcv__dT_dv__MHG(T, v) RESULT (dcv)
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
!       CALL Q_F__T_v__MHG(T, v,  QQ, FF)
! 
!       cv_T = phipp__T(T)  +  SUM( (T*QQ(:,3) + QQ(:,2))*FF(:,0) )
!       cv_v =                 SUM(  T*QQ(:,2)           *FF(:,1) )
! 
!       dcv(1) = cv_T;  dcv(2) = cv_v
!       
! 
!    END FUNCTION dcv__dT_dv__MHG
!    !============================================================ 
! 
! 
! 
 END MODULE martin_hou_gas
