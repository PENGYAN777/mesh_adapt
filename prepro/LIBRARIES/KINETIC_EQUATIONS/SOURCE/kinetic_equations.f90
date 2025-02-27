! ==============================================================================
!
!      Module: kinetic_equations
!
! Description: Computes moments of maxwellian distribution function.
!              Implementation of appendix B of Gas Kinetic Schemes for
!              unstedy compressible flow simulations. K.Xu (VKI 1998-03)
!
!              Moments of the maxwellian with respect single velocity and thermal 
!              degrees of freedom are computed. Moreover moments of the maxwellian 
!              multiplied by the vector PSI of collisional invariants and moments
!              for diadic product of PSI times itself (PSI x PSI) are computed.
!
!              For a single specie gas in 3D and in thermodynamic equilibrium, the
!              collisional invariants are:
!
!                PSI = [1, u, v, w, 0.5*(u^2 + v^2 + w^2 + xi^2)]
!
!
!      Author: Marco Fossati
!              Department of Aerospace Engineering
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
! ==============================================================================

   MODULE  kinetic_equations

   !----------------------------------------------------------------------------
   IMPLICIT NONE
   PRIVATE

   REAL(KIND=8), PARAMETER :: PI = 3.141592653589793
   
   INTEGER,      PARAMETER :: EQUILIBRIUM = 0,  &
   			      QUASI_EQUIL = 1,  &
   			      TRANS_REGIM = 2

   REAL(KIND=8), DIMENSION(0:10,3,3) :: IuT = 1.d0
   REAL(KIND=8), DIMENSION(0:10,3,3) :: IuP = 1.d0
   REAL(KIND=8), DIMENSION(0:10,3,3) :: IuN = 1.d0   
   REAL(KIND=8), DIMENSION(0:6,3)    :: Ixi = 1.d0

   REAL(KIND=8), DIMENSION(:,:,:), POINTER :: gradV
   REAL(KIND=8), DIMENSION(:,:),   POINTER :: gradD, gradT
   
   REAL(KIND=8) :: C1, C2, alpha, gamma, nsBreak
   REAL(KIND=8) :: Re_ref
   
   INTEGER :: sD, tau_fmla, phiFmla

   PUBLIC  init_kinetic_equations,          &
           mean_collision_time,             &
	   compute_scalar_moments,          &
	   P__ww, L__ww, K__ww,             &
	   nKnudsen__ww,                    &
	   iKnudsen__ww,                    &

	   IuP, IuN, &
	   I__c_PSI_g,     IP__c_PSI_g,     IN__c_PSI_g,      &
	   I__c_PSIxPSI_g, IP__c_PSIxPSI_g, IN__c_PSIxPSI_g,  &
	   
	   DPstressTensor, heatFlux, compute_DF_derivatives,  &

	   gamma, tau_fmla, gradD, gradV, gradT, alpha, phiFmla, &
	   nsBreak
   !----------------------------------------------------------------------------

   CONTAINS


   SUBROUTINE  init_kinetic_equations(ww, flow_rgm)
   !---------------------------------------------------------------------------
   USE structures,       ONLY: REF
   USE gas_properties,   ONLY: gas
   
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
   
   INTEGER, INTENT(IN), OPTIONAL :: flow_rgm
   !---------------------------------------------------------------------------

   sD = SIZE(ww,1) - 2

   gamma = ((gas % thd + 3.d0) + 2.d0)/ (gas % thd + 3.d0)

   Re_ref = REF % Re

   IF (PRESENT(flow_rgm)) THEN

     SELECT CASE (flow_rgm) 
 
       CASE (EQUILIBRIUM)   
       tau_fmla = 0	    
       C1 = 0.05	    
       C2 = 0.3 	    
 
       CASE (QUASI_EQUIL)   
       tau_fmla = 1	    

       CASE (TRANS_REGIM)   
       tau_fmla = 2

     END SELECT 	    

   ENDIF

   END SUBROUTINE  init_kinetic_equations





   FUNCTION  nKnudsen__ww(ww)   RESULT(Kn)
   !---------------------------------------------------------------------------
   USE metric_coefficients,    ONLY: cell, eta_fem => eta, &
   				     chi_b, xi_bp

   USE node_pair_structure,    ONLY: j_c_d_fem  => j_c_d, j_c_b
   USE nodes,		       ONLY: jd_jb
   USE np_quadrature   
   USE thermodynamics,         ONLY: T__P_v
   USE transport_properties,   ONLY: transport_coefficients
   USE structures,             ONLY: REF

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
   
   REAL(KIND=8), DIMENSION(SIZE(ww,2)) :: Kn
   
   REAL(KIND=8), DIMENSION(sD,SIZE(ww,2)) :: gradD, gradV, gradT
   
   REAL(KIND=8), DIMENSION(SIZE(ww,2)) :: r, P, V, T
   
   REAL(KIND=8) :: mu, kk, mFpath, KnD, KnV, KnT
   
   INTEGER :: j
   !---------------------------------------------------------------------------

   DO j = 1, SIZE(ww,2)

     r(j) = ww(1,j)
     P(j) = P__ww(ww(:,j))
     V(j) = SQRT(SUM((ww(2:sD+1,j) / r(j))**2))
     T(j) = T__P_v(P(j), 1.d0/r(j))
!     CALL transport_coefficients(T(j), 1.d0/r(j), mu, kk)

!     mFpath = 16.d0/(5*SQRT(PI)) * mu / (r(j)*SQRT(2*P(j)/r(j)))
!     Kn(j) = mFpath

   ENDDO
   
   
!    ! Computation of the gradient of the density velocity 
!    ! modulus and temperature   
   gradD = np_quadr_w_Gu(j_c_d_fem, j_c_b, jd_jb,  &
                       eta_fem, chi_b, xi_bp, r)

   gradV = np_quadr_w_Gu(j_c_d_fem, j_c_b, jd_jb,  &
                       eta_fem, chi_b, xi_bp, V)
                              
   gradT = np_quadr_w_Gu(j_c_d_fem, j_c_b, jd_jb,  &
                       eta_fem, chi_b, xi_bp, T)


   DO j = 1, SIZE(ww,2)
   
     KnD = 0.d0
     KnV = 0.d0
     KnT = 0.d0
      
     gradD(:,j) = gradD(:,j) / cell(j)
     gradV(:,j) = gradV(:,j) / cell(j)
     gradT(:,j) = gradT(:,j) / cell(j)
 
     CALL transport_coefficients(T(j), 1.d0/r(j), mu, kk)

     mFpath = 16.d0/(5*SQRT(PI)) * mu / (r(j)*SQRT(2*P(j)/r(j)))
   
     KnD = mFpath * SQRT(SUM(gradD(:,j)**2)) / r(j)
     KnT = mFpath * SQRT(SUM(gradT(:,j)**2)) / T(j)
     
     IF (V(j) /= 0.d0) THEN
       KnV = mFpath * SQRT(SUM(gradV(:,j)**2)) / V(j)
     ELSE
       KnV = 0.d0
     ENDIF  
     
     Kn(j) = MAX(KnD, KnT, KnV)

   ENDDO
   
   END FUNCTION  nKnudsen__ww





   FUNCTION  iKnudsen__ww(ww_0, ww_l, ww_r, Drr_ij)   RESULT(Kn)
   !---------------------------------------------------------------------------
   USE  thermodynamics,         ONLY: T__P_v
   USE  transport_properties,   ONLY: transport_coefficients
   
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww_0, ww_l, ww_r
   
   REAL(KIND=8) :: Drr_ij
   
   REAL(KIND=8) :: Kn
      
   REAL(KIND=8) :: r0, P0, V0, T0,  &
                   rl, Pl, Vl, Tl,  &   
		   rr, Pr, Vr, Tr,  &
		   gD, gV, gT
		   
   REAL(KIND=8) :: mu, kk, mFpath, KnD, KnV, KnT
   !---------------------------------------------------------------------------
     
   KnD = 0.d0
   KnV = 0.d0
   KnT = 0.d0
       
   r0 = ww_0(1)
   V0 = SQRT(SUM((ww_0(2:sD+1)/r0)**2))
   P0 = P__ww(ww_0)
   T0 = T__P_v(P0, 1/r0)

   CALL transport_coefficients(T0, 1.d0/r0, mu, kk)
         
   mFpath = 16.d0/(5*SQRT(PI)) * mu / (r0*SQRT(2*P0/r0))

   rl = ww_l(1)
   Vl = SQRT(SUM((ww_l(2:sD+1)/rl)**2))
   Pl = P__ww(ww_l)
   Tl = T__P_v(Pl, 1/rl)

   rr = ww_r(1)
   Vr = SQRT(SUM((ww_r(2:sD+1)/rr)**2))
   Pr = P__ww(ww_r)
   Tr = T__P_v(Pr, 1/rr)
  
   gD = ABS(rr - rl) / Drr_ij
   gV = ABS(Vr - Vl) / Drr_ij
   gT = ABS(Tr - Tl) / Drr_ij
  
   KnD = mFpath * gD / r0
   KnT = mFpath * gT / T0
   
   IF (V0 /= 0.d0) THEN
     KnV = mFpath * gV / V0
   ELSE
     KnV = 0.d0
   ENDIF  
   
   Kn = MAX(KnD, KnT)
   
   END FUNCTION  iKnudsen__ww





   FUNCTION  K__ww(ww)   RESULT(k)
   !---------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww
   REAL(KIND=8) :: k
   !---------------------------------------------------------------------------

   k = (5.0 - 3.0*gamma) / (gamma - 1.0) + (3.d0 - (SIZE(ww,1) - 2.0))

   END FUNCTION  K__ww





   FUNCTION  L__ww(ww)   RESULT(lambda)
   !---------------------------------------------------------------------------   
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww
   REAL(KIND=8) :: lambda

   REAL(KIND=8), DIMENSION(SIZE(ww,1)-2) :: u
   REAL(KIND=8) :: k
   !---------------------------------------------------------------------------

   u = ww(2:sD+1) / ww(1)
   k = K__ww(ww)
   lambda = (k+sD) * ww(1) / (4.0*(ww(sD+2) - 0.5*ww(1)*SUM(u**2)))
   
   END FUNCTION  L__ww





   FUNCTION  P__ww(ww)   RESULT(P)
   !---------------------------------------------------------------------------
   USE thermodynamics,   ONLY: P__e_r

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww
   REAL(KIND=8) :: P

   REAL(KIND=8) :: mod_u_2, rho, e

   INTEGER :: k
   !---------------------------------------------------------------------------

   k = SIZE(ww)-2

   mod_u_2 = SUM((ww(2:k+1)/ww(1))**2)

   rho = ww(1)
   e = ww(k+2) / rho - 0.5 * mod_u_2 

   P = P__e_r(e, rho)

   END FUNCTION  P__ww





   FUNCTION  mean_collision_time(dt, ww_l, ww_r, ww_0, iKN, pMI, Gx, Gt, ww_0_xx)   RESULT(tau)
   !-----------------------------------------------------------------------------------------------------
   USE  lin_algebra,            ONLY: diad_product
   USE  thermodynamics,         ONLY: T__P_v
   USE  transport_properties,   ONLY: transport_coefficients

   IMPLICIT NONE

   REAL(KIND=8),	       INTENT(IN) :: dt
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww_l, ww_r, ww_0
   
   REAL(KIND=8),                           OPTIONAL :: iKN
   REAL(KIND=8), DIMENSION(:), INTENT(IN), OPTIONAL :: pMI, Gx, Gt, ww_0_xx
     
   REAL(KIND=8) :: tau

   REAL(KIND=8), DIMENSION(sD+2,sD+2,sD+2) :: I_u0_PSI3, I_u1_PSI3, I_u2_PSI3

   REAL(KIND=8), DIMENSION(sD+2,sD+2) :: I_u0_PSI2, I_u1_PSI2, I_u2_PSI2, I_u3_PSI2, I_u4_PSI2
   REAL(KIND=8), DIMENSION(sD+2,sD+2) :: Gx_x_Gx, Gt_x_Gx, Gt_x_Gt
   
   REAL(KIND=8), DIMENSION((sD+2)*(sD+3)/2) :: sysMat
   REAL(KIND=8), DIMENSION(sD) :: u_0

   REAL(KIND=8), DIMENSION(sD+2) :: Gxx, Gxt, Gtt

   REAL(KIND=8), DIMENSION(sD+2) :: I_u0_PSI3_GxGx, I_u0_PSI3_GxGt, I_u0_PSI3_GtGt, &
                   		    I_u1_PSI3_GxGx, I_u1_PSI3_GxGt, &
		   		    I_u2_PSI3_GxGx

   REAL(KIND=8) :: I_u0_PSI2_GtGt, &
                   I_u1_PSI2_GtGt, I_u1_PSI2_GtGx, &
		   I_u2_PSI2_GtGt, I_u2_PSI2_GtGx, I_u2_PSI2_GxGx, &
		   I_u3_PSI2_GtGx, I_u3_PSI2_GxGx, &
		   I_u4_PSI2_GxGx
   
   REAL(KIND=8) :: I_Dg, I_DDg, cTerm

   REAL(KIND=8) :: lambda_l, k_l, r_l, lambda_r, k_r, r_r,  &
   		   r_0, P_0, T_0, mu_0, kk_0, C, U0

   INTEGER, DIMENSION(SIZE(ww_0)) :: IPIV
   INTEGER :: info, p, q

   ! REAL(KIND=8) :: P_l, sqL_R_l, P_r, sqL_R_r, lambda_0, k_0
   !-----------------------------------------------------------------------------------------------------

   SELECT CASE (tau_fmla)

     CASE (EQUILIBRIUM)

       r_l = ww_l(1);	k_l = K__ww(ww_l);   lambda_l = L__ww(ww_l)
       r_r = ww_r(1);	k_r = K__ww(ww_r);   lambda_r = L__ww(ww_r)

       C = C2 * ABS(r_l/lambda_l - r_r/lambda_r)  &
              / ABS(r_l/lambda_l + r_r/lambda_r)

       tau = C1*dt + dt*C


!      CASE (EQUILIBRIUM)
!      ! Alternative not used anymore ------------------------------
! 
!        r_l = ww_l(1);	k_l = K__ww(ww_l);   lambda_l = L__ww(ww_l)
!        r_r = ww_r(1);	k_r = K__ww(ww_r);   lambda_r = L__ww(ww_r)
!        r_0 = ww_0(1);	k_0 = K__ww(ww_0);   lambda_0 = L__ww(ww_0)
! 
!        ! tau = C1*SQRT(lambda_0)/r_0 + C2*SQRT(lambda_0)/r_0 * &
!        !       ABS(r_l/(2*lambda_l) - r_r/(2*lambda_r))
! 
!        P_l = P__ww(ww_l)
!        P_r = P__ww(ww_r)
!        
!        sqL_R_l = SQRT(lambda_l) / r_l
!        sqL_R_r = SQRT(lambda_r) / r_r
! 	      
! 
!        tau = C1*SQRT(lambda_0)/r_0 + C2*(ABS(sqL_R_l - sqL_R_r)/(sqL_R_l + sqL_R_r)) &
!            *(ABS(P_l - P_r)/(P_l + P_r))



     CASE (QUASI_EQUIL)

       r_l = ww_l(1);	k_l = K__ww(ww_l);   lambda_l = L__ww(ww_l)
       r_r = ww_r(1);	k_r = K__ww(ww_r);   lambda_r = L__ww(ww_r)

       r_0 = ww_0(1)
       P_0 = P__ww(ww_0)
       T_0 = T__P_v(P_0, 1/r_0)

       CALL transport_coefficients(T_0, 1.d0/r_0, mu_0, kk_0)

       tau = 1/Re_ref  *  mu_0 / P_0 + dt*(  ABS(r_l/lambda_l - r_r/lambda_r)  &
   			                   / ABS(r_l/lambda_l + r_r/lambda_r))
					   


     CASE (TRANS_REGIM)

       ! Works only for resolved regions, i.e. profiles 
       ! are assumed to be smooth
       ! Classical tau from Navier-Stokes assumption       
       r_0 = ww_0(1)
       u_0 = ww_0(2:sD+1)/r_0
       P_0 = P__ww(ww_0)
       T_0 = T__P_v(P_0, 1/r_0)

       CALL transport_coefficients(T_0, 1.d0/r_0, mu_0, kk_0)

       tau = 1/Re_ref  *  mu_0 / P_0
  
       IF (iKN < nsBreak) THEN
       
         cTerm = 0.d0
       
       ELSE
	
         ! Coefficients for Mawellian second space derivative, Gxx ---------------------------
         Gx_x_Gx = diad_product(Gx,Gx)
 
         I_u0_PSI3 = r_0 * I__c_PSI_PSIxPSI_g((/0,0,0/),'I')

         I_u0_PSI3_GxGx = 0.d0

         ! Double scalar product
         DO p = 1, SIZE(ww_0)
           DO q = 1, SIZE(ww_0)
             I_u0_PSI3_GxGx = I_u0_PSI3_GxGx + Gx_x_Gx(p,q)*I_u0_PSI3(p,q,:)
           ENDDO
         ENDDO
 
         Gxx = ww_0_xx - I_u0_PSI3_GxGx
 
         sysMat = pMI
         CALL DSPSV('U', SIZE(ww_0), 1, sysMat, IPIV, Gxx, SIZE(ww_0), info)

 
         ! Coefficients for Maxwellian mixed derivative, Gxt ---------------------------------
         Gt_x_Gx = diad_product(Gt,Gx)
 
         I_u1_PSI3 = r_0 * I__c_PSI_PSIxPSI_g((/1,0,0/),'I')
 
         I_u0_PSI3_GxGt = 0.d0
         I_u1_PSI3_GxGx = 0.d0
       
         DO p = 1, SIZE(ww_0)
           DO q = 1, SIZE(ww_0)
           
             I_u0_PSI3_GxGt = I_u0_PSI3_GxGt + Gt_x_Gx(p,q)*I_u0_PSI3(p,q,:)
             I_u1_PSI3_GxGx = I_u1_PSI3_GxGx + Gx_x_Gx(p,q)*I_u1_PSI3(p,q,:)

           ENDDO
         ENDDO      
 
         Gxt = - I_u0_PSI3_GxGt - I_u1_PSI3_GxGx - MATMUL(r_0*I__c_PSIxPSI_g((/1,0,0/),'I'),Gxx)
 
         sysMat = pMI
         CALL DSPSV('U', SIZE(ww_0), 1, sysMat, IPIV, Gxt, SIZE(ww_0), info)
  
  
	  ! Coefficients for Maxwellian second time derivative, Gtt ----------------------------
	  Gt_x_Gt = diad_product(Gt,Gt)
 
	  I_u2_PSI3 = r_0 * I__c_PSI_PSIxPSI_g((/2,0,0/),'I')
 
	  I_u0_PSI3_GtGt = 0.d0
	  I_u1_PSI3_GxGt = 0.d0
	  I_u2_PSI3_GxGx = 0.d0
	  
	  DO p = 1, SIZE(ww_0)
	    DO q = 1, SIZE(ww_0)
	    
	      I_u0_PSI3_GtGt = I_u0_PSI3_GtGt + Gt_x_Gt(p,q)*I_u0_PSI3(p,q,:)
	      I_u1_PSI3_GxGt = I_u1_PSI3_GxGt + Gt_x_Gx(p,q)*I_u1_PSI3(p,q,:)
	      I_u2_PSI3_GxGx = I_u2_PSI3_GxGx + Gx_x_Gx(p,q)*I_u2_PSI3(p,q,:)
	  
	    ENDDO
	  ENDDO
 
	  Gtt = - I_u0_PSI3_GtGt - 2*I_u1_PSI3_GxGt - I_u2_PSI3_GxGx   &
	 	- 2*MATMUL(r_0* I__c_PSIxPSI_g((/1,0,0/),'I'),Gxt)     &
	 	-   MATMUL(r_0* I__c_PSIxPSI_g((/2,0,0/),'I'),Gxx)
	  
	  sysMat = pMI
	  CALL DSPSV('U', SIZE(ww_0), 1, sysMat, IPIV, Gtt, SIZE(ww_0), info)


	  !    _		_
	  !  _/ (u-U)**2 Dg = _/ (u-U)**2 (A + au) g --------------------------------------------
	  U0 = ww_0(2)/ww_0(1)

	  !    _	    _
	  !  _/ u**2 Dg = _/ u**2 (A + au) g	   
	  I_Dg = DOT_PRODUCT(r_0 * I__c_PSI_g((/2,0,0/),'I'),Gt)  &
	       + DOT_PRODUCT(r_0 * I__c_PSI_g((/3,0,0/),'I'),Gx)

          IF (phiFmla == 1) THEN
	    !	 _	      _
	    !  _/ -2uU Dg = _/ -2uU (A + au) g
	    I_Dg = I_Dg - 2*U0 * ( DOT_PRODUCT(r_0 * I__c_PSI_g((/1,0,0/),'I'),Gt)  &
	 			 + DOT_PRODUCT(r_0 * I__c_PSI_g((/2,0,0/),'I'),Gx))

	    !	 _	     _
	    !  _/ U^2 Dg = _/ U^2 (A + au) g
	    I_Dg = I_Dg + U0**2 * ( DOT_PRODUCT(r_0 * I__c_PSI_g((/0,0,0/),'I'),Gt)  &
	 			  + DOT_PRODUCT(r_0 * I__c_PSI_g((/1,0,0/),'I'),Gx))
          ENDIF


	  !    _		   _
	  !  _/ (u-U)**2 D**2g = _/ (u-U)**2 [(A^2 + B) + 2u(C + Aa) + u^2(a^2 + b)] g ----------	  
	  !    _
	  !  _/ u^n PSI x PSI g   with n = (0,...,4)
	  I_u0_PSI2 = r_0*I__c_PSIxPSI_g((/0,0,0/),'I')
	  I_u1_PSI2 = r_0*I__c_PSIxPSI_g((/1,0,0/),'I')
	  I_u2_PSI2 = r_0*I__c_PSIxPSI_g((/2,0,0/),'I')
	  I_u3_PSI2 = r_0*I__c_PSIxPSI_g((/3,0,0/),'I')      
	  I_u4_PSI2 = r_0*I__c_PSIxPSI_g((/4,0,0/),'I')
	  
	  !		    _
	  !  (Gi x Gj) *  _/ u^n PSI x PSI g  with i,j = x,t and n = (0,...,4) 
	  I_u0_PSI2_GtGt = 0.d0
	  
	  I_u1_PSI2_GtGt = 0.d0
	  I_u1_PSI2_GtGx = 0.d0
	  
	  I_u2_PSI2_GtGt = 0.d0
	  I_u2_PSI2_GtGx = 0.d0
	  I_u2_PSI2_GxGx = 0.d0
	  
	  I_u3_PSI2_GtGx = 0.d0
	  I_u3_PSI2_GxGx = 0.d0
	  
	  I_u4_PSI2_GxGx = 0.d0

	  DO p = 1, SIZE(ww_0)
	    DO q = 1, SIZE(ww_0)

	      I_u0_PSI2_GtGt = I_u0_PSI2_GtGt + Gt_x_Gt(p,q) * I_u0_PSI2(p,q)
	      
	      I_u1_PSI2_GtGt = I_u1_PSI2_GtGt + Gt_x_Gt(p,q) * I_u1_PSI2(p,q)	      
	      I_u1_PSI2_GtGx = I_u1_PSI2_GtGx + Gt_x_Gx(p,q) * I_u1_PSI2(p,q)	      
	      
	      I_u2_PSI2_GtGt = I_u2_PSI2_GtGt + Gt_x_Gt(p,q) * I_u2_PSI2(p,q)
	      I_u2_PSI2_GtGx = I_u2_PSI2_GtGx + Gt_x_Gx(p,q) * I_u2_PSI2(p,q)
	      I_u2_PSI2_GxGx = I_u2_PSI2_GxGx + Gx_x_Gx(p,q) * I_u2_PSI2(p,q)
	      
	      I_u3_PSI2_GtGx = I_u3_PSI2_GtGx + Gt_x_Gx(p,q) * I_u3_PSI2(p,q)
	      I_u3_PSI2_GxGx = I_u3_PSI2_GxGx + Gx_x_Gx(p,q) * I_u3_PSI2(p,q)
	      
	      I_u4_PSI2_GxGx = I_u4_PSI2_GxGx + Gx_x_Gx(p,q) * I_u4_PSI2(p,q)
	    
	    ENDDO
	  ENDDO


	  !    _	       _
	  !  _/ u**2 D**2g = _/ u**2 [(A^2 + B) + 2u(C + Aa) + u^2(a^2 + b)] g
	  I_DDg = I_u2_PSI2_GtGt + 2*I_u3_PSI2_GtGx + I_u4_PSI2_GxGx	 &
	 	  +   r_0 * DOT_PRODUCT(I__c_PSI_g((/2,0,0/),'I'), Gtt)  &
	 	  + 2*r_0 * DOT_PRODUCT(I__c_PSI_g((/3,0,0/),'I'), Gxt)  &  
	 	  +   r_0 * DOT_PRODUCT(I__c_PSI_g((/4,0,0/),'I'), Gxx)

          IF (phiFmla == 1) THEN
	    !	 _		 _
	    !  _/ -2uU D**2g = _/ -2uU [(A^2 + B) + 2u(C + Aa) + u^2(a^2 + b)] g
	    I_DDg = I_DDg - 2*U0 * (I_u1_PSI2_GtGt + 2*I_u2_PSI2_GtGx + I_u3_PSI2_GxGx      &
	 			     +   r_0 * DOT_PRODUCT(I__c_PSI_g((/1,0,0/),'I'), Gtt)  &
	 			     + 2*r_0 * DOT_PRODUCT(I__c_PSI_g((/2,0,0/),'I'), Gxt)  &
	 			     +   r_0 * DOT_PRODUCT(I__c_PSI_g((/3,0,0/),'I'), Gxx))


	    !	 _		 _
	    !  _/ U**2 D**2g = _/ U**2 [(A^2 + B) + 2u(C + Aa) + u^2(a^2 + b)] g
	    I_DDg = I_DDg + U0**2 * (I_u0_PSI2_GtGt + 2*I_u1_PSI2_GtGx + I_u2_PSI2_GxGx      &
	 			      +   r_0 * DOT_PRODUCT(I__c_PSI_g((/0,0,0/),'I'), Gtt)  &
	 			      + 2*r_0 * DOT_PRODUCT(I__c_PSI_g((/1,0,0/),'I'), Gxt)  &
	 			      +   r_0 * DOT_PRODUCT(I__c_PSI_g((/2,0,0/),'I'), Gxx))
          ENDIF

         cTerm = MAX(-alpha, MIN(tau*I_DDg/(I_Dg + 1.d-15),0.d0))

       ENDIF
       !---------------------------------------------------------------------

       ! Regularized tau formula
       tau = tau / (1 + cTerm)


    END SELECT

   END FUNCTION  mean_collision_time










   !--------------------------------------------------------------------------------------------
   ! MOMENTS OF THE MAXWELLIAN: SCALAR INTEGRALS
   !--------------------------------------------------------------------------------------------

   SUBROUTINE compute_scalar_moments(ww, state)

   ! VELOCITY MOMENTS:
   !					      _+oo
   ! if U == MACROscopic velocity in X  =>  _/  u^n g dudvdwdxi  =  <u^n><v^0><w^0>
   !					     -oo
   !					      _+oo
   ! if U == MACROscopic velocity in Y  =>  _/  v^n g dudvdwdxi  =  <u^0><v^n><w^0>
   !					     -oo
   !					      _+oo
   ! if U == MACROscopic velocity in Z  =>  _/  w^n g dudvdwdxi  =  <u^0><v^0><w^n>
   !					     -oo
   !
   !					      _+oo  _+oo
   ! if U == MACROscopic velocity in X  =>  _/    _/  u^n g  du dvdwdxi  =  <u^n><v^0><w^0>
   !					     -oo   0
   !					      _+oo  _+oo
   ! if U == MACROscopic velocity in Y  =>  _/    _/  v^n g  dv dudwdxi  =  <u^0><v^n><w^0>
   !					     -oo   0
   !					      _+oo  _+oo
   ! if U == MACROscopic velocity in Z  =>  _/    _/  w^n g  dw dudvdxi  =  <u^0><v^0><w^n>
   !					     -oo   0
   !
   !					      _+oo  _0
   ! if U == MACROscopic velocity in X  =>  _/    _/  u^n g  du dvdwdxi  =  <u^n><v^0><w^0>
   !					     -oo   -oo
   !					      _+oo  _0
   ! if U == MACROscopic velocity in Y  =>  _/    _/  v^n g  dv dudwdxi  =  <u^0><v^n><w^0>
   !					     -oo   -oo
   !					      _+oo  _0
   ! if U == MACROscopic velocity in Z  =>  _/    _/  w^n g  dw dudvdxi  =  <u^0><v^0><w^n>
   !					     -oo   -oo   
   !
   ! ENERGY MOMENTS:
   ! 	_+oo
   !  _/  xi^n g dudvdwdxi  =  <u^0><v^0><w^0><xi^n>
   !   -oo
   !---------------------------------------------------------------------------------------------
   IMPLICIT NONE
   
   REAL(KIND=8),     DIMENSION(:), INTENT(IN) :: ww
   CHARACTER(LEN=1),               INTENT(IN) :: state
   
   REAL(KIND=8) :: U, K, l
   REAL(KIND=8) :: DERFC
   
   INTEGER :: i, c, n
   !---------------------------------------------------------------------------------------------

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3

   K = K__ww(ww)   
   l = L__ww(ww)

   ! Velocity Moments
   DO i = 1, sD
   
     U = ww(i+1) / ww(1)

     IuT(0,i,c) = 1.d0
     IuP(0,i,c) = 0.5*DERFC(-SQRT(l)*U)
     IuN(0,i,c) = 0.5*DERFC( SQRT(l)*U)
       
     IuT(1,i,c) = U
     IuP(1,i,c) = U*IuP(0,i,c) + 0.5*(EXP(-l*U*U))/(SQRT(PI*l))
     IuN(1,i,c) = U*IuN(0,i,c) - 0.5*(EXP(-l*U*U))/(SQRT(PI*l))

     DO n = 0,8
       IuT(n+2,i,c) = U*IuT(n+1,i,c) + (n+1)/(2*l)*IuT(n,i,c)
       IuP(n+2,i,c) = U*IuP(n+1,i,c) + (n+1)/(2*l)*IuP(n,i,c)
       IuN(n+2,i,c) = U*IuN(n+1,i,c) + (n+1)/(2*l)*IuN(n,i,c)
     ENDDO    

   ENDDO
   
   ! Energy Moments
   Ixi(2,c) = 0.5 * K / l
   Ixi(4,c) = 0.5 * Ixi(2,c) * (K+2.d0)/ l   !3*K/(4*l*l) + K*(K - 1)/(4*l*l)
   Ixi(6,c) = 0.5 * Ixi(4,c) * (K+4.d0)/ l
   
   END SUBROUTINE compute_scalar_moments







   !--------------------------------------------------------------------------------------------
   ! MOMENTS OF THE MAXWELLIAN: VECTOR INTEGERALS
   !--------------------------------------------------------------------------------------------

   FUNCTION  I__c_PSI_g(n, state)   RESULT(Int)

   !	_+oo
   !  _/  u^n PSI g dudvdwdxi => 
   !   -oo
   !
   !  component-wise  =>
   !
   !	_+oo
   !  _/  u^n g dudvdwdxi =  <u^n><v^0><w^0>
   !   -oo
   !
   !	_+oo
   !  _/  u^(n+1) g dudvdwdxi =  <u^(n+1)><v^0><w^0>
   !   -oo
   !
   !	_+oo
   !  _/  u^n v g dudvdwdxi =  <u^n><v^1><w^0>
   !   -oo
   !
   !	_+oo
   !  _/  u^n w g dudvdwdxi =  <u^n><v^0><w^1>
   !   -oo
   !
   !	_+oo
   !  _/  0.5*(u^2 + v^2 + w^2 + xi^2) u^n g dudvdwdxi = 0.5*(<u^(n+2)> <v^0> <w^0> +
   !   -oo						      <u^n>	<v^2> <w^0> +
   !							      <u^n>	<v^0> <w^2> +
   !							      <u^n> <v^0> <w^0> <xi^2>)
   !--------------------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),      INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2) :: Int

   INTEGER :: i, c, nu, nv, nw
   !--------------------------------------------------------------------------------------------
   
   nu = n(1)
   nv = n(2)
   nw = n(3)
   
   Int = 0.d0

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3

   ! Mass component
   Int(1) = IuT(nu,1,c) * IuT(nv,2,c) * IuT(nw,3,c)

   ! Velocity components
   DO i = 1, sD

     IF (i == 1)  Int(i+1) = IuT(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
     IF (i == 2)  Int(i+1) = IuT(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
     IF (i == 3)  Int(i+1) = IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)

     ! Accumulating energy components
     IF (i == 1)  Int(sD+2) =             IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
     IF (i == 2)  Int(sD+2) = Int(sD+2) + IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c)
     IF (i == 3)  Int(sD+2) = Int(sD+2) + IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c)

   ENDDO

   ! Energy component
   Int(sD+2) = 0.5*(Int(sD+2) + IuT(nu,1,c) * IuT(nv,2,c) * IuT(nw,3,c) * Ixi(2,c))

   END FUNCTION  I__c_PSI_g





   FUNCTION  IP__c_PSI_g(n, state)   RESULT(IntP)

   !	_+oo _+oo
   !  _/   _/  u^n PSI g du dvdwdxi => 
   !   -oo  0
   !
   !  component-wise  =>
   !
   !	_+oo _+oo
   !  _/   _/  u^n g du dvdwdxi =  <u^n><v^0><w^0>
   !   -oo  0
   !
   !	_+oo _+oo
   !  _/   _/  u^(n+1) g du dvdwdxi =  <u^(n+1)><v^0><w^0>
   !   -oo  0
   !
   !	_+oo _ +oo
   !  _/   _/  u^n v g du dvdwdxi =  <u^n><v^1><w^0>
   !   -oo  0
   !
   !	_+oo _+oo
   !  _/   _/  u^n w g du dvdwdxi =  <u^n><v^0><w^1>
   !   -oo  0
   !
   !	_+oo _+oo
   !  _/   _/  0.5*(u^2 + v^2 + w^2 + xi^2) u^n g du dvdwdxi =  0.5*(<u^(n+2)> <v^0> <w^0> +
   !   -oo  0							     <u^n>     <v^2> <w^0> +
   !								     <u^n>     <v^0> <w^2> +
   !								     <u^n>     <v^0> <w^0> <xi^2>)
   !--------------------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),      INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2) :: IntP

   INTEGER :: i, c, nu, nv, nw
   !--------------------------------------------------------------------------------------------
   
   nu = n(1)
   nv = n(2)
   nw = n(3)
   
   IntP = 0.d0

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3

   ! Mass component
   IntP(1) = IuP(nu,1,c) * IuT(nv,2,c) * IuT(nw,3,c)

   ! Velocity components
   DO i = 1, sD

     IF (i == 1)  IntP(i+1) = IuP(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
     IF (i == 2)  IntP(i+1) = IuP(nu  ,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
     IF (i == 3)  IntP(i+1) = IuP(nu  ,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)

     ! Accumulating energy components
     IF (i == 1)  IntP(sD+2) =              IuP(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
     IF (i == 2)  IntP(sD+2) = IntP(sD+2) + IuP(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c)
     IF (i == 3)  IntP(sD+2) = IntP(sD+2) + IuP(nu  ,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c)

   ENDDO

   ! Energy component
   IntP(sD+2) = 0.5*(IntP(sD+2) + IuP(nu,1,c) * IuT(nv,2,c) * IuT(nw,3,c) * Ixi(2,c))

   END FUNCTION  IP__c_PSI_g	





   FUNCTION  IN__c_PSI_g(n, state)   RESULT(IntN)

   !	_+oo _0
   !  _/   _/  u^n PSI g du dvdwdxi => 
   !   -oo  -oo
   !
   !  component-wise  =>
   !
   !	_+oo _0
   !  _/   _/  u^n g du dvdwdxi =  <u^n><v^0><w^0>
   !   -oo  -oo
   !
   !	_+oo _0
   !  _/   _/  u^(n+1) g du dvdwdxi =  <u^(n+1)><v^0><w^0>
   !   -oo  -oo
   !
   !	_+oo _0
   !  _/   _/  u^n v g du dvdwdxi =  <u^n><v^1><w^0>
   !   -oo  -oo
   !
   !	_+oo _0
   !  _/   _/  u^n w g du dvdwdxi =  <u^n><v^0><w^1>
   !   -oo  -oo
   !
   !	_+oo _0
   !  _/   _/  0.5*(u^2 + v^2 + w^2 + xi^2) u^n g du dvdwdxi =  0.5*(<u^(n+2)> <v^0> <w^0> +
   !   -oo  -oo 						    <u^n>     <v^2> <w^0> +
   !								    <u^n>     <v^0> <w^2> +
   !								    <u^n> <v^0> <w^0> <xi^2>)
   !--------------------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER, DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),      INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2) :: IntN

   INTEGER :: i, c, nu, nv, nw
   !--------------------------------------------------------------------------------------------

   nu = n(1)
   nv = n(2)
   nw = n(3)
   
   IntN = 0.d0

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3

   ! Mass component
   IntN(1) = IuN(nu,1,c) * IuT(nv,2,c) * IuT(nw,3,c)


   ! Velocity components
   DO i = 1, sD

     IF (i == 1)  IntN(i+1) = IuN(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
     IF (i == 2)  IntN(i+1) = IuN(nu  ,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
     IF (i == 3)  IntN(i+1) = IuN(nu  ,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)

     ! Accumulating energy components
     IF (i == 1)  IntN(sD+2) =              IuN(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
     IF (i == 2)  IntN(sD+2) = IntN(sD+2) + IuN(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c)
     IF (i == 3)  IntN(sD+2) = IntN(sD+2) + IuN(nu  ,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c)

   ENDDO

   ! Energy component
   IntN(sD+2) = 0.5*(IntN(sD+2) + IuN(nu,1,c) * IuT(nv,2,c) * IuT(nw,3,c) * Ixi(2,c))

   END FUNCTION  IN__c_PSI_g	





   FUNCTION  I__c_xi_PSI_g(n, state)   RESULT(Int)
   !--------------------------------------------------------------------------------------------
   INTEGER, DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),      INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2) :: Int

   INTEGER :: i, c, nu, nv, nw
   !--------------------------------------------------------------------------------------------
   
   nu = n(1)
   nv = n(2)
   nw = n(3)
   
   Int = 0.d0

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3
   
   DO i = 1, sD
   
     ! Mass component -------------------------------------------------------------------------------
     IF (i == 1)  Int(1) = IuT(nu+2,1,c) * IuT(nv,2,c) * IuT(nw,3,c)  &
                         + IuT(nu,  1,c) * IuT(nv,2,c) * IuT(nw,3,c) * Ixi(2,c)

     IF (i == 2)  Int(1) = Int(1) + IuT(nu,1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c)
     IF (i == 3)  Int(1) = Int(1) + IuT(nu,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) 
                  
     
     ! Velocity components --------------------------------------------------------------------------
     IF (i == 1)  Int(2) = IuT(nu+3,1,c) * IuT(nv,2,c) * IuT(nw,3,c)  &
                         + IuT(nu+1,1,c) * IuT(nv,2,c) * IuT(nw,3,c) * Ixi(2,c)
     
     IF (i == 2) THEN
     
       Int(2) = Int(2) + IuT(nu+1,1,c) * IuT(nv+2,2,c) * IuT(nw,3,c)
       
       Int(3) = IuT(nu+2,1,c) * IuT(nv+1,2,c) * IuT(nw,3,c)  &
              + IuT(nu,  1,c) * IuT(nv+3,2,c) * IuT(nw,3,c)  &
              + IuT(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,3,c) * Ixi(2,c)
     ENDIF
     
     IF (i == 3) THEN
     
       Int(2) = Int(2) + IuT(nu+1,1,c) * IuT(nv,2,c) * IuT(nw+2,3,c)

       Int(3) = Int(3) + IuT(nu,  1,c) * IuT(nv,2,c) * IuT(nw+2,3,c)

       Int(4) = IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)  &
              + IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+1,3,c)  &
	      + IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+3,3,c)  &
              + IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c) * Ixi(2,c)
     ENDIF
     
     
     ! Energy components ---------------------------------------------------------------------------
     IF (i == 1)  Int(sD+2) =   IuT(nu+4,1,c) * IuT(nv,2,c) * IuT(nw,3,c)  &
                            + 2*IuT(nu+2,1,c) * IuT(nv,2,c) * IuT(nw,3,c) * Ixi(2,c)  &
			    +   IuT(nu,  1,c) * IuT(nv,2,c) * IuT(nw,3,c) * Ixi(4,c)

     IF (i == 2)  Int(sD+2) = Int(sD+2) +   IuT(nu,  1,c) * IuT(nv+4,2,c) * IuT(nw,3,c)  &
                                        + 2*IuT(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw,3,c)  &	  
                                        + 2*IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,3,c) * Ixi(2,c)  	  

     IF (i == 3)  Int(sD+2) = Int(sD+2) +   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+4,3,c)  &
                                        + 2*IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c)  &
                                        + 2*IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c)  &
                                        + 2*IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) * Ixi(2,c)
   ENDDO

   Int = 0.5*Int
   Int(sD+2) = 0.5*Int(sD+2)

   END FUNCTION  I__c_xi_PSI_g
   
   
   
   
   
   
   
   FUNCTION  I__c_xi2_PSI_g(n, state)   RESULT(Int)
   !--------------------------------------------------------------------------------------------
   INTEGER, DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),      INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2) :: Int

   INTEGER :: i, c, nu, nv, nw
   !--------------------------------------------------------------------------------------------
   
   nu = n(1)
   nv = n(2)
   nw = n(3)
   
   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3
      
   DO i = 1, sD
   
     ! Mass component -------------------------------------------------------------------------------
     IF (i == 1)  Int(1) =   IuT(nu+4,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c)  &
                         + 2*IuT(nu+2,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  &
			 +   IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c) * Ixi(4,c)

     IF (i == 2)  Int(1) = Int(1) +   IuT(nu+0,1,c) * IuT(nv+4,2,c) * IuT(nw+0,3,c)  &
                                  + 2*IuT(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw+0,3,c)  & 
                                  + 2*IuT(nu+0,1,c) * IuT(nv+2,2,c) * IuT(nw+0,3,c) * Ixi(2,c) 

     IF (i == 3)  Int(1) = Int(1) +   IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+4,3,c)  &
                                  + 2*IuT(nu+0,1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c)  &
                                  + 2*IuT(nu+2,1,c) * IuT(nv+0,2,c) * IuT(nw+2,3,c)  &
                                  + 2*IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+2,3,c) * Ixi(2,c)

     ! Velocity components --------------------------------------------------------------------------
     IF (i == 1)  Int(2) =   IuT(nu+5,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c)  &
                         + 2*IuT(nu+3,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  &
			 +   IuT(nu+1,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c) * Ixi(4,c)

     IF (i == 2) THEN
     
       Int(2) = Int(2) +   IuT(nu+1,1,c) * IuT(nv+4,2,c) * IuT(nw+0,3,c)  &
        	       + 2*IuT(nu+3,1,c) * IuT(nv+2,2,c) * IuT(nw+0,3,c)  & 
        	       + 2*IuT(nu+1,1,c) * IuT(nv+2,2,c) * IuT(nw+0,3,c) * Ixi(2,c)
		       
       Int(3) =   IuT(nu+4,1,c) * IuT(nv+1,2,c) * IuT(nw+0,3,c)  &
              +   IuT(nu+0,1,c) * IuT(nv+5,2,c) * IuT(nw+0,3,c)  &
	      +   IuT(nu+0,1,c) * IuT(nv+1,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  &
              + 2*IuT(nu+2,1,c) * IuT(nv+1,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  & 
              + 2*IuT(nu+0,1,c) * IuT(nv+3,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  &
	      + 2*IuT(nu+2,1,c) * IuT(nv+3,2,c) * IuT(nw+0,3,c)
		       
     ENDIF


     IF (i == 3) THEN
      
       Int(2) = Int(2) +   IuT(nu+1,1,c) * IuT(nv+0,2,c) * IuT(nw+4,3,c)  &
                       + 2*IuT(nu+1,1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c)  &
                       + 2*IuT(nu+3,1,c) * IuT(nv+0,2,c) * IuT(nw+2,3,c)  &
                       + 2*IuT(nu+1,1,c) * IuT(nv+0,2,c) * IuT(nw+2,3,c) * Ixi(2,c)

       Int(3) = Int(3) +   IuT(nu+0,1,c) * IuT(nv+1,2,c) * IuT(nw+4,3,c)  &
                       + 2*IuT(nu+0,1,c) * IuT(nv+1,2,c) * IuT(nw+2,3,c) * Ixi(2,c)  & 
                       + 2*IuT(nu+0,1,c) * IuT(nv+3,2,c) * IuT(nw+2,3,c)  &
	               + 2*IuT(nu+2,1,c) * IuT(nv+1,2,c) * IuT(nw+2,3,c)

       Int(4) =   IuT(nu+4,1,c) * IuT(nv+0,2,c) * IuT(nw+1,3,c)  &
              +   IuT(nu+0,1,c) * IuT(nv+4,2,c) * IuT(nw+1,3,c)  &
	      +   IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+5,3,c)  &
	      +   IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+1,3,c) * Ixi(2,c)  &
              + 2*IuT(nu+2,1,c) * IuT(nv+0,2,c) * IuT(nw+1,3,c) * Ixi(2,c)  & 
              + 2*IuT(nu+0,1,c) * IuT(nv+2,2,c) * IuT(nw+1,3,c) * Ixi(2,c)  &
	      + 2*IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+3,3,c) * Ixi(2,c)  &
	      + 2*IuT(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw+1,3,c)  &
	      + 2*IuT(nu+0,1,c) * IuT(nv+2,2,c) * IuT(nw+3,3,c)	&
	      + 2*IuT(nu+2,1,c) * IuT(nv+0,2,c) * IuT(nw+3,3,c)
     ENDIF
     
     ! Energy components ---------------------------------------------------------------------------
     IF (i == 1)  Int(sD+2) =   IuT(nu+6,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c)  &
                            + 3*IuT(nu+2,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c) * Ixi(4,c)  &
                            + 3*IuT(nu+4,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  &
			    +   IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+0,3,c) * Ixi(6,c)

     IF (i == 2)  Int(sD+2) = Int(sD+2) +   IuT(nu+0,1,c) * IuT(nv+6,2,c) * IuT(nw+0,3,c)  &
                                        + 3*IuT(nu+0,1,c) * IuT(nv+2,2,c) * IuT(nw+0,3,c) * Ixi(4,c)  &
                                        + 3*IuT(nu+0,1,c) * IuT(nv+4,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  &
                                        + 3*IuT(nu+2,1,c) * IuT(nv+4,2,c) * IuT(nw+0,3,c)  &
					+ 3*IuT(nu+4,1,c) * IuT(nv+2,2,c) * IuT(nw+0,3,c)

     IF (i == 3)  Int(sD+2) = Int(sD+2) +   IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+6,3,c)  &
                                        + 3*IuT(nu+4,1,c) * IuT(nv+0,2,c) * IuT(nw+2,3,c)  &
                                        + 3*IuT(nu+0,1,c) * IuT(nv+4,2,c) * IuT(nw+2,3,c)  &
                                        + 3*IuT(nu+0,1,c) * IuT(nv+2,2,c) * IuT(nw+4,3,c)  &
                                        + 3*IuT(nu+2,1,c) * IuT(nv+0,2,c) * IuT(nw+4,3,c)  &
                                        + 3*IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+4,3,c) * Ixi(2,c)  &
					+ 3*IuT(nu+0,1,c) * IuT(nv+0,2,c) * IuT(nw+2,3,c) * Ixi(4,c)  &
					+ 6*IuT(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c)  &
					+ 6*IuT(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw+0,3,c) * Ixi(2,c)  &
					+ 6*IuT(nu+2,1,c) * IuT(nv+0,2,c) * IuT(nw+2,3,c) * Ixi(2,c)  &
					+ 6*IuT(nu+0,1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c) * Ixi(2,c)
   ENDDO 

   Int = 0.25*Int
   Int(sD+2) = 0.5*Int(sD+2)
   
   END FUNCTION  I__c_xi2_PSI_g
   
   
   
   
   
   
      
   !----------------------------------------------------------------------------
   ! MOMENTS OF THE MAXWELLIAN: MATRIX INTEGRALS
   ! (based on diadic product PSI^(T) x PSI)
   !----------------------------------------------------------------------------

   FUNCTION  I__c_PSIxPSI_g(n, state)   RESULT(IInt)

   ! Input vector n contains the order of moment with respect
   ! the three MICROscopic velocity components. Let us define
   !
   !  nu = n(1);   nv = n(2);	nw = n(3)
   !
   ! The moment reads
   !
   !	_+oo
   !  _/  u^nu v^nv w^nw  PSI^(T) x PSI  g dudvdwdxi
   !   -oo
   !
   ! The resulting matrix is symmetric. Below are reported
   ! only the distinct elements of the matrix.
   !
   !			 _+oo
   !       IInt(1,1) = _/  u^nu v^nv w^nw g dudvdwdxi = <u^nu><v^nv><w^nw>
   !			-oo
   !
   !			 _+oo
   !       IInt(1,2) = _/  u^(nu+1) v^nv w^nw g dudvdwdxi = <u^(nu+1)><v^nv><w^nw>
   !			-oo
   !
   !			 _+oo
   !       IInt(1,3) = _/  u^nu v^(nv+1) w^nw g dudvdwdxi = <u^nu><v^(nv+1)><w^nw>
   !			-oo
   !
   !			 _+oo
   !       IInt(1,4) = _/  u^nu v^nv w^(nw+1) g dudvdwdxi = <u^nu><v^nv><w^(nw+1)>
   !			-oo
   !
   !			 _+oo
   !       IInt(1,5) = _/  u^nu v^nv w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g dudvdwdxi = 
   !			-oo
   !		     
   !		     = 0.5 * (<u^(nu+2)><v^nv><w^nw> + <u^nu><v^(nv+2)><w^nw> + 
   !			      <u^nu><v^nv><w^(nw+2)> + <u^nu><v^nv><w^nw><xi^2>)
   !
   !
   !
   !			 _+oo
   !       IInt(2,2) = _/  u^(nu+2) v^nv w^nw g dudvdwdxi = <u^(nu+2)><v^nv><w^nw>
   !			-oo
   !
   !			 _+oo
   !       IInt(2,3) = _/  u^(nu+1) v^(nv+1) w^nw g dudvdwdxi = <u^(nu+1)><v^(nv+1)><w^nw>
   !			-oo
   !			 _+oo
   !       IInt(2,4) = _/  u^(nu+1) v^nv w^(nw+1) g dudvdwdxi = <u^(nu+1)><v^nv><w^(nw+1)>
   !			-oo
   !
   !			 _+oo
   !      IInt(2,5) = _/  u^(nu+1) v^nv w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g dudvdwdxi = 
   !			-oo
   !
   !		     = 0.5 * (<u^(nu+3)><v^nv><w^nw> + <u^(nu+1)><v^(nv+2)><w^nw> + 
   !			      <u^(nu+1)><v^nv><w^(nw+2)> + <u^(nu+1)><v^nv><w^nw><xi^2>)
   !
   !
   !
   !			 _+oo
   !       IInt(3,3) = _/  u^nu v^(nv+2) w^nw g dudvdwdxi = <u^nu><v^(nv+2)><w^nw>
   !			-oo
   !
   !			 _+oo
   !       IInt(3,4) = _/  u^nu v^(nv+1) w^(nw+1) g dudvdwdxi = <u^nu><v^(nv+1)><w^(nw+1)>
   !			-oo
   !
   !			 _+oo
   !       IInt(3,5) = _/  u^nu v^(nv+1) w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g dudvdwdxi = 
   !			-oo
   !
   !		     = 0.5 * (<u^(nu+2)><v^(nv+1)><w^nw> + <u^nu><v^(nv+3)><w^nw> + 
   !			      <u^nu><v^(nv+1)><w^(nw+2)> + <u^nu><v^(nv+1)><w^nw><xi^2>)
   !
   !
   !
   !			 _+oo
   !       IInt(4,4) = _/  u^nu v^nv w^(nw+2) g dudvdwdxi = <u^nu><v^nv><w^(nw+2)>
   !			-oo
   !
   !			 _+oo
   !       IInt(4,5) = _/  u^nu v^nv w^(nw+1) 0.5*(u^2 + v^2 + w^2 + xi^2) g dudvdwdxi = 
   !			-oo
   !
   !		     = 0.5 * (<u^(nu+2)><v^nv><w^(nw+1)> + <u^nu><v^(nv+2)><w^(nw+1)> + 
   !			      <u^nu><v^nv><w^(nw+3)> + <u^nu><v^nv><w^(nw+1)><xi^2>)
   !
   !
   !
   !			 _+oo
   !       IInt(5,5) = _/  u^nu v^nv w^nw 0.25*(u^2 + v^2 + w^2 + xi^2)^2 g dudvdwdxi =
   !			-oo
   !
   !		     = 0.25 * (<u^(nu+4)><v^nv><w^nw> + <u^nu><v^(nv+4)><w^nw>  + 
   !			      <u^nu><v^nv><w^(nw+4)> + <u^nu><v^nv><w^nw><xi^4> +
   !			      2*<u^(nu+2)><v^(nv+2)><w^nw> + 2*<u^(nu+2)><v^nv><w^(nw+2)> +
   !			      2*<u^(nu+2)><v^nv><w^nw><xi^2> + 2*<u^nu><v^(nv+2)><w^(nw+2)> +
   !			      2*<u^nu><v^(nv+2)><w^nw><xi^2> + 2*<u^nu><v^nv><w^(nw+2)><xi^2>)
   !--------------------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,          DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),               INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2,sD+2) :: IInt

   INTEGER :: c, nu, nv, nw
   !--------------------------------------------------------------------------------------------

   nu = n(1)
   nv = n(2)
   nw = n(3)

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3


   SELECT CASE (sD)

     CASE (1)	! ONE SPACE DIMENSION

       ! Assembling matrix of moments
       IInt(1,1) = IuT(nu  ,1,c)
       IInt(1,2) = IuT(nu+1,1,c)
       IInt(1,3) = 0.5 * (IuT(nu+2,1,c) + IuT(nu,  1,c)*Ixi(2,c))

       IInt(2,1) = IInt(1,2)  
       IInt(2,2) = IuT(nu+2,1,c)
       IInt(2,3) = 0.5 * (IuT(nu+3,1,c) + IuT(nu+1,1,c)*Ixi(2,c))
       
       IInt(3,1) = IInt(1,3)  
       IInt(3,2) = IInt(2,3)  
       IInt(3,3) = 0.25 * (IuT(nu+4,1,c) + IuT(nu,  1,c)*Ixi(4,c) + 2*IuT(nu+2,1,c)*Ixi(2,c))



     CASE (2)	! TWO SPACE DIMENSIONS

       IInt(1,1) = IuT(nu,  1,c) * IuT(nv,  2,c)
       IInt(1,2) = IuT(nu+1,1,c) * IuT(nv,  2,c)
       IInt(1,3) = IuT(nu,  1,c) * IuT(nv+1,2,c)
       IInt(1,4) = 0.50 * (IuT(nu+2,1,c) * IuT(nv,  2,c) &
   			  +IuT(nu,  1,c) * IuT(nv+2,2,c) &
   			  +IuT(nu,  1,c) * IuT(nv,  2,c) * Ixi(2,c))

       IInt(2,1) = IInt(1,2)
       IInt(2,2) = IuT(nu+2,1,c) * IuT(nv,  2,c)
       IInt(2,3) = IuT(nu+1,1,c) * IuT(nv+1,2,c)
       IInt(2,4) = 0.50 * (IuT(nu+3,1,c) * IuT(nv,  2,c) +  &
   			   IuT(nu+1,1,c) * IuT(nv+2,2,c) +  &
   			   IuT(nu+1,1,c) * IuT(nv,  2,c) * Ixi(2,c))

       IInt(3,1) = IInt(1,3)
       IInt(3,2) = IInt(2,3)
       IInt(3,3) = IuT(nu,  1,c) * IuT(nv+2,2,c)
       IInt(3,4) = 0.50 * (IuT(nu+2,1,c) * IuT(nv+1,2,c) +  &
   			   IuT(nu,  1,c) * IuT(nv+3,2,c) +  &
   			   IuT(nu,  1,c) * IuT(nv+1,2,c) * Ixi(2,c))

       IInt(4,1) = IInt(1,4)
       IInt(4,2) = IInt(2,4)  
       IInt(4,3) = IInt(3,4)        
       IInt(4,4) = 0.25 * (IuT(nu+4,1,c) * IuT(nv,  2,c) +  &
   		           IuT(nu,  1,c) * IuT(nv+4,2,c) +  &
   		       2 * IuT(nu+2,1,c) * IuT(nv+2,2,c) +  &
   		           IuT(nu,  1,c) * IuT(nv,  2,c) * Ixi(4,c) +  &
   		       2 * IuT(nu+2,1,c) * IuT(nv,  2,c) * Ixi(2,c) +  &
   		       2 * IuT(nu,  1,c) * IuT(nv+2,2,c) * Ixi(2,c))



     CASE (3)	! THREE SPACE DIMENSIONS

       IInt(1,1) =	   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IInt(1,2) =	   IuT(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IInt(1,3) =	   IuT(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
       IInt(1,4) =	   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)
       IInt(1,5) = 0.50 * (IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(2,c))
       
       IInt(2,1) = IInt(1,2)
       IInt(2,2) =	   IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IInt(2,3) =	   IuT(nu+1,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
       IInt(2,4) =	   IuT(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)
       IInt(2,5) = 0.50 * (IuT(nu+3,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
   	  		   IuT(nu+1,1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
   	  		   IuT(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
   	  		   IuT(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,3,c) * Ixi(2,c))
       
       IInt(3,1) = IInt(1,3)
       IInt(3,2) = IInt(2,3)
       IInt(3,3) =	   IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c)
       IInt(3,4) =	   IuT(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw+1,3,c)
       IInt(3,5) = 0.50 * (IuT(nu+2,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv+3,2,c) * IuT(nw,  3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw+2,3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c) * Ixi(2,c))
       
       IInt(4,1) = IInt(1,4)
       IInt(4,2) = IInt(2,4)
       IInt(4,3) = IInt(3,4)
       IInt(4,4) =	   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c)
       IInt(4,5) = 0.50 * (IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+1,3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+3,3,c) +  &
   	  		   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c) * Ixi(2,c))

       IInt(5,1) = IInt(1,5)
       IInt(5,2) = IInt(2,5)
       IInt(5,3) = IInt(3,5)
       IInt(5,4) = IInt(4,5)
       IInt(5,5) = 0.25 * (IuT(nu+4,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
   			   IuT(nu,  1,c) * IuT(nv+4,2,c) * IuT(nw,  3,c) +  &
   			   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+4,3,c) +  &
   			   IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(4,c) +  &
   		       2 * IuT(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
   		       2 * IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
   		       2 * IuT(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(2,c) +  &
   		       2 * IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c) +  & 
   		       2 * IuT(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) * Ixi(2,c) +  &
   		       2 * IuT(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) * Ixi(2,c))

   END SELECT
 
   END FUNCTION  I__c_PSIxPSI_g







   FUNCTION  IP__c_PSIxPSI_g(n, state)   RESULT(IIntP)

   ! Input vector n contains the order of moment with respect
   ! the three MICROscopic velocity components. Let us define
   !
   !  nu = n(1);   nv = n(2);	nw = n(3)
   !
   ! The moment reads
   !
   !	_+oo _+oo
   !  _/   _/  u^nu v^nv w^nw  PSI^(T) x PSI  g du dvdwdxi
   !   -oo  0
   !
   ! IMPORTANT: Note the extremes of integration
   !
   ! The resulting matrix is symmetric. Below are reported
   ! only the distinct elements of the matrix.
   !
   !			  _+oo _+oo
   !   IIntP(1,1) = _/   _/  u^nu v^nv w^nw g du dvdwdxi = <u^nu><v^nv><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntP(1,2) = _/   _/  u^(nu+1) v^nv w^nw g du dvdwdxi = <u^(nu+1)><v^nv><w^nw>
   !			-oo   0
   !
   !			  _+oo
   !   IIntP(1,3) = _/  u^nu v^(nv+1) w^nw g du dvdwdxi = <u^nu><v^(nv+1)><w^nw>
   !			-oo
   !
   !			  _+oo _+oo
   !   IIntP(1,4) = _/   _/  u^nu v^nv w^(nw+1) g du dvdwdxi = <u^nu><v^nv><w^(nw+1)>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntP(1,5) = _/   _/  u^nu v^nv w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !		     
   !		     = 0.5 * (<u^(nu+2)><v^nv><w^nw> + <u^nu><v^(nv+2)><w^nw> + 
   !			      <u^nu><v^nv><w^(nw+2)> + <u^nu><v^nv><w^nw><xi^2>)
   !
   !
   !
   !			  _+oo _+oo
   !   IIntP(2,2) = _/   _/  u^(nu+2) v^nv w^nw g du dvdwdxi = <u^(nu+2)><v^nv><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntP(2,3) = _/   _/  u^(nu+1) v^(nv+1) w^nw g du dvdwdxi = <u^(nu+1)><v^(nv+1)><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntP(2,4) = _/   _/  u^(nu+1) v^nv w^(nw+1) g du dvdwdxi = <u^(nu+1)><v^nv><w^(nw+1)>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntP(2,5) = _/   _/  u^(nu+1) v^nv w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !
   !		     = 0.5 * (<u^(nu+3)><v^nv><w^nw> + <u^(nu+1)><v^(nv+2)><w^nw> + 
   !			      <u^(nu+1)><v^nv><w^(nw+2)> + <u^(nu+1)><v^nv><w^nw><xi^2>)
   !
   !
   !
   !			  _+oo _+oo
   !   IIntP(3,3) = _/   _/  u^nu v^(nv+2) w^nw g du dvdwdxi = <u^nu><v^(nv+2)><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntP(3,4) = _/   _/  u^nu v^(nv+1) w^(nw+1) g du dvdwdxi = <u^nu><v^(nv+1)><w^(nw+1)>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntP(3,5) = _/   _/  u^nu v^(nv+1) w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !
   !		     = 0.5 * (<u^(nu+2)><v^(nv+1)><w^nw> + <u^nu><v^(nv+3)><w^nw> + 
   !			      <u^nu><v^(nv+1)><w^(nw+2)> + <u^nu><v^(nv+1)><w^nw><xi^2>)
   !
   !
   !
   !			  _+oo
   !   IIntP(4,4) = _/  u^nu v^nv w^(nw+2) g du dvdwdxi = <u^nu><v^nv><w^(nw+2)>
   !			-oo
   !
   !			  _+oo _+oo
   !   IIntP(4,5) = _/   _/  u^nu v^nv w^(nw+1) 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !
   !		     = 0.5 * (<u^(nu+2)><v^nv><w^(nw+1)> + <u^nu><v^(nv+2)><w^(nw+1)> + 
   !			      <u^nu><v^nv><w^(nw+3)> + <u^nu><v^nv><w^(nw+1)><xi^2>)
   !
   !
   !
   !			  _+oo _+oo
   !   IIntP(5,5) = _/   _/  u^nu v^nv w^nw 0.25*(u^2 + v^2 + w^2 + xi^2)^2 g du dvdwdxi =
   !			-oo   0
   !
   !		     = 0.25 * (<u^(nu+4)><v^nv><w^nw> + <u^nu><v^(nv+4)><w^nw>  + 
   !			      <u^nu><v^nv><w^(nw+4)> + <u^nu><v^nv><w^nw><xi^4> +
   !			      2*<u^(nu+2)><v^(nv+2)><w^nw> + 2*<u^(nu+2)><v^nv><w^(nw+2)> +
   !			      2*<u^(nu+2)><v^nv><w^nw><xi^2> + 2*<u^nu><v^(nv+2)><w^(nw+2)> +
   !			      2*<u^nu><v^(nv+2)><w^nw><xi^2> + 2*<u^nu><v^nv><w^(nw+2)><xi^2>)
   !--------------------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,          DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),               INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2,sD+2) :: IIntP

   INTEGER :: c, nu, nv, nw
   !--------------------------------------------------------------------------------------------

   nu = n(1)
   nv = n(2)
   nw = n(3)

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3


   SELECT CASE (sD)

     CASE (1)	! ONE SPACE DIMENSION

       IIntP(1,1) =	    IuP(nu  ,1,c)
       IIntP(1,2) =	    IuP(nu+1,1,c)
       IIntP(1,3) = 0.50 * (IuP(nu+2,1,c) + IuP(nu,  1,c)*Ixi(2,c))
       
       IIntP(2,1) = IIntP(1,2)
       IIntP(2,2) =	    IuP(nu+2,1,c)
       IIntP(2,3) = 0.50 * (IuP(nu+3,1,c) + IuP(nu+1,1,c)*Ixi(2,c))

       IIntP(3,1) = IIntP(1,3)  
       IIntP(3,2) = IIntP(2,3)  
       IIntP(3,3) = 0.25 * (IuP(nu+4,1,c) + IuP(nu,  1,c)*Ixi(4,c) + 2*IuP(nu+2,1,c)*Ixi(2,c))



     CASE (2)	! TWO SPACE DIMENSIONS

       IIntP(1,1) =	    IuP(nu,  1,c) * IuT(nv,  2,c)
       IIntP(1,2) =	    IuP(nu+1,1,c) * IuT(nv,  2,c)
       IIntP(1,3) =	    IuP(nu,  1,c) * IuT(nv+1,2,c)
       IIntP(1,4) = 0.50 * (IuP(nu+2,1,c) * IuT(nv,  2,c) &
   				+IuP(nu,  1,c) * IuT(nv+2,2,c) &
   				+IuP(nu,  1,c) * IuT(nv,  2,c) * Ixi(2,c))
       
       IIntP(2,1) = IIntP(1,2)
       IIntP(2,2) =	    IuP(nu+2,1,c) * IuT(nv,  2,c)
       IIntP(2,3) =	    IuP(nu+1,1,c) * IuT(nv+1,2,c)
       IIntP(2,4) = 0.50 * (IuP(nu+3,1,c) * IuT(nv,  2,c) +  &
   				 IuP(nu+1,1,c) * IuT(nv+2,2,c) +  &
   				 IuP(nu+1,1,c) * IuT(nv,  2,c) * Ixi(2,c))
       
       IIntP(3,1) = IIntP(1,3)
       IIntP(3,2) = IIntP(2,3)
       IIntP(3,3) =	    IuP(nu,  1,c) * IuT(nv+2,2,c)
       IIntP(3,4) = 0.50 * (IuP(nu+2,1,c) * IuT(nv+1,2,c) +  &
   				 IuP(nu,  1,c) * IuT(nv+3,2,c) +  &
   				 IuP(nu,  1,c) * IuT(nv+1,2,c) * Ixi(2,c))

       IIntP(4,1) = IIntP(1,4)
       IIntP(4,2) = IIntP(2,4)  
       IIntP(4,3) = IIntP(3,4)  
       IIntP(4,4) = 0.25 * (IuP(nu+4,1,c) * IuT(nv,  2,c) +  &
   				 IuP(nu,  1,c) * IuT(nv+4,2,c) +  &
   			     2 * IuP(nu+2,1,c) * IuT(nv+2,2,c) +  &
   				 IuP(nu,  1,c) * IuT(nv,  2,c) * Ixi(4,c) +  &
   			     2 * IuP(nu+2,1,c) * IuT(nv,  2,c) * Ixi(2,c) +  &
   			     2 * IuP(nu,  1,c) * IuT(nv+2,2,c) * Ixi(2,c))



     CASE (3)	! THREE SPACE DIMENSIONS

       IIntP(1,1) =	    IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IIntP(1,2) =	    IuP(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IIntP(1,3) =	    IuP(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
       IIntP(1,4) =	    IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)
       IIntP(1,5) = 0.50 * (IuP(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
 			    IuP(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
 			    IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
 			    IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(2,c))
       
       IIntP(2,1) = IIntP(1,2)
       IIntP(2,2) =	    IuP(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IIntP(2,3) =	    IuP(nu+1,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
       IIntP(2,4) =	    IuP(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)
       IIntP(2,5) = 0.50 * (IuP(nu+3,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
   	   		    IuP(nu+1,1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
   	   		    IuP(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
   	   		    IuP(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,3,c) * Ixi(2,c))
       
       IIntP(3,1) = IIntP(1,3)
       IIntP(3,2) = IIntP(2,3)
       IIntP(3,3) =	    IuP(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c)
       IIntP(3,4) =	    IuP(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw+1,3,c)
       IIntP(3,5) = 0.50 * (IuP(nu+2,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c) +  &
   	   		    IuP(nu,  1,c) * IuT(nv+3,2,c) * IuT(nw,  3,c) +  &
   	   		    IuP(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw+2,3,c) +  &
   	   		    IuP(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c) * Ixi(2,c))
       
       IIntP(4,1) = IIntP(1,4)
       IIntP(4,2) = IIntP(2,4)
       IIntP(4,3) = IIntP(3,4)
       IIntP(4,4) =	    IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c)
       IIntP(4,5) = 0.50 * (IuP(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c) +  &
   	   		    IuP(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+1,3,c) +  &
   	   		    IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+3,3,c) +  &
   	   		    IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c) * Ixi(2,c))

       IIntP(5,1) = IIntP(1,5)
       IIntP(5,2) = IIntP(2,5)
       IIntP(5,3) = IIntP(3,5)
       IIntP(5,4) = IIntP(4,5)
       IIntP(5,5) = 0.25 * (IuP(nu+4,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
   				 IuP(nu,  1,c) * IuT(nv+4,2,c) * IuT(nw,  3,c) +  &
   				 IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+4,3,c) +  &
   				 IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(4,c) +  &
   			     2 * IuP(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
   			     2 * IuP(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
   			     2 * IuP(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(2,c) +  &
   			     2 * IuP(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c) +  & 
   			     2 * IuP(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) * Ixi(2,c) +  &
   			     2 * IuP(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) * Ixi(2,c))

   END SELECT

   END FUNCTION  IP__c_PSIxPSI_g







   FUNCTION  IN__c_PSIxPSI_g(n, state)   RESULT(IIntN)

   ! Input vector n contains the order of moment with respect
   ! the three MICROscopic velocity components. Let us define
   !
   !  nu = n(1);   nv = n(2);	nw = n(3)
   !
   ! The moment reads
   !
   !	_+oo _+oo
   !  _/   _/  u^nu v^nv w^nw  PSI^(T) x PSI  g du dvdwdxi
   !   -oo  0
   !
   ! IMPORTANT: Note the extremes of integration
   !
   ! The resulting matrix is symmetric. Below are reported
   ! only the distinct elements of the matrix.
   !
   !			  _+oo _+oo
   !   IIntN(1,1) = _/   _/  u^nu v^nv w^nw g du dvdwdxi = <u^nu><v^nv><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntN(1,2) = _/   _/  u^(nu+1) v^nv w^nw g du dvdwdxi = <u^(nu+1)><v^nv><w^nw>
   !			-oo   0
   !
   !			  _+oo
   !   IIntN(1,3) = _/  u^nu v^(nv+1) w^nw g du dvdwdxi = <u^nu><v^(nv+1)><w^nw>
   !			-oo
   !
   !			  _+oo _+oo
   !   IIntN(1,4) = _/   _/  u^nu v^nv w^(nw+1) g du dvdwdxi = <u^nu><v^nv><w^(nw+1)>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntN(1,5) = _/   _/  u^nu v^nv w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !		     
   !		     = 0.5 * (<u^(nu+2)><v^nv><w^nw> + <u^nu><v^(nv+2)><w^nw> + 
   !			      <u^nu><v^nv><w^(nw+2)> + <u^nu><v^nv><w^nw><xi^2>)
   !
   !
   !
   !			  _+oo _+oo
   !   IIntN(2,2) = _/   _/  u^(nu+2) v^nv w^nw g du dvdwdxi = <u^(nu+2)><v^nv><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntN(2,3) = _/   _/  u^(nu+1) v^(nv+1) w^nw g du dvdwdxi = <u^(nu+1)><v^(nv+1)><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntN(2,4) = _/   _/  u^(nu+1) v^nv w^(nw+1) g du dvdwdxi = <u^(nu+1)><v^nv><w^(nw+1)>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntN(2,5) = _/   _/  u^(nu+1) v^nv w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !
   !		     = 0.5 * (<u^(nu+3)><v^nv><w^nw> + <u^(nu+1)><v^(nv+2)><w^nw> + 
   !			      <u^(nu+1)><v^nv><w^(nw+2)> + <u^(nu+1)><v^nv><w^nw><xi^2>)
   !
   !
   !
   !			  _+oo _+oo
   !   IIntN(3,3) = _/   _/  u^nu v^(nv+2) w^nw g du dvdwdxi = <u^nu><v^(nv+2)><w^nw>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntN(3,4) = _/   _/  u^nu v^(nv+1) w^(nw+1) g du dvdwdxi = <u^nu><v^(nv+1)><w^(nw+1)>
   !			-oo   0
   !
   !			  _+oo _+oo
   !   IIntN(3,5) = _/   _/  u^nu v^(nv+1) w^nw 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !
   !		     = 0.5 * (<u^(nu+2)><v^(nv+1)><w^nw> + <u^nu><v^(nv+3)><w^nw> + 
   !			      <u^nu><v^(nv+1)><w^(nw+2)> + <u^nu><v^(nv+1)><w^nw><xi^2>)
   !
   !
   !
   !			  _+oo
   !   IIntN(4,4) = _/  u^nu v^nv w^(nw+2) g du dvdwdxi = <u^nu><v^nv><w^(nw+2)>
   !			-oo
   !
   !			  _+oo _+oo
   !   IIntN(4,5) = _/   _/  u^nu v^nv w^(nw+1) 0.5*(u^2 + v^2 + w^2 + xi^2) g du dvdwdxi = 
   !			-oo   0
   !
   !		     = 0.5 * (<u^(nu+2)><v^nv><w^(nw+1)> + <u^nu><v^(nv+2)><w^(nw+1)> + 
   !			      <u^nu><v^nv><w^(nw+3)> + <u^nu><v^nv><w^(nw+1)><xi^2>)
   !
   !
   !
   !			  _+oo _+oo
   !   IIntN(5,5) = _/   _/  u^nu v^nv w^nw 0.25*(u^2 + v^2 + w^2 + xi^2)^2 g du dvdwdxi =
   !			-oo   0
   !
   !		     = 0.25 * (<u^(nu+4)><v^nv><w^nw> + <u^nu><v^(nv+4)><w^nw>  + 
   !			      <u^nu><v^nv><w^(nw+4)> + <u^nu><v^nv><w^nw><xi^4> +
   !			      2*<u^(nu+2)><v^(nv+2)><w^nw> + 2*<u^(nu+2)><v^nv><w^(nw+2)> +
   !			      2*<u^(nu+2)><v^nv><w^nw><xi^2> + 2*<u^nu><v^(nv+2)><w^(nw+2)> +
   !			      2*<u^nu><v^(nv+2)><w^nw><xi^2> + 2*<u^nu><v^nv><w^(nw+2)><xi^2>)
   !--------------------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,          DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),               INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2,sD+2) :: IIntN

   INTEGER :: c, nu, nv, nw
   !--------------------------------------------------------------------------------------------
   
   nu = n(1)
   nv = n(2)
   nw = n(3)

   IF (state == 'L')  c = 1
   IF (state == 'I')  c = 2
   IF (state == 'R')  c = 3


   SELECT CASE (sD)

     CASE (1)	! ONE SPACE DIMENSION

       IIntN(1,1) =	    IuN(nu  ,1,c)
       IIntN(1,2) =	    IuN(nu+1,1,c)
       IIntN(1,3) = 0.50 * (IuN(nu+2,1,c) + IuN(nu,  1,c)*Ixi(2,c))
       
       IIntN(2,1) = IIntN(1,2) 
       IIntN(2,2) =	    IuN(nu+2,1,c)
       IIntN(2,3) = 0.50 * (IuN(nu+3,1,c) + IuN(nu+1,1,c)*Ixi(2,c))

       IIntN(3,1) = IIntN(1,3)  
       IIntN(3,2) = IIntN(2,3)  
       IIntN(3,3) = 0.25 * (IuN(nu+4,1,c) + IuN(nu,  1,c)*Ixi(4,c) + 2*IuN(nu+2,1,c)*Ixi(2,c))



     CASE (2)	! TWO SPACE DIMENSIONS

       IIntN(1,1) =	    IuN(nu,  1,c) * IuT(nv,  2,c)
       IIntN(1,2) =	    IuN(nu+1,1,c) * IuT(nv,  2,c)
       IIntN(1,3) =	    IuN(nu,  1,c) * IuT(nv+1,2,c)
       IIntN(1,4) = 0.50 * (IuN(nu+2,1,c) * IuT(nv,  2,c) &
   	   		   +IuN(nu,  1,c) * IuT(nv+2,2,c) &
   	   		   +IuN(nu,  1,c) * IuT(nv,  2,c) * Ixi(2,c))
       
       IIntN(2,1) = IIntN(1,2)
       IIntN(2,2) =	    IuN(nu+2,1,c) * IuT(nv,  2,c)
       IIntN(2,3) =	    IuN(nu+1,1,c) * IuT(nv+1,2,c)
       IIntN(2,4) = 0.50 * (IuN(nu+3,1,c) * IuT(nv,  2,c) +  &
   	   		    IuN(nu+1,1,c) * IuT(nv+2,2,c) +  &
   	   		    IuN(nu+1,1,c) * IuT(nv,  2,c) * Ixi(2,c))

       IIntN(3,1) = IIntN(1,3)
       IIntN(3,2) = IIntN(2,3)
       IIntN(3,3) =	    IuN(nu,  1,c) * IuT(nv+2,2,c)
       IIntN(3,4) = 0.50 * (IuN(nu+2,1,c) * IuT(nv+1,2,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv+3,2,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv+1,2,c) * Ixi(2,c))

       IIntN(4,1) = IIntN(1,4)
       IIntN(4,2) = IIntN(2,4)  
       IIntN(4,3) = IIntN(3,4)  
       IIntN(4,4) = 0.25 * (IuN(nu+4,1,c) * IuT(nv,  2,c) +  &
   		   	    IuN(nu,  1,c) * IuT(nv+4,2,c) +  &
 			2 * IuN(nu+2,1,c) * IuT(nv+2,2,c) +  &
 			    IuN(nu,  1,c) * IuT(nv,  2,c) * Ixi(4,c) +  &
 			2 * IuN(nu+2,1,c) * IuT(nv,  2,c) * Ixi(2,c) +  &
 			2 * IuN(nu,  1,c) * IuT(nv+2,2,c) * Ixi(2,c))



     CASE (3)	! THREE SPACE DIMENSIONS

       IIntN(1,1) =	    IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IIntN(1,2) =	    IuN(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IIntN(1,3) =	    IuN(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
       IIntN(1,4) =	    IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)
       IIntN(1,5) = 0.50 * (IuN(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(2,c))
       
       IIntN(2,1) = IIntN(1,2)
       IIntN(2,2) =	    IuN(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c)
       IIntN(2,3) =	    IuN(nu+1,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c)
       IIntN(2,4) =	    IuN(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c)
       IIntN(2,5) = 0.50 * (IuN(nu+3,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
   	   		    IuN(nu+1,1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
   	   		    IuN(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
   	   		    IuN(nu+1,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(2,c))

       IIntN(3,1) = IIntN(1,3)
       IIntN(3,2) = IIntN(2,3)
       IIntN(3,3) =	    IuN(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c)
       IIntN(3,4) =	    IuN(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw+1,3,c)
       IIntN(3,5) = 0.50 * (IuN(nu+2,1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv+3,2,c) * IuT(nw,  3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw+2,3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv+1,2,c) * IuT(nw,  3,c) * Ixi(2,c))
       
       IIntN(4,1) = IIntN(1,4)
       IIntN(4,2) = IIntN(2,4)
       IIntN(4,3) = IIntN(3,4)
       IIntN(4,4) =	    IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c)
       IIntN(4,5) = 0.50 * (IuN(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+1,3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+3,3,c) +  &
   	   		    IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+1,3,c) * Ixi(2,c))

       IIntN(5,1) = IIntN(1,5)
       IIntN(5,2) = IIntN(2,5)
       IIntN(5,3) = IIntN(3,5)
       IIntN(5,4) = IIntN(4,5)
       IIntN(5,5) = 0.25 * (IuN(nu+4,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) +  &
 		            IuN(nu,  1,c) * IuT(nv+4,2,c) * IuT(nw,  3,c) +  &
 		            IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+4,3,c) +  &
 		            IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(4,c) +  &
 		        2 * IuN(nu+2,1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) +  &
 		        2 * IuN(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) +  &
 		        2 * IuN(nu+2,1,c) * IuT(nv,  2,c) * IuT(nw,  3,c) * Ixi(2,c) +  &
 		        2 * IuN(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw+2,3,c) +  &
 		        2 * IuN(nu,  1,c) * IuT(nv+2,2,c) * IuT(nw,  3,c) * Ixi(2,c) +  &
 		        2 * IuN(nu,  1,c) * IuT(nv,  2,c) * IuT(nw+2,3,c) * Ixi(2,c))

   END SELECT

   END FUNCTION  IN__c_PSIxPSI_g







   FUNCTION  I__c_PSI_PSIxPSI_g(n,state)   RESULT(IIInt)
   !--------------------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,          DIMENSION(:), INTENT(IN) :: n
   CHARACTER(LEN=1),               INTENT(IN) :: state
   
   REAL(KIND=8), DIMENSION(sD+2,sD+2,sD+2) :: IIInt

   INTEGER :: nu, nv, nw
   !--------------------------------------------------------------------------------------------   
   
   nu = n(1)
   nv = n(2)
   nw = n(3)

   SELECT CASE (sD)
   
     CASE (1)
     
       IIInt(1,1,:) = I__c_PSI_g((/nu,  nv,nw/),state)
       IIInt(1,2,:) = I__c_PSI_g((/nu+1,nv,nw/),state)
       IIInt(1,3,:) = I__c_xi_PSI_g((/nu,nv,nw/),state)

       IIInt(2,1,:) = IIInt(1,2,:)
       IIInt(2,2,:) = I__c_PSI_g((/nu+2,nv,nw/),state)
       IIInt(2,3,:) = I__c_xi_PSI_g((/nu+1,nv,nw/),state)

       IIInt(3,1,:) = IIInt(1,3,:)
       IIInt(3,2,:) = IIInt(2,3,:)
       IIInt(3,3,:) = I__c_xi2_PSI_g((/nu,nv,nw/),state)


     CASE (2)
     
       IIInt(1,1,:) = I__c_PSI_g((/nu,  nv,  nw/),state)
       IIInt(1,2,:) = I__c_PSI_g((/nu+1,nv,  nw/),state)
       IIInt(1,3,:) = I__c_PSI_g((/nu,  nv+1,nw/),state)
       IIInt(1,4,:) = I__c_xi_PSI_g((/nu,nv,nw/),state)

       IIInt(2,1,:) = IIInt(1,2,:)
       IIInt(2,2,:) = I__c_PSI_g((/nu+2,nv,  nw/),state)
       IIInt(2,3,:) = I__c_PSI_g((/nu+1,nv+1,nw/),state)       
       IIInt(2,4,:) = I__c_xi_PSI_g((/nu+1,nv,nw/),state)

       IIInt(3,1,:) = IIInt(1,3,:)
       IIInt(3,2,:) = IIInt(2,3,:)
       IIInt(3,3,:) = I__c_PSI_g((/nu,nv+2,nw/),state)       
       IIInt(3,4,:) = I__c_xi2_PSI_g((/nu,nv+1,nw/),state)

       IIInt(4,1,:) = IIInt(1,4,:)
       IIInt(4,2,:) = IIInt(2,4,:)
       IIInt(4,3,:) = IIInt(3,4,:)
       IIInt(4,4,:) = I__c_xi2_PSI_g((/nu,nv,nw/),state)


     CASE (3)
     
       IIInt(1,1,:) = I__c_PSI_g((/nu,  nv,  nw  /),state)
       IIInt(1,2,:) = I__c_PSI_g((/nu+1,nv,  nw  /),state)
       IIInt(1,3,:) = I__c_PSI_g((/nu,  nv+1,nw  /),state)
       IIInt(1,4,:) = I__c_PSI_g((/nu,  nv,  nw+1/),state)
       IIInt(1,5,:) = I__c_xi_PSI_g((/nu,nv,nw/),state)

       IIInt(2,1,:) = IIInt(1,2,:)
       IIInt(2,2,:) = I__c_PSI_g((/nu+2,nv,  nw  /),state)
       IIInt(2,3,:) = I__c_PSI_g((/nu+1,nv+1,nw  /),state)
       IIInt(2,4,:) = I__c_PSI_g((/nu+1,nv  ,nw+1/),state)
       IIInt(2,5,:) = I__c_xi_PSI_g((/nu+1,nv,nw/),state)

       IIInt(3,1,:) = IIInt(1,3,:)
       IIInt(3,2,:) = IIInt(2,3,:)
       IIInt(3,3,:) = I__c_PSI_g((/nu,nv+2,nw  /),state)
       IIInt(3,4,:) = I__c_PSI_g((/nu,nv+1,nw+1/),state)     
       IIInt(3,5,:) = I__c_xi2_PSI_g((/nu,nv+1,nw/),state)

       IIInt(4,1,:) = IIInt(1,4,:)
       IIInt(4,2,:) = IIInt(2,4,:)
       IIInt(4,3,:) = IIInt(3,4,:)
       IIInt(4,4,:) = I__c_PSI_g((/nu,nv,nw+2/),state)       
       IIInt(4,5,:) = I__c_xi2_PSI_g((/nu,nv,nw+1/),state)

       IIInt(5,1,:) = IIInt(1,5,:)
       IIInt(5,2,:) = IIInt(2,5,:)
       IIInt(5,3,:) = IIInt(3,5,:)
       IIInt(5,4,:) = IIInt(4,5,:)   
       IIInt(5,5,:) = I__c_xi2_PSI_g((/nu,nv,nw/),state)

   END SELECT   
   
   END FUNCTION  I__c_PSI_PSIxPSI_g







   SUBROUTINE compute_DF_derivatives(ww, Gww, gF, Ft)
   !------------------------------------------------------------------------------------------   
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: ww
   REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: Gww
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(Gww,1), SIZE(ww,2)), INTENT(OUT) :: gF
   REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,2)),              INTENT(OUT) :: Ft
   
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1)*(SIZE(ww,1)+1)/2) :: pM, sysMat
   REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,1)) :: M
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: Fx, Fy, Fz, F_t
   
   REAL(KIND=8) :: rho

   INTEGER, DIMENSION(SIZE(ww,1))  :: IPIV
   INTEGER, DIMENSION(3) :: order_0, order_1   

   INTEGER :: i, p_index, p, info
   !------------------------------------------------------------------------------------------

   DO i = 1, SIZE(ww,2)

     CALL compute_scalar_moments(ww(:,i), 'I')
     
     rho = ww(1,i)

     order_0 = 0
     M = rho * I__c_PSIxPSI_g(order_0,'I')

     p_index = 0
     DO p = 1, SIZE(ww,1)
       pM(p_index+1:p_index+p) = M(1:p,p)
       p_index = p_index + p
     ENDDO   

     ! Space derivatives of F ----------------------------------------------------------------
     sysMat = pM
     Fx = Gww(1,i,:)
     CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fx, SIZE(ww,1), info)       

     IF (sD >= 2) THEN  	     
       sysMat = pM
       Fy = Gww(2,i,:)
       CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fy, SIZE(ww,1), info)    
     ENDIF;IF (sD == 3) THEN 
       sysMat = pM
       Fz = Gww(3,i,:)
       CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fz, SIZE(ww,1), info)    
     ENDIF
     !----------------------------------------------------------------------------------------
   
     ! Time derivatives of F. Compatibility on Dg --------------------------------------------
     order_1 = (/1, 0, 0/)    
     F_t = - rho * (MATMUL(I__c_PSIxPSI_g(order_1,'I'), Fx))

     IF (sD >= 2) THEN
       order_1 = (/0, 1, 0/)
       F_t = F_t - rho * (MATMUL(I__c_PSIxPSI_g(order_1,'I'), Fy))       
     ENDIF;IF (sD == 3) THEN 
       order_1 = (/0, 0, 1/)
       F_t = F_t - rho * (MATMUL(I__c_PSIxPSI_g(order_1,'I'), Fz))       
     ENDIF

     sysMat = pM
     CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, F_t, SIZE(ww,1), info) 
     !----------------------------------------------------------------------------------------


     gF(:,1,i) = Fx
     Ft(:,i) = F_t
     
     IF (sD >= 2) THEN
       gF(:,2,i) = Fy
     ENDIF;IF (sD == 3) THEN 
       gF(:,3,i) = Fz
     ENDIF


   ENDDO
	  
   END SUBROUTINE  compute_DF_derivatives





   SUBROUTINE  compute_vector_moments(ww_i, gF_i, Ft_i, I_uPSI_f, I_vPSI_f, I_wPSI_f)
   !------------------------------------------------------------------------------------------   
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww_i
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: gF_i
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: Ft_i
   
   REAL(KIND=8), DIMENSION(:) :: I_uPSI_f, I_vPSI_f, I_wPSI_f 
   
   REAL(KIND=8), DIMENSION(SIZE(ww_i)) :: Fx, Fy, Fz, Ft
   REAL(KIND=8) :: rho, tau

   INTEGER, DIMENSION(3) :: order_1, order_2   
   !------------------------------------------------------------------------------------------

   CALL compute_scalar_moments(ww_i, 'I')

   ! NOTE- to be extended for the case of transitional regime ---
   tau_fmla = 1
   tau = mean_collision_time(1.d0, ww_i, ww_i, ww_i)
   !-------------------------------------------------------------
   
   rho = ww_i(1)
   Fx  = gF_i(:,1)
   Ft  = Ft_i
   
   IF (sD >= 2)  Fy = gF_i(:,2)
   IF (sD == 3)  Fz = gF_i(:,3)



   ! I__u_PSI_f --------------------------------------------------------------
   order_1 = (/1, 0, 0/)
   order_2 = (/2, 0, 0/)

   I_uPSI_f = rho * I__c_PSI_g(order_1,'I')  &
            - rho * tau *MATMUL(I__c_PSIxPSI_g(order_1,'I'),Ft) &
            - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fx)

   IF (sD >= 2) THEN
   
     order_2 = (/1, 1, 0/)
     I_uPSI_f = I_uPSI_f - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fy)
   
   ENDIF;IF (sD == 3) THEN
   
     order_2 = (/1, 0, 1/)
     I_uPSI_f = I_uPSI_f - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fz)
   
   ENDIF



   ! I__v_PSI_f --------------------------------------------------------------
   IF (sD >= 2) THEN   
   
     order_1 = (/0, 1, 0/)
     order_2 = (/0, 2, 0/)

     I_vPSI_f = rho * I__c_PSI_g(order_1,'I')  &
     	      - rho * tau *MATMUL(I__c_PSIxPSI_g(order_1,'I'),Ft) &
     	      - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fy)

     order_2 = (/1, 1, 0/)
     
     I_vPSI_f = I_vPSI_f - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fx)

     IF (sD == 3) THEN
     
       order_2 = (/0, 1, 1/)
       I_vPSI_f = I_vPSI_f - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fz)
     
     ENDIF
     
   ENDIF
   
   
   
   ! I__w_PSI_f --------------------------------------------------------------
   IF (sD == 3) THEN     
   
     order_1 = (/0, 0, 1/)
     order_2 = (/0, 0, 2/)

     I_wPSI_f = rho * I__c_PSI_g(order_1,'I')  &
     	      - rho * tau *MATMUL(I__c_PSIxPSI_g(order_1,'I'),Ft) &
     	      - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fz)

     order_2 = (/1, 0, 1/)     
     I_wPSI_f = I_wPSI_f - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fx)

     order_2 = (/0, 1, 1/)
     I_wPSI_f = I_wPSI_f - rho * tau *MATMUL(I__c_PSIxPSI_g(order_2,'I'),Fy)

   ENDIF


   END SUBROUTINE  compute_vector_moments
   
   
   
   
   
   FUNCTION  DPstressTensor(ww, gF, Ft)   RESULT(t_ij)
   
   ! Compute the nodal deviatoric part of the stress tensor.
   !------------------------------------------------------------------------------------
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww   
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: gF
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: Ft 
   
   REAL(KIND=8), DIMENSION(:), POINTER :: t_ij
   
   
   REAL(KIND=8), DIMENSION(SIZE(ww)-2, SIZE(ww)-2) :: t
   
   REAL(KIND=8), DIMENSION(SIZE(ww)) :: I_uPSI_f, I_vPSI_f, I_wPSI_f
   
   REAL(KIND=8) :: rho, U, V, W , P
   !------------------------------------------------------------------------------------
   
   I_uPSI_f = 0.d0
   I_vPSI_f = 0.d0
   I_wPSI_f = 0.d0
   
   CALL compute_vector_moments(ww, gF, Ft, I_uPSI_f, I_vPSI_f, I_wPSI_f)

   P   = P__ww(ww)
   rho = ww(1)
   U   = ww(2) / rho
   
   IF (sD >= 2)  V = ww(3) / rho
   IF (sD == 3)  W = ww(4) / rho
   
   ! Diagonal terms
   t(1,1) = - I_uPSI_f(2) + rho*U**2 + P

   IF (sD >= 2)  t(2,2) = - I_vPSI_f(3) + rho*V**2 + P 
   IF (sD == 3)  t(3,3) = - I_wPSI_f(4) + rho*W**2 + P


   ! Extra-diagonal terms
   IF (sD >= 2)  t(1,2) = - I_uPSI_f(3) + rho*U*V
   
   IF (sD == 3) THEN
     t(1,3) = - I_uPSI_f(4) + rho*W*U   
     t(2,3) = - I_vPSI_f(4) + rho*W*V
   ENDIF
   
      
   ALLOCATE(t_ij(sD + (sD*sD - sD)/2))

   t_ij(1) = t(1,1)
   
   IF (sD == 2) THEN
     t_ij(2) = t(2,2)
     t_ij(3) = t(1,2)
   ENDIF
   
   IF (sD == 3) THEN
     t_ij(2) = t(2,2)
     t_ij(3) = t(3,3)
     t_ij(4) = t(1,2)
     t_ij(5) = t(2,3)
     t_ij(6) = t(1,3)
   ENDIF  

      
   END FUNCTION  DPstressTensor
   
   
   
   
   
   FUNCTION   heatFlux(ww, gF, Ft)   RESULT(q)
   !----------------------------------------------------------------------
   USE transport_properties,  ONLY: Pr
   USE structures,            ONLY: REF
   
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: gF
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: Ft
   
   REAL(KIND=8), DIMENSION(:), POINTER :: q
   
   REAL(KIND=8), DIMENSION(SIZE(ww)) :: I_uPSI_f, I_vPSI_f, I_wPSI_f
   
   REAL(KIND=8) :: rho, U, V, W, P
   !----------------------------------------------------------------------

   I_uPSI_f = 0.d0
   I_vPSI_f = 0.d0
   I_wPSI_f = 0.d0
   
   CALL compute_vector_moments(ww, gF, Ft, I_uPSI_f, I_vPSI_f, I_wPSI_f)

   P = P__ww(ww)
   rho = ww(1)
   U   = ww(2) / rho
   
   IF (sD >= 2)  V = ww(3) / rho
   IF (sD == 3)  W = ww(4) / rho

   ALLOCATE (q(sD))

   ! x component --------------------------------------------------------------
   q(1) = I_uPSI_f(SIZE(ww)) - U*I_uPSI_f(2) + rho*U**3 - U*ww(SIZE(ww))
   
   IF (sD >= 2)  q(1) = q(1) - V*I_uPSI_f(3) + rho*U*V**2
   IF (sD == 3)  q(1) = q(1) - W*I_uPSI_f(4) + rho*U*W**2


   ! y component --------------------------------------------------------------
   IF (sD >= 2) THEN
   
     q(2) = I_vPSI_f(SIZE(ww)) - U*I_vPSI_f(2) - V*I_vPSI_f(3)  &
          + rho*V*U**2 + rho*V**3 - V*ww(SIZE(ww))

     IF(sD == 3)  q(2) = q(2) - W*I_vPSI_f(4) + rho*V*W**2
   
   ENDIF


    ! z component --------------------------------------------------------------
   IF (sD == 3) THEN
   
     q(3) = I_wPSI_f(SIZE(ww)) - U*I_wPSI_f(2) - V*I_wPSI_f(3) - W*I_wPSI_f(4)  &
          + rho*W*U**2 + rho*W*V**2 + rho*W**3 - W*ww(SIZE(ww))
   
   ENDIF

   q = q / (Pr / REF % Pr)

   END FUNCTION  heatFlux
   
   
   END MODULE  kinetic_equations
