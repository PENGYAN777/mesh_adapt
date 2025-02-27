!============================================================ 
!
!      Module: ns_num_fluxes
!
! Description: driver for the space discretizion scheme for 
!              the evaluation of Navier-Stokes numerical
!              fluxes
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!    Modified: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!============================================================ 

   MODULE ns_num_fluxes

   !----------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8)  ::  Re_ref, &  ! Reynolds reference number
   		     Pr_ref, &  ! Prandtl reference number
   		     C_visc, &  ! == 1/Re_ref
   		     C_ther	! == kk_ref/(Rgas*mu_ref*Re_ref)
   !----------------------------------------------------------

   CONTAINS


   SUBROUTINE init_ns_num_fluxes 
   !------------------------------------------------------------ 
   USE structures,   ONLY: REF
   
   IMPLICIT NONE
   !------------------------------------------------------------ 

   Re_ref = REF % Re
   Pr_ref = REF % Pr
   
   C_visc = 1 / REF % Re
   C_ther = 1 / (REF % Re * REF % Pr)

   END SUBROUTINE  init_ns_num_fluxes





   FUNCTION  ns_viscous_num_flux(vv, GG_vv, vb, GG_vb, mu, ll)   RESULT(rhs)
   !----------------------------------------------------------------------------
   USE np_quadrature,        ONLY: np_quadr_w_Gu,  &
                                   np_quadr_w_Dv,  &
                                   np_quadr_Gw_nGu
                                   
   USE node_pair_structure,  ONLY: j_c_d, j_c_b
   USE nodes,                ONLY: jd_jb
   USE metric_coefficients,  ONLY: eta, cell, chi_b, xi_bp, stiff

   !----------------------------------------------------------------------------
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: vv
   REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: GG_vv
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: vb
   REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: GG_vb
   REAL(KIND=8), DIMENSION(:),     INTENT(IN) :: mu, ll

   REAL(KIND=8), DIMENSION(SIZE(vv,1)+2, SIZE(vv,2)) :: rhs

   ! Divergence of the velocity
   REAL(KIND=8), DIMENSION(SIZE(vv,2)) :: D_vv
   
   ! Momentum equation: temporary variables
   REAL(KIND=8), DIMENSION(SIZE(GG_vv,1), SIZE(GG_vv,2), SIZE(GG_vv,3)) :: mu_GG_vv
   REAL(KIND=8), DIMENSION(SIZE(GG_vb,1), SIZE(GG_vb,2), SIZE(GG_vb,3)) :: mu_GG_vb
   REAL(KIND=8), DIMENSION(SIZE(GG_vv,1), SIZE(GG_vv,2), SIZE(GG_vv,3)) :: mu_GG_vt
   REAL(KIND=8), DIMENSION(SIZE(GG_vb,1), SIZE(GG_vb,2), SIZE(GG_vb,3)) :: mu_GG_vbt
   
   ! Energy equation: temporary variables
   REAL(KIND=8), DIMENSION(SIZE(vv,2)) :: k_e
   REAL(KIND=8), DIMENSION(SIZE(vb,2)) :: k_eb 
   
   REAL(KIND=8), DIMENSION(SIZE(vv,1), SIZE(vv,2)) :: mu_G_ke
   REAL(KIND=8), DIMENSION(SIZE(vv,1), SIZE(vv,2)) :: mu_GG_vv_vv, lvDv
   REAL(KIND=8), DIMENSION(SIZE(vb,1), SIZE(vb,2)) :: mu_GG_vb_vb, lvDvb

   ! Temporary vector and scalar quantities							      
   REAL(KIND=8), DIMENSION(SIZE(vv,1)) :: ff   

   REAL(KIND=8) :: ss

   INTEGER :: i, j, i_, j_, k, c, k_d
   !----------------------------------------------------------------------------

   k_d = SIZE(vv,1)
   rhs = 0.d0
   

   ! Initialization of D_vv and mu*GG_vv ---------------------------------------       
   ! Divergence, D_vv, of the velocity vector vv
   D_vv = (1/cell) * np_quadr_w_Dv(j_c_d, j_c_b, jd_jb, &
                                   eta, chi_b, xi_bp, vv, vb)


   ! mu_GG_vv = mu * GG_vv   and   mu_GG_vb = mu * GG_vb
   DO i = 1, SIZE(vv,2)
     mu_GG_vv(:,:,i)  =  mu(i) * GG_vv(:,:,i)        
     mu_GG_vt(:,:,i)  =  mu(i) * TRANSPOSE(GG_vv(:,:,i))
   ENDDO

   DO i = 1, SIZE(vb,2)
     mu_GG_vb (:,:,i)  =  mu(jd_jb(i)) * GG_vb(:,:,i)	   
     mu_GG_vbt(:,:,i)  =  mu(jd_jb(i)) * TRANSPOSE(GG_vb(:,:,i))   
   ENDDO
   ! ---------------------------------------------------------------------------
   

   !============================================================================
   ! MOMENTUM EQUATION
   !============================================================================
   
   ! DIV [mu GG_vv] ------------------------------------------------------------   
   ! Domain contribution
   DO k = 1, k_d
     rhs(k+1,:) = - np_quadr_Gw_nGu(j_c_d, stiff, vv(k,:), mu)
   ENDDO

   ! Boundary node-pairs contribution
   DO c = 1, SIZE(chi_b,2)

     i_ = j_c_b(1,c);	j_ = j_c_b(2,c)
     i  = jd_jb(i_);	j  = jd_jb(j_)
     
     DO k = 1, SIZE(vv,1)
       ff(k) = SUM(chi_b(:,c) * (mu_GG_vb(k,:,j_) - mu_GG_vb(k,:,i_)))
     ENDDO

     rhs(2:k_d+1,i) = rhs(2:k_d+1,i) + ff
     rhs(2:k_d+1,j) = rhs(2:k_d+1,j) - ff

   ENDDO


   ! Boundary nodes contributions
   DO i_ = 1, SIZE(xi_bp,2) 
   
     i = jd_jb(i_)
     
     DO k = 1, SIZE(vv,1)
       ff(k) = SUM(xi_bp(:,i_) * mu_GG_vb(k,:,i_))
     ENDDO

     rhs(2:k_d+1,i) = rhs(2:k_d+1,i) + ff 

   ENDDO


   ! GRAD [(mu + lambda) DIV V] ------------------------------------------------
   ! Treated using a FV-like approach. Boundary condition included in the 
   ! computation of the divergence D_vv.

   rhs(2:k_d+1,:) = rhs(2:k_d+1,:)  +  np_quadr_w_Gu(j_c_d, j_c_b, jd_jb,  &
   				                     eta, chi_b, xi_bp, (mu + ll) * D_vv)


   ! DIV [mu (GRAD V)^T] -------------------------------------------------------
   ! Identically zero if viscosity mu is constant. Treated using a FV-like approach
   ! Boundary condition is included by mu_GG_vbt(k,:,:).

   DO k = 1, SIZE(vv,1)
     rhs(k+1,:) = rhs(k+1,:)  +  np_quadr_w_Dv(j_c_d, j_c_b, jd_jb, eta, chi_b, xi_bp,  &
        				       mu_GG_vt (k,:,:), mu_GG_vbt(k,:,:))
   ENDDO


   ! - GRAD [mu DIV V] ---------------------------------------------------------
   ! Identically zero if viscosity mu is constant. Treated using a FV-like 
   ! approach. Boundary condition included in the computation of D_vv.
   
   rhs(2:k_d+1,:) = rhs(2:k_d+1,:)  -  np_quadr_w_Gu(j_c_d, j_c_b, jd_jb,   &
   				                     eta, chi_b, xi_bp, mu * D_vv)


   !============================================================================
   ! ENERGY EQUATION 
   !============================================================================
   
   ! DIV [mu GRAD |vv|^2/2] ----------------------------------------------------   
   ! Computation of the kinetic energy per unit mass 
   ! and of its gradient
   DO i = 1, SIZE(vv,2)
     k_e(i) = SUM(vv(:,i)**2)/2
   ENDDO
   
   DO i = 1, SIZE(vb,2)
     k_eb(i) = SUM(vb(:,i)**2)/2
   ENDDO
   
   mu_G_ke = np_quadr_w_Gu(j_c_d, j_c_b, jd_jb, eta, &
                           chi_b, xi_bp, k_e, k_eb)
   						   
   DO i = 1, SIZE(vv,2)
     mu_G_ke(:,i) = mu(i)*mu_G_ke(:,i)/cell(i)
   ENDDO

   ! Domain contribution
   rhs(k_d+2,:) = rhs(k_d+2,:)  -  np_quadr_Gw_nGu(j_c_d, stiff, k_e, mu)

   ! Boundary node-pairs contribution
   DO c = 1, SIZE(chi_b,2) 

     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)

     ss = SUM(chi_b(:,c) * (mu_G_ke(:,j) - mu_G_ke(:,i)))  
     
     rhs(k_d+2,i) = rhs(k_d+2,i) + ss
     rhs(k_d+2,j) = rhs(k_d+2,j) - ss

   ENDDO
   
   ! Boundary nodes contributions
   DO i_ = 1, SIZE(xi_bp,2)

     i = jd_jb(i_)
     
     ss = SUM(xi_bp(:,i_) * mu_G_ke(:,i))
     rhs(k_d+2,i) = rhs(k_d+2,i) + ss

   ENDDO



   ! DIV [(mu + lambda) (DIV V) V] ---------------------------------------------
   
   DO k = 1, SIZE(vv,1)
     lvDv(k,:) = (mu + ll) * vv(k,:) * D_vv
   ENDDO
   
   DO k = 1, SIZE(vb,1)
     lvDvb(k,:) = (mu(jd_jb) + ll(jd_jb)) * vb(k,:) * D_vv(jd_jb)
   ENDDO 

   rhs(k_d+2,:) = rhs(k_d+2,:)  +  np_quadr_w_Dv(j_c_d, j_c_b, jd_jb, &
   				                 eta, chi_b, xi_bp, lvDv, lvDvb)



   ! DIV [mu(GRAD V)^T * V  -  mu(DIV V) V]] -----------------------------------
   
   DO i = 1, SIZE(vv,2)
     DO k = 1, SIZE(vv,1)
       mu_GG_vv_vv(k,i) = SUM(vv(:,i) * mu_GG_vv(k,:,i)) - mu(i) * D_vv(i) * vv(k,i)
     ENDDO	   
   ENDDO
   
   DO i_ = 1, SIZE(vb,2)
   
     i = jd_jb(i_)
     DO k = 1, SIZE(vb,1)
       mu_GG_vb_vb(k,i_) = SUM(vb(:,i_) * mu_GG_vb(k,:,i_)) - mu(i) * D_vv(i) * vb(k,i_) 	   
     ENDDO
     	   
   ENDDO

   rhs(k_d+2,:) = rhs(k_d+2,:)  +  np_quadr_w_Dv(j_c_d, j_c_b, jd_jb, eta, chi_b, xi_bp,  &
   				                 mu_GG_vv_vv, mu_GG_vb_vb)


   ! Normalization
   rhs = C_visc * rhs  

   END FUNCTION  ns_viscous_num_flux

   
      


   FUNCTION  ns_thermal_num_flux(T, G_T, Tb, G_Tb, kk)   RESULT(phi_T)
   !----------------------------------------------------------------------------
   USE np_quadrature,        ONLY: np_quadr_Gw_nGu
   USE node_pair_structure,  ONLY: j_c_d, j_c_b
   USE nodes,                ONLY: jd_jb
   USE metric_coefficients,  ONLY: chi_b, xi_bp, stiff
   
   !----------------------------------------------------------------------------
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: T
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: Tb
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: G_T
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: G_Tb
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: kk
   
   REAL(KIND=8), DIMENSION(SIZE(T)) :: phi_T
   
   REAL(KIND=8) :: ss

   INTEGER :: c, i, j, i_, j_

   REAL(KIND=8), DIMENSION(SIZE(G_T,1),SIZE(G_T,2)) :: void_2
   REAL(KIND=8), DIMENSION(SIZE(Tb))                :: void_1
   !----------------------------------------------------------------------------

   ! Dummy Assignation
   void_1 = Tb
   void_2 = G_T


   ! DIV [kk GG_T] -------------------------------------------------------------   
   ! Domain contribution
   phi_T = - np_quadr_Gw_nGu(j_c_d, stiff, T, kk)

   ! Boundary node-pairs contribution
   DO c = 1, SIZE(chi_b,2)
   
     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)
     
     ss = SUM(chi_b(:,c) * (kk(j)*G_Tb(:,j_) - kk(i)*G_Tb(:,i_)))
     
     phi_T(i) = phi_T(i) + ss
     phi_T(j) = phi_T(j) - ss

   ENDDO
   
   ! Boundary nodes contributions
   DO i = 1, SIZE(xi_bp,2)
     
     ss = SUM(xi_bp(:,i) * kk(jd_jb(i)) * G_Tb(:,i))
     phi_T(jd_jb(i)) = phi_T(jd_jb(i)) + ss

   ENDDO


   ! Normalization
   phi_T = C_ther * phi_T      
   
   END FUNCTION  ns_thermal_num_flux

   
       
   

   FUNCTION  ns_compute_dt(ww, eos, CFL, variable_dt)   RESULT(dt)

   !  Compute the time step.  Local time stepping (steady flow
   !  only) or constant time step are allowed.
   !-----------------------------------------------------------------------
   USE structures,            ONLY: eos_type
   USE node_pair_structure,   ONLY: j_c_d
   USE metric_coefficients,   ONLY: eta, cell
   USE transport_properties,  ONLY: transport_coefficients
   USE thermodynamics,        ONLY: cv__T_v
   USE mp_interface
   
   !-----------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: ww
   TYPE(eos_type), DIMENSION(:),   INTENT(IN) :: eos
   REAL(KIND=8),		   INTENT(IN) :: CFL 
   LOGICAL,			   INTENT(IN) :: variable_dt  

   REAL(KIND=8), DIMENSION(SIZE(ww,2))  :: dt 
   REAL(KIND=8), DIMENSION(MP_wProc)    :: MP_minDt
   REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ver

   REAL(KIND=8) :: mod_eta, T, v, cv, nu, mu, kk, minDt
   
   INTEGER :: c, i, j
   !-----------------------------------------------------------------------

   ! Variable time step. (ATTN: for steady solution only)
   !
   !		  h_i^2
   !  dt_i = CFL -------
   !		  nu_i
   !
   !	     2 |C_i|
   !  h_i = ----------,   for each i in K_i
   !	     INT_dC_i  

   dt = 0
   DO c = 1, SIZE(j_c_d,2) 

     i = j_c_d(1,c);  j = j_c_d(2,c)

     mod_eta = SQRT(SUM(eta(:,c)*eta(:,c)))
     eta_ver = eta(:,c)/mod_eta

     dt(i) = dt(i) + mod_eta
     dt(j) = dt(j) + mod_eta

   ENDDO

   ! Local time step (The sum of 1.d-8
   ! is an artifact necessary only in multi-processor
   ! computation. It affects only the time step
   ! relative to the extended nodes whose
   ! time advance is not relevant)
   dt = (2 * cell) / (dt+1.d-8)

   DO i = 1, SIZE(ww,2)
   
      T = eos(i)%T
      v = 1/ww(1,i)
      
      CALL transport_coefficients(T, v,  mu, kk) 
      cv = cv__T_v(T,v)
      
      nu = 2 * MAX(mu, kk/cv) * v

      dt(i) = CFL * (dt(i)**2) / (nu/Re_ref)
   	       
   ENDDO


   ! Global time step  --->  dt = MIN dt(i),   for all i in K_i
   IF (.NOT. MP_job .AND. .NOT.variable_dt)  dt = MINVAL(dt) 
   
   IF (MP_job  .AND. .NOT.variable_dt) THEN
     
     minDt = MINVAL(dt(1:NjP))
     
     CALL MP_gather_all_minDt(minDt, MP_minDt)
     
     dt = MINVAL(MP_minDt)
     
   ENDIF  


   END FUNCTION  ns_compute_dt
   
   
   END  MODULE ns_num_fluxes
     
       

    

    
