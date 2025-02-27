!============================================================ 
!
!      Module: first_order_fluxes
!
! Description: First order upwind (ROe) schemes 
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

MODULE first_order_fluxes


   !============================================================ 
   USE euler_equations
   USE roe_average,         ONLY : entropy_fix_eos
   USE csr
   USE csr_pair_sys
   !============================================================ 



 !============================================================ 
 CONTAINS
 !============================================================ 



   !============================================================ 
   FUNCTION roe_fo_flux (ww, eos, ww_, eos_, j_c_d, eta) RESULT (phi)
   !============================================================ 


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(IN)  ::  eos_
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  eta

      REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,2))  ::  phi
      !------------------------------------------------------------
      ! Node-pair numerical flux
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  phi_ij
      ! Centred contribution
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_i,  ww_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(eta,1))  ::  ff_i, ff_j
      REAL(KIND=8), DIMENSION(SIZE(eta,1))  ::  etaij, eta_ver
      ! Roe scheme
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  RR, LL
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  lambda, lambda_ef 
      REAL(KIND=8)  ::  mod_eta
      ! Indices
      INTEGER  ::  c, i, j, p
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
            i  = j_c_d(1,c);    j  = j_c_d(2,c)
         
         ww_i  = ww(:,i);    ww_j  = ww(:,j)
        
         etaij = eta(:,c)
         
 
         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
            phi_ij(p) = 0.5*SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  

         ! Computation of the upwind correction
         ! ------------------------------------
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta
         
         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = ABS( eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) )

         lambda_ef = entropy_fix_eos(lambda,ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij

      ENDDO
      !============================================================


   END FUNCTION roe_fo_flux 
   !============================================================ 



   !============================================================ 
   FUNCTION implicit_roe_fo_flux ( ww, eos, ww_, eos_, j_c_d, eta, &
                                   MM)  RESULT (phi)
   !============================================================ 


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)     ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)     ::  eos
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)     ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(IN)     ::  eos_
      INTEGER,            DIMENSION(:,:), INTENT(IN)     ::  j_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)     ::  eta

      TYPE( CSR_matrix ),                 INTENT(INOUT)  ::  MM
      REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,2))    ::  phi
      !------------------------------------------------------------
      ! Node-pair numerical flux
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  phi_ij
      ! Centred contribution
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_i,  ww_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(eta,1))  ::  ff_i, ff_j
      REAL(KIND=8), DIMENSION(SIZE(eta,1))  ::  etaij, eta_ver
      ! Roe scheme
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  RR, LL
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  lambda, lambda_ef 
      REAL(KIND=8)  ::  mod_eta
      ! Implicit scheme
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  A_Roe, MM_ij
      TYPE(eos_ext_type)  ::  eos_i, eos_j ! extended eos
      ! Indices
      INTEGER  ::  c, i, j, p
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
            i  = j_c_d(1,c);    j  = j_c_d(2,c)
         
         ww_i  = ww(:,i);    ww_j  = ww(:,j)
        
         etaij = eta(:,c)
 
         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
            phi_ij(p) = 0.5*SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  

         ! Computation of the upwind correction
         ! ------------------------------------
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = ABS( eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) )

         lambda_ef = entropy_fix_eos(lambda,ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij
         
         ! Matrix for the implicit scheme
         A_Roe = 0.d0
         DO p = 1, SIZE(ww_i)
            A_Roe(p,p) = lambda_ef(p)* mod_eta
         ENDDO
         A_Roe = MATMUL( RR, MATMUL(A_Roe, LL) ) 
        
         eos_i = eos_ext__eos(eos(i))
         MM_ij =  (jacobian__ww_eos_nn(ww_i, eos_i, etaij) + A_Roe) / 2  
         CALL add_CSR_ij_sys(i, i,  MM_ij, MM)
         CALL add_CSR_ij_sys(j, i, -MM_ij, MM)
         
         eos_j = eos_ext__eos(eos(j))
         MM_ij =  (jacobian__ww_eos_nn(ww_j, eos_j, etaij) - A_Roe) / 2  
         CALL add_CSR_ij_sys(j, j, -MM_ij, MM)
         CALL add_CSR_ij_sys(i, j,  MM_ij, MM)
 
      ENDDO
      !============================================================


   END FUNCTION implicit_roe_fo_flux 
   !============================================================ 



   !============================================================ 
   SUBROUTINE init_fo_fluxes

   END SUBROUTINE init_fo_fluxes
   !============================================================ 


END MODULE first_order_fluxes


