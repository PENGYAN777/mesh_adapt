!============================================================ 
!
!      Module: high_resolution_fluxes
!
! Description: Implicit/explicit agorithms for 
!              high-resolution schemes in node-pair 
!              formulation. Galerkin/Roe schemes.
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

MODULE high_resolution_fluxes


   !============================================================ 
   USE euler_equations
   USE csr
   USE csr_pair_sys

   USE roe_average,   ONLY: entropy_fix_eos
   !============================================================ 

    
   !============================================================ 
   INTEGER, PARAMETER  ::  HR_NO_LIMITER           = -1, &
                           HR_FULL_LIMITER         =  0, &
                           HR_VAN_LEER_LIMITER     =  1, &
                           HR_MINMOD_LIMITER       =  2, &
                           HR_SUPERBEE_LIMITER     =  3, &
                           HR_MONOCENTR_LIMITER    =  4, &
                           HR_VAN_ALBADA_LIMITER   =  5, &                        
                           HR_YEE_MINMOD_2         =  6, &                        
                           HR_YEE_MINMOD_3         =  7                        
                           
   CHARACTER(*), DIMENSION(-1:7), PARAMETER  ::  hr_lim_names = &
                        (/ 'NO LIMITER: FULL II ORDER  ', &
                           'FULL LIMITER: I ORDER      ', &
                           'VAN LEER LIMITER           ', &
                           'MINMOD LIMITER             ', &
                           'SUPERBEE LIMITER           ', &
                           'MONOTONIZED CENTRAL LIMITER', &
                           'VAN ALBADA LIMITER         ', & 
                           'YEE MINMOD 2               ', &
                           'YEE MINMOD 3               ' /)

   INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: limiter_type
   INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: hr_limiter_type
   !============================================================ 

   CONTAINS

   !============================================================ 
   SUBROUTINE read_param_hr_fluxes(idf)
   !------------------------------------------------------------
   IMPLICIT NONE

   INTEGER :: idf
   INTEGER :: ww_size
   !------------------------------------------------------------

   READ (idf,*) ww_size

   ALLOCATE (   limiter_type(ww_size))
   ALLOCATE (hr_limiter_type(ww_size))

   READ (idf,*) hr_limiter_type

   limiter_type = hr_limiter_type

   END SUBROUTINE read_param_hr_fluxes
   !============================================================ 




   !============================================================ 
   SUBROUTINE init_hr_fluxes

   END SUBROUTINE init_hr_fluxes
   !============================================================




   !============================================================ 
   SUBROUTINE write_param_hr_fluxes(idf)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE

      INTEGER  ::  idf
      !------------------------------------------------------------
      INTEGER  ::  p
      !------------------------------------------------------------


      WRITE(idf,*) '   PARAMETERS FOR HIGH-RESOLUTION NUMERICAL FLUXES'
      
      DO p = 1, SIZE(hr_limiter_type)
         WRITE(idf,*) '   Limiter ', p,': ', hr_lim_names(hr_limiter_type(p))
      ENDDO
      
      WRITE(idf,*)
      

   END SUBROUTINE write_param_hr_fluxes
   !============================================================ 



   !============================================================ 
   FUNCTION Roe_Galerkin_HR_flux (ww, eos, ww_, eos_, &
                                  j_c_d, cs_c_d, eta, Drr) RESULT (phi)
   !============================================================ 


      ! HIGH RESOLUTION FLUX-LIMITED METHOD:
      ! ROE 1st ORDER UPWIND + (PURE)GALERKIN 2nd ORDER
      ! Limiters in Rebay's form, upwinding selection
      !============================================================
      use mp_interface
            
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(IN)  ::  eos_
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d, cs_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  eta
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  Drr

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
      ! High resolution
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  LL_isi, LL_jjs      
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_i, alpha_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_cen, alpha_upw, &
                                                alpha_lim
      REAL(KIND=8)  ::  Drr_ij, Drr_iis, Drr_jjs
      REAL(KIND=8)  ::  HR_corr
      ! Indices
      INTEGER  ::  c, i, j, is, js, p, c_isi, c_jjs
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
          !------------------------------------------------------------ 
          ! Topology, solution and metric quantities 
          ! of the c-th node-pair
          ! ----------------------------------------
          i  = j_c_d(1,c);    j  = j_c_d(2,c)
          is = j_c_d(3,c);    js = j_c_d(4,c)
         
          ww_i  = ww(:,i);    ww_j  = ww(:,j)
          ww_is = ww(:,is);   ww_js = ww(:,js)
        
        
          etaij = eta(:,c)
         
          Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)

 
         !============================================================
         ! Centred contribution to phi_ij
         !==============================
	 
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
             phi_ij(p) = 0.5d0 * SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         !=====================  
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 

         lambda_ef = entropy_fix_eos(ABS(lambda),ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta 
         ENDDO
         !============================================================


         !============================================================
         ! Higher order correction
         ! =======================  
    
         ! Computation of ``extended'' variation of is--i and j--js
         ! --------------------------------------------------------
         c_isi = ABS(cs_c_d(1,c))	 
         LL_isi =  left_eigenvectors__ww_eos_nn(ww_(:,c_isi), eos_(c_isi), eta_ver)
         alpha_i  = MATMUL( LL_isi, (ww_i - ww_is) * (drr_ij/drr_iis) )

         c_jjs = ABS(cs_c_d(2,c))
         LL_jjs =  left_eigenvectors__ww_eos_nn(ww_(:,c_jjs), eos_(c_jjs), eta_ver)
         alpha_j  = MATMUL( LL_jjs, (ww_js - ww_j)* (drr_ij/drr_jjs) )

         ! Computation of the HR correction for each characteristic field
         ! --------------------------------------------------------------
         DO p = 1, SIZE(ww_i)
         
            alpha_cen(p) = alpha(p)                           ! centred variation
            alpha_upw(p) = (alpha_i(p) + alpha_j(p)) / 2  &   ! upwind  variation
                      + (alpha_i(p) - alpha_j(p)) * SIGN(0.5d0, lambda(p))

            alpha_lim(p) = rebay_form_limiter( alpha_cen(p),   &
                                               alpha_upw(p),   &
                                               limiter_type(p) )

            HR_corr = 0.5 * lambda_ef(p) * mod_eta * alpha_lim(p)

            phi_ij  =  phi_ij  +  RR(:,p) * HR_corr 
           !^^^^^^

         ENDDO
         !============================================================

         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij

      ENDDO
      !============================================================



   END FUNCTION Roe_Galerkin_HR_flux 
   !============================================================ 



   !============================================================ 
   FUNCTION Roe_Galerkin_HR_flux_cen (ww, eos, ww_, eos_, &
                                      j_c_d, cs_c_d, eta, Drr) RESULT (phi)
   !============================================================ 


   ! HIGH RESOLUTION FLUX-LIMITED METHOD:
   ! ROE 1st ORDER UPWIND + (PURE)GALERKIN 2nd ORDER
   ! Limiters in Rebay's form, centred selection (YEE)


      !============================================================
      USE roe_average
      !============================================================


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(IN)  ::  eos_
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d, cs_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  eta
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  Drr

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
      ! High resolution
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  LL_isi, LL_jjs      
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_i, alpha_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_lim
      REAL(KIND=8)  ::  Drr_ij, Drr_iis, Drr_jjs
      REAL(KIND=8)  ::  HR_corr
      ! Indices
      INTEGER  ::  c, i, j, is, js, p, c_isi, c_jjs
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
          !------------------------------------------------------------ 
          ! Topology, solution and metric quantities 
          ! of the c-th node-pair
          ! ----------------------------------------
             i  = j_c_d(1,c);    j  = j_c_d(2,c)
             is = j_c_d(3,c);    js = j_c_d(4,c)
         
          ww_i  = ww(:,i);    ww_j  = ww(:,j)
          ww_is = ww(:,is);   ww_js = ww(:,js)
        
        
          etaij = eta(:,c)
         
          Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)

 
         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
             phi_ij(p) = 0.5d0 * SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 

         lambda_ef = entropy_fix_eos(ABS(lambda),ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         !============================================================
         ! Higher order correction
         ! =======================  

         ! Computation of ``extended'' variation of is--i and j--js
         ! --------------------------------------------------------
         c_isi = ABS(cs_c_d(1,c))
         LL_isi =  left_eigenvectors__ww_eos_nn(ww_(:,c_isi), eos_(c_isi), eta_ver)
         alpha_i  = MATMUL( LL_isi, (ww_i - ww_is) * (drr_ij/drr_iis) )

         c_jjs = ABS(cs_c_d(2,c))
         LL_jjs =  left_eigenvectors__ww_eos_nn(ww_(:,c_jjs), eos_(c_jjs), eta_ver)
         alpha_j  = MATMUL( LL_jjs, (ww_js - ww_j)* (drr_ij/drr_jjs) )

         ! Computation of the HR correction for each characteristic field
         ! --------------------------------------------------------------
         DO p = 1, SIZE(ww_i)
            alpha_lim(p) = rebay_form_centred_limiter( alpha(p),       &
                                                       alpha_i(p),     &
                                                       alpha_j(p),     &
                                                       limiter_type(p) )

            HR_corr = 0.5 * lambda_ef(p) * mod_eta * alpha_lim(p)

            phi_ij  =  phi_ij  +  RR(:,p) * HR_corr 
           !^^^^^^  

         ENDDO
         !============================================================

         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij

      ENDDO
      !============================================================


   END FUNCTION Roe_Galerkin_HR_flux_cen 
   !============================================================ 



   !============================================================ 
   FUNCTION Roe_LaxWendroff_HR_flux (ww, eos, ww_, eos_, dt,   &
                  j_c_d, cs_c_d, eta, stiff_T, Drr) RESULT (phi)
   !============================================================ 


   ! HIGH RESOLUTION FLUX-LIMITED METHOD:
   ! ROE 1st ORDER UPWIND + LAX WENDROFF 2nd ORDER


      !============================================================
      USE roe_average
      !============================================================


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),     INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),     INTENT(IN)  ::  eos_
      REAL(KIND=8),       DIMENSION(:),     INTENT(IN)  ::  dt
      INTEGER,            DIMENSION(:,:),   INTENT(IN)  ::  j_c_d, cs_c_d
      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  eta
      REAL(KIND=8),       DIMENSION(:,:,:), INTENT(IN)  ::  stiff_T
      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  Drr

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
      ! High resolution
      REAL(KIND=8), DIMENSION(SIZE(ww,1), &
                              SIZE(ww,1), &
                              SIZE(ww,1)-2) ::  JA
      REAL(KIND=8), DIMENSION(SIZE(ww,1), &
                              SIZE(ww,1))   ::  KA, AL
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  LL_isi, LL_jjs      
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_i, alpha_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_cen, alpha_upw, &
                                                alpha_lim
      REAL(KIND=8)  ::  Drr_ij, Drr_iis, Drr_jjs
      ! Indices
      INTEGER  ::  c, i, j, is, js, p, c_isi, c_jjs, k1, k2, space_dim
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
      
          !------------------------------------------------------------ 
          ! Topology, solution and metric quantities 
          ! of the c-th node-pair
          ! ----------------------------------------
            i  = j_c_d(1,c);    j  = j_c_d(2,c)
            is = j_c_d(3,c);    js = j_c_d(4,c)
         
         ww_i  = ww(:,i);    ww_j  = ww(:,j)
         ww_is = ww(:,is);   ww_js = ww(:,js)
        
         etaij = eta(:,c)
         
         Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)

 
         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
             phi_ij(p) = 0.5d0 * SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 

         lambda_ef = entropy_fix_eos(ABS(lambda),ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         !============================================================
         ! Higher order correction
         ! =======================  
    
         ! Computation of ``extended'' variation of is--i and j--js
         ! --------------------------------------------------------
         c_isi = ABS(cs_c_d(1,c))
         LL_isi =  left_eigenvectors__ww_eos_nn(ww_(:,c_isi), eos_(c_isi), eta_ver)
         alpha_i  = MATMUL( LL_isi, (ww_i - ww_is) * (drr_ij/drr_iis) )

         c_jjs = ABS(cs_c_d(2,c))
         LL_jjs =  left_eigenvectors__ww_eos_nn(ww_(:,c_jjs), eos_(c_jjs), eta_ver)
         alpha_j  = MATMUL( LL_jjs, (ww_js - ww_j) * (drr_ij/drr_jjs) )

         ! Computation of the HR correction for each characteristic field
         ! --------------------------------------------------------------
         DO p = 1, SIZE(ww_i)

            ! Centred jumps
            alpha_cen(p) = alpha(p)
            ! Upwind jumps
            alpha_upw(p) = (alpha_i(p) + alpha_j(p)) / 2  & 
                         + (alpha_i(p) - alpha_j(p)) * SIGN(0.5d0, lambda(p))
            ! Limited jumps
            alpha_lim(p) = rebay_form_limiter( alpha_cen(p),   &
                                               alpha_upw(p),   & 
                                               limiter_type(p) )

!**********
!            IF (ABS(alpha_cen(p)) > 0.e-16) THEN
!               alpha_limited(p) = &
!               standard_form_limiter(alpha_upw(p)/alpha_cen(p), &
!                                                   limiter_type(p))
!
!            ELSE 
!               alpha_limited(p) =0.d0
!            ENDIF
! ALTERNATIVE using standard form limiter
!**********

         ENDDO

         JA = jacobian__ww_eos(ww_(:,c), eos_(c))

         space_dim = SIZE(ww_i) - 2 
         KA = 0.d0
         DO k1 = 1, space_dim
            DO k2 = 1, space_dim
               KA = KA + MATMUL(JA(:,:,k1),JA(:,:,k2)) * stiff_T(k1,k2,c)
            ENDDO
         ENDDO

         AL = 0.d0
         DO p = 1, SIZE(ww_i)
             AL(p,p) = lambda_ef(p)*mod_eta
         ENDDO

         KA = MATMUL(RR,AL) + dt(i) * MATMUL( KA, RR)
         phi_ij = phi_ij  +  MATMUL(KA,alpha_lim) / 2
        !^^^^^^
!         phi_ij = phi_ij  +  MATMUL(KA,alpha_limited*alpha) / 2
! ALTERNATIVE using standard form limiter
         !============================================================


         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij

      ENDDO
      !============================================================


   END FUNCTION Roe_LaxWendroff_HR_flux 
   !============================================================ 




   !============================================================ 
   FUNCTION Roe_LaxWendroff_HR_flux_rot (ww, eos, ww_, eos_, dt, &
                     j_c_d, cs_c_d, eta, stiff, Drr) RESULT (phi)
   !============================================================ 


   ! HIGH RESOLUTION FLUX-LIMITED METHOD:
   ! ROE 1st ORDER UPWIND + LAX WENDROFF 2nd ORDER


      !============================================================
      USE roe_average
      !============================================================


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(IN)  ::  eos_
      REAL(KIND=8),       DIMENSION(:),   INTENT(IN)  ::  dt
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d, cs_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  eta
      REAL(KIND=8),       DIMENSION(:),   INTENT(IN)  ::  stiff
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  Drr

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
      ! High resolution
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  LL_isi, LL_jjs      
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_i, alpha_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_cen, alpha_upw, &
                                                alpha_lim
      REAL(KIND=8)  ::  Drr_ij, Drr_iis, Drr_jjs
      REAL(KIND=8)  ::  HR_corr
      ! Indices
      INTEGER  ::  c, i, j, is, js, p, c_isi, c_jjs
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
          !------------------------------------------------------------ 
          ! Topology, solution and metric quantities 
          ! of the c-th node-pair
          ! ----------------------------------------
      
            i  = j_c_d(1,c);    j  = j_c_d(2,c)
            is = j_c_d(3,c);    js = j_c_d(4,c)
         
         ww_i  = ww(:,i);    ww_j  = ww(:,j)
         ww_is = ww(:,is);   ww_js = ww(:,js)
        
         etaij = eta(:,c)
         
         Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)


         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
             phi_ij(p) = 0.5d0 * SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 

         lambda_ef = entropy_fix_eos(ABS(lambda),ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         !============================================================
         ! Higher order correction
         ! =======================  
    
         ! Computation of ``extended'' variation of is--i and j--js
         ! --------------------------------------------------------
         c_isi = ABS(cs_c_d(1,c))
         LL_isi =  left_eigenvectors__ww_eos_nn(ww_(:,c_isi), eos_(c_isi), eta_ver)
         alpha_i  = MATMUL( LL_isi, (ww_i - ww_is) * (drr_ij/drr_iis) )

         c_jjs = ABS(cs_c_d(2,c))
         LL_jjs =  left_eigenvectors__ww_eos_nn(ww_(:,c_jjs), eos_(c_jjs), eta_ver)
         alpha_j  = MATMUL( LL_jjs, (ww_js - ww_j)* (drr_ij/drr_jjs) )

         ! Computation of the HR correction for each characteristic field
         ! --------------------------------------------------------------
         DO p = 1, SIZE(ww_i)

            alpha_cen(p) = alpha(p)                           ! centred variation
            alpha_upw(p) = (alpha_i(p) + alpha_j(p)) / 2  &   ! upwind  variation
              + (alpha_i(p) - alpha_j(p)) * SIGN(0.5d0, lambda(p))

            alpha_lim(p) = rebay_form_limiter( alpha_cen(p),   &
                                               alpha_upw(p),   & 
                                               limiter_type(p) )

            HR_corr =  lambda_ef(p)* mod_eta &
!                    +  dt(i)*(lambda_ef(p)**2) * stiff(c)  
                    +  dt(i)*(lambda(p)**2) * stiff(c)  
                    
            HR_corr = (HR_corr / 2 ) * alpha_lim(p)

            phi_ij  =  phi_ij  +  RR(:,p) * HR_corr   

         ENDDO
         !============================================================


         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij

      ENDDO
      !============================================================


   END FUNCTION Roe_LaxWendroff_HR_flux_rot 
   !============================================================ 



!============================================================ 
!************************************************************
!
!   IMPLICIT SCHEMES
!
!************************************************************
!============================================================ 



   !============================================================ 
   FUNCTION imp_Roe_Galerkin_HR_flux (ww, eos, ww_, eos_, &
                      j_c_d, cs_c_d, eta, Drr,  MM) RESULT (phi)
   !============================================================ 


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(IN)  ::  eos_
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d, cs_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  eta
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  Drr

      TYPE( CSR_matrix ),           INTENT(INOUT)      ::  MM
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
      ! High resolution
      REAL(KIND=8), DIMENSION(SIZE(ww,1), &
                              SIZE(ww,1))   ::  A_Roe, MM_ij
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  LL_isi, LL_jjs
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_i, alpha_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_cen, alpha_upw, &
                                                alpha_lim
      TYPE(eos_ext_type)  ::  eos_i, eos_j ! extended eos
      REAL(KIND=8)  ::  Drr_ij, Drr_iis, Drr_jjs
      REAL(KIND=8)  ::  HR_corr
      ! Indices
      INTEGER  ::  c, i, j, is, js, p, c_isi, c_jjs
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
      
          !------------------------------------------------------------ 
          ! Topology, solution and metric quantities 
          ! of the c-th node-pair
          ! ----------------------------------------
            i  = j_c_d(1,c);    j  = j_c_d(2,c)
            is = j_c_d(3,c);    js = j_c_d(4,c)
         
         ww_i  = ww(:,i);    ww_j  = ww(:,j)
         ww_is = ww(:,is);   ww_js = ww(:,js)
        
         etaij = eta(:,c)
         
         Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)

 
         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
             phi_ij(p) = 0.5d0 * SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 

         lambda_ef = entropy_fix_eos(ABS(lambda),ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         !============================================================
         ! Higher order correction
         ! =======================  
    
         ! Computation of ``extended'' variation of is--i and j--js
         ! --------------------------------------------------------
         c_isi = ABS(cs_c_d(1,c))
         LL_isi =  left_eigenvectors__ww_eos_nn(ww_(:,c_isi), eos_(c_isi), eta_ver)
         alpha_i  = MATMUL( LL_isi, (ww_i - ww_is) * (drr_ij/drr_iis) )

         c_jjs = ABS(cs_c_d(2,c))
         LL_jjs =  left_eigenvectors__ww_eos_nn(ww_(:,c_jjs), eos_(c_jjs), eta_ver)
         alpha_j  = MATMUL( LL_jjs, (ww_js - ww_j)* (drr_ij/drr_jjs) )

          
         ! Computation of the HR correction for each characteristic field
         ! --------------------------------------------------------------
         DO p = 1, SIZE(ww_i)
         
            alpha_cen(p) = alpha(p)                           ! centred variation
            alpha_upw(p) = (alpha_i(p) + alpha_j(p)) / 2  &   ! upwind  variation
                         + (alpha_i(p) - alpha_j(p)) * SIGN(0.5d0, lambda(p))

            alpha_lim(p) = rebay_form_limiter( alpha_cen(p),   &
                                               alpha_upw(p),   &
                                               limiter_type(p) )

            HR_corr = 0.5 * lambda_ef(p) * mod_eta * alpha_lim(p)

            phi_ij  =  phi_ij  +  RR(:,p) * HR_corr   
           !^^^^^^

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


   END FUNCTION imp_Roe_Galerkin_HR_flux 
   !============================================================ 



   !============================================================ 
   FUNCTION imp_Roe_LaxWendroff_HR_flux (ww, eos, ww_, eos_, dt, &
                 j_c_d, cs_c_d, eta, stiff_T,  Drr,  MM) RESULT (phi)
   !============================================================ 


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),     INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),     INTENT(IN)  ::  eos_
      REAL(KIND=8),       DIMENSION(:),     INTENT(IN)  ::  dt
      INTEGER,            DIMENSION(:,:),   INTENT(IN)  ::  j_c_d, cs_c_d
      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  eta
      REAL(KIND=8),       DIMENSION(:,:,:), INTENT(IN)  ::  stiff_T
      REAL(KIND=8),       DIMENSION(:,:),   INTENT(IN)  ::  Drr
 
      TYPE( CSR_matrix ),           INTENT(INOUT)      ::  MM
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
      ! High resolution
      REAL(KIND=8), DIMENSION(SIZE(ww,1), &
                              SIZE(ww,1), &
                              SIZE(ww,1)-2) ::  JA
      REAL(KIND=8), DIMENSION(SIZE(ww,1), &
                              SIZE(ww,1))   ::  KA, AL, A_Roe, MM_ij
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  LL_isi, LL_jjs
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_i, alpha_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_cen, alpha_upw, &
                                                alpha_lim
      TYPE(eos_ext_type)  ::  eos_i, eos_j ! extended eos
      REAL(KIND=8)  ::  Drr_ij, Drr_iis, Drr_jjs
      ! Indices
      INTEGER  ::  c, i, j, is, js, p, c_isi, c_jjs, k1, k2, space_dim
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
      
          !------------------------------------------------------------ 
          ! Topology, solution and metric quantities 
          ! of the c-th node-pair
          ! ----------------------------------------
            i  = j_c_d(1,c);    j  = j_c_d(2,c)
            is = j_c_d(3,c);    js = j_c_d(4,c)
         
         ww_i  = ww(:,i);    ww_j  = ww(:,j)
         ww_is = ww(:,is);   ww_js = ww(:,js)
        
         etaij = eta(:,c)
         
         Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)


         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
             phi_ij(p) = 0.5d0 * SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 

         lambda_ef = entropy_fix_eos(ABS(lambda),ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         !============================================================
         ! Higher order correction
         ! =======================  
    
         ! Computation of ``extended'' variation of is--i and j--js
         ! --------------------------------------------------------
         c_isi = ABS(cs_c_d(1,c))
         LL_isi =  left_eigenvectors__ww_eos_nn(ww_(:,c_isi), eos_(c_isi), eta_ver)
         alpha_i  = MATMUL( LL_isi, (ww_i - ww_is) * (drr_ij/drr_iis) )

         c_jjs = ABS(cs_c_d(2,c))
         LL_jjs =  left_eigenvectors__ww_eos_nn(ww_(:,c_jjs), eos_(c_jjs), eta_ver)
         alpha_j  = MATMUL( LL_jjs, (ww_js - ww_j)* (drr_ij/drr_jjs) )


         ! Computation of the HR correction for each characteristic field
         ! --------------------------------------------------------------
         DO p = 1, SIZE(ww_i)

            alpha_cen(p) = alpha(p)                           ! centred variation
            alpha_upw(p) = (alpha_i(p) + alpha_j(p)) / 2  &   ! upwind  variation
                      + (alpha_i(p) - alpha_j(p)) * SIGN(0.5d0, lambda(p))

            alpha_lim(p) = rebay_form_limiter( alpha_cen(p),   &
                                               alpha_upw(p),   &
                                               limiter_type(p) )

!            IF (ABS(alpha_cen(p)) > 0.e-16) THEN
!               alpha_limited(p) = &
!               standard_form_limiter(alpha_upw(p)/alpha_cen(p),  &
!                                                 limiter_type(p) )
!
!            ELSE
!               alpha_limited(p) =0.d0
!            ENDIF
            
         ENDDO

!         JA = jacobian(ww_)
!
!         space_dim = SIZE(ww,1) - 2 
!         KA = 0.d0
!         DO k1 = 1, space_dim
!            DO k2 = 1, space_dim
!               KA = KA + MATMUL(JA(:,:,k1),JA(:,:,k2)) * stiff_T(k1,k2,c)
!            ENDDO
!         ENDDO
!
!         AL = 0.d0
!         DO p = 1, SIZE(ww_i)
!             AL(p,p) = lambda_ef(p)*mod_eta
!         ENDDO
!                            ! Warning: variable time step not allowed
!         KA = MATMUL(RR,AL) + dt(i) * MATMUL( KA, RR)
!
!         AL = 0.d0
!         DO p = 1, SIZE(ww_i)
!             AL(p,p) = alpha_limited(p)
!         ENDDO
!         phi_ij = phi_ij  &
!                +  (MATMUL(MATMUL(KA, MATMUL(AL,LL)),(ww_j-ww_i))) / 2

         JA = jacobian__ww_eos(ww_(:,c), eos_(c))

         space_dim = SIZE(ww_i) - 2 
         KA = 0.d0
         DO k1 = 1, space_dim
            DO k2 = 1, space_dim
               KA = KA + MATMUL(JA(:,:,k1),JA(:,:,k2)) * stiff_T(k1,k2,c)
            ENDDO
         ENDDO

         AL = 0.d0
         DO p = 1, SIZE(ww_i)
             AL(p,p) = lambda_ef(p)*mod_eta
         ENDDO

         KA = MATMUL(RR,AL) + dt(i) * MATMUL( KA, RR)
         phi_ij = phi_ij  +  MATMUL(KA,alpha_lim) / 2
         !============================================================


         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij

         ! Matrix for the implicit scheme
         A_Roe = 0.d0
         DO p = 1, SIZE(ww_i)
            A_Roe(p,p) = lambda_ef(p)* mod_eta
         ENDDO
         A_Roe = MATMUL( RR, MATMUL(A_Roe, LL) ) 
!         A_Roe = MATMUL( RR, MATMUL(A_Roe, LL) ) &
!               + MATMUL( KA, MATMUL(AL,    LL) )
! WARNING PLUS OR MINUS???
        
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


   END FUNCTION imp_Roe_LaxWendroff_HR_flux 
   !============================================================ 

   !============================================================ 
   FUNCTION imp_Roe_Galerkin_HR_flux_cen (ww, eos, ww_, eos_, &
                      j_c_d, cs_c_d, eta, Drr,  MM) RESULT (phi)
   !============================================================ 


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(IN)  ::  eos_
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d, cs_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  eta
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  Drr

      TYPE( CSR_matrix ),           INTENT(INOUT)      ::  MM
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
      ! High resolution
      REAL(KIND=8), DIMENSION(SIZE(ww,1), &
                              SIZE(ww,1))   ::  A_Roe, MM_ij
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  LL_isi, LL_jjs
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_i, alpha_j
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha_lim !, alpha_cen, alpha_upw

      TYPE(eos_ext_type)  ::  eos_i, eos_j ! extended eos
      REAL(KIND=8)  ::  Drr_ij, Drr_iis, Drr_jjs
      REAL(KIND=8)  ::  HR_corr
      ! Indices
      INTEGER  ::  c, i, j, is, js, p, c_isi, c_jjs
      !============================================================


      phi = 0.d0

      !============================================================
      ! Loop on the node-pairs
      ! ====================== 
      DO c = 1, SIZE(j_c_d,2)
      
      
          !------------------------------------------------------------ 
          ! Topology, solution and metric quantities 
          ! of the c-th node-pair
          ! ----------------------------------------
            i  = j_c_d(1,c);    j  = j_c_d(2,c)
            is = j_c_d(3,c);    js = j_c_d(4,c)
         
         ww_i  = ww(:,i);    ww_j  = ww(:,j)
         ww_is = ww(:,is);   ww_js = ww(:,js)
        
         etaij = eta(:,c)
         
         Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)

 
         !============================================================
         ! Centred contribution to phi_ij
         ! ==============================
         ff_i = flux__ww_eos(ww_i, eos(i))  
         ff_j = flux__ww_eos(ww_j, eos(j))

         DO p = 1, SIZE(ww_i)
             phi_ij(p) = 0.5d0 * SUM( (ff_i(p,:) + ff_j(p,:)) * etaij )
         ENDDO
         !============================================================


         !============================================================
         ! Roe upwind correction
         ! =====================  
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver) 

         lambda_ef = entropy_fix_eos(ABS(lambda),ww_(:,c), eos_(c), eta_ver)  
                    !^^^^^^^^^^^ 

         DO p = 1, SIZE(ww_i)
            phi_ij  =  phi_ij &
                    -  0.5 * RR(:,p) * lambda_ef(p) * alpha(p) * mod_eta
         ENDDO
         !============================================================


         !============================================================
         ! Higher order correction
         ! =======================  
    
         ! Computation of ``extended'' variation of is--i and j--js
         ! --------------------------------------------------------
         c_isi = ABS(cs_c_d(1,c))
         LL_isi =  left_eigenvectors__ww_eos_nn(ww_(:,c_isi), eos_(c_isi), eta_ver)
         alpha_i  = MATMUL( LL_isi, (ww_i - ww_is) * (drr_ij/drr_iis) )

         c_jjs = ABS(cs_c_d(2,c))
         LL_jjs =  left_eigenvectors__ww_eos_nn(ww_(:,c_jjs), eos_(c_jjs), eta_ver)
         alpha_j  = MATMUL( LL_jjs, (ww_js - ww_j)* (drr_ij/drr_jjs) )

          
         ! Computation of the HR correction for each characteristic field
         ! --------------------------------------------------------------
         DO p = 1, SIZE(ww_i)
            alpha_lim(p) = rebay_form_centred_limiter( alpha(p),       &
                                                       alpha_i(p),     &
                                                       alpha_j(p),     &
                                                       limiter_type(p) )

            HR_corr = 0.5 * lambda_ef(p) * mod_eta * alpha_lim(p)

            phi_ij  =  phi_ij  +  RR(:,p) * HR_corr 
           !^^^^^^  

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


   END FUNCTION imp_Roe_Galerkin_HR_flux_cen 
   !============================================================ 
 

!============================================================ 
!************************************************************
!
!   LIMITERS
!
!************************************************************
!============================================================ 



   !============================================================ 
   FUNCTION rebay_form_limiter(a_cen, a_upw, limiter_type)     &
                                                  RESULT (a_lim)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN)  ::  a_cen, a_upw
      INTEGER,      INTENT(IN)  ::  limiter_type

      REAL(KIND=8)              ::  a_lim
      !------------------------------------------------------------
      REAL(KIND=8), PARAMETER   ::  zero = 0,  half = 0.5,  one = 1, &
                                    eps = 1.e-16
      !------------------------------------------------------------
 

      SELECT CASE(limiter_type)

         !------------------------------------------------------------
         CASE(HR_NO_LIMITER);  ! second order
            a_lim = a_cen
         !------------------------------------------------------------
         CASE(HR_FULL_LIMITER);  ! no limiter, first-order upwind
            a_lim = 0
         !------------------------------------------------------------
         CASE(HR_VAN_LEER_LIMITER);  ! van Leer
            a_lim = (a_cen*ABS(a_upw) + ABS(a_cen)*a_upw) &
                  / (ABS(a_cen) + ABS(a_upw) + eps)
         !------------------------------------------------------------
         CASE(HR_MINMOD_LIMITER);  ! minmod
            a_lim = (SIGN(half,a_cen) + SIGN(half,a_upw)) &
                  * MIN(ABS(a_upw), ABS(a_cen))
         !------------------------------------------------------------
         CASE(HR_SUPERBEE_LIMITER);  ! superbee
            a_lim = (SIGN(half,a_cen) + SIGN(half,a_upw))   &
                  *  MAX( MIN(  ABS(a_cen), 2*ABS(a_upw)),  &
                         MIN(2*ABS(a_cen),   ABS(a_upw)) )
         !------------------------------------------------------------
         CASE(HR_MONOCENTR_LIMITER);  ! Monotonized Central
            a_lim = MAX( zero,  MIN((a_cen+a_upw)/2, 2*a_cen, 2*a_upw) ) &
                  + MIN( zero,  MAX((a_cen+a_upw)/2, 2*a_cen, 2*a_upw) )
         !------------------------------------------------------------
         CASE(HR_VAN_ALBADA_LIMITER);  ! van Albada
            a_lim = a_cen * a_upw * (a_cen + a_upw) &
                  /(a_cen**2 + a_upw**2 + eps)
         !------------------------------------------------------------
         CASE DEFAULT
            WRITE(*,*) ' Limiter of unknown type. STOP'
            STOP
         !------------------------------------------------------------

      END SELECT


   END FUNCTION rebay_form_limiter
   !============================================================ 



   !============================================================ 
   FUNCTION rebay_form_centred_limiter(a_cen, a_bkw, a_fwd,    &
                                    limiter_type) RESULT (a_lim)
   !============================================================ 

 
      !------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN)  ::  a_bkw, a_cen, a_fwd
      INTEGER,      INTENT(IN)  ::  limiter_type
      
      REAL(KIND=8)              ::  a_lim
      !------------------------------------------------------------
      REAL(KIND=8)              ::  psi_fwd, psi_bkw, psi_cen
      REAL(KIND=8), PARAMETER   ::  zero = 0,  half = 0.5,  one = 1, &
                                    eps = 1.e-16
      !------------------------------------------------------------
  

      SELECT CASE(limiter_type)

         !------------------------------------------------------------
         CASE(HR_NO_LIMITER); 
            ! second order = a_cen
            a_lim = a_cen
         !------------------------------------------------------------
         CASE(HR_FULL_LIMITER);  
            ! no limiter, first-order upwind = 0.d0
            a_lim = 0.d0
         !------------------------------------------------------------
         CASE(HR_VAN_LEER_LIMITER);  
            ! van Leer = vanLeer(a_cen, a_bkw) 
            !          + vanLeer(a_cen, a_fwd)
            !          - a_cen 
            psi_bkw = (a_cen*ABS(a_bkw) + ABS(a_cen)*a_bkw) &
                    / (ABS(a_cen) + ABS(a_bkw) + eps)  
            psi_fwd = (a_cen*ABS(a_fwd) + ABS(a_cen)*a_fwd) &
                    / (ABS(a_cen) + ABS(a_fwd) + eps)  

            a_lim = psi_bkw + psi_fwd - a_cen
         !------------------------------------------------------------
         CASE(HR_MINMOD_LIMITER);  
            ! minmod = minmod(a_cen, b_bkw)
            !        + minmod(a_cen, b_fwd)   
            !        - a_cen       
            psi_bkw = (SIGN(half,a_bkw) + SIGN(half,a_cen)) &
                      * MIN(ABS(a_bkw), ABS(a_cen))         
            psi_fwd = (SIGN(half,a_fwd) + SIGN(half,a_cen)) & 
                      * MIN(ABS(a_fwd), ABS(a_cen))         

            a_lim = psi_bkw + psi_fwd - a_cen
         !------------------------------------------------------------
         CASE(HR_SUPERBEE_LIMITER);  
            ! superbee = superbee(a_cen, a_bkw)
            !          + superbee(a_cen, a_fwd)
            !          - a_cen
            psi_bkw = (SIGN(half,a_cen) + SIGN(half,a_bkw))   &
                    *  MAX( MIN(  ABS(a_cen), 2*ABS(a_bkw)),  &
                            MIN(2*ABS(a_cen),   ABS(a_bkw)) )       
            psi_fwd = (SIGN(half,a_cen) + SIGN(half,a_fwd))   & 
                    *  MAX( MIN(  ABS(a_cen), 2*ABS(a_fwd)),  &
                            MIN(2*ABS(a_cen),   ABS(a_fwd)) )

            a_lim = psi_bkw + psi_fwd - a_cen
         !------------------------------------------------------------
         CASE(HR_MONOCENTR_LIMITER);  
            ! Monotonized Central (MC) = MC(a_cen, a_bkw)
            !                          + MC(a_cen, a_fwd)
            !                          - a_cen
            psi_bkw = MAX( zero, MIN((a_cen+a_bkw)/2, 2*a_cen, 2*a_bkw) ) &
                    + MIN( zero, MAX((a_cen+a_bkw)/2, 2*a_cen, 2*a_bkw) ) 
            psi_fwd = MAX( zero, MIN((a_cen+a_fwd)/2, 2*a_cen, 2*a_fwd) ) &
                    + MIN( zero, MAX((a_cen+a_fwd)/2, 2*a_cen, 2*a_fwd) ) 

            a_lim = psi_bkw + psi_fwd - a_cen
         !------------------------------------------------------------
         CASE(HR_VAN_ALBADA_LIMITER);  
            ! van Albada = vanAlbada(a_cen, a_bkw) 
            !            + vanAlbada(a_cen, a_fwd) 
            !            - a_cen
            psi_bkw = a_cen * a_bkw * (a_cen + a_bkw) &
                    / (a_cen**2 + a_bkw**2 + eps) 
            psi_fwd = a_cen * a_fwd * (a_cen + a_fwd) &
                    / (a_cen**2 + a_fwd**2 + eps) 
                         
            a_lim = psi_bkw + psi_fwd - a_cen
         !------------------------------------------------------------
         CASE(HR_YEE_MINMOD_2);  
            ! minmod 2 = MIN( ABS(minmod(a_cen, a_bkw)), 
            !                 ABS(minmod(a_cen, a_fwd)) 
            !               ) * SIGN(a_cen) 
            psi_bkw = (SIGN(half,a_bkw) + SIGN(half,a_cen)) &
                      * MIN(ABS(a_bkw), ABS(a_cen))         
            psi_fwd = (SIGN(half,a_fwd) + SIGN(half,a_cen)) &
                      * MIN(ABS(a_fwd), ABS(a_cen))         
            
            a_lim = MIN( ABS(psi_bkw) , ABS(psi_fwd) ) &
                  * SIGN(1.d0, a_cen)
         !------------------------------------------------------------
         CASE(HR_YEE_MINMOD_3);  
            ! minmod 3 = 2*MIN( ABS(minmod(a_cen, a_bkw)), 
            !                   ABS(minmod(a_cen, a_fwd))  
            !                   ABS(minmod(a_cen, (a_bkw+a_fwd)/4)) 
            !                 ) * SIGN(a_cen) 
            psi_bkw = (SIGN(half,a_bkw) + SIGN(half,a_cen)) &
                      * MIN(ABS(a_bkw), ABS(a_cen))         
            psi_fwd = (SIGN(half,a_fwd) + SIGN(half,a_cen)) &
                      * MIN(ABS(a_fwd), ABS(a_cen))         
            psi_cen = (SIGN(half,(a_bkw+a_fwd)/4) + SIGN(half,a_cen)) &
                      * MIN(ABS((a_bkw+a_fwd)/4), ABS(a_cen))

            a_lim = 2 * MIN( ABS(psi_bkw) , ABS(psi_fwd), ABS(psi_cen) ) &
                  * SIGN(1.d0, a_cen)
         !------------------------------------------------------------
         CASE DEFAULT
            WRITE(*,*) ' Limiter of unknown type. STOP'
            STOP
         !------------------------------------------------------------


      END SELECT


   END FUNCTION rebay_form_centred_limiter
   !============================================================ 



   !============================================================ 
   FUNCTION standard_form_limiter(t, limiter_type) RESULT (lim)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN)  ::  t
      INTEGER,      INTENT(IN)  ::  limiter_type
      
      REAL(KIND=8)              ::  lim
      !------------------------------------------------------------ 
      REAL(KIND=8)              ::  c
      REAL(KIND=8), PARAMETER   ::  zero = 0,  half = 0.5,  one = 1
      !------------------------------------------------------------
  

      SELECT CASE(limiter_type)
      
         !------------------------------------------------------------
         CASE(HR_NO_LIMITER); ! second order
            lim = 1
         !------------------------------------------------------------
         CASE(HR_FULL_LIMITER);  ! no limiter, first-order upwind
            lim = 0
         !------------------------------------------------------------
         CASE(HR_VAN_LEER_LIMITER);  ! van Leer
            lim = (t + ABS(t))/(1.d0 + ABS(t)) 
         !------------------------------------------------------------
         CASE(HR_MINMOD_LIMITER);  ! minmod
            lim = MAX(0.d0, MIN(1.d0,t)) 
         !------------------------------------------------------------
         CASE(HR_SUPERBEE_LIMITER);  ! superbee
            lim = MAX(0.d0, MAX(MIN(1.d0, 2*t),MIN(2.d0, t)) ) 
         !------------------------------------------------------------
         CASE(HR_MONOCENTR_LIMITER);  ! Monotonized Central
            c = (1 + t) / 2 
            lim = MIN(ABS(c), 2.d0, 2.d0*ABS(t)) * SIGN(1.d0,c)
         !------------------------------------------------------------
         CASE(HR_VAN_ALBADA_LIMITER);  ! van Albada
            lim = (t + ABS(t)) / (1 + t**2) 
         !------------------------------------------------------------
         CASE DEFAULT
            WRITE(*,*) ' Limiter of unknown type. STOP'
            STOP
         !------------------------------------------------------------


      END SELECT


   END FUNCTION standard_form_limiter
   !============================================================ 
   


   !============================================================ 
   FUNCTION centred_limiter(t_bkw, t_fwd, limiter_type)   &
                                                    RESULT (lim)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN)  ::  t_bkw, t_fwd
      INTEGER,      INTENT(IN)  ::  limiter_type
           
      REAL(KIND=8)              ::  lim
      !------------------------------------------------------------
      REAL(KIND=8)              ::  c_t_bkw, c_t_fwd
      REAL(KIND=8), PARAMETER  ::  zero = 0,  half = 0.5,  one = 1
      !------------------------------------------------------------
  

      SELECT CASE(limiter_type)

         !------------------------------------------------------------
         CASE(HR_NO_LIMITER); ! second order
             lim = 1.d0
         !------------------------------------------------------------
         CASE(HR_FULL_LIMITER);  ! no limiter, first-order upwind
             lim = 0.d0
         !------------------------------------------------------------
         CASE(HR_VAN_LEER_LIMITER);  ! van Leer
             lim = (t_bkw+ABS(t_bkw))/(1.d0+ABS(t_bkw)) & 
                 + (t_fwd+ABS(t_fwd))/(1.d0+ABS(t_fwd)) - 1.d0
         !------------------------------------------------------------
         CASE(HR_MINMOD_LIMITER);  ! minmod
             lim = MAX(0.d0, MIN(1.d0,t_bkw)) &
                 + MAX(0.d0, MIN(1.d0,t_fwd)) -1.d0         
         !------------------------------------------------------------
         CASE(HR_SUPERBEE_LIMITER);  ! superbee
             lim =  MAX(0.d0,MIN(1.d0,2.d0*t_bkw),MIN(2.d0,t_bkw)) &
                 +  MAX(0.d0,MIN(1.d0,2.d0*t_fwd),MIN(2.d0,t_fwd)) &
                 -  1.d0
         !------------------------------------------------------------
         CASE(HR_MONOCENTR_LIMITER);  ! Monotonized Central
             c_t_bkw = .5d0*(1.d0+t_bkw)
             c_t_fwd = .5d0*(1.d0+t_fwd)
             lim = MIN(ABS(c_t_bkw),2.d0,2.d0*ABS(t_bkw))*SIGN(1.d0,c_t_bkw) &
                 + MIN(ABS(c_t_fwd),2.d0,2.d0*ABS(t_fwd))*SIGN(1.d0,c_t_fwd) &
                 - 1.d0
         !------------------------------------------------------------
         CASE(HR_VAN_ALBADA_LIMITER);  ! van Albada
         !------------------------------------------------------------
         CASE(HR_YEE_MINMOD_2);  ! minmod 2
             lim = MIN( MAX(0.d0, MIN(1.d0,t_bkw)) &
                      , MAX(0.d0, MIN(1.d0,t_fwd)) ) 
         !------------------------------------------------------------
         CASE(HR_YEE_MINMOD_3);  ! minmod 3
             lim = 2.d0*MIN( MAX( 0.d0, MIN(1.d0,t_bkw) ),           &
                             MAX( 0.d0, MIN(1.d0,t_fwd) ),           &
                             MAX( 0.d0, MIN(1.d0,(t_bkw+t_fwd)/4) )  )
         !------------------------------------------------------------
         CASE DEFAULT
            WRITE(*,*) ' Limiter of unknown type. STOP'
            STOP
         !------------------------------------------------------------


      END SELECT


   END FUNCTION centred_limiter
   !============================================================  


END MODULE high_resolution_fluxes

