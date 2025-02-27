!============================================================ 
!
!      Module: muscl_fluxes
!
! Description: MUSCL schemes 
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

MODULE muscl_fluxes


   !============================================================ 
   USE euler_equations
   USE thermodynamics
   USE roe_average,   ONLY: intermediate_ww, &
                      entropy_fix_eos

   USE high_resolution_fluxes,   ONLY: limiter_type
   !============================================================ 

   INTEGER, PARAMETER  ::  MUSCL_NO_LIMITER           = -1, &
                           MUSCL_FULL_LIMITER         =  0, &
                           MUSCL_VAN_LEER_LIMITER     =  1, &
                           MUSCL_MINMOD_LIMITER       =  2, &
                           MUSCL_SUPERBEE_LIMITER     =  3, &
                           MUSCL_MONOCENTR_LIMITER    =  4, &
                           MUSCL_VAN_ALBADA_LIMITER   =  5                        
                           
   CHARACTER(*), DIMENSION(-1:5), PARAMETER  ::  m_lim_names = &
                        (/ 'NO LIMITER: FULL II ORDER  ', &
                           'FULL LIMITER: I ORDER      ', &
                           'VAN LEER LIMITER           ', &
                           'MINMOD LIMITER             ', &
                           'SUPERBEE LIMITER           ', &
                           'MONOTONIZED CENTRAL LIMITER', &
                           'VAN ALBADA LIMITER         ' /)

   !============================================================ 
   ! Parameter of MUSCL reinterpolation
   ! ----------------------------------
   REAL(KIND=8), PARAMETER :: kappa   = 1.d0/3,    &
                              onemk_2 = 0.5*(1-kappa)
   !============================================================ 


   !============================================================ 
   ! Variables: private/public policy
   ! --------------------------------
   PUBLIC  ::  kappa, onemk_2
   !============================================================ 


   !============================================================ 
   ! Procedures: private/public policy
   ! ---------------------------------
   PUBLIC  ::  roe_MUSCL_flux, s_limiter
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   FUNCTION roe_muscl_flux (ww, eos, j_c_d, eta, Drr) RESULT (phi)
   !============================================================ 


      !------------------------------------------------------------
      ! MUSCL scheme
      ! Extrapolation of nonconservative variables:
      ! density, velocity, total enthalpy per unit volume
      !------------------------------------------------------------


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d
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
      TYPE(eos_type)                        ::  eos_i, eos_j 
      ! Roe scheme
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_
      TYPE(eos_ext_type)                    ::  eos_
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  RR, LL
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  lambda, lambda_ef 
      REAL(KIND=8)  ::  mod_eta
      ! MUSCL reconstruction
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js, &
                                                ww_ip, ww_jm
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  Dvv_ij, Dvv_iis, Dvv_jjs
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  vv_i,  vv_j,  &
                                                vv_is, vv_js, &
                                                vv_ip, vv_jm
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::   s_i,   s_j
      REAL(KIND=8) ::  Drr_ij, Drr_iis, Drr_jjs
      ! Indices
      INTEGER  ::  c, i, j, is, js, p
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
         !------------------------------------------------------------ 
         

         !============================================================
         ! MUSCL reconstruction
         ! ====================
         !
         !   M        U                S          C            L 
         !   Monotone Upstream-centred Scheme for Conservation Laws 
         !
         !          
         !                     |
         !                     |
         !                ww_i+|
         !                  ___|
         !              ___/   |
         !          ___/ ww_i  |ww_j-
         !      ___/           |____  ww_j
         !     /               |    \____
         !   ww_i*             |         \____   ww_j*
         !                     |              \____
         !     O---------O------------O-----------O
         !    i*         i   i+|j-    j          j*
         !                     |
         !                     |
         !                     


         !------------------------------------------------------------
         ! Conservative to primitive variables
         ! (density,velocity,total enthalpy per unit volume)
         ! abuse of notation: vv(p) = rho htot instead of htot
         ! ---------------------------------------------------
         p = SIZE(ww_i)

         vv_i(1)      = ww_i(1)
         vv_i(2:p-1)  = ww_i(2:p-1)/ww_i(1)
         vv_i(p)      = ww_i(p) + eos(i)%P

         vv_j(1)      = ww_j(1)
         vv_j(2:p-1)  = ww_j(2:p-1)/ww_j(1)
         vv_j(p)      = ww_j(p) + eos(j)%P

         vv_is(1)     = ww_is(1)
         vv_is(2:p-1) = ww_is(2:p-1)/ww_is(1)
         vv_is(p)     = ww_is(p) + eos(is)%P

         vv_js(1)     = ww_js(1)
         vv_js(2:p-1) = ww_js(2:p-1)/ww_js(1)
         vv_js(p)     = ww_js(p) + eos(js)%P
         !------------------------------------------------------------


         !------------------------------------------------------------
         ! Left interface
         ! --------------

         ! Variations of the primitive variables 
         ! -------------------------------------
         Dvv_ij  =                       vv_j - vv_i  
         Dvv_iis = (Drr_ij/Drr_iis) * (  vv_i - vv_is )
                   !^^^^^^^^^^^^^^ Normalization factor

         ! Computation of the limiter 
         ! --------------------------
         DO p = 1, SIZE(vv_i)
            s_i(p) = s_limiter(Dvv_iis(p), Dvv_ij(p), limiter_type(p))
         ENDDO 

         ! Conservative variables at the (left) interface
         ! ----------------------------------------------
         vv_ip = vv_i  +  0.5 * s_i * ( Dvv_ij - onemk_2*(Dvv_ij - Dvv_iis) ) 

         p = SIZE(vv_j)

         vv_ip(p) = vv_ip(p)/vv_ip(1) ! <<< rho htot

         ww_ip(1)     = vv_ip(1)
         ww_ip(2:p-1) = vv_ip(2:p-1)*vv_ip(1)
         ww_ip(p)     = vv_ip(p)*vv_ip(1) - Pi__vv(vv_ip)
         !------------------------------------------------------------


         !------------------------------------------------------------
         ! Right interface
         ! ---------------

         ! Variations of the conservative variables 
         ! ----------------------------------------
         Dvv_ij  =                       vv_j - vv_i  
         Dvv_jjs = (Drr_ij/Drr_jjs) * ( vv_js - vv_j )
                   !^^^^^^^^^^^^^^ Normalization factor

         ! Computation of the limiter 
         ! --------------------------
         DO p = 1, SIZE(vv_j)
            s_j(p) = s_limiter(Dvv_ij(p), Dvv_jjs(p), limiter_type(p))
         ENDDO 

         ! Conservative variables at the (right) interface
         ! -----------------------------------------------
         vv_jm = vv_j  -  0.5 * s_j * ( Dvv_ij + onemk_2*(Dvv_jjs - Dvv_ij) ) 

         p = SIZE(vv_j)

         vv_jm(p) = vv_jm(p)/vv_jm(1) ! <<< rho htot

         ww_jm(1)     = vv_jm(1)
         ww_jm(2:p-1) = vv_jm(2:p-1)*vv_jm(1)
         ww_jm(p)     = vv_jm(p)*vv_jm(1) - Pi__vv(vv_jm)
         !============================================================

         ! For ease of notation...
         ww_i = ww_ip;  ww_j = ww_jm
         eos_i = eos__e_v( ww_i(p) / ww_i(1) &
                           - 0.5*SUM( (ww_i(2:p-1)/ww_i(1))**2 ),&
                           1/ww_i(1))
         eos_j = eos__e_v( ww_j(p) / ww_j(1) &
                           - 0.5*SUM( (ww_j(2:p-1)/ww_j(1))**2 ),&
                           1/ww_j(1))

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

         ! Intermediate state of Roe in conservative variables
         ! ---------------------------------------------------
         CALL intermediate_ww (ww_i, eos_i, ww_j, eos_j, &
                               ww_, eos_) 
         ! Computation of the upwind correction
         ! ------------------------------------
         mod_eta = SQRT( SUM(etaij*etaij) );  eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_, eos_, eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_, eos_, eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = ABS( eigenvalues__ww_eos_nn(ww_, eos_, eta_ver) )

         lambda_ef = entropy_fix_eos(lambda,ww_, eos_, eta_ver)  
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


   END FUNCTION roe_muscl_flux 
   !============================================================ 



   !============================================================ 
   FUNCTION roe_muscl_flux_cons (ww, eos, j_c_d, eta, Drr) RESULT (phi)
   !============================================================ 


      !------------------------------------------------------------
      ! MUSCL scheme
      ! Extrapolation of conservative variable
      !------------------------------------------------------------


      !============================================================
      IMPLICIT NONE

      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)  ::  ww
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)  ::  eos
      INTEGER,            DIMENSION(:,:), INTENT(IN)  ::  j_c_d
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
      TYPE(eos_type)                        ::  eos_i, eos_j 
      ! Roe scheme
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_
      TYPE(eos_ext_type)                    ::  eos_
      REAL(KIND=8), DIMENSION(SIZE(ww,1),&
                              SIZE(ww,1))   ::  RR, LL
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  alpha
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  lambda, lambda_ef 
      REAL(KIND=8)  ::  mod_eta
      ! MUSCL reconstruction
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  ww_is, ww_js, &
                                                ww_ip, ww_jm
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::  Dww_ij, Dww_iis, Dww_jjs
      REAL(KIND=8), DIMENSION(SIZE(ww,1))   ::   s_i,   s_j
      REAL(KIND=8) ::  Drr_ij, Drr_iis, Drr_jjs
      ! Indices
      INTEGER  ::  c, i, j, is, js, p
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
          !------------------------------------------------------------ 
         

         !============================================================
         ! MUSCL reconstruction
         ! ====================
         !
         !   M        U                S          C            L 
         !   Monotone Upstream-centred Scheme for Conservation Laws 
         !
         !          
         !                     |
         !                     |
         !                ww_i+|
         !                  ___|
         !              ___/   |
         !          ___/ ww_i  |ww_j-
         !      ___/           |____  ww_j
         !     /               |    \____
         !   ww_i*             |         \____   ww_j*
         !                     |              \____
         !     O---------O------------O-----------O
         !    i*         i   i+|j-    j          j*
         !                     |
         !                     |
         !                     


         p = SIZE(ww_i)



         !------------------------------------------------------------
         ! Left interface
         ! --------------

         ! Variations of the conservative variables 
         ! -------------------------------------
         Dww_ij  =                       ww_j - ww_i  
         Dww_iis = (Drr_ij/Drr_iis) * (  ww_i - ww_is )
                   !^^^^^^^^^^^^^^ Normalization factor

         ! Computation of the limiter 
         ! --------------------------
         DO p = 1, SIZE(ww_i)
            s_i(p) = s_limiter(Dww_iis(p), Dww_ij(p), limiter_type(p))
         ENDDO 

         ! Conservative variables at the (left) interface
         ! ----------------------------------------------
         ww_ip = ww_i  +  0.5 * s_i * ( Dww_ij - onemk_2*(Dww_ij - Dww_iis) ) 
 
         !------------------------------------------------------------


         !------------------------------------------------------------
         ! Right interface
         ! ---------------

         ! Variations of the conservative variables 
         ! ----------------------------------------
         Dww_ij  =                       ww_j - ww_i  
         Dww_jjs = (Drr_ij/Drr_jjs) * ( ww_js - ww_j )
                   !^^^^^^^^^^^^^^ Normalization factor

         ! Computation of the limiter 
         ! --------------------------
         DO p = 1, SIZE(ww_j)
            s_j(p) = s_limiter(Dww_ij(p), Dww_jjs(p), limiter_type(p))
         ENDDO 

         ! Conservative variables at the (right) interface
         ! -----------------------------------------------
         ww_jm = ww_j  -  0.5 * s_j * ( Dww_ij + onemk_2*(Dww_jjs - Dww_ij) ) 

         !============================================================

         ! For ease of notation...
         ww_i = ww_ip;  ww_j = ww_jm
         eos_i = eos__e_v( ww_i(p) / ww_i(1) &
                           - 0.5*SUM( (ww_i(2:p-1)/ww_i(1))**2 ),&
                           1/ww_i(1))
         eos_j = eos__e_v( ww_j(p) / ww_j(1) &
                           - 0.5*SUM( (ww_j(2:p-1)/ww_j(1))**2 ),&
                           1/ww_j(1))

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

         ! Intermediate state of Roe in conservative variables
         ! ---------------------------------------------------
         CALL intermediate_ww (ww_i, eos_i, ww_j, eos_j, &
                               ww_, eos_) 
         ! Computation of the upwind correction
         ! ------------------------------------
         mod_eta = SQRT( SUM(etaij*etaij) ); eta_ver = etaij/mod_eta

         RR = right_eigenvectors__ww_eos_nn(ww_, eos_, eta_ver)
         LL =  left_eigenvectors__ww_eos_nn(ww_, eos_, eta_ver)
         alpha  = MATMUL(LL, ww_j - ww_i)
         lambda = ABS( eigenvalues__ww_eos_nn(ww_, eos_, eta_ver) )

         lambda_ef = entropy_fix_eos(lambda,ww_, eos_, eta_ver)  
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


   END FUNCTION roe_muscl_flux_cons 
   !============================================================ 


   !============================================================ 
   !============================================================ 
   !============================================================ 

   

   !============================================================ 
   FUNCTION s_limiter(a,b, limiter_type) RESULT (s)
   !============================================================ 


      ! Limiters in standard form


      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN)  ::  a,b
      INTEGER,      INTENT(IN)  ::  limiter_type
      
      REAL(KIND=8)              ::  s


      SELECT CASE(limiter_type)

         CASE(-1);  ! no limiter, linear reconstruction
                    ! ---------------------------------
             s = 1
 
         CASE(0);   ! full limiter, no reconstruction
                    ! -------------------------------
             s = 0
 
         CASE(1);   ! van Leer
                    ! --------
             s = 2 * ( a*b + ABS(a*b) ) / ( (a+b)**2 + 1.0d-8 )

         CASE(2);   ! minmod
                    ! ------
             s = 0.5 * (SIGN(1.d0,a) + SIGN(1.d0,b) ) * MIN( ABS(a), ABS(b) )

         CASE(3);   ! superbee
                    ! --------
             WRITE(*,*) ' Superbee limiter not implemented. STOP'
             STOP

         CASE(4);   ! Monotonized Central
                    ! -------------------
             WRITE(*,*) ' Monotonized Central limiter not implemented. STOP'
             STOP

         CASE(5);   ! van Albada
                    ! ----------
             s =  ( a*b + ABS(a*b) ) / ( a**2 + b**2 + 1.0d-12 )

         CASE DEFAULT
            WRITE (*,*) ''
            WRITE (*,*) 'ERROR. S_LIMITER:'
            WRITE (*,*) 'Limiter of unknown type'
            WRITE (*,*) ''
            STOP

      END SELECT


   END FUNCTION s_limiter
   !============================================================ 



   !============================================================ 
   !============================================================ 
   !============================================================ 



!    !============================================================ 
!    SUBROUTINE read_param_muscl_fluxes(idf)
!    !============================================================ 
! 
! 
!       !------------------------------------------------------------
!       IMPLICIT NONE
! 
!       INTEGER  ::  idf
!       !------------------------------------------------------------
!       INTEGER  ::  ww_size
!       !------------------------------------------------------------
! 
! 
!           READ (idf,*) ww_size
!           ALLOCATE (limiter_type(ww_size))         
!           READ (idf,*) limiter_type
! 
!    END SUBROUTINE read_param_muscl_fluxes
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE write_param_muscl_fluxes(idf)
!    !============================================================ 
! 
! 
!       !------------------------------------------------------------
!       IMPLICIT NONE
! 
!       INTEGER  ::  idf
!       !------------------------------------------------------------
!       INTEGER  ::  p
!       !------------------------------------------------------------
! 
! 
!       WRITE(idf,*) '   PARAMETERS FOR MUSCL NUMERICAL FLUXES'
!       
!       DO p = 1, SIZE(limiter_type)
!          WRITE(idf,*) '   Limiter ', p,': ', m_lim_names(limiter_type(p))
!       ENDDO
!       
!       WRITE(idf,*)
!       
! 
!    END SUBROUTINE write_param_muscl_fluxes
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE init_muscl_fluxes
!    !============================================================ 
! 
! 
!       !------------------------------------------------------------
!       IMPLICIT NONE
!       !------------------------------------------------------------
! 
! 
! 
!    END SUBROUTINE init_muscl_fluxes
!    !============================================================ 


END MODULE muscl_fluxes

