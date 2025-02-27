!============================================================ 
!
!      Module: euler_num_fluxes
!
! Description: driver for the different space discretizion
!              schemes for the evaluation of Euler numerical
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

  MODULE euler_num_fluxes

  !----------------------------------------------------------------------------
  USE structures
  
  !----------------------------------------------------------------------------
  IMPLICIT NONE

  PRIVATE

  ! Parameters 
  INTEGER :: flux_type    
  
  ! Node-pair connectivity and metric quantities 
  INTEGER,      DIMENSION(:,:), POINTER :: j_c_d
  INTEGER,      DIMENSION(:,:), POINTER :: cs_c_d
  REAL(KIND=8), DIMENSION(:,:), POINTER :: Drr
  REAL(KIND=8), DIMENSION(:,:), POINTER :: eta
  
  ! Roe average
  REAL(KIND=8),       DIMENSION(:,:), ALLOCATABLE :: ww_
  TYPE(eos_ext_type), DIMENSION(:),   ALLOCATABLE :: eos_
  
  ! HR and MUSCL methods
  REAL(KIND=8),       DIMENSION(:,:), ALLOCATABLE :: phie
  REAL(KIND=8),       DIMENSION(:,:), ALLOCATABLE :: wwe
  TYPE(eos_type),     DIMENSION(:),   ALLOCATABLE :: eose
  REAL(KIND=8),       DIMENSION(:,:), ALLOCATABLE :: wwe_
  TYPE(eos_ext_type), DIMENSION(:),   ALLOCATABLE :: eose_

  INTEGER :: ext_size    

  INTEGER, PARAMETER :: FO_ROE           = 1, &
                        MUSCL_ROE_CONS   = 2, &
                        MUSCL_ROE        = 3, &
                        HR_ROE_GALERKIN  = 4, &
                        HR_ROE_GAL_CEN   = 5, &
                        HR_ROE_LW        = 6, &
                        HR_ROE_LW_ROT    = 7

  CHARACTER(*), DIMENSION(1:7), PARAMETER  ::  ft_names = &
                       (/ 'ROE SCHEME, FIRST ORDER UPWIND            ', &
                          'ROE SCHEME, MUSCL FOR CONSERVATIVE VARS   ', &
                          'ROE SCHEME, MUSCL FOR rho, u, rho htot    ', &
                          'ROE SCHEME + GALERKIN, HIGH RESOLUTION    ', &
                          'ROE SCHEME + GALERKIN, HR, CENTRED LIMITER', &
                          'ROE SCHEME + LAX-WENDROFF, HIGH RESOLUTION', &
                          'ROE SCHEME + LW ROTATED, HIGH RESOLUTION  ' /)

  !----------------------------------------------------------------------------
  PUBLIC   ::  euler_num_flux, euler_implicit_num_flux, &
               euler_compute_dt,                        &
               euler_compute_CFL,                       &
               read_param_euler_num_fluxes,             &
               write_param_euler_num_fluxes,            &
               init_euler_num_fluxes,                   &
               set_eul_full_limiters,                   &
               set_eul_hr_limiters, flux_type
  !----------------------------------------------------------------------------

  CONTAINS

 
   SUBROUTINE  read_param_euler_num_fluxes(idf)
   !----------------------------------------------------------------------------
   USE roe_average,              ONLY: read_param_roe_average
   USE high_resolution_fluxes,   ONLY: read_param_hr_fluxes
   USE extrap_boundary,          ONLY: read_param_extrap_boundary
   
   IMPLICIT NONE
   INTEGER :: idf
   !----------------------------------------------------------------------------

   READ(idf,*) flux_type
         
   ! Read parameters for Roe linearization
   CALL read_param_roe_average(idf)
   
   ! Read parameters for numerical fluxes
   SELECT CASE (flux_type)
   
      CASE (MUSCL_ROE_CONS,  &
            MUSCL_ROE,       &
            HR_ROE_GALERKIN, &
            HR_ROE_GAL_CEN,  &
            HR_ROE_LW,       &
            HR_ROE_LW_ROT)
            
         CALL read_param_hr_fluxes(idf)

   END SELECT
   
   ! Read parameters for extrapolation at boundaries
   SELECT CASE (flux_type)

      CASE (MUSCL_ROE_CONS,  &
            MUSCL_ROE,       &
            HR_ROE_GALERKIN, &
            HR_ROE_GAL_CEN,  &
            HR_ROE_LW,       &
            HR_ROE_LW_ROT)
         
         CALL read_param_extrap_boundary(idf)
   
   END SELECT

   END SUBROUTINE  read_param_euler_num_fluxes





   SUBROUTINE  init_euler_num_fluxes(fem_metric, ww)
   !------------------------------------------------------------
   USE node_pair_structure,  ONLY: j_c_d_fem => j_c_d,    &
                                   j_c_d_fvm => j_c_fv,   & 
                                   cs_c_d_fem => cs_c_d,  &
                                   cs_c_d_fvm => cs_c_fv

   USE metric_coefficients,  ONLY: eta_fem => eta,     &   
                                   eta_fvm => eta_fv,  &
                                   Drr_fem => Drr,     &                                    
                                   Drr_fvm => Drr_fv
   
   USE extrap_boundary,      ONLY: init_extrap_boundary, ww_ext_size
   USE nodes,                ONLY: jd_jb, bound_p
   USE roe_average,          ONLY: init_roe_average
   USE mp_interface
   
   !------------------------------------------------------------   
   IMPLICIT NONE
   
   LOGICAL,                      INTENT(IN) :: fem_metric
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
   !------------------------------------------------------------


   ! Selects between FEM and FVM topology and metrics
   IF (fem_metric) THEN 
      
      j_c_d  =>  j_c_d_fem
      cs_c_d => cs_c_d_fem
      eta    =>    eta_fem
      Drr    =>    Drr_fem
      
   ELSE
      
      j_c_d  =>  j_c_d_fvm
      cs_c_d => cs_c_d_fvm
      eta    =>    eta_fvm
      Drr    =>    Drr_fvm

      IF( (flux_type == HR_ROE_LW) .OR. &
          (flux_type == HR_ROE_LW_ROT)  ) THEN

         WRITE(*,*)
         WRITE(*,*) 'ERROR. INIT_EULER_NUM_FLUXES:'
         WRITE(*,*) 'Lax Wendroff scheme and FVM metrics not'
         WRITE(*,*) 'compatible, yet.'
         WRITE(*,*)
         STOP
      
      ENDIF

   ENDIF      
     
   ! Init extrapolation at boundaries
   SELECT CASE (flux_type)
   
      CASE (MUSCL_ROE_CONS,  &
            MUSCL_ROE,       &
            HR_ROE_GALERKIN, &
            HR_ROE_GAL_CEN,  &
            HR_ROE_LW,       &
            HR_ROE_LW_ROT)
         
         IF (.NOT. MP_job  .OR.  NjP_b /= 0) THEN
           CALL init_extrap_boundary(ww, j_c_d, cs_c_d, jd_jb, bound_p, Drr, fem_metric)	 
	   ext_size = ww_ext_size()
	 ELSE
	   ext_size = 0
	 ENDIF  
         
         ALLOCATE (phie(SIZE(ww,1), SIZE(ww,2) + ext_size),    &
                   wwe(SIZE(ww,1), SIZE(ww,2) + ext_size),     &
                   eose(SIZE(ww,2) + ext_size),                &
                   wwe_(SIZE(ww,1), SIZE(j_c_d,2) + ext_size), &
                   eose_(SIZE(j_c_d,2) + ext_size))           

   END SELECT


   ! Init Roe average
   ALLOCATE (ww_(SIZE(ww,1), SIZE(j_c_d,2)), eos_(SIZE(j_c_d,2)))

   CALL init_roe_average

   END SUBROUTINE  init_euler_num_fluxes





   SUBROUTINE  set_eul_full_limiters
   !-----------------------------------------------------------------------
   USE high_resolution_fluxes,   ONLY: limiter_type, HR_FULL_LIMITER

   IMPLICIT NONE 
   !-----------------------------------------------------------------------

   SELECT CASE (flux_type)

      CASE (MUSCL_ROE_CONS,  &
            MUSCL_ROE,       &
            HR_ROE_GALERKIN, &
            HR_ROE_GAL_CEN,  &
            HR_ROE_LW,       &
            HR_ROE_LW_ROT)

        limiter_type = HR_FULL_LIMITER

   END SELECT

   END SUBROUTINE  set_eul_full_limiters





   SUBROUTINE  set_eul_hr_limiters
   !-----------------------------------------------------------------------
   USE high_resolution_fluxes,   ONLY: limiter_type, hr_limiter_type

   IMPLICIT NONE 
   !-----------------------------------------------------------------------

   SELECT CASE (flux_type)

      CASE (MUSCL_ROE_CONS,  &
            MUSCL_ROE,       &
            HR_ROE_GALERKIN, &
            HR_ROE_GAL_CEN,  &
            HR_ROE_LW,       &
            HR_ROE_LW_ROT)

        limiter_type = hr_limiter_type

   END SELECT

   END SUBROUTINE  set_eul_hr_limiters





   SUBROUTINE  write_param_euler_num_fluxes(idf)
   !------------------------------------------------------------
   USE high_resolution_fluxes,  ONLY: write_param_hr_fluxes
   USE extrap_boundary,         ONLY: write_param_extrap_boundary
   USE roe_average,             ONLY: write_param_roe_average
   
   IMPLICIT NONE
   INTEGER  ::  idf
   !------------------------------------------------------------


   WRITE(idf,*) '   PARAMETERS FOR EULER NUMERICAL FLUXES'
   WRITE(idf,*) '   Numerical flux type: ', ft_names(flux_type)
   WRITE(idf,*)

   
   ! Write parameters for numerical fluxes
   SELECT CASE (flux_type)
   
      CASE (MUSCL_ROE_CONS,  &
            MUSCL_ROE,       &
            HR_ROE_GALERKIN, &
            HR_ROE_GAL_CEN,  &
            HR_ROE_LW,       &
            HR_ROE_LW_ROT)
            
         CALL write_param_hr_fluxes(idf)
       
   END SELECT
   
   ! Write parameters for extrapolation at boudnaries
   SELECT CASE (flux_type)

      CASE (MUSCL_ROE_CONS,  &
            MUSCL_ROE,       &
            HR_ROE_GALERKIN, &
            HR_ROE_GAL_CEN,  &
            HR_ROE_LW,       &
            HR_ROE_LW_ROT)
         
         CALL write_param_extrap_boundary(idf)
   
   END SELECT
         
   ! Write parameters for Roe linearization
   CALL write_param_roe_average(idf)

   END SUBROUTINE  write_param_euler_num_fluxes

 



   FUNCTION  euler_num_flux(fem_bc, ww, eos, wb, dt, stage)   RESULT(phi)
   !------------------------------------------------------------------
   USE roe_average,              ONLY: intermediate_ww
   USE extrap_boundary,          ONLY: ww_extrap, ww_extrap_roe
   USE euler_equations,          ONLY: eos__ww, flux__ww
   USE first_order_fluxes,       ONLY: Roe_FO_flux
   USE high_resolution_fluxes,   ONLY: Roe_Galerkin_HR_flux,     &
                                       Roe_Galerkin_HR_flux_cen, &
                                       Roe_LaxWendroff_HR_flux,  &
                                       Roe_LaxWendroff_HR_flux_rot
   USE muscl_fluxes,             ONLY: Roe_MUSCL_flux,  &
                                       Roe_MUSCL_flux_cons
   
   USE nodes,                    ONLY: jd_jb
   USE node_pair_structure,      ONLY: j_c_b
   USE metric_coefficients,      ONLY: chi_b, xi_bp, stiff, stiff_T
   use mp_interface
   !------------------------------------------------------------------   
   IMPLICIT NONE

   LOGICAL,                        INTENT(IN) :: fem_bc
   REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: ww
   TYPE(eos_type), DIMENSION(:),   INTENT(IN) :: eos
   REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: wb
   REAL(KIND=8),   DIMENSION(:),   INTENT(IN) :: dt
   INTEGER,                        INTENT(IN) :: stage

   REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(ww,2)) :: phi
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: phi_ij, phi_i

   REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(xi_bp,1)) :: ff_i, ff_j

   REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ver
   REAL(KIND=8)  ::  mod_eta
      
   INTEGER :: c, i, j, p
   !------------------------------------------------------------------
      
   phi = 0

   ! Compute Roe's intermediate state
   IF (stage > 1) THEN
     
     ! Already computed in compute_dt
     DO c = 1, SIZE(j_c_d,2)                                                      

       i = j_c_d(1,c)
       j = j_c_d(2,c)
       
       mod_eta = SQRT(SUM(eta(:,c)**2))
       eta_ver = eta(:,c)/mod_eta

       ! Intermediate state for each couple of nodes
       CALL intermediate_ww (ww(:,i), eos(i), ww(:,j), eos(j), ww_(:,c), eos_(c))
   
     ENDDO                                                                        
   
   ENDIF


   
   ! Extrapolation at boundaries      
   IF (ext_size > 0) THEN 
      
      phie = 0 
      
      wwe(:,            1:SIZE(ww,2) ) = ww
      wwe(:, SIZE(ww,2)+1:SIZE(wwe,2)) = ww_extrap(ww, j_c_d, Drr)
               
      eose(1:SIZE(ww,2)) = eos
      DO i = SIZE(ww,2)+1, SIZE(wwe,2)
        eose(i) = eos__ww(wwe(:,i))
      ENDDO
      
       wwe_(:, 1:SIZE(ww_,2)) = ww_
      eose_(   1:SIZE(ww_,2)) = eos_
      
      CALL ww_extrap_roe(j_c_d, cs_c_d, wwe, eose,  wwe_, eose_)

   ELSE
   
      wwe   = ww
      eose  = eos
      wwe_  = ww_
      eose_ = eos_

   ENDIF   
   
   
   ! DOMAIN INTEGRAL

   SELECT CASE (flux_type)

      ! Roe's first order upwind
      CASE (FO_ROE)
      
         phi = Roe_FO_flux(ww, eos, ww_, eos_, j_c_d, eta) 


      ! Roe's method with MUSCL
      CASE (MUSCL_ROE)
                  
         phie = Roe_MUSCL_flux(wwe, eose, j_c_d, eta, Drr) 
         
         DO i = 1, SIZE(ww,2)
            phi(:,i) = phie(:,i) 
         ENDDO


      ! Roe's method with MUSCL on conservative variables
      CASE (MUSCL_ROE_CONS)
                  
         phie = Roe_MUSCL_flux_cons (wwe, eose, j_c_d, eta, Drr)
         
         DO i = 1, SIZE(ww,2)
             phi(:,i) = phie(:,i) 
         ENDDO


      ! Roe + Galerkin HR method
      CASE (HR_ROE_GALERKIN)
      
         phie = Roe_Galerkin_HR_flux(wwe, eose, wwe_, eose_, &
                                     j_c_d, cs_c_d, eta, Drr) 
         
         DO i = 1, SIZE(ww,2)
             phi(:,i) = phie(:,i) 
         ENDDO


      ! Roe + Galerkin HR method centred scheme         
      CASE (HR_ROE_GAL_CEN)
      
         phie = Roe_Galerkin_HR_flux_cen(wwe, eose, wwe_, eose_, &
                                   j_c_d, cs_c_d, eta, Drr)
         DO i = 1, SIZE(ww,2)
             phi(:,i) = phie(:,i) 
         ENDDO


      ! Roe + Lax-Wendroff HR method
      CASE (HR_ROE_LW)

         phie = Roe_LaxWendroff_HR_flux(wwe, eose, wwe_, eose_, dt,   &
                                        j_c_d, cs_c_d, eta, stiff_T, Drr)
         DO i = 1, SIZE(ww,2)
             phi(:,i) = phie(:,i) 
         ENDDO


      ! Roe + Lax-Wendroff HR method (rotated?)
      CASE (HR_ROE_LW_ROT)

         phie = Roe_LaxWendroff_HR_flux_rot(wwe, eose, wwe_, eose_, dt, &
                                            j_c_d, cs_c_d, eta, stiff, Drr)
         DO i = 1, SIZE(ww,2)
             phi(:,i) = phie(:,i) 
         ENDDO


      ! Unknown method
      CASE DEFAULT
      
         WRITE(*,*)
         WRITE(*,*) 'ERROR. Unknown numerical flux.'
         WRITE(*,*);   STOP

   END SELECT


   ! BOUNDARY INTEGRALS:
   
   ! Node-pair contributions
   IF (fem_bc) THEN

      DO c = 1, SIZE(chi_b,2) 

         i = j_c_b(1,c);  j = j_c_b(2,c)
      
         ff_i = flux__ww(wb(:,i)); ff_j = flux__ww(wb(:,j))
      
         DO p = 1, SIZE(wb,1)
           phi_ij(p) = SUM( chi_b(:,c) * (ff_j(p,:) - ff_i(p,:)) ) / 2
         ENDDO

         phi(:,jd_jb(i)) = phi(:,jd_jb(i)) + phi_ij
         phi(:,jd_jb(j)) = phi(:,jd_jb(j)) - phi_ij
 
      ENDDO

   ENDIF

   ! Nodal contributions   
   DO i = 1, SIZE(xi_bp,2)
   
      ff_i = flux__ww(wb(:,i))

      DO p = 1, SIZE(wb,1)
        phi_i(p) = SUM(xi_bp(:,i) * ff_i(p,:))
      ENDDO

      phi(:,jd_jb(i)) = phi(:,jd_jb(i)) + phi_i

   ENDDO


   END FUNCTION  euler_num_flux
      
      
      


   FUNCTION  euler_implicit_num_flux(fem_bc, ww, eos, wb, dwb_dwi, dt, stage, MM)   RESULT(phi)
   !-----------------------------------------------------------------------------------------
   USE csr
   USE csr_pair_sys   
   USE roe_average,              ONLY: intermediate_ww
   USE extrap_boundary,          ONLY: ww_extrap, ww_extrap_roe
   USE euler_equations,          ONLY: eos__ww, flux__ww
   USE first_order_fluxes,       ONLY: implicit_Roe_FO_flux
   USE high_resolution_fluxes,   ONLY: imp_Roe_Galerkin_HR_flux,     &
                                       imp_Roe_Galerkin_HR_flux_cen, &
                                       imp_Roe_LaxWendroff_HR_flux

   USE nodes,                    ONLY: jd_jb
   USE node_pair_structure,      ONLY: j_c_b
   USE metric_coefficients,      ONLY: chi_b, xi_bp, stiff_T                                       
   USE euler_equations,          ONLY: jacobian__ww_nn
   
   !-----------------------------------------------------------------------------------------
   IMPLICIT NONE

   LOGICAL,                          INTENT(IN) :: fem_bc
   REAL(KIND=8),   DIMENSION(:,:),   INTENT(IN) :: ww
   TYPE(eos_type), DIMENSION(:),     INTENT(IN) :: eos
   REAL(KIND=8),   DIMENSION(:,:),   INTENT(IN) :: wb
   REAL(KIND=8),   DIMENSION(:,:,:), INTENT(IN) :: dwb_dwi 
   REAL(KIND=8),   DIMENSION(:),     INTENT(IN) :: dt
   INTEGER,                          INTENT(IN) :: stage
   
   TYPE(CSR_matrix), INTENT(INOUT) :: MM

   REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(ww,2)) :: phi
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: phi_ij, phi_i
   
   REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ver
   REAL(KIND=8) :: mod_eta
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(xi_bp,1)) :: ff_i, ff_j
   REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,1)) :: AA_xi, MM_ij

   INTEGER :: c, i, j, p
   !-----------------------------------------------------------------------------------------
    
   !============================================================ 
   ! Domain contribution
   ! ===================      
   phi = 0.d0
         
   ! Compute Roe's intermediate state
   IF (stage > 1) THEN  ! Already computed in compute_dt
     DO c = 1, SIZE(j_c_d,2)
   
        i = j_c_d(1,c)
        j = j_c_d(2,c)

        mod_eta = SQRT(SUM(eta(:,c)**2))
        eta_ver = eta(:,c)/mod_eta

        CALL intermediate_ww (ww(:,i), eos(i), ww(:,j), eos(j), &
                              ww_(:,c), eos_(c) ) 
   
     ENDDO
   ENDIF

   IF (ext_size > 0) THEN 
      
      phie = 0.d0 
      
      wwe(:, 1:SIZE(ww,2)) =  ww
      wwe(:, SIZE(ww,2)+1:SIZE(wwe,2)) = ww_extrap(ww, j_c_d, Drr)
 
      eose(1:SIZE(ww,2)) = eos
      
      DO i = SIZE(ww,2)+1, SIZE(wwe,2)
        eose(i) = eos__ww(wwe(:,i))
      ENDDO
      
      wwe_(:, 1:SIZE(ww_,2)) = ww_
      eose_(1:SIZE(ww_,2)) = eos_
      CALL ww_extrap_roe(j_c_d, cs_c_d, wwe, eose,  wwe_, eose_)

   ENDIF

   SELECT CASE (flux_type)


      !------------------------------------------------------------
      ! Roe's first order upwind
      ! ------------------------
      CASE ( FO_ROE )
      
         phi = implicit_Roe_FO_flux ( ww, eos, ww_, eos_, j_c_d, eta, &
                                      MM)
      !------------------------------------------------------------



      !------------------------------------------------------------
      ! Roe + Galerkin HR method
      ! ------------------------
      CASE ( HR_ROE_GALERKIN )
      
         phie = imp_Roe_Galerkin_HR_flux (wwe, eose, wwe_, eose_, &
                   j_c_d, cs_c_d, eta, Drr,  MM)
         DO i = 1, SIZE(ww,2)
             phi(:,i) = phie(:,i) 
         ENDDO
        
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! Roe + Galerkin HR method
      ! ------------------------
      CASE ( HR_ROE_LW )

         phie = imp_Roe_LaxWendroff_HR_flux (wwe, eose, wwe_, eose_, dt, &
              j_c_d, cs_c_d, eta, stiff_T,  Drr,  MM)
         DO i = 1, SIZE(ww,2)
             phi(:,i) = phie(:,i) 
         ENDDO
 
        
      !------------------------------------------------------------
      ! Roe + Galerkin HR method centred limiter
      ! ----------------------------------------
      CASE ( HR_ROE_GAL_CEN )
       
         phie = imp_Roe_Galerkin_HR_flux_cen(wwe, eose, wwe_, eose_, &
                                    j_c_d, cs_c_d, eta, Drr,  MM)
         DO i = 1, SIZE(ww,2)
              phi(:,i) = phie(:,i) 
         ENDDO

      !------------------------------------------------------------
      ! Unknown method
      ! --------------
      CASE DEFAULT
      
         WRITE(*,*) ' euler_num_fluxes: method of unknown type. STOP'
         STOP
        
      !------------------------------------------------------------

   END SELECT
   !============================================================ 


   !============================================================ 
   ! Boundary contribution to the numerical flux
   ! ===========================================

   !------------------------------------------------------------
   ! Node-pair contributions
   ! -----------------------
   IF( fem_bc ) THEN

      DO c = 1, SIZE(chi_b,2) 

         i = j_c_b(1,c);  j = j_c_b(2,c)
      
         ff_i = flux__ww(wb(:,i)); ff_j = flux__ww(wb(:,j))
      
         DO p = 1, SIZE(wb,1)
            phi_ij(p) = SUM( chi_b(:,c) * (ff_j(p,:) - ff_i(p,:)) ) / 2
         ENDDO

         phi(:,jd_jb(i)) = phi(:,jd_jb(i)) + phi_ij
         phi(:,jd_jb(j)) = phi(:,jd_jb(j)) - phi_ij
 
      ENDDO

   ENDIF      
   !------------------------------------------------------------
   
   !------------------------------------------------------------
   ! Nodal contributions
   ! -------------------
   DO i = 1, SIZE(xi_bp,2)
   
      ff_i = flux__ww(wb(:,i))
   
      DO p = 1, SIZE(wb,1)
         phi_i(p) = SUM( xi_bp(:,i) * ff_i(p,:) )
      ENDDO

      phi(:,jd_jb(i)) = phi(:,jd_jb(i)) + phi_i

      ! Implicit bc
      AA_xi = jacobian__ww_nn(wb(:,i),xi_bp(:,i))         
      MM_ij  = MATMUL(AA_xi, dwb_dwi(:,:,i))
      
      CALL add_CSR_ij_sys(jd_jb(i), jd_jb(i), MM_ij, MM)


   ENDDO
   
   END FUNCTION  euler_implicit_num_flux





   FUNCTION euler_compute_dt(ww, eos, CFL, variable_dt)  RESULT(dt)
   !
   !  Compute the time step.  Local time stepping (steady flow
   !  only) or constant time step are allowed.
   !------------------------------------------------------------
   USE roe_average,          ONLY: intermediate_ww
   USE euler_equations,      ONLY: eigenvalues__ww_eos_nn
   USE metric_coefficients,  ONLY: cell   
   USE mp_interface
   
   !------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: ww
   TYPE(eos_type), DIMENSION(:),   INTENT(IN) :: eos
   REAL(KIND=8),                   INTENT(IN) :: CFL 
   LOGICAL,                        INTENT(IN) :: variable_dt  

   REAL(KIND=8), DIMENSION(SIZE(ww,2))  :: dt    
   REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ver
   REAL(KIND=8), DIMENSION(MP_wProc)    :: MP_minDt
   
   REAL(KIND=8) :: l_max, mod_eta, minDt
   
   INTEGER :: c, i, j
   !------------------------------------------------------------ 

   ! [Selmin: Bridge...]
   !
   !                     2 |C_i|
   !  dt_i = CFL -----------------------,   for each i in K_i
   !              INT_dC_i [lambda_max]

   dt = 0.d0
   DO c = 1, SIZE(j_c_d,2) 

     i = j_c_d(1,c);  j = j_c_d(2,c)

     mod_eta = SQRT(SUM(eta(:,c)**2 ))
     eta_ver = eta(:,c)/mod_eta

     CALL intermediate_ww(ww(:,i), eos(i), ww(:,j), eos(j), ww_(:,c), eos_(c)) 
    			    
     l_max = MAXVAL(ABS(eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)))

     l_max = l_max * mod_eta

     dt(i) = dt(i) + l_max
     dt(j) = dt(j) + l_max

   ENDDO


   ! Local time step (The sum of 1.d-8
   ! is an artifact necessary only in multi-processor
   ! computation. It affects only the time step
   ! relative to the extended nodes whose
   ! time advance is not relevant)
   IF (.NOT. MP_job) THEN
     dt = CFL * 2*cell / dt
   ELSE
     dt = CFL * 2*cell / (dt+1.d-10)   
   ENDIF  
   

   ! Global time step  --->  dt = MIN dt(i),   for all i in K_i
   IF (.NOT. MP_job  .AND. .NOT.variable_dt)  dt = MINVAL(dt)
   
   IF (MP_job  .AND. .NOT.variable_dt) THEN
     
     minDt = MINVAL(dt(1:NjP))
     
     CALL MP_gather_all_minDt(minDt, MP_minDt)
     
     dt = MINVAL(MP_minDt)
     
   ENDIF  

   END FUNCTION euler_compute_dt





   FUNCTION euler_compute_CFL(ww, eos, dt)  RESULT(iCFL)
   !------------------------------------------------------------
   USE roe_average,          ONLY: intermediate_ww
   USE euler_equations,      ONLY: eigenvalues__ww_eos_nn
   USE metric_coefficients,  ONLY: cell   
      
   !------------------------------------------------------------ 
   IMPLICIT NONE

   REAL(KIND=8),   DIMENSION(:,:),        INTENT(IN)  ::  ww   
   TYPE(eos_type), DIMENSION(:),          INTENT(IN)  ::  eos  
   REAL(KIND=8),   DIMENSION(SIZE(ww,2)), INTENT(IN)  ::  dt

   REAL(KIND=8) :: iCFL
   
   REAL(KIND=8), DIMENSION(SIZE(ww,2))  :: vCFL
   REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ver
   
   REAL(KIND=8) :: l_max, mod_eta
   INTEGER :: c, i, j
   !------------------------------------------------------------

   DO c = 1, SIZE(j_c_d,2) 

     i = j_c_d(1,c)
     j = j_c_d(2,c)

     mod_eta = SQRT(SUM(eta(:,c)**2 ))
     eta_ver = eta(:,c)/mod_eta
     
     CALL intermediate_ww(ww(:,i), eos(i), ww(:,j), eos(j), ww_(:,c), eos_(c)) 
    			    
     l_max = MAXVAL(ABS(eigenvalues__ww_eos_nn(ww_(:,c), eos_(c), eta_ver)))

     l_max = l_max * mod_eta
     
     vCFL(i) = vCFL(i) + l_max
     vCFL(j) = vCFL(j) + l_max
    
   ENDDO

   vCFL = vCFL * dt / (2*cell)

   iCFL = MAXVAL(vCFL)

   END FUNCTION euler_compute_CFL


  END MODULE euler_num_fluxes
