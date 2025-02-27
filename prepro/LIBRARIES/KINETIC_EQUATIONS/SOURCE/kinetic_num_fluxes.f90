!=============================================================================== 
!
!      Module: kinetic_num_fluxes
!
! Description: Compute interface fluxes on the basis of the 
!              kinetic theory of gases.
!
!      Author: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!===============================================================================

  MODULE kinetic_num_fluxes

  !-----------------------------------------------------------------------------
  USE nodes,		    ONLY: bound_p
  USE node_pair_structure,  ONLY:  j_c_d_fvm =>  j_c_fv,  &
        			  cs_c_d_fvm => cs_c_fv
  USE metric_coefficients,  ONLY: eta_fvm => eta_fv,  &
                                  Drr_fvm => Drr_fv, cell

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: KFVS_FO = 1, &
                        BGK_FO  = 2, &
                        KFVS_HR = 3, &
                        BGK_HR  = 4
			
  CHARACTER(*), DIMENSION(1:4), PARAMETER :: ft_names =  &
                     (/ 'KFVS SCHEME, FIRST ORDER', &
                        'BGK SCHEME, FIRST ORDER ', &
                        'KFVS SCHEME, HI-RES     ', &
                        'BGK SCHEME, HI-RES      '/)

  INTEGER :: num_flux, varNumFlux

  ! Node-pairs connectivity and metric quantities 
  INTEGER,      DIMENSION(:,:), POINTER :: j_c_d
  INTEGER,      DIMENSION(:,:), POINTER :: cs_c_d
  REAL(KIND=8), DIMENSION(:,:), POINTER :: Drr
  REAL(KIND=8), DIMENSION(:,:), POINTER :: eta

  ! HR methods
  INTEGER  ::  ext_size
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: phie
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: wwe


  PUBLIC :: kinetic_domain_flux,             &
            read_param_kinetic_num_fluxes,   &
            write_param_kinetic_num_fluxes,  &
            init_kinetic_num_fluxes,         &
            kinetic_compute_dt,              &
            kinetic_compute_CFL,             &
            switch_first_order,              &
            switch_second_order
  !---------------------------------------------------------------------------------

  CONTAINS


  SUBROUTINE  read_param_kinetic_num_fluxes(flow_rgme, idf)
  !---------------------------------------------------------------------------------
  USE  kinetic_extrap_boundary, 	 ONLY: read_param_kin_extrap_boundary
  USE  high_resolution_kinetic_fluxes,   ONLY: read_param_hr_kin_fluxes
  USE  kinetic_equations,                ONLY: nsBreak, alpha, phiFmla

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: flow_rgme, idf
  !----------------------------------------------------------------------------------

  READ(idf,*) num_flux
  
  IF (flow_rgme == 2) THEN
    READ(idf,*) nsBreak
    READ(idf,*) phiFmla
    READ(idf,*) alpha
  ENDIF  


  SELECT CASE (num_flux)

    ! High order spatial and temporal accuracy parameteres
    CASE (KFVS_HR, BGK_HR)
    CALL read_param_hr_kin_fluxes(idf)
    CALL read_param_kin_extrap_boundary(idf)

  END SELECT

  END SUBROUTINE  read_param_kinetic_num_fluxes





  SUBROUTINE  init_kinetic_num_fluxes(flow_rgme, ww)
  !----------------------------------------------------------------------
  USE mp_interface
  USE kinetic_equations,	 ONLY: init_kinetic_equations
  USE kinetic_extrap_boundary,   ONLY: init_kin_extrap_boundary,  &
  				       ww_kin_ext_size
  IMPLICIT NONE
  
  INTEGER,                      INTENT(IN) :: flow_rgme
  REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
  !----------------------------------------------------------------------

  ! Finite Volumes metric
  j_c_d  => j_c_d_fvm
  cs_c_d => cs_c_d_fvm
  eta	 => eta_fvm
  Drr	 => Drr_fvm

  CALL init_kinetic_equations(ww, flow_rgme)

  varNumFlux = num_flux
  
  ! Init extrapolation at boundaries
  SELECT CASE (num_flux)
  
    CASE (KFVS_HR, BGK_HR)

      IF (.NOT. MP_job  .OR.  NjP_b /= 0) THEN
        CALL init_kin_extrap_boundary(ww, j_c_d, cs_c_d, bound_p, Drr)
        ext_size = ww_kin_ext_size()
      ELSE
        ext_size = 0
      ENDIF

      ALLOCATE (phie(SIZE(ww,1), SIZE(ww,2) + ext_size),   &
  		 wwe(SIZE(ww,1), SIZE(ww,2) + ext_size))
  END SELECT

  END SUBROUTINE  init_kinetic_num_fluxes





  SUBROUTINE  switch_first_order
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  !--------------------------------------------------------------------------

  IF (num_flux > 2)  varNumFlux = 2

  END SUBROUTINE  switch_first_order





  SUBROUTINE  switch_second_order
  !-------------------------------------------------------------------------				       
  IMPLICIT NONE
  !-------------------------------------------------------------------------

  IF (num_flux > 2)  varNumFlux = num_flux

  END SUBROUTINE  switch_second_order





  SUBROUTINE write_param_kinetic_num_fluxes(idf)
  !---------------------------------------------------------------------
  USE  high_resolution_kinetic_fluxes,   ONLY: write_param_hr_fluxes
    
  IMPLICIT NONE
  INTEGER :: idf
  !---------------------------------------------------------------------

  WRITE(idf,*) '   PARAMETERS FOR KINETIC NUMERICAL FLUXES'
  WRITE(idf,*) '   Numerical flux type: 	', ft_names(num_flux)
  WRITE(idf,*)
  
  ! Write parameters for numerical fluxes
  SELECT CASE (num_flux)
  
    CASE (KFVS_HR, BGK_HR)
    CALL write_param_hr_fluxes(idf)
      
  END SELECT

  END SUBROUTINE  write_param_kinetic_num_fluxes





  FUNCTION  kinetic_domain_flux(flow_rgme, ww, dt)   RESULT (phi)
  !---------------------------------------------------------------------------
  USE structures  
  USE first_order_kinetic_fluxes,	ONLY: KFVS_FO_flux,  &
  					      BGK_FO_flux
  USE high_resolution_kinetic_fluxes,   ONLY: BGK_HR_flux
  USE mp_interface					      
  USE kinetic_extrap_boundary,          ONLY: ww_kin_extrap
  !USE kinetic_boundary_cond, 	        ONLY: ww_kin_extrap
  
  IMPLICIT NONE

  INTEGER,                        INTENT(IN) :: flow_rgme
  REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: ww
  REAL(KIND=8),   DIMENSION(:),   INTENT(IN) :: dt

  REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(ww,2)) :: phi

  INTEGER :: i
  !---------------------------------------------------------------------------

  ! Extrapolation at boundaries      
  IF (ext_size > 0  .AND.  varNumFlux == num_flux) THEN 
    phie = 0 
    wwe(:,1:SIZE(ww,2)) = ww
    wwe(:,SIZE(ww,2)+1:SIZE(wwe,2)) = ww_kin_extrap(ww, j_c_d, Drr)
  ELSE
    phie = 0
    wwe = ww  
  ENDIF


  SELECT CASE (varNumFlux)

    CASE (KFVS_FO)
    phi = KFVS_FO_flux(ww, j_c_d, eta, dt)

    CASE (BGK_FO)
    phi = BGK_FO_flux(ww, j_c_d, eta, dt)

    !CASE (KFVS_HR)
    !phie = KFVS_HR_flux(wwe, j_c_d, eta, Drr, dt)

    CASE (BGK_HR)
    phie = BGK_HR_flux(flow_rgme, wwe, j_c_d, eta, Drr, dt)

  END SELECT


  IF (ext_size > 0  .AND.  varNumFlux == num_flux) THEN

    DO i = 1, SIZE(ww,2)
      phi(:,i) = phie(:,i) 
    ENDDO

  ELSEIF (MP_job  .AND.  MP_worker  .AND.  &
          ext_size == 0  .AND.  varNumFlux == num_flux) THEN

    DO i = 1, SIZE(ww,2)
      phi(:,i) = phie(:,i) 
    ENDDO

  ENDIF


  END FUNCTION  kinetic_domain_flux





  FUNCTION  kinetic_compute_dt(ww, CFL, variable_dt)  RESULT(dt)
  
  ! Actually time step is computed considering only the 
  ! spectral radius of Jacobian of inviscid equations
  !-----------------------------------------------------------------------
  USE thermodynamics,	      ONLY: c__P_r
  USE kinetic_equations,      ONLY: compute_scalar_moments, P__ww, &
                                    IP__c_PSI_g, IN__c_PSI_g, L__ww
  USE messages,               ONLY: terminate
  USE mp_interface			    
  
  IMPLICIT NONE

  REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: ww
  REAL(KIND=8), 		  INTENT(IN) :: CFL 
  LOGICAL,			  INTENT(IN) :: variable_dt  

  REAL(KIND=8), DIMENSION(SIZE(ww,2))  :: dt 
  REAL(KIND=8), DIMENSION(MP_wProc)    :: MP_minDt
  REAL(KIND=8), DIMENSION(SIZE(ww,1))  :: ww_ij
  REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ver, un

  REAL(KIND=8) :: mod_eta, rho_ij, P_ij, c_ij, mod_un, minDt

  INTEGER :: c, i, j, sD
  !----------------------------------------------------------------------- 

  !			   2 C_i
  !  dt_i = CFL -----_------------------,   for each i in K_i
  !		   _/ (|u.n| + c) dCi 
  !               dCi
  !
  ! see: Energy stability analysis of multistep methods
  !      on unstructured meshes. M. Giles

  dt = 0.d0

  sD = SIZE(ww_ij) - 2

  DO c = 1, SIZE(j_c_d,2) 

    i = j_c_d(1,c);  j = j_c_d(2,c)


    ! Error check
    IF (L__ww(ww(:,i)) <= 0.d0) THEN
      OPEN(1000, FILE='npCode.err')
      WRITE(1000,*) 'node pair', c, ' node', i
      WRITE(1000,*) 'ww', ww(:,i)
      CLOSE(1000)
      CALL terminate(6,'KINETIC_COMPUTE_DT','Negative temperature')
    ELSEIF (L__ww(ww(:,j)) <= 0.d0) THEN
      OPEN(1000, FILE='npCode.err')
      WRITE(1000,*) 'node pair', c, ' node', j
      WRITE(1000,*) 'ww', ww(:,j)
      CLOSE(1000)
      CALL terminate(6,'KINETIC_COMPUTE_DT','Negative temperature')
    ENDIF  


    CALL compute_scalar_moments(ww(:,i), 'L')
    CALL compute_scalar_moments(ww(:,j), 'R')
      
    ww_ij = ww(1,i) * IP__c_PSI_g((/0,0,0/),'L')  &
          + ww(1,j) * IN__c_PSI_g((/0,0,0/),'R')

    rho_ij = ww_ij(1)

    mod_eta = SQRT(SUM(eta(:,c)*eta(:,c)))    
    eta_ver = eta(:,c) / mod_eta

    un = SUM((ww_ij(2:sD+1)/ww_ij(1)) * eta_ver)
    mod_un = SQRT(SUM(un*un))

    P_ij = P__ww(ww_ij)    
    c_ij = c__P_r(P_ij, rho_ij)

    dt(i) = dt(i) + (mod_un + c_ij)*mod_eta
    dt(j) = dt(j) + (mod_un + c_ij)*mod_eta

  ENDDO


  ! Local time step  (The sum of 1.d-8
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

  END FUNCTION  kinetic_compute_dt





  FUNCTION  kinetic_compute_CFL(ww, dt)  RESULT(iCFL)  
  !-----------------------------------------------------------------------
  USE thermodynamics,	      ONLY: c__P_r
  USE kinetic_equations,      ONLY: compute_scalar_moments, P__ww, &
                                    IP__c_PSI_g, IN__c_PSI_g, L__ww
  USE messages,               ONLY: terminate				    
  
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
  REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: dt

  REAL(KIND=8) :: iCFL

  REAL(KIND=8), DIMENSION(SIZE(ww,2))  :: vCFL
  REAL(KIND=8), DIMENSION(SIZE(ww,1))  :: ww_ij
  REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ver

  REAL(KIND=8) :: mod_eta, rho_ij, P_ij, c_ij, mod_un, un

  INTEGER :: c, i, j, sD
  !----------------------------------------------------------------------- 

  sD = SIZE(ww_ij) - 2

  vCFL = 0.d0

  DO c = 1, SIZE(j_c_d,2) 

    i = j_c_d(1,c);  j = j_c_d(2,c)


    ! Error check
    IF (L__ww(ww(:,i)) <= 0.d0) THEN
      OPEN(1000, FILE='npCode.err')
      WRITE(1000,*) 'node pair', c, ' node', i
      WRITE(1000,*) 'ww', ww(:,i)
      CALL terminate(6,'KINETIC_COMPUTE_CFL','Negative temperature')
    ELSEIF (L__ww(ww(:,j)) <= 0.d0) THEN
      OPEN(1000, FILE='npCode.err')
      WRITE(1000,*) 'node pair', c, ' node', j
      WRITE(1000,*) 'ww', ww(:,j)
      CALL terminate(6,'KINETIC_COMPUTE_CFL','Negative temperature')
    ENDIF  


    mod_eta = SQRT(SUM(eta(:,c)*eta(:,c)))
    eta_ver = eta(:,c) / mod_eta

    CALL compute_scalar_moments(ww(:,i), 'L')
    CALL compute_scalar_moments(ww(:,j), 'R')
      
    ww_ij = ww(1,i) * IP__c_PSI_g((/0,0,0/),'L') + ww(1,j) * IN__c_PSI_g((/0,0,0/),'R')
      
    rho_ij = ww_ij(1)

    un = SUM((ww_ij(2:sD+1)/ww_ij(1)) * eta_ver)
    mod_un = ABS(un)

    P_ij = P__ww(ww_ij)
    c_ij = c__P_r(P_ij, rho_ij)
    
    vCFL(i) = vCFL(i) + (mod_un + c_ij)*mod_eta
    vCFL(j) = vCFL(j) + (mod_un + c_ij)*mod_eta
    
  ENDDO

  vCFL = vCFL * dt / (2*cell)

  iCFL = MAXVAL(vCFL)

  END FUNCTION  kinetic_compute_CFL


END MODULE  kinetic_num_fluxes
