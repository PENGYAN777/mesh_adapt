!=============================================================================== 
!
!      Module: Turbulence Models
!
! Description: Computation of the turbulence models.
!
!      Author: Alberto Guardone, Stefano Cinquina
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
!===============================================================================
   
  MODULE  turbulence_models
   
  !-----------------------------------------------------------------------------
  USE commons
  USE np_quadrature 
  !USE gas_properties,       ONLY : gas
  !USE transport_properties, ONLY : transport_coefficients
  !USE nodes,                ONLY : jd_jb, rr, bound_p, Np_d, Np_b  
  !USE node_pair_structure,  ONLY : j_c_d_fem  => j_c_d, j_c_b
  !USE metric_coefficients,  ONLY : cell, eta_fem => eta, chi_b, xi_bp, stiff
  !USE euler_equations,      ONLY : T__ww
  !-----------------------------------------------------------------------------

  IMPLICIT NONE;  PRIVATE

  TYPE distance_data
    ! Number of node in domain numeration
    INTEGER :: domain_node
    ! Distance from the nearest wall 
    REAL(KIND=8) :: distance  
  END TYPE distance_data


  !REAL(KIND=8),        DIMENSION(:),   ALLOCATABLE ::  dist_wall
  !REAL(KIND=8),        DIMENSION(:,:), ALLOCATABLE ::  dist_trip
   
  !TYPE(distance_data), DIMENSION(:),   ALLOCATABLE ::  wall_distance

  LOGICAL :: trip
     
  !REAL(KIND=8) :: Re
    
  ! Parameter to impose the trip point on both 
  ! the surface of the profile
  INTEGER, PARAMETER  :: UP   = 1,  &
                         LOW  = 2

  ! Grid spacing along the wall at the trip, along x and y, weights,
  ! % of the chord for the trip point, chord supposed = 1
  !REAL(KIND=8), DIMENSION(2) :: delta_wall, delta_x, delta_y, weight_l
  REAL(KIND=8) :: trip_abscissa

  ! Coordinates of the trip point
  REAL(KIND=8), PARAMETER     :: Prandtl_t = 0.9d0   
  !REAL(KIND=8), DIMENSION(2,2):: trip_point
  !INTEGER,      DIMENSION(2)  :: trip_point_l, trip_point_r
   
  ! Parameter of the Spalart-Allmaras model
  REAL(KIND=8), PARAMETER :: k_kar = 0.41d0
  REAL(KIND=8), PARAMETER :: sigma = (2.d0/3.d0)
  REAL(KIND=8), PARAMETER :: c_v1 = 7.1d0, c_v2 = 5.d0
  REAL(KIND=8), PARAMETER :: c_b1 = 0.1355d0, c_b2 = 0.622d0
  REAL(KIND=8), PARAMETER :: c_w1 = (c_b1/k_kar**2) + ((1+c_b2)/sigma)
  REAL(KIND=8), PARAMETER :: c_w2 = 0.3d0, c_w3 = 2.d0
  REAL(KIND=8), PARAMETER :: c_t1 = 1.0d0, c_t2 = 2.0d0
  REAL(KIND=8), PARAMETER :: c_t3 = 1.2d0, c_t4 = 0.5d0

  INTEGER :: turboModel

  INTEGER, PARAMETER :: SPALART_ALLMARAS = 1
  CHARACTER(*), DIMENSION(1:1), PARAMETER :: tm_names = (/'SPALART-ALLMARAS'/)

  !-----------------------------------------------------------------------------
  PUBLIC :: turbulent_flow, read_param_turbulence!, &
  	    !write_param_ns_turbulent,	      &
            !ns_turbulent_bc, compute_rhs_SA, &
            !compute_mu_t, init_ns_turbulent, &
            !compute_kk_t, write_resid_tur,   &
            !init_mu_tilde, compute_resid_tur
  !-----------------------------------------------------------------------------

  CONTAINS


  SUBROUTINE  read_param_turbulence(idf)
  !------------------------------------------------------------ 
  IMPLICIT NONE
  INTEGER, INTENT(IN)  ::  idf
  !------------------------------------------------------------ 

  READ(idf,*) turboModel    
  READ(idf,*) trip
  
  IF (trip) THEN
    READ(idf,*) trip_abscissa
    trip_abscissa = trip_abscissa / 100.0
  ENDIF

  END SUBROUTINE  read_param_turbulence




!    !============================================================ 
!    SUBROUTINE write_param_ns_turbulent (idf) 
!    !============================================================ 
! 
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       INTEGER,                 INTENT(IN)  ::  idf
!       !------------------------------------------------------------ 
! 
!       WRITE(idf,*) '-----------------------------------------'
!       WRITE(idf,*) ' SPALART-ALLMARAS turbulence model?  ', turbulent_flow
! 
!       IF ( trip ) THEN
!          WRITE(*,*) ' trip point is located at ', (trip_abscissa * 100), &
! 	            '% of the chord'
!       ENDIF
! 
!       WRITE(idf,*) '-----------------------------------------'
!       WRITE(idf,*) ' '
!     
! 
!    END SUBROUTINE write_param_ns_turbulent
!    !============================================================ 
! 
! 
! 
!    !============================================================== 
!    SUBROUTINE init_ns_turbulent (grid_name, name_length) 
!    !============================================================== 
! 
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       ! Strings for the problem name
!       CHARACTER(len=64), INTENT(IN)  ::  grid_name
!       INTEGER,           INTENT(IN)  ::  name_length
! 
!       TYPE boundary_w
!          ! Number of th boundary
!          INTEGER                                     ::  number
!          ! is it a wall? (TRUE/FALSE)
!          LOGICAL                                     ::  t_f
!       END TYPE boundary_w
! 
!       TYPE(boundary_w), DIMENSION(:),   ALLOCATABLE  :: walls
!        
!       INTEGER                        :: i, N_wall
!       REAL(KIND=8)                   :: dummy
!       REAL(KIND=8),     DIMENSION(2) :: d_l, d_r
!       !------------------------------------------------------------ 
! 
! 
!       IF (turbulent_flow) THEN
! 
!          ! read distance from the nearest wall
!          ALLOCATE(wall_distance(Np_d), dist_wall(Np_d))
! 
!          OPEN(5,file='distances.'//grid_name(1:name_length))
! 
! 	 DO i = 1, Np_d
! 	    READ(5,*) wall_distance(i) %domain_node, wall_distance(i) %distance 
! 	 ENDDO
! 	 
! 	 dist_wall = wall_distance %distance
! 	 WHERE(dist_wall <= 1.e-6) dist_wall = 1.e-6
! 
! 	 CLOSE(5)
!  
! 	 ! quantities for the trip term are computed
!          IF ( trip ) THEN
! 
! 	    OPEN(5,file='distances.param')
! 
! 	       READ(5,*)
! 	       READ(5,*) N_wall
! 	       ALLOCATE(walls(N_wall))
! 
! 	       DO i = 1, N_wall
! 
! 		  READ(5,*) walls(i) % number, walls(i) % t_f
! 
! 	       ENDDO
! 
! 	    CLOSE(5) 
! 
! 	    d_l  = - HUGE(d_l); d_r  = HUGE(d_r);
! 
! 	    ! loop to find the nearest points to the trips on both surface
! 	    ! we impose the abscissa of the trip and we look for the ordinate
! 	    ! and the number (domain numeration) of the nearest nodes
! 	    DO i = 1, SIZE(jd_jb)
! 
! 
!   	       IF ( walls (bound_p(i)) % t_f ) THEN
! 
! 
! 		  dummy = rr(1,jd_jb(i)) - trip_abscissa
! 
! 		  IF ( rr(2,jd_jb(i)) > 0 ) THEN ! to consider upper surface
! 
! 		     IF ( dummy >= 0 ) THEN
! 
! 	        	IF ( dummy < d_r(UP) ) THEN
! 	                   d_r(UP) = dummy
! 		           trip_point_r(UP) = jd_jb(i)
!   	        	ENDIF
! 
! 	             ELSE
! 
! 	        	IF ( dummy > d_l(UP) ) THEN
! 	                   d_l(UP) = dummy
! 		           trip_point_l(UP) = jd_jb(i)
! 	        	ENDIF
! 
! 		     ENDIF
! 
! 		  ELSE                           ! to consider lower surface
! 
! 		     IF ( dummy >= 0 ) THEN
! 
! 	        	IF ( dummy < d_r(LOW) ) THEN
! 	                   d_r(LOW) = dummy
! 		           trip_point_r(LOW) = jd_jb(i)
!   	        	ENDIF
! 
! 	             ELSE
! 
! 	        	IF ( dummy > d_l(LOW) ) THEN
! 	                   d_l(LOW) = dummy
! 		           trip_point_l(LOW) = jd_jb(i)
! 	        	ENDIF
! 
! 		     ENDIF
! 
!         	  ENDIF
! 
! 	       ENDIF
! 
! 	    ENDDO
! 
! 	    delta_x    = rr(1,trip_point_r) - rr(1,trip_point_l)
! 	    delta_y    = rr(2,trip_point_r) - rr(2,trip_point_l)
! 	    delta_wall = SQRT (delta_x **2 + delta_y **2)
! 
! 	    weight_l(UP)  = (trip_abscissa - rr(1,trip_point_l(UP)))  / delta_x(UP)
! 	    weight_l(LOW) = (trip_abscissa - rr(1,trip_point_l(LOW))) / delta_x(LOW)
! 
! 	    trip_point(:,1) = trip_abscissa
! 	    trip_point(:,2) = rr(2,trip_point_l) + weight_l * delta_y
! 
! 	    ! distance from trip points is computed
! 	    CALL distance_from_trip_point(trip_point)
!         
!          ENDIF
! 
!       ENDIF
! 
!    END SUBROUTINE init_ns_turbulent
!    !============================================================ 
! 
! 
! 
!    !============================================================== 
!    SUBROUTINE ns_turbulent_ref(Re_ref)
!    !============================================================== 
! 
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8),   INTENT(IN)  ::  Re_ref
!       !------------------------------------------------------------ 
! 
!       IF (turbulent_flow) THEN
!       
!          Re = Re_ref
! 	 
!       ENDIF
! 
!    END SUBROUTINE ns_turbulent_ref
!    !============================================================ 
! 
! 
! 
!    !============================================================== 
!    SUBROUTINE init_mu_tilde(ww, mu_tilde)
!    !============================================================== 
! 
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8),  DIMENSION(:,:), INTENT(IN)  ::  ww
!       REAL(KIND=8),  DIMENSION(:),   INTENT(OUT) ::  mu_tilde
! 
!       !------------------------------------------------------------ 
!       REAL(KIND=8) :: T, mu, dummy
!       INTEGER :: i
!       !------------------------------------------------------------ 
! 
!       IF (turbulent_flow) THEN
!         
! 	 DO i = 1, SIZE(ww,2)
! 
!             T = T__ww(ww(:,i))
!             CALL transport_coefficients(T, 1.d0/ww(1,i), mu, dummy)
! 
!             mu_tilde(i) = mu / 10
! 
! 	 ENDDO	 
! 
!       ELSE
! 	 
! 	 mu_tilde = 0.d0
!      
!       ENDIF
! 
!    END SUBROUTINE init_mu_tilde
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION ns_turbulent_bc( mu_tilde_b, mu, value_type, bpoints ) &
!                                        RESULT( mu_tildeb )
!    !============================================================ 
!      
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8),       DIMENSION(:),    INTENT(IN) ::  mu_tilde_b
!       REAL(KIND=8),       DIMENSION(:),    INTENT(IN) ::  mu
!       INTEGER,            DIMENSION(:),    INTENT(IN) ::  bpoints
!       INTEGER,                             INTENT(IN) ::  value_type
! 
!       REAL(KIND=8),       DIMENSION(SIZE(mu_tilde_b)) ::  mu_tildeb
!       !------------------------------------------------------------ 
!       INTEGER,      PARAMETER  ::  BV__VV_NOSLIP = 1
!       !------------------------------------------------------------ 
! 
!       IF (turbulent_flow) THEN
!   
!             SELECT CASE ( value_type )            
!               
!               !------------------
!               CASE(BV__VV_NOSLIP)
!               !------------------
!               mu_tildeb(bpoints)  = 0.d0 
! 	      
!               !--------------
!               CASE DEFAULT
!               !--------------
!               mu_tildeb(bpoints)  = mu(bpoints) / 10
!            
!             END SELECT
!        
!       ELSE
!       
!       mu_tildeb = 0.d0
!       
!       ENDIF
! 
!    END FUNCTION ns_turbulent_bc
!    !============================================================ 
! 
! 
! 
!    !=========================================================
!    FUNCTION compute_rhs_SA(ni_tilde, ni_tildeb, ni, rho,    &
!                                             vv, vb, R_vv )  &
!                                             RESULT(rhs_t)
!    !=========================================================
!   
!       !--------------------------------------------------------------- 
!       IMPLICIT NONE
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)   ::  ni_tilde
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)   ::  ni_tildeb
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)   ::  ni
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)   ::  rho
!       REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   ::  vv
!       REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   ::  vb
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)   ::  R_vv
! 
!       REAL(KIND=8), DIMENSION(SIZE(ni_tilde))    ::  rhs_t
!       !---------------------------------------------------------------
!       
!       REAL(KIND=8), DIMENSION(SIZE(vv,2))            ::  mod_vv
!       REAL(KIND=8), DIMENSION(SIZE(vv,1),SIZE(vv,2)) ::  cc, G_ni_tilde
!       REAL(KIND=8), DIMENSION(SIZE(vb,1),SIZE(vb,2)) ::  ccb, G_ni_tildeb
!       ! auxiliar function for the model
!       REAL(KIND=8), DIMENSION(SIZE(ni_tilde))        ::  f_v1, f_v2, f_v3,  &
!       						         f_t1, f_t2, f_w,   &
! 						         source,            &
! 							 r, g, e, q,        &
! 							 chi, S_tilde
!       REAL(KIND=8), DIMENSION(2,SIZE(ni_tilde))      ::  power, dummy, g_t
!       ! RHS and numerical fluxes						 
!       REAL(KIND=8), DIMENSION(SIZE(ni_tilde))        ::  phi_c, phi_d, phi_s,  &
! 							 phi_p, phi_trip
!       ! Vorticity at trip points
!       REAL(KIND=8), DIMENSION(2)                     ::  omega_t
!       REAL(KIND=8)    ::   ss
!       ! Indices
!       INTEGER         ::   c, i_, j_, i, j, k
!       !---------------------------------------------------------------
! 			   
! 
!       rhs_t = 0.d0
!   
!       !========================
!       IF ( turbulent_flow ) THEN
!       !========================
!       
!          ! auxiliar functions
!          chi     = ni_tilde / ni
! 	 WHERE( chi <= 0.001 ) chi = 0.001
! 	 
!          f_v1 = chi**3 / (chi**3 + c_v1**3)
! 	 
!          f_v2 = ( 1 + chi/c_v2 )**(-3)
! 	 
! 	 f_v3 = ( 1 + chi*f_v1) * ( 1 - f_v2 ) / chi
! 
! 	 S_tilde =    R_vv * f_v3 &
! 	            + ( ni_tilde*f_v2 ) / ( k_kar**2 * dist_wall**2 * Re )
!          WHERE( S_tilde <= 1.e-10 ) S_tilde = 1.e-10
! 
!          r =  ni_tilde /  ( S_tilde * k_kar**2 * dist_wall**2 * Re )
!          
!          g =  r + c_w2 * ( r**6 - r )
! 
!          f_w =  g * (  ( 1    + c_w3**6 ) / &
!                        ( g**6 + c_w3**6 )   ) **(1.d0/6.d0)
! 
! ! old        f_v2    = 1 - chi / ( 1+c_v1*f_v1 )
! ! version    S_tilde  = R_vv + (ni_tilde/Re) * &
! ! of the SA  f_v2 / (k_kar**2 * dist_wall**2)
! ! model      WHERE( S_tilde <= 0 ) S_tilde = 0
! 
!          ! nodal gradient of ni_tilde in weak form
!          G_ni_tilde = np_quadr_w_Gu (j_c_d_fem, j_c_b, jd_jb,  &
!    		        	     eta_fem, chi_b, xi_bp,    &
! 				     ni_tilde, ni_tildeb )
! 
!          DO i = 1, SIZE(vv,2)
!              G_ni_tilde (:,i) = G_ni_tilde (:,i) / cell(i)
!          ENDDO
! 
!          G_ni_tildeb(:,:) = G_ni_tilde(:,jd_jb)
!          
! 	 
! 	 ! functions for transition
! 	 !---------------------------------------------------
! 
!          f_t2 = 0.d0; phi_trip = 0.d0
! 
!          IF ( trip ) THEN
! 
!             ! module of the velocity 
!             DO i = 1, SIZE(vv,2)
!         	mod_vv(i) = SQRT( SUM( vv(:,i)**2 ) )
!             ENDDO
! 
! 	    WHERE( mod_vv <= 1.e-10) mod_vv = 1.e-10  ! to avoid blow-up, 
! 	    ! anyway Dirichlet condition is imposed on wall nodes
! 
!             ! vorticity at the trip points, linear interpolation
!             omega_t =    R_vv(trip_point_l)  &
!                        + weight_l * ( R_vv(trip_point_r) - R_vv(trip_point_l) )
! 
!             f_t2 = c_t3 * exp ( - c_t4 * chi**2 ) 
! 
!             g_t(UP, :) = mod_vv / ( omega_t(UP)  * delta_wall(UP) )
! 	    g_t(LOW,:) = mod_vv / ( omega_t(LOW) * delta_wall(LOW) )
! 
!             WHERE (g_t >= 0.1) g_t = 0.1
! 
!             power(UP, :) =  c_t2 * omega_t(UP)**2   &
!                             * ( dist_wall**2 + (dist_trip(UP, :)*g_t(UP, :))**2 ) / mod_vv**2
!             power(LOW,:) =  c_t2 * omega_t(LOW)**2  &
!                             * ( dist_wall**2 + (dist_trip(LOW,:)*g_t(LOW,:))**2 ) / mod_vv**2
! 
!             dummy(UP, :) = g_t(UP, :) * exp( - power(UP, :))
!             dummy(LOW,:) = g_t(LOW,:) * exp( - power(LOW,:))
! 
!             f_t1 = c_t1 * MAX(dummy(UP,:), dummy(LOW,:)) 
! 
!          ENDIF    
!     
!          !==================================================!
! 	 ! Computation of the RHS of the transport equation !
!          !==================================================!
! 
! 
!          ! Advective Term      !
!          !====================!
! 	 
!          ! - DIV( mu_tilde * vv )
! 	 
! 	 DO k = 1, SIZE(vv,1)
!             cc (k,:) = ni_tilde * vv(k,:) * rho
!             ccb(k,:) = ni_tilde(jd_jb) * vb(k,:) * rho(jd_jb)
!          ENDDO
!         
!      
!          phi_c = np_quadr_w_Dv_uw( j_c_d_fem, j_c_b, jd_jb, &
!                                    eta_fem, chi_b, xi_bp,   &
!                                    (ni_tilde*rho), cc, ccb )
! 
!  
!          ! Diffusive Term    !
!          !===================!
! 	 
!          ! DIV[ ( mu + mu_tilde ) G_ni_tilde ]
! 	 
! 	 e = (ni + ni_tilde) * rho
! 
!          ! Domain contribution 
!          phi_d = - np_quadr_Gw_nGu( j_c_d_fem, stiff , &
! 	                            ni_tilde, e )
! 
!          ! Boundary contribution      
!             ! Boundary node-pairs contribution
!          DO c = 1, SIZE( chi_b,2 ) 
!             i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
!             i  = jd_jb(i_);   j  = jd_jb(j_)
!          
!             ss =   SUM( chi_b(:,c) * ( e(j) * G_ni_tildeb(:,j_) &
! 	                             - e(i) * G_ni_tildeb(:,i_) ) )
!          
!             phi_d(i) = phi_d(i) + ss
!             phi_d(j) = phi_d(j) - ss
! 
!          ENDDO
!       
!          ! Boundary nodes contributions
!          DO i_ = 1, SIZE( xi_bp,2 );  i = jd_jb(i_)
!          
!             ss =  SUM( xi_bp(:,i_) * e(i) * G_ni_tildeb(:,i_) )
!          
! 	    phi_d(i) = phi_d(i) + ss
! 
!          ENDDO
!   
! 
!          ! Production Term    !
! 	 !====================!
! 
!          ! G_ni_tilde * G_ni_tilde
! 	 			  
!          DO i = 1, SIZE(vv,2)
!              q(i) = SUM( G_ni_tilde(:,i)**2 )
!          ENDDO
! 	      
!          phi_p = c_b2 * rho * cell * q
! 
! 
!          ! Source Term      !
! 	 !==================!
! 	 
!          source =   c_b1 * S_tilde * ni_tilde * ( 1 - f_t2 )    &
! 	         -  ( c_w1 * f_w - c_b1 * f_t2 / k_kar**2 )     &
! 		 * (ni_tilde/dist_wall)**2 / Re
! 		 
!          phi_s = cell * source * rho
! 
! 
!          ! Trip Term        !
! 	 !==================!
! 
!          IF ( trip ) phi_trip = Re * rho * f_t1 * mod_vv**2 * cell
! 
! 
! 	 ! RHS !
! 	 !=====!
! 
! 	 rhs_t = - phi_c + (phi_d  +  phi_p) / (sigma*Re) + phi_s + phi_trip
! 	       
!       ENDIF
!       
!       
!    END FUNCTION compute_rhs_SA
!    !==========================================
! 
! 
! 
!    !==========================================
!    FUNCTION compute_mu_t( mu_tilde, mu ) &
!                                  RESULT(mu_t)
!    !==========================================
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  mu_tilde
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  mu
!             
!       REAL(KIND=8), DIMENSION(SIZE(mu_tilde))   ::  mu_t
!       !------------------------------------------------------------ 
!       
!       REAL(KIND=8), DIMENSION(SIZE(mu_tilde))   ::  f_v1
!       REAL(KIND=8), DIMENSION(SIZE(mu_tilde))   ::  chi
!       !------------------------------------------------------------ 
!   
!       IF(turbulent_flow) THEN
!          
!          ! we don't need to divide by the reference viscosity
!          ! becaouse all the quantities are dimensionless,
!          ! chi is a number and mu_tilde is dimensionless
! 
!          chi   = mu_tilde / mu
!          f_v1 = chi**3 / (chi**3 + c_v1**3)
! 
!          mu_t  = mu_tilde * f_v1
! 
!       ELSE
!          
!  	 mu_t = 0.d0
! 	
!       ENDIF 
! 
! 
!    END FUNCTION compute_mu_t
!    !==========================================
! 
!    !==========================================
!    FUNCTION compute_kk_t( mu_t ) &
!                                  RESULT(kk_t)
!    !==========================================
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8), DIMENSION(:),   INTENT(IN)     ::  mu_t
! 
!       REAL(KIND=8), DIMENSION(SIZE(mu_t))          ::  kk_t
!       !------------------------------------------------------------
! 
!       REAL(KIND=8)    :: Cp
!       !------------------------------------------------------------ 
!   
!       IF(turbulent_flow) THEN
!       
! 	 ! in MODULE transport_properties FUNCTION transport_coefficients
! 	 ! we have 
!          !
! 	 !    kk = (mu * Rgas*cp__T_v(T_adim,v_adim)) / Pr
! 	 !    kk_adim = kk / kk_ref
! 	 !
!          ! where the mu is the dimensional one, kk_adim is the coefficient
! 	 ! we put in NS equations, so we look 
! 	 ! for the dimensional turbulent viscosity (mu_t * gas%mu_ref)
! 	 ! and after we divide by (gas%kk_ref) to have the kk_adim.
! 	 ! => the coefficient we put in NS equations is adimensional
! 
! 	 Cp     = (7.d0/2.d0) * gas%Rgas
! 
! 	 kk_t   = ( mu_t * gas%mu_ref ) * Cp / Prandtl_t
! 
! 	 kk_t   = kk_t / gas%kk_ref
! 
!       ELSE
!          
!  	 kk_t = 0.d0
! 	
!       ENDIF 
! 
! 
!    END FUNCTION compute_kk_t
!    !==========================================
! 
! 
! 
!    !============================================================== 
!    SUBROUTINE distance_from_trip_point(trip_point)
!    !============================================================== 
! 
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8), DIMENSION(:,:),  INTENT(IN)  ::  trip_point
!       
!       INTEGER     :: i
!       !------------------------------------------------------------ 
! 
!       ALLOCATE(dist_trip(2,Np_d))
!       
!       DO i = 1, Np_d
!          
!          dist_trip(UP,i)  = SQRT((rr(1,i) - trip_point(UP,1))**2  + &
!                                  (rr(2,i) - trip_point(UP,2))**2    )
!          dist_trip(LOW,i) = SQRT((rr(1,i) - trip_point(LOW,1))**2 + &
!                                  (rr(2,i) - trip_point(LOW,2))**2   )
! 			   
!       ENDDO
! 
!    END SUBROUTINE distance_from_trip_point
!    !============================================================= 
! 
! 
! 
!    !============================================================== 
!    SUBROUTINE write_resid_tur(mu_tilde_resid)
!    !============================================================== 
! 
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8),  INTENT(IN)  :: mu_tilde_resid
!       !------------------------------------------------------------ 
! 
!       
!       IF ( turbulent_flow ) THEN
!          
!          WRITE(*,*) '   tur_resid =', mu_tilde_resid
! 	 
!       ENDIF
! 
! 
!    END SUBROUTINE write_resid_tur
!    !============================================================= 
! 
! 
! 
!    !============================================================== 
!    FUNCTION compute_resid_tur(mu_tilde, mu_tilde_old, dt, cell)  &
!                                             RESULT (mu_tilde_resid)
!    !============================================================== 
! 
! 
!       !------------------------------------------------------------ 
!       IMPLICIT NONE
! 
!       REAL(KIND=8),  DIMENSION(:),  INTENT(IN)   :: mu_tilde
!       REAL(KIND=8),  DIMENSION(:),  INTENT(IN)   :: mu_tilde_old
!       REAL(KIND=8),  DIMENSION(:),  INTENT(IN)   :: dt
!       REAL(KIND=8),  DIMENSION(:),  INTENT(IN)   :: cell
! 
!       REAL(KIND=8)   :: mu_tilde_resid
!       !------------------------------------------------------------ 
! 
!       
!       mu_tilde_resid = 1.d0
! 
!       IF ( turbulent_flow ) THEN
! 
!          mu_tilde_resid = SQRT( SUM( (((mu_tilde - mu_tilde_old)/dt)**2) * cell ))
! 	 
!       ENDIF
!       
! 
!    END FUNCTION compute_resid_tur
!    !============================================================= 

  END MODULE  turbulence_models
