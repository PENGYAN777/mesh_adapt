!=============================================================================== 
!
!      Module: kinetic_boundary_cond
!
! Description: Enforces boundary conditions according to the
!              kinetic formulation
!
!      Author: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!===============================================================================

    MODULE  kinetic_boundary_cond

    USE commons

    IMPLICIT NONE

    !-----------------------------------------------------------------------------
    TYPE boundary_data
  
       ! Type of boundary condition to be imposed
       INTEGER  			      ::  cond_type
       ! Type of boundary value 
       INTEGER  			      ::  value_type
       ! Size of and vector of boundary data 
       INTEGER  			      ::  size_bdata
       REAL(KIND=8), DIMENSION(:,:), POINTER  ::  bdata
       ! Number of and list of boundary nodes 
       ! (boundary indices)
       INTEGER  			      ::  Np_bound    
       INTEGER,      DIMENSION(:),   POINTER  ::  bpoints
       ! Vector of imposed boundary values 
       ! (conservative variables)
       REAL(KIND=8), DIMENSION(:,:), POINTER  ::  bb
  
    END TYPE boundary_data


    ! Boundary condition table
    TYPE(boundary_data), DIMENSION(:), ALLOCATABLE  ::  bc_table
    ! List of boundary points belonging to solid bodies
    INTEGER, DIMENSION(:), ALLOCATABLE  ::  solid_wall_points

    ! Constants: 
    !	Boundary value types						 
    INTEGER, PARAMETER  ::   BV__NONE		 = 0, &
    			     BV__WW_OO  	 = 1, &
    			     BV__WW_0		 = 3, &
    			     BV__CONST_WW	 = 4, &
    			     BV__CONST_P	 = 5, &
    			     BV__CONST_RHO_P	 = 6, &
    			     BV__CONST_RHO_U_P   = 7

    CHARACTER(*), DIMENSION(0:7), PARAMETER  ::  bv_names = &
    			 (/ 'NO BOUNDARY VALUE IMPOSED              ', &
    			    'WW_OO (CONS. VARIABLE AT INFINITY)     ', &
    			    'not used for kinetic                   ', &
                            'WW_0 (INITIAL VALUE)                   ', &
    			    'CONSTANT CONSERVATIVE VARIABLE WW_C    ', &
    			    'CONSTANT PRESSURE                      ', &
    			    'CONSTANT DENSITY AND PRESSURE          ', &
    			    'CONSTANT DENSITY, VELOCITY AND PRESSURE' /)

    !	Boundary conditions types					      
    INTEGER, PARAMETER  ::   BC__NONE		 = 0, &
    			     BC__SOLID_WALL	 = 1, &
    			     BC__RIEMANN	 = 2, &
    			     BC__SIMM_PLANE	 = 3

    CHARACTER(*), DIMENSION(0:3), PARAMETER  ::  bc_names = &
    			 (/ 'NO BOUNDARY CONDITION',  &
    			    'SLIP CONDITION       ',  &
    			    'RIEMANN              ',  &
    			    'SIMMETRY PLANE       '  /)

    INTEGER, PARAMETER :: EQUILIBRIUM = 0, &
    			  QUASI_EQUIL = 1, &
    			  TRANS_REGIM = 2
    			     
    INTEGER :: equilType

    REAL(KIND=8) :: zeta, T_wall, Rgas

    PRIVATE  ::  bc_table, solid_wall_points,			&
    		 BV__NONE,		 BV__WW_OO,		&
    		 BV__CONST_WW,  	 BV__CONST_P,		&
    		 BV__CONST_RHO_P,	 BV__CONST_RHO_U_P,	&
    		 BC__SOLID_WALL,	 BC__RIEMANN,		&
    		 BC__SIMM_PLANE
    		 !BV__WW_OO_FFCORRECT,	 BV__WW_0,		&
		 !BC__NON_REFLECT  

    PUBLIC   ::  kinetic_boundary_flux, 			&
    		 init_kinetic_boundary_cond,			&
    		 read_param_kin_boundary_cond,  		&
    		 write_param_kin_boundary_cond, 		&
    		 kinetic_type_of_boundary
   		 !wb__kinetic_boundary_cond,			&
                 !ww_b__bc, S_kin_impose_strong,                &
  		 !kbody_pnts, kN_body_pnts, kCL_body, kCD_body, &
                 !ww_kin_extrap

    !PRIVATE  ::  bb__boundary_data,    impose_kinetic_bc,	&
    		 !impose_wall_slip_bc,  impose_free_stream_bc,	&
    		 !impose_simm_plane_bc
    !-----------------------------------------------------------------------------

    CONTAINS


    SUBROUTINE  kinetic_boundary_flux(equil_type, ww, GG_ww_b, HH_ww_b, normal, dt_, rhs, v_slip)
    !----------------------------------------------------------------------------------------------
    USE kinetic_equations
    USE messages,               ONLY: terminate
    USE euler_equations,        ONLY: flux__ww
    USE nodes,                  ONLY: jd_jb
    USE metric_coefficients,    ONLY: cosines_XI_BP
    USE axial_symmetry,         ONLY: xi_bp_2d
    USE thermodynamics,         ONLY: T__P_v
    USE transport_properties,   ONLY: transport_coefficients
    !USE nodes,                  ONLY: jd_jb, bound_p, jb_jd, jdm_jb
    
    IMPLICIT NONE
    
    INTEGER,                          INTENT(IN)    :: equil_type
    REAL(KIND=8), DIMENSION(:,:),     INTENT(IN)    :: ww
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)    :: GG_ww_b
    REAL(KIND=8), DIMENSION(:,:),     INTENT(IN)    :: normal
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)    :: HH_ww_b
    REAL(KIND=8), DIMENSION(:),       INTENT(IN)    :: dt_
    REAL(KIND=8), DIMENSION(:,:),     INTENT(INOUT) :: rhs
    
    REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,1)) :: RR, Mb
    
    REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(jd_jb,1)) :: ww_b
    REAL(KIND=8), DIMENSION(SIZE(ww,1),sDim) :: Gw_b
    REAL(KIND=8), DIMENSION(sDim,sDim,sDim+2) :: Hw_b
        
    REAL(KIND=8), DIMENSION(:,:), POINTER :: bdata
    
    REAL(KIND=8), DIMENSION(SIZE(ww,1)*(SIZE(ww,1)+1)/2) :: pM_b, sysMat

    REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(jd_jb,1)) :: wb
    REAL(KIND=8), DIMENSION(:,:), POINTER :: bb
    
    REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: ww_INC, ww_REF, Fb, ww_xx
    REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: Gx, Gy, Gz, Gt, v1

    REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(normal,1)) :: ff_i
    
    REAL(KIND=8), DIMENSION(SIZE(ww,2)) :: v_slip
    
    REAL(KIND=8), DIMENSION(sDim) :: nVer

    REAL(KIND=8), PARAMETER :: PI_G = 3.141592653589793

    REAL(KIND=8) :: rb, r_wall, l_wall, r_INC, r_REF, l_REF, dt, &
                    tau, PHI, iKn

    INTEGER, DIMENSION(:), POINTER :: bpoints
       
    INTEGER, DIMENSION(SIZE(ww_b,1)) :: IPIV
    INTEGER, DIMENSION(3) :: order_0, order_1, order_2
    
    INTEGER :: cond_type, value_type
    INTEGER :: j, j_, b, p, q, pi, info !, jd

     REAL(KIND=8) :: mu_b, kk_b, mfp, Tb, Pb

    ! Gocken model
    !REAL(KIND=8) :: v_n, coeff_A, coeff_S
    !-------------------------------------------------------------------------------------
        
    v_slip = 0.d0    
    
    ! The conserved values at boundary is taken element-wise, in the sense
    ! that instead of the node variables a weighted average of the
    ! values at the two nodes of the boundary element is considered
    ! (the choice of the weigths 1/4 and 3/4 is VALID in two dimensions only)
    !DO p = 1, SIZE(jd_jb)
    !  
    !  jd = jd_jb(p)      
    !  
    !  IF (bound_p(jb_jd(1,jd))  /=  bound_p(jb_jd(2,jd))) THEN
    !    ww_b(:,p) = 3.d0/4.d0 * ww(:,jdm_jb(1,p)) + 1.d0/4.d0 * ww(:,jdm_jb(2,p))      
    !  ELSE
    !    ww_b(:,p) = ww(:,jd_jb(p))
    !  ENDIF 
    !   
    !ENDDO

    ww_b = ww(:,jd_jb)


    DO b = 1, SIZE(bc_table)

       cond_type  =  bc_table(b) % cond_type
       value_type =  bc_table(b) % value_type
       bpoints    => bc_table(b) % bpoints
       bdata	  => bc_table(b) % bdata
       bb	  => bc_table(b) % bb

       ! Conserved variables on the basis of the boundary data imposed 
       bb = imposed_boundary_values(ww_b(:,bpoints), bdata, value_type)

       ! Conserved variables after the imposition of the boundary
       ! conditions. These are used to compute boundary fluxes.
       IF (axisym) THEN
         wb(:,bpoints) = impose_boundary_cond(bb, ww_b(:,bpoints), xi_bp_2d(:,bpoints), &
                                              cond_type, 0, equil_type)
       ELSE
         wb(:,bpoints) = impose_boundary_cond(bb, ww_b(:,bpoints), normal(:,bpoints), &
                                              cond_type, 0, equil_type)
       ENDIF				    

       DO j_ = 1, SIZE(bpoints)
         
	 j = bpoints(j_)

         IF ( cond_type == 2  .or. equil_type == 0 ) THEN
	 
	        ff_i = flux__ww(wb(:,j))
	 	 
	        DO p = 1, SIZE(wb,1)
	          rhs(p,jd_jb(j)) = rhs(p,jd_jb(j)) - SUM(normal(:,j) * ff_i(p,:))
	        ENDDO
	 
	 ELSE
       
                dt = dt_(jd_jb(j))
                nVer = normal(:,j) / SQRT(SUM(normal(:,j)**2))

                ! Rotation matrix to rotate from original frame
                ! into one aligned with normal vector xi_bp
                RR = 0.d0
                RR(1,1) = 1.d0
                RR(SIZE(wb,1),SIZE(wb,1)) = 1.d0

                DO p = 1, sDim
                  RR(p+1,2:1+sDim) = cosines_XI_BP(p,:,j)
                ENDDO

	        ! Conserved variables, gradient and Hessian 
                ! in the local rotated frame
	        wb(:,j) = MATMUL(RR, wb(:,j))

                IF (L__ww(wb(:,j)) <= 0.d0) THEN
                  OPEN(1000, FILE='npCode.err')
                  WRITE(1000,*) 'boundary node', j
                  WRITE(1000,*) 'corresponding domain node', jd_jb(j)
	          WRITE(1000,*) 'boundary state', wb(:,j)
   	          CLOSE(1000)		
	          CALL terminate(6,'KINETIC_BOUNDARY_FLUX:','Negative temperature.')
	        ENDIF  
	 
	        CALL compute_scalar_moments(wb(:,j), 'L')
	 

		! Gradient at the boundary node j ---------------------------------------
                IF (sDim == 1) THEN
         
                  Gw_b = MATMUL(RR,GG_ww_b(:,:,j))
         
                ELSEIF (sDim == 2) THEN
	          
                  Gw_b(1,:)   = MATMUL(GG_ww_b(1,:,j),TRANSPOSE(cosines_XI_BP(:,:,j)))
                  Gw_b(2:3,:) = MATMUL(cosines_XI_BP(:,:,j), GG_ww_b(2:3,:,j))
                  Gw_b(2:3,:) = MATMUL(Gw_b(2:3,:),   TRANSPOSE(cosines_XI_BP(:,:,j)))
                  Gw_b(4,:)   = MATMUL(GG_ww_b(4,:,j),TRANSPOSE(cosines_XI_BP(:,:,j)))

                ENDIF

	        ! Hessian Matrix at the boundary node j ---------------------------------
	        IF (equilType == TRANS_REGIM) THEN
	        
		  DO  p = 1, SIZE(Hw_b,2)
		    
		    Hw_b(:,p,1)   = MATMUL(HH_ww_b(:,p,j,1),   TRANSPOSE(cosines_XI_BP(:,:,j)))
		    Hw_b(:,p,2:3) = MATMUL(cosines_XI_BP(:,:,j), HH_ww_b(:,p,j,2:3))
		    Hw_b(:,p,2:3) = MATMUL(HH_ww_b(:,p,j,2:3), TRANSPOSE(cosines_XI_BP(:,:,j)))
		    Hw_b(:,p,4)   = MATMUL(HH_ww_b(:,p,j,4),   TRANSPOSE(cosines_XI_BP(:,:,j)))

		  ENDDO
		
		ENDIF


                IF (cond_type == BC__SIMM_PLANE) THEN
		  Gw_b(:,1)   = 0.d0
                  Hw_b(1,1,:) = 0.d0
                ENDIF

                order_0 = 0
	        rb = wb(1,j)										
                Mb = rb * I__c_PSIxPSI_g(order_0,'L')							       
                												     
                pi = 0  											     
                DO p = 1, SIZE(ww_b,1)  									     
                  pM_b(pi+1:pi+p) = Mb(1:p,p)									
                  pi = pi + p											     
                ENDDO												     

                ! Space derivatives of distribution function --------------------------------------------
                sysMat = pM_b								  
                Gx = Gw_b(:,1)
                CALL DSPSV('U', SIZE(ww_b,1), 1, sysMat, IPIV, Gx, SIZE(ww_b,1), info)      

                IF (sDim >= 2) THEN							      
                  sysMat = pM_b 							  
                  Gy = Gw_b(:,2)						  
                  CALL DSPSV('U', SIZE(ww_b,1), 1, sysMat, IPIV, Gy, SIZE(ww_b,1), info)    
                ENDIF;IF (sDim == 3) THEN						      
                  sysMat = pM_b 							  
                  Gz = Gw_b(:,3)		  
                  CALL DSPSV('U', SIZE(ww_b,1), 1, sysMat, IPIV, Gz, SIZE(ww_b,1), info)    
                ENDIF 
                									    
                										       
                ! Time derivatives of distribution function. Compatibility on Dg ------------------------
                order_1 = (/1, 0, 0/)							    
                Gt = - rb * (MATMUL(I__c_PSIxPSI_g(order_1,'L'),Gx))			 
                									    
                IF (sDim >= 2) THEN							      
                  order_1 = (/0, 1, 0/) 						    
                  Gt = Gt - rb * (MATMUL(I__c_PSIxPSI_g(order_1,'L'),Gy))		 
                ENDIF;IF (sDim == 3) THEN						      
                  order_1 = (/0, 0, 1/) 						    
                  Gt = Gt - rb * (MATMUL(I__c_PSIxPSI_g(order_1,'L'),Gz))		 
                ENDIF


                sysMat = pM_b
                CALL DSPSV('U', SIZE(ww_b,1), 1, sysMat, IPIV, Gt, SIZE(ww_b,1), info)


	        IF (equil_type == TRANS_REGIM) THEN
	 
                    ww_xx = 0.d0 !Hw_b(1,1,:)
                    		      
                    Pb = P__ww(wb(:,j))
                    Tb = T__P_v(Pb, 1.d0/wb(1,j))
                    		      
                    CALL transport_coefficients(Tb, 1.d0/wb(1,j), mu_b, kk_b)

                    mfp = 16.d0/(5*SQRT(PI_G)) * mu_b / (wb(1,j)*SQRT(2*Pb/wb(1,j)))
                    		      
                    ! Only density is considered here instead of V and T
                    iKn = 0.d0 !mfp * SQRT(SUM(Gw_b(:,1))**2) / wb(1,j)
                    	     
                    	     
                    tau = mean_collision_time(dt, wb(:,j), wb(:,j), wb(:,j), &
                    				 iKN, pM_b, Gx, Gt, ww_xx)
                  
                ELSE

                    tau = mean_collision_time(dt, wb(:,j), wb(:,j), wb(:,j))	 

                ENDIF


                IF (cond_type == BC__SOLID_WALL) THEN

                     ! Computation of the slip velocity as specified by Gokcen
                     ! model. Macroscopic approach used for comparison -----------------------------------------
                     
                     ! Derivative in the normal direction of the tangential 
                     ! component of the velocity
                     !coeff_A = SQRT(2.0/PI_G)
                     !coeff_S = 1.0 !0.858d0
                      
                     !v_n = 1.0/wb(1,j)*(Gw_b(3,1) - wb(3,j)/wb(1,j)*Gw_b(1,1)) 

	             !Pb = P__ww(wb(:,j))
	             !Tb = T__P_v(Pb, 1.d0/wb(1,j))
	          
	             !CALL transport_coefficients(Tb, 1.d0/wb(1,j), mu_b, kk_b)

                     !mfp = 16.d0/(5*SQRT(PI_G)) * mu_b / (wb(1,j)*SQRT(2*Pb/wb(1,j)))  
                     
                     !v_slip(jd_jb(j)) = v_slip(jd_jb(j)) + coeff_A*((2.0-coeff_S)/coeff_S * mfp * v_n) / 2.0
                     ! -----------------------------------------------------------------------------------------
                     
                     
                     ! INCOMING GAS
                     ww_INC = wb(:,j)
                     r_INC  = wb(1,j)

                     order_1 = (/1,0,0/)
                     order_2 = (/2,0,0/)

                     v1 = 0.d0;   v1(1) = 1.d0
                     
                     PHI = 2.0*SQRT(PI_G) * r_INC*(		   DOT_PRODUCT(IP__c_PSI_g(order_1,'L'),v1)  & 
                				  + (0.5*dt-tau) * DOT_PRODUCT(IP__c_PSI_g(order_1,'L'),Gt)  & 
                				  -	    tau  * DOT_PRODUCT(IP__c_PSI_g(order_2,'L'),Gx))
                     IF (sDim >= 2) THEN							   
                       order_2 = (/1, 1, 0/)							 
                       PHI = PHI - 2.0*SQRT(PI_G) * r_INC*tau * DOT_PRODUCT(IP__c_PSI_g(order_2,'L'),Gy)     
                     ENDIF;IF (sDim == 3) THEN  						   
                       order_2 = (/1, 0, 1/)							 
                       PHI = PHI - 2.0*SQRT(PI_G) * r_INC*tau * DOT_PRODUCT(IP__c_PSI_g(order_2,'L'),Gz)     
                     ENDIF
 
 
                     ! REFLECTED GAS 
                     IF (zeta == 0.d0) THEN
                     
	               l_REF = PHI**2 / (2*(gamma-1)*ww_INC(SIZE(ww_INC)))**2
                       r_REF = SQRT(l_REF) * PHI
                     
                     ELSEIF (zeta == 1.d0) THEN

                       l_wall = 1.d0 / (2.0*Rgas*T_wall)
                       r_wall = SQRT(l_wall) * PHI
                       
                       l_REF = l_wall
                       r_REF = r_wall
                     
                     ELSE
                     
                       l_wall = 1.d0 / (2.0*Rgas*T_wall)	      
                       r_wall = SQRT(l_wall) * PHI	    
	               
                       l_REF = PHI**2 / (2*(1-zeta)*(gamma-1)*ww_INC(SIZE(ww_INC)) + zeta*r_wall/l_wall)**2
                       r_REF = SQRT(l_REF) * PHI

                     ENDIF


                     ww_REF(1) = r_REF
	             ww_REF(2:sDim+1) = 0.d0
                     ww_REF(SIZE(ww_REF)) = r_REF / (2.0*(gamma-1)*l_REF)

                     CALL compute_scalar_moments(ww_REF, 'R')


                     ! Wall Boundary Flux --------------------------------------------------------------------------------
	             order_1 = (/1,0,0/)
                     order_2 = (/2,0,0/)
                     
                     Fb = r_INC*(IP__c_PSI_g(order_1,'L') + (0.5*dt-tau)*MATMUL(IP__c_PSIxPSI_g(order_1,'L'),Gt)  &
	        	
	        	- tau*MATMUL(IP__c_PSIxPSI_g(order_2,'L'),Gx)) + r_REF*IN__c_PSI_g(order_1,'R')
	
                     IF (sDim >= 2) THEN
	               order_2 = (/1,1,0/)
	               Fb = Fb - tau*r_INC*MATMUL(IP__c_PSIxPSI_g(order_2, 'L'),Gy)
	             ENDIF; IF (sDim == 3) THEN
	               order_2 = (/1,0,1/)
	               Fb = Fb - tau*r_INC*MATMUL(IP__c_PSIxPSI_g(order_2, 'L'),Gz)
	             ENDIF    

                ELSE

                     ! Boundary Flux -------------------------------------------------------------------------------------
	             order_1 = (/1,0,0/)
                     order_2 = (/2,0,0/)
	             
                     Fb = rb*(I__c_PSI_g(order_1,'L') + (0.5*dt-tau)*MATMUL(I__c_PSIxPSI_g(order_1,'L'),Gt)  &
	        	
	        	- tau*MATMUL(I__c_PSIxPSI_g(order_2,'L'),Gx))

                     IF (sDim >= 2) THEN
	               order_2 = (/1,1,0/)
	               Fb = Fb - tau*rb*MATMUL(I__c_PSIxPSI_g(order_2,'L'),Gy)
	             ENDIF; IF (sDim == 3) THEN
	               order_2 = (/1,0,1/)
	               Fb = Fb - tau*rb*MATMUL(I__c_PSIxPSI_g(order_2,'L'),Gz)
	             ENDIF	   

                ENDIF

	        Fb = MATMUL(TRANSPOSE(RR),Fb)

                rhs(:,jd_jb(j)) = rhs(:,jd_jb(j)) - Fb * SQRT(SUM(normal(:,j)**2))
	   
	 ENDIF

       ENDDO
       
    ENDDO

    END SUBROUTINE  kinetic_boundary_flux





    FUNCTION  imposed_boundary_values(ww_b, bdata, value_type)   RESULT (bb)   

    ! Computes the value of the conserved variables corresponding to the imposed      
    ! boundary data                                                                   
    !---------------------------------------------------------------------------------
    USE euler_equations                                                               
                                                                                      
    IMPLICIT NONE                                                                     

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: bdata
    INTEGER                     , INTENT(IN) :: value_type

    REAL(KIND=8), DIMENSION(SIZE(ww_b,1), SIZE(ww_b,2)) :: bb

    ! Variables needed for primitive variables manipulation                           
    REAL(KIND=8), DIMENSION(SIZE(ww_b,1) - 2) :: uu
    REAL(KIND=8) :: rho, pres

    INTEGER :: i, p                                                                 
    !---------------------------------------------------------------------------------

    p = SIZE(ww_b,1)                                                                  

    SELECT CASE (value_type)
 
       ! No imposed boundary values                                                   
       CASE(BV__NONE)                                                                 

          DO i = 1, SIZE(bb,2)                                                        
            bb(:,i) = ww_b(:,i)                                                
          ENDDO                                                                       


       ! Asymptotic constant conservative variable                                                    
       CASE(BV__WW_OO)                                                 

          DO i = 1, SIZE(bb,2)                                                        
            bb(:,i) = bdata(:,1)                                               
          ENDDO


       ! Conservative variable distribution is                                        
       ! set equal to the initial condition                                           
       CASE(BV__WW_0)                                                                 

          DO i = 1, SIZE(bb,2)                                                        
            bb(:,i) = bdata(:,i)                                               
          ENDDO


       ! Assigned constant conservative variable
       CASE (BV__CONST_WW)

          DO i = 1, SIZE(bb,2)        
            bb(:,i) = bdata(:,1) 
          ENDDO


       ! Constant pressure imposed                                                    
       CASE (BV__CONST_P)                                                             

          pres = bdata(1,1)                                                           

          DO i = 1, SIZE(bb,2)                                                        

             rho = ww_b(1,i)                                                          
             uu  = ww_b(2:p-1,i)/ww_b(1,i)                                            

             bb(1,    i) = rho                                                        
             bb(2:p-1,i) = rho * uu                                                   
             bb(p,    i) = Etot__rho_m_P(rho, rho*uu, pres)                           

          ENDDO


       ! Constant pressure and density imposed                                        
       CASE (BV__CONST_RHO_P)                                                         

          rho  = bdata(1,1)                                                           
          pres = bdata(2,1)                                                           

          DO i = 1, SIZE(bb,2)                                                        

             uu  = ww_b(2:p-1,i)/ww_b(1,i)                                            

             bb(1,    i) = rho                                                        
             bb(2:p-1,i) = rho * uu                                                   
             bb(p,    i) = Etot__rho_m_P(rho, rho*uu, pres)                           

          ENDDO                                                                       


       ! Conservative variables fully determined by                                   
       ! means of a (constant) primitive variable                                     
       ! set: density, velocity, pressure                                             
       CASE (BV__CONST_RHO_U_P)                                                       

          rho  = bdata(1,1)                                                           
          uu   = bdata(2:p-1,1)                                                       
          pres = bdata(p,1)                                                           

          DO i = 1, SIZE(bb,2)                                                        

             bb(1,    i) = rho                                                        
             bb(2:p-1,i) = rho * uu                                                   
             bb(p,    i) = Etot__rho_m_P(rho, rho*uu, pres)                           

          ENDDO


       CASE DEFAULT                                                                   

         PRINT*, '';   PRINT*, 'ERROR. IMPOSED_BOUNDARY_VALUES:'                      
         PRINT*, 'Unknown boundary value.';   PRINT*, ''                              
         STOP                                                                         

    END SELECT                                                                        

    END FUNCTION  imposed_boundary_values            





    FUNCTION  imposed_boundary_values_GG(ww_b)   RESULT (bb)   

    ! Computes the value of the conserved variables corresponding to the imposed      
    ! boundary data                                                                   
    !---------------------------------------------------------------------------------
    USE euler_equations                                                               
                                                                                      
    IMPLICIT NONE                                                                     

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b

    REAL(KIND=8), DIMENSION(SIZE(ww_b,1), SIZE(ww_b,2)) :: bb

    REAL(KIND=8), DIMENSION(:,:), POINTER :: bdata

    INTEGER, DIMENSION(:), POINTER :: bpoints

    ! Variables needed for primitive variables manipulation                           
    REAL(KIND=8), DIMENSION(SIZE(ww_b,1) - 2) :: uu
    REAL(KIND=8) :: rho, pres

    INTEGER :: value_type, Np_bound
    INTEGER :: i, p, b
    !---------------------------------------------------------------------------------

    p = SIZE(ww_b,1)                                                                  

    DO b = 1, SIZE(bc_table)
 
       value_type =  bc_table(b) % value_type
       bdata      => bc_table(b) % bdata
       Np_bound   =  bc_table(b) % Np_bound
       bpoints    => bc_table(b) % bpoints

       SELECT CASE (value_type)
 
          ! No imposed boundary values                                                
          CASE(BV__NONE)                                                              

             DO i = 1, SIZE(bb,2)                                                     
               bb(:,i) = ww_b(:,i)                                             
             ENDDO                                                                    


          ! Asymptotic constant conservative variable                                                 
          CASE(BV__WW_OO)                                              

             DO i = 1, SIZE(bb,2)                                                     
               bb(:,i) = bdata(:,1)                                            
             ENDDO


          ! Conservative variable distribution is                                     
          ! set equal to the initial condition                                        
          CASE(BV__WW_0)                                                              

             DO i = 1, SIZE(bb,2)                                                     
               bb(:,i) = bdata(:,i)                                            
             ENDDO


          ! Assigned constant conservative variable
          CASE (BV__CONST_WW)

             DO i = 1, SIZE(bb,2)        
               bb(:,i) = bdata(:,1) 
             ENDDO


          ! Constant pressure imposed                                                 
          CASE (BV__CONST_P)                                                          

             pres = bdata(1,1)                                                        

             DO i = 1, SIZE(bb,2)                                                     

                rho = ww_b(1,i)                                                       
                uu  = ww_b(2:p-1,i)/ww_b(1,i)                                         

                bb(1,    i) = rho                                                     
                bb(2:p-1,i) = rho * uu                                                
                bb(p,    i) = Etot__rho_m_P(rho, rho*uu, pres)                        

             ENDDO


          ! Constant pressure and density imposed                                     
          CASE (BV__CONST_RHO_P)                                                      

             rho  = bdata(1,1)                                                        
             pres = bdata(2,1)                                                        

             DO i = 1, SIZE(bb,2)                                                     

                uu  = ww_b(2:p-1,i)/ww_b(1,i)                                         

                bb(1,    i) = rho                                                     
                bb(2:p-1,i) = rho * uu                                                
                bb(p,    i) = Etot__rho_m_P(rho, rho*uu, pres)                        

             ENDDO                                                                    


          ! Conservative variables fully determined by                                
          ! means of a (constant) primitive variable                                  
          ! set: density, velocity, pressure                                          
          CASE (BV__CONST_RHO_U_P)                                                    

             rho  = bdata(1,1)                                                        
             uu   = bdata(2:p-1,1)                                                    
             pres = bdata(p,1)                                                        

             DO i = 1, SIZE(bb,2)                                                     

                bb(1,    i) = rho                                                     
                bb(2:p-1,i) = rho * uu                                                
                bb(p,    i) = Etot__rho_m_P(rho, rho*uu, pres)                        

             ENDDO


          CASE DEFAULT                                                                

            PRINT*, '';   PRINT*, 'ERROR. IMPOSED_BOUNDARY_VALUES:'                   
            PRINT*, 'Unknown boundary value.';   PRINT*, ''                           
            STOP                                                                      

       END SELECT                                                                     
    
    DEALLOCATE (bpoints)
    
    ENDDO

    END FUNCTION  imposed_boundary_values_GG



    FUNCTION  impose_boundary_cond(bb, ww_b, normal, cond_type, cal, eqT)   RESULT(wb)     
   
    ! Explicit time integration scheme                                                 
    !--------------------------------------------------------------------------------- 
    IMPLICIT NONE                                                                      

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: bb, ww_b                               
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: normal                                 
    INTEGER,                      INTENT(IN) :: cond_type
    INTEGER,                      INTENT(IN) :: eqT                              

    REAL(KIND=8), DIMENSION(SIZE(ww_b,1), SIZE(ww_b,2)) :: wb                          
    REAL(KIND=8), DIMENSION(SIZE(ww_b,1)-2) :: nn
    
    INTEGER :: i, cal, p                    
    !--------------------------------------------------------------------------------- 
   
    SELECT CASE (cond_type)                                                            

       ! No boundary condition ---> wb == bb                                           
       CASE (BC__NONE)

          DO i = 1, SIZE(ww_b,2) 
            wb(:,i) = bb(:,i)
          ENDDO


       ! Solid Wall                                                                    
       CASE (BC__SOLID_WALL)                                                           

          IF (eqT == EQUILIBRIUM) THEN
	  
   	    DO i = 1, SIZE(ww_b,2)
   	      wb(:,i) = impose_symmetry_bc(ww_b(:,i), normal(:,i))
            ENDDO	  
	  
	  ELSE

            IF (cal == 0) THEN
          
              DO i = 1, SIZE(ww_b,2)						         
            	wb(:,i) = bb(:,i)						         
              ENDDO
              
            ELSE

              p = SIZE(ww_b,1)

              DO i = 1, SIZE(ww_b,2)						         
            		  
            	nn = normal(:,i) / SQRT(SUM(normal(:,i)**2))
   
            	wb(1,i)     = ww_b(1,i)
            	wb(2:p-1,i) = ww_b(2:p-1,i) - SUM(ww_b(2:p-1,i)*nn)*nn
            	wb(p,i)     = ww_b(p,i) - 0.5*(SUM(ww_b(2:p-1,i)**2) / ww_b(1,i)) &
            				+ 0.5*(SUM(wb  (2:p-1,i)**2) / wb  (1,i))
              ENDDO

            ENDIF
	    
	  ENDIF  


       ! Riemann-type boundary conditions                                              
       CASE (BC__RIEMANN)                                                              

          DO i = 1, SIZE(ww_b,2)                                                       
            wb(:,i) = impose_riemann_bc(ww_b(:,i), bb(:,i), normal(:,i))               
          ENDDO                                                                        


       ! Simmetry plane                                                                
       CASE (BC__SIMM_PLANE)                                                           

          DO i = 1, SIZE(ww_b,2)                                                       
            wb(:,i) = impose_symmetry_bc(ww_b(:,i), normal(:,i))                      
          ENDDO                                                                        

    END SELECT                                                                         
   
    END FUNCTION  impose_boundary_cond
    
    
    
    
    
    SUBROUTINE  read_param_kin_boundary_cond(idf) 
    !---------------------------------------------------------------------------
    USE mp_interface
    
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: idf
    INTEGER :: b, N_bound, dType, dVType, dSize
    !---------------------------------------------------------------------------

    ! Read the boundary conditions table
    IF (.NOT. MP_job  .OR.  MP_master) THEN
      READ(idf,*) N_bound
      Nb_G = N_bound
    ELSE
      READ(idf,*)
      N_bound = NbP
    ENDIF
      
    ALLOCATE (bc_table(N_bound))

    DO b = 1, Nb_G
   
       IF (.NOT. MP_job  .OR.  MP_master) THEN
   
    	 READ (idf,*) bc_table(b) % cond_type,   &
    		      bc_table(b) % value_type,  &
    		      bc_table(b) % size_bdata
   
   
    	 IF (bc_table(b) % size_bdata > 0) THEN
       
    	   ALLOCATE (bc_table(b) % bdata(bc_table(b) % size_bdata,1))
    	   READ(idf,*) bc_table(b) % bdata(:,1)
    	   
    	 ENDIF 
   
      ELSE
      
    	 IF (bP_bG(b) /= 0) THEN

    	   READ (idf,*) bc_table(bP_bG(b)) % cond_type,   &
    			bc_table(bP_bG(b)) % value_type,  &
    			bc_table(bP_bG(b)) % size_bdata
   
   
    	   IF (bc_table(bP_bG(b)) % size_bdata > 0) THEN
       
    	     ALLOCATE (bc_table(bP_bG(b)) % bdata(bc_table(bP_bG(b)) % size_bdata,1))
    	     READ(idf,*) bc_table(bP_bG(b)) % bdata(:,1)
    	     
    	   ENDIF 
    	   
        ELSE
	
	  READ(idf,*) dType, dVType, dSize
	  IF (dSize > 0) READ(idf,*)
	
	ENDIF

      
      ENDIF
   
    ENDDO
   
    END SUBROUTINE  read_param_kin_boundary_cond 


 


    SUBROUTINE  init_kinetic_boundary_cond(flow_rgme, ww_oo, ww_b, bound_p)
    !----------------------------------------------------------------------------
    USE dynamic_vector
    USE structures,       ONLY: REF
    USE gas_properties,   ONLY: gas
       
    IMPLICIT NONE
   
    INTEGER,                      INTENT(IN) :: flow_rgme
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww_oo
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b
    INTEGER,      DIMENSION(:),   INTENT(IN) :: bound_p

    TYPE(D_I_V), DIMENSION(:),  ALLOCATABLE :: bound_p_DIV, p_bound_DIV

    INTEGER :: i_, i, b, Np_bound, Np_wall
    !----------------------------------------------------------------------------

    ! Dimensionless Gas constant   
    Rgas = gas % Rgas / REF % Rgas
    equilType = flow_rgme

    ! Computes the connectivity boundary / nodes
    ! and allocates the vector bb
    ALLOCATE (bound_p_DIV(SIZE(bound_p)))
   
    DO i = 1, SIZE(bound_p)
      ALLOCATE (bound_p_DIV(i) % vec(1))
      bound_p_DIV(i) % vec(1) = bound_p(i)
    ENDDO

    ! p_bound = from boundary index returns 
    ! boundary nodes index
    ALLOCATE (p_bound_DIV(size_DIV(bound_p_DIV,3)))
    p_bound_DIV = invert_DIV(bound_p_DIV)


    DO b = 1, SIZE(bc_table)

       ! Number of boundary nodes for each boundary line
       bc_table(b) % Np_bound = SIZE(p_bound_DIV(b) % vec)

       ! Indices of the boundary nodes that belong to boundary line 'b'
       ALLOCATE (bc_table(b) % bpoints(bc_table(b)%Np_bound))
       bc_table(b) % bpoints = p_bound_DIV(b) % vec

       ! bb = conserved variables for each boundary node
       ALLOCATE (bc_table(b) % bb(SIZE(ww_oo), bc_table(b)%Np_bound))

    ENDDO
   
    DEALLOCATE (bound_p_DIV, p_bound_DIV)


    ! Initializes the boundary condition table
    DO b = 1, SIZE(bc_table)
   
       SELECT CASE (bc_table(b) % value_type)            

         CASE (BV__NONE)
            
            IF (bc_table(b) % cond_type /= BC__SOLID_WALL  .OR. &
	        flow_rgme == 0) THEN
              ALLOCATE (bc_table(b) % bdata(1,1))
              bc_table(b) % bdata(:,1) = 0.d0
            ENDIF        

         CASE (BV__WW_OO)
        
            ALLOCATE (bc_table(b) % bdata(SIZE(ww_oo),1))
            bc_table(b) % bdata(:,1) = ww_oo    
            

         CASE (BV__WW_0)
                     
            Np_bound = SIZE(bc_table(b) % bpoints)
            ALLOCATE (bc_table(b) % bdata(SIZE(ww_oo), Np_bound))
            
            DO i_ = 1, Np_bound
              i = bc_table(b) % bpoints(i_)
              bc_table(b) % bdata(:,i_) = ww_b(:,i)
            ENDDO          


         CASE(BV__CONST_WW)
         ! Read in read_param_kinetic_boundary_cond
         
	 CASE(BV__CONST_P)
         ! Read in read_param_kinetic_boundary_cond
         
	 CASE(BV__CONST_RHO_P)
         ! Read in read_param_kinetic_boundary_cond
         
	 CASE(BV__CONST_RHO_U_P)
         ! Read in read_param_kinetic_boundary_cond
	 
	 CASE DEFAULT

            PRINT*, '';   PRINT*, 'ERROR. INIT_BC_DATA'
            PRINT*, 'Unknown boundary value.';   PRINT*, ''
            STOP       
       
       END SELECT
   
    ENDDO

   
   
    ! The number of wall points is computed
    Np_wall = 0
    DO b = 1, SIZE(bc_table)
    
      IF (bc_table(b) % cond_type == BC__SOLID_WALL) THEN
        
        Np_wall = Np_wall + SIZE(bc_table(b) % bpoints)
        
        IF (SIZE(bc_table(b) % bdata,1) == 1) THEN        
          zeta   = bc_table(b) % bdata(1,1)
        ELSEIF (SIZE(bc_table(b) % bdata,1) == 2) THEN          
          zeta   = bc_table(b) % bdata(1,1)
          T_wall = bc_table(b) % bdata(2,1)
        ENDIF
         
      ENDIF
       
    ENDDO

    ALLOCATE (solid_wall_points(Np_wall))
    
    
    ! The indices of wall points is stored
    Np_wall = 0
    DO b = 1, SIZE(bc_table)
    
      IF (bc_table(b) % cond_type == BC__SOLID_WALL) THEN
      
     	 Np_bound = SIZE(bc_table(b) % bpoints)
     	 solid_wall_points(Np_wall+1: Np_wall+Np_bound) = bc_table(b) % bpoints
     	 Np_wall = Np_wall + Np_bound
	 
      ENDIF

    ENDDO

    END SUBROUTINE init_kinetic_boundary_cond





    SUBROUTINE  write_param_kin_boundary_cond(idf) 
    !----------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: idf
    INTEGER  ::  b
    !----------------------------------------------------------------------------

    ! Write the boundary conditions table
    WRITE(idf,*) ' Parameters for kinetic boundary conditions'
    WRITE(idf,*) ' ----------------------------------------'
    WRITE(idf,*) ' Number of boundaries:', SIZE(bc_table)
    WRITE(idf,*)

    DO b = 1, SIZE(bc_table)

       WRITE(idf,*) '  Boundary', b
       WRITE(idf,*) '     Boundary condition type:  ', &
                    bc_names(bc_table(b)%cond_type)
       WRITE(idf,*) '     Boundary value type:      ', &
                    bv_names(bc_table(b)%value_type)
       WRITE(idf,*)

    ENDDO

    END SUBROUTINE  write_param_kin_boundary_cond





   !============================================================
   !  Procedures for the computation of the conservative 
   !  variables for each type of boundary condition:  
   !============================================================

   FUNCTION  impose_symmetry_bc(ww_i, normal)   RESULT(ww_)

   ! Slip condition (null normal velocity) 
   ! [ uu_ = uu - uu_normal ]
   !------------------------------------------------------------ 
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww_i
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: normal

   REAL(KIND=8), DIMENSION(SIZE(ww_i)) :: ww_

   REAL(KIND=8), DIMENSION(SIZE(normal)) :: nn

   INTEGER :: p
   !------------------------------------------------------------

   ! Slip condition (null normal velocity) 
   ! is imposed [ uu_ = uu - uu_normal ]
   p = SIZE(ww_i)

   nn = normal / SQRT( SUM(normal**2) )
   
   ww_(1)     = ww_i(1)
   ww_(2:p-1) = ww_i(2:p-1) - SUM(ww_i(2:p-1)*nn)*nn
   ww_(p)     = ww_i(p) - 0.5*( SUM(ww_i(2:p-1)**2) / ww_i(1)) &
   			+ 0.5*( SUM(ww_ (2:p-1)**2) / ww_ (1)) 

   END FUNCTION  impose_symmetry_bc





   FUNCTION  impose_riemann_bc(ww_i, ww_out, normal)   RESULT (ww_)

   ! Imposes Riemann boundary conditions. 
   ! Characteristic reconstruction of the unknown
   ! The system is linearized around the state ww_i
   !------------------------------------------------------------ 
   USE euler_equations

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww_i, ww_out
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: normal

   REAL(KIND=8), DIMENSION(SIZE(ww_i)) :: ww_

   REAL(KIND=8), DIMENSION(SIZE(ww_i),SIZE(ww_i)) :: RR, LL
   
   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: alpha, alpha_Neg
   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: lambda
   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: ww
   REAL(KIND=8), DIMENSION(SIZE(normal)) :: nn
   !------------------------------------------------------------

   nn = normal / SQRT( SUM(normal**2) )

   ww = ww_i  ! <<< State used for the linearization


   ! Eigenstructure of the system linearized in ww 
       RR = right_eigenvectors__ww_nn(ww, nn)
       LL =  left_eigenvectors__ww_nn(ww, nn)
   lambda =	   eigenvalues__ww_nn(ww, nn)


   ! Jump between the initial values of the characteristic 
   ! variables. alpha(p) = L^p Dw
   alpha = MATMUL(LL, ww_out - ww_i)


   !------------------------------------------------------------ 
   !		       _
   !		      | alpha(p)  if lambda(p) < 0
   ! alpha^Neg (p) == |
   !		      |_   0	  if lambda(p) > 0
   !
   ! -----------------------------------------------------------
   
   WHERE (lambda < 0)  
     alpha_Neg = alpha
   ELSEWHERE
     alpha_Neg = 0
   END WHERE


   !------------------------------------------------------------ 
   ! Reconstruction of the (vector) unknown ww_
   !
   !   ww_  =	ww_i  +  SUM_p alpha^Neg (p)  RR(:,p)
   !
   ! alternatively (not used here),
   !
   !   ww_  = ww_out  -  SUM_p alpha^Pos (p)  RR(:,p)
   ! -----------------------------------------------------------
   
   ww_ = ww_i + MATMUL(RR, alpha_Neg)

   END FUNCTION  impose_riemann_bc





   FUNCTION  kinetic_type_of_boundary(bound)   RESULT(btype)
   !------------------------------------------------------------ 
   IMPLICIT NONE

   INTEGER, INTENT(IN)  ::  bound     
   INTEGER              ::  btype
   !------------------------------------------------------------ 

   btype = bc_table(bound) % cond_type
   
   END FUNCTION  kinetic_type_of_boundary


  END MODULE  kinetic_boundary_cond




!     FUNCTION  ww_kin_extrap(ww, j_c_d)   RESULT (ww_Bext)
!     !----------------------------------------------------------------------------
!     USE nodes,   ONLY: k_d
!     USE kinetic_extrap_boundary
! 
!     IMPLICIT NONE
! 
!     REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
!     INTEGER,      DIMENSION(:,:), INTENT(IN) :: j_c_d
! 
!     REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(c_ext)) :: ww_Bext
! 
!     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: bb
!     REAL(KIND=8), DIMENSION(:,:), POINTER :: bdata
!      
!     INTEGER :: cond_type, value_type
!     INTEGER :: i, j, is, js, itp, jtp, itm, jtm, &
!                ibp, jbp, ibm, jbm, c, c_e, j_ext, b
!     !----------------------------------------------------------------------------
! 
!     ALLOCATE (bb(SIZE(ww,1),1))
! 
!     DO c_e = 1, SIZE(c_ext)
!    
!        c = c_ext(c_e)
!        b = ext_bc(c_e)
! 
!        j_ext = j_c__ext(c_e)
! 
! 
!        ! Stencil nodes ---------------------------
!        i  = j_c_d(1,c);  j  = j_c_d(2,c)
!        is = j_c_d(3,c);  js = j_c_d(4,c)
! 
!        IF (k_d >= 2) THEN
!          itp = j_c_d(5,c);  itm = j_c_d(7,c)
!          jtp = j_c_d(6,c);  jtm = j_c_d(8,c)
!        ENDIF
! 
!        IF (k_d == 3) THEN
!          ibp = j_c_d(9,c);   ibm = j_c_d(11,c)
!          jbp = j_c_d(10,c);  jbm = j_c_d(12,c)
!        ENDIF
!        !-------------------------------------------
! 
! 
!        IF (                  js  == j_ext  .OR. &
!            (k_d >= 2  .AND.  jtp == j_ext) .OR. &
!            (k_d >= 2  .AND.  jtm == j_ext) .OR. &
!            (k_d == 3  .AND.  jbp == j_ext) .OR. &
!            (k_d == 3  .AND.  jbm == j_ext)) THEN
! 
!            IF (b == 0) THEN
! 
!              ww_Bext(:,c_e) = ww(:,j)
! 
!            ELSE
! 
!              cond_type  =  bc_table(b) % cond_type
!              value_type =  bc_table(b) % value_type
!              bdata      => bc_table(b) % bdata
! 
!              ! Conserved variables on the basis of the boundary data imposed 
!              bb = imposed_boundary_values(ww(:,j:j), bdata, value_type)
! 
!              ! Conserved variables after the imposition of the boundary conditions.
!              ww_Bext(:,c_e:c_e) = impose_boundary_cond(bb, ww(:,j:j), normD(:,j:j), cond_type, 1)       
!          
!            ENDIF
! 
! 
!        ELSEIF ((k_d >= 2  .AND.  itp == j_ext) .OR. &
!                (k_d >= 2  .AND.  itm == j_ext) .OR. &
!                (k_d == 3  .AND.  ibp == j_ext) .OR. &
!                (k_d == 3  .AND.  ibm == j_ext)) THEN
! 
!            IF (b == 0) THEN
! 
!              ww_Bext(:,c_e) = ww(:,i)
! 
!            ELSE
! 
!              cond_type  =  bc_table(b) % cond_type
!              value_type =  bc_table(b) % value_type
!              bdata      => bc_table(b) % bdata
! 
!              ! Conserved variables on the basis of the boundary data imposed 
!              bb = imposed_boundary_values(ww(:,i:i), bdata, value_type)
! 
!              ! Conserved variables after the imposition of the boundary conditions.
!              ww_Bext(:,c_e:c_e) = impose_boundary_cond(bb, ww(:,i:i), normD(:,i:i), cond_type, 1)       
!          
!            ENDIF
! 
!        ENDIF
! 
!     ENDDO
! 
!     DEALLOCATE (bb)
! 
!     END FUNCTION  ww_kin_extrap







! REFLECTED GAS - LOCAL EQUILIBRIUM (=EULER) MODEL
!l_INC = L__ww(ww_INC)
!
!r_REF = -2.d0*SQRT(PI_G*l_INC)* (r_INC*(IuN(1,1,1) + 0.5*dt*DOT_PRODUCT(IN__c_PSI_g(order_1,'L'),Gt)))
!l_REF = l_INC
!
!ww_REF(1) = r_REF
!ww_REF(2) = -r_REF*ww_INC(2)/r_INC
!
!IF (sDim >=2) ww_REF(2+1:sDim-1) = ww_INC(2+1:sDim-1)
!
!T_REF = 1.d0 / (2.d0 * l_REF)
!e_REF = e__T_v(T_REF, 1.d0/r_REF)
!
!ww_REF(SIZE(ww_REF)) = r_REF*e_REF + 0.5*r_REF*SUM(ww_REF(2:sDim-1)**2)/r_REF
!
!CALL compute_scalar_moments(ww_REF, 'R')
