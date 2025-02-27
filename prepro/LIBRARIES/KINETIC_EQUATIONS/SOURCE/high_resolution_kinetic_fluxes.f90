!=============================================================================== 
!
!      Module: high_resolution_kinetic_fluxes
!
! Description: High order kinetic fluxes.
!
!              In a local reference frame with x-axis aligned with the node-pair,
!              the flux results from the following Riemann problem for the BGK
!              equation for the particle distribution function f:
!
!              ( f_t + f_x = (g - f) / tau
!              ( 
!              ( f(x,0) = f_0 
!
!              The initial state f_0 is a suitable expansion in space, time, and
!              tau of the distribution function f. tau is the mean collision time
!              between molecules.
!              The function g is assumed to be a suitable expansion in space and 
!              time of the Maxwellian distribution function at the interface.
!
!              The macroscopic state at the interface is computed by means
!              of the compatibility constraint for t = 0 and x = x_0. [See
!              VKI_LS 1998-03 K.Xu]
! 
!              Finally the flux at the interface is the result of the fol-
!              lowing integral:
!
!	          _dt  _+oo
!               _/   _/  u * PSI * f(x_0,t,u,v,w,xi) dudvdwdxi
!	         0    -oo
!
!              where PSI is the vector of collisional invariants
!
!
!      Author: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!===============================================================================

   MODULE  high_resolution_kinetic_fluxes

   !-----------------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER, PARAMETER :: SL_NO_LIMITER      = -1, &
   			 SL_FULL_LIMITER    =  0, &
   			 SL_VAN_LEER        =  1, &
		         SL_VAN_LEER_MOD    =  2, &	 
   			 SL_MINMOD	    =  3, &
   			 SL_SUPERBEE 	    =  4, &
   			 SL_MONOCENTR 	    =  5, &
   			 SL_VAN_ALBADA      =  6, &
			 SL_VENKATAKRISHNAN =  7
   			   
   CHARACTER(*), DIMENSION(-1:7), PARAMETER  ::  sl_lim_names = &
   			      (/ 'NO LIMITER: FULL II ORDER   ', &
   				 'FULL LIMITER: I ORDER       ', &
   				 'VAN LEER LIMITER            ', &
				 'VAN LEER LIMITER MODIFIED   ', &
   				 'MINMOD LIMITER              ', &
   				 'SUPERBEE LIMITER            ', &
   				 'MONOTONIZED CENTRAL LIMITER ', &
   				 'VAN ALBADA LIMITER          ', &   				 
				 'VENKATAKRISHNAN             '/)

   INTEGER, PARAMETER :: ONE_SIDED      = -1, &
                         LINEAR_UPWIND  =  0, &
                         THREEP_UPWIND  =  1, &
                         CENTRAL_SCHEME =  2
			 
   INTEGER, PARAMETER :: ZERO = 0, &
                         ONE  = 1
			 
   INTEGER, PARAMETER :: EQUILIBRIUM = 0, &
                         QUASI_EQUIL = 1, &
			 TRANS_REGIM = 2
			 
   INTEGER, PARAMETER :: BASE_VERSION  = 0, &
                         SIMPLIFIED_VERSION = 1
                         
   INTEGER :: kFmla

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PRIVATE :: ww_i,   ww_j,   ww_ip,  ww_jp,  ww_im,  ww_jm,   &
   					               ww_is,  ww_js,  ww_itp, ww_jtp, ww_itm, ww_jtm,  &
   					               ww_ibp, ww_jbp, ww_ibm, ww_jbm

   REAL(KIND=8), PRIVATE :: Drr_ij,   Drr_iis,  Drr_jjs,  Drr_iitp, Drr_jjtp, Drr_iitm,   &
   		            Drr_iibp, Drr_jjbp, Drr_iibm, Drr_jjbm, Drr_jjtm

   REAL(KIND=8) :: kMUSCL

   INTEGER, DIMENSION(:), ALLOCATABLE ::    sl_limiter_type

   INTEGER, PRIVATE :: i, j, is, itp, itm, ibp, ibm, js, jtp, jtm, jbp, jbm
   !-----------------------------------------------------------------------------
    
   CONTAINS


   FUNCTION  BGK_HR_flux(flow_rgme, ww, j_c_d, eta, Drr, timeStep)   RESULT (phi)
   !------------------------------------------------------------------------------------------------
   USE messages,              ONLY: terminate
   USE kinetic_equations
   USE dynamic_vector
   USE transport_properties,  ONLY: Pr
   USE structures
   USE metric_coefficients,   ONLY: cosines_ETA_FV
   USE mp_interface
   
   IMPLICIT NONE
   
   INTEGER,        		   INTENT(IN) :: flow_rgme
   REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: ww
   INTEGER,	   DIMENSION(:,:), INTENT(IN) :: j_c_d
   REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: eta
   REAL(KIND=8),   DIMENSION(:,:), INTENT(IN) :: Drr
   REAL(KIND=8),   DIMENSION(:),   INTENT(IN) :: timeStep

   REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,2)) :: phi      


   REAL(KIND=8), DIMENSION(SIZE(ww,1),  SIZE(ww,1))  :: M_l, M_r, M_0, MP_l, MN_r, RR
   REAL(KIND=8), DIMENSION(SIZE(ww,1),  SIZE(eta,1)) :: s_wi, s_wj
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1)*(SIZE(ww,1)+1)/2) :: pM_l, pM_r, pM_0, sysMat

   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: Fx_L, Fy_L, Fz_L
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: Fx_R, Fy_R, Fz_R   
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: Ft_L, Ft_R
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: Gx_L, Gx_R
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: Gx, Gy, Gz, Gt
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1), 4) :: ww_CONS, wwCR
   
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: ww_0, ww_0_x, ww_0_xx
   REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: phi_ij, ML, MR

   REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ij

   REAL(KIND=8) :: mod_eta_ij, expDt, tau, dt, Prandtl, iKn,  &
   		   r_i, r_j, r_0, G_1, G_2, G_3, G_4, G_5, G_6

   INTEGER, DIMENSION(SIZE(ww,1))  :: IPIV
   INTEGER, DIMENSION(3) :: order_0, order_1, order_2

   INTEGER :: c, p, sD, info, pi   
   
   !REAL(KIND=8), DIMENSION(SIZE(ww,1)) :: wwINT, wwINTx, wwINTxx
   !REAL(KIND=8) :: epsilon, P_is, P_i, P_j, P_js, dP_I, dP_J, Dp_ij, Dp_iis, Dp_jjs
   !------------------------------------------------------------------------------------------------
       
   IF (.NOT. ALL(timeStep == timeStep(1)))     &
   CALL terminate(6,'BGK_HR_FLUX:','Local time step not allowed.')

   dt = timeStep(1)
   sD = SIZE(eta,1)

   phi    = 0.d0
   phi_ij = 0.d0
   Gt     = 0.d0

   order_0 = 0

   ! Loop on the node-pairs
   DO c = 1, SIZE(j_c_d, 2)

      ! Integrated normal vector
      eta_ij = eta(:,c)
      mod_eta_ij = SQRT(SUM(eta_ij*eta_ij))

      ! Rotation matrix to rotate from original frame
      ! into one aligned with normal vector eta_ij
      RR = 0.d0
      RR(1,1) = 1.d0
      RR(SIZE(ww,1),SIZE(ww,1)) = 1.d0

      DO p = 1, sD
   	RR(p+1, 2:1+sD) = cosines_ETA_FV(p,:,c)
      ENDDO
      
      ! Defines the extended nodes for second order accuracy, the lenght
      ! of the extended node-pairs and the values of the conserved varia
      ! bles in the nodes in the rotated frame
      CALL set_stencil(c, j_c_d, ww, Drr, RR)


      ! Due to the will to rename the reconstructed varibles at cell interface
      ! as the values at nodes, I save the values in the nodes is, i, j, js
      DO p = 1, SIZE(ww,1)
        ww_CONS(p,:) = (/ww_is(p), ww_i(p), ww_j(p), ww_js(p)/)
      ENDDO

      ! Reconstruction at LEFT and RIGHT sides of the interface
      CALL subcell_reconstruction(ww_i, ww_j, s_wi, s_wj)

      r_i  = ww_i(1)   
      CALL compute_scalar_moments(ww_i, 'L')

      r_j  = ww_j(1)
      CALL compute_scalar_moments(ww_j, 'R')


      ! Computation of the INTERFACE state ww_0
      IF (flow_rgme == TRANS_REGIM) THEN
      
	DO p = 1, SIZE(ww,1)
          wwCR(p,:) = ww_CONS(p,:)
        ENDDO
	
	! Interface state from polynomial cubic reconstruction
        !CALL cubic_reconstruction(wwCR, Drr(:,c), wwINT, wwINTx, wwINTxx)	
	!ww_0	= wwINT
	!ww_0_x  = wwINTx
	!ww_0_xx = wwINTxx		
	

        ww_0 = r_i*IP__c_PSI_g(order_0,'L')  +  r_j*IN__c_PSI_g(order_0,'R')
	
	! Finite differences derivatives:
	ww_0_x  = ( ww_CONS(:,3) - ww_CONS(:,2)) / Drr(1,c)
	ww_0_xx = ((ww_CONS(:,4) - ww_CONS(:,3)) / Drr(3,c)  - &
	           (ww_CONS(:,2) - ww_CONS(:,1)) / Drr(2,c)) / &
	           (Drr(1,c) + 0.5*Drr(2,c) + 0.5*Drr(3,c))


      ELSE
        ! Interface state from compatibility
        ww_0 = r_i*IP__c_PSI_g(order_0,'L')  +  r_j*IN__c_PSI_g(order_0,'R')
      ENDIF


      IF (L__ww(ww_0) <= 0.d0) THEN
        OPEN(1000, FILE='npCode.err')
        WRITE(1000,*) 'node pair', c
	WRITE(1000,*) 'Left state ', ww_i
	WRITE(1000,*) 'Right state', ww_j
        WRITE(1000,*) 'Intermediate state', ww_0
	CLOSE(1000)
        CALL terminate(6,'BGK_HR_Flux:','Negative temperature at interface.')
      ENDIF
      
      r_0  = ww_0(1)
      
      CALL compute_scalar_moments(ww_0,'I')



      IF (flow_rgme == TRANS_REGIM) THEN
	     	      
             ! To be used for viscous computations in the case of well resolved 			           
	     ! regions, i.e. smooth profiles (Left state = Interface = Right state)			           
             order_0 = 0										             
             M_0 = r_0 * I__c_PSIxPSI_g(order_0,'I')							             
             												             
             pi = 0											             
             DO p = 1, SIZE(ww,1)									             
               pM_0(pi+1:pi+p) = M_0(1:p,p)								             
               pi = pi + p										             
             ENDDO											             

             ! Space derivatives of Maxwellian G --------------------------------------------		           
             sysMat = pM_0										           
             Gx = ww_0_x										           
             CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gx, SIZE(ww,1), info) 			           

             IF (sD >= 2) THEN  									           
               sysMat = pM_0										           
               Gy = s_wi(:,2)										           
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gy, SIZE(ww,1), info)			           
             ENDIF;IF (sD == 3) THEN									           
               sysMat = pM_0										           
               Gz = s_wi(:,3)										           
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gz, SIZE(ww,1), info)			           
             ENDIF !-------------------------------------------------------------------------		           
             												           
             												           
             ! Time derivatives of Maxwellian G. Compatibility on Dg ------------------------		           
             order_1 = (/1, 0, 0/)									           
             Gt = - r_0 * (MATMUL(I__c_PSIxPSI_g(order_1,'I'),Gx))					           
             												           
             IF (sD >= 2) THEN  									           
               order_1 = (/0, 1, 0/)									           
               Gt = Gt - r_0 * (MATMUL(I__c_PSIxPSI_g(order_1,'I'),Gy)) 				           
             ENDIF;IF (sD == 3) THEN									           
               order_1 = (/0, 0, 1/)									           
               Gt = Gt - r_0 * (MATMUL(I__c_PSIxPSI_g(order_1,'I'),Gz)) 				           
             ENDIF											           

             sysMat = pM_0										           
             CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gt, SIZE(ww,1), info) 			           
             !-------------------------------------------------------------------------------		           


             ! Interface Knudsen Number
	     iKN = iKnudsen__ww(ww_0, ww_CONS(:,2), ww_CONS(:,3), Drr(1,c))				           
	     												           
             ! Computing mean collision time (up to second order C-E)					           
             tau = mean_collision_time(dt, ww_i, ww_j, ww_0, iKN, pM_0, Gx, Gt, ww_0_xx)		           


             ! INTERFACE FLUX ----------------------------------------------------------------  	           
             G_1 = dt											           
             G_2 = 0.5 * dt*dt - tau*dt 								        		    
             G_3 = tau*dt										           
             												           
             ! ij Node-pair numerical flux								           
             order_1 = (/1, 0, 0/)									           

             phi_ij = G_1 * r_0*I__c_PSI_g(order_1,'I')  &						           
             	    + G_2 * r_0*MATMUL(I__c_PSIxPSI_g(order_1,'I'),Gt)  				           

             order_2 = (/2, 0, 0/)									           
             phi_ij = phi_ij - G_3 * r_0*MATMUL(I__c_PSIxPSI_g(order_2,'I'),Gx) 			           

             IF (sD >= 2) THEN  									           
               order_2 = (/1, 1, 0/)									           
               phi_ij = phi_ij - G_3 * r_0*MATMUL(I__c_PSIxPSI_g(order_2,'I'),Gy)			           
             ENDIF;IF (sD == 3) THEN									           
               order_2 = (/1, 0, 1/)									           
               phi_ij = phi_ij - G_3 * r_0*MATMUL(I__c_PSIxPSI_g(order_2,'I'),Gz)			           
             ENDIF											           


             ! Prandtl Fix -------------------------------------------------------------------  	           
	     ! (Note: variable "Pr" is the dimensional value,						           
             !        variable "Prandtl" is the nondimensional one)					           
             Prandtl = Pr / REF % Pr									           

             ! Local assignment to use the same function for Prandtl					           
             ! number correction of the full flux version (To FIX)					           
             r_i = r_0;   r_j = r_0									           
             												           
             Fx_L = Gx;   Fx_R = Gx;   Fy_L = Gy;   Fy_R = Gy						           
             Fz_L = Gz;   Fz_R = Gz;   Ft_L = Gt;   Ft_R = Gt						           
             Gx_L = Gx;   Gx_R = Gx									           

             IF (Pr /= 1.d0) THEN									           
	     												           
	       CALL fix_Prandtl(phi_ij, Prandtl, r_i, r_j, r_0, dt, tau,       &			           
             			Fx_L, Fy_L, Fz_L, Fx_R, Fy_R, Fz_R, Ft_L, Ft_R, &			           
             			Gx_L, Gx_R, Gy, Gz, Gt) 						           
	     ENDIF											           


      ELSE


             ! Computing matrices M_alpha.beta_L,R,0  and  packing symmetric					  
             ! matrices as required by subroutine DSPSV 							  
             order_0 = 0											  
             M_l = r_i * I__c_PSIxPSI_g(order_0,'L')								  
             M_r = r_j * I__c_PSIxPSI_g(order_0,'R')								  
             M_0 = r_0 * I__c_PSIxPSI_g(order_0,'I')								  
      
             pi = 0												  
             DO p = 1, SIZE(ww,1)										  
               pM_l(pi+1:pi+p) = M_l(1:p,p)									  
               pM_r(pi+1:pi+p) = M_r(1:p,p)									  
               pM_0(pi+1:pi+p) = M_0(1:p,p)									  
               pi = pi + p											  
             ENDDO												  
 

             ! Spatial derivatives of F (Left/Right) ------------------------------------------------------------ 
             sysMat = pM_l											  
             Fx_L = s_wi(:,1)											  
             CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fx_L, SIZE(ww,1), info)				  

             sysMat = pM_r											  
             Fx_R = s_wj(:,1)											  
             CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fx_R, SIZE(ww,1), info)				  
             	    												  
             IF (sD >= 2) THEN  										  
             	    												  
               sysMat = pM_l											  
               Fy_L = s_wi(:,2) 										  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fy_L, SIZE(ww,1), info)				  

               sysMat = pM_r											  
               Fy_R = s_wj(:,2) 										  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fy_R, SIZE(ww,1), info)				  
             													  
             ENDIF;IF (sD == 3) THEN										  
             	    												  
               sysMat = pM_l											  
               Fz_L = s_wi(:,3) 										  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fz_L, SIZE(ww,1), info)				  

               sysMat = pM_r											  
               Fz_R = s_wj(:,3) 										  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Fz_R, SIZE(ww,1), info)				  
             													  
             ENDIF !--------------------------------------------------------------------------------------------- 

             ! Time derivative of F. Only for first nonequilibrium term in CE ----------------------------------- 
             IF (flow_rgme == QUASI_EQUIL) THEN										  
             													  
               order_1 = (/1, 0, 0/)										  

               Ft_L = - MATMUL(r_i*I__c_PSIxPSI_g(order_1,'L'), Fx_L)						  
               Ft_R = - MATMUL(r_j*I__c_PSIxPSI_g(order_1,'R'), Fx_R)						  

               IF (sD >= 2) THEN										  

             	 order_1 = (/0, 1, 0/)  									  

             	 Ft_L = Ft_L - MATMUL(r_i*I__c_PSIxPSI_g(order_1,'L'), Fy_L)					  
             	 Ft_R = Ft_R - MATMUL(r_j*I__c_PSIxPSI_g(order_1,'R'), Fy_R)					  

               ENDIF;IF (sD == 3) THEN  									  

             	 order_1 = (/0, 0, 1/)  									  

             	 Ft_L = Ft_L - MATMUL(r_i*I__c_PSIxPSI_g(order_1,'L'), Fz_L)					  
             	 Ft_R = Ft_R - MATMUL(r_j*I__c_PSIxPSI_g(order_1,'R'), Fz_R)					  

               ENDIF												  

               sysMat = pM_l											  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Ft_L, SIZE(ww,1), info)				  
               sysMat = pM_r											  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Ft_R, SIZE(ww,1), info)				  
             													  
             ENDIF !--------------------------------------------------------------------------------------------- 
             													  
             													  
             													  
             ! Spatial derivatives in x of Maxwellian G (Left/Right) -------------------------------------------- 
             Gx_L = (ww_0 - MATMUL(RR,ww(:,i))) / (0.5*Drr_ij)  						  
             Gx_R = (MATMUL(RR,ww(:,j)) - ww_0) / (0.5*Drr_ij)  						  
             													  
             sysMat = pM_0											  
             CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gx_L, SIZE(ww,1), info)				  
             sysMat = pM_0											  
             CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gx_R, SIZE(ww,1), info)				  


             ! Spatial derivatives in y and z of Maxwellian G							  
             IF (sD >= 2) THEN  										  

               order_0 = 0											  
               MP_l = r_i * IP__c_PSIxPSI_g(order_0,'L')							  
               MN_r = r_j * IN__c_PSIxPSI_g(order_0,'R')							  


               ML = MATMUL(MP_l, Fy_L)  									  
               MR = MATMUL(MN_r, Fy_R)  									  

               sysMat = pM_0											  
               Gy = (ML + MR)											  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gy, SIZE(ww,1), info)				  

             ENDIF;IF(sD == 3) THEN										  

               ML = MATMUL(MP_l, Fz_L)  									  
               MR = MATMUL(MN_r, Fz_R)  									  

               sysMat = pM_0											  
               Gz = (ML + MR)											  
               CALL DSPSV('U', SIZE(ww,1), 1, sysMat, IPIV, Gz, SIZE(ww,1), info)				  

             ENDIF !--------------------------------------------------------------------------------------------- 


             ! Computing mean collision time									      
             tau = mean_collision_time(dt, ww_i, ww_j, ww_0)

             ! Time derivative of Maxwellian G  								
             Gt = time_derivative(r_i, r_j, r_0, pM_0, dt, tau, Fx_L, Fy_L, Fz_L,  &				
                		  Fx_R, Fy_R, Fz_R, Gx_L, Gx_R, Gy, Gz, Ft_L, Ft_R, flow_rgme)


             ! INTERFACE FLUX ----------------------------------------------------------------------------------- 
             expDt = EXP(-dt/tau)										  
     
             G_1 = dt + tau*expDt - tau 									  
             G_2 = 2*tau*tau - dt*tau - expDt*(2*tau*tau + dt*tau)						  
             G_3 = tau*tau + 0.5*dt*dt - dt*tau - expDt*tau*tau 						  
             G_4 = tau*(1 - expDt)										  
     
             IF (flow_rgme == EQUILIBRIUM) THEN
	     									  
               G_5 = tau*tau - tau*expDt*(tau + dt)
	       								  
             ELSEIF (flow_rgme == QUASI_EQUIL) THEN
	     										  
               G_5 = tau*tau - tau*expDt*(tau + dt) + tau*tau*(1 - expDt)					  
               G_6 = tau*tau*(1 - expDt)
	       									  
             ENDIF												  

									  
             order_1 = (/1, 0, 0/)										  
             order_2 = (/2, 0, 0/)										  

             ! 1-D term 											  
             phi_ij = G_1 * r_0*I__c_PSI_g(order_1,'I')   &							  
  
             	    + G_2 * r_0*( MATMUL(IP__c_PSIxPSI_g(order_2,'I'), Gx_L)   &				  
             	    		+ MATMUL(IN__c_PSIxPSI_g(order_2,'I'), Gx_R))  &				  
             	    												  
             	    + G_3 * r_0*(MATMUL(I__c_PSIxPSI_g(order_1,'I'), Gt))  &					  
             	    												  
             	    + G_4 * ( r_i*IP__c_PSI_g(order_1,'L')    & 						  
             	    	    + r_j*IN__c_PSI_g(order_1,'R'))   & 						  

             	    - G_5 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_2,'L'), Fx_L)   &				  
             	    	    + MATMUL(r_j*IN__c_PSIxPSI_g(order_2,'R'), Fx_R))					  

             ! Multidimensional terms										  
             IF (sD >= 2) THEN  										  

               order_2 = (/1, 1, 0/)										  

               phi_ij = phi_ij + G_2 * r_0*(MATMUL(I__c_PSIxPSI_g(order_2,'I'), Gy))  & 			  

             	    	       - G_5 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_2,'L'), Fy_L)  &			  
             	    		       + MATMUL(r_j*IN__c_PSIxPSI_g(order_2,'R'), Fy_R))			  

             ENDIF;IF (sD == 3) THEN										  
             	    												  
               order_2 = (/1, 0, 1/)										  

               phi_ij = phi_ij + G_2 * r_0*(MATMUL(I__c_PSIxPSI_g(order_2,'I'), Gz))  & 			  

             	    	       - G_5 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_2,'L'), Fz_L)  &			  
             	    		       + MATMUL(r_j*IN__c_PSIxPSI_g(order_2,'R'), Fz_R))			  
             ENDIF												  


             IF (flow_rgme == QUASI_EQUIL) THEN										  
             													  
               order_1 = (/1, 0, 0/)										  
             													  
               phi_ij = phi_ij - G_6 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_1,'L'), Ft_L)   &			  
             	    		       + MATMUL(r_j*IN__c_PSIxPSI_g(order_1,'R'), Ft_R))			  
             	    												  
               ! Prandtl Fix (Note: Pr is the dimensional value,						  
               !	      Prandtl is the nondimensional one)						  
               Prandtl = Pr / REF % Pr  									  
             	    												  
               IF (Pr /= 1.d0)  CALL fix_Prandtl(phi_ij, Prandtl, r_i, r_j, r_0, dt,  tau,  &			  
             	    				 Fx_L, Fy_L, Fz_L, Fx_R, Fy_R, Fz_R, Ft_L,  &		  
             	    				 Ft_R, Gx_L, Gx_R, Gy, Gz, Gt) 					  
             ENDIF
	     	     
      ENDIF



      ! Averaging Flux over time step
      phi_ij = phi_ij / dt

      ! Normal flux vector in original reference
      phi_ij = MATMUL(TRANSPOSE(RR),phi_ij)

      ! Multiply by area
      phi_ij = mod_eta_ij * phi_ij

      ! Accumulate fluxes contribute in the two nodes
      phi(:,i) = phi(:,i) + phi_ij
      phi(:,j) = phi(:,j) - phi_ij

    ENDDO

    CALL release_arrays

    CONTAINS


    SUBROUTINE  cubic_reconstruction(wwCR, Drr, w, wx, wxx)

    ! It determines the coefficients of the following cubic reconstruction
    ! for the conserved variables at interface:
    ! 
    !  w(x) = w0  +  w1*x  +  w2*x^2  +  w3*x^3
    !
    ! enforcing the following conditions:
    ! 
    !  w(-Diis - Dij/2) = w_is
    !  w(-Dij/2) = w_i
    !  w(+Dij/2) = w_j
    !  w(Dij/2 + Djjs) = w_js
    !
    ! having assumed the origin of the reference system at the interface. That results
    ! in the following linear system of equations:
    !
    !  | 1   (-Diis - Dij/2)   (-Diis - Dij/2)^2   (-Diis - Dij/2)^3 |   | w0 |     | wis |
    !  | 1   (-Dij/2)	       (-Dij/2)^2	   (-Diis - Dij/2)^3 |   | w1 |  =  | wi  |
    !  | 1   (Dij/2)	       (Dij/2)^2	   (-Diis - Dij/2)^3 |   | w2 |     | wj  |
    !  | 1   (Dij/2 + Djjs)    (Dij/2 + Djjs)*2    (Dij/2 + Djjs)*3  |   | w3 |     | wjs |
    !
    !----------------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: wwCR
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: Drr

    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: w, wx, wxx   

    REAL(KIND=8), DIMENSION(4,4) :: A
    REAL(KIND=8), DIMENSION(4)   :: B

    REAL(KIND=8) :: Dij, Diis, Djjs, x_iis, x_i, x_j, x_jjs
   
    INTEGER, DIMENSION(4) :: IPIV
   
    INTEGER :: p, INFO
    !----------------------------------------------------------------------------------------

    Dij  = Drr(1)
    Diis = Drr(2)
    Djjs = Drr(3)

    x_iis = -0.5*Dij - Diis
    x_i   = -0.5*Dij
    x_j   = +0.5*Dij
    x_jjs = +0.5*Dij + Djjs

   
    DO p = 1, SIZE(w,1)

      A(1,:) = (/ 1.d0, x_iis, x_iis**2, x_iis**3 /)
      A(2,:) = (/ 1.d0, x_i,   x_i**2,   x_i**3   /)
      A(3,:) = (/ 1.d0, x_j,   x_j**2,   x_j**3   /)
      A(4,:) = (/ 1.d0, x_jjs, x_jjs**2, x_jjs**3 /)  

      B = wwCR(p,:)

      ! Solution is obtained by means of numerical subroutine. In this
      ! case a symbolic solution may speed up computation.
      CALL DGESV(4, 1, A, 4, IPIV, B, 4, INFO)
      
    	w(p) = B(1)
       wx(p) = B(2)
      wxx(p) = B(3)*2


      ! Analytical solution for the 1D case with uniform spacing --------------------------
      !   w(p) = (9*(wwCR(p,2) + wwCR(p,3)) - (wwCR(p,1) + wwCR(p,4))) / 16.d0
      !  wx(p) = (9*(wwCR(p,3) - wwCR(p,2))/8.d0 - (wwCR(p,4) - wwCR(p,1))/24.d0) / Dij
      ! wxx(p) = 2*(wwCR(p,1) - wwCR(p,2) - wwCR(p,3) + wwCR(p,4)) / (4.d0*Dij*Dij)

    ENDDO

    END SUBROUTINE  cubic_reconstruction





    FUNCTION  time_derivative(r_i, r_j, r_0, pM_0, dt, tau, Fx_L, Fy_L, Fz_L, &
                              Fx_R, Fy_R, Fz_R, Gx_L, Gx_R, Gy, Gz, Ft_L, Ft_R, flow_rgme)   RESULT (Gt)
    !---------------------------------------------------------------------------------
    USE kinetic_equations

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r_i, r_j, r_0

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: pM_0
    REAL(KIND=8),		INTENT(IN) :: tau, dt
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Fx_L, Fy_L, Fz_L, Ft_L  
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Fx_R, Fy_R, Fz_R, Ft_R
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Gx_L, Gx_R, Gy, Gz
    INTEGER,                    INTENT(IN) :: flow_rgme
    
    REAL(KIND=8), DIMENSION(SIZE(Gz)) :: Gt

    REAL(KIND=8) :: expDt, B_0, B_1, B_2, B_3, B_4, B_5

    REAL(KIND=8), DIMENSION(SIZE(pM_0)) :: sysMat
    
    INTEGER, DIMENSION(SIZE(Gz))   :: IPIV
    INTEGER, DIMENSION(3) :: order_0, order_1

    INTEGER :: sD, info
    !---------------------------------------------------------------------------------
    
    sD = SIZE(Fx_L)-2
    
    expDt = EXP(-dt/tau)
    
    B_0 = tau*(dt - tau*(1 - expDt))
    B_1 = -tau*(1 - expDt) / B_0
    B_2 = (2*tau*tau - dt*tau - expDt*(2*tau*tau + dt*tau)) / B_0
    B_3 = tau*(1 - expDt) / B_0
    
    IF (flow_rgme == EQUILIBRIUM) THEN 
    
      B_4 = (tau*(expDt*(tau + dt) - tau)) / B_0
      
    ELSEIF (flow_rgme == QUASI_EQUIL) THEN
    
      B_4 = (tau*(expDt*(tau + dt) - tau) + tau*tau*(1 - expDt)) / B_0
      B_5 = tau*tau*(1 - expDt) / B_0
      
    ENDIF
   
   order_0 = 0
   order_1 = (/1,0,0/)
   
   ! One dimensional term
    Gt = B_1 * r_0*I__c_PSI_g(order_0,'I')   &

       + B_2 * r_0*( MATMUL(IP__c_PSIxPSI_g(order_1,'I'), Gx_L)   &
                   + MATMUL(IN__c_PSIxPSI_g(order_1,'I'), Gx_R))  &
                                                                                
       + B_3 * ( r_i*IP__c_PSI_g(order_0,'L')  &
               + r_j*IN__c_PSI_g(order_0,'R')) &


       + B_4 * ( r_i*MATMUL(IP__c_PSIxPSI_g(order_1,'L'), Fx_L)  &
               + r_j*MATMUL(IN__c_PSIxPSI_g(order_1,'R'), Fx_R))

                                                        
    IF (sD >= 2) THEN                                                           

      order_1 = (/0, 1, 0/)                                                     

      Gt = Gt + B_2 * r_0*(MATMUL(I__c_PSIxPSI_g(order_1,'I'), Gy))  &  

              + B_4 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_1,'L'), Fy_L)  &  
                      + MATMUL(r_j*IN__c_PSIxPSI_g(order_1,'R'), Fy_R))  

    ENDIF;IF (sD == 3) THEN                                                     

      order_1 = (/0, 0, 1/)                                                     

      Gt = Gt + B_2 * r_0*(MATMUL(I__c_PSIxPSI_g(order_1,'I'), Gz))  &  

              + B_4 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_1,'L'), Fz_L)  &  
                      + MATMUL(r_j*IN__c_PSIxPSI_g(order_1,'R'), Fz_R))  
    ENDIF                                                                       


    IF (flow_rgme == QUASI_EQUIL) THEN                                                                          
                                                                                                        
      order_1 = (/1, 0, 0/)                                                                             
                                                                                                        
      Gt = Gt + B_5 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_1,'L'), Ft_L)   &                      
                      + MATMUL(r_j*IN__c_PSIxPSI_g(order_1,'R'), Ft_R))
    ENDIF

    sysMat = pM_0
    CALL DSPSV('U', SIZE(Gt,1), 1, sysMat, IPIV, Gt, SIZE(Gt), info)
					 
    END FUNCTION  time_derivative





    SUBROUTINE  fix_Prandtl(phi_ij, PrandtlNo, r_i, r_j, r_0, dt, tau,      &
                            Fx_L, Fy_L, Fz_L, Fx_R, Fy_R, Fz_R, Ft_L, Ft_R, &
			    Gx_L, Gx_R, Gy, Gz, Gt)
    !-----------------------------------------------------------------------------------------       
    USE kinetic_equations		   
    
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: phi_ij
    
    REAL(KIND=8),               INTENT(IN) :: PrandtlNo
    REAL(KIND=8),               INTENT(IN) :: r_i, r_j, r_0
    REAL(KIND=8),		INTENT(IN) :: tau, dt
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Fx_L, Fy_L, Fz_L  
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Fx_R, Fy_R, Fz_R
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Ft_L, Ft_R
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Gx_L, Gx_R, Gy, Gz, Gt

    REAL(KIND=8), DIMENSION(SIZE(phi_ij)) :: psi_ij

    REAL(KIND=8) :: expDt, G_1, G_2, G_3, G_4, G_5, G_6, U, V, W, q_ij

    INTEGER, DIMENSION(3) :: order_0, order_1
    !-----------------------------------------------------------------------------------------
    
    expDt = EXP(-dt/tau)
     
    G_1 = dt + tau*expDt - tau
    G_2 = 2*tau*tau - dt*tau - expDt*(2*tau*tau + dt*tau)
    G_3 = tau*tau + 0.5*dt*dt - dt*tau - expDt*tau*tau
    G_4 = tau*(1 - expDt)
    G_5 = tau*tau - tau*expDt*(tau + dt) + tau*tau*(1 - expDt)
    G_6 = tau*tau*(1 - expDt)
   
    order_0 = 0
    order_1 = (/1, 0, 0/)

    ! 1-D term
    psi_ij = G_1 * r_0*I__c_PSI_g(order_0,'I')   &

    	   + G_2 * r_0*( MATMUL(IP__c_PSIxPSI_g(order_1,'I'), Gx_L)   &
    		       + MATMUL(IN__c_PSIxPSI_g(order_1,'I'), Gx_R))  &
           
	   + G_3 * r_0*(MATMUL(I__c_PSIxPSI_g(order_0,'I'), Gt))  &
    	   
	   + G_4 * ( r_i*IP__c_PSI_g(order_0,'L')   &
    		   + r_j*IN__c_PSI_g(order_0,'R'))  &

    	   - G_5 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_1,'L'), Fx_L)   &
    		   + MATMUL(r_j*IN__c_PSIxPSI_g(order_1,'R'), Fx_R))  &

    	   - G_6 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_0,'L'), Ft_L)  &
    		   + MATMUL(r_j*IN__c_PSIxPSI_g(order_0,'R'), Ft_R))

    IF (sD >= 2) THEN									    

      order_1 = (/0, 1, 0/)									    

      psi_ij = psi_ij + G_2 * r_0*(MATMUL(I__c_PSIxPSI_g(order_1,'I'), Gy))  &  		    

    		      - G_5 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_1,'L'), Fy_L)  &		    
    			      + MATMUL(r_j*IN__c_PSIxPSI_g(order_1,'R'), Fy_R)) 		    

    ENDIF;IF (sD == 3) THEN								    
    												    
      order_1 = (/0, 0, 1/)									    

      psi_ij = psi_ij + G_2 * r_0*(MATMUL(I__c_PSIxPSI_g(order_1,'I'), Gz))   & 		    

    		      - G_5 * ( MATMUL(r_i*IP__c_PSIxPSI_g(order_1,'L'), Fz_L)  &		    
    			      + MATMUL(r_j*IN__c_PSIxPSI_g(order_1,'R'), Fz_R)) 		    
    ENDIF											    


    ! 1-D term (quasi, since the first and fourth term at the right hand side			    
    !		accounts for multidimensionality directly)					    
    U = psi_ij(2) / psi_ij(1)
    								    
    q_ij =   phi_ij(SIZE(ww,1)) + 0.5*U*U*  phi_ij(1) - U*  phi_ij(2)  &			    
         - U*psi_ij(SIZE(ww,1)) - 0.5*U*U*U*psi_ij(1) + U*U*psi_ij(2)

    IF (sD >= 2) THEN									    
    
      V = psi_ij(3) / psi_ij(1)
      
      q_ij = q_ij + 0.5*V*V*  phi_ij(1) - V*  phi_ij(3)  &
                  - 0.5*U*V*V*psi_ij(1) + U*V*psi_ij(3)
    
    ENDIF;IF (sD == 3) THEN					    

      W = psi_ij(4) / psi_ij(1)
      
      q_ij = q_ij + 0.5*W*W*  phi_ij(1) - W*  phi_ij(4)  &
                  - 0.5*U*W*W*psi_ij(1) + U*W*psi_ij(4)							   
    
    ENDIF											    

    phi_ij(SIZE(phi_ij)) = phi_ij(SIZE(phi_ij)) + (1.d0/PrandtlNo - 1) * q_ij

    END SUBROUTINE  fix_Prandtl


    END FUNCTION  BGK_HR_flux





   SUBROUTINE  set_stencil(node_pair, j_c_d, ww, Drr, RR)
   !---------------------------------------------------------------------------
   USE messages,            ONLY: terminate
   USE kinetic_equations,   ONLY: L__ww
   USE mp_interface
   
   IMPLICIT NONE
   
   INTEGER,                      INTENT(IN) :: node_pair
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: j_c_d
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Drr
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: RR

   REAL(KIND=8) :: L_i, L_j
   INTEGER :: c, sD
   CHARACTER(LEN=2) :: fex
   !---------------------------------------------------------------------------

   sD = SIZE(ww,1) - 2

   IF (.NOT. ALLOCATED(ww_i)) THEN

     ALLOCATE (ww_i(SIZE(ww,1)),  ww_j(SIZE(ww,1)),   &
               ww_is(SIZE(ww,1)), ww_js(SIZE(ww,1)),  &
	       ww_ip(SIZE(ww,1)), ww_jp(SIZE(ww,1)),  &
	       ww_im(SIZE(ww,1)), ww_jm(SIZE(ww,1)))

     IF (sD >= 2)  ALLOCATE (ww_itp(SIZE(ww,1)), ww_jtp(SIZE(ww,1)),  &
                             ww_itm(SIZE(ww,1)), ww_jtm(SIZE(ww,1)))
   
     IF (sD == 3)  ALLOCATE (ww_ibp(SIZE(ww,1)), ww_jbp(SIZE(ww,1)),  &
                             ww_ibm(SIZE(ww,1)), ww_jbm(SIZE(ww,1)))
   ENDIF

   c = node_pair

   ! Node-pair and extended nodes
   i = j_c_d(1,c);   is = j_c_d(3,c)
   j = j_c_d(2,c);   js = j_c_d(4,c)

   ww_i = MATMUL(RR, ww(:,i));   ww_is = MATMUL(RR, ww(:,is))
   ww_j = MATMUL(RR, ww(:,j));   ww_js = MATMUL(RR, ww(:,js))

   Drr_ij  = Drr(1,c)
   Drr_iis = Drr(2,c)   
   Drr_jjs = Drr(3,c)

   IF (sD >= 2) THEN

     itp = j_c_d(5,c);   itm = j_c_d(7,c)
     jtp = j_c_d(6,c);   jtm = j_c_d(8,c)

     ww_itp = MATMUL(RR, ww(:,itp));   ww_itm = MATMUL(RR, ww(:,itm))
     ww_jtp = MATMUL(RR, ww(:,jtp));   ww_jtm = MATMUL(RR, ww(:,jtm))

     Drr_iitp = Drr(4,c);   Drr_iitm = Drr(6,c)
     Drr_jjtp = Drr(5,c);   Drr_jjtm = Drr(7,c)

   ENDIF;IF (sD == 3) THEN

     ibp = j_c_d(9,c);   ibm = j_c_d(11,c)
     jbp = j_c_d(10,c);  jbm = j_c_d(12,c)

     ww_ibp = MATMUL(RR, ww(:,ibp));   ww_ibm = MATMUL(RR, ww(:,ibm))
     ww_jbp = MATMUL(RR, ww(:,jbp));   ww_jbm = MATMUL(RR, ww(:,jbm))

     Drr_iibp = Drr(8,c);   Drr_iibm = Drr(10,c)
     Drr_jjbp = Drr(9,c);   Drr_jjbm = Drr(11,c)

   ENDIF

   L_i = L__ww(ww_i)
   L_j = L__ww(ww_j)

   IF (L_i <= 0.d0) THEN
     IF (.NOT. MP_job) THEN   
       OPEN (1000, FILE='npCode.err')
       WRITE(1000,*) ' pair', c, ' node i', i, 'extensions', is, j, js
       WRITE(1000,*) 'ww', ww(:,i)
       CLOSE(1000)
     ELSE
       WRITE(fex,'(i2)') MP_pRank
       OPEN (1000, FILE='npCode_P'//fex//'.err')
       WRITE(1000,*) ' pair', cG_cP_FV(c), ' node i', jG_jP(i),  &
                              'extensions', jG_jP(is), & 
			                    jG_jP(j),  &
					    jG_jP(js)
       WRITE(1000,*) 'ww', ww(:,i)
       CLOSE(1000)     
     ENDIF  
     CALL terminate(6, 'SET_STENCIL:', 'Negative temperature.')
   ENDIF  

   IF (L_j <= 0.d0)  THEN
     IF (.NOT. MP_job) THEN   
       OPEN (1000, FILE='npCode.err')
       WRITE(1000,*) ' pair', cG_cP_FV(c), ' node j', j, 'extensions', js, i, is
       WRITE(1000,*) 'ww', ww(:,j)
       CLOSE(1000)
     ELSE
       WRITE(fex,'(i2)') MP_pRank
       OPEN (1000, FILE='npCode_P'//fex//'.err')
       WRITE(1000,*) ' pair', c, ' node j', jG_jP(j),  &
                              'extensions', jG_jP(js), & 
			                    jG_jP(i),  &
					    jG_jP(is)
       WRITE(1000,*) 'ww', ww(:,j)
       CLOSE(1000)
     ENDIF  
     CALL terminate(6, 'SET_STENCIL:', 'Negative temperature.')
   ENDIF  

   END SUBROUTINE  set_stencil







   SUBROUTINE  subcell_reconstruction(ww_L, ww_R, slope_L, slope_R)   
   !--------------------------------------------------------------------------------
   USE kinetic_equations,   ONLY: L__ww
      
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: ww_L, ww_R
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: slope_L, slope_R

   REAL(KIND=8), DIMENSION(SIZE(ww_L), SIZE(ww_L)-2) :: ww_ipH, ww_imH, ww_jpH, ww_jmH

   REAL(KIND=8), DIMENSION(SIZE(ww_L)) :: Dww_ij, Dww_iis, Dww_jjs,  &
                                          Dww_iitp, Dww_iitm, Dww_iibp, Dww_iibm,  &
                                          Dww_jjtp, Dww_jjtm, Dww_jjbp, Dww_jjbm
   REAL(KIND=8) :: a, b
   					     
   INTEGER :: sD, p, limType  
   !--------------------------------------------------------------------------------
   
   sD = SIZE(ww_i) - 2
   
   Dww_jjs = ww_js - ww_j
   Dww_ij  = ww_j  - ww_i
   Dww_iis = ww_i  - ww_is
   
   IF (sD >= 2) THEN
     
     Dww_iitm = ww_itm - ww_i
     Dww_iitp = ww_i   - ww_itp       
     
     Dww_jjtp = ww_jtp - ww_j
     Dww_jjtm = ww_j   - ww_jtm

   ENDIF;IF (sD == 3) THEN
   
     Dww_iibp = ww_ibp - ww_i
     Dww_iibm = ww_i   - ww_ibm       
     
     Dww_jjbp = ww_jbp - ww_j
     Dww_jjbm = ww_j   - ww_jbm

   ENDIF


   ! LEFT side of interface     
   ! Computing Limited derivative as the product of a limiter function
   ! and an approximation of the derivative of the variables in the cell
   DO p = 1, SIZE(ww_i)
   
     a = Dww_iis(p)*(Drr_ij/Drr_iis)
     b = Dww_ij(p)
     
     slope_L(p,1) = limFun(a, b, sl_limiter_type(p), 'L', p)  &
                  * 0.5*((1+kMUSCL)*Dww_ij(p)/Drr_ij + (1-kMUSCL)*Dww_iis(p)/Drr_iis)

     IF (sD >= 2) THEN
     
       a = Dww_iitp(p)*(Drr_iitm/Drr_iitp)
       b = Dww_iitm(p)

       IF (sl_limiter_type(p) == SL_VENKATAKRISHNAN) THEN
         limType = SL_VAN_LEER
       ELSE
         limType = sl_limiter_type(p)
       ENDIF

       slope_L(p,2) = limFun(a, b, limType)  &
                    * 0.5*(Dww_iitm(p)/Drr_iitm + Dww_iitp(p)/Drr_iitp)
		     
     ENDIF;IF (sD == 3) THEN
     
       a = Dww_iibm(p)*(Drr_iibp/Drr_iibm)
       b = Dww_iibp(p)
       
       IF (sl_limiter_type(p) == SL_VENKATAKRISHNAN) THEN
         limType = SL_VAN_LEER
       ELSE
         limType = sl_limiter_type(p)
       ENDIF
     
       slope_L(p,3) = limFun(a, b, limType)  &
                    * 0.5*(Dww_iibm(p)/Drr_iibm + Dww_iibp(p)/Drr_iibp)
     ENDIF
     
   ENDDO


   ! RIGHT side of interface     
   ! Computing Limited derivative as the product of a limiter function
   ! and an approximation of the derivative of the variables in the cell
   DO p = 1, SIZE(ww_j)     
     
     a = Dww_ij(p)
     b = Dww_jjs(p)*(Drr_ij/Drr_jjs)

     slope_R(p,1) = limFun(a, b, sl_limiter_type(p), 'R', p)  &
                  * 0.5*((1+kMUSCL)*Dww_ij(p)/Drr_ij + (1-kMUSCL)*Dww_jjs(p)/Drr_jjs)

     IF (sD >= 2) THEN
     
       a = Dww_jjtm(p)*(Drr_jjtp/Drr_jjtm)
       b = Dww_jjtp(p)
       
       IF (sl_limiter_type(p) == SL_VENKATAKRISHNAN) THEN
         limType = SL_VAN_LEER
       ELSE
         limType = sl_limiter_type(p)
       ENDIF
     
       slope_R(p,2) = limFun(a, b, limType)  &
                    * 0.5*(Dww_jjtm(p)/Drr_jjtm + Dww_jjtp(p)/Drr_jjtp)
		     
     ENDIF;IF (sD == 3) THEN
     
       a = Dww_jjbm(p)*(Drr_jjbp/Drr_jjbm)
       b = Dww_jjbp(p)
       
       IF (sl_limiter_type(p) == SL_VENKATAKRISHNAN) THEN
         limType = SL_VAN_LEER
       ELSE
         limType = sl_limiter_type(p)
       ENDIF
     
       slope_R(p,3) = limFun(a, b, limType)  &
                    * 0.5*(Dww_jjbm(p)/Drr_jjbm + Dww_jjbp(p)/Drr_jjbp)
     ENDIF

   ENDDO


   ! LEFT reconstruction check
   ww_ipH(:,1) = ww_i + slope_L(:,1) * Drr_ij/2.d0
   ww_imH(:,1) = ww_i - slope_L(:,1) * Drr_iis/2.d0
      
   IF (L__ww(ww_ipH(:,1)) <= 0.d0 .OR. &
       L__ww(ww_imH(:,1)) <= 0.d0) THEN	
     slope_L(:,1) = 0.d0
     ww_L = ww_i
   ELSE
     ww_L = ww_ipH(:,1)
   ENDIF
   
   ! RIGHT reconstruction check
   ww_jpH(:,1) = ww_j + slope_R(:,1) * Drr_jjs/2.d0
   ww_jmH(:,1) = ww_j - slope_R(:,1) * Drr_ij/2.d0
   
   IF (L__ww(ww_jpH(:,1)) <= 0.d0 .OR. &
       L__ww(ww_jmH(:,1)) <= 0.d0) THEN	
     slope_R(:,1) = 0.d0
     ww_R = ww_j
   ELSE
     ww_R = ww_jmH(:,1)
   ENDIF 
   
   
   IF (sD >= 2) THEN   
   
     ! LEFT slope check
     ww_ipH(:,2) = ww_L + slope_L(:,2) * Drr_iitm/2.d0
     ww_imH(:,2) = ww_L - slope_L(:,2) * Drr_iitp/2.d0
     	
     IF (L__ww(ww_ipH(:,2)) <= 0.d0 .OR. &
     	 L__ww(ww_imH(:,2)) <= 0.d0)  slope_L(:,2) = 0.d0
   
     ! RIGHT slope check
     ww_jpH(:,2) = ww_R + slope_R(:,2) * Drr_jjtp/2.d0
     ww_jmH(:,2) = ww_R - slope_R(:,2) * Drr_jjtm/2.d0
   
     IF (L__ww(ww_jpH(:,2)) <= 0.d0 .OR. &
     	 L__ww(ww_jmH(:,2)) <= 0.d0)  slope_R(:,2) = 0.d0

   ENDIF;IF (sD == 3) THEN     
   
     ! LEFT slope check
     ww_ipH(:,3) = ww_L + slope_L(:,3) * Drr_iibp/2.d0
     ww_imH(:,3) = ww_L - slope_L(:,3) * Drr_iibm/2.d0
     	
     IF (L__ww(ww_ipH(:,3)) <= 0.d0 .OR. &
     	 L__ww(ww_imH(:,3)) <= 0.d0)  slope_L(:,3) = 0.d0
   
     ! RIGHT slope check
     ww_jpH(:,3) = ww_R + slope_R(:,3) * Drr_jjbp/2.d0
     ww_jmH(:,3) = ww_R - slope_R(:,3) * Drr_jjbm/2.d0
   
     IF (L__ww(ww_jpH(:,3)) <= 0.d0 .OR. &
     	 L__ww(ww_jmH(:,3)) <= 0.d0)  slope_R(:,3) = 0.d0

   ENDIF
   
   END SUBROUTINE  subcell_reconstruction







   FUNCTION  limFun(a, b, limiter_type, state, eq)   RESULT(s)

   ! Limiters in standard form
   !------------------------------------------------------------
   USE messages,   ONLY: terminate
   
   IMPLICIT NONE

   REAL(KIND=8),     INTENT(IN) :: a, b
   INTEGER,	     INTENT(IN) :: limiter_type
   
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: state
   INTEGER,          INTENT(IN), OPTIONAL :: eq
   
   REAL(KIND=8) :: s
   !------------------------------------------------------------

   SELECT CASE (limiter_type)
   
     CASE (SL_NO_LIMITER)
     s = 1
 
     CASE (SL_FULL_LIMITER)
     s = 0
 
     CASE (SL_VAN_LEER)
     s = 2 * (a*b + ABS(a*b)) / ((a+b)**2 + 1.0d-8)
     
     CASE (SL_VAN_LEER_MOD)
     s = (SIGN(1.d0,(a+1.e-10)) + SIGN(1.d0,(b+1.e-10)))*(ABS(a)*ABS(b))/(ABS(a)+ABS(b)+1.e-10)

     CASE (SL_MINMOD)
     s = 0.5 * (SIGN(1.d0,a) + SIGN(1.d0,b)) * MIN(ABS(a), ABS(b))

     CASE (SL_SUPERBEE)
     s = 0.5 * MAX(0.d0, MIN(2*(a/(b+1.e-8)),1.d0), MIN(2.d0,a/(b+1.e-8) ))

     CASE (SL_MONOCENTR)
     s = 0.5 * MAX(0.d0, MIN(2.d0, 2*(a/(b+1.e-8)), (1+a/(b+1.e-8))/2.d0))

     CASE (SL_VAN_ALBADA)
     s = (a*b + ABS(a*b)) / (a**2 + b**2 + 1.0d-12)

     CASE (SL_VENKATAKRISHNAN)
     s = venkaFun(eq, state)     

     CASE DEFAULT
     CALL terminate(6,'BGK-HR', 'Unknown limiter.')

   END SELECT


   CONTAINS
   
   
   FUNCTION venkaFun(eq, state)   RESULT(psi)
   !--------------------------------------------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER,          INTENT(IN) :: eq
   CHARACTER(LEN=1), INTENT(IN) :: state
   
   REAL(KIND=8) :: psi

   REAL(KIND=8), DIMENSION(3) :: wMax_A, wMin_A

   REAL(KIND=8) :: wwMax, wwMin, DrMax, DrMin, wpH, wmH, wwI 
   REAL(KIND=8) :: PHI_pH, PHI_mH, Dm, Dp, sqEps, kVENK = 0.3, DrpH, DrmH, Dx
   !-------------------------------------------------------------------------------------------------------- 
     
   IF (state == 'L') THEN
	  
     wMax_A = (/ww_is(eq), ww_i(eq), ww_j(eq)/)
     wMin_A = (/ww_is(eq), ww_i(eq), ww_j(eq)/)
     	  
     wwMax = MAXVAL(wMax_A)
     wwMin = MINVAL(wMin_A)
 
     IF (ALL(MAXLOC(wMax_A) == 1)) DrMax = Drr_iis
     IF (ALL(MAXLOC(wMax_A) == 2)) DrMax = 1.d0
     IF (ALL(MAXLOC(wMax_A) == 3)) DrMax = Drr_ij
        		
     IF (ALL(MINLOC(wMax_A) == 1)) DrMin = Drr_iis
     IF (ALL(MINLOC(wMax_A) == 2)) DrMin = 1.d0
     IF (ALL(MINLOC(wMax_A) == 3)) DrMin = Drr_ij    
      
     wpH = ww_i(eq) + 0.5*((1+kMUSCL)*(ww_j(eq) - ww_i(eq))/Drr_ij + (1-kMUSCL)*(ww_i(eq) - ww_is(eq))/Drr_iis)*Drr_ij/2.d0
     wmH = ww_i(eq) - 0.5*((1+kMUSCL)*(ww_j(eq) - ww_i(eq))/Drr_ij + (1-kMUSCL)*(ww_i(eq) - ww_is(eq))/Drr_iis)*Drr_iis/2.d0
     wwI = ww_i(eq)

     DrpH = Drr_ij / 2.d0
     DrmH = Drr_iis / 2.d0

     Dx = 0.5*(Drr_iis + Drr_ij)

   ENDIF
   
   IF (state == 'R') THEN
	  
     wMax_A = (/ww_i(eq), ww_j(eq), ww_js(eq)/)
     wMin_A = (/ww_i(eq), ww_j(eq), ww_js(eq)/)
     	  
     wwMax = MAXVAL(wMax_A)
     wwMin = MINVAL(wMin_A)
 
     IF (ALL(MAXLOC(wMax_A) == 1)) DrMax = Drr_ij
     IF (ALL(MAXLOC(wMax_A) == 2)) DrMax = 1.d0
     IF (ALL(MAXLOC(wMax_A) == 3)) DrMax = Drr_jjs
        		
     IF (ALL(MINLOC(wMax_A) == 1)) DrMin = Drr_ij
     IF (ALL(MINLOC(wMax_A) == 2)) DrMin = 1.d0
     IF (ALL(MINLOC(wMax_A) == 3)) DrMin = Drr_jjs    
      
     wpH = ww_j(eq) + 0.5*((1+kMUSCL)*(ww_j(eq) - ww_i(eq))/Drr_ij + (1-kMUSCL)*(ww_js(eq) - ww_j(eq))/Drr_jjs)*Drr_jjs/2.d0
     wmH = ww_j(eq) - 0.5*((1+kMUSCL)*(ww_j(eq) - ww_i(eq))/Drr_ij + (1-kMUSCL)*(ww_js(eq) - ww_j(eq))/Drr_jjs)*Drr_ij/2.d0
     wwI = ww_j(eq)

     DrpH = Drr_jjs / 2.d0
     DrmH = Drr_ij / 2.d0

     Dx = 0.5*(Drr_jjs + Drr_ij)

   ENDIF
   
   
   ! PHI i-1/2
   IF ((wmH - wwI) == 0.d0) THEN
     PHI_mH = 1.d0
   ELSE
   
     Dm = wmH - wwI
     
     IF ((wmH - wwI) < 0.d0)  Dp = (wwMin - wwI) * DrmH/DrMin
     IF ((wmH - wwI) > 0.d0)  Dp = (wwMax - wwI) * DrmH/DrMax

     sqEps = (kVENK*Dx)**3
     
     PHI_mH = 1.d0/(SIGN(1.d0,Dm)*(ABS(Dm)+1.d-12))  &
	    * ((Dp**2 + sqEps)*Dm + 2*Dp*Dm**2)      &
	    / (Dp**2 + 2*Dm**2 + Dm*Dp + sqEps)
   ENDIF

   ! PHI i+1/2
   IF ((wpH - wwI) == 0.d0) THEN
     PHI_pH = 1.d0
   ELSE
   
     Dm = wpH - wwI

     IF ((wpH - wwI) < 0.d0)  Dp = (wwMin - wwI) * DrpH/DrMin
     IF ((wpH - wwI) > 0.d0)  Dp = (wwMax - wwI) * DrpH/DrMax
     
     sqEps = (kVENK*Dx)**3
     
     PHI_pH = 1.d0/(SIGN(1.d0,Dm)*(ABS(Dm)+1.d-12))  &
	    * ((Dp**2 + sqEps)*Dm + 2*Dp*Dm**2)      &
	    / (Dp**2 + 2*Dm**2 + Dm*Dp + sqEps)
   ENDIF       


   ! Resulting Limiter function
   psi = MIN(PHI_mH, PHI_pH)
   
   END FUNCTION venkaFun

   END FUNCTION  limFun





   SUBROUTINE  read_param_hr_kin_fluxes(idf)
   !------------------------------------------------------------
   IMPLICIT NONE
   INTEGER :: idf, ww_size
   !------------------------------------------------------------

   READ(idf,*) kFmla

   ! Constant k in MUSCL reconstruction 
   SELECT CASE (kFmla)
   
     CASE (ONE_SIDED)
     kMUSCL = -1.d0
     
     CASE (LINEAR_UPWIND)
     kMUSCL = 0.d0
      
     CASE (THREEP_UPWIND)
     kMUSCL = 1.d0/3.d0
     
     CASE (CENTRAL_SCHEME)
     kMUSCL = 1.d0
        
   END SELECT

   READ(idf,*) ww_size
   
   ALLOCATE (sl_limiter_type(ww_size))
   READ(idf,*) sl_limiter_type
   
   END SUBROUTINE  read_param_hr_kin_fluxes





   SUBROUTINE  write_param_hr_fluxes(idf)
   !------------------------------------------------------------
   IMPLICIT NONE
   INTEGER :: idf, p
   !------------------------------------------------------------

   WRITE(idf,*) '   PARAMETERS FOR HIGH-RESOLUTION NUMERICAL FLUXES'

   DO p = 1, SIZE(sl_limiter_type)
      WRITE(idf,*) '   Limiter ', p,': ', sl_lim_names(sl_limiter_type(p))
   ENDDO
   
   WRITE(idf,*)
   
   END SUBROUTINE  write_param_hr_fluxes
   
   
   
   
   
   SUBROUTINE  release_arrays
   !------------------------------------------------------------
   IMPLICIT NONE
   !------------------------------------------------------------

   DEALLOCATE (ww_i, ww_j, ww_is, ww_js, ww_ip, ww_jp, ww_im, ww_jm)

   IF (ALLOCATED(ww_itp))  DEALLOCATE (ww_itp, ww_jtp, ww_itm, ww_jtm)
   IF (ALLOCATED(ww_ibp))  DEALLOCATE (ww_ibp, ww_jbp, ww_ibm, ww_jbm)

   END SUBROUTINE  release_arrays

   
   END MODULE  high_resolution_kinetic_fluxes
   
   
   
   
! TENTATIVE OF RECONSTRUCTION BY MEANS OF PRIMITIVE VARIABLES

! Primitive variables (rho, u, P) ---------------------------------------------------
!DO p = 1, 4
!  wwCR(1,p)	  = ww_CONS(1,p)
!  wwCR(2:sD+1,p) = ww_CONS(2:sD+1,p)/ww_CONS(1,p)
!  wwCR(sD+2,p)   = P__ww(ww_CONS(:,p))
!ENDDO   ----------------------------------------------------------------------------
   
! Conservative variables from primitives --------------------------------------------
!ww_0(1)      = wwINT(1)
!ww_0(2:sD+1) = wwINT(2:sD+1)*wwINT(1)
!ww_0(sD+2)   = wwINT(sD+2)/(gamma-1) + 0.5*wwINT(1)*SUM(wwINT(2:sD+1)**2)

! First Derivatives from primitives
!ww_0_x(1)	= wwINTx(1)
!ww_0_x(2:sD+1) = wwINT(1)*wwINTx(2:sD+1) + wwINT(2:sD+1)*wwINTx(1)
!ww_0_x(sD+2)	= wwINTx(sD+2)/(gamma-1) + 0.5*wwINTx(1)*SUM(wwINT(2:sD+1)**2) &
!		+ wwINT(1)*wwINT(2)*wwINTx(2)  ! OneD

! Second Derivatives from primitives
!ww_0_xx(1)	 = wwINTxx(1)
!ww_0_xx(2:sD+1) = wwINT(1)*wwINTxx(2:sD+1) + wwINT(2:sD+1)*wwINTxx(1) &
!		 + 2*wwINTx(1)*wwINTx(2)  ! OneD
!ww_0_xx(sD+2)   = wwINTxx(sD+2)/(gamma-1)	 &
!		 + wwINT(1)*wwINT(2)*wwINTxx(2)  &  ! OneD
!		 + 0.5*wwINTxx(1)*wwINT(2)**2	 &  ! OneD
!		 + ww_0_x(2)*wwINTx(2)  	 &  ! OneD
!		 + wwINT(2)*wwINTx(1)*wwINTx(2)     ! OneD
! ------------------------------------------------------------------------------------

   
   
   
! SOME COLLECTION OF USEFUL DEBUG PRINTS

! print*, ''
! print*, 'interface', c
! print*, i, ww_i
! print*, '0', ww_0
! print*, j, ww_j
! print*, ''
! print*, 'Left G source', g_uL
! print*, 'Left G coeff', g_uL
! print*, ''
! print*, 'Right G source', g_uR
! print*, 'Right G coeff', g_uR
! call print_sM('L')
! call print_sM('R')
! call print_sM('I')
! print*, 'Dt', dt, ' - tau', tau
! print*, 'A            =', G_1*r_0
! print*, 'I__u_PSI_g0  =', I__c_PSI_g(1,'I')
! print*, 'B            =', G_2
! print*, 'guL          =', g_uL
! print*, 'MatrixP      =', r_0 * IP__c_PSIxPSI_g(order_2,'I')
! print*, 'guR          =', g_uR
! print*, 'MatrixN      =', IN__c_PSIxPSI_g(order_2,'I')
! print*, 'I__u2_PSI_g0 =', r_0*(MATMUL(IP__c_PSIxPSI_g(order_2,'I'), g_uL) + MATMUL(IN__c_PSIxPSI_g(order_2,'I'), g_uR))
! print*, 'D            =', G_3
! print*, 'A1, A2, A3   =', At
! print*, 'I__u_PSIxPSI =', r_0*(MATMUL(I__c_PSIxPSI_g(order_1,'I'), At))
! print*, ''
! print*, 'A Source term', At
! print*, 'A coeff', At
! print*, c, phi_ij
