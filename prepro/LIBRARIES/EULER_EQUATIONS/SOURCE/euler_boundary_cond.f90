!===============================================================================
!
!      Module: euler_boundary_cond
!
! Description: Procedure for the evaluation of the boundary
!              conditions and boundary numerical fluxes for
!              the Euler equations (implicit/explicit schemes) 
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
!   Copyright: 1998-2008 Alberto Guardone, Marco Fossati
!              See COPYING file for copyright notice
!
!===============================================================================

   MODULE  euler_boundary_cond

   !============================================================================
   USE csr
   USE csr_pair_sys
   
   IMPLICIT NONE
   !============================================================================

   !============================================================================
   TYPE boundary_data
   
      ! Type of boundary condition to be imposed
      INTEGER                               :: cond_type
      ! Type of boundary value 
      INTEGER                               :: value_type
      ! Size of and vector of boundary data 
      INTEGER                               :: size_bdata
      REAL(KIND=8), DIMENSION(:,:), POINTER :: bdata
      ! Number of and list of boundary nodes
      ! (boundary indices)
      INTEGER                               :: Np_bound    
      INTEGER,      DIMENSION(:),   POINTER :: bpoints
      ! Vector of imposed boundary values 
      ! (conservative variables)
      REAL(KIND=8), DIMENSION(:,:), POINTER :: bb
   
   END TYPE boundary_data
   !============================================================================


   ! Boundary condition table
   TYPE(boundary_data), DIMENSION(:), ALLOCATABLE  ::  bc_table
   ! List of boundary points belonging to solid bodies
   INTEGER, DIMENSION(:), ALLOCATABLE  ::  solid_wall_points

   ! Boundary value types                                             
   INTEGER, PARAMETER  ::   BV__NONE            = 0, &
                            BV__WW_OO           = 1, &
                            BV__WW_OO_FFCORRECT = 2, &
                            BV__WW_0            = 3, &
                            BV__CONST_WW        = 4, &
                            BV__CONST_P         = 5, &
                            BV__CONST_RHO_P     = 6, &
                            BV__CONST_RHO_U_P   = 7

   CHARACTER(*), DIMENSION(0:7), PARAMETER  ::  bv_names = &
                        (/ 'NO BOUNDARY VALUE IMPOSED              ', &
                           'WW_OO (CONS. VARIABLE AT INFINITY)     ', &
                           'WW_OO + FAR FIELD CORRECTION           ', &
                           'WW_0 (INITIAL VALUE)                   ', &
                           'CONSTANT CONSERVATIVE VARIABLE WW_C    ', &
                           'CONSTANT PRESSURE                      ', &
                           'CONSTANT DENSITY AND PRESSURE          ', &
                           'CONSTANT DENSITY, VELOCITY AND PRESSURE' /)

   ! Boundary conditions types                                             
   INTEGER, PARAMETER  ::   BC__NONE            = 0, &
                            BC__SOLID_WALL      = 1, &
                            BC__RIEMANN         = 2, &
                            BC__SIMM_PLANE      = 3

   CHARACTER(*), DIMENSION(0:3), PARAMETER  ::  bc_names = &
                        (/ 'NO BOUNDARY CONDITION',  &
                           'SLIP CONDITION       ',  &
                           'RIEMANN              ',  &
                           'SIMMETRY PLANE       '  /)

   ! Variables: private/public policy
   PRIVATE  ::  bc_table, solid_wall_points,                  &
                BV__NONE,               BV__WW_OO,            &
                BV__WW_OO_FFCORRECT,    BV__WW_0,             &
                BV__CONST_WW,           BV__CONST_P,          &
                BV__CONST_RHO_P,        BV__CONST_RHO_U_P,    &
                BC__SOLID_WALL,         BC__RIEMANN,          &
                BC__SIMM_PLANE!,         BC__NON_REFLECT  

   ! Procedures: private/public policy
   PUBLIC   ::  euler_boundary_flux,                          &
                wb__euler_boundary_cond,                      &
                wb__euler_imp_boundary_cond,                  &
                init_euler_boundary_cond,                     &
                read_param_euler_boundary_cond,               &
                write_param_euler_boundary_cond,              &
                euler_type_of_boundary,                       &
                body_pnts, N_body_pnts

   PRIVATE  ::  impose_wall_slip_bc
                !bb__boundary_data,      & 
		!impose_euler_bc,        &
                !impose_free_stream_bc,  &
                !impose_simm_plane_bc
   !============================================================================
   
   CONTAINS


   SUBROUTINE  euler_boundary_flux(wb, xi_bp, jd_jb, rhs)
   
   ! Only nodal contribution (FVM)
   !---------------------------------------------------------------
   USE euler_equations,   ONLY: flux__ww
   
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: wb
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: xi_bp
   INTEGER,      DIMENSION(:),   INTENT(IN) :: jd_jb

   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: rhs
   
   REAL(KIND=8), DIMENSION(SIZE(wb,1), SIZE(xi_bp,1)) :: ff_i
   
   INTEGER :: i, p
   !---------------------------------------------------------------
   
   DO i = 1, SIZE(wb,2)
   
     ff_i = flux__ww(wb(:,i))
   
     DO p = 1, SIZE(wb,1)
       rhs(p,jd_jb(i)) = rhs(p,jd_jb(i)) - SUM(xi_bp(:,i) * ff_i(p,:))
     ENDDO

   ENDDO

   END SUBROUTINE  euler_boundary_flux





   SUBROUTINE  read_param_euler_boundary_cond(idf) 
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
   
        READ (idf,*) bc_table(b) % cond_type,	&
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
    
   END SUBROUTINE  read_param_euler_boundary_cond 


 


   SUBROUTINE  write_param_euler_boundary_cond(idf) 
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf
   INTEGER :: b
   !----------------------------------------------------------------------------

   ! Write the boundary conditions table
   WRITE(idf,*) ' Parameters for euler boundary conditions'
   WRITE(idf,*) ' ----------------------------------------'
   WRITE(idf,*) ' Number of boundaries:', SIZE(bc_table)
   WRITE(idf,*)

   DO b = 1, SIZE(bc_table)

     WRITE(idf,*) '  Boundary', b
     WRITE(idf,*) '	Boundary condition type:  ', &
        	  bc_names(bc_table(b)%cond_type)
     WRITE(idf,*) '	Boundary value type:	  ', &
        	  bv_names(bc_table(b)%value_type)
     WRITE(idf,*)

   ENDDO

   END SUBROUTINE  write_param_euler_boundary_cond 




 
   SUBROUTINE  init_euler_boundary_cond(ww_oo, ww_b, bound_p) 
   !----------------------------------------------------------------------------
   USE dynamic_vector
   USE euler_equations,   ONLY: Mach__ww, mod_u__ww

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww_oo
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b
   INTEGER,	 DIMENSION(:),   INTENT(IN) :: bound_p

   TYPE(D_I_V),  DIMENSION(:),	 ALLOCATABLE :: bound_p_DIV, p_bound_DIV
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: tmp_array
   
   REAL(KIND=8) :: beta, Mach_oo, mod_u_oo, alpha_oo
   
   INTEGER :: i_, i, b, Np_bound, Np_wall
   !----------------------------------------------------------------------------

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
   	   
   	   ALLOCATE (bc_table(b) % bdata(1,1))
   	   bc_table(b) % bdata(:,1) = 0.d0
   	   

   	CASE (BV__WW_OO)
       
   	   ALLOCATE (bc_table(b) % bdata(SIZE(ww_oo),1))
   	   bc_table(b) % bdata(:,1) = ww_oo    


   	CASE (BV__WW_OO_FFCORRECT)

   	   IF (SIZE(ww_oo) .NE. 4) THEN
   	     PRINT*, '';   PRINT*, 'ERROR. INIT_EULER_BOUNDARY_COND:'
   	     PRINT*, 'Far field correction allowed for 2d flows only.'
   	     PRINT*, '';   STOP
   	   ENDIF

   	   IF (bc_table(b) % size_bdata == 0) THEN
   	     ALLOCATE (bc_table(b) % bdata(3,1))
   	     bc_table(b) % bdata(:,1) = (/ 1.d0, 0.25d0, 0.d0 /)
   	   ENDIF
   	   
   	   ALLOCATE (tmp_array(SIZE(bc_table(b) % bdata)))
   	   tmp_array = bc_table(b) % bdata(:,1)
   	   
   	   DEALLOCATE (bc_table(b) % bdata) 
   	   
   	   Mach_oo  =  Mach__ww(ww_oo)
   	   mod_u_oo = mod_u__ww(ww_oo)
   	   
   	   beta = 0.d0
   	   IF (Mach_oo < 1.d0) THEN
   	      beta = SQRT(1.d0 - Mach_oo**2)
   	   ENDIF

   	   alpha_oo = ATAN2(ww_oo(3),ww_oo(2))

   	   ALLOCATE (bc_table(b) % bdata(SIZE(ww_oo),3))
   	   bc_table(b) % bdata(:,1) = ww_oo
   	   bc_table(b) % bdata(:,2) = (/ tmp_array, beta /) 
   	   bc_table(b) % bdata(:,3) = (/ Mach_oo, mod_u_oo, alpha_oo, 0.d0/)
   	   
   	   DEALLOCATE (tmp_array)	     
   	   

   	CASE (BV__WW_0)
   		    
   	   Np_bound = SIZE(bc_table(b) % bpoints)
   	   ALLOCATE (bc_table(b) % bdata(SIZE(ww_oo), Np_bound))
   	   
   	   DO i_ = 1, Np_bound
   	     i = bc_table(b) % bpoints(i_)
   	     bc_table(b) % bdata(:,i_) = ww_b(:,i)
   	   ENDDO	  


   	CASE(BV__CONST_WW)
   	! Read in read_param_euler_boundary_cond
   	
        CASE(BV__CONST_P)
   	! Read in read_param_euler_boundary_cond
   	
        CASE(BV__CONST_RHO_P)
   	! Read in read_param_euler_boundary_cond
   	
        CASE(BV__CONST_RHO_U_P)
   	! Read in read_param_euler_boundary_cond
        
        CASE DEFAULT

   	   PRINT*, '';   PRINT*, 'ERROR. READ_BC_DATA'
   	   PRINT*, 'Unknown boundary value.';	PRINT*, ''
   	   STOP       
      
      END SELECT
   
   ENDDO

   
   
   ! The number of wall points is computed
   Np_wall = 0
   DO b = 1, SIZE(bc_table)
   
     IF (bc_table(b) % cond_type == BC__SOLID_WALL) &
     Np_wall = Np_wall + SIZE(bc_table(b) % bpoints)
      
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

   END SUBROUTINE  init_euler_boundary_cond 





   FUNCTION  wb__euler_boundary_cond(ww_b, normal, rr_b)   RESULT(wb)
   !------------------------------------------------------------ 
   USE euler_equations

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: normal
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr_b

   REAL(KIND=8), DIMENSION(SIZE(ww_b,1),SIZE(ww_b,2)) :: wb
   
   REAL(KIND=8), DIMENSION(:,:), POINTER :: bdata, bb
   INTEGER,	 DIMENSION(:),   POINTER :: bpoints

   INTEGER :: cond_type, value_type
   INTEGER :: b
   !------------------------------------------------------------ 
   
   ! Computation of the imposed values (bb) and
   ! of the conservative variable (wb) after
   ! the boundary condition have been imposed
   DO b = 1, SIZE(bc_table), 1

! print*, ''
! print*, ''
! print*, b
    
      cond_type  =  bc_table(b) % cond_type
      value_type =  bc_table(b) % value_type
      bpoints	 => bc_table(b) % bpoints
      bb	 => bc_table(b) % bb
      bdata	 => bc_table(b) % bdata

      bb = imposed_boundary_values(ww_b(:,bpoints),   &
   				   rr_b(:,bpoints),   &
   				   normal(:,bpoints), &
   				   bdata, value_type, & 
   				   ww_b(:,solid_wall_points), & 
   				   normal(:,solid_wall_points))

      wb(:,bpoints) = impose_euler_boundary_cond(bb, ww_b(:,bpoints),  &
                                                 normal(:,bpoints), cond_type)
   ENDDO


   END FUNCTION  wb__euler_boundary_cond





   FUNCTION  wb__euler_imp_boundary_cond(ww_b, normal, rr_b, dwb_dwi)   RESULT(wb)
   !------------------------------------------------------------------------------
   USE euler_equations

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: normal
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr_b

   REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: dwb_dwi
   
   REAL(KIND=8), DIMENSION(SIZE(ww_b,1), SIZE(ww_b,2)) :: wb

   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Dwb
   REAL(KIND=8), DIMENSION(:,:),   POINTER     :: bdata, bb
   INTEGER,	 DIMENSION(:),     POINTER     :: bpoints
   
   INTEGER :: cond_type, value_type
   INTEGER :: b, i
   !------------------------------------------------------------------------------ 
   
   ! Computation of the imposed values (bb) and
   ! of the conservative variable (wb) after
   ! the boundary condition have been imposed
   DO b = 1, SIZE(bc_table)
 

      cond_type  =  bc_table(b) % cond_type
      value_type =  bc_table(b) % value_type
      bpoints	 => bc_table(b) % bpoints
      bb	 => bc_table(b) % bb
      bdata	 => bc_table(b) % bdata


      bb = imposed_boundary_values(ww_b(:,bpoints), rr_b(:,bpoints),     &
   				   normal(:,bpoints), bdata, value_type, &
   				   ww_b(:,solid_wall_points),            & 
   				   normal(:,solid_wall_points))

      ALLOCATE (Dwb(SIZE(ww_b,1), SIZE(ww_b,1), SIZE(bpoints))) 
 
      wb(:,bpoints) = impose_euler_imp_boundary_cond(bb, ww_b(:,bpoints), &
   						     normal(:,bpoints),	  &
   						     cond_type, Dwb)
       
      DO i = 1, SIZE(bpoints)
   	dwb_dwi(:,:,bpoints(i)) = Dwb(:,:,i)
      ENDDO
      
      DEALLOCATE (Dwb) 

   ENDDO

   END FUNCTION  wb__euler_imp_boundary_cond





   FUNCTION  imposed_boundary_values(ww_b, rr_b, normal, bdata, &		      
   				     value_type, wwbody, normalbody)   RESULT (bb)   

   ! Computes the value of the conserved variables corresponding to the imposed       
   ! boundary data								      
   !--------------------------------------------------------------------------------- 
   USE euler_equations  							      
   										      
   IMPLICIT NONE								      

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b				    
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr_b				    
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: normal				    
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: bdata				    
   INTEGER		       , INTENT(IN) :: value_type			    
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: wwbody, normalbody		    

   REAL(KIND=8), DIMENSION(SIZE(ww_b,1), SIZE(ww_b,2)) :: bb

   ! Variables needed for primitive variables manipulation			      
   REAL(KIND=8)  ::  rho, pres  						      
   REAL(KIND=8), DIMENSION(SIZE(ww_b,1) - 2)  ::  uu				      

   ! Variables needed for the vortex correction 				      
   REAL(KIND=8), DIMENSION(SIZE(normal,1))  ::  r, X_cp 			      
   REAL(KIND=8), DIMENSION(SIZE(ww_b,1))    ::  ww_oo				      
   REAL(KIND=8)  ::  mod_u_oo, Mach_oo, alpha_oo, &				      
   		     theta, F, beta, chord, CL, dr				      
   REAL(KIND=8)  ::  pi ! pi greek						      

   INTEGER  ::  i, p								      
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
   										      
   										      
      ! Conservative variables at infinity with vortex correction  						      
      ! (2D subsonic flows only)
      CASE(BV__WW_OO_FFCORRECT) 						      

   	 ww_oo    = bdata(:,  1)  ! Reference conserv. variables		      
   	 chord    = bdata(1,  2)  ! Geom. quantities: choord			      
   	 X_cp	  = bdata(2:3,2)  !		      X_cp			      
   	 beta	  = bdata(4,  2)  ! Flow quantities:  Beta			      
   	 Mach_oo  = bdata(1,  3)  !		      Mach_oo			      
   	 mod_u_oo = bdata(2,  3)  !		      Module of the velocity	      
   	 alpha_oo = bdata(3,  3)  !		      flow direction		      
   	 pi = 4.d0 * ATAN(1.d0)   ! Pi greek					      
   										      
   	 CL = 1.d0 !CL_body( wwbody, ww_oo, normalbody, chord )   ! Lift coefficient        

   	 ! Computation of the vortex correction to the  			      
   	 ! reference value ww_oo for each node  				      
   	 DO i = 1, SIZE(bb,2)							      

   	    r	  = rr_b(:,i) - X_cp						      
   	    dr    = SQRT( SUM(r**2) )						      
   	    theta = ATAN2(r(2),r(1))						      

   	    F = (0.25*CL*chord) * (beta/(pi*dr)) &				      
   	      * (1/(1 - (Mach_oo * SIN(theta - alpha_oo))**2))  		      

   	    bb(1,i) = ww_oo(1)  						      
   	    bb(2,i) = ww_oo(1) * mod_u_oo * (COS(alpha_oo) + F*SIN(theta))	      
   	    bb(3,i) = ww_oo(1) * mod_u_oo * (SIN(alpha_oo) - F*COS(theta))	      
   	    bb(4,i) = ww_oo(4)  						      

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





   FUNCTION  impose_euler_boundary_cond(bb, ww_b, normal, cond_type)   RESULT(wb)
   
   ! Explicit time integration scheme
   !---------------------------------------------------------------------------------
   use mp_interface
   
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: bb, ww_b
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: normal
   INTEGER,			 INTENT(IN) :: cond_type   

   REAL(KIND=8), DIMENSION(SIZE(ww_b,1), SIZE(ww_b,2)) :: wb

   INTEGER  ::  i
   !---------------------------------------------------------------------------------
   
   SELECT CASE (cond_type)

      ! No boundary condition ---> wb == bb
      CASE (BC__NONE)							  
   									  
   	 DO i = 1, SIZE(ww_b,2)
   	   wb(:,i) = bb(:,i)
   	 ENDDO


      ! Solid Wall
      CASE (BC__SOLID_WALL)						  
   									  
   	 DO i = 1, SIZE(ww_b,2)	 
   	   wb(:,i) = impose_wall_slip_bc(ww_b(:,i), normal(:,i))
         ENDDO

      ! Riemann-type boundary conditions
      CASE (BC__RIEMANN)

   	 DO i = 1, SIZE(ww_b,2)
   	   wb(:,i) = impose_riemann_bc(ww_b(:,i), bb(:,i), normal(:,i))
   	 ENDDO


      ! Simmetry plane
      CASE (BC__SIMM_PLANE)

   	 DO i = 1, SIZE(ww_b,2)
   	   wb(:,i) = impose_wall_slip_bc(ww_b(:,i), normal(:,i))
   	 ENDDO

   END SELECT
   
   END FUNCTION  impose_euler_boundary_cond





   FUNCTION  impose_euler_imp_boundary_cond(bb, ww_b, normal, cond_type, Dwb)   RESULT(wb)
   
   ! Implicit time integration scheme
   !-------------------------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  bb, ww_b
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  normal
   INTEGER,			 INTENT(IN)  ::  cond_type   
   REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)  :: Dwb 
   
   REAL(KIND=8), DIMENSION(SIZE(ww_b,1), SIZE(ww_b,2))  ::  wb

   INTEGER  ::  i, p
   !-------------------------------------------------------------------------------------------

      SELECT CASE (cond_type)

         ! No boundary condition ---> wb == bb
         CASE (BC__NONE)
         
            DO i = 1, SIZE(ww_b,2)
     
               wb(:,i) = bb(:,i)

               Dwb = 0.d0
               DO p = 1, SIZE(ww_b,1)
                 Dwb(p,p,i) = 1.d0
               ENDDO

            ENDDO


         ! Solid wall
         CASE (BC__SOLID_WALL)
         
            DO i = 1, SIZE(ww_b,2)
              wb(:,i) = impose_imp_wall_slip_bc(ww_b(:,i), normal(:,i), Dwb(:,:,i))
            ENDDO


         ! Riemann-type boundary conditions
         CASE (BC__RIEMANN)

            DO i = 1, SIZE(ww_b,2)
              wb(:,i) = impose_imp_riemann_bc(ww_b(:,i), bb(:,i), normal(:,i), Dwb(:,:,i))
            ENDDO


         ! Simmetry plane
         CASE (BC__SIMM_PLANE)
 
            DO i = 1, SIZE(ww_b,2)
              wb(:,i) = impose_imp_wall_slip_bc(ww_b(:,i), normal(:,i), Dwb(:,:,i))
            ENDDO

      END SELECT

   END FUNCTION  impose_euler_imp_boundary_cond





   !============================================================
   !  Procedures for the computation of the conservative
   !  variables for each type of boundary condition:
   !============================================================

   FUNCTION  impose_wall_slip_bc(ww_i, normal)   RESULT(ww_)

   ! Slip condition (null normal velocity) 
   ! [ uu_ = uu - uu_normal ]
   !------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ww_i
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: normal

   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: ww_
   REAL(KIND=8), DIMENSION(SIZE(normal)) :: nn

   INTEGER :: p
   !------------------------------------------------------------ 

   ! Slip condition (null normal velocity) 
   ! is imposed [ uu_ = uu - uu_normal ]
   ! ------------------------------------- 
   p = SIZE(ww_i)

   nn = normal / SQRT(SUM(normal**2))
   
   ww_(1)     = ww_i(1)
   ww_(2:p-1) = ww_i(2:p-1) - SUM(ww_i(2:p-1)*nn)*nn
   ww_(p)     = ww_i(p) - 0.5*(SUM(ww_i(2:p-1)**2) / ww_i(1)) &
   			+ 0.5*(SUM(ww_ (2:p-1)**2) / ww_ (1)) 

   END FUNCTION  impose_wall_slip_bc





   FUNCTION  impose_riemann_bc(ww_i, ww_out, normal)   RESULT(ww_)

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
   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: lambda, ww
   REAL(KIND=8), DIMENSION(SIZE(normal)) :: nn
   !------------------------------------------------------------ 

   nn = normal / SQRT(SUM(normal**2))

   ww = ww_i  ! <<< State used for the linearization

! print*, 'boundary value', ww_i
! print*, 'imposed value', ww_out
! print*, 'pivot state', ww

   ! Eigenstructure of the system linearized in ww 
       RR = right_eigenvectors__ww_nn(ww, nn)
       LL =  left_eigenvectors__ww_nn(ww, nn)
   lambda =	   eigenvalues__ww_nn(ww, nn)


   ! Jump between the initial values of the characteristic 
   ! variables. alpha(p) = L^p Dw
   alpha = MATMUL(LL, ww_out - ww_i)

!print*, 'alpha', alpha

   !------------------------------------------------------------ 
   !		      /
   !		      | alpha(p)  if lambda(p) < 0
   ! alpha^Neg (p) == |
   !		      |   0	  if lambda(p) > 0
   !		      \
   !------------------------------------------------------------
   
   WHERE (lambda < 0)  
     alpha_Neg = alpha
   ELSEWHERE
     alpha_Neg = 0.d0
   END WHERE

!print*, 'alphaNeg', alpha_Neg

   !------------------------------------------------------------
   ! Reconstruction of the (vector) unknown ww_
   !
   !   ww_  =	ww_i  +  SUM_p alpha^Neg (p)  RR(:,p)
   !
   ! alternatively (not used here),
   !
   !   ww_  = ww_out  -  SUM_p alpha^Pos (p)  RR(:,p)
   !------------------------------------------------------------
   
   ww_ = ww_i + MATMUL(RR, alpha_Neg)

! print*, 'imposed BC', ww_
! print*, ''

   END FUNCTION  impose_riemann_bc





   !=====================================================================
   !  Procedures for the computation of the conservative
   !  variables for each type of boundary condition: IMPLICIT SCHEME
   !=====================================================================

   FUNCTION  impose_imp_wall_slip_bc(ww_i, normal, dwb_dwi)   RESULT(ww_)

   ! Slip condition (null normal velocity) 
   ! [ uu_ = uu - uu_normal ]
   !---------------------------------------------------------------------
   USE euler_equations

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ww_i
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: normal
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: dwb_dwi

   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: ww_
   REAL(KIND=8), DIMENSION(SIZE(normal)) :: nn
   
   REAL(KIND=8) :: m_n
   
   INTEGER :: p, k, l
   !--------------------------------------------------------------------- 

   ! Slip condition (null normal velocity) 
   ! is imposed [ uu_ = uu - uu_normal ]
   p = SIZE(ww_i)

   nn = normal / SQRT(SUM(normal**2))
   m_n = SUM(ww_i(2:p-1)*nn)   
   
   ww_(1)     = ww_i(1)
   ww_(2:p-1) = ww_i(2:p-1) - m_n*nn
   ww_(p)     = ww_i(p) - 0.5*(SUM(ww_i(2:p-1)**2) / ww_i(1)) &
   			+ 0.5*(SUM(ww_ (2:p-1)**2) / ww_ (1)) 

   dwb_dwi = 0.d0
   
   ! Mass
   dwb_dwi(1,1) = 1.d0         
   
   ! Momentum
   DO k = 1, SIZE(ww_i) - 2
      dwb_dwi(k+1,k+1) = 1.d0
      DO l = 1, SIZE(ww_i) - 2
      dwb_dwi(k+1,l+1) =  dwb_dwi(k+1,l+1) - nn(k)*nn(l)
      ENDDO
   ENDDO

   ! Total energy
   dwb_dwi(p,1)     = 0.5*(m_n/ww_i(1))**2
   dwb_dwi(p,2:p-1) = - (m_n/ww_i(1))*nn
   dwb_dwi(p,p)     = 1.d0

   END FUNCTION  impose_imp_wall_slip_bc





   FUNCTION  impose_imp_riemann_bc(ww_i, ww_out, normal, dwb_dwi)   RESULT(ww_)

   ! Imposes Riemann boundary conditions. 
   ! Characteristic reconstruction of the unknown
   ! The system is linearized around the state ww_i
   !----------------------------------------------------------------------------
   USE euler_equations

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ww_i, ww_out
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: normal
   REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: dwb_dwi

   REAL(KIND=8), DIMENSION(SIZE(ww_i)) :: ww_

   REAL(KIND=8), DIMENSION(SIZE(ww_i),SIZE(ww_i)) :: RR, LL
   REAL(KIND=8), DIMENSION(SIZE(ww_i),SIZE(ww_i)) :: Neg
   
   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: alpha, alpha_Neg
   REAL(KIND=8), DIMENSION(SIZE(ww_i))   :: lambda, ww
   REAL(KIND=8), DIMENSION(SIZE(normal)) :: nn

   INTEGER  ::  p
   !----------------------------------------------------------------------------
   
   nn = normal / SQRT(SUM(normal**2))

   ww = ww_i  ! <<< State used for the linearization

   ! Eigenstructure of the system linearized in ww 
       RR = right_eigenvectors__ww_nn(ww, nn)
       LL =  left_eigenvectors__ww_nn(ww, nn)
   lambda =	   eigenvalues__ww_nn(ww, nn)


   ! alpha(p) = L^p Dw
   alpha = MATMUL(LL, ww_out - ww_i)


   !------------------------------------------------------------ 
   !		      /
   !		      | alpha(p)  if lambda(p) < 0
   ! alpha^Neg (p) == |
   !		      |   0	  if lambda(p) > 0
   !		      \
   !------------------------------------------------------------ 

   WHERE (lambda < 0)  
     alpha_Neg = alpha
   ELSEWHERE
     alpha_Neg = 0.d0
   END WHERE


   !------------------------------------------------------------ 
   ! Reconstruction of the (vector) unknown ww_
   !
   !   ww_  =	ww_i  +  SUM_p alpha^Neg (p)  RR(:,p)
   !
   ! alternatively (not used here),
   !
   !   ww_  = ww_out  -  SUM_p alpha^Pos (p)  RR(:,p)
   !
   !------------------------------------------------------------

   ww_ = ww_i + MATMUL(RR, alpha_Neg)
   
   Neg = 0.d0
   DO  p = 1, SIZE(ww_i)
     Neg(p,p) = - MIN(0.d0, SIGN(1.d0,lambda(p)))
   ENDDO
   
   dwb_dwi = - MATMUL(RR, MATMUL(Neg,LL))
   DO  p = 1, SIZE(ww_i)
     dwb_dwi(p,p) = dwb_dwi(p,p) + 1.d0
   ENDDO

   END FUNCTION  impose_imp_riemann_bc





   !============================================================
   !  Procedures for the computation of integral quantities
   !  related to the body:
   !
   !	 -  Lift Coefficient CL and drag coefficient CD
   !
   !	 -  List of nodes belonging to the body(ies)
   !============================================================

   FUNCTION  euler_type_of_boundary(bound)   RESULT(btype)
   !------------------------------------------------------------ 
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: bound
   INTEGER :: btype
   !------------------------------------------------------------ 

   btype = bc_table(bound) % cond_type
   
   END FUNCTION  euler_type_of_boundary
 




   FUNCTION  body_pnts()   RESULT(jb_jbody)
   !------------------------------------------------------------ 
   IMPLICIT NONE
   INTEGER, DIMENSION(:), POINTER :: jb_jbody
   !------------------------------------------------------------ 

   ALLOCATE (jb_jbody(SIZE(solid_wall_points)))

   jb_jbody = solid_wall_points

   END FUNCTION  body_pnts





   FUNCTION  N_body_pnts()  RESULT(N_jbody)
   !------------------------------------------------------------
   IMPLICIT NONE
   INTEGER :: N_jbody
   !------------------------------------------------------------ 

   N_jbody = SIZE(solid_wall_points)
      
   END FUNCTION  N_body_pnts
   

   END MODULE  euler_boundary_cond


