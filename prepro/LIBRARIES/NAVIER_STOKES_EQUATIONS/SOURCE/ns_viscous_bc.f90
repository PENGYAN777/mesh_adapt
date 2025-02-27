!===============================================================================
!
!      Module: ns_viscous_bc
!
! Description: Procedure for the evaluation of the boundary
!              conditions and boundary numerical fluxes for
!              the velocity in the Navier-Stokes equations 
!              (explicit schemes only) 
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
!===============================================================================

   MODULE  ns_viscous_bc

   !============================================================================
   USE dynamic_vector
   
   IMPLICIT NONE
   !============================================================================

   !============================================================================
   TYPE boundary_data
   
      ! Type of boundary condition (weak/strong) 
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
  
   END TYPE boundary_data
   !============================================================================


   PRIVATE
   
   ! Boundary condition table
   TYPE(boundary_data), DIMENSION(:), ALLOCATABLE  ::  bc_table
   ! List of boundary points belonging to solid bodies
   INTEGER, DIMENSION(:), ALLOCATABLE  ::  solid_wall_points
   ! Constants: 
   !   Boundary value types                        
   INTEGER, PARAMETER  ::  BV__VV_NONE         =  0, &
                           BV__VV_NOSLIP       =  1, &
                           BV__VV_OO           =  2, &
                           BV__VV_0            =  3, &
                           BV__VV_CONST        =  4, &
                           BV__VV_SIMM         =  5, &
                           BV__VV_SLIP         =  6, &
                           BV__GGVV_NONE       = 10, &
                           BV__GGVV_NEUM_HOMOG = 11

   CHARACTER(*), DIMENSION(0:11), PARAMETER  ::  bv_names = &
                        (/ 'VELOCITY: NO BOUNDARY VALUE IMPOSED    ', &
                           'VELOCITY: NO SLIP CONDITION            ', &
                           'VELOCITY: VV_OO (VELOCITY AT INFINITY) ', &
                           'VELOCITY: VV_0 (INITIAL VALUE)         ', &
                           'VELOCITY: CONSTANT VELOCITY            ', &
                           'VELOCITY: SIMMETRY PLANE               ', &
                           'VELOCITY: SLIP CONDITION (GOKCEN)      ', &
                           'UNKNOWN TYPE                           ', &
                           'UNKNOWN TYPE                           ', &
                           'UNKNOWN TYPE                           ', &
                           'GG VV: NONE                            ', &
                           'GG VV: HOMOGENEOUS                     ' /)

   ! Boundary conditions types                                             
   INTEGER, PARAMETER  ::   BC__WEAK   = 0, &
                            BC__STRONG = 1

   CHARACTER(*), DIMENSION(0:1), PARAMETER  ::  bc_names = &
                        (/ 'WEAK FORM  ',  &
                           'STRONG FORM' /)

   REAL(KIND=8) :: C_visc  
		    
   PUBLIC   ::  ns_viscous_boundary_flux,  &
                vb__ns_viscous_bc,         &
		GG_vb__ns_viscous_bc,      &
                ns_impose_strong_vv_bc,    &
                      init_ns_viscous_bc,  &
                read_param_ns_viscous_bc,  &
               write_param_ns_viscous_bc
   !============================================================================

   CONTAINS
 
 
   SUBROUTINE  ns_viscous_boundary_flux(vv, vb, GG_vb, mu, ll, rhs)
   !------------------------------------------------------------------------------------------------
   USE np_quadrature
   USE node_pair_structure,  ONLY: j_c_d, j_c_b, jd_jb
   USE metric_coefficients,  ONLY: eta, cell, chi_b, xi_bp

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: vv
   REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: vb
   REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: GG_vb
   REAL(KIND=8), DIMENSION(:),     INTENT(IN) :: mu, ll
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: rhs
   
   REAL(KIND=8), DIMENSION(SIZE(GG_vb,1), SIZE(GG_vb,2), SIZE(GG_vb,3)) :: mu_GG_vb
   REAL(KIND=8), DIMENSION(SIZE(GG_vb,1), SIZE(GG_vb,2), SIZE(GG_vb,3)) :: mu_GG_vbt

   REAL(KIND=8), DIMENSION(SIZE(vb,1), SIZE(vb,2)) :: mu_GG_vb_vb, lvDvb   
   REAL(KIND=8), DIMENSION(SIZE(vv,1), SIZE(vv,2)) :: mu_G_ke
   
   REAL(KIND=8), DIMENSION(SIZE(vv,2)) :: k_e
   REAL(KIND=8), DIMENSION(SIZE(vb,2)) :: k_eb 
   REAL(KIND=8), DIMENSION(SIZE(vv,2)) :: D_vv
   REAL(KIND=8), DIMENSION(SIZE(vv,1)) :: ff
   REAL(KIND=8), DIMENSION(SIZE(vv,1)) :: vv_ij, pp_ij
   
   REAL(KIND=8) :: uu_ij, ss
   
   INTEGER :: i, j, c, i_, j_, k, k_d
   !------------------------------------------------------------------------------------------------
   
   k_d = SIZE(rhs,1) - 2
   
   ! Divergence, D_vv, of the velocity vector vv
   D_vv = (1/cell) * np_quadr_w_Dv(j_c_d, j_c_b, jd_jb, &
                                   eta, chi_b, xi_bp, vv, vb)
   
   DO i = 1, SIZE(vb,2)
     mu_GG_vb (:,:,i)  =  mu(jd_jb(i)) * GG_vb(:,:,i)	   
     mu_GG_vbt(:,:,i)  =  mu(jd_jb(i)) * TRANSPOSE(GG_vb(:,:,i))   
   ENDDO
   
      
   !============================================================================
   ! MOMENTUM EQUATION
   !============================================================================
   
   ! DIV [mu GG_vv] ------------------------------------------------------------
   
   ! Node-pair contribution
   DO c = 1, SIZE(chi_b,2)

     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)
    
     DO k = 1, k_d
       ff(k) = SUM(chi_b(:,c) * (mu_GG_vb(k,:,j_) - mu_GG_vb(k,:,i_)))
     ENDDO

     rhs(2:k_d+1,i) = rhs(2:k_d+1,i) + ff * C_visc
     rhs(2:k_d+1,j) = rhs(2:k_d+1,j) - ff * C_visc

   ENDDO
   
   ! Node contribution
   DO i_ = 1, SIZE(xi_bp,2) 
   
     i = jd_jb(i_)
     
     DO k = 1, k_d
       ff(k) = SUM(xi_bp(:,i_) * mu_GG_vb(k,:,i_))
     ENDDO

     rhs(2:k_d+1,i) = rhs(2:k_d+1,i) + ff * C_visc 

   ENDDO

       
   ! GRAD [(mu + lambda) DIV V] ------------------------------------------------

   ! Boundary condition included in the computation of the divergence D_vv.
   ! Node-pair contribution
   DO c = 1, SIZE(chi_b,2)

     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)

     pp_ij = chi_b(:,c) * 0.5*((mu(j) + ll(j))*D_vv(j) - (mu(i) + ll(i))*D_vv(i))
   
     rhs(2:k_d+1,i) = rhs(2:k_d+1,i) + pp_ij * C_visc
     rhs(2:k_d+1,j) = rhs(2:k_d+1,j) - pp_ij * C_visc
     
   ENDDO

   ! Node contribution
   DO i_ = 1, SIZE(xi_bp,2)

     i = jd_jb(i_)
     rhs(2:k_d+1,i) = rhs(2:k_d+1,i)  +  (mu(i) + ll(i))*D_vv(i) * xi_bp(:,i_) * C_visc 

   ENDDO



   ! DIV [mu (GRAD V)^T] -------------------------------------------------------
   
   DO k = 1, SIZE(vv,1)
   
     ! Node-pair contribution
     DO c = 1, SIZE(chi_b,2)

       i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
       i  = jd_jb(i_);   j  = jd_jb(j_)

       vv_ij = 0.5 * (mu_GG_vbt(k,:,j_)  -  mu_GG_vbt(k,:,i_))
    
       uu_ij = SUM(vv_ij * chi_b(:,c))

       rhs(k+1,i) = rhs(k+1,i) + uu_ij * C_visc
       rhs(k+1,j) = rhs(k+1,j) - uu_ij * C_visc

     ENDDO

     ! Nodal contribution
     DO i_ = 1, SIZE(xi_bp,2)

       i = jd_jb(i_)
       rhs(k+1,i) = rhs(k+1,i) + SUM(mu_GG_vbt(k,:,i_) * xi_bp(:,i_)) * C_visc

     ENDDO

   ENDDO


      
   ! - GRAD [mu DIV V] ---------------------------------------------------------

   ! Boundary condition included in the computation of the divergence D_vv.
   ! Node-pair contribution
   DO c = 1, SIZE(chi_b,2)

     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)

     pp_ij = - chi_b(:,c) * 0.5*(mu(j)*D_vv(j) - mu(i)*D_vv(i))
   
     rhs(2:k_d+1,i) = rhs(2:k_d+1,i) + pp_ij * C_visc
     rhs(2:k_d+1,j) = rhs(2:k_d+1,j) - pp_ij * C_visc
     
   ENDDO

   ! Node contribution
   DO i_ = 1, SIZE(xi_bp,2)

     i = jd_jb(i_)
     rhs(2:k_d+1,i) = rhs(2:k_d+1,i)  -  mu(i)*D_vv(i) * xi_bp(:,i_) * C_visc 

   ENDDO


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

   ! Node-pair contribution
   DO c = 1, SIZE(chi_b,2) 

     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)

     ss = SUM(chi_b(:,c) * (mu_G_ke(:,j) - mu_G_ke(:,i)))  
    
     rhs(k_d+2,i) = rhs(k_d+2,i) + ss * C_visc
     rhs(k_d+2,j) = rhs(k_d+2,j) - ss * C_visc

   ENDDO
   
   ! Node contributions
   DO i_ = 1, SIZE(xi_bp,2)

     i = jd_jb(i_)
     
     ss = SUM(xi_bp(:,i_) * mu_G_ke(:,i))
     rhs(k_d+2,i) = rhs(k_d+2,i) + ss * C_visc

   ENDDO
   

   
   ! DIV [(mu + lambda) (DIV V) V] ---------------------------------------------
   
   DO k = 1, SIZE(vb,1)
     lvDvb(k,:) = (mu(jd_jb) + ll(jd_jb)) * vb(k,:) * D_vv(jd_jb)
   ENDDO    
   
   ! Node-pair contribution
   DO c = 1, SIZE(chi_b,2)

     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)

     vv_ij = 0.5 * (lvDvb(:,j_) - lvDvb(:,i_))
    
     uu_ij = SUM(vv_ij * chi_b(:,c))

     rhs(k_d+2,i) = rhs(k_d+2,i) + uu_ij * C_visc
     rhs(k_d+2,j) = rhs(k_d+2,j) - uu_ij * C_visc

   ENDDO

   ! Nodal contribution
   DO i_ = 1, SIZE(xi_bp,2)

     i = jd_jb(i_)
     rhs(k_d+2,i) = rhs(k_d+2,i) + SUM(lvDvb(:,i_) * xi_bp(:,i_)) * C_visc

   ENDDO



   ! DIV [mu(GRAD V)^T * V  -  mu(DIV V) V]] -----------------------------------	 
	 
   DO i_ = 1, SIZE(vb,2)
   
     i = jd_jb(i_)
     DO k = 1, SIZE(vb,1)
       mu_GG_vb_vb(k,i_) = SUM(vb(:,i_) * mu_GG_vb(k,:,i_)) - mu(i) * D_vv(i) * vb(k,i_) 	   
     ENDDO
     	   
   ENDDO	 

   ! Node-pair contribution
   DO c = 1, SIZE(chi_b,2)

     i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
     i  = jd_jb(i_);   j  = jd_jb(j_)

     vv_ij = 0.5 * (mu_GG_vb_vb(:,j_)  -  mu_GG_vb_vb(:,i_))
    
     uu_ij = SUM(vv_ij * chi_b(:,c))

     rhs(k_d+2,i) = rhs(k_d+2,i) + uu_ij * C_visc
     rhs(k_d+2,j) = rhs(k_d+2,j) - uu_ij * C_visc

   ENDDO

   ! Nodal contribution
   DO i_ = 1, SIZE(xi_bp,2)

     i = jd_jb(i_)
     rhs(k_d+2,i) = rhs(k_d+2,i) + SUM(mu_GG_vb_vb(:,i_) * xi_bp(:,i_)) * C_visc

   ENDDO   
   
   END SUBROUTINE  ns_viscous_boundary_flux
   
   
   
   
   
   SUBROUTINE  read_param_ns_viscous_bc(idf)
   !------------------------------------------------------------ 
   USE mp_interface
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf
   INTEGER  :: b, N_bound, dType, dVType, dSize
   !------------------------------------------------------------ 

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

   END SUBROUTINE  read_param_ns_viscous_bc 

    
 


   SUBROUTINE write_param_ns_viscous_bc (idf) 
   !------------------------------------------------------------ 
   IMPLICIT NONE

   INTEGER,		    INTENT(IN)  ::  idf
   INTEGER  ::  b
   !------------------------------------------------------------ 

   ! Write the boundary conditions table
   WRITE(idf,*) '   Parameters for ns viscous boundary conditions'
   WRITE(idf,*) ' -----------------------------------------------'
   WRITE(idf,*) '   Number of boundaries:', SIZE(bc_table)
   WRITE(idf,*)

   DO b = 1, SIZE(bc_table)

      WRITE(idf,*) '  Boundary', b
      WRITE(idf,*) '     Boundary value type:	 ', &
   		   bv_names(bc_table(b)%value_type)
      WRITE(idf,*) '     Boundary condition is imposed in a: ', &
   		   bc_names(bc_table(b)%cond_type)
      WRITE(idf,*)

   ENDDO

   END SUBROUTINE  write_param_ns_viscous_bc 




   SUBROUTINE  init_ns_viscous_bc(ww_oo, ww_b, bound_p) 
   !-------------------------------------------------------------------------
   USE dynamic_vector
   USE structures,   ONLY: REF

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww_oo
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww_b
   INTEGER,      DIMENSION(:),   INTENT(IN) :: bound_p

   TYPE(D_I_V), DIMENSION(:), ALLOCATABLE :: bound_p_DIV, p_bound_DIV
						  
   INTEGER :: i_, i, b, Np_bound,  &
              value_type, Np_wall, k_d
   !-------------------------------------------------------------------------

   k_d = SIZE(ww_oo)-2

   ! Boundary to node connectivity
   ALLOCATE (bound_p_DIV(SIZE(bound_p)))
   DO i = 1, SIZE(bound_p)
     ALLOCATE (bound_p_DIV(i) % vec(1))
     bound_p_DIV(i) % vec(1) = bound_p(i)
   ENDDO

   ALLOCATE (p_bound_DIV(size_DIV(bound_p_DIV,3)))
   p_bound_DIV = invert_DIV(bound_p_DIV)
   
   DO b = 1, SIZE(bc_table)

     bc_table(b) % Np_bound = SIZE(p_bound_DIV(b)% vec)

     ALLOCATE (bc_table(b) % bpoints(bc_table(b) % Np_bound))
     bc_table(b) % bpoints = p_bound_DIV(b) % vec

   ENDDO
   
   DEALLOCATE (bound_p_DIV, p_bound_DIV)


   ! Initializes the boundary condition values
   DO b = 1, SIZE(bc_table)

      Np_bound   = bc_table(b) % Np_bound
      value_type = bc_table(b) % value_type
      
      SELECT CASE (value_type)            

         ! Velocity boundary values --------------------------------------------
         CASE (BV__VV_NONE)
	 
           ALLOCATE (bc_table(b) % bdata(1,1))
           bc_table(b) % bdata(1,1) = 0.d0


         CASE (BV__VV_NOSLIP)

           ALLOCATE (bc_table(b) % bdata(1,1))
           bc_table(b) % bdata(1,1) = 0.d0


         CASE (BV__VV_OO)
    
           ALLOCATE (bc_table(b) % bdata(k_d,1))
           bc_table(b) % bdata(:,1) = ww_oo(2:k_d+1)/ww_oo(1)


         CASE (BV__VV_CONST)
         ! Read in read_param_ns_boundary_cond


         CASE (BV__VV_0)

           ALLOCATE (bc_table(b) % bdata(k_d, Np_bound))
           DO i_ = 1, Np_bound
             i = bc_table(b) % bpoints(i_)
             bc_table(b) % bdata(:,i_) = ww_b(2:k_d+1,i)/ww_b(1,i)
           ENDDO


         CASE(BV__VV_SIMM)

           ALLOCATE (bc_table(b) % bdata(1,1))
           bc_table(b) % bdata(1,1) = 0.d0



         CASE (BV__VV_SLIP)

           ALLOCATE (bc_table(b) % bdata(1,1))
           bc_table(b) % bdata(1,1) = 0.d0
	 !----------------------------------------------------------------------  


         ! Normal velocity gradient boundary values ----------------------------
         CASE (BV__GGVV_NONE)

           ALLOCATE (bc_table(b) % bdata(1,1))
           bc_table(b) % bdata(1,1) = 0.d0

         CASE (BV__GGVV_NEUM_HOMOG)

           ALLOCATE (bc_table(b) % bdata(1,1))
           bc_table(b) % bdata(1,1) = 0.d0
         !----------------------------------------------------------------------


         CASE DEFAULT
         WRITE(*,*) ' Boundary value of unknown type. STOP.'
         WRITE(*,*) '   In MODULE ns_viscous_bc'
         WRITE(*,*) '      SUBROUTINE init_ns_viscous_bc'
         STOP
               
      END SELECT
   
   ENDDO
   
   
   ! The number of solid wall points is computed
   ! and the list of wall points is stored
   Np_wall = 0
   DO b = 1, SIZE(bc_table)
   
      IF (bc_table(b) % value_type == BV__VV_NOSLIP) THEN
        Np_wall = Np_wall + SIZE(bc_table(b) % bpoints) 
      ENDIF

   ENDDO

   ALLOCATE (solid_wall_points(Np_wall))
   
   Np_wall = 0
   DO b = 1, SIZE(bc_table)
   
      IF (bc_table(b) % value_type == BV__VV_NOSLIP) THEN
      
         Np_bound = SIZE(bc_table(b) % bpoints)
         solid_wall_points(Np_wall+1: Np_wall+Np_bound) = bc_table(b) % bpoints 
         Np_wall = Np_wall + Np_bound
      
      ENDIF

   ENDDO
   
   
   ! Used only in the subroutine NS_VISCOUS_BOUNDARY_FLUX 
   C_visc = 1 / REF % Re
   
   END SUBROUTINE  init_ns_viscous_bc 




   !============================================================ 
   FUNCTION vb__ns_viscous_bc( vv_b, normal, rr_b )    &
                                       RESULT( vb )
   !============================================================ 
     

      !------------------------------------------------------------ 
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  vv_b
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  normal
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  rr_b

      REAL(KIND=8), DIMENSION(SIZE(vv_b,1), SIZE(vv_b,2))  ::  vb
      !------------------------------------------------------------ 
      INTEGER                                  ::  value_type
      REAL(KIND=8), DIMENSION(:,:), POINTER    ::  bdata
      INTEGER                                  ::  Np_bound
      INTEGER,      DIMENSION(:),   POINTER    ::  bpoints
      REAL(KIND=8), DIMENSION(SIZE(normal,1))  ::  nn  
      INTEGER  ::  i, i_, b

      REAL(KIND=8), DIMENSION(SIZE(rr_b,1),  SIZE(rr_b,2)) :: void_1
      !------------------------------------------------------------ 

      ! Dummy assignation
      void_1 = rr_b

      DO b = 1, SIZE(bc_table)
 
         value_type =  bc_table(b) % value_type
         bdata      => bc_table(b) % bdata
         Np_bound   =  bc_table(b) % Np_bound
         bpoints    => bc_table(b) % bpoints

         SELECT CASE (value_type)            

           !--------------
           CASE(BV__VV_NONE)
           !--------------
           vb(:,bpoints) = vv_b(:,bpoints)
              
           !------------------
           CASE(BV__VV_NOSLIP)
           !------------------
           vb(:,bpoints) = 0.d0
              
           !--------------
           CASE(BV__VV_OO)
           !--------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              vb(:,i) = bdata(:,1)              
           ENDDO

           !-----------------
           CASE(BV__VV_CONST)
           !-----------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              vb(:,i) = bdata(:,1)              
           ENDDO
           
           !--------------------------
           CASE(BV__VV_0)
           !--------------------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              vb(:,i) = bdata(:,i_) 
           ENDDO

           !--------------------------
           CASE(BV__VV_SIMM)
           !--------------------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              nn = normal(:,i) / SQRT( SUM(normal(:,i)**2) )
              vb(:,i) = vv_b(:,i) - SUM( vv_b(:,i)*nn)*nn
           ENDDO              

           !------------------
           CASE(BV__VV_SLIP)
           !------------------
           IF (SIZE(vb,1) == 3) THEN
             PRINT*, ''
             PRINT*, 'ERROR. Slip condition not implemented yet'
             PRINT*, ''
           ENDIF
           
!           DO i_ = 1, Np_bound
!              i = bpoints(i_)
!              nn = normal(:,i) / SQRT( SUM(normal(:,i)**2) )              
!              Vs = Gokcen_Slip_Formula('V',)              
!           ENDDO
           
           !--------------
           CASE DEFAULT
           !--------------
           vb(:,bpoints) = vv_b(:,bpoints)

         END SELECT
            
      ENDDO
      

   END FUNCTION vb__ns_viscous_bc
   !============================================================ 



   !============================================================ 
   FUNCTION GG_vb__ns_viscous_bc( GG_vv_b, normal, rr_b )    &
                                          RESULT( GG_vb )
   !============================================================ 
     

      !------------------------------------------------------------ 
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  ::  GG_vv_b
      REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  ::  normal
      REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  ::  rr_b

      REAL(KIND=8), DIMENSION( SIZE(GG_vv_b,1), &
                               SIZE(GG_vv_b,2), &
                               SIZE(GG_vv_b,3)  )  ::  GG_vb
      !------------------------------------------------------------ 
      INTEGER                                 :: value_type
      !REAL(KIND=8), DIMENSION(:,:), POINTER  :: bdata
      INTEGER                                 :: Np_bound
      INTEGER,      DIMENSION(:),   POINTER   :: bpoints
      REAL(KIND=8), DIMENSION(SIZE(normal,1)) :: nn  
      INTEGER  ::  i, i_, k, b

      REAL(KIND=8), DIMENSION(SIZE(rr_b,1),  SIZE(rr_b,2)) :: void_1
      !------------------------------------------------------------ 

      ! Dummy assignation
      void_1 = rr_b

     
      DO b = 1, SIZE(bc_table)
 
         value_type =  bc_table(b) % value_type
         bpoints    => bc_table(b) % bpoints
         Np_bound   =  bc_table(b) % Np_bound
         !bdata      => bc_table(b) % bdata

         SELECT CASE ( value_type )            
           
           !------------------
           CASE(BV__GGVV_NONE)
           !------------------
           GG_vb(:,:,bpoints) = GG_vv_b(:,:,bpoints)


           !------------------------
           CASE(BV__GGVV_NEUM_HOMOG)
           !------------------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              nn = normal(:,i) / SQRT( SUM(normal(:,i)**2) )
              DO k = 1, SIZE(GG_vb,1)
                GG_vb(k,:,i) = GG_vv_b(k,:,i) - SUM(GG_vv_b(k,:,i)*nn)*nn
              ENDDO
           ENDDO
	   
	   
           !------------------
           CASE DEFAULT
           !------------------
           GG_vb(:,:,bpoints) = GG_vv_b(:,:,bpoints)
                 
         END SELECT

      ENDDO
      

   END FUNCTION GG_vb__ns_viscous_bc
   !============================================================ 



   !============================================================ 
   SUBROUTINE ns_impose_strong_vv_bc(ww, jd_jb, normal, rr_b)
   !============================================================ 
     

      !------------------------------------------------------------ 
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT)  ::  ww
      INTEGER,      DIMENSION(:),   INTENT(IN)     ::  jd_jb
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  normal
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  rr_b
      !------------------------------------------------------------ 
      INTEGER                                  ::  cond_type, value_type
      REAL(KIND=8), DIMENSION(:,:), POINTER    ::  bdata
      INTEGER                                  ::  Np_bound
      INTEGER,      DIMENSION(:),   POINTER    ::  bpoints
      REAL(KIND=8), DIMENSION(SIZE(normal,1))  ::  nn  
      INTEGER  ::  i, i_, b, p

      REAL(KIND=8), DIMENSION(SIZE(rr_b,1),  SIZE(rr_b,2))   :: void_1
      !------------------------------------------------------------ 

      ! Dummy assignation
      void_1 = rr_b
     
      p = SIZE(ww,1)
      
      DO b = 1, SIZE(bc_table)
 
         cond_type  =  bc_table(b) % cond_type
         value_type =  bc_table(b) % value_type
         bdata      => bc_table(b) % bdata
         Np_bound   =  bc_table(b) % Np_bound
         bpoints    => bc_table(b) % bpoints

         IF (cond_type == BC__STRONG) THEN

            SELECT CASE (value_type)            

              !--------------
              CASE(BV__VV_NONE)
              !--------------

              !------------------
              CASE(BV__VV_NOSLIP)
              !------------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 ww(p,i) = ww(p,i) - 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
                 ww(2:p-1,i) = 0.d0
              ENDDO

              !--------------
              CASE(BV__VV_OO)
              !--------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 ww(p,i) = ww(p,i) - 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
                 ww(2:p-1,i) = bdata(:,1) * ww(1,i)              
                 ww(p,i) = ww(p,i) + 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
              ENDDO

              !-----------------
              CASE(BV__VV_CONST)
              !-----------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 ww(p,i) = ww(p,i) - 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
                 ww(2:p-1,i) = bdata(:,1) * ww(1,i)              
                 ww(p,i) = ww(p,i) + 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
              ENDDO

              !--------------------------
              CASE(BV__VV_0)
              !--------------------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 ww(p,i) = ww(p,i) - 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
                 ww(2:p-1,i) = bdata(:,i_) * ww(1,i)              
                 ww(p,i) = ww(p,i) + 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
              ENDDO

              !--------------------------
              CASE(BV__VV_SIMM)
              !--------------------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 ww(p,i) = ww(p,i) - 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
                 nn = normal(:,bpoints(i_)) &
                    / SQRT( SUM(normal(:,bpoints(i_))**2) )
                 ww(2:p-1,i) = ww(2:p-1,i) - SUM( ww(2:p-1,i)*nn)*nn
                 ww(p,i) = ww(p,i) + 0.5*SUM( ww(2:p-1,i)**2 ) / ww(1,i)
              ENDDO              

            END SELECT
            
         ENDIF
            
      ENDDO

   END SUBROUTINE ns_impose_strong_vv_bc


END MODULE ns_viscous_bc

