!============================================================ 
!
!      Module: np_metric_gen
!
! Description: Computes all metrics quantities pertaining to
!              the domain and its boundary 
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

!=====================================================================
!
!   ---------------------------
!   Domain node-pair quantities
!   ---------------------------
!
!            mass(c)   Mass associated to the c-th node-pair
!
!           eta(k,c)   Metric vector (of the physical space)
!                      associated to the c-th node-pair
!
!           stiff(c)   Stiffnes associated to the c-th node-pair
!
!         stiff_T(c)   Stiffness tensor associate to the c-th node-pair
!
!           Drr(1,c)   Length of the node-pair  i--j
!           Drr(2,c)   Length of the extension i*--j
!           Drr(3,c)   Length of the extension  i--j*
!
!   ----------------------
!   Domain node quantities
!   ----------------------
!
!            cell(i)   Cell size == diagonally lumped mass matrix
!
!         mass_ii(i)   Diagonal element of the mass matrix
!
!        stiff_ii(i)   Diagonal element of the stiffness matrix
!
!
!   -----------------------------
!   Boundary node-pair quantities
!   -----------------------------
!
!          mass_b(c)   Mass (surface integral) associated to
!                      the c-th surface node-pair
!
!         chi_b(k,c)   Boundary metric vector (of the lhysical space)
!                      associated to the c-th surface node-pair
!
!
!   ------------------------
!   Boundary node quantities
!   ------------------------
!
!         xi_bp(k,i)   Boundary metric vector (of the physical space)
!                      associated to the i-th surface point
!
!   ------------------------
!   Finite Volume quantities
!   ------------------------
!
!        eta_fv(k,c)   Metric vector (of the physical space)
!                      associated to the c-th node-pair
!
!=====================================================================

   MODULE  np_metric_gen

   USE metric_coefficients
   USE fem_ele_types
   USE fem_ref_elements
   USE fem_gauss_points
   USE nodes
   USE mesh_structure
   USE node_pair_structure

   CONTAINS


   SUBROUTINE  metrics_gen
   !--------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(element_domain)   :: ele_d
   TYPE(element_boundary) :: ele_b
   
   TYPE(reference_element), DIMENSION(:), POINTER :: ref_ele_dom,  &
   						     ref_ele_bou
   INTEGER :: m, c, i, j, i_b
   !--------------------------------------------------------------------------

   ! Allocation of metric quantities
   ALLOCATE (mass(Nc_d), stiff(Nc_d), eta(k_d, Nc_d), stiff_T(k_d, k_d, Nc_d))
   ALLOCATE (mass_ii(Np_d), stiff_ii(Np_d), cell(Np_d))
   ALLOCATE (mass_b(Nc_b), chi_b(k_d, Nc_b), kappa_b(2,k_d, Nc_b), xi_bp(k_d, Np_b))

   ALLOCATE (cosines_ETA_fv(k_d, k_d, Nc_d), cosines_XI_BP(k_d, k_d, Np_b))
   	
   ! Axisymmetric metric quantities
   ALLOCATE (mass_y(Nc_d))
   ALLOCATE (mass_y_ii(Np_d))
   ALLOCATE (eta_y(k_d, Nc_d))     
   ALLOCATE (xi_y_bp(k_d, Np_b))
   ALLOCATE (cell_y(Np_d))
   
   ! ALLOCATE (Drr(4*k_d-1, Nc_d)) --> Moved to np_topology_gen: extend_node_pair
   ! WARNING NP_D == NN_D HERE !!! (ONLY P1 ELEMENTS)



   ! FINITE ELEMENTS DOMAIN METRIC QUANTITIES -------------------------------------
   WRITE(*,*) '   Domain finite element metrics...'

   ! Initialization of the metric quantities
   eta     = 0.d0;    stiff    = 0.d0
   mass    = 0.d0;    stiff_T  = 0.d0
   mass_ii = 0.d0;    stiff_ii = 0.d0
 
   !rot_stiff = 0.d0
   	
   ! Axisymmetric metric quantities
   eta_y     = 0.d0
   mass_y    = 0.d0
   mass_y_ii = 0.d0

   ! Generation of Gauss points
   CALL init_ref_elements(k_d, ref_ele_dom, ref_ele_bou)

   DO m = 1, Ne_d
   
       ele_d % ele_type = ele_type_d(m)
   
       ALLOCATE (ele_d % nn(SIZE(j_m_d(m) % vec)))
       ele_d % nn = j_m_d(m) % vec
       
       CALL gen_gauss_ele_dom(ele_d, rr(:,j_m_d(m) % vec),    &
   			      ref_ele_dom(ele_d % ele_type))
       
       CALL metric_dom(ele_d, j_m_d(m) % vec, c_m_d(m) % vec, &
   		       j_c_d(:,c_m_d(m) % vec))
       
       DEALLOCATE (ele_d%nn,   ele_d%dw, ele_d%rr_G, &
   		   ele_d%dw_C, ele_d%rjp)

   ENDDO 


   ! Diagonal elements of the stiffness matrix
   !
   ! stiff_ii = - SUM_(k in K_i,!=) stiff_ik

   DO c = 1, Nc_d
     
     i = j_c_d(1,c)
     j = j_c_d(2,c)	   
     
     stiff_ii(i) = stiff_ii(i) - stiff(c)
     stiff_ii(j) = stiff_ii(j) - stiff(c)
   
   ENDDO


   ! Cell size === diagonally lumped mass matrix
   !
   ! cell_i = SUM_(k in K_i) mass_ik

   cell = mass_ii
   DO c = 1, Nc_d
   
     i = j_c_d(1,c)
     j = j_c_d(2,c)

     cell(i) = cell(i) + mass(c)
     cell(j) = cell(j) + mass(c)

   ENDDO

   ! Axisymmetric
   cell_y = mass_y_ii
   DO c = 1, Nc_d

      i = j_c_d(1,c) 
      j = j_c_d(2,c)

      cell_y(i) = cell_y(i) + mass_y(c) 
      cell_y(j) = cell_y(j) + mass_y(c) 

   ENDDO


   ! Computation of the length of each  
   ! node-pair and its extensions
   !
   ! Marco Fossati: Now, for convenience, this computation 
   !		    is performed into module 'np_topology_gen', 
   !		    subroutine 'extend_node_pair'.
   !
   !DO c = 1,  Nc_d

   !  !i*---i---j---j*
   !  i  = j_c_d(1,c);  j  = j_c_d(2,c)
   !  is = j_c_d(3,c);  js = j_c_d(4,c)

     ! lenght of i--j
   !  Drr(1,c) = SQRT(SUM((rr(:,j) - rr(:,i))**2))

     ! lenght of i--i*
   ! Drr(2,c) = SQRT(SUM((rr(:,i) - rr(:,is))**2))

     !lenght of j--j*
   !  IF (js .NE. i) THEN
   !	Drr(3,c) = SQRT(SUM((rr(:,js) - rr(:,j))**2))
   !  ELSE
   !	Drr(3,c) =  - Drr(1,c) 
   !  ENDIF

   !ENDDO





   ! FINITE VOLUME DOMAIN METRIC QUANTITIES ---------------------------------------
   WRITE(*,*) '   Domain finite volumes metrics...'

   ALLOCATE (eta_fv(k_d, Nc_fv), eta_y_fv(k_d, Nc_fv))

   !ALLOCATE (Drr_fv(3, Nc_fv)) --> Moved to np_topology_gen: extend_node_pair

   eta_fv   = 0.d0
   eta_y_fv = 0.d0

   CALL compute_eta_fv(j_m_d, ele_type_d, c_m_fv, j_c_fv, rr) !Modifica
   

   ! Computation of the length of each  
   ! node-pair and its two extensions
   !
   ! Marco Fossati: Now, for convenience, this computation 
   !		    is performed into module 'np_topology_gen', 
   !		    subroutine 'extend_node_pair'.
   !
   !DO c = 1,  Nc_fv
   !
   !   !  i*---i---j---j*
   !   i  = j_c_fv(1,c);  j  = j_c_fv(2,c)
   !   is = j_c_fv(3,c);  js = j_c_fv(4,c)
   !
   !   ! lenght of i--j
   !   Drr_fv(1,c) = SQRT( SUM( (rr(:,j)  - rr(:,i))**2 ) )
   !
   !   ! lenght of i--i*
   !   Drr_fv(2,c) = SQRT( SUM( (rr(:,i)  - rr(:,is))**2 ) )
   !
   !   ! lenght of j--j*
   !   IF ( js .NE. i ) THEN
   !	  Drr_fv(3,c) = SQRT( SUM( (rr(:,js) - rr(:,j))**2 ) )
   !   ELSE
   !	  Drr_fv(3,c) =  - Drr_fv(1,c) 
   !   ENDIF
   !
   !ENDDO





   ! BOUNDARY METRIC QUANTITIES ---------------------------------------------------
   WRITE(*,*) '   Boundary metrics...'

   ! Initialization of the metric quantities
   mass_b  = 0.d0
   chi_b   = 0.d0  
   kappa_b = 0.d0   
   xi_bp   = 0.d0
   
   ! Axisymmetric metric quantities
   xi_y_bp = 0.d0

   DO m = 1, Ne_b

       ele_b % ele_type = ele_type_b(m)
       
       ALLOCATE (ele_b % nn(SIZE(j_m_b(m)%vec)), &
   		 ele_b % nb(SIZE(j_m_b(m)%vec)))

       ele_b % nn = jd_jb(j_m_b(m)%vec)
       ele_b % nb = j_m_b(m) % vec
   	
       CALL gen_gauss_ele_bou(ele_b, rr(:,ele_b%nn), ref_ele_bou(ele_b%ele_type))
       
       CALL metric_bou(ele_b, j_m_b(m)%vec, c_m_b(m)%vec)
 
       DEALLOCATE (ele_b % nn,  ele_b % nb, ele_b % rr_G, &
   		   ele_b % rjp, ele_b % rnorms)
       
   ENDDO

   ! One dimensional case
   IF (k_d == 1) THEN

      DO c = 1, Nc_d

   	 i = j_c_d(1,c)
         j = j_c_d(2,c)
 
   	 DO i_b = 1, Np_b
     
   	    IF(jd_jb(i_b) == i)  xi_bp(:,i_b) = - eta(:,c)
   	    IF(jd_jb(i_b) == j)  xi_bp(:,i_b) =   eta(:,c)

   	 ENDDO

      ENDDO
   ENDIF
 

   END SUBROUTINE  metrics_gen





   SUBROUTINE metric_dom(ele_d, j_m, c_m, j_c)
   !------------------------------------------------------------ 
   IMPLICIT NONE

   INTEGER, DIMENSION(:),   INTENT(IN) :: j_m, c_m
   INTEGER, DIMENSION(:,:), INTENT(IN) :: j_c
   TYPE(element_domain),    INTENT(IN) :: ele_d
   !------------------------------------------------------------ 

   ! Indices
   INTEGER :: c       ! Node-pair     (global index)
   INTEGER :: c_      ! Node-pair     (local index)
   INTEGER :: i       ! Node i        (global index)
   INTEGER :: i_, j_  ! Nodes i and j (local element index)

   ! Elemental quantities
   INTEGER :: n_w, l_G, k, h   
   REAL(KIND=8), DIMENSION(:,:),   POINTER :: w
   REAL(KIND=8), DIMENSION(:,:,:), POINTER :: dw
   REAL(KIND=8), DIMENSION(:),     POINTER :: rjp
   REAL(KIND=8), DIMENSION(:,:),   POINTER :: rr_G

   ! Others
   REAL(KIND=8) :: s ! Used to accumulate the stiffness vector

   ! Local orientation of the node-pair with respect 
   ! to the globally defined one.
   INTEGER :: cn_orientation
   !------------------------------------------------------------ 

   !------------------------------------------------------------ 
   ! Retrieve the metric of the element
   ! ----------------------------------
   n_w  =  ele_d % n_w 
   l_G  =  ele_d % l_G

   w	=>  ele_d % ww 
   dw	=>  ele_d % dw
   rjp  =>  ele_d % rjp
   
   rr_G =>  ele_d % rr_G
   !------------------------------------------------------------ 

   
   !------------------------------------------------------------ 
   ! Computes metric quantities for every node-pair 
   ! and node of the element
   ! ----------------------------------------------
   c_ = 0
   DO i_ = 1, n_w;  i = j_m(i_)


      !------------------------------------------------------------ 
      ! Nodal quantities
      ! ----------------

      ! Mass = N_i  N_i
      mass_ii(i)   = mass_ii(i)   + SUM(	    w(i_,:)*w(i_,:) * rjp)
      
      ! Mass_y = y N_i N_i  (Axisymmetric)
      IF (k_d == 2)  &
      mass_y_ii(i) = mass_y_ii(i) + SUM(rr_G(2,:) * w(i_,:)*w(i_,:) * rjp)


      DO j_ = i_+1, n_w
      
   	 c_ = c_ + 1;  c = c_m(c_)

   	 !------------------------------------------------------------ 
   	 ! Node-pair quantities
   	 ! --------------------

   	 ! Mass = N_i N_j
   	 mass(c)    =  mass(c)   +  SUM(	    w(i_,:)*w(j_,:) * rjp)

   	 ! Mass_y = y N_i N_j  (Axisymmetric)
   	 IF (k_d == 2)  &
     	 mass_y(c)  =  mass_y(c) +  SUM(rr_G(2,:) * w(i_,:)*w(j_,:) * rjp)



   	 ! Eta = N_i G(N_j) - N_j G(N_i)
   	 cn_orientation = 1
   	 IF (i .NE. j_c(1,c_)) cn_orientation = -1

   	 DO k = 1, k_d
   	   eta(k,c)	= eta(k,c)   + DBLE(cn_orientation) * SUM(	      (  w(i_,:)*dw(k,j_,:)  &
   									       - w(j_,:)*dw(k,i_,:)) * rjp)
   	 ENDDO

   	 ! Eta_y = y (N_i G(N_j) - N_j G(N_i))  (Axisymmetric)
   	 IF (k_d == 2) THEN
   	   DO k = 1, k_d
   	     eta_y(k,c) = eta_y(k,c) + DBLE(cn_orientation) * SUM(rr_G(2,:) * (  w(i_,:)*dw(k,j_,:)  &
   									       - w(j_,:)*dw(k,i_,:)) * rjp)
   	   ENDDO
   	 ENDIF     


   		     
   	 ! stiff = G(N_i).G(N_j)
   	 ! stiff tensor_kh  = G_k(N_i) G_h(N_j) 
   	 ! ------------------------------------
   	 s = 0.d0
   	 DO k = 1, k_d

   	    s = s + SUM( dw(k,i_,:) * dw(k,j_,:) * rjp )
   			    
   	    DO h = 1, k_d

   	       stiff_T(k,h,c)  =  stiff_T(k,h,c)  &
   			       +  SUM( dw(k,i_,:) * dw(h,j_,:) * rjp )

   	    ENDDO

   	 ENDDO

   	 stiff(c) = stiff(c) + s
   	 !------------------------------------------------------------ 

   	 ! Rotated stiffness = G(N_i) x z . G(N_j)
   	 !------------------------------------------------------------ 
   	 !IF (k_d == 2) THEN
   	 !
   	 !   rot_stiff(c) = rot_stiff(c)  &
   	 !		  + SUM( (dw(2,i_,:) * dw(1,j_,:)   & 
   	 !			- dw(1,i_,:) * dw(2,j_,:)) * rjp )
   	 !
   	 !ENDIF        
   	 !------------------------------------------------------------ 


      ENDDO

   ENDDO


   END SUBROUTINE metric_dom





   SUBROUTINE metric_bou(ele_b, j_m, c_m)
   !------------------------------------------------------------ 
   IMPLICIT NONE

   TYPE(element_boundary), INTENT(IN) :: ele_b
   INTEGER, DIMENSION(:),  INTENT(IN) :: j_m, c_m

   ! Indices
   INTEGER :: c	! Node-pair	(global index)
   INTEGER :: c_	! Node-pair	(local index)
   INTEGER :: i	! Node i	(global index)
   INTEGER :: i_, j_  ! Nodes i and j (local element index)
   
   ! Elemental quantities
   INTEGER  ::  n_w, l_G, k 
   REAL(KIND=8), DIMENSION(:,:), POINTER :: w
   REAL(KIND=8), DIMENSION(:),   POINTER :: rjp
   REAL(KIND=8), DIMENSION(:,:), POINTER :: rnorms  
   REAL(KIND=8), DIMENSION(:,:), POINTER :: rr_G
   !------------------------------------------------------------ 


   !------------------------------------------------------------ 
   ! Retrieve the metric of the element
   ! ----------------------------------
   n_w  =  ele_b % n_w 
   l_G  =  ele_b % l_G

   w	   =>  ele_b % ww 
   rjp     =>  ele_b % rjp
   rnorms  =>  ele_b % rnorms
   
   rr_G    =>  ele_b % rr_G
   !------------------------------------------------------------ 


   !------------------------------------------------------------ 
   ! Computes metric quantities for every node-pair 
   ! and node of the element
   ! ----------------------------------------------
   c_ = 0
   DO i_ = 1, n_w;  i = j_m(i_)

      !------------------------------------------------------------ 
      ! ----------------
      ! Nodal quantities
      ! ----------------

      ! xi_bp(1:k_d, i) = Ns_i  n dS   
      ! (n is the unit normal outward vector)
      ! -------------------------------------
      DO k = 1, k_d
   	 xi_bp(k,i) = xi_bp(k,j_m(i_)) & 
   		    + SUM( rnorms(k,:) * w(i_,:) * rjp )
      ENDDO

      ! Axisymmetric metric quantity
      DO k = 1, k_d
   	 xi_y_bp(k,i) = xi_y_bp(k,j_m(i_)) & 
   		    + SUM( rr_G(2,:) * rnorms(k,:) * w(i_,:) * rjp )
      ENDDO
      !------------------------------------------------------------ 


      DO j_ = i_+1, n_w
   	 c_ = c_ + 1;  c  = c_m(c_)

   	 !------------------------------------------------------------ 
   	 ! Node-pair quantities
   	 ! --------------------

   	 ! chi_b(1:k_d, cs) = Ns_i  Ns_j  n dS   
   	 ! (n is the unit normal outward vector)
   	 ! -------------------------------------
   	 DO k = 1, k_d
   	    chi_b(k,c) = chi_b(k,c) & 
   	       + SUM( rnorms(k,:) * w(i_,:) * w(j_,:) * rjp ) 
   	 ENDDO
   	 
   	 ! kappa_b(1, 1:k_d, cs) = Ns_i  Ns_i  n dS   
   	 ! kappa_b(2, 1:k_d, cs) = Ns_j  Ns_j  n dS   
   	 ! (n is the unit normal outward vector)
   	 ! -------------------------------------
   	 DO k = 1, k_d
   	    kappa_b(1,k,c) = kappa_b(1,k,c) & 
   	       + SUM( rnorms(k,:) * w(i_,:)**2 * rjp ) 
   	    kappa_b(2,k,c) = kappa_b(2,k,c) & 
   	       + SUM( rnorms(k,:) * w(j_,:)**2 * rjp ) 
   	 ENDDO

   	 ! mass_b(c) = Ns_i  Ns_j  dS  !! SCALAR quantity
   	 ! ----------------------------------------------
   	 mass_b(c) = mass_b(c) + SUM( w(i_,:) * w(j_,:) * rjp )
   	 !------------------------------------------------------------ 
   
      ENDDO

   ENDDO
   !------------------------------------------------------------ 

   END SUBROUTINE metric_bou





   SUBROUTINE metric_check
   !------------------------------------------------------------ 
   IMPLICIT NONE

   INTEGER  :: c, i, i_max
   REAL(KIND=8)  :: resid_i, maxresid, domain_area
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  ::  resid
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  ::  resid_norm
   LOGICAL,	 DIMENSION(:),   ALLOCATABLE  ::  boundary
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  ::  boundary_closure
   !------------------------------------------------------------ 

   ALLOCATE (resid(k_d, Np_d), resid_norm(Np_d), boundary(Np_d))
   
   maxresid = -HUGE(maxresid)
   resid = 0.d0    
   i_max = 1

   ! Set the flag for boundary nodes
   boundary = .FALSE.
   boundary(jd_jb) = .TRUE.


   ! FINITE ELEMENT ETA RESIDUAL ----------------------------------------------------
   
   ! Residual on INTERNAL POINTS
   DO c = 1, SIZE(j_c_d,2)
     resid(:,j_c_d(1,c)) = resid(:,j_c_d(1,c)) + eta(:,c)
     resid(:,j_c_d(2,c)) = resid(:,j_c_d(2,c)) - eta(:,c)
   END DO

   i_max = 1
   DO i = 1, Np_d

      IF (.NOT. boundary(i)) THEN

   	 resid_i = DOT_PRODUCT(resid(:,i),resid(:,i))
   	 resid_i = resid_i / cell(i)

   	 IF (resid_i > maxresid) THEN
   	   maxresid = resid_i
   	   i_max = i
   	 ENDIF

      ENDIF

   ENDDO

   WRITE (*,1000) maxresid
   WRITE (*,*) '     - on domain node ', i_max
   WRITE (*,*) '     - of coordinates ', rr(:,i_max)



   ! Residual on BOUNDARY POINTS
   DO i = 1, Np_b   
     resid(:,jd_jb(i)) = resid(:,jd_jb(i)) + xi_bp(:,i)
   ENDDO

   i_max = 1
   maxresid = -HUGE(maxresid)
   DO i = 1, Np_b   

     resid_i = DOT_PRODUCT(resid(:,jd_jb(i)),resid(:,jd_jb(i)))
     resid_i = resid_i / cell(jd_jb(i))

     IF (resid_i > maxresid) THEN
     	maxresid = resid_i
     	i_max = i
     ENDIF

   END DO

   WRITE (*,1001) maxresid
   WRITE (*,*) '     - on boundary node ', i_max, '   (domain  ', jd_jb(i_max), ')'
   WRITE (*,*) '     - of coordinates ', rr(:,jd_jb(i_max))
   WRITE (*,*) 
   !--------------------------------------------------------------------------------- 



   ! DOMAIN AREA(VOLUME) CHECK:
   ! Domain area(volume) check on cell
   domain_area = 0.d0
   DO i = 1, Np_d
     domain_area = domain_area + cell(i) 
   ENDDO
   WRITE(*,*) '     Domain area (cell) =', domain_area

   ! Domain area(volume) check on mass
   domain_area = 0.d0
   DO c = 1, Nc_d
     domain_area = domain_area + 2.d0*mass(c) 
   ENDDO
   DO i = 1, Np_d
     domain_area = domain_area + mass_ii(i) 
   ENDDO
   WRITE(*,*) '     Domain area (mass,mass_ii) =', domain_area




   ! BOUNDARY CHECK:
   ! Boundary check on xi_bp
   ALLOCATE (boundary_closure(k_d))
   boundary_closure = 0.d0
   
   DO i = 1, Np_b
     boundary_closure = boundary_closure + xi_bp(:,i)
   ENDDO

   maxresid = SQRT(DOT_PRODUCT(boundary_closure,boundary_closure))
   WRITE(*,*) '     Boundary closure (xi) =', maxresid

   ! Boundary check on chi_b
   boundary_closure = 0.d0
   
   DO i = 1, Nc_b
     boundary_closure = boundary_closure + chi_b(:,i)
   ENDDO
   
   maxresid = SQRT(DOT_PRODUCT(boundary_closure,boundary_closure))
   WRITE(*,*) '     Boundary closure (chi) =', maxresid
   WRITE(*,*) 



   ! FINITE VOLUME ETA RESIDUAL -----------------------------------------------------
   resid = 0.d0
   i_max = 1
   
   ! Residual on INTERNAL POINTS
   DO c = 1, SIZE(j_c_fv,2)  
     resid(:,j_c_fv(1,c)) = resid(:,j_c_fv(1,c)) + eta_fv(:,c)
     resid(:,j_c_fv(2,c)) = resid(:,j_c_fv(2,c)) - eta_fv(:,c)
   ENDDO

   maxresid = -HUGE(maxresid)
   DO i = 1, Np_d

      IF (.NOT. boundary(i)) THEN

   	 resid_i = DOT_PRODUCT(resid(:,i),resid(:,i))
   	 resid_i = resid_i / cell(i)

   	 IF (resid_i > maxresid) THEN
   	   maxresid = resid_i
   	   i_max = i
   	 ENDIF

      ENDIF

   ENDDO

   WRITE (*,1002) maxresid
   WRITE (*,*) '     - on domain node ', i_max
   WRITE (*,*) '     - of coordinates ', rr(:,i_max)



   ! Residual on BOUNDARY POINTS
   DO i = 1, Np_b   
     resid(:,jd_jb(i)) = resid(:,jd_jb(i)) + xi_bp(:,i)
   END DO

   i_max = 1
   maxresid = -HUGE(maxresid)
   DO i = 1, Np_b   

     resid_i = DOT_PRODUCT(resid(:,jd_jb(i)),resid(:,jd_jb(i)))
     resid_i = resid_i / cell(jd_jb(i))

     IF (resid_i > maxresid) THEN
       maxresid = resid_i
       i_max = i
     ENDIF

   ENDDO

   WRITE (*,1003) maxresid
   WRITE (*,*) '     - on boundary node ', i_max, '   (domain  ', jd_jb(i_max), ')'
   WRITE (*,*) '     - of coordinates ', rr(:,jd_jb(i_max))
   WRITE (*,*) 
   !------------------------------------------------------------ 

1000 FORMAT ('      FEM maximum eta residual on DOMAIN points --->', es24.16)
1001 FORMAT ('      FEM maximum eta residual on BOUNDARY points --->', es24.16)
1002 FORMAT ('      FVM maximum eta residual on DOMAIN points --->', es24.16)
1003 FORMAT ('      FVM maximum eta residual on BOUNDARY points --->', es24.16)

   END SUBROUTINE metric_check





   SUBROUTINE metric_check_axi
   !------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER  :: c, i, i_max
   REAL(KIND=8)  :: resid_i, maxresid, domain_area
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  ::  resid
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  ::  resid_norm
   LOGICAL,	 DIMENSION(:),   ALLOCATABLE  ::  boundary
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE  ::  boundary_closure
   !------------------------------------------------------------------

   ALLOCATE (resid(k_d, Np_d), resid_norm(Np_d), boundary(Np_d))
   
   maxresid = -HUGE(maxresid)
   resid = 0.d0
   i_max = 1

   ! Set the flag for boundary nodes
   boundary = .FALSE.
   boundary(jd_jb) = .TRUE.


   ! FINITE ELEMENT ETA RESIDUAL ----------------------------------------------------

   ! Residual on INTERNAL POINTS
   DO c = 1, SIZE(j_c_d,2)  

      ! Node i of c
      resid(:,j_c_d(1,c)) = resid(:,j_c_d(1,c)) + eta_y(:,c)
      resid(2,j_c_d(1,c)) = resid(2,j_c_d(1,c)) - mass(c)

      ! Node j of c
      resid(:,j_c_d(2,c)) = resid(:,j_c_d(2,c)) - eta_y(:,c)
      resid(2,j_c_d(2,c)) = resid(2,j_c_d(2,c)) - mass(c)

   ENDDO
   
   resid(2,:) = resid(2,:) - mass_ii
   ! Or equivalently:
   !resid(2,:) = resid(2,:) - cell

   i_max = 1
   DO i = 1, Np_d

      IF (.NOT. boundary(i)) THEN

   	 resid_i = DOT_PRODUCT(resid(:,i),resid(:,i))
   	 resid_i = resid_i / cell(i)

   	 IF (resid_i > maxresid) THEN
   	   maxresid = resid_i
   	   i_max = i
   	 ENDIF

      ENDIF

   ENDDO

   WRITE (*,1004) maxresid
   WRITE (*,*) '     - on domain node ', i_max
   WRITE (*,*) '     - of coordinates ', rr(:,i_max)



   ! Residual on BOUNDARY POINTS
   DO i = 1, Np_b   
     resid(:,jd_jb(i)) = resid(:,jd_jb(i)) + xi_y_bp(:,i)
   ENDDO

   i_max = 1
   maxresid = -HUGE(maxresid)
   DO i = 1, Np_b   

     resid_i = DOT_PRODUCT(resid(:,jd_jb(i)),resid(:,jd_jb(i)))
     resid_i = resid_i / cell(jd_jb(i))

     IF (resid_i > maxresid) THEN
     	maxresid = resid_i
     	i_max = i
     ENDIF

   ENDDO

   WRITE (*,1005) maxresid
   WRITE (*,*) '     - on boundary node ', i_max, '   (domain  ', jd_jb(i_max), ')'
   WRITE (*,*) '     - of coordinates ', rr(:,jd_jb(i_max))
   WRITE (*,*) 
   !---------------------------------------------------------------------------------



   ! DOMAIN AREA(VOLUME) CHECK:
   ! Domain area(volume) check on cell
   domain_area = 0.d0
   DO i = 1, Np_d
     domain_area = domain_area + cell(i) 
   ENDDO
   WRITE(*,*) '     Domain area (cell) = ', domain_area

   ! Domain area(volume) check on mass
   domain_area = 0.d0
   DO c = 1, Nc_d
     domain_area = domain_area + 2.d0*mass(c) 
   ENDDO
   DO i = 1, Np_d
     domain_area = domain_area + mass_ii(i) 
   ENDDO
   WRITE(*,*) '     Domain area (mass,mass_ii) = ', domain_area




   ! BOUNDARY CHECK:
   ! Boundary check on xi_y_bp
   ALLOCATE (boundary_closure(k_d))
   boundary_closure = 0.d0
   
   DO i = 1, Np_b
     boundary_closure = boundary_closure + xi_y_bp(:,i)
   ENDDO
   boundary_closure(2) = boundary_closure(2) - domain_area

   maxresid = SQRT(DOT_PRODUCT(boundary_closure,boundary_closure))
   WRITE(*,*) '     Boundary closure (xi) = ', maxresid

   ! Boundary check on chi_y_b
   !boundary_closure = 0.d0
   
   !DO i = 1, Nc_b
   !   boundary_closure = boundary_closure + chi_y_b(:,i)
   !ENDDO
   
   !maxresid = SQRT(DOT_PRODUCT(boundary_closure,boundary_closure))
   !WRITE(*,*) '  Boundary closure (chi) =  ', maxresid
   WRITE(*,*)


   ! FINITE VOLUME ETA RESIDUAL -----------------------------------------------------
   resid = 0.d0
   i_max = 1
   
   ! Residual on INTERNAL POINTS
   DO c = 1, SIZE(j_c_fv,2)
     resid(:,j_c_fv(1,c)) = resid(:,j_c_fv(1,c)) + eta_y_fv(:,c)
     resid(:,j_c_fv(2,c)) = resid(:,j_c_fv(2,c)) - eta_y_fv(:,c)
   ENDDO
   resid(2,:) = resid(2,:) - cell

   maxresid = -HUGE(maxresid)
   DO i = 1, Np_d

      IF (.NOT. boundary(i)) THEN

   	 resid_i = DOT_PRODUCT(resid(:,i),resid(:,i))
   	 resid_i = resid_i / cell(i)

   	 IF (resid_i > maxresid) THEN
   	   maxresid = resid_i
   	   i_max = i
   	 ENDIF

      ENDIF

   ENDDO

   WRITE (*,1006) maxresid
   WRITE (*,*) '     - on domain node ', i_max
   WRITE (*,*) '     - of coordinates ', rr(:,i_max)



   ! Residual on BOUNDARY POINTS
   DO i = 1, Np_b   
     resid(:,jd_jb(i)) = resid(:,jd_jb(i)) + xi_y_bp(:,i)
   END DO

   i_max = 1
   maxresid = -HUGE(maxresid)
   DO i = 1, Np_b   

     resid_i = DOT_PRODUCT(resid(:,jd_jb(i)),resid(:,jd_jb(i)))
     resid_i = resid_i / cell(jd_jb(i))

     IF (resid_i > maxresid) THEN
       maxresid = resid_i
       i_max = i
     ENDIF

   ENDDO

   WRITE (*,1000) maxresid
   WRITE (*,*) '     - on boundary node ', i_max, '   (domain  ', jd_jb(i_max), ')'
   WRITE (*,*) '     - of coordinates ', rr(:,jd_jb(i_max))
   WRITE (*,*) '     - resid --> ', resid(:,jd_jb(i_max))
   WRITE (*,*) 
   !---------------------------------------------------------------------------------

   DEALLOCATE (resid)
   ALLOCATE (resid(k_d, Nc_d))

1004 FORMAT ('      FEM maximum eta residual on DOMAIN points --->', es24.16)
1005 FORMAT ('      FEM maximum eta residual on BOUNDARY points --->', es24.16)
1006 FORMAT ('      FVM maximum eta residual on DOMAIN points --->', es24.16)
1000 FORMAT ('      FVM maximum eta residual on BOUNDARY points --->', es24.16)

   END SUBROUTINE metric_check_axi





   SUBROUTINE compute_eta_fv(j_m, ele_type, c_m_fv, j_c_fv, rr)

      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE(D_I_V),  DIMENSION(:),   INTENT(IN)  ::  j_m, c_m_fv
      INTEGER,      DIMENSION(:),   INTENT(IN)  ::  ele_type
      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c_fv
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  rr
      !------------------------------------------------------------ 
      ! Indices
      ! -------
      INTEGER ::  c       ! Node-pair     (global index)
      INTEGER ::  c_      ! Node-pair     (local index)
      INTEGER ::  i,j , i_f, j_f, i_f_, j_f_
      TYPE(D_I_V), DIMENSION(:), POINTER ::  faces

      INTEGER  ::  m, k, f

      REAL(KIND=8)  ::  cn_orientation   ! Local orientation of the 
                         ! node-pair with respect to the global one
      REAL(KIND=8), DIMENSION(SIZE(rr,1)) ::  xg, xf, xm , norm, norm_y
      REAL(KIND=8), DIMENSION(SIZE(rr,1)-1, SIZE(rr,1)) ::  Dx
      INTEGER, DIMENSION(:), POINTER  ::  jm, cm
      INTEGER, DIMENSION(:), ALLOCATABLE  ::  jf
      !------------------------------------------------------------ 

      k_d = SIZE(eta_fv,1)

      SELECT CASE (k_d)

         ! ONE SPACE DIMENSION -------------------------------------------------
	 CASE (1)
         eta_fv = 1.d0



         ! TWO SPACE DIMENSIONS ------------------------------------------------
         CASE (2)

         DO m = 1, SIZE(j_m)

            jm =>    j_m(m) % vec
            cm => c_m_fv(m) % vec

            ! Baricenter of the element m
            DO k = 1, k_d
              xg(k) = SUM(rr(k,jm)) / SIZE(jm)
            ENDDO

            DO c_ = 1, SIZE(cm)

	       c = cm(c_)
	       
               i = j_c_fv(1,c)
	       j = j_c_fv(2,c)
  
               ! Baricenter of the node-pair (midpoint)
               xm = 0.5 * (rr(:,i) + rr(:,j)) 

               Dx(1,:) = xm - xg
               norm = vec_prod(Dx)
               cn_orientation = SIGN(1.d0, SUM(norm * (xm - rr(:,i))))

               ! CARTESIAN - eta_fv
               eta_fv(:,c) = eta_fv(:,c) + cn_orientation * norm

               ! AXISYMMETRIC - eta_y_fv
               norm_y = norm * ABS((xg(2) + xm(2))) / 2.d0
               eta_y_fv(:,c) = eta_y_fv(:,c) + cn_orientation * norm_y
               
            ENDDO

         ENDDO



         ! THREE SPACE DIMENSIONS ----------------------------------------------
         CASE (3)

         DO m = 1, SIZE(j_m)

            jm =>    j_m(m)%vec
            cm => c_m_fv(m)%vec

            ! Baricenter of the element
            ! -------------------------
            DO k = 1, k_d
               xg(k) = SUM( rr(k,jm) ) / SIZE(jm)
            ENDDO

            faces => ele_faces(ele_type(m))

            DO f = 1, SIZE(faces)

               IF ( ALLOCATED(jf) ) DEALLOCATE (jf)
               ALLOCATE( jf(SIZE(faces(f)%vec)) ) 
               jf = jm(faces(f)%vec)

               ! Baricenter of the face
               ! ----------------------
               DO k = 1, k_d
         	  xf(k) = SUM( rr(k,jf) ) / SIZE(jf)
               ENDDO


               Dx(1,:) = xf - xg
            
               DO c_ = 1, SIZE(cm);  c = cm(c_)

         	  i = j_c_fv(1,c);  j = j_c_fv(2,c)
  
         	  DO i_f_ = 1, SIZE(jf); i_f = jf(i_f_)
         	     DO j_f_ = 1, SIZE(jf); j_f = jf(j_f_)

         		!IF (( (i == i_f) .AND. (j == j_f) ) .OR. &
         		!    ( (j == i_f) .AND. (i == j_f) )	  ) THEN
         		IF ( (i == i_f) .AND. (j == j_f) ) THEN
         		    
         		   ! Baricenter of the node-pair (midpoint)
         		   ! --------------------------------------
         		   xm = 0.5 * ( rr(:,i) + rr(:,j) ) 

         		   Dx(2,:) = xm - xf

         		   norm = 0.5*vec_prod(Dx)

         		   cn_orientation = SIGN(1.d0, &
         			     SUM( norm*(rr(:,j) - rr(:,i) ) ) )
 
         		   eta_fv(:,c) = eta_fv(:,c)  &
         			       +  cn_orientation * norm
         		ENDIF
         	     ENDDO
         	  ENDDO

               ENDDO

            ENDDO
         
	    DEALLOCATE(faces)
         
	 ENDDO

      END SELECT

   END SUBROUTINE compute_eta_fv
  
   
   
   
   
   SUBROUTINE  rotation_matrices_gen

   ! cosMat is the cosine direction matrix for the directions aligned 
   ! with normal vector (n) and two orthogonal directions tangent to n.
   ! 1st row contains cosine directors in the direction 'n', the other 
   ! two rows the components of the versors in the two tangential directions
   !---------------------------------------------------------------------------
   USE lin_algebra,   ONLY: vec_prod

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(2, SIZE(cosines_ETA_fv,1)) :: VV
   REAL(KIND=8), DIMENSION(SIZE(cosines_ETA_fv,1))    :: eta_ij, xi_j, t, b, temp

   REAL(KIND=8) :: mod_eta_ij, mod_xi_j, mod_t, mod_b

   INTEGER :: c, j, sD, z_pos
   !---------------------------------------------------------------------------

   sD = SIZE(cosines_ETA_fv,1)

   ! Domain matrices, referred to integrated normal eta_fv
   DO c = 1, Nc_fv
     
     eta_ij = eta_fv(:,c)
     mod_eta_ij = SQRT(SUM(eta_ij**2))
     IF (mod_eta_ij == 0.0) THEN
       eta_ij = 0.0
     ELSE
       eta_ij = eta_ij / mod_eta_ij
     ENDIF  
     
     
     ! Direction along the normal to the interface
     cosines_ETA_fv(1,:,c) = eta_ij


     IF (sD == 2) THEN
     
       ! Direction along the tangent to the interface
       cosines_ETA_fv(2,:,c) = (/ -cosines_ETA_fv(1,2,c), &
                                   cosines_ETA_fv(1,1,c) /)

     ELSEIF (sD == 3) THEN
       
       VV(1,:) = cosines_ETA_fv(1,:,c)					       

       IF (COUNT(cosines_ETA_fv(1,:,c) == 0.d0) == 2) THEN		       

         VV(2,:) = cosines_ETA_fv(1,:,c)			       

         temp = cosines_ETA_fv(1,:,c)*cosines_ETA_fv(1,:,c)
 
         z_pos = MINLOC(temp,1) 				       

         VV(2,z_pos) = VV(2,z_pos) + 1.d0			       

       ELSE							       

         VV(2,:) = (/ cosines_ETA_fv(1,1,c),  &
	              cosines_ETA_fv(1,2,c),  &
		      cosines_ETA_fv(1,3,c) + 1.d0 /) 

       ENDIF							       

       t = vec_prod(VV)

       VV(1,:) = t
       VV(2,:) = cosines_ETA_fv(1,:,c)

       b = vec_prod(VV)

       mod_t = SQRT(SUM(t*t))
       mod_b = SQRT(SUM(b*b))


       ! Directions along the tangent (row 2) and
       ! binormal (row 3) to the interface
       cosines_ETA_fv(2,:,c) = t / mod_t
       cosines_ETA_fv(3,:,c) = b / mod_b

     ENDIF

   ENDDO
   
   
   
   
   ! Boundary matrices, referred to the integrated normals xi_bp
   DO j = 1, Np_b
   
     xi_j = xi_bp(:,j)
     mod_xi_j = SQRT(SUM(xi_j**2))
     IF (mod_xi_j == 0.0) THEN
       xi_j = 0.0
     ELSE
       xi_j = xi_j / mod_xi_j
     ENDIF
     
     ! Direction along the normal to the interface
     cosines_XI_BP(1,:,j) = xi_j


     IF (sD == 2) THEN
     
       ! Direction along the tangent to the interface
       cosines_XI_BP(2,:,j) = (/ -cosines_XI_BP(1,2,j), &
                                  cosines_XI_BP(1,1,j) /)

     ELSEIF (sD == 3) THEN
       
       VV(1,:) = cosines_XI_BP(1,:,j)					       

       IF (COUNT(cosines_XI_BP(1,:,j) == 0.d0) == 2) THEN		       

         VV(2,:) = cosines_XI_BP(1,:,j)			       

         temp = cosines_XI_BP(1,:,j)*cosines_XI_BP(1,:,j)
 
         z_pos = MINLOC(temp,1) 				       

         VV(2,z_pos) = VV(2,z_pos) + 1.d0			       

       ELSE							       

         VV(2,:) = (/ cosines_XI_BP(1,1,j),  &
	              cosines_XI_BP(1,2,j),  &
		      cosines_XI_BP(1,3,j) + 1.d0 /) 

       ENDIF							       

       t = vec_prod(VV)

       VV(1,:) = t
       VV(2,:) = cosines_XI_BP(1,:,j)

       b = vec_prod(VV)

       mod_t = SQRT(SUM(t*t))
       mod_b = SQRT(SUM(b*b))


       ! Directions along the tangent (row 2) and
       ! binormal (row 3) to the interface
       cosines_XI_BP(2,:,j) = t / mod_t
       cosines_XI_BP(3,:,j) = b / mod_b

     ENDIF

   
   ENDDO

   END SUBROUTINE  rotation_matrices_gen


END MODULE  np_metric_gen

