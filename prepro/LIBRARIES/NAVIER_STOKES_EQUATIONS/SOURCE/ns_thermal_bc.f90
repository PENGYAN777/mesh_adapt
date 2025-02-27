!============================================================ 
!
!      Module: ns_viscous_bc
!
! Description: Procedure for the evaluation of the boundary
!              conditions and boundary numerical fluxes for
!              the temperature in the Navier-Stokes equations 
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
!============================================================ 

   MODULE ns_thermal_bc

   !============================================================ 
   USE dynamic_vector
   USE euler_equations

   IMPLICIT NONE
   !============================================================ 

   !============================================================ 
   TYPE boundary_data
   
      ! Type of boundary condition (weak/strong) 
      INTEGER                                ::  cond_type
      ! Type of boundary value 
      INTEGER                                ::  value_type
      ! Size of and vector of boundary data 
      INTEGER                                ::  size_bdata
      REAL(KIND=8), DIMENSION(:,:), POINTER  ::  bdata
      ! Number of and list of boundary nodes 
      ! (boundary indices)
      INTEGER                                ::  Np_bound    
      INTEGER,      DIMENSION(:),   POINTER  ::  bpoints
  
   END TYPE boundary_data
   !============================================================ 


   PRIVATE
   
   ! Boundary condition table
   TYPE(boundary_data), DIMENSION(:), ALLOCATABLE :: bc_table
   
   ! Boundary value types                        
   INTEGER, PARAMETER  ::  BV__T_NONE       = 0,  &
                           BV__T_OO         = 1,  &
                           BV__T_CONST      = 2,  &
                           BV__T_0          = 3,  &
                           BV__GT_NONE      = 10, &
                           BV__GT_ADIABATIC = 11
			   
   CHARACTER(*), DIMENSION(0:11), PARAMETER  :: bv_names = &
                        (/ 'TEMPERATURE: NO BOUNDARY VALUE IMPOSED      ', &
                           'TEMPERATURE: T_OO (TEMPERATURE AT INFINITY) ', &
                           'TEMPERATURE: CONSTANT TEMPERATURE           ', &
                           'TEMPERATURE: INITIAL VALUE                  ', &
                           'UNKNOWN TYPE                                ', &
                           'UNKNOWN TYPE                                ', &
                           'UNKNOWN TYPE                                ', &
                           'UNKNOWN TYPE                                ', &
                           'UNKNOWN TYPE                                ', &
                           'UNKNOWN TYPE                                ', &
                           'GRAD T: NO BOUNDARY VALUE IMPOSED           ', &
                           'GRAD T: ADIABATIC WALL                      ' /) 

   ! Boundary conditions types                                             
   INTEGER, PARAMETER  ::   BC__WEAK   = 0, &
                            BC__STRONG = 1

   CHARACTER(*), DIMENSION(0:1), PARAMETER  ::  bc_names = &
                        (/ 'WEAK FORM  ',  &
                           'STRONG FORM' /)

   REAL(KIND=8) :: C_ther

   PUBLIC   ::   ns_thermal_boundary_flux, &
                 Tb__ns_thermal_bc,        &
               G_Tb__ns_thermal_bc,        &
                ns_impose_strong_T_bc,     &
                      init_ns_thermal_bc,  &
                read_param_ns_thermal_bc,  &
               write_param_ns_thermal_bc    
   !============================================================ 

   CONTAINS

   
   SUBROUTINE  ns_thermal_boundary_flux(G_Tb, kk, xi_bp, chi_b, j_c_b, jd_jb, rhs)
   !-----------------------------------------------------------------------------------
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: G_Tb
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: kk
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: xi_bp
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: chi_b
   INTEGER,      DIMENSION(:),   INTENT(IN) :: jd_jb
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: j_c_b
   
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: rhs
   
   REAL(KIND=8) :: ss
   
   INTEGER :: i, j, c
   !-----------------------------------------------------------------------------------
   
   ! Node-pair contribution
   DO c = 1, SIZE(chi_b,2)
  
     i = j_c_b(1,c)  
     j = j_c_b(2,c)

     ss = SUM(chi_b(:,c) * (kk(jd_jb(j))*G_Tb(:,j) - kk(jd_jb(i))*G_Tb(:,i))) * C_ther
    
     rhs(SIZE(rhs,1),jd_jb(i)) = rhs(SIZE(rhs,1),jd_jb(i)) + ss
     rhs(SIZE(rhs,1),jd_jb(j)) = rhs(SIZE(rhs,1),jd_jb(j)) - ss

   ENDDO   
   
   ! Nodal contribution
   DO i = 1, SIZE(xi_bp,2)
     
     ss = SUM(xi_bp(:,i) * kk(jd_jb(i)) * G_Tb(:,i)) * C_ther
     
     rhs(SIZE(rhs,1),jd_jb(i)) = rhs(SIZE(rhs,1),jd_jb(i)) + ss

   ENDDO       

   END SUBROUTINE  ns_thermal_boundary_flux





   SUBROUTINE  read_param_ns_thermal_bc(idf) 
   !------------------------------------------------------------ 
   USE mp_interface
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf
   INTEGER :: b, N_bound, dType, dVType, dSize
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

   END SUBROUTINE  read_param_ns_thermal_bc


 


   SUBROUTINE  write_param_ns_thermal_bc(idf) 
   !------------------------------------------------------------ 
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf
   INTEGER :: b
   !------------------------------------------------------------ 

   ! Write the boundary conditions table
   WRITE(idf,*) '   Parameters for ns thermal boundary conditions'
   WRITE(idf,*) ' -----------------------------------------------'
   WRITE(idf,*) '   Number of boundaries:', SIZE(bc_table)
   WRITE(idf,*)

   DO b = 1, SIZE(bc_table)

      WRITE(idf,*) '  Boundary', b
      WRITE(idf,*) '	 Boundary value type:	   ', &
   		   bv_names(bc_table(b)%value_type)
      WRITE(idf,*) '	 Boundary condition is imposed in a: ', &
   		   bc_names(bc_table(b)%cond_type)
      WRITE(idf,*)

   ENDDO

   END SUBROUTINE  write_param_ns_thermal_bc 





   SUBROUTINE  init_ns_thermal_bc(ww_oo, ww_b, bound_p) 
   !------------------------------------------------------------ 
   USE dynamic_vector
   USE structures,   ONLY: REF
   
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  ww_oo
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  ww_b
   INTEGER,	 DIMENSION(:),   INTENT(IN)  ::  bound_p

   TYPE(D_I_V), DIMENSION(:),	 ALLOCATABLE  ::  bound_p_DIV, &
   						  p_bound_DIV

   INTEGER ::  i_, i, b, Np_bound, value_type
   !------------------------------------------------------------ 

   ! Boundary to node connectivity
   ALLOCATE (bound_p_DIV(SIZE(bound_p)))
   
   DO i = 1, SIZE(bound_p)
     ALLOCATE (bound_p_DIV(i)%vec(1))
     bound_p_DIV(i)%vec(1) = bound_p(i)
   ENDDO

   ALLOCATE (p_bound_DIV(size_DIV(bound_p_DIV,3)))
   p_bound_DIV = invert_DIV(bound_p_DIV)
   
   DO b = 1, SIZE(bc_table)

     bc_table(b) % Np_bound = SIZE(p_bound_DIV(b)% vec)

     ALLOCATE (bc_table(b) % bpoints(bc_table(b) % Np_bound)) 
     bc_table(b) % bpoints = p_bound_DIV(b) % vec

   ENDDO
   
   DEALLOCATE (bound_p_DIV, p_bound_DIV)


   ! Initializes the boundary condition table
   DO b = 1, SIZE(bc_table)

      Np_bound   = bc_table(b) % Np_bound
      value_type = bc_table(b) % value_type
      
      SELECT CASE (value_type)	    

   	! Temperature boundary values ------------------------------------------
   	CASE (BV__T_NONE)
   	  
	  ALLOCATE (bc_table(b) % bdata(1,1))
   	  bc_table(b) % bdata(1,1) = 0.d0 	     
   	    

   	CASE (BV__T_OO)
 
     	  ALLOCATE (bc_table(b) % bdata(1,1))
   	  bc_table(b) % bdata(1,1) = T__ww(ww_oo)
 

   	CASE (BV__T_CONST)
   	! Read in read_param_ns_boundary_cond
   	   
   	CASE (BV__T_0)

     	 ALLOCATE (bc_table(b) % bdata(1, Np_bound))	    
   	 DO i_ = 1, Np_bound
   	   i = bc_table(b) % bpoints(i_)
   	   bc_table(b) % bdata(1,i_) = T__ww(ww_b(:,i))
   	 ENDDO
	 !----------------------------------------------------------------------


   	! Temperature normal gradient boundary values --------------------------
   	CASE (BV__GT_NONE)

     	  ALLOCATE (bc_table(b) % bdata(1,1))
   	  bc_table(b) % bdata(1,1) = 0.d0

   	CASE (BV__GT_ADIABATIC)

   	  ALLOCATE (bc_table(b) % bdata(1,1))
   	  bc_table(b) % bdata(1,1) = 0.d0
	 !----------------------------------------------------------------------


   	CASE DEFAULT
   	WRITE(*,*) ' Boundary value of unknown type. STOP.'
   	WRITE(*,*) '   In MODULE ns_thermal_bc'
   	WRITE(*,*) '	  SUBROUTINE init_ns_thermal_bc'
   	STOP
      
      END SELECT
   	 
   ENDDO


   ! Used only in the subroutine NS_THERMAL_BOUNDARY_FLUX 
   C_ther = 1 / (REF % Re * REF % Pr)
   
   END SUBROUTINE  init_ns_thermal_bc 



   !============================================================ 
   FUNCTION Tb__ns_thermal_bc( T_b, normal, rr_b ) RESULT( Tb )
   !============================================================ 
     

      !------------------------------------------------------------ 
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  T_b
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  normal
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  rr_b

      REAL(KIND=8), DIMENSION(SIZE(T_b))        ::  Tb
      !------------------------------------------------------------ 
      INTEGER                                ::  value_type
      REAL(KIND=8), DIMENSION(:,:), POINTER  ::  bdata
      INTEGER                                ::  Np_bound
      INTEGER,      DIMENSION(:),   POINTER  ::  bpoints
      INTEGER  ::  i, i_, b

      REAL(KIND=8), DIMENSION(SIZE(rr_b,1),  SIZE(rr_b,2))   :: void_1
      REAL(KIND=8), DIMENSION(SIZE(normal,1),SIZE(normal,2)) :: void_2
      !------------------------------------------------------------ 

      ! Dummy assignation
      void_1 = rr_b
      void_2 = normal
     

      DO b = 1, SIZE(bc_table)
 
         value_type =  bc_table(b) % value_type
         bdata      => bc_table(b) % bdata
         Np_bound   =  bc_table(b) % Np_bound
         bpoints    => bc_table(b) % bpoints

         SELECT CASE (value_type)            
           
           !--------------
           CASE (BV__T_NONE)
           !--------------
           Tb(bpoints) = T_b(bpoints)
              
           !--------------
           CASE (BV__T_OO)
           !--------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              Tb(i) = bdata(1,1)              
           ENDDO

           !-----------------
           CASE (BV__T_CONST)
           !-----------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              Tb(i) = bdata(1,1)              
           ENDDO
 
           !--------------------------
           CASE (BV__T_0)
           !--------------------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              Tb(i) = bdata(1,i_) 
           ENDDO
	   
           !--------------
           CASE DEFAULT
           !--------------
           Tb(bpoints) = T_b(bpoints)

         END SELECT
      
      ENDDO
      

   END FUNCTION Tb__ns_thermal_bc
   !============================================================ 



   !============================================================ 
   FUNCTION G_Tb__ns_thermal_bc( G_T_b, normal, rr_b ) RESULT( G_Tb )
   !============================================================ 
     

      !------------------------------------------------------------ 
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  G_T_b
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  normal
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  rr_b

      REAL(KIND=8), DIMENSION( SIZE(G_T_b,1), &
                               SIZE(G_T_b,2)  )  ::  G_Tb
      !------------------------------------------------------------ 
      INTEGER                                ::  value_type
!      REAL(KIND=8), DIMENSION(:,:), POINTER  ::  bdata
      INTEGER                                ::  Np_bound
      INTEGER,      DIMENSION(:),   POINTER  ::  bpoints
      REAL(KIND=8), DIMENSION(SIZE(normal,1))  ::  nn
      INTEGER  ::  i, i_, b

      REAL(KIND=8), DIMENSION(SIZE(rr_b,1),  SIZE(rr_b,2))   :: void_1
      REAL(KIND=8), DIMENSION(SIZE(normal,1),SIZE(normal,2)) :: void_2
      !------------------------------------------------------------ 

      ! Dummy assignation
      void_1 = rr_b
      void_2 = normal

    
      DO b = 1, SIZE(bc_table)

         value_type =  bc_table(b) % value_type
!         bdata      => bc_table(b) % bdata
         Np_bound   =  bc_table(b) % Np_bound
         bpoints    => bc_table(b) % bpoints

         SELECT CASE (value_type)            

           !--------------
           CASE (BV__GT_NONE)
           !--------------
           G_Tb(:,bpoints) = G_T_b(:,bpoints)

           !--------------
           CASE (BV__GT_ADIABATIC)
           !--------------
           DO i_ = 1, Np_bound
              i = bpoints(i_)
              nn = normal(:,i) / SQRT( SUM(normal(:,i)**2) )
              G_Tb(:,i) = G_T_b(:,i) - SUM(G_T_b(:,i)*nn)*nn
           ENDDO
	   
           !--------------
           CASE DEFAULT
           !--------------
           G_Tb(:,bpoints) = G_T_b(:,bpoints)
                       
         END SELECT

      ENDDO


   END FUNCTION G_Tb__ns_thermal_bc
   !============================================================ 



   !============================================================ 
   SUBROUTINE ns_impose_strong_T_bc( ww, jd_jb, normal, rr_b )
   !============================================================ 
     

      !------------------------------------------------------------ 
      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT)  ::  ww
      INTEGER,      DIMENSION(:),   INTENT(IN)     ::  jd_jb
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  normal
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  rr_b
      !------------------------------------------------------------ 
      INTEGER                                ::  cond_type, value_type
      REAL(KIND=8), DIMENSION(:,:), POINTER  ::  bdata
      INTEGER                                ::  Np_bound
      INTEGER,      DIMENSION(:),   POINTER  ::  bpoints
      REAL(KIND=8)  ::  T_b
      INTEGER  ::  i, i_, b, p

      REAL(KIND=8), DIMENSION(SIZE(rr_b,1),  SIZE(rr_b,2))   :: void_1
      REAL(KIND=8), DIMENSION(SIZE(normal,1),SIZE(normal,2)) :: void_2
      !------------------------------------------------------------ 

      ! Dummy assignation
      void_1 = rr_b
      void_2 = normal
     
      p = SIZE(ww,1)
      
      DO b = 1, SIZE(bc_table)
 
         cond_type  =  bc_table(b) % cond_type
         value_type =  bc_table(b) % value_type
         bdata      => bc_table(b) % bdata
         Np_bound   =  bc_table(b) % Np_bound
         bpoints    => bc_table(b) % bpoints

         IF (cond_type == BC__STRONG) THEN

            SELECT CASE ( value_type )            

              !--------------
              CASE(BV__T_NONE)
              !--------------

              !--------------
              CASE(BV__T_OO)
              !--------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 T_b = bdata(1,1)              
                 ww(p,i) = Etot__rho_m_T (ww(1,i), ww(2:p-1,i), T_b) 
              ENDDO

              !-----------------
              CASE(BV__T_CONST)
              !-----------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 T_b = bdata(1,1)              
                 ww(p,i) = Etot__rho_m_T (ww(1,i), ww(2:p-1,i), T_b) 
              ENDDO

              !--------------------------
              CASE(BV__T_0)
              !--------------------------
              DO i_ = 1, Np_bound
                 i = jd_jb(bpoints(i_))
                 T_b = bdata(1,i_)              
                 ww(p,i) = Etot__rho_m_T (ww(1,i), ww(2:p-1,i), T_b) 
              ENDDO

            END SELECT

         ENDIF
      
      ENDDO
      

   END SUBROUTINE ns_impose_strong_T_bc
   !============================================================ 



END MODULE ns_thermal_bc

