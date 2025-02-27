MODULE error_estimator
 
 USE structures   
 USE derivatives
 USE thermodynamic
 USE grid_utils
 USE io_proc
 
 !-------------------------------------------------------------------------------
 IMPLICIT NONE
 
 INTEGER, PARAMETER   ::  ERR_VARIABLE_RHO         = 1, &
                          ERR_VARIABLE_MACH        = 2, &
                          ERR_VARIABLE_PRESSURE    = 3, &
                          ERR_VARIABLE_TEMPERATURE = 4, &
                          ERR_VARIABLE_VORTICITY   = 5, &
                          ERR_VARIABLE_KNUDSEN     = 6, &
                          ERR_VARIABLE_XVELOCITY   = 7, &
                          ERR_VARIABLE_YVELOCITY   = 8, &
                          ERR_VARIABLE_TURBVISC    = 9, &
                          ERR_VARIABLE_LWC         = 10, &
                          ERR_VARIABLE_ICEGROWTH   = 11
                          
                                                  
 INTEGER, PARAMETER, PUBLIC  :: ERR_ESTIMATOR_GRADIENT      = 1, &  ! Scalar
                                ERR_ESTIMATOR_II_DERIVATIVE = 2, &  ! Vector
                                ERR_ESTIMATOR_WEBSTER       = 3, &  ! Vector
                                ERR_ESTIMATOR_ANISOTROPIC   = 4, &  ! Scalar
                                ERR_ESTIMATOR_METRIC_BASED  = 5     ! Scalar
 !-------------------------------------------------------------------------------

 CONTAINS
   
 
 
 SUBROUTINE  Compute_Error_Estimator( grid, grid_name, solution, estimator, mode )
 !------------------------------------------------------------------------------
 USE np_quadrature, ONLY: np_quadr_w_ROTu

 IMPLICIT NONE

 TYPE(grid_type),     INTENT(IN) :: grid
 CHARACTER(LEN=30) :: grid_name
 TYPE(solution_type), INTENT(IN) :: solution

 TYPE(estimator_type), DIMENSION(:), INTENT(INOUT) :: estimator
 integer, intent(in) :: mode

 REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: key_variable

 REAL(KIND=8), DIMENSION(grid % k_d, grid % Nj_d) :: V

 REAL(KIND=8), DIMENSION(3, grid % Nj_d) :: omega

 INTEGER, DIMENSION(:), ALLOCATABLE :: iB, mB
 INTEGER :: i, j, j_, p, i_, Npb, Neb
 !------------------------------------------------------------------------------

  ALLOCATE (key_variable(grid % Nj_d))
  
  DO j = 1, SIZE(estimator)

   ! Compute variable that will drive adaption
   SELECT CASE (estimator(j)%error_variable)

    CASE (ERR_VARIABLE_RHO)

       key_variable = solution%ww(1,:)

    CASE (ERR_VARIABLE_MACH)

       SELECT CASE (solution%sol_fmt) 
          ! CASE(SU2_EU, SU2_NSL, SU2_NSSA)
          !    key_variable = solution % ww(10,:)
          CASE(SU2_EU)
             key_variable = solution % ww(6+grid%k_d,:)
          CASE(SU2_NSL)
             key_variable = solution % ww(6+grid%k_d,:)
          CASE(SU2_NSSA)
             key_variable = solution % ww(7+grid%k_d,:)
          CASE(SU2_NSSST)
             key_variable = solution % ww(8+grid%k_d,:)
          CASE DEFAULT
             DO i = 1, grid % Nj_d  
               key_variable(i) = Mach__ww(grid % k_d+2, solution % ww(:,i)) 
             ENDDO
       END SELECT

    CASE (ERR_VARIABLE_PRESSURE)

       SELECT CASE (solution%sol_fmt) 
          ! CASE(SU2_EU, SU2_NSL, SU2_NSSA)
          !    key_variable = solution % ww(7,:)
          CASE(SU2_EU)
             key_variable = solution % ww(3+grid%k_d,:)
          CASE(SU2_NSL)
             key_variable = solution % ww(3+grid%k_d,:)
          CASE(SU2_NSSA)
             key_variable = solution % ww(4+grid%k_d,:)
          CASE(SU2_NSSST)
             key_variable = solution % ww(5+grid%k_d,:)
          CASE DEFAULT
             DO i = 1, grid % Nj_d  
               key_variable(i) = Pressure__ww(grid % k_d+2, solution % ww(:,i)) 
             ENDDO    
       END SELECT


    CASE (ERR_VARIABLE_TEMPERATURE)

       SELECT CASE (solution%sol_fmt) 
          ! CASE(SU2_EU, SU2_NSL, SU2_NSSA)
          !    key_variable = solution % ww(8,:)
          CASE(SU2_EU)
             key_variable = solution % ww(4+grid%k_d,:)
          CASE(SU2_NSL)
             key_variable = solution % ww(4+grid%k_d,:)
          CASE(SU2_NSSA)
             key_variable = solution % ww(5+grid%k_d,:)
          CASE(SU2_NSSST)
             key_variable = solution % ww(6+grid%k_d,:)
          CASE DEFAULT
             DO i = 1, grid % Nj_d  
               key_variable(i) = Temperature__ww(grid % k_d+2, solution % ww(:,i)) 
             ENDDO    
       END SELECT

    CASE (ERR_VARIABLE_VORTICITY)

       ! Velocity vector
       DO i = 1, grid % Nj_d
         V(:,i) = solution % ww(2:1+grid % k_d,i) / solution % ww(1,i)
       ENDDO

       ! Vorticity omega = curl(V)
       omega = np_quadr_w_ROTu(grid % jcd_fem, grid % jcb_fem, grid % jd_jb,  &
                               grid % eta,     grid % chi_b,   grid % xi_bp,  V)

       DO i = 1, grid % Nj_d
         omega(:,i) = omega(:,i) / grid % cell(i)
       ENDDO


       !-----------------------------------------------------------------------------------------------
       OPEN (UNIT = 11, FILE = 'vorticity.plt', STATUS = 'unknown')
  
       WRITE(11,*) 'TITLE = "Vorticity"'
       WRITE(11,*) 'VARIABLES = "X", "Y", "u", "v", "omega_x", "omega_y", "omega_z" "ABS(omega)"'
  
       WRITE(11,*) 'ZONE T = "Domain", N=', grid % Nj_d, 'E=', grid % Nm_d,  &
                   'DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
  
       DO i = 1, grid % Nj_d
         WRITE(11,*) grid % rr(:,i), V(:,i), omega(:, i),  SQRT(SUM(omega(:,i)**2))
       ENDDO
  
       DO i = 1, grid % Nm_d
         IF (grid % ele_type_d(i) == 2)  WRITE(11,*) grid % j_m_d(i) % vec, grid % j_m_d(i) % vec(3) 
         IF (grid % ele_type_d(i) == 3)  WRITE(11,*) grid % j_m_d(i) % vec
       ENDDO 
         
       ! Boundaries
!        DO i = 1, MAXVAL(grid % bound_p)
!        
!          Npb = COUNT(bound_p == i)
!        Neb = COUNT(bound_m == i)
!        
!          ALLOCATE (iB(Npb), mB(Neb))
!        
!          p = 1
!          DO j_ = 1, grid % Nj_b
!         IF (bound_p(j_) == i) THEN
!           iB(p) = j_
!           p = p+1
!         ENDIF
!        ENDDO
!        
!          p = 1
!          DO j_ = 1, grid % Nm_b
!         IF (bound_m(j_) == i) THEN
!           mB(p) = j_
!           p = p+1
!         ENDIF
!        ENDDO
!        
!        WRITE(11,*)
!          WRITE(11,*) 'ZONE T = "Boundary", N=', Npb, 'E=', Neb,  &
!                    'DATAPACKING=POINT, ZONETYPE=FELINESEG'
!   
!          DO i_ = 1, Npb
!          p = jd_jb(iB(i_))
!            WRITE(11,*) grid % rr(:,p), V(:,i), omega(:, p),  SQRT(SUM(omega(:,p)**2))
!          ENDDO
!   
!          p = 1
!          DO i_ = 1, Neb
!            WRITE(11,*) p, p+1
!          p = p+1
!          ENDDO   
!        
!          DEALLOCATE (iB, mB)
!        
!        ENDDO 
       
       CLOSE(11)
       !-----------------------------------------------------------------------------------------------
  



       DO i = 1, grid % Nj_d
         key_variable(i) = SQRT(SUM(omega(:,i)**2))
       ENDDO    

    
    CASE (ERR_VARIABLE_KNUDSEN)
      
       key_variable = Knudsen__ww(grid, solution % ww) 
       
    CASE (ERR_VARIABLE_XVELOCITY)
       key_variable = solution%ww(2,:) / solution%ww(1,:)
    
    CASE (ERR_VARIABLE_YVELOCITY)
       key_variable = solution%ww(3,:) / solution%ww(1,:)
    
    CASE (ERR_VARIABLE_TURBVISC )
       SELECT CASE (solution%sol_fmt) 
          ! CASE(SU2_EU, SU2_NSL, SU2_NSSA)
          !    key_variable = solution % ww(15,:)
          CASE(SU2_NSSA)
             key_variable = solution % ww(13+grid%k_d,:)
          CASE(SU2_NSSST)
             key_variable = solution % ww(14+grid%k_d,:)
          CASE DEFAULT
             key_variable = solution%ww(5,:)
       END SELECT
    
    CASE (ERR_VARIABLE_LWC      )
       key_variable = solution%ww(6,:)
    
    CASE (ERR_VARIABLE_ICEGROWTH)
       key_variable = solution%ww(8,:)


   END SELECT 


   ! Compute error
   SELECT CASE (estimator(j) % estimator_function)

    CASE (ERR_ESTIMATOR_GRADIENT)

       ALLOCATE ( estimator(j)%errors(1,grid%Nc_d) )
       estimator(j)%errors = Gradient_estimator( grid, key_variable )          

    CASE (ERR_ESTIMATOR_II_DERIVATIVE)

       ALLOCATE ( estimator(j)%errors(2,grid%Nc_d) )
       estimator(j)%errors = II_derivative_estimator( grid, key_variable, solution )

    CASE (ERR_ESTIMATOR_WEBSTER)

       ALLOCATE ( estimator(j)%errors(2,grid%Nc_d) )
       estimator(j)%errors = Webster_estimator( grid, key_variable, solution )

    CASE (ERR_ESTIMATOR_ANISOTROPIC)

       ALLOCATE ( estimator(j)%errors(1,grid%Nc_d) )
       estimator(j)%errors = Anisotropic_estimator( grid, key_variable, solution )        

    CASE (ERR_ESTIMATOR_METRIC_BASED)

       ALLOCATE ( estimator(j)%errors(1,grid%Nc_d) )
       estimator(j)%errors = Metric_based_estimator( grid, grid_name, key_variable, mode )        

   END SELECT 
  
  ENDDO   

  DEALLOCATE (key_variable)

 END SUBROUTINE  compute_error_estimator 






 
 FUNCTION  Gradient_estimator(grid, key_variable)   RESULT(edge_error)
 !---------------------------------------------------------------------------
  IMPLICIT NONE 
  
  REAL(KIND=8),    DIMENSION(:),   INTENT(IN)  :: key_variable
  TYPE(grid_type),                 INTENT(IN)  :: grid   
  
  INTEGER                                      :: i, nd1, nd2
  REAL(KIND=8)                                 :: length    
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: nodal_error 
  REAL(KIND=8),    DIMENSION(:,:), POINTER     :: edge_error
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: G_Kv
 !---------------------------------------------------------------------------
 
  ALLOCATE (nodal_error(1, grid % Nj_d), &
             edge_error(1, grid % Nc_d), &
                   G_Kv(2, grid % Nj_d))

  CALL Gradient(grid, key_variable, G_Kv)

  ! Nodal error  
  DO i = 1, grid % Nj_d 
   nodal_error(1,i) = SQRT(G_Kv(1,i)**2 + G_Kv(2,i)**2) 
  ENDDO   
  
  ! Edge error
  DO i = 1, grid % Nc_d
  
    nd1 = grid % j_c(1,i)
    nd2 = grid % j_c(2,i)    

    length = SQRT( (grid % rr(2, nd2) - grid % rr(2, nd1))**2  +  (grid % rr(1, nd2) - grid % rr(1, nd1))**2 )

    edge_error(1,i) = length * MAX(nodal_error(1, nd1), nodal_error(1, nd2))

  ENDDO   

  DEALLOCATE (nodal_error, G_Kv)
 
 END FUNCTION Gradient_estimator
 





 FUNCTION II_derivative_estimator( grid, key_variable, solution ) &
 RESULT(edge_error)
 !---------------------------------------------------------------------------
 IMPLICIT NONE 
 
 REAL(KIND=8),    DIMENSION(:),     INTENT(IN)  :: key_variable
 TYPE(grid_type),                   INTENT(IN)  :: grid 
 TYPE(solution_type),               INTENT(IN)  :: solution    
 
 INTEGER                                        :: i, j, k 
 INTEGER                                        :: nd1, nd2
 REAL(KIND=8)                                   :: u, v, w, G_Kv_mod
 REAL(KIND=8)                                   :: equiv_rad
 REAL(KIND=8),    DIMENSION(2)                  :: Tg, Nm
 REAL(KIND=8)                                   :: ddMT, ddMN  
 REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: nodal_error 
 REAL(KIND=8),    DIMENSION(:,:),   POINTER     :: edge_error
 REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: G_Kv
 REAL(KIND=8),    DIMENSION(:,:,:), ALLOCATABLE :: H_Kv
 !---------------------------------------------------------------------------

  ALLOCATE( nodal_error(2,grid%Nj_d) )
  ALLOCATE( edge_error(2,grid%Nc_d) )
  
  ALLOCATE( G_Kv(2,grid%Nj_d) )
  CALL gradient( grid, key_variable, G_Kv )        
  
  ALLOCATE( H_Kv(grid%k_d,2,grid%Nj_d) )
  CALL hessian( grid, key_variable, H_Kv )

  ! Nodal error  
  DO j = 1, grid % Nj_d

    u = solution%ww(2,j)/solution%ww(1,j)
    v = solution%ww(3,j)/solution%ww(1,j)
    w = SQRT( u**2 + v**2 )
   
   IF ( (w .EQ. 0.d0) .AND. (grid%jb_jd(1,j) .NE. 0) ) THEN   ! Wall Boundary Node
    
    Tg = Compute_Tangent_Vector(j, grid) ! Vector Tangent to boundary wall in node 'j'    
        
    Nm(1) = -Tg(2)         
    Nm(2) =  Tg(1)     
   
   ELSE  ! Node not belonging to a boundary wall     
    
    IF ( w .EQ. 0.d0 ) THEN
     
     !PRINT*, 'WARNING. ERROR ESTIMATOR:'
     !PRINT*, 'Velocity  is  zero somewhere inside the domain.'
     !PRINT*, 'Considering GRADIENT VECTOR instead of VELOCITY'
     
     G_Kv_mod = SQRT( G_Kv(1,j)**2 + G_Kv(2,j)**2 )
     
     Tg(1) = G_Kv(1,j) / G_Kv_mod
     Tg(2) = G_Kv(2,j) / G_Kv_mod
          
     Nm(1) = -Tg(2)
     Nm(2) =  Tg(1)
     
    ELSE

     Tg(1) =  u / w
     Tg(2) =  v / w
         
     Nm(1) = -Tg(2)   
     Nm(2) =  Tg(1)        
    ENDIF
     
   ENDIF
   
    equiv_rad = SQRT( grid%cell(j) )
    
    ddMT = 0.d0;  ddMN = 0.d0   
    DO i = 1, 2
      DO k = 1, 2
       ddMT = ddMT + Tg(i)*Tg(k)*H_Kv(i,k,j)
       ddMN = ddMN + Nm(i)*Nm(k)*H_Kv(i,k,j)
      ENDDO
    ENDDO

    nodal_error(1,j) = (equiv_rad**2)*ABS(ddMT)
    nodal_error(2,j) = (equiv_rad**2)*ABS(ddMN)
  
  ENDDO


  ! Edge error
  DO i = 1, grid%Nc_d
    nd1 = grid%j_c(1,i)
    nd2 = grid%j_c(2,i)
    DO j = 1, SIZE(edge_error,1)
      edge_error(j,i) = MAX( nodal_error(j,nd1), nodal_error(j,nd2) )    
    ENDDO
  ENDDO   

  DEALLOCATE (nodal_error, G_Kv, H_Kv)
  
 END FUNCTION II_derivative_estimator







 FUNCTION  Webster_estimator( grid, key_variable, solution )  RESULT(edge_error)
 !-----------------------------------------------------------------------------------
 IMPLICIT NONE 
 
 REAL(KIND=8),    DIMENSION(:),     INTENT(IN)  :: key_variable
 TYPE(grid_type),                   INTENT(IN)  :: grid 
 TYPE(solution_type),               INTENT(IN)  :: solution    
 
 INTEGER                                        :: i, j, k 
 INTEGER                                        :: nd1, nd2
 REAL(KIND=8)                                   :: u, v, w 
 REAL(KIND=8)                                   :: equiv_rad, eps, mean 
 REAL(KIND=8),    DIMENSION(2)                  :: Tg, Nm
 REAL(KIND=8)                                   :: dMT, dMN  
 REAL(KIND=8)                                   :: ddMT, ddMN  
 REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: nodal_error 
 REAL(KIND=8),    DIMENSION(:,:),   POINTER     :: edge_error
 REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: G_Kv
 REAL(KIND=8),    DIMENSION(:,:,:), ALLOCATABLE :: H_Kv
 !------------------------------------------------------------------------------------

  ALLOCATE( nodal_error(2,grid%Nj_d) )
  ALLOCATE( edge_error(2,grid%Nc_d) )
  
  ALLOCATE( G_Kv(2,grid%Nj_d) )
  CALL gradient( grid, key_variable, G_Kv )
  
  ALLOCATE( H_Kv(grid%k_d,2,grid%Nj_d) )
  CALL hessian( grid, key_variable, H_Kv )
  
  ! Nodal error
  DO j = 1, grid%Nj_d
  
   u = solution%ww(2,j)/solution%ww(1,j)
   v = solution%ww(3,j)/solution%ww(1,j)
   w = SQRT( u**2 + v**2 )  

   Tg(1) =  u / (w+1.E-6);   Tg(2) = v / (w+1.E-6)    
   Nm(1) = -Tg(2);   Nm(2) =  Tg(1)

   eps = 0.12
   equiv_rad = SQRT( grid%cell(j) )   
   mean = Average(key_variable)

   dMT = ( Tg(1)*G_Kv(1,j) + Tg(2)*G_Kv(2,j) )
   dMN = ( Nm(1)*G_Kv(1,j) + Nm(2)*G_Kv(2,j) )

   ddMT = 0.d0;  ddMN = 0.d0    
   DO i = 1, 2
    DO k = 1, 2
      ddMT = ddMT + Tg(i)*Tg(k)*H_Kv(i,k,j)
      ddMN = ddMN + Nm(i)*Nm(k)*H_Kv(i,k,j)
    ENDDO
   ENDDO

   nodal_error(1,j) = (equiv_rad**2)*ABS(ddMT) / ( equiv_rad*ABS(dMT) + eps*ABS(mean) )
   nodal_error(2,j) = (equiv_rad**2)*ABS(ddMN) / ( equiv_rad*ABS(dMN) + eps*ABS(mean) )
  
  ENDDO

  ! Edge error
  DO i = 1, grid%Nc_d
    nd1 = grid%j_c(1,i)
    nd2 = grid%j_c(2,i)      
    DO j = 1, SIZE(edge_error,1)
      edge_error(j,i) = MAX( nodal_error(j,nd1), nodal_error(j,nd2) )    
    ENDDO
  ENDDO   

  DEALLOCATE (nodal_error, G_Kv, H_Kv)

 END FUNCTION Webster_estimator







 FUNCTION Anisotropic_estimator( grid, key_variable, solution ) &
                                                         RESULT(edge_error)
 !---------------------------------------------------------------------------
  IMPLICIT NONE 
  
  REAL(KIND=8),    DIMENSION(:),   INTENT(IN)  :: key_variable
  TYPE(grid_type),                 INTENT(IN)  :: grid
  TYPE(solution_type),             INTENT(IN)  :: solution     
  
  INTEGER                                      :: i, nd1, nd2, j, k
  REAL(KIND=8),    DIMENSION(2)                :: E_ij, V_ij
  REAL(KIND=8)                                 :: u, v, w 
  REAL(KIND=8)                                 :: equiv_rad
  REAL(KIND=8),    DIMENSION(2)                :: Tg, Nm
  REAL(KIND=8)                                 :: ddMT, ddMN    
  REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: nodal_error 
  REAL(KIND=8),    DIMENSION(:,:),   POINTER     :: edge_error
  REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: G_Kv
  REAL(KIND=8),    DIMENSION(:,:,:), ALLOCATABLE :: H_Kv  
 !---------------------------------------------------------------------------
 
  ALLOCATE( nodal_error(2,grid%Nj_d) )
  ALLOCATE( edge_error(1,grid%Nc_d) )
  
  ALLOCATE( G_Kv(2,grid%Nj_d) )  
  CALL gradient( grid, key_variable, G_Kv )

  ALLOCATE( H_Kv(grid%k_d,2,grid%Nj_d) )
  CALL hessian( grid, key_variable, H_Kv )

  ! Nodal error  
  DO j = 1, grid%Nj_d
  
    u = solution%ww(2,j)/solution%ww(1,j)
    v = solution%ww(3,j)/solution%ww(1,j)
    w = SQRT( u**2 + v**2 )  

    Tg(1) =  u / w;   Tg(2) = v / w    
    Nm(1) = -Tg(2);   Nm(2) =  Tg(1)

    equiv_rad = SQRT( grid%cell(j) )

    ddMT = 0.d0;  ddMN = 0.d0   
    DO i = 1, 2
      DO k = 1, 2
       ddMT = ddMT + Tg(i)*Tg(k)*H_Kv(i,k,j)
       ddMN = ddMN + Nm(i)*Nm(k)*H_Kv(i,k,j)
      ENDDO
    ENDDO

    nodal_error(1,j) = (equiv_rad**2)*ABS(ddMT)
    nodal_error(2,j) = (equiv_rad**2)*ABS(ddMN)
  
  ENDDO
  
  
  ! Edge error
  DO i = 1, grid%Nc_d
  
    nd1 = grid%j_c(1,i)
    nd2 = grid%j_c(2,i)

    V_ij(1) = grid%rr(1,nd1) - grid%rr(1,nd2)  ! x component of vector linking nodes (edge)
    V_ij(2) = grid%rr(2,nd1) - grid%rr(2,nd2)  ! y component of vector linking nodes (edge)
        
    E_ij(1) = MAX( nodal_error(1,nd1), nodal_error(1,nd2) )
    E_ij(2) = MAX( nodal_error(2,nd1), nodal_error(2,nd2) )

    edge_error(1,i) = ABS( E_ij(1)*V_ij(1) + E_ij(2)*V_ij(2) ) 
        
  ENDDO   

  DEALLOCATE (nodal_error, G_Kv, H_Kv)
 
 END FUNCTION Anisotropic_estimator





  FUNCTION  Metric_based_estimator ( grid, grid_name, key_variable, mode ) &
    RESULT(edge_error)
  !---------------------------------------------------------------------------
  IMPLICIT NONE 
 
  TYPE(grid_type),               INTENT(IN)  :: grid   
  CHARACTER(LEN=30) :: grid_name 
  REAL(KIND=8),    DIMENSION(:), INTENT(IN)  :: key_variable
  integer, intent(in) :: mode
  
  REAL(KIND=8),    DIMENSION(:),     ALLOCATABLE :: lambda, vector_1, vector_2
  REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: R, D, inv_R
  REAL(KIND=8),    DIMENSION(:,:),   POINTER     :: edge_error
!   REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: G_Kv
  REAL(KIND=8),    DIMENSION(:,:,:), ALLOCATABLE :: H_Kv
  REAL(KIND=8),    DIMENSION(:,:,:), ALLOCATABLE :: M  

  REAL(KIND=8),    DIMENSION(2) :: gamma_p

  REAL(KIND=8) :: h11, h12, h21, h22, l1, l2
  REAL(KIND=8) :: gamma_pT__M__gamma_p 

  INTEGER :: i, j, k
  INTEGER :: nd1, nd2
  !---------------------------------------------------------------------------
 
  ALLOCATE (edge_error(1,grid % Nc_d))

!   ALLOCATE (G_Kv(2,grid % Nj_d))
!   CALL gradient(grid, key_variable, G_Kv)

  ALLOCATE (H_Kv(grid % k_d, grid % k_d, grid % Nj_d))
  CALL hessian(grid, key_variable, H_Kv)

  ! Computing Metric tensor in each node
  ALLOCATE (M(2,2,grid % Nj_d))
  ALLOCATE (D(2,2), R(2,2), inv_R(2,2))
  ALLOCATE (lambda(2), vector_1(2), vector_2(2))

  DO k = 1, grid % Nj_d
   
    h11 = H_Kv(1,1,k);   h12 = H_Kv(1,2,k)
    h21 = H_Kv(2,1,k);   h22 = H_Kv(2,2,k)         
   
    lambda(1) = ( ( h11+h22 ) + SQRT( (h11+h22)**2 - 4*(h11*h22-h12*h21) ) ) / 2     
    vector_1(1) = 1.d0
    vector_1(2) = ( lambda(1) - h11 )*vector_1(1) / h12
   
    lambda(2) = ( ( h11+h22 ) - SQRT( (h11+h22)**2 - 4*(h11*h22-h12*h21) ) ) / 2
    vector_2(1) = 1.d0
    vector_2(2) = ( lambda(2) - h11 )*vector_2(1)/h12

    R(:,1) = vector_1
    R(:,2) = vector_2


    IF (ALLOCATED(inv_R))   DEALLOCATE(inv_R)

    ALLOCATE (inv_R(SIZE(R,1), SIZE(R,2)))

    CALL Invert_Matrix(R, inv_R)

   
    D(1,1) = ABS( lambda(1) );  D(1,2) = 0.d0
    D(2,1) = 0.d0;              D(2,2) = ABS( lambda(2) )
   
    M(:,:,k) = MATMUL( R,D )
    M(:,:,k) = MATMUL( M(:,:,k),inv_R )
   
  ENDDO
  
 
  ! Edge length in the specified tensor metric 
  DO i = 1, grid%Nc_d
  
    nd1 = grid%j_c(1,i)
    nd2 = grid%j_c(2,i)

    gamma_p(1) = grid%rr(1,nd2) - grid%rr(1,nd1)
    gamma_p(2) = grid%rr(2,nd2) - grid%rr(2,nd1)



    gamma_pT__M__gamma_p = (SUM( gamma_p*MATMUL(M(:,:,nd1), gamma_p) ) + 1.e-7)
    l1 = SQRT( gamma_pT__M__gamma_p )
        
    gamma_pT__M__gamma_p = (SUM( gamma_p*MATMUL(M(:,:,nd2), gamma_p) ) + 1.e-7)
    l2 = SQRT( gamma_pT__M__gamma_p )

    edge_error(1,i) = (2.d0/3.d0)*ABS( (l1**2+l1*l2+l2**2) / (l1+l2) )
        
  ENDDO   

  DEALLOCATE (M, H_Kv)


 CONTAINS 
 
 
   SUBROUTINE Invert_Matrix(matrix,  inverse) 
   !---------------------------------------------------------------------------
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: matrix

   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: inverse
    
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Id, matrix_i 

   REAL(KIND=8) :: det_m, det_im
   !---------------------------------------------------------------------------
      
    ALLOCATE (       Id( SIZE(matrix,1),SIZE(matrix,2)),  &
               matrix_i( SIZE(matrix,1),SIZE(matrix,2)) )
  
    Id = 0.d0
    DO i = 1, SIZE(matrix,1)
     Id(i,i) = 1.d0
    ENDDO
    
    det_m = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
     
    DO i = 1, SIZE(matrix,2)    
     DO j = 1, SIZE(matrix,1)

      matrix_i = matrix
      matrix_i(:,j) = Id(:,i)
              
      det_im = matrix_i(1,1)*matrix_i(2,2) - matrix_i(1,2)*matrix_i(2,1)
        
      inverse(j,i) = det_im / det_m
     
     ENDDO
    ENDDO
 
    DEALLOCATE (Id, matrix_i)
 
   END SUBROUTINE Invert_Matrix
 
 END FUNCTION  Metric_based_estimator
 




 FUNCTION Average(vector)   RESULT(mean_value) 
 !---------------------------------------------------------------------------
 IMPLICIT NONE

 REAL(KIND=8), DIMENSION(:), INTENT(IN) :: vector
 REAL(KIND=8)                           :: mean_value
 !---------------------------------------------------------------------------
 
  mean_value = SUM(vector)
  mean_value = mean_value / (1.d0*SIZE(vector)) 
 
 END FUNCTION Average

 

END MODULE error_estimator
