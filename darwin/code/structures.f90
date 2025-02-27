MODULE structures

 USE dynamic_vector
   
! STRUCTURE TYPES DEFINITIONS:
    
  
  TYPE curve      
    INTEGER                               :: nd, ns
    REAL(KIND=8), DIMENSION(:),   POINTER :: s, h 
    REAL(KIND=8), DIMENSION(:,:), POINTER :: x, xs      
  END TYPE curve
  
  
        
  TYPE grid_type

    INTEGER  ::  adaptation_Level
    INTEGER  ::  k_d
    INTEGER  ::  Nb
    ! Grid nodes
    INTEGER                                ::  Nj_d, Nj_b
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  ::  rr
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  jd_jb
    INTEGER,      DIMENSION(:,:), ALLOCATABLE  ::  jb_jd
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  bound_p
    ! Domain elements
    INTEGER                                ::  Nm_d
    TYPE (D_I_V), DIMENSION(:),   ALLOCATABLE  ::  j_m_d
    TYPE (D_I_V), DIMENSION(:),   POINTER  ::  m_j_d
    TYPE (D_I_V), DIMENSION(:),   POINTER  ::  ma_m_d
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  ele_type_d
    ! Boundary elements
    INTEGER                                ::  Nm_b
    TYPE (D_I_V), DIMENSION(:),   ALLOCATABLE  ::  j_m_b
    TYPE (D_I_V), DIMENSION(:),   POINTER  ::  m_j_b
    TYPE (D_I_V), DIMENSION(:),   POINTER  ::  ma_m_b
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  ele_type_b
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  bound_m
    ! Grid egdes
    INTEGER                                ::  Nc_d, Nc_b
    INTEGER,      DIMENSION(:,:), ALLOCATABLE  ::  j_c
    TYPE (D_I_V), DIMENSION(:),   POINTER  ::  c_j
    TYPE (D_I_V), DIMENSION(:),   POINTER  ::  m_c
    TYPE (D_I_V), DIMENSION(:),   POINTER  ::  c_m
    INTEGER,      DIMENSION(:,:), ALLOCATABLE  ::  j_c_b
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  cd_cb
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  cb_cd
    INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  bound_c
    ! Metrics
    INTEGER,      DIMENSION(:,:),   ALLOCATABLE  ::  jcd_fem
    INTEGER,      DIMENSION(:,:),   ALLOCATABLE  ::  jcb_fem
    TYPE (D_I_V), DIMENSION(:),     ALLOCATABLE  ::  cjd_fem
    REAL(KIND=8), DIMENSION(:),     ALLOCATABLE  ::  cell
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  ::  eta
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  ::  chi_b
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  ::  xi_bp
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  ::  stiff_T
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  ::  Kb_ijs, Kb_ija

    ! Boundary Curve Coordinates
    REAL(KIND=8), DIMENSION(:),   POINTER  ::  curve_coordinates
    ! Boundary Geometry
    TYPE(curve),  DIMENSION(:),   POINTER  ::  line

  END TYPE grid_type  
 

 
 
 
  TYPE E_R_P   ! Edge_Refinement_Pattern_Type

    ! Refinement pattern from grid m to grid m+1
    INTEGER N_refined_edges   ! Number of inserted nodes
 
    ! Indices i of node inserted in the c-th refined edge
    ! (according to grid m+1 numbering), with 1 < c < N_refined_edges
    !                      c
    INTEGER,     DIMENSION(:),   POINTER  :: inserted_nodes

    ! Indices j (1) and k (2) of the two nodes belonging to the refined
    ! edge c
    !                      2 c
    INTEGER,     DIMENSION(:,:), POINTER  :: edge_nodes 
    TYPE(D_I_V), DIMENSION(:),   POINTER  :: adjacent_nodes

  END TYPE E_R_P  





  TYPE :: R_P
 
    INTEGER                                 :: e_type
    LOGICAL                                 :: ref_flag
    INTEGER                                 :: ref_np_nmr
    INTEGER,          DIMENSION(:), POINTER :: np_idx
    CHARACTER(LEN=3)                        :: strategy
 
  END TYPE
   


  TYPE estimator_type

   INTEGER :: error_variable
   INTEGER :: estimator_function
   INTEGER :: passages
   
   REAL(KIND=8), DIMENSION(:,:), POINTER :: errors
   REAL(KIND=8), DIMENSION(:,:), POINTER :: R_Thres
   REAL(KIND=8), DIMENSION(:,:), POINTER :: C_Thres
    
  END TYPE estimator_type



  TYPE solution_type

    ! Number of spatial dimensions and number of points
    INTEGER :: k_d, Nj_d
    ! Solution format
    INTEGER           :: sol_fmt
    ! Solution
    REAL(KIND=8), DIMENSION(:,:), POINTER :: ww
    
  END TYPE solution_type
  

  

  TYPE crv_crd
    REAL(KIND=8), DIMENSION(:),   POINTER :: points
  END TYPE crv_crd 

  
   
  TYPE cnv_idx
    INTEGER,      DIMENSION(:),   POINTER :: index
  END TYPE cnv_idx 
 
 

  TYPE adaption_param

    ! Problem Parameters
    CHARACTER(LEN=30) :: prob_name 
    CHARACTER(LEN=30) :: grid_name
    INTEGER           :: sol_fmt

    ! Type of adaption and boxes
    INTEGER :: ref_type, nbbox
    TYPE(D_I_V_R), DIMENSION(:), POINTER :: box

    ! Adaption Levels
    INTEGER :: number_adaptation_levels
    INTEGER :: entry_level

    ! Estimator List
    INTEGER                                     :: N_estimators
    TYPE(estimator_type), DIMENSION(:), POINTER :: estimator

    ! Boundary Treatment
    INTEGER                        :: Nb_lines
    LOGICAL, DIMENSION(:), POINTER :: solid_bou
    INTEGER                        :: bou_tre_typ
    INTEGER                        :: jou_mod_fml
    REAL(KIND=8)                   :: jmf_const
    REAL(KIND=8)                   :: jmf_alpha
    REAL(KIND=8)                   :: Poi_Coe



    ! Quality Parameters    
    LOGICAL  :: swap, smooth
    INTEGER  :: l_its

    INTEGER :: max_number_QL
    REAL(KIND=8) :: min_angle, min_length, relax_1, relax_2

  END TYPE adaption_param


 TYPE(grid_type) :: grid_zero 

 INTEGER, PARAMETER :: PIG_EU     = 0, &  ! PIG         Euler
                       PIG_NSL    = 1, &  ! PIG         Navier Stokes Laminar
                       PIG_NST    = 2, &  ! PIG         Navier Stokes Turbulent
                       VdW_EU     = 3, &  ! Vd Waals    Euler
                       VdW_NSL    = 4, &  ! Vd Waals    Navier Stokes Laminar
                       VdW_NST    = 5, &  ! Vd Waals    Navier Stokes Turbulent
                       FENSAP = 6, &
                       FENSAPDROP = 7, &
                       FENSAPDROPICE = 8, &
                       SU2_EU    =  9, &  ! SU2 Euler solution
                       SU2_NSL   = 10, &  ! SU2 Laminar Navier-Stokes solution 
                       SU2_NSSA  = 11, &  ! SU2 Turbulent Navier-Stokes solution, SA model
                       SU2_NSSST = 12     ! SU2 Turbulent Navier-Stokes solution, SST model
 
END MODULE structures
