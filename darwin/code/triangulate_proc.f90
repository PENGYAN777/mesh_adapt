MODULE triangulate_proc
 
 USE structures
 USE nodes
 USE mesh_structure
 USE grid_utils
 USE marking_proc
 USE boundary_geometry
 USE new_nodes
 USE new_elements
 
 IMPLICIT NONE

 CONTAINS
 
 ! Triangulate accepts as input 'mesh', then it modifies 'mesh' according to the
 ! refinement informations provided by np_refine and bnp_refine.
 
 SUBROUTINE Triangulate(np_refine,    bnp_refine,       &
                        nmr_p_to_add, bnmr_p_to_add,    &
                        mesh, ref_history, bou_tre_typ, &
                        displ, sol_ww, add_ww)

  !-------------------------------------------------------------------------------------------
  IMPLICIT NONE

  LOGICAL,           DIMENSION(:),   INTENT(IN)    :: np_refine, bnp_refine
  INTEGER                                            :: nmr_p_to_add, bnmr_p_to_add
  TYPE(grid_type),                   INTENT(INOUT) :: mesh
  TYPE(grid_type)                                  :: old_mesh
  TYPE(E_R_P),                       INTENT(IN)    :: ref_history
  INTEGER,                           INTENT(IN)    :: bou_tre_typ
  TYPE(D_I_V_R),     DIMENSION(:),   POINTER       :: displ
  REAL(KIND=8),      DIMENSION(:,:), INTENT(IN)    :: sol_ww
  REAL(KIND=8),      DIMENSION(:,:), POINTER       :: add_ww

  
  LOGICAL,          DIMENSION(:),   ALLOCATABLE          :: el_refine, bel_refine

  INTEGER                                                :: i
  INTEGER                                                :: dnmr_e_to_add, &
                                                            bnmr_e_to_add
  INTEGER                                                :: crop, fNp_d 
  INTEGER                                                :: rNp_d, rNe_d 
  INTEGER                                                :: rNp_b, rNe_b
  INTEGER                                                :: nds_nn   
  INTEGER,          DIMENSION(:),   ALLOCATABLE          :: rj_b
  INTEGER,          DIMENSION(:,:), ALLOCATABLE          :: rj_c_d
  INTEGER,          DIMENSION(:,:), ALLOCATABLE          :: rj_c_b 

  REAL(KIND=8),     DIMENSION(:,:), ALLOCATABLE          :: updt_cc
  REAL(KIND=8),     DIMENSION(:,:), ALLOCATABLE          :: aux_rr
  REAL(KIND=8),     DIMENSION(:,:), ALLOCATABLE          :: COPY_add_ww

  TYPE(crv_crd),    DIMENSION(:),   ALLOCATABLE          :: s_new

  TYPE(cnv_idx),    DIMENSION(:),   ALLOCATABLE          :: jd_jsb

  TYPE(R_P),        DIMENSION(:),   ALLOCATABLE          ::   pattern
  TYPE(R_P),        DIMENSION(:),   ALLOCATABLE          :: b_pattern 
  
  CHARACTER(LEN=1)                                       :: zone
  !-------------------------------------------------------------------------------------------


  ! Arrays allocation
  ALLOCATE ( el_refine(mesh%Nm_d))  ! Domain Refine Elements flag
  ALLOCATE (bel_refine(mesh%Nm_b))  ! Boundary Refine Elements flag
  
  ALLOCATE (rj_c_d(mesh%Nc_d,2))    ! For each domain edge assign inserted node indx
  ALLOCATE (rj_c_b(mesh%Nc_b,2))    ! For each boundary edge assign inserted node indx
  
  ALLOCATE (  pattern(mesh%Nm_d))   ! Domain Elements refinement pattern
  ALLOCATE (b_pattern(mesh%Nm_b))   ! Boundary Elements refinement pattern


  CALL Mark_d_Elements(mesh, np_refine,  el_refine,  nmr_p_to_add, &
                       dnmr_e_to_add, pattern)

  CALL Mark_b_Elements(mesh, bnp_refine, bel_refine, nmr_p_to_add, &
                       bnmr_e_to_add, b_pattern, bnmr_p_to_add)
  
  bnmr_p_to_add = 2 * bnmr_p_to_add


  !--Result--------------------------------------------------------------------------!
  !  nmr_p_to_add           -> Points to add to the grid
  !  bnmr_p_to_add          -> Points to add to the BOUNDARY grid 
  !  dnmr_e_to_add          -> Elements to add to the grid   
  !  bnmr_e_to_add          -> Elements to add to the BOUNDARY grid 
  !
  !  pattern(Ne_d)          -> (DIV) Refinement informations for domain elements
  !  b_pattern(Ne_b)        -> (DIV) Refinement informations for boundary elements
  !                            -see patterns.f90 for more details-
  !----------------------------------------------------------------------------------!
  
  ! Copying 'mesh' in 'old_mesh'
  old_mesh = Duplicate_Grid(mesh)  
  
  rNp_d = old_mesh % Nj_d + nmr_p_to_add
  rNp_b = old_mesh % Nj_b + bnmr_p_to_add  

  rNe_d = old_mesh % Nm_d + dnmr_e_to_add
  rNe_b = old_mesh % Nm_b + bnmr_e_to_add

  ALLOCATE (updt_cc(bnmr_p_to_add,2))

  DEALLOCATE (mesh%rr);              DEALLOCATE (mesh%jd_jb)
  ALLOCATE (mesh%rr(2,rNp_d));       ALLOCATE (mesh%jd_jb(rNp_b))

  DEALLOCATE (mesh%ele_type_d);      DEALLOCATE (mesh%ele_type_b)
  ALLOCATE (mesh%ele_type_d(rNe_d)); ALLOCATE (mesh%ele_type_b(rNe_b))

  DEALLOCATE (mesh%j_m_d);           DEALLOCATE (mesh%ma_m_d)
  ALLOCATE (mesh%j_m_d(rNe_d));      ALLOCATE (mesh%ma_m_d(rNe_d))
  
  DEALLOCATE (mesh%j_m_b);           DEALLOCATE (mesh%ma_m_b)
  ALLOCATE (mesh%j_m_b(rNe_b));      ALLOCATE (mesh%ma_m_b(rNe_b))

  DEALLOCATE (mesh%bound_p);         DEALLOCATE (mesh%bound_m)
  ALLOCATE (mesh%bound_p(rNp_b));    ALLOCATE (mesh%bound_m(rNe_b))




  ! Saving old grid informations and reallocating variables for new grid.
  mesh % rr      = 0.d0;  mesh %   rr(:, 1:old_mesh % Nj_d) = old_mesh % rr 
  mesh % jd_jb   = 0;     mesh %   jd_jb(1:old_mesh % Nj_b) = old_mesh % jd_jb
  mesh % bound_p = 0;     mesh % bound_p(1:old_mesh % Nj_b) = old_mesh % bound_p

  rj_c_d = 0

  ALLOCATE ( s_new(old_mesh % Nb)) !-For each boundary line spec new node curve coord
  ALLOCATE (jd_jsb(old_mesh % Nb)) !-For each boundary line spec domain index vs cc index
  ALLOCATE (  rj_b(old_mesh % Nb)) !-For each boundary line spec the number of nodes added

  rj_b = 0

  DO i = 1, old_mesh % Nb

   ALLOCATE (jd_jsb(i) %  index( bnmr_p_to_add/2 ))
   ALLOCATE ( s_new(i) % points( bnmr_p_to_add/2 ))

   jd_jsb(i) % index  = 0
   s_new(i) % points  = 2.d0

  ENDDO  
  

  CALL Add_Nodes(mesh, old_mesh, rj_c_d, rj_c_b, nds_nn, s_new, &
                 jd_jsb, updt_cc, rj_b, ref_history)


  ! Updating grid%curve_coordinates
  IF (ASSOCIATED(mesh % curve_coordinates))  DEALLOCATE(mesh % curve_coordinates)   
  ALLOCATE (mesh % curve_coordinates(old_mesh % Nj_b + bnmr_p_to_add))

  mesh%curve_coordinates(1:old_mesh%Nj_b) = old_mesh%curve_coordinates
  mesh%curve_coordinates(old_mesh%Nj_b+1:old_mesh%Nj_b+bnmr_p_to_add) = updt_cc(:,1)
  
  mesh%ele_type_d = 0
  mesh%ele_type_b = 0   
  mesh%bound_m    = 0


  ! Solution Interpolation
  ALLOCATE (add_ww(SIZE(sol_ww,1), nmr_p_to_add))   !rNp_d
  add_ww = 0.d0

  ! New Domain Elements                                                                    
  zone = 'D'                                                                               
  crop = 0                                                                                 
  CALL Build_Elements(    mesh%ele_type_d,     mesh%j_m_d,     mesh%ma_m_d, pattern,    & 
                      old_mesh%ele_type_d, old_mesh%j_m_d, old_mesh%ma_m_d, rj_c_d,     & 
                      np_refine, mesh%bound_m, old_mesh%bound_m, nds_nn, crop, mesh%rr, & 
                      mesh%c_m, zone, sol_ww, add_ww, old_mesh % Nj_d)

  ! New Boundary Elements                                                                  
  zone = 'B'                                                                               
  CALL Build_Elements(    mesh%ele_type_b,     mesh%j_m_b,     mesh%ma_m_b, b_pattern,  & 
                      old_mesh%ele_type_b, old_mesh%j_m_b, old_mesh%ma_m_b, rj_c_b,     & 
                      np_refine, mesh%bound_m, old_mesh%bound_m, nds_nn, crop, mesh%rr, & 
                      mesh%c_m, zone, sol_ww, add_ww, old_mesh % Nj_d)



  ! Updating nodes number as consequence of Q1 refinement: degeneration of two
  ! quadrilaters in triangles. This has effect only on the number of nodes
  IF (crop == 0) THEN 

    fNp_d = rNp_d

  ELSE

    fNp_d = rNp_d - crop
    
    IF(ALLOCATED(aux_rr)) DEALLOCATE(aux_rr)
    ALLOCATE(aux_rr(2,fNp_d))                
    
    aux_rr = mesh % rr(:,1:fNp_d)

    DEALLOCATE (mesh % rr) 
    ALLOCATE (mesh % rr(2,fNp_d))
    
    mesh % rr = aux_rr

    DEALLOCATE(aux_rr)


    ! Adjusting add_sol for solution interpolation
    ALLOCATE (COPY_add_ww(SIZE(add_ww,1), fNp_d - old_mesh % Nj_d))
    COPY_add_ww = add_ww(:, 1:fNp_d - old_mesh % Nj_d)

    DEALLOCATE (add_ww)
    ALLOCATE (add_ww(SIZE(sol_ww,1), fNp_d - old_mesh % Nj_d))
    add_ww = COPY_add_ww

    DEALLOCATE (COPY_add_ww)


  ENDIF

  mesh % Nj_d = fNp_d
  mesh % Nj_b = rNp_b

  mesh % Nm_d = rNe_d
  mesh % Nm_b = rNe_b

  mesh % ma_m_d = ma_m_connectivity(mesh % j_m_d, mesh % ele_type_d)
  mesh % ma_m_b = ma_m_connectivity(mesh % j_m_b, mesh % ele_type_b)


  IF (bou_tre_typ /= 0) THEN

    ! Computing nodes displacement towards real boundary
    ALLOCATE (displ(mesh % Nj_d))

    DO i = 1, mesh % Nj_d
       ALLOCATE (displ(i) % vec(2))
    ENDDO   

    displ = Bnd_Displ(mesh, s_new, jd_jsb, rj_b)

  ENDIF


  IF (ALLOCATED(rj_b))       DEALLOCATE(rj_b)
  IF (ALLOCATED(rj_c_d))     DEALLOCATE(rj_c_d)
  IF (ALLOCATED(rj_c_b))     DEALLOCATE(rj_c_b) 
  IF (ALLOCATED(updt_cc))    DEALLOCATE(updt_cc)
  IF (ALLOCATED(aux_rr))     DEALLOCATE(aux_rr)
  IF (ALLOCATED(s_new))      DEALLOCATE(s_new)
  IF (ALLOCATED(jd_jsb))     DEALLOCATE(jd_jsb)
  IF (ALLOCATED(pattern))    DEALLOCATE(pattern)
  IF (ALLOCATED(b_pattern))  DEALLOCATE(b_pattern)
 
 END SUBROUTINE Triangulate 
 
 

END MODULE triangulate_proc
