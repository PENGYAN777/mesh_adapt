MODULE  grid_utils

 !-------------------------------------------------------------------------
 !   Description: Contains subroutines for grid management, as the sub-
 !                routine for a full duplication of a grid, or routines
 !                for quality enforcement.
 !
 !
 !        Author: Marco Fossati
 !                Department of Aerospace Engineering
 !                Politecnico di Milano
 !                Via La Masa 34, 20156 Milano, ITALY
 !                e-mail: fossati@aero.polimi.it
 !
 !          Year:  2006, September
 !-------------------------------------------------------------------------

 USE  dynamic_vector
 USE  lib_geo
 USE  mesh_structure
 USE  metric_coefficients
 USE  node_pair_structure 
 USE  nodes
 USE  np_topology_gen 
 USE  structures

 USE  new_elements,   ONLY: angle
 
 ! Notation
 !
 !  - j  -->  nodes 
 !  - m  -->  elements
 !  - c  -->  edges
 ! 
 ! Nodes are ordered ANTICLOCKWISE around elements
 ! Edges are ordered     CLOCKWISE around TRIANGULAR elements
 ! Edges for quadrilateral elements are locally ordered to
 ! recover clockwise sense.

 !-------------------------------------------------------------------------------------
 IMPLICIT NONE

 TYPE  m_points
   REAL(KIND=8), DIMENSION(:,:), POINTER :: points
 END TYPE  m_points

 INTEGER, PARAMETER ::    TRIANGLE = 3, &
                       QUADRILATER = 4
 !-------------------------------------------------------------------------------------

 CONTAINS          

 
   FUNCTION  Compute_Tangent_Vector (node, grid)  RESULT (tangent_vector)
   !----------------------------------------------------------------------------------
   IMPLICIT NONE  
 
   TYPE(grid_type), INTENT(IN) :: grid
   INTEGER,         INTENT(IN) :: node

   REAL(KIND=8), DIMENSION(2) :: tangent_vector


   REAL(KIND=8) :: x1, x2, y1, y2

   INTEGER :: i, j, node1, node2 

   INTEGER, DIMENSION(:), ALLOCATABLE :: edges

   INTEGER, DIMENSION(2) :: boundary_edges
   INTEGER, DIMENSION(2) :: nodes_c1, nodes_c2
   !----------------------------------------------------------------------------------
   
   ALLOCATE (edges(SIZE(grid % c_j(node) % vec)))

   edges = grid % c_j(node) % vec
   
   j = 1
   
   DO i = 1, SIZE(edges)
    
     IF (grid % cb_cd(edges(i)) .NE. 0) THEN

       boundary_edges(j) = edges(i)
       j = j + 1

     ENDIF
      
   ENDDO
   
   nodes_c1 = grid % j_c(:, boundary_edges(1))
   nodes_c2 = grid % j_c(:, boundary_edges(2))
   
   IF (nodes_c1(1) .NE. node) THEN
     node1 = nodes_c1(1)
   ELSE
     node1 = nodes_c1(2)
   ENDIF       
      
   IF ( nodes_c2(1) .NE. node ) THEN
     node2 = nodes_c2(1)
   ELSE
     node2 = nodes_c2(2)
   ENDIF             
   
   x1 = grid % rr(1,node1);     y1 = grid % rr(2,node1);
   x2 = grid % rr(1,node2);     y2 = grid % rr(2,node2);
   
   tangent_vector(1) = (x2 - x1) / SQRT((y2 - y1)**2 + (x2 - x1)**2)
   tangent_vector(2) = (y2 - y1) / SQRT((y2 - y1)**2 + (x2 - x1)**2)

   DEALLOCATE (edges)

   END FUNCTION  Compute_Tangent_Vector




 
   SUBROUTINE  Invert_jd_jb (grid)  
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE  
  
   TYPE(grid_type), INTENT(INOUT) :: grid


   INTEGER :: dmn_idx, i
   !-------------------------------------------------------------------------------------

! Guardone
     IF (ALLOCATED(grid % jb_jd)) DEALLOCATE (grid % jb_jd)
! Guardone
   ALLOCATE (grid % jb_jd(2, grid % Nj_d)) 

   grid % jb_jd = 0
  
   DO i = 1, grid % Nj_b

     dmn_idx = grid % jd_jb(i)

     IF ( grid % jb_jd(1,dmn_idx) .EQ. 0 ) THEN

       grid % jb_jd(1,dmn_idx) = i

     ELSE

       grid % jb_jd(2,dmn_idx) = i 

     ENDIF   

   ENDDO 

   END SUBROUTINE  Invert_jd_jb
 
   
   

 
   SUBROUTINE  Invert_cd_cb (grid) 
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
 
   TYPE(grid_type), INTENT(INOUT) :: grid

  
   INTEGER :: np_d_idx, i
   !-------------------------------------------------------------------------------------

! Guardone
     IF (ALLOCATED(grid % cb_cd)) DEALLOCATE (grid % cb_cd)
! Guardone
   ALLOCATE (grid % cb_cd(grid % Nc_d))
  
   grid % cb_cd = 0 
  
   DO i = 1, grid % Nc_b
  
     np_d_idx = grid % cd_cb(i)
     grid % cb_cd(np_d_idx) = i
     
   ENDDO

   END SUBROUTINE  Invert_cd_cb




   
   SUBROUTINE  Edge_Structure_Gen (grid, action)

   ! Build node_pairs structure on the basis
   ! of dynamic vectors j_m_d and j_m_b.
 
   ! action = 'P'  partial edge structure generation  = generates only FV connectivity
   ! action = 'C'  complete edge structure generation = generates FE and FV connectivity 
   !---------------------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type),  INTENT(INOUT) :: grid
   CHARACTER(LEN=1), INTENT(IN)    :: action


   INTEGER                                  :: i, Ncd

   INTEGER,     DIMENSION(:),   ALLOCATABLE :: fv_fe
   INTEGER,     DIMENSION(:,:), ALLOCATABLE :: j_c, jc_fv, jcd, jcb

   TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE :: cmd, cmb, j_c_DIV, c_j, m_c
   !---------------------------------------------------------------------------------------
  
   IF (ALLOCATED(cmd)) DEALLOCATE(cmd, jcd, cmb, jcb, cd_cb, jc_fv)
                                                
   ! DOMAIN edges:
   
   ! c_m_d (FE)
   ALLOCATE (cmd(grid % Nm_d))
   cmd = c_m_connectivity(grid % j_m_d)
   
   ! Number of node pairs (FE)
   Ncd = size_DIV(cmd,3)

   ! j_c_d (FE)
   ALLOCATE (jcd(4,Ncd))
   jcd = j_c_connectivity(grid % j_m_d, cmd)


   ! JCD_FEM, CJD_FEM
   IF (action .EQ. 'C') THEN

! Guardone
     IF (ALLOCATED(grid % jcd_fem)) DEALLOCATE (grid % jcd_fem)
! Guardone
     ALLOCATE (grid % jcd_fem(4, Ncd))

     grid % jcd_fem = 0
     grid % jcd_fem = jcd

   ENDIF 
   
   ! FV-FE connectivity
   IF (ALLOCATED(fv_fe)) DEALLOCATE(fv_fe)
   ALLOCATE (fv_fe(Ncd))

   fv_fe = fv_fe_np_connectivity(jcd, grid % j_m_d, cmd, grid % ele_type_d, grid % Nc_d)
   
   
   ! J_C
   ALLOCATE (jc_fv(4,grid % Nc_d)) 
   jc_fv = fv_node_pair(jcd, fv_fe, grid % Nc_d)
   
   IF (action .EQ. 'C') THEN
          
! Guardone
     IF (ALLOCATED(grid % j_c)) DEALLOCATE (grid % j_c)
! Guardone
     ALLOCATE (grid % j_c(2, grid % Nc_d))
     grid % j_c = jc_fv(1:2,:)

   ENDIF 

   ! C_J
   ALLOCATE (j_c(4, grid % Nc_d))  
   j_c(1:2,:) = grid % j_c
   j_c(3:4,:) = 0  ! dummy

   ALLOCATE (j_c_DIV(SIZE(j_c,2)))
   j_c_DIV = convert_matrix_to_DIV(j_c)

   ALLOCATE (c_j(size_DIV(j_c_DIV, 3)))
   c_j = invert_DIV(j_c_DIV)

   DEALLOCATE (j_c_DIV, j_c)    
      
! Guardone
     IF (ASSOCIATED(grid % c_j)) DEALLOCATE (grid % c_j)
! Guardone
   ALLOCATE (grid % c_j(grid % Nj_d))
    
   DO i = 1, grid % Nj_d       
     ALLOCATE (grid % c_j(i) % vec(SIZE(c_j(i) % vec)))
     grid % c_j(i) % vec = c_j(i) % vec
   ENDDO 
   
   DEALLOCATE (c_j)   


   ! C_M 
! Guardone
     IF (ASSOCIATED(grid % c_m)) DEALLOCATE (grid % c_m)
! Guardone
   ALLOCATE (grid % c_m(grid % Nm_d))
   grid % c_m = fv_c_m_connectivity(cmd, fv_fe)


   ! M_C 
   ALLOCATE (m_c(size_DIV(grid % c_m, 3)))
   m_c = invert_DIV(grid % c_m)
   
! Guardone
     IF (ASSOCIATED(grid % m_c)) DEALLOCATE (grid % m_c)
! Guardone
   ALLOCATE (grid % m_c(grid % Nc_d))
   
   DO i = 1, grid % Nc_d  

     ALLOCATE(grid % m_c(i) % vec(SIZE(m_c(i) % vec)))
     grid % m_c(i) % vec = m_c(i) % vec

   ENDDO



  ! BOUNDARY edges:

   ALLOCATE (cmb(grid % Nm_b))
   cmb = c_m_connectivity(grid % j_m_b)

   grid % Nc_b = size_DIV(cmb, 3)

   ALLOCATE (jcb(4, grid % Nc_b))
   jcb = j_c_connectivity(grid % j_m_b, cmb)

   ! J_C_B 
   IF (action .EQ. 'C') THEN

! Guardone
     IF (ALLOCATED(grid % j_c_b)) DEALLOCATE (grid % j_c_b)
! Guardone
     ALLOCATE (grid % j_c_b(2,grid % Nc_b))
     grid % j_c_b = jcb(1:2,:)

   ENDIF
   
   ! JCB_FEM
   IF (action .EQ. 'C') THEN

! Guardone
     IF (ALLOCATED(grid % jcb_fem)) DEALLOCATE (grid % jcb_fem)
! Guardone
     ALLOCATE(grid % jcb_fem(4, grid % Nc_b))
     grid % jcb_fem = jcb(:,:)   

   ENDIF 
   
   
   ! CD_CB 
   IF (action .EQ. 'C') THEN

     IF (ALLOCATED(cd_cb)) DEALLOCATE(cd_cb)
     ALLOCATE (cd_cb(grid % Nc_b))

     cd_cb = cd_cb_connectivity(grid % jcd_fem, grid % jd_jb, jcb)

! Guardone
     IF (ALLOCATED(grid % cd_cb)) DEALLOCATE (grid % cd_cb)
! Guardone
     ALLOCATE (grid % cd_cb(grid % Nc_b))

     DO i = 1, grid % Nc_b
       grid % cd_cb(i) = fv_fe(cd_cb(i))
     ENDDO 

   ENDIF
   
   IF (action .EQ. 'P') THEN

     DO i = 1, SIZE(grid % cd_cb)
       grid % cd_cb(i) = fv_fe(grid % cd_cb(i))
     ENDDO  

   ENDIF 
   
   ! CB_CD 
   CALL Invert_cd_cb(grid)
   
   
   ! BOUND_C
   IF (action .EQ. 'C') THEN

! Guardone
     IF (ALLOCATED(grid % bound_c)) DEALLOCATE (grid % bound_c)
! Guardone
     ALLOCATE (grid % bound_c(grid % Nc_b))
     grid % bound_c = bound_connectivity(cmb, grid % bound_m, grid % Nc_b)

   ENDIF        

   ! Deallocating
   IF (ALLOCATED(fv_fe))  DEALLOCATE (fv_fe)
   IF (ALLOCATED(jcd))    DEALLOCATE (jcd)
   IF (ALLOCATED(jcb))    DEALLOCATE (jcb)
   IF (ALLOCATED(jc_fv))  DEALLOCATE (jc_fv)
   IF (ALLOCATED(m_c))    DEALLOCATE (m_c)
   IF (ALLOCATED(cmd))    DEALLOCATE (cmd)
   IF (ALLOCATED(cmb))    DEALLOCATE (cmb)

   END SUBROUTINE  Edge_structure_gen 




   
   FUNCTION  Extract_np_index(cj, n1, n2)  RESULT (np_index)
   !--------------------------------------------------------------------
   IMPLICIT NONE
 
   TYPE(D_I_V), DIMENSION(:), INTENT(IN) :: cj
   INTEGER,                   INTENT(IN) :: n1, n2

   INTEGER :: np_index, np_found
 

   INTEGER, DIMENSION(:), ALLOCATABLE :: cj_n1, cj_n2
 
   INTEGER :: j, k, dim1, dim2
   !--------------------------------------------------------------------

   np_index = 0 
   np_found = 1 

   dim1 = SIZE(cj(n1) % vec)
   dim2 = SIZE(cj(n2) % vec)

   ALLOCATE (cj_n1(dim1))
   ALLOCATE (cj_n2(dim2))

   cj_n1 = cj(n1) % vec
   cj_n2 = cj(n2) % vec
!WRITE(*,*) 'cj_n1=', cj_n1
!WRITE(*,*) 'cj_n2=', cj_n2
   DO j = 1, dim1

     DO k = 1, dim2

       IF (cj_n1(j) == cj_n2(k)) THEN

         np_index = cj_n1(j)
         np_found = 1
 
         DEALLOCATE (cj_n1, cj_n2)
       
         RETURN

       ENDIF

     ENDDO

   ENDDO


   IF (np_index == 0) THEN

     PRINT*, ''
     PRINT*, 'ERROR. EXTRACT_NP_INDEX:' 
     PRINT*, 'Edge index not found.'
     PRINT*, 'np_found = ', np_found
     PRINT*, ''

     STOP

   ENDIF

   END FUNCTION  Extract_np_index   




 
   FUNCTION  Duplicate_Grid (grid_in)   RESULT (grid_out)
   !------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type), INTENT(IN) :: grid_in
   
   TYPE(grid_type):: grid_out
 

   INTEGER :: i, ns, nd
   !------------------------------------------------------------------- 

   grid_out%adaptation_level = grid_in%adaptation_level
   grid_out%k_d              = grid_in%k_d
   grid_out%Nb               = grid_in%Nb
   
   ! Nodes
   grid_out%Nj_d = grid_in%Nj_d;  grid_out%Nj_b = grid_in%Nj_b

                 
   ALLOCATE( grid_out%rr(grid_in%k_d,grid_in%Nj_d)  )
   ALLOCATE( grid_out%jd_jb(grid_in%Nj_b)   )
   ALLOCATE( grid_out%jb_jd(2,grid_in%Nj_d) )
   ALLOCATE( grid_out%bound_p(grid_in%Nj_b) )
   
   grid_out%rr      = grid_in%rr
   grid_out%jd_jb   = grid_in%jd_jb
   grid_out%jb_jd   = grid_in%jb_jd
   grid_out%bound_p = grid_in%bound_p
   
   ! Domain grid
   grid_out%Nm_d = grid_in%Nm_d  
   
   ALLOCATE( grid_out%m_j_d(grid_in%Nj_d)      )
   ALLOCATE( grid_out%j_m_d(grid_in%Nm_d)      )
   ALLOCATE( grid_out%ma_m_d(grid_in%Nm_d)     )

   DO i = 1, grid_out%Nj_d        
     ALLOCATE( grid_out%m_j_d(i)%vec(SIZE(grid_in%m_j_d(i)%vec)) )      
     grid_out%m_j_d(i)%vec = grid_in%m_j_d(i)%vec   
   ENDDO    
   
   DO i = 1, grid_out%Nm_d    
     
     ALLOCATE( grid_out%j_m_d(i)%vec(SIZE(grid_in%j_m_d(i)%vec)) )
     grid_out%j_m_d(i)%vec = grid_in%j_m_d(i)%vec
     
     ALLOCATE( grid_out%ma_m_d(i)%vec(SIZE(grid_in%ma_m_d(i)%vec)) )
     grid_out%ma_m_d(i)%vec = grid_in%ma_m_d(i)%vec
     
   ENDDO
   
   ALLOCATE( grid_out%ele_type_d(grid_in%Nm_d) )    
   grid_out%ele_type_d = grid_in%ele_type_d


   ! Boundary grid
   grid_out%Nm_b = grid_in%Nm_b  
   
   ALLOCATE( grid_out%m_j_b(grid_in%Nj_b)      )
   ALLOCATE( grid_out%j_m_b(grid_in%Nm_b)      )
   ALLOCATE( grid_out%ma_m_b(grid_in%Nm_b)     )

   DO i = 1, grid_out%Nj_b        
     ALLOCATE( grid_out%m_j_b(i)%vec(SIZE(grid_in%m_j_b(i)%vec)) )      
     grid_out%m_j_b(i)%vec = grid_in%m_j_b(i)%vec   
   ENDDO    
   
   DO i = 1, grid_out%Nm_b    
     
     ALLOCATE( grid_out%j_m_b(i)%vec(SIZE(grid_in%j_m_b(i)%vec)) )
     grid_out%j_m_b(i)%vec = grid_in%j_m_b(i)%vec
     
     ALLOCATE( grid_out%ma_m_b(i)%vec(SIZE(grid_in%ma_m_b(i)%vec)) )
     grid_out%ma_m_b(i)%vec = grid_in%ma_m_b(i)%vec
     
   ENDDO

   ALLOCATE( grid_out%ele_type_b(grid_in%Nm_b) )    
   grid_out%ele_type_b = grid_in%ele_type_b
   
   ALLOCATE( grid_out%bound_m(grid_in%Nm_b) )
   grid_out%bound_m = grid_in%bound_m
   
   
   ! Edges
   grid_out%Nc_d = grid_in%Nc_d
   grid_out%Nc_b = grid_in%Nc_b
   
   ALLOCATE( grid_out%j_c(2,grid_in%Nc_d) ) 
   grid_out%j_c = grid_in%j_c
   
   ALLOCATE( grid_out%jcd_fem(4,SIZE(grid_in%jcd_fem,2)) )    
   grid_out%jcd_fem = grid_in%jcd_fem    
   
   ALLOCATE( grid_out%j_c_b(2,grid_in%Nc_b) )
   grid_out%j_c_b = grid_in%j_c_b    
   
   ALLOCATE( grid_out%jcb_fem(4,SIZE(grid_in%jcb_fem,2)) )    
   grid_out%jcb_fem = grid_in%jcb_fem        
   
   ALLOCATE( grid_out%c_j(grid_in%Nj_d) )        
   DO i = 1, grid_in%Nj_d
     ALLOCATE( grid_out%c_j(i)%vec(SIZE(grid_in%c_j(i)%vec)) )
     grid_out%c_j(i)%vec = grid_in%c_j(i)%vec
   ENDDO    
   
   ALLOCATE( grid_out%m_c(grid_in%Nc_d) )    
   DO i = 1, grid_in%Nc_d  
     ALLOCATE( grid_out%m_c(i)%vec(SIZE(grid_in%m_c(i)%vec)) )
     grid_out%m_c(i)%vec = grid_in%m_c(i)%vec
   ENDDO
   
   ALLOCATE( grid_out%c_m(grid_in%Nm_d) )     
   DO i = 1, grid_in%Nm_d
     ALLOCATE( grid_out%c_m(i)%vec(SIZE(grid_in%c_m(i)%vec)) )
     grid_out%c_m(i)%vec = grid_in%c_m(i)%vec
   ENDDO
           
   ALLOCATE( grid_out%cd_cb(grid_in%Nc_b) )    
   grid_out%cd_cb = grid_in%cd_cb 

   ALLOCATE( grid_out%cb_cd(grid_in%Nc_d) )    
   grid_out%cb_cd = grid_in%cb_cd     
   
   ALLOCATE( grid_out%bound_c(grid_in%Nc_b) )    
   grid_out%bound_c = grid_in%bound_c
   

WRITE(*,*) SIZE(grid_in%cell), grid_in%Nj_d
   ALLOCATE( grid_out%cell(grid_in%Nj_d) )     
!   ALLOCATE( grid_out%cell(SIZE(grid_in%cell)) )     
   grid_out%cell = grid_in%cell 
   
   ALLOCATE( grid_out%curve_coordinates(grid_in%Nj_b) )
   grid_out%curve_coordinates = grid_in%curve_coordinates
   
   ALLOCATE( grid_out%line(grid_in%Nb) )
    
   DO i = 1, grid_in % Nb
     
     grid_out % line(i) % nd = grid_in % line(i) % nd
     grid_out % line(i) % ns = grid_in % line(i) % ns
     
     nd = grid_in % line(i) % nd
     ns = grid_in % line(i) % ns
     
     ALLOCATE (grid_out % line(i) % s(ns),    grid_out % line(i) % h(ns))
     ALLOCATE (grid_out % line(i) % x(nd,ns), grid_out % line(i) % xs(nd,ns))
     
     grid_out % line(i) % s  = grid_in % line(i) % s
     grid_out % line(i) % h  = grid_in % line(i) % h 
     grid_out % line(i) % x  = grid_in % line(i) % x
     grid_out % line(i) % xs = grid_in % line(i) % xs
     
   ENDDO
    
   END FUNCTION  Duplicate_grid



 

 
  
   SUBROUTINE  Mark_for_Quality (grid, angle_treshold, maxN_loop, refined_edges_flag)
   !---------------------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type),       INTENT(IN)    :: grid
   REAL(KIND=8),          INTENT(IN)    :: angle_treshold
   INTEGER,               INTENT(IN)    :: maxN_loop
   LOGICAL, DIMENSION(:), INTENT(INOUT) :: refined_edges_flag

   TYPE(m_points), DIMENSION(:), POINTER :: v__rr

   TYPE(D_I_V) :: c_m__loc

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: nodes_m

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: length, alpha

   REAL(KIND=8), DIMENSION(2) :: rr_1, rr_2, rr_3, rr_4

   REAL(KIND=8) :: min_alpha, diagonal

   INTEGER, DIMENSION(:,:), ALLOCATABLE :: edges_m

   CHARACTER(LEN=3) :: ref_type

   INTEGER :: c, j, j_, m, m_, Nref_np, np_1, np_2, np_3, np_4, &
              added_edges, quality_loop

   LOGICAL :: start

   REAL(KIND=8), PARAMETER :: PI = 3.141592654
   !---------------------------------------------------------------------------------------

   start = .TRUE.
   added_edges = 0
   quality_loop = 0

   PRINT*, '   Quality enforcement...'

   DO WHILE ((start  .OR.  added_edges > 0)  .AND.  quality_loop < maxN_loop)

     start = .FALSE.
     added_edges = 0

     DO m = 1, grid % Nm_d

       ! Element nodes coordinates
       ALLOCATE (nodes_m(SIZE(grid % j_m_d(m) % vec),2))

       DO j_ = 1, SIZE(nodes_m,1)

         j = grid % j_m_d(m) % vec(j_)

         nodes_m(j_,1) = grid % rr(1,j)
         nodes_m(j_,2) = grid % rr(2,j)

       ENDDO


       ! Element edges and relative refinement intention
       Nref_np = 0

       ALLOCATE (c_m__loc % vec(SIZE(grid % c_m(m) % vec)))

       ! Force clockwise order of edges in quadrilateral elements
       IF (grid % ele_type_d(m) == 2) THEN
 
         c_m__loc % vec = grid % c_m(m) % vec 

       ELSEIF (grid % ele_type_d(m) == 3) THEN

         c_m__loc % vec(1:2) = grid % c_m(m) % vec(1:2)

         c_m__loc % vec(3) = grid % c_m(m) % vec(4)
         c_m__loc % vec(4) = grid % c_m(m) % vec(3)

       ENDIF

       ! write(*,*) 'Error 672:',m, SIZE(grid % c_m(m) % vec), c_m__loc % vec

       ALLOCATE (edges_m(SIZE(grid % c_m(m) % vec),2))      
       edges_m(:,1) = c_m__loc % vec

       DO c = 1, SIZE(edges_m, 1)
 
         IF (refined_edges_flag(edges_m(c,1))) THEN
           edges_m(c,2) = 1
         ELSE
           edges_m(c,2) = 0
         ENDIF

       ENDDO

       Nref_np = SUM(edges_m(:,2))


       ! Element refinement intention
       ref_type = '000'

       ! write(*,*) 'Error 695:',m, ele_type_d(m)

       IF (grid % ele_type_d(m) == 2  .AND.  Nref_np == 2)  ref_type = 'T2 '
       IF (grid % ele_type_d(m) == 3  .AND.  Nref_np == 1)  ref_type = 'Q1 '
       IF (grid % ele_type_d(m) == 3  .AND.  Nref_np == 3)  ref_type = 'Q3 '

       IF ((grid % ele_type_d(m) == 3  .AND.  Nref_np == 2)  .AND.  &

           (ALL(edges_m(:,2) == (/1,1,0,0/))  .OR. &  
            ALL(edges_m(:,2) == (/0,1,1,0/))  .OR. &
            ALL(edges_m(:,2) == (/0,0,1,1/))  .OR. &
            ALL(edges_m(:,2) == (/1,0,0,1/))))  ref_type = 'Q2 '

       IF (ref_type /= '000') THEN

         min_alpha = HUGE(min_alpha)

         CALL Virtual_Elements(ref_type, nodes_m, edges_m, v__rr)


         DO m_ = 1, SIZE(v__rr)

           ALLOCATE (length(SIZE(v__rr(m_) % points,1)),  &
                      alpha(SIZE(v__rr(m_) % points,1)))

           SELECT CASE (SIZE(length))

             CASE (TRIANGLE)

               rr_1 = v__rr(m_) % points(1,:)
               rr_2 = v__rr(m_) % points(2,:)
               rr_3 = v__rr(m_) % points(3,:)

               length(1) = SQRT((rr_2(2) - rr_1(2))**2 + (rr_2(1) - rr_1(1))**2)
               length(2) = SQRT((rr_3(2) - rr_1(2))**2 + (rr_3(1) - rr_1(1))**2)
               length(3) = SQRT((rr_2(2) - rr_3(2))**2 + (rr_2(1) - rr_3(1))**2)

               alpha(1) = ACOS((length(1)**2 + length(2)**2 - length(3)**2) / (2*length(1)*length(2)))
               alpha(2) = ACOS((length(1)**2 + length(3)**2 - length(2)**2) / (2*length(1)*length(3)))
               alpha(3) = ACOS((length(3)**2 + length(2)**2 - length(1)**2) / (2*length(3)*length(2)))

               min_alpha = MIN(min_alpha, MINVAL(alpha,1))


             CASE (QUADRILATER)

               rr_1 = v__rr(m_) % points(1,:)
               rr_2 = v__rr(m_) % points(2,:)
               rr_3 = v__rr(m_) % points(3,:)
               rr_4 = v__rr(m_) % points(4,:)

               length(1) = SQRT((rr_2(2) - rr_1(2))**2 + (rr_2(1) - rr_1(1))**2)
               length(2) = SQRT((rr_4(2) - rr_1(2))**2 + (rr_4(1) - rr_1(1))**2)
               length(3) = SQRT((rr_4(2) - rr_3(2))**2 + (rr_4(1) - rr_3(1))**2)
               length(4) = SQRT((rr_2(2) - rr_3(2))**2 + (rr_2(1) - rr_3(1))**2)

               diagonal = SQRT((rr_2(2) - rr_4(2))**2 + (rr_2(1) - rr_4(1))**2)

               alpha(1) = ACOS((length(1)**2 + length(2)**2 - diagonal**2) / (2*length(1)*length(2)))
               alpha(3) = ACOS((length(3)**2 + length(4)**2 - diagonal**2) / (2*length(3)*length(4)))

               diagonal = SQRT((rr_3(2) - rr_1(2))**2 + (rr_3(1) - rr_1(1))**2)

               alpha(2) = ACOS((length(1)**2 + length(4)**2 - diagonal**2) / (2*length(1)*length(4)))
               alpha(4) = ACOS((length(3)**2 + length(2)**2 - diagonal**2) / (2*length(3)*length(2)))

               min_alpha = MIN(min_alpha, MINVAL(alpha,1))

           END SELECT

           DEALLOCATE (length, alpha)

         ENDDO


         IF (min_alpha < angle_treshold*(PI/180.d0)) THEN

           SELECT CASE (ref_type)

             CASE ('T2 ')

               DO c = 1, SIZE(edges_m,1)

                 IF (edges_m(c,2) == 0) THEN

                     np_1 = edges_m(c,1) 
                     refined_edges_flag(np_1) = .TRUE.

                     EXIT
 
                 ENDIF

               ENDDO 

               added_edges = added_edges + 1


             CASE ('Q1 ')

               DO c = 1, SIZE(edges_m,1)

                 IF (edges_m(c,2) == 1  .AND. (c == 1  .OR. c == 3)) THEN

                     np_1 = edges_m(1,1)
                     np_3 = edges_m(3,1)
 
                     refined_edges_flag(np_1) = .TRUE.
                     refined_edges_flag(np_3) = .TRUE.

                     EXIT

                 ELSEIF (edges_m(c,2) == 1  .AND. (c == 2  .OR. c == 4)) THEN

                     np_2 = edges_m(2,1)
                     np_4 = edges_m(4,1)

                     refined_edges_flag(np_2) = .TRUE.
                     refined_edges_flag(np_4) = .TRUE.

                     EXIT

                 ENDIF

               ENDDO 

               added_edges = added_edges + 1


             CASE ('Q2 ')

               DO c = 1, SIZE(edges_m,1)

                 IF (edges_m(c,2) == 0) THEN

                     np_1 = edges_m(c,1) 
                     refined_edges_flag(np_1) = .TRUE.

                 ENDIF

               ENDDO 

               added_edges = added_edges + 2


             CASE ('Q3 ')

               DO c = 1, SIZE(edges_m,1)

                 IF (edges_m(c,2) == 0) THEN

                     np_1 = edges_m(c,1) 
                     refined_edges_flag(np_1) = .TRUE.

                 ENDIF

               ENDDO 

               added_edges = added_edges + 1

           END SELECT

         ENDIF


         DO m_ = 1, SIZE(v__rr)
           DEALLOCATE (v__rr(m_) % points)
         ENDDO

         DEALLOCATE (v__rr)

       ENDIF     

       DEALLOCATE (edges_m, nodes_m, c_m__loc % vec)

     ENDDO

     quality_loop = quality_loop + 1

   ENDDO

!   PRINT*, ' - number of quality loops performed:', quality_loop



   CONTAINS


     SUBROUTINE Virtual_Elements(m_ref_type, nodes_m, edges_m, v__rr)
     !-------------------------------------------------------------------------------
     IMPLICIT NONE

     CHARACTER(LEN=3), INTENT(IN) :: m_ref_type
     REAL(KIND=8), DIMENSION(:,:) :: nodes_m
     INTEGER,      DIMENSION(:,:) :: edges_m

     TYPE(m_points), DIMENSION(:), POINTER :: v__rr


     REAL(KIND=8), DIMENSION(2) :: rnd1, rnd2, rnd3, rnd4
     REAL(KIND=8), DIMENSION(2) :: aux_R

     CHARACTER(LEN=1) :: config

     INTEGER, DIMENSION(2) :: aux_I

     INTEGER :: l
     !-------------------------------------------------------------------------------

     SELECT CASE (m_ref_type)

       CASE('T2 ')

1         l = 1
         
          IF (edges_m(l,2) == 0) THEN   
          
            ! Node pairs
            aux_I          = edges_m(l,  :)
            edges_m(l,  :) = edges_m(l+1,:)
            edges_m(l+1,:) = edges_m(l+2,:)
            edges_m(l+2,:) = aux_I
            
            ! Nodes
            aux_R          = nodes_m(l+2,:)
            nodes_m(l+2,:) = nodes_m(l+1,:)
            nodes_m(l+1,:) = nodes_m(l,  :)
            nodes_m(l,  :) = aux_R

            GOTO 1

          ENDIF

          l = l + 1

          IF (edges_m(l,2) == 1) THEN

            ! Node pairs
            aux_I          = edges_m(l-1,:)
            edges_m(l-1,:) = edges_m(l,  :)
            edges_m(l,  :) = edges_m(l+1,:)
            edges_m(l+1,:) = aux_I

            ! Nodes
            aux_R          = nodes_m(l+1,:)
            nodes_m(l+1,:) = nodes_m(l,  :)
            nodes_m(l,  :) = nodes_m(l-1,:)
            nodes_m(l-1,:) = aux_R

            GOTO 1

           ENDIF            


           rnd1 = (nodes_m(1,:) + nodes_m(2,:)) / 2.d0
           rnd2 = (nodes_m(2,:) + nodes_m(3,:)) / 2.d0

           ! Coordinates of nodes of virtual new triangles
           ALLOCATE (v__rr(3))

           ALLOCATE (v__rr(1) % points(3,2))
           v__rr(1) % points(1,:) = nodes_m(1,:)
           v__rr(1) % points(2,:) = rnd1
           v__rr(1) % points(3,:) = nodes_m(3,:)

           ALLOCATE (v__rr(2) % points(3,2))
           v__rr(2) % points(1,:) = rnd1
           v__rr(2) % points(2,:) = nodes_m(2,:)
           v__rr(2) % points(3,:) = rnd2

           ALLOCATE (v__rr(3) % points(3,2))
           v__rr(3) % points(1,:) = rnd1
           v__rr(3) % points(2,:) = rnd2
           v__rr(3) % points(3,:) = nodes_m(3,:)



        CASE('Q1 ')

2         l = 1  

          IF (edges_m(l,2) == 0) THEN   

            ! Node pairs
            aux_I          = edges_m(l,  :)
            edges_m(l,  :) = edges_m(l+1,:)
            edges_m(l+1,:) = edges_m(l+2,:)
            edges_m(l+2,:) = edges_m(l+3,:)
            edges_m(l+3,:) = aux_I
            
            ! Nodes
            aux_R          = nodes_m(l+3,:)
            nodes_m(l+3,:) = nodes_m(l+2,:)
            nodes_m(l+2,:) = nodes_m(l+1,:)
            nodes_m(l+1,:) = nodes_m(l,  :)
            nodes_m(l,  :) = aux_R

            GOTO 2

          ENDIF 


          rnd1 = (nodes_m(1,:) + nodes_m(2,:)) / 2.d0
          rnd2 = Q1_node_coords(nodes_m, rnd1, config)

          IF (config == 'T') THEN

            ! Coordinates of nodes of virtual new triangles
            ALLOCATE (v__rr(3))

            ALLOCATE (v__rr(1) % points(3,2))
            v__rr(1) % points(1,:) = nodes_m(1,:)
            v__rr(1) % points(2,:) = rnd1
            v__rr(1) % points(3,:) = nodes_m(4,:)

            ALLOCATE (v__rr(2) % points(3,2))
            v__rr(2) % points(1,:) = rnd1
            v__rr(2) % points(2,:) = nodes_m(2,:)
            v__rr(2) % points(3,:) = nodes_m(3,:)

            ALLOCATE (v__rr(3) % points(3,2))
            v__rr(3) % points(1,:) = rnd1
            v__rr(3) % points(2,:) = nodes_m(3,:)
            v__rr(3) % points(3,:) = nodes_m(4,:)

          ELSE              
                 
            ! Coordinates of nodes of virtual new elements
            ALLOCATE (v__rr(3))

            ALLOCATE (v__rr(1) % points(4,2))
            v__rr(1) % points(1,:) = nodes_m(1,:)
            v__rr(1) % points(2,:) = rnd1
            v__rr(1) % points(3,:) = rnd2
            v__rr(1) % points(4,:) = nodes_m(4,:)

            ALLOCATE (v__rr(2) % points(4,2))
            v__rr(2) % points(1,:) = rnd1
            v__rr(2) % points(2,:) = nodes_m(2,:)
            v__rr(2) % points(3,:) = nodes_m(3,:)
            v__rr(2) % points(4,:) = rnd2

            ALLOCATE (v__rr(3) % points(3,2))
            v__rr(3) % points(1,:) = rnd2
            v__rr(3) % points(2,:) = nodes_m(3,:)
            v__rr(3) % points(3,:) = nodes_m(4,:)

          ENDIF


        CASE('Q2 ')

3         l = 1  

          IF (edges_m(l,2) == 0) THEN   

            ! Node pairs
            aux_I          = edges_m(l,  :)
            edges_m(l,  :) = edges_m(l+1,:)
            edges_m(l+1,:) = edges_m(l+2,:)
            edges_m(l+2,:) = edges_m(l+3,:)
            edges_m(l+3,:) = aux_I
            
           ! Nodes
            aux_R          = nodes_m(l+3,:)
            nodes_m(l+3,:) = nodes_m(l+2,:)
            nodes_m(l+2,:) = nodes_m(l+1,:)
            nodes_m(l+1,:) = nodes_m(l,  :)
            nodes_m(l,  :) = aux_R

            GOTO 3

          ENDIF

          l = l + 1

          IF (edges_m(l,2) == 1) THEN   
                    
            ! Node pairs
            aux_I          = edges_m(l-1,:)
            edges_m(l-1,:) = edges_m(l,  :)
            edges_m(l,  :) = edges_m(l+1,:)
            edges_m(l+1,:) = edges_m(l+2,:)
            edges_m(l+2,:) = aux_I
            
            ! Nodes
            aux_R          = nodes_m(l+2,:)
            nodes_m(l+2,:) = nodes_m(l+1,:)
            nodes_m(l+1,:) = nodes_m(l,  :)
            nodes_m(l,  :) = nodes_m(l-1,:)
            nodes_m(l-1,:) = aux_R
            
            GOTO 3      
           
          ENDIF     


          rnd1 = (nodes_m(1,:) + nodes_m(2,:)) / 2.d0
          rnd2 = (nodes_m(2,:) + nodes_m(3,:)) / 2.d0
          rnd3 = SUM(nodes_m,1) / 4.d0

          ! Coordinates of nodes of virtual new elements
          ALLOCATE (v__rr(3))

          ALLOCATE (v__rr(1) % points(4,2))
          v__rr(1) % points(1,:) = nodes_m(1,:)
          v__rr(1) % points(2,:) = rnd1
          v__rr(1) % points(3,:) = rnd3
          v__rr(1) % points(4,:) = nodes_m(4,:)

          ALLOCATE (v__rr(2) % points(4,2))
          v__rr(2) % points(1,:) = rnd1
          v__rr(2) % points(2,:) = nodes_m(2,:)
          v__rr(2) % points(3,:) = rnd2
          v__rr(2) % points(4,:) = rnd3

          ALLOCATE (v__rr(3) % points(4,2))
          v__rr(3) % points(1,:) = rnd2
          v__rr(3) % points(2,:) = nodes_m(3,:)
          v__rr(3) % points(3,:) = nodes_m(4,:)
          v__rr(3) % points(4,:) = rnd3



        CASE('Q3 ')
        
4         l = 1

          IF (edges_m(l,2) == 0) THEN   

            ! Node pairs
            aux_I          = edges_m(l,  :)
            edges_m(l,  :) = edges_m(l+1,:)
            edges_m(l+1,:) = edges_m(l+2,:)
            edges_m(l+2,:) = edges_m(l+3,:)
            edges_m(l+3,:) = aux_I
            
            ! Nodes
            aux_R          = nodes_m(l+3,:)
            nodes_m(l+3,:) = nodes_m(l+2,:)
            nodes_m(l+2,:) = nodes_m(l+1,:)
            nodes_m(l+1,:) = nodes_m(l,  :)
            nodes_m(l,  :) = aux_R

            GOTO 4

          ENDIF

          l = l + 1

          IF (edges_m(l,2) == 1) THEN
                    
            ! Node pairs
            aux_I          = edges_m(l-1,:)
            edges_m(l-1,:) = edges_m(l,  :)
            edges_m(l,  :) = edges_m(l+1,:)
            edges_m(l+1,:) = edges_m(l+2,:)
            edges_m(l+2,:) = aux_I
            
            ! Nodes
            aux_R          = nodes_m(l+2,:)
            nodes_m(l+2,:) = nodes_m(l+1,:)
            nodes_m(l+1,:) = nodes_m(l,  :)
            nodes_m(l,  :) = nodes_m(l-1,:)
            nodes_m(l-1,:) = aux_R

            GOTO 4      
           
           ENDIF            


          rnd1 = (nodes_m(1,:) + nodes_m(2,:)) / 2.d0
          rnd2 = (nodes_m(2,:) + nodes_m(3,:)) / 2.d0
          rnd3 = (nodes_m(3,:) + nodes_m(4,:)) / 2.d0
          rnd4 = SUM(nodes_m,1) / 4.d0

          ! Coordinates of nodes of virtual new elements
          ALLOCATE (v__rr(5))

          ALLOCATE (v__rr(1) % points(3,2))
          v__rr(1) % points(1,:) = nodes_m(1,:)
          v__rr(1) % points(2,:) = rnd1
          v__rr(1) % points(3,:) = rnd4

          ALLOCATE (v__rr(2) % points(4,2))
          v__rr(2) % points(1,:) = rnd1
          v__rr(2) % points(2,:) = nodes_m(2,:)
          v__rr(2) % points(3,:) = rnd2
          v__rr(2) % points(4,:) = rnd4

          ALLOCATE (v__rr(3) % points(4,2))
          v__rr(3) % points(1,:) = rnd4
          v__rr(3) % points(2,:) = rnd2
          v__rr(3) % points(3,:) = nodes_m(3,:)
          v__rr(3) % points(4,:) = rnd3

          ALLOCATE (v__rr(4) % points(3,2))
          v__rr(4) % points(1,:) = rnd4
          v__rr(4) % points(2,:) = rnd3
          v__rr(4) % points(3,:) = nodes_m(4,:)

          ALLOCATE (v__rr(5) % points(3,2))
          v__rr(5) % points(1,:) = rnd4
          v__rr(5) % points(2,:) = nodes_m(4,:)
          v__rr(5) % points(3,:) = nodes_m(1,:)


     END SELECT

     END SUBROUTINE  Virtual_Elements



     FUNCTION  Q1_node_coords(rr_nodes, rr_rnd1, config)   RESULT(node_coords)
     !------------------------------------------------------------------------------- 
     IMPLICIT NONE

     REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: rr_nodes
     REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: rr_rnd1 
     CHARACTER(LEN=1),             INTENT(INOUT) :: config

     REAL(KIND=8), DIMENSION(2) :: node_coords

     REAL(KIND=8) :: alpha
     REAL(KIND=8) :: beta_3, beta_4, beta, beta_half
     REAL(KIND=8) :: gamma, delta, t
     REAL(KIND=8) :: num, den, pi = 3.141592654
     REAL(KIND=8) :: x_m43, y_m43, lenght
     REAL(KIND=8) :: xx1, xx2, xx3, xx4, xx5, xxm
     REAL(KIND=8) :: yy1, yy2, yy3, yy4, yy5, yym
     REAL(KIND=8) :: xx_1, xx_2, xx_3, &
                     xx_4, xx_5, xx_m
     REAL(KIND=8) :: yy_1, yy_2, yy_3, &
                     yy_4, yy_5, yy_m

     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: points
     REAL(KIND=8), DIMENSION(2)                  :: coords
     !-------------------------------------------------------------------------------

     x_m43 = (rr_nodes(4,1) + rr_nodes(3,1)) / 2.d0
     y_m43 = (rr_nodes(4,2) + rr_nodes(3,2)) / 2.d0

     num = rr_rnd1(2) - y_m43
     den = rr_rnd1(1) - x_m43
     
     IF (ABS(den) < 1.E-8)  den = 1.e-8
       
     alpha = ATAN2(num,den) 
       
     lenght = SQRT((rr_rnd1(2) - y_m43)**2 + (rr_rnd1(1) - x_m43)**2)


     ! Transforming coordinates in local reference system
     ! Origin of local system is point '5'. 'X' axis passes for points '5' 
     ! and 'm43'
     xx1 = rr_nodes(1,1) - rr_rnd1(1);  yy1 = rr_nodes(1,2) - rr_rnd1(2)
     xx2 = rr_nodes(2,1) - rr_rnd1(1);  yy2 = rr_nodes(2,2) - rr_rnd1(2)
     xx3 = rr_nodes(3,1) - rr_rnd1(1);  yy3 = rr_nodes(3,2) - rr_rnd1(2)
     xx4 = rr_nodes(4,1) - rr_rnd1(1);  yy4 = rr_nodes(4,2) - rr_rnd1(2)
     xxm =         x_m43 - rr_rnd1(1);  yym =         y_m43 - rr_rnd1(2)
     xx5 = 0.0;                         yy5 = 0.0
             
     xx_1 =  xx1*COS( pi-alpha ) - yy1*SIN( pi-alpha )
     yy_1 =  xx1*SIN( pi-alpha ) + yy1*COS( pi-alpha )
     
     xx_2 =  xx1*COS( pi-alpha ) - yy2*SIN( pi-alpha )
     yy_2 =  xx2*SIN( pi-alpha ) + yy2*COS( pi-alpha )
     
     xx_3 =  xx3*COS( pi-alpha ) - yy3*SIN( pi-alpha )
     yy_3 =  xx3*SIN( pi-alpha ) + yy3*COS( pi-alpha )
     
     xx_4 =  xx4*COS( pi-alpha ) - yy4*SIN( pi-alpha )
     yy_4 =  xx4*SIN( pi-alpha ) + yy4*COS( pi-alpha )
     
     xx_m =  xxm*COS( pi-alpha ) - yym*SIN( pi-alpha )       
     yy_m =  xxm*SIN( pi-alpha ) + yym*COS( pi-alpha )
     
     xx_5 = xx5
     yy_5 = yy5


     ! Computing quadrilater's angles in order to fix configuration:
     ! Configuration is determined by the intersection  of  two  straight  lines. 
     ! Line one passes for points '5' and 'm', line two passes  for  the  vertex
     ! associated to the minumun angle between '4' and '3', with slope  equal to
     ! half that angle
     ALLOCATE ( points(2,4) )
     points(1,1) = xx_1;         points(2,1) = yy_1
     points(1,2) = xx_2;         points(2,2) = yy_2
     points(1,3) = xx_3;         points(2,3) = yy_3
     points(1,4) = xx_4;         points(2,4) = yy_4
                     
     beta_3 = angle( points, 3 )
     beta_4 = angle( points, 4 )
     
     beta      = MIN( beta_3, beta_4 )
     beta_half = beta*0.5 
     
     IF ( beta .EQ. beta_3 ) THEN    
       points(1,1) = xx_5;         points(2,1) = yy_5
       points(1,2) = xx_2;         points(2,2) = yy_2
       points(1,3) = xx_3;         points(2,3) = yy_3
       points(1,4) = xx_m;         points(2,4) = yy_m        
        
       gamma = angle( points, 4)
       delta = pi - ( beta_half + gamma )
          
       t = ( TAN(-delta)*(xx_5 - xx_3) + yy_3 - yy_5 ) / &
           ( yy_m - yy_5 - TAN(-delta)*(xx_m - xx_5) )                                                 
     ELSE
       points(1,1) = xx_1;         points(2,1) = yy_1
       points(1,2) = xx_5;         points(2,2) = yy_5
       points(1,3) = xx_m;         points(2,3) = yy_m
       points(1,4) = xx_4;         points(2,4) = yy_4
       
       gamma = angle( points, 3)
       delta = pi - ( beta_half + gamma )
       
       t = ( TAN(delta)*(xx_5 - xx_4) + yy_4 - yy_5 ) / &
           ( yy_m - yy_5 - TAN(delta)*(xx_m - xx_5) )          
           
     ENDIF  

     IF ( t .GT. 1.0 ) THEN
     
       node_coords(1) = xx_5  
       node_coords(2) = yy_5                
     
     ELSEIF ( t .LT. 0.0 ) THEN
     
       node_coords(1) = xx_5
       node_coords(2) = yy_5
     
     ELSEIF ((t .GT. 0.d0) .AND. (t .LT. 1.d0)) THEN
     
       node_coords(1) = xx_5 + t*(xx_m - xx_5)  
       node_coords(2) = yy_5 + t*(yy_m - yy_5)
     
     ENDIF                           
     
     IF (SQRT((node_coords(2) - yy_5)**2 +           &
              (node_coords(1) - xx_5)**2 ) .LT.  1.E-1*lenght) THEN
       config = 'T' 
     ELSE
       config = 'Q'
     ENDIF  


     ! Transforming coordinates in global reference system
     coords(1) =  node_coords(1)*COS( pi-alpha ) + &
                  node_coords(2)*SIN( pi-alpha )
                
     coords(2) = -node_coords(1)*SIN( pi-alpha ) + &
                  node_coords(2)*COS( pi-alpha )

     node_coords(1) = coords(1) + rr_rnd1(1)
     node_coords(2) = coords(2) + rr_rnd1(2)


     END FUNCTION  Q1_node_coords

   END SUBROUTINE  Mark_for_Quality



 
  
   SUBROUTINE  Smooth_Grid (grid, l_its, relax_1, relax_2)

   ! This function is very badly programmed. This is because is has been modifyied
   ! adding pieces without reprogramming from scratch to insert new features.
   ! It wuold be useful a complete redesign.
   !---------------------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type), INTENT(INOUT) :: grid
   INTEGER,         INTENT(IN)    :: l_its
   REAL(KIND=8),    INTENT(IN)    :: relax_1, relax_2


   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dr
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: cm

   REAL(KIND=8), DIMENSION(3) :: x, y
   REAL(KIND=8), DIMENSION(2) :: x0

   REAL(KIND=8) :: t, u, delta,       &
                   delta_t, delta_u,  &
                   x1, x2, x3, x4,    &
                   y1, y2, y3, y4


   REAL(KIND=8), PARAMETER :: tol = 1.d-6

   INTEGER, DIMENSION(:), ALLOCATABLE :: p_edges, p_nodes

   INTEGER, DIMENSION(3) :: m_nodes

   INTEGER :: i, j, j_, k, m, &
              element_type,   &
              element_index,  &
              Nm_j, c_, c,    &
              count, l_out,   &
              m_out, Nj_j, n
   !---------------------------------------------------------------------------------------

   ALLOCATE (dr(grid % k_d, grid % Nj_d),  &
             cm(grid % k_d))   


   DO k = 1, l_its

      dr = 0.d0

      DO m = 1, grid % Nm_d

        ! Vector cm contains the coordinates of the center
        ! of mass for element m

        cm = 0.d0

        DO j_ = 1, SIZE(grid % j_m_d(m) % vec)

          j = grid % j_m_d(m) % vec(j_)
         
          cm = cm + grid % rr(:,j) / SIZE(grid % j_m_d(m) % vec)

        ENDDO


        ! For each node compute the sum of vectors joining 
        ! the node and the center of mass of adjacent elements
        ! The following loop determines the nodes displacement
        ! for the current laplacian iteration that will be
        ! adopted with smoothing formula 1.

        DO j_ = 1, SIZE(grid % j_m_d(m) % vec)

          j = grid % j_m_d(m) % vec(j_)

          dr(:,j) = dr(:,j) + cm(:) - grid % rr(:,j)

        ENDDO

      ENDDO


      ! Loop over nodes
      DO j = 1, grid % Nj_d

        IF (ANY(grid % jd_jb == SPREAD(j,1,grid % Nj_b))) THEN

          ! For boundary nodes displacement is set to zero.
          ! Set dr(:,j) = 0 for smoothing formula 1 and 
          ! omega = 0 for smoothing formula 2

          dr(:,j) = 0.d0

        ENDIF


        ! Number of elements sharing node j 
        Nm_j = SIZE(grid % m_j_d(j) % vec)

        ! Edges that are sides of the polygon. This
        ! vector will be used only after the test
        ! for quadrilaters.

        ALLOCATE (p_edges(Nm_j))


        ! Number of nodes of the patch of node j
        Nj_j = SIZE(grid % c_j(j) % vec)

        ALLOCATE (p_nodes(Nj_j))

        DO i = 1, Nj_j

          c = grid % c_j(j) % vec(i)

          IF (grid % j_c(1,c) /= j) THEN

            p_nodes(i) = grid % j_c(1,c)

          ELSE

            p_nodes(i) = grid % j_c(2,c)

          ENDIF 

        ENDDO


        ! This loop determines if one of the elements adjacent the 
        ! j-th node is a quadrilateral, if the prescribed 
        ! displacement (dr(:,j)) brings the node outside the 
        ! patch of elements sharing node j, and if the movement
        ! of node j inside the patch produces elements of negative
        ! area.
        ! If one of these three events occurs the movement of 
        ! node j is forbidden. Avoiding the second event is an artificial 
        ! trick to preserve positivity of elements' area, but
        ! it is necessary due to the fact that laplacian smoothing
        ! is not so reliable. Moreover this trick solves the problem
        ! of nonconvex patches with high skewness.

        loc_elm: DO count = 1, Nm_j
       
          element_index = grid % m_j_d(j) % vec(count)       
          element_type  = grid % ele_type_d(element_index)


          ! Quadrilater occurrence

          IF (element_type == 3) THEN

            dr(:,j) = 0.d0

            EXIT loc_elm

          ENDIF

          
          ! Veryfies that moving node j doesn't cause
          ! that a node of the patch finds itself inside
          ! element element_index of the patch as a consequence of
          ! the movement of internal rays (the lines that)
          ! join node j with nodes of the patch. Note
          ! that at this points elements are only triangles.

          ! Nodes of element element_index
          m_nodes = grid % j_m_d(element_index) % vec
          
          DO i = 1, SIZE(m_nodes)

            IF (m_nodes(i) == j) THEN

              x(i) = grid % rr(1,j) + dr(1,j)/relax_1
              y(i) = grid % rr(2,j) + dr(2,j)/relax_1              

            ELSE

              x(i) = grid % rr(1, m_nodes(i))
              y(i) = grid % rr(2, m_nodes(i))

            ENDIF 

          ENDDO 

          n = SIZE(m_nodes) 

          DO i = 1, SIZE(p_nodes)

            IF (ANY(m_nodes == SPREAD(p_nodes(i),1,SIZE(m_nodes)))) CYCLE

            x0 = grid % rr(:,p_nodes(i))

            CALL locpt (x0(1), x0(2), x, y, n, l_out, m_out)

            IF (l_out >= 0) THEN
  
              dr(:,j) = 0.d0

              EXIT loc_elm

            ENDIF 

          ENDDO


          ! Extract the edges that constitute the sides of
          ! the polygon centered in j. Element is necessarily
          ! a triangle.

          DO c_ = 1, SIZE(grid % c_m(element_index) % vec)
          
            c = grid % c_m(element_index) % vec(c_)
            
            IF (grid % j_c(1,c) /= j  .AND.  &
                grid % j_c(2,c) /= j) THEN

              p_edges(count) = c
              
              EXIT
                        
            ENDIF 
          
          ENDDO

        ENDDO loc_elm


        ! This condition is adopted to treat the case in which
        ! one quadrilater is found in the patch. In this case 
        ! no more processing for node j is needed, and the 
        ! following instructions are jumped. (ONLY 2D damned !)

        IF (ANY(dr(:,j) /= 0.d0)) THEN

          x1 = grid % rr(1,j)
          y1 = grid % rr(2,j)

          x3 = grid % rr(1,j) + dr(1,j)/relax_1
          y3 = grid % rr(2,j) + dr(2,j)/relax_1


          ! Loop over edges that are sides of the patch. The
          ! goal is to determine if the line through 1-3 
          ! intersects the line through the two nodes of the
          ! sides of the patch in such a way that corresponds
          ! to the location of point 3 outside the patch.
          ! Node 3 is node j displaced with dr(:,j)

          DO c_ = 1, SIZE(p_edges)
        
            c = p_edges(c_)
        
            x2 = grid % rr(1,grid % j_c(1,c))
            y2 = grid % rr(2,grid % j_c(1,c))
            
            x4 = grid % rr(1,grid % j_c(2,c))
            y4 = grid % rr(2,grid % j_c(2,c))


            ! Solving the system for the intersection of the lines
            ! through 1-3 and 2-4 in parametric form. t is the parameter
            ! for line 1-3, while s is the parameter for line 2-4.
            
            delta = (y4 - y2)*(x3 - x1) - (x4 -x2)*(y3 - y1)
            
            delta_t = (x2 - x1)*(y4 - y2) - (y2 - y1)*(x4 - x2)
            delta_u = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)
            
            t = delta_t / (delta + 1.e-15)
            u = delta_u / (delta + 1.e-15)


            ! If the movement of the node j brings node outside
            ! the polygon centered in j, then the displacement
            ! is forced to zero.
            
            IF ((u >= 0.d0  .AND.  u <= 1.d0)  .AND.  &
                (t >= 0.d0  .AND.  t <= 1.d0*relax_2)) THEN
              
              dr(:,j) = 0.d0

              EXIT
                
            ENDIF

          ENDDO
        
        ENDIF

        grid % rr(:,j) = grid % rr(:,j) + dr(:,j)/relax_1

        DEALLOCATE (p_edges, p_nodes)

      ENDDO              

!       e = 0.d0
!       DO j = 1, grid % Nj_d
!         e = e + SUM(dr(:,j)**2)
!       ENDDO
! 
!       IF (e .LT. tol)   EXIT
     
   ENDDO         

   DEALLOCATE (dr, cm)

   END SUBROUTINE  Smooth_Grid


 
 
 !____________________________________________________________________________

   SUBROUTINE  Swap_Edges (grid)

 ! Valid only in 2D and for triangles 
 ! (Please forgive me!!!)
 !____________________________________________________________________________

   IMPLICIT NONE

   TYPE(grid_type), INTENT(INOUT) :: grid

   REAL(KIND=8), DIMENSION(2) :: v1, v2, v3

   REAL(KIND=8) :: inR_m1, inR_m2, inR_m3, inR_m4, t, s, &
                   Delta, Delta_t, Delta_s, theta2, theta3

   INTEGER, DIMENSION(:), ALLOCATABLE :: n_, local_edges

   INTEGER, DIMENSION(4) :: n, e
   INTEGER, DIMENSION(3) :: j_m_d_cc_order

   INTEGER :: i, j, k, c, m1, m2, n1_ei, n2_ei, m_j, check

 !----------------------------------------------------------------------------


   CALL Edge_Structure_Gen(grid, 'C') 

 
   DO c = 1, grid % Nc_d

     IF (ANY((grid % cd_cb) == SPREAD(c, 1, grid % Nc_b))) CYCLE

     m1 = grid % m_c(c) % vec(1)
     m2 = grid % m_c(c) % vec(2)

     IF (grid % ele_type_d(m1) == 2  .AND.  grid % ele_type_d(m2) == 2) THEN

       ! Define n1, n2, n3, n4 as the vertices  of  the  quadrilater 
       ! formed by the two triangles.  The four nodes are ordered in
       ! counterclockwise sense over the quadrilater.

       v1(1) = grid % rr(1,grid % j_m_d(m1) % vec(1)); v1(2) = grid % rr(2,grid % j_m_d(m1) % vec(1))
       v2(1) = grid % rr(1,grid % j_m_d(m1) % vec(2)); v2(2) = grid % rr(2,grid % j_m_d(m1) % vec(2))
       v3(1) = grid % rr(1,grid % j_m_d(m1) % vec(3)); v3(2) = grid % rr(2,grid % j_m_d(m1) % vec(3))
       v2 = v2-v1; v3 = v3-v1

       j_m_d_cc_order(1) = grid % j_m_d(m1) % vec(1)
       theta3 = ATAN2(v3(2),v3(1)); IF(theta3<0d0) theta3 = 2*3.14159265359d0+theta3
       theta2 = ATAN2(v2(2),v2(1)); IF(theta2<0d0) theta2 = 2*3.14159265359d0+theta2
       IF( theta3 - theta2 >= 0d0) THEN
         ! Counterclockwise order: 1 - 2 - 3
         j_m_d_cc_order(2) = grid % j_m_d(m1) % vec(2)
         j_m_d_cc_order(3) = grid % j_m_d(m1) % vec(3)
       ELSE
         ! Counterclockwise order: 1 - 3 - 2
         j_m_d_cc_order(2) = grid % j_m_d(m1) % vec(3)
         j_m_d_cc_order(3) = grid % j_m_d(m1) % vec(2)
!	WRITE(*,*) 'c order'
       ENDIF

       IF (j_m_d_cc_order(1) /= grid % j_c(1,c)  .AND.  &
           j_m_d_cc_order(1) /= grid % j_c(2,c)) THEN

         n(1) = j_m_d_cc_order(1)
         n(2) = j_m_d_cc_order(2)
         n(4) = j_m_d_cc_order(3)

       ELSEIF (j_m_d_cc_order(2) /= grid % j_c(1,c)  .AND.  &
               j_m_d_cc_order(2) /= grid % j_c(2,c)) THEN
 
         n(1) = j_m_d_cc_order(2)
         n(2) = j_m_d_cc_order(3)
         n(4) = j_m_d_cc_order(1)

       ELSE

         n(1) = j_m_d_cc_order(3)
         n(2) = j_m_d_cc_order(1)
         n(4) = j_m_d_cc_order(2)

       ENDIF


       IF (grid % j_m_d(m2) % vec(1) /= grid % j_c(1,c)  .AND. &
           grid % j_m_d(m2) % vec(1) /= grid % j_c(2,c)) THEN

           n(3) = grid % j_m_d(m2) % vec(1)

       ELSEIF (grid % j_m_d(m2) % vec(2) /= grid % j_c(1,c)  .AND. &
               grid % j_m_d(m2) % vec(2) /= grid % j_c(2,c)) THEN
 
           n(3) = grid % j_m_d(m2) % vec(2)

       ELSE 
      
           n(3) = grid % j_m_d(m2) % vec(3)

       ENDIF


       ! Verify the convexity of the quadrilater
       Delta   = (grid % rr(1,n(1)) - grid % rr(1,n(3))) * (grid % rr(2,n(4)) - grid % rr(2,n(2))) -  &
                 (grid % rr(1,n(4)) - grid % rr(1,n(2))) * (grid % rr(2,n(1)) - grid % rr(2,n(3)))
               
       Delta_t = (grid % rr(1,n(4)) - grid % rr(1,n(3))) * (grid % rr(2,n(4)) - grid % rr(2,n(2))) -  &
                 (grid % rr(1,n(4)) - grid % rr(1,n(2))) * (grid % rr(2,n(4)) - grid % rr(2,n(3)))
               
       Delta_s = (grid % rr(1,n(1)) - grid % rr(1,n(3))) * (grid % rr(2,n(4)) - grid % rr(2,n(3))) -  &
                 (grid % rr(1,n(4)) - grid % rr(1,n(3))) * (grid % rr(2,n(1)) - grid % rr(2,n(3)))    
                 
       t = Delta_t / (Delta + 1.e-10)
       s = Delta_s / (Delta + 1.e-10)


       ! Quadrilater is not convex. Swapping is not allowed.
       IF (t .LT. 0.d0  .OR.  t .GT. 1.d0  .OR.  &
           s .LT. 0.d0  .OR.  s .GT. 1.d0) CYCLE


       ! Define vector  e  as the  egdes  of the quadrilater 
       ! formed by the two triangles. The four edges are or-
       ! dered in clockwise sense over the quadrilater.
       !
       ! - e(1) -> nodes n1, n4
       ! - e(2) -> nodes n4, n3
       ! - e(3) -> nodes n3, n2
       ! - e(4) -> nodes n2, n1
       !
       ! Moreover the index of the loop c, identifies the edge
       ! corresponding to one diagonal of the quadrilater.


       e(1) = Extract_np_index(grid % c_j, n(1), n(4))
       e(2) = Extract_np_index(grid % c_j, n(4), n(3))
       e(3) = Extract_np_index(grid % c_j, n(3), n(2))
       e(4) = Extract_np_index(grid % c_j, n(2), n(1))

       ! Up to now the four possible triangles are:
       ! 
       ! - m1 ->   nodes n1, n2, n4;   edges c, e(4), e(1)  (original)
       ! - m2 ->   nodes n3, n4, n2;   edges c, e(2), e(3)  (original)
       ! - m3 ->   nodes n1, n3, n4;   edges c, e(1), e(2)  
       ! - m4 ->   nodes n1, n2, n3;   edges c, e(3), e(4)


       ! Computing for each triangle the inradius

       ! Triangle m1     
       v1(1) = grid % rr(1,n(1)); v1(2) = grid % rr(2,n(1))
       v2(1) = grid % rr(1,n(2)); v2(2) = grid % rr(2,n(2))
       v3(1) = grid % rr(1,n(4)); v3(2) = grid % rr(2,n(4))

       inR_m1 = InRadius(v1, v2, v3)


       ! Triangle m2
       v1(1) = grid % rr(1,n(3)); v1(2) = grid % rr(2,n(3))
       v2(1) = grid % rr(1,n(4)); v2(2) = grid % rr(2,n(4))
       v3(1) = grid % rr(1,n(2)); v3(2) = grid % rr(2,n(2))

       inR_m2 = InRadius(v1, v2, v3)


       ! Triangle m3     
       v1(1) = grid % rr(1,n(1)); v1(2) = grid % rr(2,n(1))
       v2(1) = grid % rr(1,n(3)); v2(2) = grid % rr(2,n(3))
       v3(1) = grid % rr(1,n(4)); v3(2) = grid % rr(2,n(4))

       inR_m3 = InRadius(v1, v2, v3)


       ! Triangle m4     
       v1(1) = grid % rr(1,n(1)); v1(2) = grid % rr(2,n(1))
       v2(1) = grid % rr(1,n(2)); v2(2) = grid % rr(2,n(2))
       v3(1) = grid % rr(1,n(3)); v3(2) = grid % rr(2,n(3))

       inR_m4 = InRadius(v1, v2, v3)


       ! The couple of triangles to which is associated the minimum
       ! radius is removed by swapping the edge and the new elements
       ! m3 and m4 are created.
       IF (MINVAL((/inR_m1, inR_m2, inR_m3, inR_m4/)) == inR_m1  .OR.  &
           MINVAL((/inR_m1, inR_m2, inR_m3, inR_m4/)) == inR_m2) THEN

         ! Element m1 is transformed in element m3
         grid % j_m_d(m1) % vec(1) = n(1)
         grid % j_m_d(m1) % vec(2) = n(3)
         grid % j_m_d(m1) % vec(3) = n(4)

         ! Element m2 is transformed in element m4
         grid % j_m_d(m2) % vec(1) = n(1)
         grid % j_m_d(m2) % vec(2) = n(2)
         grid % j_m_d(m2) % vec(3) = n(3)


         ! Modifying j_c
         grid % j_c(1,c) = n(1)
         grid % j_c(2,c) = n(3)


         ! Modifying c_j
         
         ! Considering nodes n1 and n3. Adding edge c 
         ! to connecitvity
         DO i = 1, 4, 2
         
           ALLOCATE (local_edges(SIZE(grid % c_j(n(i)) % vec)))
           
           local_edges = grid % c_j(n(i)) % vec

           DEALLOCATE (grid % c_j(n(i)) % vec)
           ALLOCATE (grid % c_j(n(i)) % vec(SIZE(local_edges) + 1))

           grid % c_j(n(i)) % vec(1) = c
           grid % c_j(n(i)) % vec(2:SIZE(grid % c_j(n(i)) % vec)) = local_edges
           
           DEALLOCATE (local_edges)
           
         ENDDO

         ! Considering nodes n2 and n4. Removing edge c 
         ! from connecitvity
         DO i = 2, 4, 2

           ALLOCATE (local_edges(SIZE(grid % c_j(n(i)) % vec)))
           
           local_edges = grid % c_j(n(i)) % vec

           DEALLOCATE (grid % c_j(n(i)) % vec)
           ALLOCATE (grid % c_j(n(i)) % vec(SIZE(local_edges) - 1))

           k = 1
           DO j = 1, SIZE(local_edges)
           
             IF (local_edges(j) /= c) THEN
             
               grid % c_j(n(i)) % vec(k) = local_edges(j)
               
               k = k + 1
             
             ENDIF
           
           ENDDO
           
           DEALLOCATE (local_edges)
           
         ENDDO


         ! Modifying m_c
         DO i = 1, 4

           ! Nodes defining edge e(i)
           n1_ei = grid % j_c(1,e(i))
           n2_ei = grid % j_c(2,e(i))

           ! Loop over elements adjacent edge e(i) (not yet updated to swap)
           DO j = 1, SIZE(grid % m_c(e(i)) % vec)

             m_j = grid % m_c(e(i)) % vec(j)
             
             ! Nodes of the element m_j (updated to swap)
             ALLOCATE (n_(SIZE(grid % j_m_d(m_j) % vec)))

             n_ = grid % j_m_d(m_j) % vec

             check = 0

             ! Check that two nodes of the elements coincide with
             ! the nodes of the selected edge. So to confirm adj-
             ! acency variable check must score 2.
             DO k = 1, SIZE(n_)

               IF (n_(k) == n1_ei  .OR. &
                   n_(k) == n2_ei) check = check + 1

             ENDDO

             IF (check /= 2  .AND.  m1 == m_j) grid % m_c(e(i)) % vec(j) = m2
             IF (check /= 2  .AND.  m2 == m_j) grid % m_c(e(i)) % vec(j) = m1

             DEALLOCATE (n_)

           ENDDO

         ENDDO

       ENDIF

     ENDIF

   ENDDO 

   END SUBROUTINE  Swap_Edges




 !_____________________________________________________________________
 
   FUNCTION  InRadius (v1, v2, v3)  RESULT (radius)
 
 ! Computes the radius of the circle inscribed in a triangle.
 ! The triangle is specified by its vertices' coordinates
 !_____________________________________________________________________

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(2), INTENT(IN) :: v1, v2, v3

   REAL(KIND=8) :: radius


   REAL(KIND=8) :: l12, l23, l31, s1, s2, s3
 
 !---------------------------------------------------------------------


   l12 = SQRT((v1(1) - v2(1))**2 + (v1(2) - v2(2))**2)
   l23 = SQRT((v2(1) - v3(1))**2 + (v2(2) - v3(2))**2)
   l31 = SQRT((v3(1) - v1(1))**2 + (v3(2) - v1(2))**2)

   s1 =  l12 + l23 - l31
   s2 =  l12 - l23 + l31
   s3 = -l12 + l23 + l31

   ! Verifies that the  triangle is not degenerate (i.e.
   ! tre vertices aligned  along the same line). If the
   ! triangle is degenerate then radius is assigned the
   ! irrealistic value of -1.
   IF (s1 .LT. 1.e-10  .OR.  &
       s2 .LT. 1.e-10  .OR.  &
       s3 .LT. 1.e-10) THEN

       radius = -1

       RETURN

   ENDIF

   radius = 0.5 * SQRT(( l12 + l23 - l31) *  &
                       ( l12 - l23 + l31) *  &
                       (-l12 + l23 + l31) / (l12 + l23 + l31))

   END FUNCTION InRadius





 !__________________________________________________________________________
 
   SUBROUTINE  Write_Edges_Length (length)
 !__________________________________________________________________________
 
   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:) :: length
 
 
   INTEGER :: k, n
   INTEGER, DIMENSION(1) :: max_pos

 !--------------------------------------------------------------------------

 
   OPEN (UNIT=1000, FILE='EDGES_LENGTH.dat', FORM='formatted')
    
     n = SIZE(length)
    
     DO k = 1, n
    
       WRITE(1000,*) k, MAXVAL(length)
       max_pos = MAXLOC(length)
       length(max_pos) = 0.0
      
     ENDDO  
    
   CLOSE(1000) 

   END SUBROUTINE  Write_Edges_Length

      
END MODULE grid_utils
