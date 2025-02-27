MODULE  reconstruction_proc

 !-------------------------------------------------------------------------
 !   Description: Contains subroutines for driving the reconstruction of
 !                of grid topology. The reconstruction at level M is ob-
 !                tained starting from level  0  through  M intermediate 
 !                steps.
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

 USE  csr
 USE  csr_pair_sys
 USE  solve_skit
 USE  dynamic_vector
 USE  grid_utils
 USE  io_proc
 USE  mesh_structure
 USE  rh_setting
 USE  structures 
 USE  triangulate_proc

 INTEGER, PARAMETER :: NO_TREAT = 0, &
                       MOVE_NDS = 1, &
                   MOVE_NDS_MSH = 2 

 CONTAINS
 
   FUNCTION  Reconstruct_Mesh (grid_zero, refinement_history, params, &
                               sol, add_sol)   RESULT (adapted_grid)
   !-----------------------------------------------------------------------------------------------
   IMPLICIT NONE
  
   TYPE(grid_type),                   INTENT(IN)    :: grid_zero
   TYPE(adaption_param),              INTENT(IN)    :: params
   TYPE(solution_type),               INTENT(IN)    :: sol ! solution on grid 0
   TYPE(solution_type), DIMENSION(:), INTENT(INOUT) :: add_sol

   TYPE(grid_type) :: adapted_grid

 
   TYPE(E_R_P), DIMENSION(0:params % entry_level) :: refinement_history

   TYPE(D_I_V_R), DIMENSION(:), POINTER :: displ
  
   TYPE(D_I_V), DIMENSION(:), ALLOCATABLE :: T_j_m_d, T_m_j_d

   TYPE(CSR_matrix) :: KK

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tmp_ww
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: d_u, F 

   REAL(KIND=8), DIMENSION(3) :: v1, v2, v3

   REAL(KIND=8) :: inR_m1, inR_m2, inR_m3, inR_m4

   INTEGER, DIMENSION(:), ALLOCATABLE :: ins_nodes_index

   INTEGER, DIMENSION(4) :: n

   INTEGER  :: i, j, k, p, nds_to_add, bnds_to_add,  &
               node, node1, node2, edge, T_Nm_d, t, s_ind

   LOGICAL, DIMENSION(:), ALLOCATABLE :: np_refine, bnp_refine
   !-----------------------------------------------------------------------------------------------

   adapted_grid = Duplicate_Grid(grid_zero)

   IF (refinement_history(0) % N_refined_edges == 0) THEN
     
     PRINT*, '____________________________________'
     PRINT*, ''
     PRINT*, 'OPTIME.EXE - No adaption required.  '
     PRINT*, '____________________________________'    
   
   ELSE    
  
     ! Dummy allocation
     ALLOCATE (T_j_m_d(1), T_m_j_d(1))
  
     PRINT*, ''

     DO i = 0, params % entry_level 

       adapted_grid % adaptation_Level = i+1
      
       ALLOCATE (np_refine(adapted_grid % Nc_d), bnp_refine(adapted_grid % Nc_b))
       ALLOCATE (ins_nodes_index(refinement_history(i) % N_refined_edges))
      
       np_refine  = .FALSE.
       bnp_refine = .FALSE.
    
    
       CALL Build_Refine_vec(adapted_grid, refinement_history(i),  &
                             nds_to_add, bnds_to_add, np_refine,  bnp_refine)

       ! Operations required for solution interpolation
       s_ind = SIZE(sol % ww,2)

       DO j = 0, i - 1, +1
         s_ind = s_ind + SIZE(add_sol(j+1) % ww,2)
       ENDDO

       ALLOCATE (tmp_ww(SIZE(sol % ww,1), s_ind))
       tmp_ww(:, 1:SIZE(sol % ww,2)) = sol % ww

       ! Reset s_ind to the value equal to the 
       ! number of nodes in grid 0
       s_ind = SIZE(sol % ww,2)

       DO j = 0, i - 1, +1

         IF (j > 0)  s_ind = s_ind + SIZE(add_sol(j) % ww,2)
         tmp_ww(:,s_ind + 1 : s_ind + SIZE(add_sol(j+1) % ww,2)) = add_sol(j+1) % ww

       ENDDO


       ! Adding nodes and reconstructiong topology
       CALL Triangulate(np_refine,  bnp_refine,               &
                        nds_to_add, bnds_to_add,              &
                        adapted_grid, refinement_history(i),  &
                        params % bou_tre_typ, displ, tmp_ww,  &
                        add_sol(i+1) % ww)

       DEALLOCATE (tmp_ww)

       IF (ALLOCATED(np_refine))   DEALLOCATE(np_refine, bnp_refine, ins_nodes_index)

       write (*,'(4x,a4,1x,i1,1x,a5,1x,i8,a8,1x,i8)') 'Grid', i+1, 'nodes', adapted_grid % Nj_d, ' elements', adapted_grid % Nm_d


       ! Boundary Treatment
       SELECT CASE (params % bou_tre_typ)

         CASE (NO_TREAT)

           ! Nodes are left at the middle point of each refined edge. Such an
           ! operation has already been performed in subroutine ADD_NODES  of
           ! module NEW_NODES.



         CASE (MOVE_NDS)

           ! Move ONLY new BOUNDARY NODES towards real boundary
           DO p = 1, SIZE(adapted_grid % rr, 2)
             adapted_grid % rr(1, p) = adapted_grid % rr(1, p) + displ(p) % vec(1)
             adapted_grid % rr(2, p) = adapted_grid % rr(2, p) + displ(p) % vec(2)
           ENDDO



         CASE (MOVE_NDS_MSH)

           ! Relaxing MESH according to BOUNDARY NODES movement
           ! towards real boundary

           ! Transforming quadrilaters in triangles
           T_Nm_d = adapted_grid % Nm_d + COUNT(adapted_grid % ele_type_d == 3)

           IF (ALLOCATED(T_j_m_d))  THEN
!              DO j = 1, SIZE(T_j_m_d)
!                IF (ASSOCIATED(T_j_m_d(j) % vec))  DEALLOCATE (T_j_m_d(j) % vec)
!              ENDDO
             DEALLOCATE (T_j_m_d)
           ENDIF

           ALLOCATE (T_j_m_d(T_Nm_d))

           DO p = 1, T_Nm_d
             ALLOCATE (T_j_m_d(p) % vec(3))
           ENDDO

           t = 1

           DO p = 1, adapted_grid % Nm_d

            IF (adapted_grid % ele_type_d(p) == 2) THEN

              T_j_m_d(t) % vec(1) = adapted_grid % j_m_d(p) % vec(1)
              T_j_m_d(t) % vec(2) = adapted_grid % j_m_d(p) % vec(2)
              T_j_m_d(t) % vec(3) = adapted_grid % j_m_d(p) % vec(3)

              t = t + 1      

            ELSEIF (adapted_grid % ele_type_d(p) == 3) THEN

              ! Select the couple of triangles characterized by
              ! the lower skweness. The triangles are selected 
              ! on the basis of the radius of the inscribed
              ! circumpherence (same criterion adopted to swap edges)

              n(1) =  adapted_grid % j_m_d(p) % vec(1)
              n(2) =  adapted_grid % j_m_d(p) % vec(2)
              n(3) =  adapted_grid % j_m_d(p) % vec(3)
              n(4) =  adapted_grid % j_m_d(p) % vec(4)

              
              v1(1) = adapted_grid % rr(1,n(1)); v1(2) = adapted_grid % rr(2,n(1))
              v2(1) = adapted_grid % rr(1,n(2)); v2(2) = adapted_grid % rr(2,n(2))
              v3(1) = adapted_grid % rr(1,n(4)); v3(2) = adapted_grid % rr(2,n(4))

              inR_m1 = InRadius(v1, v2, v3)


              v1(1) = adapted_grid % rr(1,n(3)); v1(2) = adapted_grid % rr(2,n(3))
              v2(1) = adapted_grid % rr(1,n(4)); v2(2) = adapted_grid % rr(2,n(4))
              v3(1) = adapted_grid % rr(1,n(2)); v3(2) = adapted_grid % rr(2,n(2))

              inR_m2 = InRadius(v1, v2, v3)


              v1(1) = adapted_grid % rr(1,n(1)); v1(2) = adapted_grid % rr(2,n(1))
              v2(1) = adapted_grid % rr(1,n(3)); v2(2) = adapted_grid % rr(2,n(3))
              v3(1) = adapted_grid % rr(1,n(4)); v3(2) = adapted_grid % rr(2,n(4))

              inR_m3 = InRadius(v1, v2, v3)


              v1(1) = adapted_grid % rr(1,n(1)); v1(2) = adapted_grid % rr(2,n(1))
              v2(1) = adapted_grid % rr(1,n(2)); v2(2) = adapted_grid % rr(2,n(2))
              v3(1) = adapted_grid % rr(1,n(3)); v3(2) = adapted_grid % rr(2,n(3))

              inR_m4 = InRadius(v1, v2, v3)


              ! The couple of triangles to which is associated the minimum
              ! radius is discarded and the altenative coulpe is selected
              
              IF (MINVAL((/inR_m1, inR_m2, inR_m3, inR_m4/)) == inR_m1  .OR.  &
                  MINVAL((/inR_m1, inR_m2, inR_m3, inR_m4/)) == inR_m2) THEN


                T_j_m_d(t) % vec(1) = adapted_grid % j_m_d(p) % vec(1)
                T_j_m_d(t) % vec(2) = adapted_grid % j_m_d(p) % vec(2)
                T_j_m_d(t) % vec(3) = adapted_grid % j_m_d(p) % vec(3)

                t = t + 1

                T_j_m_d(t) % vec(1) = adapted_grid % j_m_d(p) % vec(1)
                T_j_m_d(t) % vec(2) = adapted_grid % j_m_d(p) % vec(3)
                T_j_m_d(t) % vec(3) = adapted_grid % j_m_d(p) % vec(4)

                t = t + 1

              ELSE
              
                T_j_m_d(t) % vec(1) = adapted_grid % j_m_d(p) % vec(1)
                T_j_m_d(t) % vec(2) = adapted_grid % j_m_d(p) % vec(2)
                T_j_m_d(t) % vec(3) = adapted_grid % j_m_d(p) % vec(4)

                t = t + 1

                T_j_m_d(t) % vec(1) = adapted_grid % j_m_d(p) % vec(2)
                T_j_m_d(t) % vec(2) = adapted_grid % j_m_d(p) % vec(3)
                T_j_m_d(t) % vec(3) = adapted_grid % j_m_d(p) % vec(4)

                t = t + 1

              ENDIF

            ENDIF

           ENDDO

           IF (ALLOCATED(T_m_j_d)) THEN
!              DO j = 1, SIZE(T_m_j_d)
!                IF (ASSOCIATED(T_m_j_d(j) % vec))  DEALLOCATE (T_m_j_d(j) % vec)
!              ENDDO
             DEALLOCATE (T_m_j_d)
           ENDIF

           ALLOCATE (T_m_j_d(size_DIV(T_j_m_d, 3))) 

           T_m_j_d = invert_DIV(T_j_m_d)

           ALLOCATE (d_u(2*adapted_grid % Nj_d), F(2*adapted_grid % Nj_d))

           ! No force is applied to nodes
           F = 0.d0


           ! WARNING. Matrice rigidezza associata alla griglia
           ! ridotta a soli triangoli
           KK = Stiffness_Matrix(adapted_grid, T_j_m_d, T_m_j_d, params)

           DO j = 1, SIZE(T_j_m_d)
             IF (ASSOCIATED(T_j_m_d(j) % vec))  DEALLOCATE (T_j_m_d(j) % vec)
           ENDDO
           DEALLOCATE (T_j_m_d)
 
           DO j = 1, SIZE(T_m_j_d)
             IF (ASSOCIATED(T_m_j_d(j) % vec))  DEALLOCATE (T_m_j_d(j) % vec)
           ENDDO
           DEALLOCATE (T_m_j_d)



           ! Dirichlet B.C.
           CALL Dirichlet_BoCo(adapted_grid, displ, KK, F)

           ! Nodes Displacement:                                        
           d_u = 0                                                        


           OPEN(UNIT = 101, FILE = 'spkit.cfg')
             CALL init_spkit(101)
           CLOSE(101)

           CALL solve_spkit(KK, F, d_u)                                   

           DEALLOCATE (F)

           DO p = 1, SIZE(adapted_grid % rr, 2)
             adapted_grid % rr(1, p) = adapted_grid % rr(1, p) + d_u(2*(p-1) + 1)
             adapted_grid % rr(2, p) = adapted_grid % rr(2, p) + d_u(2*(p-1) + 2)
           ENDDO

           DEALLOCATE (d_u)           

       END SELECT


       ! Swapping edges (ONLY for 2D meshes)
       IF (params % swap)   CALL Swap_Edges(adapted_grid)


       CALL Edge_structure_gen(adapted_grid, 'C')
       
       CALL Invert_jd_jb(adapted_grid)

       ! Inverting j_m_d
       IF (ASSOCIATED(adapted_grid % m_j_d))  DEALLOCATE(adapted_grid % m_j_d)
       ALLOCATE (adapted_grid % m_j_d(size_DIV(adapted_grid % j_m_d, 3)))

       adapted_grid % m_j_d = invert_DIV(adapted_grid % j_m_d)

       ! Inverting j_m_b
       IF (ASSOCIATED(adapted_grid % m_j_b))  DEALLOCATE(adapted_grid % m_j_b)
       ALLOCATE (adapted_grid % m_j_b(size_DIV(adapted_grid % j_m_b,3)))

       adapted_grid % m_j_b = invert_DIV(adapted_grid % j_m_b)

  
       ! Updating Refinement_history%adjacent_nodes:
       IF (ASSOCIATED (refinement_history(i) % adjacent_nodes)) &
           DEALLOCATE (refinement_history(i) % adjacent_nodes)

       ALLOCATE (refinement_history(i) % adjacent_nodes(refinement_history(i) % N_refined_edges))
        
       DO j = 1, refinement_history(i) % N_refined_edges
    
         node = refinement_history(i) % inserted_nodes(j)
         ALLOCATE (refinement_history(i) % adjacent_nodes(j) % vec(SIZE(adapted_grid % c_j(node) % vec)))
              
         DO k = 1, SIZE(adapted_grid%c_j(node) % vec)
        
           edge  = adapted_grid % c_j(node) % vec(k)

           node1 = adapted_grid % j_c(1,edge)
           node2 = adapted_grid % j_c(2,edge)

           IF (node1 /= node) refinement_history(i) % adjacent_nodes(j) % vec(k) = node1
           IF (node2 /= node) refinement_history(i) % adjacent_nodes(j) % vec(k) = node2      
           
         ENDDO
    
       ENDDO


       ! Smoothing grid
       IF (params % smooth)   CALL Smooth_grid(adapted_grid, params % l_its,  &
                                               params % relax_1, params % relax_2)

     ENDDO  
  
     CALL Write_Adaption_History(refinement_history, adapted_grid % adaptation_level,  &
                                 params % entry_level)

   ENDIF


   CONTAINS

  
     FUNCTION  Stiffness_Matrix (mesh, T_j_m_d, T_m_j_d, params)  RESULT (KK)
     !--------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(grid_type),           INTENT(INOUT) :: mesh 
     TYPE(D_I_V), DIMENSION(:), INTENT(IN) :: T_j_m_d, T_m_j_d
     TYPE(adaption_param),      INTENT(IN) :: params

     TYPE(CSR_matrix) :: KK


     REAL(KIND=8), DIMENSION(6,6) :: KK_e                                        
     REAL(KIND=8), DIMENSION(36)  :: h                                           
     REAL(KIND=8), DIMENSION(3)   :: x, y, l                                     

     REAL(KIND=8) :: Jou_Mod, Poi_Coe, ke_i

     INTEGER, DIMENSION(mesh % Nj_d) :: dim
     INTEGER, DIMENSION(36) :: r, c
     INTEGER, DIMENSION(6)  :: n
                                                                                 
     INTEGER  ::  e, i, nnz,     &
                  j, nrow,       &
                  idf = 11,      &
                  row_i, col_i,  & 
                  ib, idx_, idx

     REAL(KIND=8) :: x_m12, x_m23, x_m31, &
                     y_m12, y_m23, y_m31, &
                     x_g, y_g,            &
                     G_M12, G_M23, G_M31

     INTEGER, PARAMETER :: CONSTANT   = 0, &
                           INV_LENGHT = 1
     !--------------------------------------------------------------------------------

     nrow =  2 * mesh % Nj_d
     dim  = 0                                                                   

     ! Loop over domain nodes number to extract for                             
     ! each node the number of adjacent elements + 1                            
     DO i = 1, SIZE(T_m_j_d)
          dim(i) = SIZE(T_m_j_d(i) % vec) + 1
     ENDDO

     ! Added by N2-DM suggestion
     DO i = 1, mesh % Nj_b
        ib = mesh % jd_jb(i)
        dim(ib) = SIZE(T_m_j_d(ib) % vec) + 2
     ENDDO

     ! Number of non zero elements of matrix KK
     nnz = 4*(SUM(dim))


     ! Stiffness Matrix
     ALLOCATE (KK % i(nrow + 1), KK % j(nnz), KK % e(nnz))
!     ALLOCATE (iw(nrow))

     CALL set_CSR_sys(0.d0, KK)


     ! Setting stiffness matrix nonzero elements
     DO i = 1, mesh % Nj_d
       KK % i(mesh % k_d*(i - 1) + 1) = mesh % k_d * dim(i)
       KK % i(mesh % k_d*(i - 1) + 2) = mesh % k_d * dim(i)
     ENDDO
    
     CALL setia(KK % i)
     CALL set_CSR_col(mesh % k_d, T_j_m_d, T_m_j_d,   KK)



     OPEN (UNIT = idf, FILE = 'stiff_mat.out')

       KK % e = 0.d0

       ! Loop over elements                                                     
       DO e = 1, SIZE(T_j_m_d)

         h = 0;   r = 0;   c = 0

         DO i = 1, 3                                              
           n((i-1)*2 + 1) = 2*((T_j_m_d(e) % vec(i)) - 1) + 1                   
           n((i-1)*2 + 2) = 2*((T_j_m_d(e) % vec(i)) - 1) + 2                   
         ENDDO


         ! Joung's Modulus
         SELECT CASE (params % jou_mod_fml)
         
           CASE (CONSTANT)

              Jou_Mod = params % jmf_const


           CASE (INV_LENGHT)

              x = mesh % rr(1, T_j_m_d(e) % vec)
              y = mesh % rr(2, T_j_m_d(e) % vec)

              ! Centroid coordinates                                        
              x_g = SUM(x) / (SIZE(x) * 1.d0)
              y_g = SUM(y) / (SIZE(y) * 1.d0)

              ! Middle points coordinates                                   
              x_m12 = (x(1) + x(2)) / 2.d0;   y_m12 = (y(1) + y(2)) / 2.d0  
              x_m23 = (x(2) + x(3)) / 2.d0;   y_m23 = (y(2) + y(3)) / 2.d0  
              x_m31 = (x(3) + x(1)) / 2.d0;   y_m31 = (y(3) + y(1)) / 2.d0  

              ! Distances of centroid from triangles edge middle points     
              G_M12 = SQRT((x_g - x_m12)**2  +  (y_g - y_m12)**2)           
              G_M23 = SQRT((x_g - x_m23)**2  +  (y_g - y_m23)**2)           
              G_M31 = SQRT((x_g - x_m31)**2  +  (y_g - y_m31)**2)           

              l = (/G_M12, G_M23, G_M31/)                                   

              Jou_Mod = 1 / MINVAL(l)**(params % jmf_alpha)                 


         END SELECT

         Poi_Coe = params % Poi_Coe

         ! Stiffness Matrix for element 'e'                                     
         CALL Mat_Tria (x, y, Poi_Coe, Jou_Mod,   KK_e)                         
                                                                                  
         DO i = 1, 6                                                              
                                                                                  
            r((i-1)*6+1:(i-1)*6+6) = n(i)                                       
            h((i-1)*6+1:(i-1)*6+6) = KK_e(i,:)                                      
                                                                                  
            DO j = 1, 6                                                         
               c((i-1)*6+j) = n(j)                                                
            ENDDO                                                                 
                                                                                  
         ENDDO                                                                  


         DO idx_ = 1, SIZE(r)
             
           row_i = r(idx_)
           col_i = c(idx_)
           ke_i  = h(idx_)

           idx = csridx_srt(row_i, col_i, KK%i, KK%j)
             
           KK%e(idx) = KK%e(idx) + ke_i
             
         ENDDO
       
       CLOSE (idf)
    
     ENDDO

     END FUNCTION  Stiffness_Matrix





     SUBROUTINE  set_CSR_col (k_d, T_j_m_d, T_m_j_d, KK)
     !--------------------------------------------------------------------------------
     IMPLICIT NONE

     INTEGER,                   INTENT(IN)    :: k_d 
     TYPE(D_I_V), DIMENSION(:), INTENT(IN)    :: T_j_m_d, T_m_j_d
     TYPE(CSR_matrix),          INTENT(INOUT) :: KK

     INTEGER,      DIMENSION(:),    ALLOCATABLE :: n
     INTEGER,      DIMENSION(:),    ALLOCATABLE :: row, col
  
     INTEGER  ::  m, i, idx_, n_j, row_i, col_i, k  
     !--------------------------------------------------------------------------------
  
     KK % j = k_d * SIZE(T_m_j_d)

     DO m = 1, SIZE(T_j_m_d)
              
        n_j = SIZE(T_j_m_d(m)%vec)   ! number of nodes
        
        ALLOCATE( n (k_d*n_j),                      &
                  row ( (k_d*n_j)**2 ),             &
                  col ( (k_d*n_j)**2 ) )
        
        row = 0; col = 0
  
        DO i = 1, 3
           DO k = 1, k_d
              n((i-1)*k_d + k) = ( (T_j_m_d(m)%vec(i))-1 )*k_d + k
           ENDDO        
        ENDDO
  
        DO i = 1, k_d*n_j
           
           row ((i-1) * k_d * n_j + 1:i * k_d * n_j) = n(i)
           col ((i-1) * k_d*n_j + 1:i * k_d*n_j ) = n
           
        ENDDO          

        DO idx_ = 1, SIZE(row)
           
           row_i = row (idx_)
           col_i = col (idx_)

           CALL setja(row_i, col_i, KK % i, KK % j)
        
        ENDDO

        DEALLOCATE (n, row, col)
              
     ENDDO   

     END SUBROUTINE  set_CSR_col




     SUBROUTINE  Mat_Tria (x, y, Poi_Coe, Jou_Mod,  KKe)
     !----------------------------------------------------------------- 
     IMPLICIT NONE                                                              
                                                                                
     REAL(KIND=8), DIMENSION(3), INTENT(IN) :: x, y                             
     REAL(KIND=8),               INTENT(IN) :: Poi_Coe, &                       
                                               Jou_Mod                          

     REAL(KIND=8), DIMENSION(6,6) :: KKe                                        

     REAL(KIND=8) :: kn, knu, ks,   &                                           
                     x21, x32, x31, &                                           
                     y21, y32, y31, &                                           
                     A123
     !------------------------------------------------------------------

                                                                                
     x21 = x(2) - x(1)                                                        
     x32 = x(3) - x(2)                                                        
     x31 = x(3) - x(1)                                                        
                                                                              
     y21 = y(2) - y(1)                                                        
     y32 = y(3) - y(2)                                                        
     y31 = y(3) - y(1)                                                        
                                                                              
     A123 = -0.5 * ( x32*y21 - x21*y32 )                                      
                                                                              
     kn  = 0.25* Jou_Mod / (A123 * (1 - Poi_Coe*Poi_Coe))                     
     knu =  kn * Poi_Coe                                                      
     ks  = 0.25* Jou_Mod / (A123 * (1 + Poi_Coe))                             
                                                                              
     KKe(1,1) =   kn*y32*y32  + ks*x32*x32                                    
     KKe(2,1) = - knu*y32*x32 - ks*x32*y32                                    
     KKe(3,1) = - kn*y32*y31  - ks*x32*x31                                    
     KKe(4,1) =   knu*y32*x31 + ks*x32*y31                                    
     KKe(5,1) =   kn*y32*y21  + ks*x32*x21                                    
     KKe(6,1) = - knu*y32*x21 - ks*x32*y21                                    
                                                                              
     KKe(2,2) =   kn*x32*x32  + ks*y32*y32                                    
     KKe(3,2) =   knu*x32*y31 + ks*y32*x31                                    
     KKe(4,2) = - kn*x32*x31  - ks*y32*y31                                    
     KKe(5,2) = - knu*x32*y21 - ks*y32*x21                                    
     KKe(6,2) =   kn*x32*x21  + ks*y32*y21                                    
                                                                              
     KKe(3,3) =   kn*y31*y31  + ks*x31*x31                                    
     KKe(4,3) = - knu*y31*x31 - ks*x31*y31                                    
     KKe(5,3) = - kn*y31*y21  - ks*x31*x21                                    
     KKe(6,3) =   knu*y31*x21 + ks*x31*y21                                    
                                                                              
     KKe(4,4) =   kn*x31*x31  + ks*y31*y31                                    
     KKe(5,4) =   knu*x31*y21 + ks*y31*x21                                    
     KKe(6,4) = - kn*x31*x21  - ks*y31*y21                                    
                                                                              
     KKe(5,5) =   kn*y21*y21  + ks*x21*x21                                    
     KKe(6,5) = - knu*y21*x21 - ks*x21*y21                                    
                                                                              
     KKe(6,6) =   kn*x21*x21  + ks*y21*y21                                    
                                                                              
                                                                              
     KKe(1,2) =  KKe(2,1)                                                     
     KKe(1,3) =  KKe(3,1)                                                     
     KKe(1,4) =  KKe(4,1)                                                     
     KKe(1,5) =  KKe(5,1)                                                     
     KKe(1,6) =  KKe(6,1)                                                     
                                                                              
     KKe(2,3) =  KKe(3,2)                                                     
     KKe(2,4) =  KKe(4,2)                                                     
     KKe(2,5) =  KKe(5,2)                                                     
     KKe(2,6) =  KKe(6,2)                                                     

     KKe(3,4) =  KKe(4,3)                                                     
     KKe(3,5) =  KKe(5,3)                                                     
     KKe(3,6) =  KKe(6,3)                                                     
                                                                              
     KKe(4,5) =  KKe(5,4)                                                     
     KKe(4,6) =  KKe(6,4)                                                     

     KKe(5,6) =  KKe(6,5)                                                     

     END SUBROUTINE  Mat_Tria                                                  




     SUBROUTINE  Dirichlet_BoCo (mesh, displ, AA, B)

     ! Enforces Dirichlet boundary conditions, The vector 'displ'
     ! contains the  value  of the nodes displacements to enforce
     ! at boundaries.
     !-----------------------------------------------------------------
     IMPLICIT NONE

     TYPE(grid_type),                INTENT(IN)    :: mesh
     TYPE(D_I_V_R),    DIMENSION(:), INTENT(IN)    :: displ
     TYPE(CSR_Matrix),               INTENT(INOUT) :: AA
     REAL(KIND=8),     DIMENSION(:), INTENT(INOUT) :: B


     INTEGER, DIMENSION(:), ALLOCATABLE :: js
     INTEGER :: i, n, p, k
     !-----------------------------------------------------------------

     ! js contains the indices of the AA rows affected by 
     ! boundary conditions. Only nodes at boundaries are
     ! affected by BC, so vector js has a dimension equal
     ! to twice the number of nodes at boundaries. Twice
     ! because for each node you have x-displacement and
     ! y-displacement. Since the variable  mesh % Nj_b 
     ! already contains the double of the number of nodes
     ! for the fact that the mesh at boundarie is not 
     ! connected the allocation is made with mesh % Nj_b.
     
     ALLOCATE (js(mesh % Nj_b))

     k = 0

     DO n = 1, mesh % Nj_d

       ! Check if the node considered (by the loop over
       ! the domain nodes indices) belongs to the boundary
       IF (ANY(mesh % jd_jb == n)) THEN

         js(k + 1) = 2*(n-1) + 1
         js(k + 2) = 2*(n-1) + 2

         ! Enforces the RHS to the value specifed by
         ! dynamic vector 'displ'
         B(2*(n-1) + 1) = displ(n) % vec(1)
         B(2*(n-1) + 2) = displ(n) % vec(2)

         k = k + 2 

       ENDIF

     ENDDO

  
     DO n = 1, SIZE(js)

        i = js(n)

        ! Only for the rows interested by B.C. the following
        ! loop places 1 in the diagonal element and 0 in the
        ! extra-diagonal elements

        DO p = AA % i(i), AA % i(i+1) - 1
          
           IF (AA % j(p) == i) THEN
              AA % e(p) = 1
           ELSE
              AA % e(p) = 0
           ENDIF
       
        ENDDO

     ENDDO

     IF (ALLOCATED(js))  DEALLOCATE (js)

     END SUBROUTINE  Dirichlet_BoCo

   END FUNCTION  Reconstruct_Mesh




   SUBROUTINE Build_Refine_vec(grid, refinement_layout, nds_to_add, bnds_to_add, np_refine, bnp_refine)
   !--------------------------------------------------------------------------------------------------------
   IMPLICIT NONE
     
   TYPE(grid_type), INTENT(IN)    :: grid
   TYPE(E_R_P),     INTENT(IN)    :: refinement_layout     

   LOGICAL, DIMENSION(:), INTENT(INOUT) :: np_refine, bnp_refine 
   INTEGER,               INTENT(INOUT) :: nds_to_add, bnds_to_add
  
  
   INTEGER :: i, node_pair_index, nd1,  nd2  
   !--------------------------------------------------------------------------------------------------------
      
   nds_to_add  = 0
   bnds_to_add = 0

   DO i = 1, refinement_layout % N_refined_edges  
  
     nd1 = refinement_layout % edge_nodes(1,i)                
     nd2 = refinement_layout % edge_nodes(2,i)
!WRITE(*,*) 'reconstruction_proc start Extract_np_index'
     node_pair_index = Extract_np_index(grid % c_j, nd1, nd2)
!WRITE(*,*) 'reconstruction_proc end Extract_np_index'

     np_refine(node_pair_index) = .TRUE.
     nds_to_add = nds_to_add + 1

     IF (grid % cb_cd(node_pair_index) .NE. 0) THEN
    
       bnp_refine(grid % cb_cd(node_pair_index)) = .TRUE.
       bnds_to_add = bnds_to_add + 1
    
     ENDIF    
      
   ENDDO
     
   END SUBROUTINE  Build_Refine_vec


END MODULE  reconstruction_proc
