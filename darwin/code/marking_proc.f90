MODULE marking_proc
 
 USE dynamic_vector
 USE error_estimator
 USE grid_utils
 USE lib_geo
 USE select_strategy
 USE structures

 CONTAINS

   FUNCTION Mark_Edges_for_Refinement(grid, ref_type, boxes, estimator, &
                                      minimum_angle, maxN_loop,         &
                                      minimum_length)   RESULT (refined_edges_flag)                           
   !----------------------------------------------------------------------------------------
   IMPLICIT NONE
 
   TYPE(grid_type),             INTENT(IN) :: grid    
   INTEGER,                     INTENT(IN) :: ref_type
   TYPE(D_I_V_R), DIMENSION(:), INTENT(IN) :: boxes
   REAL(KIND=8),                INTENT(IN) :: minimum_angle
   INTEGER,                     INTENT(IN) :: maxN_loop
   REAL(KIND=8),                INTENT(IN) :: minimum_length

   LOGICAL, DIMENSION(grid%Nc_d) :: refined_edges_flag

   LOGICAL, DIMENSION(grid%Nc_d) :: c__refined_edges_flag

   TYPE(estimator_type), DIMENSION(:) :: estimator

   REAL(KIND=8)                            :: x1, x2, y1, y2
 
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: limit_R
 
   REAL(KIND=8), DIMENSION(grid%Nc_d) :: length 
   REAL(KIND=8), DIMENSION(grid%Nc_d) :: length_copy

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x, y
   INTEGER :: l, m, n, nd1_in, nd2_in

   LOGICAL, DIMENSION(grid%Nc_d) :: skip_edge 
   LOGICAL, DIMENSION(grid%Nc_d) :: skip_edge_multivariables 
 
   INTEGER, DIMENSION(grid%Nc_d) :: edges_copy
   INTEGER, DIMENSION(:), ALLOCATABLE :: active_edges

   INTEGER :: pass, est, edge_skipped_counter, &
              k, nd1, nd2, i, g, j
   !----------------------------------------------------------------------------------------

   PRINT*, '   Averaging error...'

   refined_edges_flag = .FALSE.
   skip_edge          = .FALSE.
   skip_edge_multivariables = .FALSE.
   edge_skipped_counter = 0


   ! Computing length for minimum length criterion 
   DO k = 1, grid % Nc_d
 
     nd1 = grid % j_c(1, k);   nd2 = grid % j_c(2, k)

     x1 = grid % rr(1, nd1);   y1 = grid % rr(2, nd1)
     x2 = grid % rr(1, nd2);   y2 = grid % rr(2, nd2)

     length(k) = SQRT((y2-y1)**2 + (x2-x1)**2)
      
   ENDDO
 
   length_copy = length

   CALL Write_Edges_Length(length_copy) 


   ! Select refinement type. 
   !
   !  0 = whole domain
   !  1 = boxed
 
   IF (ref_type == 0) THEN

     ALLOCATE (active_edges(grid % Nc_d))
     active_edges = (/ (i, i = 1,grid % Nc_d) /)

   ELSE

     edges_copy = 0
     j = 1

     DO i = 1, grid % Nc_d

       nd1 = grid%j_c(1, i);   nd2 = grid%j_c(2, i)

       x1 = grid%rr(1, nd1);   y1 = grid%rr(2, nd1)
       x2 = grid%rr(1, nd2);   y2 = grid%rr(2, nd2)

       DO k = 1, SIZE(boxes)

         ALLOCATE (x(SIZE(boxes(k) % vec)/2), &
                   y(SIZE(boxes(k) % vec)/2)) 

         x = (/ boxes(k) % vec(1:SIZE(boxes(k) % vec):2) /)
         y = (/ boxes(k) % vec(2:SIZE(boxes(k) % vec):2) /)

         n = SIZE(boxes(k) % vec) / 2

         CALL locpt (x1, y1, x, y, n, l, m)
         nd1_in = l

         CALL locpt (x2, y2, x, y, n, l, m)
         nd2_in = l

         IF (nd1_in /= -1  .OR.  nd2_in /= -1) THEN

           edges_copy(j) = i
           j = j + 1

         ENDIF 

         DEALLOCATE (x, y)

       ENDDO

     ENDDO 

     ALLOCATE (active_edges(COUNT(edges_copy /= 0)))
     active_edges = edges_copy(1:SIZE(active_edges))

   ENDIF



   DO est = 1, SIZE(estimator), +1
    
    ALLOCATE (limit_R(SIZE(estimator(est) % errors, 1)))
    ALLOCATE (estimator(est) % R_Thres(     estimator(est) % passages, &
                                       SIZE(estimator(est) % errors, 1)))

    ! Loop for multiple passages strategy
    DO pass = 1, estimator(est)%passages   

      IF ( pass .eq. 1 ) THEN
        CALL Compute_Limit('R', estimator(est), skip_edge_multivariables,  limit_R)
        estimator(est)%R_Thres(pass,:) = limit_R
      ELSE
        CALL Compute_Limit('R', estimator(est), skip_edge,  limit_R)
        estimator(est)%R_Thres(pass,:) = limit_R
      END IF


      ! Loop over active edges
      DO g = 1, SIZE(active_edges)             

        k = active_edges(g)

        IF (.NOT. skip_edge(k)) THEN 

          SELECT CASE ( estimator(est)%estimator_function )

           CASE (ERR_ESTIMATOR_GRADIENT)  

            IF (estimator(est)%errors(1,k) .GT. limit_R(1)) THEN

              IF (length(k)*0.5 .LT. minimum_length) THEN
                skip_edge(k) = .TRUE.
                edge_skipped_counter = edge_skipped_counter + 1          
              ELSE    
                refined_edges_flag(k) = .TRUE.
                skip_edge(k)          = .TRUE.
              ENDIF
 
            ENDIF    



           CASE (ERR_ESTIMATOR_ANISOTROPIC)  

            IF( estimator(est)%errors(1,k) .GT. limit_R(1) ) THEN

              IF ( length(k)*0.5 .LT. minimum_length ) THEN
                skip_edge(k) = .TRUE.
                edge_skipped_counter = edge_skipped_counter + 1                          
              ELSE    
                refined_edges_flag(k) = .TRUE.
                skip_edge(k)          = .TRUE.
              ENDIF
              
            ENDIF                


           CASE (ERR_ESTIMATOR_METRIC_BASED)  

            IF( estimator(est)%errors(1,k) .GT. limit_R(1) ) THEN

              IF ( length(k)*0.5 .LT. minimum_length ) THEN
                skip_edge(k) = .TRUE.
                edge_skipped_counter = edge_skipped_counter + 1                          
              ELSE    
                refined_edges_flag(k) = .TRUE.
                skip_edge(k)          = .TRUE.
              ENDIF

            ENDIF                


           CASE DEFAULT

            IF( ( estimator(est)%errors(1,k) .GT. limit_R(1) )    .OR. & 
                ( estimator(est)%errors(2,k) .GT. limit_R(2) ) )  THEN 

              IF ( length(k)*0.5 .LT. minimum_length ) THEN
                skip_edge(k) = .TRUE.
                edge_skipped_counter = edge_skipped_counter + 1                          
              ELSE    
                refined_edges_flag(k) = .TRUE.
                skip_edge(k)          = .TRUE.
              ENDIF

            ENDIF

          END SELECT

        ENDIF
      ENDDO     
     ENDDO 
     
     DEALLOCATE (limit_R)
            
    ENDDO

    DEALLOCATE (active_edges)


    c__refined_edges_flag = refined_edges_flag

    CALL Mark_for_Quality(grid, minimum_angle, maxN_loop, refined_edges_flag)

    write (*,'(4x,a39,1x,i8)') '- edges marked for quality enforcement:', COUNT(refined_edges_flag) - &
                                                        COUNT(c__refined_edges_flag)

   END FUNCTION Mark_Edges_for_Refinement





   FUNCTION  Mark_Edges_for_Derefinement(grid, estimator)   RESULT(coarsen_edges_flag) 
   !-----------------------------------------------------------------------------------------------
   IMPLICIT NONE
 
   TYPE(grid_type), INTENT(IN)  :: grid
 
   INTEGER :: k  
   LOGICAL,              DIMENSION(:), POINTER     :: coarsen_edges_flag
   LOGICAL,              DIMENSION(:), ALLOCATABLE :: skip_edge 
   INTEGER                                         :: pass, est
   REAL(KIND=8),         DIMENSION(:), ALLOCATABLE :: limit_D
   TYPE(estimator_type), DIMENSION(:)              :: estimator
   !-----------------------------------------------------------------------------------------------
 
   ALLOCATE (coarsen_edges_flag(grid%Nc_d), &
                      skip_edge(grid%Nc_d))
    
   coarsen_edges_flag = .TRUE.
   skip_edge          = .FALSE.
 

   DO est = 1, SIZE(estimator)
    
     ALLOCATE (limit_D(SIZE(estimator(est) % errors, 1)))
     ALLOCATE (estimator(est)% C_Thres(     estimator(est) % passages, &
                                       SIZE(estimator(est) % errors, 1)))


     ! Loop for multiple passages strategy
     DO pass = 1, estimator(est)%passages   

      CALL Compute_Limit('D', estimator(est), skip_edge,  limit_D)

      estimator(est)%C_Thres(pass,:) = limit_D


      ! Loop over edges
      DO k = 1, grid % Nc_d

        IF (.NOT. skip_edge(k)) THEN


          SELECT CASE ( estimator(est)%estimator_function )

          ! coarsen_edges_flag(k) = TRUE if the error associated to edge 'k' is 
          !                         LOWER than the thresold
          ! coarsen_edges_flag(k) = FALSE if the error associated to edge 'k' is 
          !                         GREATER than the thresold

           CASE (ERR_ESTIMATOR_GRADIENT)  

            IF( estimator(est)%errors(1,k) .GT. limit_D(1) )  THEN 

                coarsen_edges_flag(k) = .FALSE.
                skip_edge(k)          = .TRUE.

            ENDIF    



           CASE (ERR_ESTIMATOR_ANISOTROPIC)  

            IF( estimator(est)%errors(1,k) .GT. limit_D(1) ) THEN

                coarsen_edges_flag(k) = .FALSE.
                skip_edge(k)          = .TRUE.

            ENDIF                


           CASE (ERR_ESTIMATOR_METRIC_BASED)  

            IF( estimator(est)%errors(1,k) .GT. limit_D(1) ) THEN

                coarsen_edges_flag(k) = .FALSE.
                skip_edge(k)          = .TRUE.

            ENDIF                


           CASE DEFAULT

            IF( ( estimator(est)%errors(1,k) .GT. limit_D(1) )   .OR.  & 
                ( estimator(est)%errors(2,k) .GT. limit_D(2) ) )  THEN

                coarsen_edges_flag(k) = .FALSE.      
                skip_edge(k)          = .TRUE.     

            ENDIF

          END SELECT

        ENDIF
      ENDDO     
     ENDDO
     
     DEALLOCATE (limit_D)          
    
    ENDDO

    DEALLOCATE (skip_edge)
 
   END FUNCTION Mark_Edges_for_Derefinement





   SUBROUTINE Compute_Limit(action, estimator, skip,   limit)
   !----------------------------------------------------------------------------
   IMPLICIT NONE
 
   CHARACTER(LEN=1),                   INTENT(IN) :: action 
   TYPE(estimator_type),               INTENT(IN) :: estimator
   LOGICAL,              DIMENSION(:), INTENT(IN) :: skip

   REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: limit
 
   INTEGER                                   :: i, j, dim
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: valid_edges
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: mean_error, deviation
   !----------------------------------------------------------------------------
 
    dim = COUNT(.NOT. skip)

    IF (ALLOCATED(valid_edges))   DEALLOCATE (valid_edges)
    ALLOCATE (valid_edges(SIZE(estimator % errors,1), dim))


    ! Selecting edges for limit computation only if edges haven't been
    ! considered in a previous passage
    i = 1

    DO j = 1, SIZE(skip)

      IF (.NOT. skip(j)) THEN

        valid_edges(:,i) = estimator % errors(:,j)   

        i = i + 1

      ENDIF 

    ENDDO   

 
    ! Computing Limit
    ALLOCATE (mean_error(SIZE(estimator%errors, 1)),  &
               deviation(SIZE(estimator%errors, 1)))
 
    IF (SIZE(valid_edges, 1) .EQ. 1) THEN   
      
      mean_error(1) = Average(valid_edges(1,:))

      deviation(1) = 0.d0

      DO j = 1, SIZE(valid_edges,2) 
        deviation(1) = deviation(1) + (valid_edges(1,j) - mean_error(1))**2
      ENDDO

      deviation(1) = SQRT(deviation(1) / (1.d0*dim))

      IF (action .EQ. 'R') limit(1) = mean_error(1) + deviation(1)
      IF (action .EQ. 'D') limit(1) = mean_error(1)       

      
    ELSEIF (SIZE(valid_edges, 1) .EQ. 2) THEN  

      DO i = 1, 2

        mean_error(i) = Average(valid_edges(i,:))

        deviation(i) = 0.d0
 
        DO j = 1, SIZE(valid_edges,2) 
          deviation(i) = deviation(i) + (valid_edges(i,j) - mean_error(i))**2
        ENDDO

        deviation(i) = SQRT(deviation(i) / (1.d0*dim))

        IF (action .EQ. 'R') limit(i) = mean_error(i) + deviation(i)
        IF (action .EQ. 'D') limit(i) = mean_error(i)       
      
      ENDDO      

    ENDIF
    
    DEALLOCATE (valid_edges, mean_error, deviation)

   END SUBROUTINE Compute_Limit
 




   FUNCTION  Mark_Nodes_for_Deletion(grid, refined_edges_flag, coarsen_edges_flag, &
                                     refinement_history)   RESULT (valid_nodes_flag)
   !===================================================================================!
   IMPLICIT NONE
 
   TYPE(grid_type),                               INTENT(IN) :: grid
   TYPE(E_R_P), DIMENSION(0:grid%adaptation_level), INTENT(IN) :: refinement_history
   LOGICAL,     DIMENSION(:),                     INTENT(IN) :: refined_edges_flag
   LOGICAL,     DIMENSION(:),                     INTENT(IN) :: coarsen_edges_flag
   
   LOGICAL,     DIMENSION(:), POINTER       :: valid_nodes_flag
   INTEGER                                  :: A_Level, edge_idx
   INTEGER                                  :: i, j, k, l, m
   INTEGER                                  :: central_node
   INTEGER                                  :: local_dim
   INTEGER                                  :: added_nodes
   INTEGER                                  :: new_nodes_dim
   INTEGER,     DIMENSION(:),   ALLOCATABLE :: adj_nodes
   INTEGER,     DIMENSION(2)                :: current_np
   INTEGER,     DIMENSION(2)                :: refined_np
   !===================================================================================!     
 
   A_Level = grid%adaptation_Level
 
   ! Number of nodes added to pass from current level to successive 
   added_nodes = COUNT(refined_edges_flag)
   new_nodes_dim = grid%Nj_d + added_nodes
    
   ALLOCATE( valid_nodes_flag(new_nodes_dim) )
   valid_nodes_flag = .FALSE.
 
   !  Error for connected edges:
   !  At the beginning  nodes  are  all  considered false. If while examining a node
   !  one  of  the  connected  edges  has an  error  greater  than the thresold (i.e. 
   !  coarsen_edge_flag = FALSE) or has been marked for refinement (angle condition),
   !  this  node is marked true and the loop goes to the following node in the list.
 
   DO i = 1, grid%Nj_d
    DO j = 1, SIZE( grid%c_j(i)%vec )  
      
      edge_idx = grid%c_j(i)%vec(j)    
      
      IF ( (.NOT. coarsen_edges_flag(edge_idx)) .OR. &  ! error condition
           (      refined_edges_flag(edge_idx)) ) THEN  ! error & angle condition 
           
        valid_nodes_flag(i) = .TRUE.      
        EXIT
        
      ENDIF
      
    ENDDO 
   ENDDO 

 
   !--Original grid nodes and added nodes are untouchable
   valid_nodes_flag(1:grid_zero % Nj_d) = .TRUE.

   valid_nodes_flag(grid % Nj_d+1:new_nodes_dim) = .TRUE.


   !--Verify, only for nodes marked for deletion, if a valid grid can be obtained. 
   !--In this case deletion is confirmed, otherwise deletion is aborted. 
 
   Levels:  DO i = 1, A_Level
     
     Nodes:  DO j = 1, SIZE(refinement_history(A_Level-i) % inserted_nodes)

      central_node = refinement_history(A_Level-i) % inserted_nodes(j)
      
      Skip_if_valid: IF (.NOT. valid_nodes_flag(central_node)) THEN
       
        local_dim = SIZE(refinement_history(A_Level-i) % adjacent_nodes(j) % vec)
        
        IF (ALLOCATED(adj_nodes)) DEALLOCATE(adj_nodes)
        ALLOCATE (adj_nodes(local_dim))

        adj_nodes = refinement_history(A_Level-i) % adjacent_nodes(j) % vec

        Ad_nodes: DO k = 1, SIZE(adj_nodes)

           current_np(1) = MIN(central_node, adj_nodes(k))
           current_np(2) = MAX(central_node, adj_nodes(k))     

      ! Loop to understand if the couple 'central_node,adj_nodes(k)' has ever been 
      ! refined. The loop is over each level in which the  couple  of nodes can be
      ! refined. For example a node inserted to pass from level 2 to 3  cannot  be
      ! part of couple of nodes refined passing  from  level 0 to 1 or from 1 to 2 
      ! or even from 2 to 3, it can only be part of a copule of nodes refined pas-
      ! sing from 3 to 4 or 4 to 5 and so on.
      
           DO l = A_Level, A_Level-i+1, -1
            DO m = 1, SIZE(refinement_history(l) % inserted_nodes)

             refined_np(1) = MINVAL(refinement_history(l) % edge_nodes(:,m)) 
             refined_np(2) = MAXVAL(refinement_history(l) % edge_nodes(:,m)) 

             IF (.NOT. ANY(current_np .NE. refined_np)) THEN

               valid_nodes_flag(central_node) = .TRUE.
               CYCLE Nodes

             ENDIF

            ENDDO
           ENDDO              
           
        ENDDO Ad_nodes
          
      ENDIF Skip_if_valid
     ENDDO Nodes 
              
   ENDDO Levels

   IF (ALLOCATED(adj_nodes)) DEALLOCATE(adj_nodes)

   END FUNCTION Mark_Nodes_for_Deletion
 
     
     
 
     
   SUBROUTINE  Mark_d_Elements(grid, np_flag, el_flag, pts_to_add, elem_to_add, patt)
   !=======================================================================================!     
   IMPLICIT NONE
   
   TYPE(grid_type),                             INTENT(IN)    :: grid
   LOGICAL,          DIMENSION(:),              INTENT(IN)    :: np_flag 
   LOGICAL,          DIMENSION(:),              INTENT(INOUT) :: el_flag
   TYPE(R_P),        DIMENSION(:),              INTENT(INOUT) :: patt
   INTEGER,                                     INTENT(OUT)   :: elem_to_add
     
   INTEGER                                                    :: k, l  
   INTEGER                                                    :: np_idx, k_vec_dim, pts_to_add, dum
   INTEGER,          DIMENSION(:), ALLOCATABLE                :: ref_np_for_elem  
   CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE                :: strategy        
   TYPE(D_I_V),      DIMENSION(:), ALLOCATABLE                :: np_idx_for_elem
   !=======================================================================================!
 
   ALLOCATE (ref_np_for_elem(grid % Nm_d))
   ALLOCATE (np_idx_for_elem(grid % Nm_d))
   ALLOCATE (strategy(grid % Nm_d))

   el_flag = .FALSE.
   elem_to_add = 0
   ref_np_for_elem = 0
   strategy = '   '

   DO k = 1, grid % Nm_d

     k_vec_dim = SIZE(grid%c_m(k) % vec)

     ALLOCATE (np_idx_for_elem(k) % vec(k_vec_dim))
     ALLOCATE (patt(k) % np_idx(k_vec_dim))

     np_idx_for_elem(k) % vec = 0

     DO l = 1, k_vec_dim

        np_idx = grid % c_m(k) % vec(l)

        IF (np_flag(np_idx)) THEN

          el_flag(k) = .TRUE.
          np_idx_for_elem(k) % vec(l) = np_idx
          ref_np_for_elem(k) = ref_np_for_elem(k) + 1

        ENDIF

     ENDDO

     CALL select_refinement(k, el_flag, grid%ele_type_d, ref_np_for_elem, elem_to_add, &
                            pts_to_add, dum, strategy, np_idx_for_elem(k)%vec)

     patt(k)%e_type     = grid%ele_type_d(k)
     patt(k)%ref_flag   = el_flag(k)
     patt(k)%ref_np_nmr = ref_np_for_elem(k)
     patt(k)%np_idx     = np_idx_for_elem(k)%vec
     patt(k)%strategy   = strategy(k)

   ENDDO

   DEALLOCATE (ref_np_for_elem, strategy, np_idx_for_elem)

   !--Result-------------------------------------------------------------------------------!
   !  elem_to_add                  ->  Number of elements selected for refinement
   !  el_flag(Ne_d)                ->  Flag (True or False) for refine elements
   !  ref_np_for_elem(Ne_d)        ->  Number of node pairs refined for each element
   !  np_idx_for_elem(Ne_d, Nc_d)  ->  (DIV) Index of node pairs refined for each 
   !                                   element:
   !                                    - if 0 = node pair not to be refined 
   !                                    - if 'idx' = node pair index to be refined
   !  strategy(Ne_d)               ->  Refinement strategy for each element
   !---------------------------------------------------------------------------------------!
 
   END SUBROUTINE Mark_d_Elements
 
 
 
 
 
 
 
   SUBROUTINE Mark_b_Elements(grid, bnp_flag, belem_flag, pts_to_add, belem_to_add,  &
                              b_patt, b_pts_to_add)

!   This routine  flag boundary elements to be refined observing which node pairs have to 
!   be refined  for  each element. Refined node pairs are identyfied starting from  nodes  
!   belonging to  the  element. It  determines  also  all  the  information required  for
!   refinement.

!  =======================================================================================!
    TYPE(grid_type),                             INTENT(IN)    :: grid
    LOGICAL,          DIMENSION(:),              INTENT(IN)    :: bnp_flag
    LOGICAL,          DIMENSION(:),              INTENT(INOUT) :: belem_flag
    INTEGER,                                     INTENT(INOUT) :: belem_to_add, &
                                                                  pts_to_add, b_pts_to_add
    TYPE(R_P),        DIMENSION(:),              INTENT(INOUT) :: b_patt
                                                                          
    INTEGER                                                    :: k, l, m, nmr, npmr
    INTEGER,          DIMENSION(2)                             :: couple
    INTEGER,          DIMENSION(:), ALLOCATABLE                :: loc_nodes, ref_bnp_for_el  
    CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE                :: b_strategy 
    TYPE(D_I_V),      DIMENSION(:), POINTER                    :: np_idx_for_elem
!  =======================================================================================!
 
   ALLOCATE ( ref_bnp_for_el(grid%Nm_b) )
   ALLOCATE ( np_idx_for_elem(grid%Nm_b) )
   ALLOCATE ( b_strategy(grid%Nm_b) ) 
 
   belem_flag     = .false.
   belem_to_add   = 0
   ref_bnp_for_el = 0
   b_strategy     = '   '
     
   DO k = 1, grid%Nm_b
     
     nmr  = N_nodes_ele(grid%ele_type_b(k))
     npmr = SIZE( ele_edges(grid%ele_type_b(k)),1 )
     
     ALLOCATE ( loc_nodes(nmr) )
     ALLOCATE ( np_idx_for_elem(k)%vec(npmr) )    
     ALLOCATE ( b_patt(k)%np_idx(npmr) )

     DO l = 1, nmr
      loc_nodes(l) = grid%j_m_b(k)%vec(l)  
     ENDDO
 
     ! Result:   loc_nodes(nmr) -> Indices of boundary nodes for 'k'-th element 
     DO l = 1, nmr
      
      couple(1) = loc_nodes(l)
      IF ( l+1 .GT. nmr ) THEN
        couple(2) = loc_nodes(1)
        GOTO 1 
      ENDIF         
      couple(2) = loc_nodes(l+1)

     ENDDO 

1    DO m = 1, npmr
      DO l = 1, grid%Nc_b
        IF ( (grid%j_c_b(1,l) == couple(1) .AND. grid%j_c_b(2,l) == couple(2)) .OR. &
             (grid%j_c_b(1,l) == couple(2) .AND. grid%j_c_b(2,l) == couple(1)) ) THEN
           
            np_idx_for_elem(k)%vec(m) = l
            IF ( bnp_flag(l) ) THEN
              belem_flag(k) = .true.
              ref_bnp_for_el(k) = ref_bnp_for_el(k) + 1
            ENDIF
            GOTO 2 
         ENDIF            
        ENDDO
2     ENDDO
     
   !--Result-------------------------------------------------------------------------------!
   !  np_idx_for_elem                     -> Boundary np indexes for 'k'-th element
   !---------------------------------------------------------------------------------------!      

        CALL select_refinement( k, belem_flag, grid%ele_type_b, ref_bnp_for_el, & 
                                belem_to_add, pts_to_add, b_pts_to_add,    &
                                b_strategy, np_idx_for_elem(k)%vec )
                                
     DEALLOCATE (loc_nodes)

     b_patt(k) % e_type     = grid % ele_type_b(k)
     b_patt(k) % ref_flag   = belem_flag(k)
     b_patt(k) % ref_np_nmr = ref_bnp_for_el(k)
     b_patt(k) % np_idx     = np_idx_for_elem(k) % vec
     b_patt(k) % strategy   = b_strategy(k)

   ENDDO
 
   !--Result-------------------------------------------------------------------------------!
   !  belem_to_add                 ->  Number of elements selected for refinement
   !  elem_flag(Ne_b)              ->  Flag (True or False) for refine elements
   !  ref_bnp_for_elem(Ne_b)       ->  Number of node pairs refined for each element
   !  np_idx_for_elem(Ne_d, Nc_d)  ->  (DIV) Index of node pairs refined for each 
   !                                   element:
   !                                    - if 0 = node pair not to be refined 
   !                                    - if 'idx' = node pair index to be refined
   !  b_strategy(Ne_b)             ->  Refinement strategy for each element
   !---------------------------------------------------------------------------------------!

   END SUBROUTINE Mark_b_Elements
  
   
END MODULE marking_proc
