MODULE rh_setting
 
 USE structures 
 USE error_estimator
 USE marking_proc

 CONTAINS


  SUBROUTINE  Set_Adaption_History(grid, grid_name, solution, params, refinement_history)
  !-------------------------------------------------------------------------------
  IMPLICIT NONE                                         

  TYPE(grid_type),      INTENT(IN) :: grid
  CHARACTER(LEN=30) :: grid_name
  TYPE(solution_type),  INTENT(IN) :: solution


  LOGICAL, DIMENSION(grid%Nc_d) :: refined_edges_flag
  LOGICAL, DIMENSION(grid%Nc_d) :: coarsen_edges_flag
  LOGICAL, DIMENSION(:), POINTER :: valid_nodes_flag
  INTEGER :: entry_level
 
  TYPE(adaption_param)                         :: params  
  TYPE(E_R_P), DIMENSION(0:params%entry_level) :: refinement_history  
  LOGICAL :: cnst_exist
  INTEGER :: Ncc, i
  LOGICAL, DIMENSION(:), ALLOCATABLE :: constraints
  !-------------------------------------------------------------------------------

  entry_level = params % entry_Level        
  
  ! Building refinement_history for new adaptation_level.
  CALL Compute_Error_Estimator( grid, grid_name, solution, params % estimator, &
    params % ref_type )
  
  ! The two following subroutines produce two lists of logical:
  ! - refined_edges_flag
  ! - coarsen_edges_flag.
  ! The  first  one  specify  if  an  edge has to be refined, the second one, 
  ! instead specify if an edge has an error level smaller than the thresold.

  refined_edges_flag = Mark_Edges_for_Refinement(grid, params%ref_type, params%box,  &
                                                 params%estimator, params%min_angle, &
                                                 params%max_number_QL, params%min_length)
  
  INQUIRE ( file='constraints.'//trim(grid_name), exist=cnst_exist )
  IF ( cnst_exist ) THEN
    open ( unit=123, file='constraints.'//trim(grid_name) )
    read (123,*) Ncc
    allocate ( constraints(Ncc) )
    do i = 1, Ncc
      read (123,*) constraints(i)
    end do
    close (123)
    
    refined_edges_flag = refined_edges_flag .and. (.not. constraints)
    
  END IF
  
  IF (entry_Level .GT. 0) THEN

   IF (ANY(grid % ele_type_d == 3)) THEN
     coarsen_edges_flag = .FALSE.
   ELSE 
     coarsen_edges_flag = Mark_Edges_for_Derefinement(grid, params % estimator)
   ENDIF 
  
  ENDIF    


  ! This subroutine specifies how to refine grid from level 'entry_Level' to level 
  !'entry_Level+1' according to edges to be refined. 

  CALL Extend_Refinement_History(grid, refined_edges_flag, refinement_history(entry_Level))
 
   
  ! Nodes are marked for deletion according to three main tests:
  ! - Error level over each edge connected to a node (coarsen_edges_flag)
  ! - Refinement status of the connected edges (refined_edges_flag)
  ! - Grid validity considerations (refinement_history)
  ! The output is a list of logical (valid_nodes_flag) that specifies if a node
  ! has to be removed (.FALSE.) or not (.TRUE.).  

  ALLOCATE( valid_nodes_flag(grid%Nj_d + COUNT(refined_edges_flag)) )
  
  IF (entry_Level == 0) THEN
  
    valid_nodes_flag = .TRUE.
    
  ELSE
   
   IF (ANY(grid % ele_type_d == 3)) THEN
   
    valid_nodes_flag = .TRUE.
   
   ELSE 
   
    valid_nodes_flag = Mark_Nodes_for_Deletion (grid, refined_edges_flag, coarsen_edges_flag, &
                                                refinement_history)
   ENDIF
                                                
  ENDIF
  

  ! This subroutine specifies how to refine grids from level '0' to level 
  !'A_Level' according to nodes to be deleted. 
  ! ( EDIT refinement_history(0:A_Level) )  

  IF (entry_Level > 0) CALL Update_Refinement_History(grid, valid_nodes_flag, refinement_history)


  IF (ASSOCIATED(valid_nodes_flag))  DEALLOCATE(valid_nodes_flag)
 
 END SUBROUTINE Set_Adaption_History
 
 
 




 SUBROUTINE  Extend_Refinement_History(grid, refined_edges_flag, &
                                       edge_refinement_pattern)
                                               
 !==========================================================================!
  IMPLICIT NONE                                         

  TYPE(grid_type),       INTENT(IN)  :: grid
  LOGICAL, DIMENSION(:), INTENT(IN)  :: refined_edges_flag


  INTEGER :: Nc_refined, idx, c


  TYPE(E_R_P), INTENT(INOUT) :: edge_refinement_pattern
 !==========================================================================!

   Nc_refined = COUNT(refined_edges_flag)
   edge_refinement_pattern%N_refined_edges = Nc_refined
   
   ALLOCATE(edge_refinement_pattern%inserted_nodes(  Nc_refined), &
            edge_refinement_pattern%edge_nodes    (2,Nc_refined))

   idx = 1 
   DO c = 1, grid % Nc_d, +1
   
    IF (refined_edges_flag(c)) THEN 

      edge_refinement_pattern%inserted_nodes(idx) = grid%Nj_d + idx
      edge_refinement_pattern%edge_nodes(1,idx) = MINVAL(grid%j_c(:,c))
      edge_refinement_pattern%edge_nodes(2,idx) = MAXVAL(grid%j_c(:,c))
      
      idx = idx + 1

    ENDIF 
      
   ENDDO
        
   ALLOCATE(edge_refinement_pattern%adjacent_nodes(1))


 END SUBROUTINE Extend_Refinement_History

   





 SUBROUTINE Update_Refinement_History( grid, valid_nodes_flag, &
                                       edge_refinement_pattern )

 !==========================================================================!
  IMPLICIT NONE
  
  TYPE(grid_type),           INTENT(IN)    :: grid
  LOGICAL,     DIMENSION(:), INTENT(IN)    :: valid_nodes_flag
  TYPE(E_R_P), DIMENSION(0:grid%adaptation_level), INTENT(INOUT) :: edge_refinement_pattern
  
  INTEGER                              :: i, j, k
  INTEGER                              :: min_,      max_  
  INTEGER                              :: min_index, max_index
  INTEGER                              :: new_inserted_nodes
  INTEGER, DIMENSION(:),   ALLOCATABLE :: Inserted_nodes
  INTEGER, DIMENSION(:),   ALLOCATABLE :: nodes_mask
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Edge_nodes
 !==========================================================================!
    
  ALLOCATE (nodes_mask(SIZE(valid_nodes_flag)))
  nodes_mask = Build_Nodes_Mask(valid_nodes_flag)
  
  DO i = 0, grid % adaptation_Level
   
   ! Copying inserted_nodes and edge_nodes
   IF (ALLOCATED(Inserted_nodes)) DEALLOCATE(Inserted_nodes)
   ALLOCATE (Inserted_nodes(edge_refinement_pattern(i) % N_refined_edges))
   Inserted_nodes = edge_refinement_pattern(i) % inserted_nodes
   
   IF (ALLOCATED(Edge_nodes)) DEALLOCATE(Edge_nodes)   
   ALLOCATE (Edge_nodes(2,edge_refinement_pattern(i) % N_refined_edges))
   Edge_nodes = edge_refinement_pattern(i) % edge_nodes   
   
   
   min_ = 1
   max_ = edge_refinement_pattern(i) % N_refined_edges
   
   IF (max_ > 0) THEN
   
     min_index = edge_refinement_pattern(i) % inserted_nodes(min_)
     max_index = edge_refinement_pattern(i) % inserted_nodes(max_)
   
     new_inserted_nodes = COUNT( valid_nodes_flag(min_index:max_index) )
     edge_refinement_pattern(i) % N_refined_edges = new_inserted_nodes
   
     DEALLOCATE (edge_refinement_pattern(i) % inserted_nodes)
     DEALLOCATE (edge_refinement_pattern(i) % edge_nodes)
   
     ALLOCATE (edge_refinement_pattern(i) % inserted_nodes(new_inserted_nodes))
     ALLOCATE (edge_refinement_pattern(i) % edge_nodes(2,new_inserted_nodes))
   
   ENDIF
   
   k = 1
   DO j = 1, SIZE(Inserted_nodes)
   
    IF (valid_nodes_flag(Inserted_nodes(j))) THEN
    
      edge_refinement_pattern(i) % inserted_nodes(k) = nodes_mask(Inserted_nodes(j))
      edge_refinement_pattern(i) % edge_nodes(1,k) = nodes_mask(Edge_nodes(1,j))
      edge_refinement_pattern(i) % edge_nodes(2,k) = nodes_mask(Edge_nodes(2,j))
      
      k = k + 1
    
    ENDIF
   
   ENDDO
   
   DEALLOCATE (edge_refinement_pattern(i) % adjacent_nodes)
      
  ENDDO 


  IF (ALLOCATED(Inserted_nodes))  DEALLOCATE (Inserted_nodes)
  IF (ALLOCATED(nodes_mask))      DEALLOCATE (nodes_mask)
  IF (ALLOCATED(Edge_nodes))      DEALLOCATE (Edge_nodes)
  
 END SUBROUTINE Update_Refinement_History
 



   
 FUNCTION  Build_Nodes_Mask(valid_nodes_flag) RESULT(nodes_mask)
 
 !===================================================================================!
  LOGICAL, DIMENSION(:), INTENT(IN)  :: valid_nodes_flag  
  INTEGER, DIMENSION(:), POINTER     :: nodes_mask  
 
  INTEGER                            :: i, j
 !===================================================================================!     
  
  ALLOCATE( nodes_mask(SIZE(valid_nodes_flag)) )
  nodes_mask = 0
  
  j = 1
  DO i = 1, SIZE(valid_nodes_flag)
          
    IF ( valid_nodes_flag(i) ) THEN  
      nodes_mask(i) = j 
      j = j + 1          
    ENDIF    
          
  ENDDO   
   
 END FUNCTION Build_Nodes_Mask


END MODULE rh_setting
