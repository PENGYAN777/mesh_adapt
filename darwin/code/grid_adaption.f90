MODULE grid_adaption 

 USE io_proc
 USE reconstruction_proc
 USE rh_setting 
 USE structures   

 INTEGER, PUBLIC :: nodes_zero
  
 CONTAINS
   
      
 FUNCTION Adapt_Grid(grid, solution, params, add_solution)   RESULT(adapted_grid)

 ! Adaptation step from grid entry_Level to grid exit_level = entry_Level+1 
 ! based on solution on grid entry_Level.  Refinement history is updated. 
 !------------------------------------------------------------------------------------
 IMPLICIT NONE

 TYPE(grid_type),                    INTENT(IN)    :: grid    
 TYPE(solution_type),                INTENT(IN)    :: solution
 TYPE(adaption_param),               INTENT(IN)    :: params
 TYPE(solution_type),  DIMENSION(:), INTENT(INOUT) :: add_solution 
    
 INTEGER :: entry_Level 
 CHARACTER(LEN=30) :: grid_name  
 TYPE(E_R_P), DIMENSION(:), ALLOCATABLE :: adaption_history
 TYPE(solution_type) :: sol

 TYPE(grid_type) :: adapted_grid 
 !------------------------------------------------------------------------------------
  
  entry_level = params % entry_level
  grid_name   = params % grid_name
  
  if ( params % ref_type .ne. 2 ) then
  ! As a consequence of the adaption strategy we have the necessity to load
  ! the original grid, while the input grid is the mesh at entry_level state.
  ! While adapting from level 0 to 1 the two grids are the same. They are
  ! different for all the other cases.
  grid_zero = Load_Grid(0, grid_name, '0      ')  
  nodes_zero = grid_zero % Nj_d
  
  ALLOCATE ( adaption_history(0:entry_level) )
  
  ! If it is not the 0 to 1 step you have to load the previous adaption history.
  IF (entry_level > 0) adaption_history = Load_Adaption_History(entry_level)
  end if

  ! These are the two basic steps of the adaption:
  ! - Set up adaption history
  ! - Mesh adaption according to adaption history
  
  CALL Set_Adaption_History(grid, grid_name, solution, params, adaption_history)

  ALLOCATE (sol % ww(SIZE(solution % ww,1), nodes_zero))
  sol % k_d = solution % k_d
  sol % Nj_d = nodes_zero
  sol % ww = solution % ww(:, 1:nodes_zero)

  adapted_grid = Reconstruct_Mesh(grid_zero, adaption_history, params, sol, add_solution)


 END FUNCTION Adapt_Grid                        
       

END MODULE grid_adaption
