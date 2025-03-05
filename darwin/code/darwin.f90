!-----------------------------------------------------------------------------
!   Description: Program for the automatic adaption of computational 
!                grids.  An  edge-based approach and a multiple pas-
!                sages strategy are adopted.  Adapted  grid  quality 
!                is obtained by enforcing techniques of grid smooth-
!                ing and edge swapping.
!
!                An input file named optiMe.cfg is required.
!
!
!        Author: Marco Fossati
!                Department of Aerospace Engineering
!                Politecnico di Milano
!                Via La Masa 34, 20156 Milano, ITALY
!                e-mail: fossati@aero.polimi.it
!
!          Year: 2003-2007
!-----------------------------------------------------------------------------

   PROGRAM  darwin
   
   !--------------------------------------------------------------------------
   USE  grid_adaption
   USE  grid_utils
   USE  io_proc
   USE  solution_interpolation
   USE  structures
 
 
   IMPLICIT NONE
 
   INTEGER           :: entry_level, exit_level
   CHARACTER(LEN=7)  :: ext
   CHARACTER(LEN=30) :: prob_name
   CHARACTER(LEN=30) :: grid_name
 
   TYPE(adaption_param) :: ad_params
   TYPE(grid_type)      :: entry_grid, exit_grid
   TYPE(solution_type)  :: solution, interpolated_solution
   TYPE(solution_type), DIMENSION(:), ALLOCATABLE  :: add_solution
   !--------------------------------------------------------------------------

   CALL write_optiMe_logo(6)

   ad_params = Load_Adaption_Param()
         
   entry_level  = ad_params % entry_level
   exit_level   = ad_params % entry_level + 1 
   prob_name    = ad_params % prob_name
   grid_name    = ad_params % grid_name
       
   WRITE(ext,1000) entry_level
   1000 FORMAT (i7)
   ext = ADJUSTL(ext)
   
   WRITE (*,'(4x,a12,1x,i1,1x,a3)') 'Reading grid', entry_level, '...'
   entry_grid = Load_Grid(entry_level, grid_name, ext)

   solution = Load_Solution ( ad_params % sol_fmt, entry_grid % k_d,  &
     entry_grid % Nj_d, entry_grid % Nb, prob_name )

   ALLOCATE ( add_solution(exit_level) )

   write (*,*)
   write (*,'(4x,a21,1x,i1,1x,a8,1x,i1)') 'Adapt grid from level', entry_level, 'to level', exit_level
   exit_grid = Adapt_Grid(entry_grid, solution, ad_params, add_solution)
!   interpolated_solution = Interp_Solution(solution, add_solution, nodes_zero, exit_grid)
!   interpolated_solution = Interp_Solution_nn(exit_grid, entry_grid, solution, add_solution, nodes_zero)
   interpolated_solution = Interp_Solution_nn_kdtree(exit_grid, entry_grid, solution, add_solution, nodes_zero)


   WRITE(ext,1000) exit_level
   ext = ADJUSTL(ext)

   PRINT*, ''  
   PRINT*, '   Postprocessing'
        
   CALL  Post_Error(entry_grid, ad_params % estimator, grid_name)
   CALL  Write_Grid(exit_grid, grid_name, ext)
   CALL  Post_Grid(exit_grid, grid_name, ad_params % ref_type, ad_params % box)  

!   CALL  Post_Interpolated_Solution(exit_grid, ad_params % sol_fmt, interpolated_solution, prob_name)
   CALL  Write_Interpolated_Solution(exit_grid, ad_params % sol_fmt, exit_grid % Nb, interpolated_solution, prob_name)

   PRINT*, ''   

   CONTAINS

     SUBROUTINE  write_optiMe_logo(idf)

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: idf

     WRITE(idf,*) ''
     WRITE(idf,*) ' -------------------------------------------------------------------'
     WRITE(idf,*) '   DARWIN'
     WRITE(idf,*) ' -------------------------------------------------------------------'
     WRITE(idf,*) '   CFD Laboratory - McGill University'     
     WRITE(idf,*) '   2009-2012 - Marco Fossati'
     WRITE(idf,*) '   Department of Aerospace Engineering - Politecnico di Milano'
     WRITE(idf,*) '   2003-2009 - Alberto Guardone, Marco Fossati'
     WRITE(idf,*) ''

     END SUBROUTINE  write_optiMe_logo

 END PROGRAM  darwin
