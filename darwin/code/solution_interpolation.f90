MODULE solution_interpolation

 USE structures
 USE kdtree2_module
 USE kdtree2_precision_module

 CONTAINS


   FUNCTION  Interp_Solution(solution, loc_sol, nodes_zero, new_grid)   RESULT(interpolated_solution)
   !-----------------------------------------------------------------------------------------------------
   IMPLICIT NONE
 
   TYPE(solution_type), INTENT(IN) :: solution
   TYPE(solution_type), DIMENSION(:), INTENT(IN) :: loc_sol
   INTEGER,             INTENT(IN) :: nodes_zero
   TYPE(grid_type),     INTENT(IN) :: new_grid
 
   TYPE(solution_type) :: interpolated_solution

   INTEGER :: i, j, n
   !-----------------------------------------------------------------------------------------------------

   interpolated_solution % k_d  = solution % k_d 
   interpolated_solution % Nj_d = new_grid % Nj_d
  
   ALLOCATE (interpolated_solution % ww(SIZE(solution % ww, 1), new_grid % Nj_d))
  
   ! Copying solution on grid zero
   interpolated_solution % ww = 0.d0
   interpolated_solution % ww(:, 1:nodes_zero) = solution % ww(:, 1:nodes_zero)

   n = nodes_zero 

   DO i = 1, SIZE(loc_sol)

      DO j = 1, COUNT(loc_sol(i) % ww(1,:) /= 0)

         n = n + 1
         interpolated_solution % ww(:,n) = loc_sol(i) % ww(:,j)

      ENDDO

   ENDDO

   END FUNCTION  Interp_Solution

	

	FUNCTION Interp_Solution_nn(new_grid, grid, solution, loc_sol, nodes_zero) RESULT(interpolated_solution)
	
	IMPLICIT NONE

	   TYPE(grid_type),     INTENT(IN) :: new_grid, grid
	   TYPE(solution_type), INTENT(IN) :: solution
	   TYPE(solution_type), DIMENSION(:), INTENT(IN) :: loc_sol
	   INTEGER,             INTENT(IN) :: nodes_zero

	   TYPE(solution_type) :: interpolated_solution
	   INTEGER :: i, k
	   REAL(8), DIMENSION(SIZE(solution%ww,2)) :: dist

	   interpolated_solution % k_d  = solution % k_d 
	   interpolated_solution % Nj_d = new_grid % Nj_d
	  
	   ALLOCATE (interpolated_solution % ww(SIZE(solution % ww, 1), new_grid % Nj_d))
	  
	   interpolated_solution % ww = 0.d0
	
		IF(SIZE(loc_sol)==1) THEN

			interpolated_solution % ww(:, 1:nodes_zero) = solution % ww(:, 1:nodes_zero)	
			k = nodes_zero 

			DO i = 1, COUNT(loc_sol(1) % ww(1,:) /= 0)
				 k = k + 1
				 interpolated_solution % ww(:,k) = loc_sol(1) % ww(:,i)
			ENDDO

		ELSE

			DO i = 1, new_grid % Nj_d	
	
				DO k = 1, SIZE(solution%ww,2)
					dist(k) = NORM2(new_grid%rr(:,i)-grid%rr(:,k))
				ENDDO
				k = MINLOC(dist,1)
		
				interpolated_solution % ww(:,i) = solution % ww(:,k)

			ENDDO

		ENDIF
	

	END FUNCTION Interp_Solution_nn



	FUNCTION Interp_Solution_nn_kdtree(new_grid, grid, solution, loc_sol, nodes_zero) RESULT(interpolated_solution)
	
	IMPLICIT NONE

	   TYPE(grid_type),     INTENT(IN) :: new_grid, grid
	   TYPE(solution_type), INTENT(IN) :: solution
	   TYPE(solution_type), DIMENSION(:), INTENT(IN) :: loc_sol
	   INTEGER,             INTENT(IN) :: nodes_zero

	   TYPE(solution_type) :: interpolated_solution
	   TYPE(kdtree2), POINTER :: tree
	   TYPE(kdtree2_result), ALLOCATABLE :: results(:)
	   INTEGER :: i, k
	   REAL(8), DIMENSION(SIZE(solution%ww,2)) :: dist

	   interpolated_solution % k_d  = solution % k_d 
	   interpolated_solution % Nj_d = new_grid % Nj_d
	  
	   ALLOCATE (interpolated_solution % ww(SIZE(solution % ww, 1), new_grid % Nj_d))
	  
	   interpolated_solution % ww = 0.d0
	
		IF(SIZE(loc_sol)==1) THEN

			interpolated_solution % ww(:, 1:nodes_zero) = solution % ww(:, 1:nodes_zero)	
			k = nodes_zero 

			DO i = 1, COUNT(loc_sol(1) % ww(1,:) /= 0)
				 k = k + 1
				 interpolated_solution % ww(:,k) = loc_sol(1) % ww(:,i)
			ENDDO

		ELSE
			
			tree => kdtree2_create(grid%rr)
			ALLOCATE(results(1))			

			DO i = 1, new_grid % Nj_d	
	
				CALL kdtree2_n_nearest(tree, new_grid%rr(:,i), 1, results)

				interpolated_solution % ww(:,i) = solution % ww(:,results(1)%idx)
	
			ENDDO

		ENDIF
	

	END FUNCTION Interp_Solution_nn_kdtree


END MODULE solution_interpolation
