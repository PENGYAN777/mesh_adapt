MODULE derivatives  
  
 USE metric_coefficients
 USE node_pair_structure
 USE nodes
 USE np_quadrature,       ONLY : np_quadr_w_Gu, np_quadr_w_Hu  
 USE structures
  
 CONTAINS  

   SUBROUTINE  gradient (grid, v, G_v)
   !-----------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type),            INTENT(IN) :: grid
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: v

   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: G_v 

   INTEGER :: j
   !-----------------------------------------------------------------------

   G_v  = np_quadr_w_Gu(grid % jcd_fem, grid % jcb_fem, grid % jd_jb, &
                        grid % eta, grid % chi_b, grid % xi_bp, v)

   DO j = 1, grid % Nj_d
     G_v(:,j) = G_v(:,j) / grid % cell(j) 
   ENDDO

   END SUBROUTINE  gradient





   SUBROUTINE  hessian(grid, v, H_v)
   !-------------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type),            INTENT(IN) :: grid
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: v

   REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: H_v 

   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: G_v 

   REAL(KIND=8) :: mixed_deriv 

   INTEGER :: i, k  
   !-------------------------------------------------------------------------------
 
   ALLOCATE (G_v(grid % k_d, grid % Nj_d)) 
 
   G_v  = np_quadr_w_Gu(grid % jcd_fem, grid % jcb_fem, grid % jd_jb, &
                        grid % eta, grid % chi_b, grid % xi_bp, v)

   DO i = 1, grid % Nj_d
     G_v(:,i) = G_v(:,i) / grid % cell(i)
   ENDDO

   DO k = 1, grid % k_d
     H_v(k,:,:) = np_quadr_w_Gu(grid % jcd_fem, grid % jcb_fem, grid % jd_jb,   &
                                grid % eta, grid % chi_b, grid % xi_bp, G_v(k,:))
   ENDDO
   
   DO i = 1, grid % Nj_d
     H_v(:,:,i) = H_v(:,:,i) / grid % cell(i)
   ENDDO

   ! Forcing equality of mixed derivatives
   
   DO i = 1, grid % Nj_d

     mixed_deriv = (H_v(1,2,i) + H_v(2,1,i)) / 2.d0
   
     H_v(1,2,i) = mixed_deriv
     H_v(2,1,i) = mixed_deriv

   ENDDO    

   END SUBROUTINE  hessian





   SUBROUTINE  hessian_(grid, v, H_v)
   !--------------------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type),            INTENT(IN) :: grid
   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: v

   REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT) :: H_v

   INTEGER :: j
   !--------------------------------------------------------------------------------------

   H_v = np_quadr_w_Hu(grid % jcd_fem, grid % jcb_fem, grid % jd_jb, &
                       grid % stiff_T, grid % Kb_ijs, grid % Kb_ija, v)

   DO j = 1, grid % Nj_d
     H_v(:,:,j) = H_v(:,:,j) / grid % cell(j)
   ENDDO

! do j = 1, grid % Nj_d
! print*, j
! print*, H_v(1,:,j)
! print*, H_v(2,:,j)
! enddo

   END SUBROUTINE  hessian_

END MODULE derivatives   
