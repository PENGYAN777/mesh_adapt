!=============================================================================== 
!
!      Module: Axial Symmetry
!
! Description: Handle axisymmetric problems
!
!      Author: Alberto Guardone, Federico Muccioli
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!    Modified: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!===============================================================================

  MODULE  axial_symmetry

  !-----------------------------------------------------------------
  USE metric_coefficients
  USE euler_equations
  USE csr
  USE csr_pair_sys,   ONLY: add_CSR_ij_sys
  
  IMPLICIT NONE;  PRIVATE
  
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: cell_2d
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: xi_bp_2d
  
  PUBLIC :: rhs_axial, rhs_implicit_axial, &
            np_metric_axial_change, xi_bp_2d 
  !-----------------------------------------------------------------

  CONTAINS


  SUBROUTINE  np_metric_axial_change
  !------------------------------------------------------------ 
  ! Switch from 2D metric quantities to Axysimmetric metric
  ! quantities.  Saves the values of cell (needed for the 
  ! evaluation of the domain pressure term) and of xi_bp 
  ! (needed for the evaluation of the boundary normals)
  !------------------------------------------------------------ 
  IMPLICIT NONE      
  !------------------------------------------------------------ 

  ALLOCATE (cell_2d(SIZE(cell,1)))
  ALLOCATE (xi_bp_2d(SIZE(xi_bp,1), SIZE(xi_bp,2)))
    
  ! Cartesian metric
  cell_2d  = cell
  xi_bp_2d = xi_bp
  
  ! Axisymmetric metric
  cell     = cell_y
  mass     = mass_y 
  eta      = eta_y  
  eta_fv   = eta_y_fv  
  mass_ii  = mass_y_ii
  xi_bp    = xi_y_bp
        
  END SUBROUTINE  np_metric_axial_change





  FUNCTION  rhs_axial(ww)  RESULT(rhs_a)
  !------------------------------------------------------------ 
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:),   INTENT (IN) :: ww
  REAL(KIND=8), DIMENSION(SIZE(ww,2)) :: rhs_a

  INTEGER :: i 
  !------------------------------------------------------------

  ! RHS Contribution
  rhs_a = 0.d0
  DO i = 1, SIZE(ww,2)
    rhs_a(i) = cell_2d(i) * Pi__ww(ww(:,i))
  ENDDO
  
  END FUNCTION  rhs_axial





  FUNCTION rhs_implicit_axial(ww, MM) RESULT (rhs_a) 
  !------------------------------------------------------------ 
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:),   INTENT (IN)    :: ww
  TYPE( CSR_matrix ),             INTENT(INOUT)  :: MM
  REAL(KIND=8), DIMENSION(SIZE(ww,2))            :: rhs_a
  REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(ww,1)) :: MM_ij
  INTEGER :: i
  !------------------------------------------------------------

  ! RHS Contribution
  rhs_a = 0.d0
  DO i = 1 , SIZE(ww,2)
     rhs_a(i) = cell_2d(i) * Pi__ww(ww(:,i))
  ENDDO
  
  ! LHS Contribution
  MM_ij = 0.d0
  DO i = 1, SIZE(ww,2)
     MM_ij(3,:) =  - cell_2d(i) * dPi_dww(ww(:,i)) 
     CALL add_CSR_ij_sys(i, i, MM_ij,  MM)         
  ENDDO
        
  END FUNCTION rhs_implicit_axial


  END MODULE  axial_symmetry
