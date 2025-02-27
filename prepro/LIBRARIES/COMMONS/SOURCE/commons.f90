!=============================================================================== 
!
!      Module: Commons
!
! Description: Collection of global variables
!
!      Author: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!===============================================================================

  MODULE commons
  
  IMPLICIT NONE
  
  ! Problem Definitions
  CHARACTER(len=64) :: pb_name
  CHARACTER(len=64) :: gd_name 
  
  INTEGER :: pb_name_l
  INTEGER :: gd_name_l 
  INTEGER :: nEqs
  INTEGER :: sDim
  
  LOGICAL :: steady_flow
  LOGICAL :: restart
  LOGICAL :: axisym
  LOGICAL :: viscous_flow
  LOGICAL :: turbulent_flow
  LOGICAL :: MP_job = .FALSE.  
  LOGICAL :: MP_master
  LOGICAL :: MP_worker
    
  END MODULE commons  
