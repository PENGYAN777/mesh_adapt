!=============================================================================== 
!
!      Module: structures
!
! Description: Variables type definitions
!
!      Author: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!
!   Copyright: 2003-2006 Marco Fossati
!              See COPYING file for copyright notice
!
!===============================================================================

  MODULE  structures


  IMPLICIT NONE


  TYPE reference_values
    REAL(KIND=8) :: L,  P,  v,    Rgas
    REAL(KIND=8) :: T,  u,  cp,   cv,  Z
    REAL(KIND=8) :: mu, k,  Re,   Pr
  END TYPE reference_values



  TYPE  eos_type
    REAL(KIND=8) :: T, v, e, P
  END TYPE  eos_type


  TYPE  eos_ext_type  
    REAL(KIND=8) :: T, v, e, P, c
    REAL(KIND=8), DIMENSION(2) :: dP    
  END TYPE  eos_ext_type


  ! VARIABLES DEFINITIONS
  TYPE(reference_values) :: REF
  
  END MODULE  structures
