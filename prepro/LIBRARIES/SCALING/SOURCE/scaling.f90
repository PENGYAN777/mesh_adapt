!=============================================================================== 
!
!      Module: scaling
!
! Description: Reference Values for scaling equations
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

  MODULE  scaling
  
  !----------------------------------------------------------------------------
  USE structures,   ONLY: REF
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_param_scaling,  &
	    write_param_scaling
  !----------------------------------------------------------------------------

  CONTAINS

  
  SUBROUTINE  read_param_scaling(idf)  
  !-----------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: idf
  !-----------------------------------------------------------------------------
  
  READ(idf,*);  READ(idf,*);  READ(idf,*)

  READ(idf,*) REF % L
  READ(idf,*) REF % P
  READ(idf,*) REF % T
  READ(idf,*) REF % Rgas

  END SUBROUTINE  read_param_scaling
  
  
  
  SUBROUTINE  write_param_scaling(idf)
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: idf
  !----------------------------------------------------------------------------

  WRITE(idf,*) REF % L, ' - [m] length'
  WRITE(idf,*) REF % P, ' - [Pa] pressure'
  WRITE(idf,*) REF % T, ' - [K] temperature'
  WRITE(idf,*) REF % Rgas, ' - [J/(kg K) gas constant]'
  
  END SUBROUTINE  write_param_scaling
  

  END MODULE  scaling
