!============================================================ 
!
!      Module: csr_dirichlet
!
! Description: Dirichlet boundary conditions for sparse 
!              matrices in CSR format. 
!
!     Credits: initial version by Stefano Rebay
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!   Copyright: 1998-2003 Alberto Guardone
!              See COPYING file for copyright notice
!
!============================================================ 

MODULE csr_dirichlet

   !============================================================ 
   USE csr
   !============================================================ 

!============================================================ 
CONTAINS
!============================================================ 


   !============================================================ 
   SUBROUTINE dirichlet_CSR( dir_nodes, AA ) 
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER,      DIMENSION(:), INTENT(IN)     ::  dir_nodes
      TYPE(CSR_matrix) :: AA

      !------------------------------------------------------------ 
      INTEGER  ::  i, j, n, d_n
      !------------------------------------------------------------ 


      DO n = 1, SIZE(dir_nodes); d_n = dir_nodes(n)

         DO i = AA%i(d_n), AA%i(d_n + 1) - 1

            j = AA%j(i)

            IF (j == d_n) THEN 
               AA%e(i) = 1
            ELSE 
               AA%e(i) = 0
            ENDIF

         ENDDO 

      ENDDO


   END SUBROUTINE dirichlet_CSR
   !============================================================ 


END MODULE csr_dirichlet
