!============================================================ 
!
!    Author:  Alberto Guardone
!             Dipartimento di Ingegneria Aerospaziale
!             Politecnico di Milano
!             Via La Msa 34, 20156 Milano, ITALY
!             e-mail: guardone@aero.polimi.it
!
! Copyright: 1998-2003 Alberto Guardone
!            See COPYING file for copyright notice
!
!============================================================ 
MODULE coo


   USE csr

   IMPLICIT NONE


   !============================================================ 
   TYPE COO_matrix
   
      INTEGER,      DIMENSION(:), POINTER  ::  i
      INTEGER,      DIMENSION(:), POINTER  ::  j
      REAL(KIND=8), DIMENSION(:), POINTER  ::  e
   
   END TYPE COO_matrix
   !============================================================ 


   !============================================================ 
   PUBLIC :: convert_CSR_to_COO
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   FUNCTION convert_CSR_to_COO(CSR_M) RESULT(COO_M)
   !============================================================ 



      IMPLICIT NONE

      TYPE(CSR_matrix),  INTENT(IN)  ::  CSR_M

      TYPE(COO_matrix)               ::  COO_M

      INTEGER  ::  Nz, i, j , e
      

      Nz = SIZE(CSR_M%j)
      ALLOCATE( COO_M%i(Nz), COO_M%j(Nz), COO_M%e(Nz) )

      DO i = 1, SIZE(CSR_M%i) - 1   
         DO e = CSR_M%i(i), CSR_M%i(i+1) - 1;  j = CSR_M%j(e) 

             COO_M%i(e) = i  
             COO_M%j(e) = j  

         ENDDO
      ENDDO

      COO_M%e = CSR_M%e
   
   END FUNCTION convert_CSR_to_COO
   !============================================================ 



END MODULE coo
