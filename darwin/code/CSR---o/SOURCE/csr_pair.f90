!============================================================ 
!
!      Module: csr_pair
!
! Description: procedure to manipulate a sparse matrix in CSR
!              format based on node-pair topology. 
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

MODULE csr_pair


   !============================================================ 
   USE dynamic_vector
   USE csr
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   SUBROUTINE alloc_CSR_pair(j_c_matrix, MM)
   !============================================================ 


      IMPLICIT NONE

      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c_matrix
      TYPE( CSR_matrix )  ::  MM      

      INTEGER,     DIMENSION(:,:), ALLOCATABLE  ::  j_c_mat_SWP
      TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE  ::  j_c, c_j
      INTEGER  :: i, j , c, c_


      ALLOCATE(j_c_mat_SWP(SIZE(j_c_matrix,1), SIZE(j_c_matrix,2)) )
      j_c_mat_SWP = j_c_matrix
      j_c_mat_SWP(3:4,:) = 0 

      ALLOCATE( j_c(SIZE(j_c_mat_SWP,2)) ) 
      j_c = convert_matrix_to_DIV(j_c_mat_SWP) 
      DEALLOCATE (j_c_mat_SWP)

      ALLOCATE( c_j(size_DIV(j_c,3)) ) 
      c_j = invert_DIV(j_c)
   

      ! Set the index vector ia
      ! -----------------------
      ALLOCATE( MM%i(SIZE(c_j) + 1) )
      DO i = 1, SIZE(c_j)
         MM%i(i) = SIZE(c_j(i)%vec) + 1
      ENDDO

      CALL setia(MM%i)


      ! Set the index vector ja
      ! -----------------------
      ALLOCATE( MM%j(MM%i(SIZE(MM%i))-1) )
      MM%j = HUGE(i)

      ! Loop on all nodes
      DO i = 1, SIZE(c_j)

         ! Loop on the node-pair bubble of i
         DO c_ = 1, SIZE(c_j(i)%vec) 

            c = c_j(i)%vec(c_)

               IF (j_c(c)%vec(1) == i) THEN
                  j = j_c(c)%vec(2)
               ELSE
                  j = j_c(c)%vec(1)
               ENDIF

               CALL setja(i, j , MM%i, MM%j)

         ENDDO

         ! Add the diagonal element ii 
         CALL setja(i,i,MM%i,MM%j)


      ENDDO

            
      ! Allocate the data vector a
      ! --------------------------

      ALLOCATE( MM%e(SIZE(MM%j)) ) 
      MM%e = 0.d0

      DEALLOCATE (j_c, c_j)

      
   END SUBROUTINE alloc_CSR_pair
   !============================================================ 



   !============================================================ 
   SUBROUTINE store_CSR_pair(j_c, a_ii, a_ij, MM)
   !============================================================ 

      IMPLICIT NONE

      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  a_ii, a_ij
      TYPE( CSR_matrix )  ::  MM
      
      INTEGER  :: i, j , c

      !--- Set the data vector a
      !-------------------------
      
      ! Loop on all node-pair
      DO c = 1, SIZE(j_c,2)

            i = j_c(1,c);  j = j_c(2,c)
            MM%e(csridx_srt(i,j,MM%i,MM%j)) = a_ij(c) 
            MM%e(csridx_srt(j,i,MM%i,MM%j)) = a_ij(c) 

      ENDDO

      DO i = 1, SIZE(a_ii)

         ! Add the diagonal element ii 
         MM%e(csridx_srt(i,i,MM%i,MM%j)) = a_ii(i) 

      ENDDO


      
   END SUBROUTINE store_CSR_pair
   !============================================================ 



END MODULE csr_pair 
