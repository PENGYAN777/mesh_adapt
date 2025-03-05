!============================================================ 
!
!      Module: csr_pair_sys
!
! Description: procedure to manipulate a sparse matrix in CSR
!              format based on node-pair topology for a 
!              system of equations.
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

MODULE csr_pair_sys


   !============================================================ 
   USE dynamic_vector
   USE csr
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   SUBROUTINE alloc_CSR_pair_sys(j_c_matrix, sys_dim,  MM)
   !============================================================ 


      IMPLICIT NONE

      INTEGER,      DIMENSION(:,:), INTENT(IN)     ::  j_c_matrix
      INTEGER,                      INTENT(IN)     ::  sys_dim
      TYPE(CSR_matrix),             INTENT(INOUT)  ::  MM

      INTEGER,     DIMENSION(2,SIZE(j_c_matrix,2))  ::  j_c_mat_SWP
      TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE  ::  j_c, c_j
      INTEGER  :: i, j, n, m, c, c_, p, iM, jM


      j_c_mat_SWP = j_c_matrix(1:2,:)

      ALLOCATE( j_c(SIZE(j_c_mat_SWP,2)) ) 
      j_c = convert_matrix_to_DIV(j_c_mat_SWP) 
 
      ALLOCATE( c_j(size_DIV(j_c,3)) ) 
      c_j = invert_DIV(j_c)
   
      ! Set the index vector ia
      ! -----------------------
      ALLOCATE( MM%i(SIZE(c_j)*sys_dim + 1) )
      DO i = 1, SIZE(c_j)
         DO p = 1, sys_dim
            iM = (i-1)*sys_dim + p
            MM%i(iM) = (SIZE(c_j(i)%vec) + 1)*sys_dim  
        ENDDO
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
 
               DO n = 1, sys_dim 
                  DO m = 1, sys_dim
                     iM = (i-1)*sys_dim + n  
                     jM = (j-1)*sys_dim + m  
                     CALL setja(iM , jM, MM%i, MM%j)
                  ENDDO
               ENDDO
               
         ENDDO

         ! Add the diagonal element ii 
        
         DO n = 1, sys_dim 
            DO m = 1, sys_dim
       
            iM = (i-1)*sys_dim + n  
            jM = (i-1)*sys_dim + m  
            CALL setja(iM , jM, MM%i, MM%j)
       
            ENDDO
         ENDDO
       

      ENDDO

            
      ! Allocate the data vector a
      ! --------------------------
      ALLOCATE( MM%e(SIZE(MM%j)) ) 
      MM%e = 0.d0

      DEALLOCATE (j_c, c_j)

      
   END SUBROUTINE alloc_CSR_pair_sys
   !============================================================ 



   !============================================================ 
   SUBROUTINE add_CSR_ij_sys(i, j, M_ij, MM)
   !============================================================ 

      IMPLICIT NONE

      INTEGER,                      INTENT(IN)     ::  i, j
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  M_ij
      TYPE(CSR_matrix),             INTENT(INOUT)  ::  MM
      
      INTEGER  :: n, m, idx, iM, jM, sys_dim



            
      !--- Set the data vector a
      !-------------------------
      
      sys_dim = SIZE(M_ij,1)
      DO n = 1, sys_dim
         DO m = 1, sys_dim 
            iM = (i-1)*sys_dim + n  
            jM = (j-1)*sys_dim + m  
            idx = csridx_srt(iM, jM, MM%i, MM%j)
            MM%e(idx) = MM%e(idx) + M_ij(n,m)
         ENDDO
      ENDDO


      
   END SUBROUTINE add_CSR_ij_sys
   !============================================================ 



   !============================================================ 
   SUBROUTINE set_CSR_ij_sys(i, j, M_ij, MM)
   !============================================================ 

      IMPLICIT NONE

      INTEGER,                      INTENT(IN)     ::  i, j
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  M_ij
      TYPE(CSR_matrix),             INTENT(INOUT)  ::  MM
      
      INTEGER  :: n, m, idx, iM, jM, sys_dim



            
      !--- Set the data vector a
      !-------------------------
      

      sys_dim = SIZE(M_ij,1)
      DO n = 1, sys_dim
         DO m = 1, sys_dim 
            iM = (i-1)*sys_dim + n  
            jM = (j-1)*sys_dim + m 
            idx = csridx_srt(iM, jM, MM%i, MM%j)
            MM%e(idx) = M_ij(n,m)
         ENDDO
      ENDDO


      
   END SUBROUTINE set_CSR_ij_sys
   !============================================================ 



   !============================================================ 
   SUBROUTINE set_CSR_sys(alpha, MM)
   !============================================================ 

      IMPLICIT NONE

      REAL(KIND=8),     INTENT(IN)     ::  alpha
      TYPE(CSR_matrix), INTENT(INOUT)  ::  MM
      
      MM%e = alpha
      
   END SUBROUTINE set_CSR_sys
   !============================================================ 


END MODULE csr_pair_sys 
