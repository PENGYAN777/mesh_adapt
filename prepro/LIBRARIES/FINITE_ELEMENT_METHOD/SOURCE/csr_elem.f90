!============================================================ 
!
!      Module: csr_elem
!
! Description: Initialize the sparse matrix in CSR format
!              from a finite element discretization  
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


MODULE csr_elem


   !============================================================ 
   USE dynamic_vector
   USE csr
   USE fem_ele_types
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   SUBROUTINE alloc_CSR_elem(ele_d, AA)
   !============================================================ 


      IMPLICIT NONE

      TYPE(element_domain), DIMENSION(:), INTENT(IN)  ::  ele_d 
      TYPE(CSR_matrix) :: AA

      TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE  ::  n_m, m_n
      LOGICAL, DIMENSION(:),  ALLOCATABLE  ::  visited
      INTEGER  :: n_tot
      INTEGER  :: i, j, j_ , m, m_


      ALLOCATE( n_m(SIZE(ele_d)) )
      DO m = 1, SIZE( ele_d )

         ALLOCATE ( n_m(m)%vec(SIZE(ele_d(m)%nn)) )
         n_m(m)%vec = ele_d(m)%nn

      ENDDO

      ALLOCATE( m_n(size_DIV(n_m,3)) ) 
      m_n = invert_DIV(n_m)
   

      ALLOCATE( visited(SIZE(m_n)) ) 

      ! Set the index vector ia
      ! -----------------------
      ALLOCATE( AA%i(SIZE(m_n) + 1) )
      DO i = 1, SIZE(m_n)

         visited = .FALSE.

         n_tot = 1;  visited(i) = .TRUE.  ! Matrix diagonal element ii 

         DO m_ = 1, SIZE(m_n(i)%vec);    m = m_n(i)%vec(m_)

            DO j_ = 1, SIZE(n_m(m)%vec); j = n_m(m)%vec(j_)

               IF (.NOT. visited(j)) THEN
                  n_tot = n_tot + 1 ! Matrix element ij
                  visited(j) = .TRUE.
               ENDIF

            ENDDO

         ENDDO

         AA%i(i) = n_tot

      ENDDO

      CALL setia(AA%i)


      ! Set the index vector ja
      ! -----------------------
      ALLOCATE( AA%j(AA%i(SIZE(AA%i))-1) )
      AA%j = HUGE(i)

      ! Loop on all nodes
      DO i = 1, SIZE(m_n)

         ! Loop on the node-pair bubble of i
         DO m_ = 1, SIZE(m_n(i)%vec);    m = m_n(i)%vec(m_)

            DO j_ = 1, SIZE(n_m(m)%vec); j = n_m(m)%vec(j_)

               CALL setja(i, j , AA%i, AA%j)
              
            ENDDO

         ENDDO

         ! Add the diagonal element ii 
         CALL setja(i,i,AA%i,AA%j)


      ENDDO

            
      ! Allocate the data vector a
      ! --------------------------

      ALLOCATE( AA%e(SIZE(AA%j)) ) 
      AA%e = 0.d0

      DEALLOCATE(n_m, m_n, visited)
      
   END SUBROUTINE alloc_CSR_elem
   !============================================================ 




END MODULE csr_elem
