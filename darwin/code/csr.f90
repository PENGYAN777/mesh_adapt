!============================================================ 
!
!      Module: csr
!
! Description: Definition of the CSR matrix type and 
!              procedures for manipulating sparse matrices in
!              CSR format. 
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

!============================================================ 
!
!  Procedures to store (retrieve) the nonzero elements
!  of a sparse matrix in (from) a global vector.
!  CSR data structure employed.
!
!  --- USAGE ---
!
!  1) Allocate ia(ni+1) and ja(nz), where ni is the number of
!     rows of the sparse matrix and nz is the number of 
!     nonzero entries of the matrix.
!
!  2) To set up the csr data structure:
!
!     2a) Set ia(1:ni) equal to the number of nonzero entries
!         for each row of the sparse matrix (the element 
!         ia(ni+1) is not used at this stage), then CALL 
!         subroutine setia.
!     2b) Initialize ja(:) = HUGE(integer), then CALL 
!         subroutine setja for each nonzero entry (i,j) of 
!         the sparse matrix.
!
!  3) Function csridx_XXX returns the index k of the (i,j) 
!     entry of the sparse matrix in the global data vector.
!     If the coefficients are stored in a vector data(1:nz),
!     data(csridx_XXX(i,j,ia,ja)) is the coefficient (i,j) of 
!     the sparse matrix.
!
!============================================================ 


MODULE csr


   IMPLICIT NONE


   !============================================================ 
   TYPE CSR_matrix
   
      INTEGER,      DIMENSION(:), POINTER  ::  i
      INTEGER,      DIMENSION(:), POINTER  ::  j
      REAL(KIND=8), DIMENSION(:), POINTER  ::  e
   
   END TYPE CSR_matrix
   !============================================================ 





   !============================================================ 
   PUBLIC :: CSR_matrix, setia, setja, csridx_nos, csridx_srt
   !============================================================ 



 !============================================================ 
 CONTAINS
 !============================================================ 



   !============================================================ 
   SUBROUTINE setia ( ia )
   !============================================================ 


      !  Input:  ia(1:n)    number of nonzero coefficients of each row
      !  Output: ia(1:n+1)  ia vector of the csr matrix data structure

      IMPLICIT NONE

      INTEGER(4), DIMENSION(:), INTENT(INOUT) :: ia

      INTEGER(4) :: i, m, n


      n = 1

      DO i = 1, SIZE(ia)
         m = ia(i);  ia(i) = n;  n = n + m
      ENDDO


   END SUBROUTINE setia
   !============================================================ 



   !============================================================ 
   SUBROUTINE setja  ( i, j, ia, ja )
   !============================================================ 


      !  Add to the column index vector ja the coefficient i,j
      !  ja sorted in increasing order


      IMPLICIT NONE

      INTEGER(4),               INTENT(IN)  :: i, j
      INTEGER(4), DIMENSION(:), INTENT(IN)  :: ia
      INTEGER(4), DIMENSION(:), INTENT(OUT) :: ja

      INTEGER(4) :: beg, end, m, n


      beg = ia(i)
      end = ia(i+1)-1

      DO m = beg, end
         IF ( j == ja(m) ) EXIT
         IF ( j < ja(m) ) THEN
            DO n = end, m+1, -1
               ja(n) = ja(n-1)
            ENDDO
            ja(m) = j
            EXIT
         ENDIF
      ENDDO


   END SUBROUTINE setja
   !============================================================ 



   !============================================================ 
   FUNCTION csridx_nos ( i, j, ia, ja ) RESULT ( k )
   !============================================================ 


      !  ja sorting not required


      IMPLICIT NONE

      INTEGER(4),               INTENT(IN) :: i, j
      INTEGER(4), DIMENSION(:), INTENT(IN) :: ia, ja
      INTEGER(4)                           :: k

      INTEGER(4) :: m


      k = 0

      DO m = ia(i), ia(i+1)-1
         IF ( j == ja(m) ) THEN
            k = m
            EXIT
         ENDIF
      ENDDO

      IF ( k == 0 ) THEN
         WRITE (*,*) 'csridx_nos error, column not found. STOP'
         STOP
      ENDIF

   END FUNCTION csridx_nos
   !============================================================ 



   !============================================================ 
   FUNCTION csridx_srt  ( i, j, ia, ja ) RESULT ( m )
   !============================================================ 


      !  ja must be sorted in increasing order


      IMPLICIT NONE

      INTEGER(4),               INTENT(IN) :: i, j
      INTEGER(4), DIMENSION(:), INTENT(IN) :: ia, ja
      INTEGER(4)                           :: m

      INTEGER :: n, k, kk, kl, kr


      m = 0;  kk = 1;  kl = ia(i);  kr = ia(i+1)-1
      DO n = ia(i), ia(i+1)-1
         kk = ABS(kk-1); k = (kl+kr)/2 + kk
         IF ( ja(k) > j ) THEN
            kr = k
         ELSEIF ( ja(k) < j ) THEN
            kl = k
         ELSE
            m = k;  EXIT
         ENDIF
      ENDDO

      IF ( m == 0 ) THEN
         WRITE (*,*) 'csridx_srt error, column not found. STOP'
         STOP
      ENDIF

   END FUNCTION csridx_srt
   !============================================================ 



!   !============================================================ 
!   SUBROUTINE prtija ( ia, ja )
!   !============================================================ 
!
!
!      IMPLICIT NONE
!
!      INTEGER(4), DIMENSION(:), INTENT(IN) :: ia, ja
!      INTEGER(4) :: i, j
!
!
!      PRINT '(2a5,5x,a4)', '    I','   IA','  JA'
!
!      DO i = 1, SIZE(ia)-1
!         PRINT '(2i5,5x,20i4)', i, ia(i), &
!                               (ja(j),j=ia(i),ia(i+1)-1)
!      ENDDO
!
!
!   END SUBROUTINE prtija
!   !============================================================ 




END MODULE csr
