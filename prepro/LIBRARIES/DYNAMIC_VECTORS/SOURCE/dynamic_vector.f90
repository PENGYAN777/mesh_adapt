!============================================================ 
!      Module: dynamic_vector
!
! Description: The module provide the definition of the 
!              dynamic integer vector format to be used when 
!              arrays of variable length are to be stored 
!              (for example when storing the connectivity 
!              matrix of a grid made of elements of different
!              kind) 
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

MODULE dynamic_vector

   !============================================================ 
   TYPE D_I_V
   !============================================================ 

      INTEGER, DIMENSION(:), POINTER :: vec

   END TYPE D_I_V
   !============================================================ 

   !============================================================ 
   TYPE D_I_V_R
   !============================================================ 

      REAL(KIND=8), DIMENSION(:), POINTER :: vec

   END TYPE D_I_V_R
   !============================================================ 

   !============================================================ 
   TYPE D_I_M
   !============================================================ 

      TYPE(D_I_V), DIMENSION(:), POINTER :: row

   END TYPE D_I_M
   !============================================================ 

!============================================================ 
CONTAINS
!============================================================

   !============================================================ 
   FUNCTION invert_DIV ( div ) RESULT( inv_div )
   !============================================================ 

   !------------------------------------------------------------ 
   !   This function is used to invert the connectivity relation 
   !   between objects (in DIV format)
   ! 
   ! INPUT:    
   !   div      vector in dynamic integer vector format
   ! 
   ! OUTPUT:   
   !   inv_div  vector in dynamic integer vector format which 
   !            contains the inverse relation between the 
   !            elements of div
   ! 
   ! NOTES: 
   !   The vector in input and output are normalized so has to 
   !   have indices going from one to their size, and hence has 
   !   to be renormalized if different numeration is needed.  
   !   Moreover, if holes are present in the numeration of the 
   !   objects in input, that will be dealt with by generatin an 
   !   output vector with the same holes.
   !------------------------------------------------------------ 


      IMPLICIT NONE

      TYPE(D_I_V), DIMENSION(:), INTENT(IN)  ::  div
      TYPE(D_I_V), DIMENSION(:), POINTER     ::  inv_div

      INTEGER  ::  max_idx_div, min_idx_div
      INTEGER, DIMENSION(:), ALLOCATABLE  ::  inv_div_vec_size 
      INTEGER  ::  i, j

      ! maximum and minimum integer appearing in div 
      ! --------------------------------------------
      max_idx_div = 0
      DO i = 1, size_DIV(div,1)
          max_idx_div = MAX( max_idx_div, MAXVAL(div(i)%vec) )
      ENDDO

      min_idx_div =  max_idx_div 
      DO i = 1, size_DIV(div,1)
         min_idx_div = MIN(min_idx_div, MINVAL(div(i)%vec) )
      ENDDO

      IF ( min_idx_div == 0) THEN
         min_idx_div = 1
      ENDIF
      
      ! Allocation of the inverse div
      ! -----------------------------
      ALLOCATE ( inv_div(max_idx_div-min_idx_div+1) )  

      ! Computation of the size of each dynamic vector in inv_div(:)
      ! ------------------------------------------------------------
      ALLOCATE ( inv_div_vec_size(max_idx_div-min_idx_div+1) )
      inv_div_vec_size = 0
      
      DO i = 1, size_DIV(div,1)
         DO j = 1, SIZE(div(i)%vec)
            IF( div(i)%vec(j) .NE. 0 ) THEN
               inv_div_vec_size( div(i)%vec(j) ) & 
                  = inv_div_vec_size( div(i)%vec(j) ) + 1
            ENDIF 
         ENDDO
      ENDDO

      ! Allocation and initialization of each dynamic 
      ! vector in inv_dir(:) 
      ! ---------------------------------------------
      DO i =1, size_DIV(inv_div,1)
          ALLOCATE ( inv_div(i)%vec(inv_div_vec_size(i)) )  
          inv_div(i)%vec   = 0
      ENDDO

      ! The connectivity relation is now inverted
      ! The vector inv_div_vec_size is used as index 
      ! --------------------------------------------
      DO i = 1, size_DIV(div,1)  
         DO j = 1, SIZE(div(i)%vec)
            IF ( div(i)%vec(j) .NE. 0) THEN  
               inv_div(div(i)%vec(j)) % vec(inv_div_vec_size(div(i)%vec(j))) &
                  = i
               inv_div_vec_size(div(i)%vec(j)) &
                  = inv_div_vec_size(div(i)%vec(j)) - 1

               IF ( inv_div_vec_size(div(i)%vec(j)) < 0 ) THEN
                  WRITE(*,*) 'invert_DIV error in inverting relations'   
                  WRITE(*,*) 'STOP'
                  STOP
               ENDIF 
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE( inv_div_vec_size )
   
   
   END FUNCTION invert_DIV 
   !============================================================ 



   !============================================================ 
   FUNCTION convert_matrix_to_DIV ( matrix ) RESULT( div )
   !============================================================ 

   !------------------------------------------------------------ 
   ! 
   !   This function is used to convert a matrix in a vector in 
   !   DVI format.  Null elements of the matrix ARE stored.
   ! 
   ! INPUT:    
   !   matrix   
   ! 
   ! OUTPUT:   
   !   div      vector in dynamic integer vector format
   ! 
   !------------------------------------------------------------ 


      IMPLICIT NONE

      INTEGER,     DIMENSION(:,:), INTENT(IN)  ::  matrix
      TYPE(D_I_V), DIMENSION(:),   POINTER     ::  div

      INTEGER  ::  i


      ! Allocation of the div vector
      ! ----------------------------
      ALLOCATE ( div(SIZE(matrix,2)) )  

      ! Initialization of the div vector
      ! --------------------------------
      DO i = 1, size_DIV(div,1)
         ALLOCATE ( div(i) % vec(SIZE(matrix,1)) )
         div(i) % vec = matrix(:,i)
      ENDDO


   END FUNCTION convert_matrix_to_DIV
   !============================================================ 



   !============================================================ 
   FUNCTION convert_pack_matrix_to_DIV ( matrix ) RESULT( div )
   !============================================================ 

   !------------------------------------------------------------ 
   !   This function is used to convert a matrix in a vector in 
   !   DVI format.  Null elements of the matrix ARE NOT stored.
   ! 
   ! INPUT:    
   !   matrix   matrix with (hopefully) some null element
   ! 
   ! OUTPUT:   
   !   div      vector in dynamic integer vector format
   ! 
   !------------------------------------------------------------ 


      IMPLICIT NONE

      INTEGER    , DIMENSION(:,:), INTENT(IN)  ::  matrix
      TYPE(D_I_V), DIMENSION(:),   POINTER     ::  div

      INTEGER  ::  i

      ! Allocation of the div vector
      ! ----------------------------
      ALLOCATE ( div(SIZE(matrix,2)) )  

      ! Initialization of the div vector
      ! --------------------------------
      DO i = 1, size_DIV(div,1)
         ALLOCATE ( div(i) % vec(COUNT(matrix(:,i) .NE. 0)) )
         WHERE ( matrix(:,i) .NE. 0 )
            div(i) % vec = matrix(:,i)
         END WHERE 
      ENDDO


   END FUNCTION convert_pack_matrix_to_DIV
   !============================================================ 



   !============================================================ 
   FUNCTION convert_DIV_to_matrix ( div ) RESULT( matrix )
   !============================================================ 

   !------------------------------------------------------------ 
   !
   !   This function is used to convert a vector in DVI format 
   !   in a matrix with fixed dimension. 
   ! 
   ! INPUT:    
   !   div      vector in dynamic integer vector format
   ! 
   ! OUTPUT:   
   !   matrix   matrix with the following dimension
   !               rows:    size_DVI(dvi,1)
   !               columns: size_DVI(dvi,2)
   ! 
   ! NOTES: 
   !   Since the number of columns of the matrix is equal
   !   to the maximum dimension of the vectors of dvi, the 
   !   memory occupied by the matrix is in general greater then 
   !   the one occupied by the starting vector.  Notice that 
   !   indices that are not present in the dvi vector are 
   !   initialized to zero.
   ! 
   !------------------------------------------------------------ 


      IMPLICIT NONE

      TYPE(D_I_V), DIMENSION(:),   INTENT(IN)  ::  div
      INTEGER    , DIMENSION(:,:), POINTER     ::  matrix

      INTEGER  ::  i, j


      ! Allocation of the matrix
      ! ------------------------
      ALLOCATE ( matrix(size_DIV(div,2),size_DIV(div,1)) )  
      matrix = 0

      ! Initialization of the matrix 
      ! ----------------------------
      DO i = 1, size_DIV(div,1)
         DO j = 1, SIZE(div(i)%vec)
            matrix(j,i) = div(i)%vec(j)
         ENDDO
      ENDDO


   END FUNCTION convert_DIV_to_matrix 
   !============================================================ 



   !============================================================ 
   FUNCTION size_DIV ( div,d ) RESULT(s)
   !============================================================ 

   !------------------------------------------------------------ 
   !
   !   Return the size (in the sense specified below), of a 
   !   vector in DIV format 
   ! 
   ! INPUT:    
   !   div      vector in dynamic integer vector format
   !   d        type of size desired:
   !               1 --> size of the array div
   !            2 --> max size of the arrays div % vec
   ! 
   ! OUTPUT:   s        size of div 
   ! 
   !------------------------------------------------------------ 


      IMPLICIT NONE

      TYPE(D_I_V), DIMENSION(:), INTENT(IN)  :: div
      INTEGER  :: d     
      INTEGER  :: s

      INTEGER  :: i, max_idx_div, min_idx_div


      SELECT CASE(d)

         ! Length of the array div of DIV elements
         ! ---------------------------------------
         CASE (1)
         s = SIZE(div)
 
         ! Maximum length of the vectors of the array div
         ! ----------------------------------------------
         CASE (2)
         s = 0
         DO i = 1, SIZE(div)
            s = MAX( s, SIZE( div(i) % vec )  )
         ENDDO

         ! Number of independent indices in div
         ! Zero is not an index: means NONE
         ! ------------------------------------
         CASE (3)
         max_idx_div = 0
         DO i = 1, SIZE(div)
            max_idx_div = MAX( max_idx_div, MAXVAL(div(i)%vec) )
         ENDDO

         min_idx_div =  max_idx_div 
         DO i = 1, SIZE(div)
            min_idx_div = MIN(min_idx_div, MINVAL(div(i)%vec) )
         ENDDO

         IF ( min_idx_div == 0 ) THEN
            min_idx_div = 1
         ENDIF 

         s = max_idx_div-min_idx_div + 1  
 
         CASE DEFAULT
         WRITE(*,*) 'size_DIV error: unknown size type'
         WRITE(*,*) 'STOP'
         STOP

      END SELECT


   END FUNCTION size_DIV
   !============================================================ 


END MODULE dynamic_vector

