!============================================================ 
!
!      Module: lin_algebra
!
! Description: procedures for (very simple) linear algebra
!              (matrices and vectors)
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

MODULE lin_algebra


!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   FUNCTION determinant (A) RESULT (det)
   !============================================================ 


      IMPLICIT NONE
     
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  A
      REAL(KIND=8)  ::  det 


      IF ( SIZE(A,1) .NE. SIZE(A,2) ) THEN
         WRITE(*,*) ' Matrix A is not square.' 
         WRITE(*,*) ' in FUNCTION determinant, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
      ENDIF

      SELECT CASE ( SIZE(A,1) )

      CASE (1)
         det = A(1,1)

      CASE (2)
         det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      CASE (3)
         det = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3))  &
             - A(1,2)*(A(2,1)*A(3,3) - A(3,1)*A(2,3))  &
             + A(1,3)*(A(2,1)*A(3,2) - A(3,1)*A(2,2))

      CASE DEFAULT 
         WRITE(*,*) ' Matrix of size ', SIZE(A,1), ' not implemented' 
         WRITE(*,*) ' in FUNCTION determinant, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
         ! CALL ......

      END SELECT


   END FUNCTION determinant 
   !============================================================ 



   !============================================================ 
   FUNCTION inverse (A) RESULT (B)
   !============================================================ 


      IMPLICIT NONE
     
      REAL(KIND=8), DIMENSION(:,:),     INTENT(IN)  ::  A
      REAL(KIND=8), DIMENSION(SIZE(A,1),SIZE(A,2))  ::  B

      REAL(KIND=8)  ::  det


      IF ( SIZE(A,1) .NE. SIZE(A,2) ) THEN
         WRITE(*,*) ' Matrix A is not square.' 
         WRITE(*,*) ' in FUNCTION inverse, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
      ENDIF

      det = determinant(A)

      IF ( det == 0 ) THEN
         WRITE(*,*) ' Matrix A is singular.' 
         WRITE(*,*) ' in FUNCTION inverse, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
      ENDIF

      SELECT CASE ( SIZE(A,1) )

      CASE (1)
      ! ------
         B(1,1) = 1.d0/det

      CASE (2) 
      ! ------
         B(1,1) =   A(2,2);  B(1,2) = - A(1,2) 
         B(2,1) = - A(2,1);  B(2,2) =   A(1,1) 
         B = B/det

      CASE (3)
      ! ------
         B(1,1) =   A(2,2)*A(3,3) - A(3,2)*A(2,3)
         B(1,2) = - A(1,2)*A(3,3) + A(3,2)*A(1,3)
         B(1,3) =   A(1,2)*A(2,3) - A(2,2)*A(1,3)

         B(2,1) = - A(2,1)*A(3,3) + A(3,1)*A(2,3)
         B(2,2) =   A(1,1)*A(3,3) - A(3,1)*A(1,3)
         B(2,3) = - A(1,1)*A(2,3) + A(2,1)*A(1,3)

         B(3,1) =   A(2,1)*A(3,2) - A(3,1)*A(2,2)
         B(3,2) = - A(1,1)*A(3,2) + A(3,1)*A(1,2)
         B(3,3) =   A(1,1)*A(2,2) - A(2,1)*A(1,2)

         B = B/det

      CASE DEFAULT
         WRITE(*,*) ' Matrix of size ', SIZE(A,1), ' not implemented' 
         WRITE(*,*) ' in FUNCTION inverse, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
         ! CALL ......

      END SELECT


   END FUNCTION inverse
   !============================================================ 



   !============================================================ 
   FUNCTION vec_prod (VV) RESULT (VN)
   !============================================================ 


      !------------------------------------------------------------
      !
      ! Generalized vector product in N, N = SIZE(VV,2) 
      ! dimensions.  For example, in 3D, 
      !
      !        [  i    j    k  ]  
      !        [               ]
      !  AA =  [ V_x  V_y  V_z ] 
      !        [               ]
      !        [ U_x  U_y  U_z ]
      !
      ! and one retrieve the standard vector product
      !
      !------------------------------------------------------------


      !------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------
      REAL(KIND=8), DIMENSION(:,:),       INTENT(IN)  ::  VV

      REAL(KIND=8), DIMENSION(SIZE(VV,2))             ::  VN
      !------------------------------------------------------------
      REAL(KIND=8), DIMENSION(SIZE(VV,1)+1,SIZE(VV,2))  ::  AA 
      INTEGER  ::  k
      !------------------------------------------------------------
      
      
      IF ( SIZE(VV,1)+1 .NE. SIZE(VV,2) ) THEN      
         WRITE(*,*) ' Inconsistent size of the vectors.' 
         WRITE(*,*) ' in FUNCTION vec_prod, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP     
      ENDIF
      
      
      AA(1,:) = 1.d0
      AA(2:SIZE(VV,1)+1,:) = VV 

      DO k = 1, SIZE(VN)  

         VN(k) = ((-1.d0)**(1+k)) * minor(AA, 1, k)

      ENDDO


   END FUNCTION vec_prod
   !============================================================ 

 
 
   !============================================================ 
   FUNCTION minor(A, i, j) RESULT (MNR)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  A
      INTEGER,                      INTENT(IN)  ::  i, j

      REAL(KIND=8)                              ::  MNR 
      !------------------------------------------------------------
      REAL(KIND=8), DIMENSION(SIZE(A,1)-1,SIZE(A,2)-1)  ::  MN
      INTEGER  ::  m, n
      !------------------------------------------------------------


      ! Check consistency
      ! -----------------
      IF ( SIZE(A,1) .NE. SIZE(A,2) ) THEN
         WRITE(*,*) ' Matrix A is not square.' 
         WRITE(*,*) ' in FUNCTION minor, MODULE lin_algebra' 
         WRITE(*,*) ' STOP. ' 
         STOP
      ENDIF

      DO m = 1, i-1
         DO n = 1, j-1
            MN(m,n) = A(m,n)
         ENDDO
         DO n = j+1, SIZE(A,2)
            MN(m,n-1) = A(m,n)
         ENDDO
      ENDDO

      DO m = i+1, SIZE(A,1)
         DO n = 1, j-1
            MN(m-1,n) = A(m,n)
         ENDDO
         DO n = j+1, SIZE(A,2)
            MN(m-1,n-1) = A(m,n)
         ENDDO
      ENDDO

      MNR = determinant(MN)

   END FUNCTION minor
   !============================================================ 



END MODULE lin_algebra
