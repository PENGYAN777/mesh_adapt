!============================================================ 
!
!      Module: fem_inlet
!
! Description: impose strong inlet boundary conditions 
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

MODULE fem_inlet

   !============================================================ 
   USE fem_ele_types
   USE csr
   !============================================================ 

!============================================================ 
CONTAINS
!============================================================ 


   !============================================================ 
   FUNCTION inlet_nodes( jd_jb, vv, normal) RESULT (in_nodes) 
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER,      DIMENSION(:),   INTENT(IN)     ::  jd_jb
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  vv
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  normal

      INTEGER,      DIMENSION(:), POINTER        :: in_nodes
      !------------------------------------------------------------ 
      INTEGER  ::  i, n, N_inlet
      LOGICAL, DIMENSION(SIZE(jd_jb))  ::  inlet
      !------------------------------------------------------------ 

      inlet = .FALSE.
      
      N_inlet = 0
      DO i = 1, SIZE(jd_jb)
      
         IF ( SUM( vv(:,jd_jb(i)) * normal(:,i) ) < 1.e-12 ) THEN
             inlet(i) = .TRUE.
             N_inlet = N_inlet + 1
         ENDIF
      
      ENDDO

      ALLOCATE ( in_nodes(N_inlet) )
      
      n = 1
      DO i = 1, SIZE(jd_jb)
      
         IF (inlet(i)) THEN
             in_nodes(n) = jd_jb(i)
             n = n + 1
         ENDIF
         
      ENDDO
      

   END FUNCTION inlet_nodes
   !============================================================ 


   !============================================================ 
   FUNCTION nodal_normals( jd_jb, ele_b ) RESULT (normal) 
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER,                 DIMENSION(:),   INTENT(IN)  ::  jd_jb
      TYPE( element_boundary), DIMENSION(:),   INTENT(IN)  :: ele_b

      REAL(KIND=8),            DIMENSION(:,:), POINTER     ::  normal
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  rnorms
      INTEGER,      DIMENSION(:),     POINTER  ::  nb
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      INTEGER  ::  m, n, l, k_d
      !------------------------------------------------------------ 


      k_d = SIZE( ele_b(1)%rnorms, 1)

      ALLOCATE( normal(k_d, SIZE(jd_jb)) )
      normal = 0

      DO m = 1, SIZE(ele_b)

         rnorms => ele_b(m)%rnorms
         nb     => ele_b(m)%nb
         rjp    => ele_b(m)%rjp
         
         DO l = 1, SIZE(rjp)

            DO n = 1, SIZE(nb)
               normal(:,nb(n)) = normal(:,nb(n))  +  rnorms(:,l) * rjp(l) 
            ENDDO
            
         ENDDO
 

      ENDDO


   END FUNCTION nodal_normals
   !============================================================ 


END MODULE fem_inlet
