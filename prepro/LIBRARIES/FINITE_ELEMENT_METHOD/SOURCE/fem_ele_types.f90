!============================================================ 
!
!      Module: fem_ele_types
!
! Description: Definition of the finite element type and of
!              (scalar and vector) functions evaluated at
!              Gauss points.
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
! Notation:
! 
!   X_k : coordinate of the k-th Gauss point 
!         in the physical space
!   S_k : coordinate of the k-th Gauss point 
!         in the reference space
!
!
!    w(n,l)  =  N_n(S_l)     [ == N_n(X_l) ]
!
!                d N_n(S)
!  d(k,n,l)  =  ----------|
!                  d S    |S = S_l
!
!                d X_k
!   dr(k,h)  =  -------
!                d S_h
!
!                d N_n(X) 
! dw(k,n,l)  =  ----------|
!                 d X_k   |X = X_l
!
!============================================================ 

!============================================================ 
!
! Remarks:
!
! To save memory, shape function are not effectively 
! associated to every element. Instead, the pointer ele % ww 
! is associated to a global variable ww, which is a component 
! of a vector whose length is given by the number of total 
! independent element types. See the Gauss points generation 
! module
!
!============================================================
 
MODULE fem_ele_types


   !============================================================ 
   TYPE element_domain
   !============================================================ 

      !------------------------------------------------------------ 
      ! Type of the element
      INTEGER  ::  ele_type
      !------------------------------------------------------------ 
      ! Number of shape functions and Gauss points
      INTEGER  ::  n_w, l_G
      !------------------------------------------------------------ 
      ! Index of the n-th degree in the global domain numbering 
      !                  n
      INTEGER, DIMENSION(:), POINTER  ::  nn 
      !------------------------------------------------------------ 
      ! k-th coordinate of l-th Gauss point
      !                       k,l
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  rr_G 
      !------------------------------------------------------------ 
      ! Value of n-th shape function evaluated at l-th Gauss point
      !                       n,l
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww 
      !------------------------------------------------------------ 
      ! Value of the k-th component of the gradient of the n-th 
      ! shape function evaluated at the l-th Gauss point
      !                       k,n,l 
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw 
      !------------------------------------------------------------ 
      ! Value of the product of the Jacobian for the weight of the 
      ! l-th Gauss point
      !                       l
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp 
      !------------------------------------------------------------ 
      ! Baricentrical quantities
      ! ------------------------
      ! k-th coordinate of baricenter of the element
      !                       k
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rr_C 
      !------------------------------------------------------------ 
      ! Value of n-th shape function evaluated at baricenter of the 
      ! element
      !                       n
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  ww_C 
      !------------------------------------------------------------ 
      ! Value of the k-th component of the gradient of the n-th 
      ! shape function evaluated at the baricenter of the element
      !                       k,n
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  dw_C 
      !------------------------------------------------------------ 

   END TYPE element_domain
   !============================================================ 



   !============================================================ 
   TYPE element_boundary
   !============================================================ 

      !------------------------------------------------------------ 
      ! Type of the element
      INTEGER  ::  ele_type
      !------------------------------------------------------------ 
      ! Number of shape function and Gauss points
      INTEGER  ::  n_w, l_G
      !------------------------------------------------------------ 
      ! Index of the n-th degree in the global domain numbering 
      !                  n
      INTEGER, DIMENSION(:), POINTER  ::  nn 
      !------------------------------------------------------------ 
      ! Index of the n-th degree in the global boundary numbering 
      !                  n
      INTEGER, DIMENSION(:), POINTER  ::  nb 
      !------------------------------------------------------------ 
      ! k-th coordinate of l-th Gauss point
      !                       k,l
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  rr_G 
      !------------------------------------------------------------ 
      ! Value of n-th shape function evaluated at l-th Gauss point
      !                       n,l
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww 
      !------------------------------------------------------------ 
      ! Value of the k-th component of the gradient of the n-th 
      ! shape function evaluated at the l-th Gauss point
      !                       k,n,l 
!      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw 
      !------------------------------------------------------------ 
      ! Value of product of the Jacobian for the weight for the 
      ! l-th Gauss point
      !                       l
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp 
      !------------------------------------------------------------ 
      ! Value of the k-th component of the outward-directed
      ! normal evaluated at the l-th Gauss point
      !                       k,l
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  rnorms 
      !------------------------------------------------------------ 
      ! Baricentrical quantities
      ! ------------------------
      ! k-th coordinate of the baricenter of the element
      !                       k
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rr_C 
      !------------------------------------------------------------ 
      ! Value of n-th shape function evaluated at baricenter of the
      ! element
      !                       n
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  ww_C 
      !------------------------------------------------------------ 
      ! Value of the k-th component of the gradient of the n-th 
      ! shape function evaluated at the baricenter of the element
      !                       k,n 
!      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  dw_C 
      !------------------------------------------------------------ 

   END TYPE element_boundary
   !============================================================ 



   !============================================================ 
   TYPE scalar_func_gauss
   !============================================================ 

      !------------------------------------------------------------ 
      ! Global index of the associated element
      INTEGER  ::  mm
      !------------------------------------------------------------ 
      ! Number of Gauss points
      INTEGER  ::  l_G
      !------------------------------------------------------------ 
      ! Value of the function evaluated at l-th Gauss point
      !                       l
      REAL(KIND=8), DIMENSION(:),   POINTER  ::  ss 
      !------------------------------------------------------------ 

   END TYPE scalar_func_gauss
   !============================================================ 



   !============================================================ 
   TYPE vector_func_gauss
   !============================================================ 

      !------------------------------------------------------------ 
      ! Global index of the associated element
      INTEGER  ::  mm
      !------------------------------------------------------------ 
      ! Number of Gauss points, number of components
      INTEGER  ::  l_G, n_p
      !------------------------------------------------------------ 
      ! Value of the p-th function evaluated at l-th Gauss point
      !                       p,l
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  vv 
      !------------------------------------------------------------ 

   END TYPE vector_func_gauss
   !============================================================ 



   !============================================================ 
   ! NOTATION
   ! ========
   !
   !   X_k : coordinate of the k-th Gauss point 
   !         in the physical space
   !   S_k : coordinate of the k-th Gauss point 
   !         in the reference space
   !
   !
   !    w(n,l)  =  N_n(S_l)     [ == N_n(X_l) ]
   !
   !                d N_n(S)
   !  d(k,n,l)  =  ----------|
   !                  d S    |S = S_l
   !
   !                d X_k
   !   dr(k,h)  =  -------
   !                d S_h
   !
   !                d N_n(X) 
   ! dw(k,n,l)  =  ----------|
   !                 d X_k   |X = X_l
   !
   !============================================================ 


   ! To save memory, shape function are not effectively 
   ! associated to every element. Instead, the pointer ele % ww 
   ! is associated to a global variable ww, which is a component 
   ! of a vector whose length is given by the number of total 
   ! independent element types. See the Gauss points generation 
   ! module


!============================================================ 
CONTAINS
!============================================================ 

   !============================================================ 
   FUNCTION u_to_f(ele_d, uu)  RESULT(ele_f)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),            DIMENSION(:),   INTENT(IN)  ::  uu

      TYPE(scalar_func_gauss), DIMENSION(SIZE(ele_d))      ::  ele_f

      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      !------------------------------------------------------------ 
      INTEGER       ::  m, l
      !------------------------------------------------------------ 
      

      DO m = 1, SIZE(ele_d)
     
         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
     
         ALLOCATE( ele_f(m)%ss(ele_d(m)%l_G))
         
         DO l = 1, ele_d(m) % l_G
 
            ele_f(m)%ss(l) = SUM( uu(nn) * ww(:,l) )
     
         ENDDO
         
      ENDDO



   END FUNCTION u_to_f
   !============================================================ 






END MODULE fem_ele_types
