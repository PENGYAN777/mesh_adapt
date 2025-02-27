!============================================================ 
!
!      Module: fem_ref_elements
!
! Description: Definition of the reference element type for 
!              the finite element method.
!              One-, two-, thre-dimensional P1 elements.
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

MODULE fem_ref_elements

   
   !============================================================ 
   TYPE reference_element     
   !============================================================ 

      ! ------------------------------------------------------------ 
      ! Number of shape functions and Gauss points
      INTEGER(4)  ::  n_w,  l_G
      ! ------------------------------------------------------------ 
      ! k-th coordinate of l-th Gauss point
      !                       k,l
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  rr_G
      ! ------------------------------------------------------------ 
      ! Value of n-th shape function at l-th Gauss point
      !                       n,l
      REAL(KIND=8), DIMENSION(:,: ),  POINTER  ::  ww
      ! ------------------------------------------------------------ 
      ! Value of the k-th component of the gradient of the n-th 
      ! shape function at the l-th Gauss point
      !                       k,n,l 
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dd
      ! ------------------------------------------------------------ 
      ! Value of the weight of the l-th Gauss point
      !                       l
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  pp
      ! ------------------------------------------------------------ 
      ! k-th coordinate of the baricenter of the element
      !                       k
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rr_C
      ! ------------------------------------------------------------ 
      ! Value of the k-th component of the gradient of the n-th 
      ! shape function at the baricenter of the element
      !                       k,n
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  dd_C
      ! ------------------------------------------------------------ 

   END TYPE reference_element
   !============================================================ 



   !============================================================ 
!   TYPE (reference_element), DIMENSION(:), ALLOCATABLE  ::  ref_ele_dom
!   TYPE (reference_element), DIMENSION(:), ALLOCATABLE  ::  ref_ele_bou
   !============================================================ 


   !============================================================ 
!   PUBLIC   ::  ref_ele_dom, ref_ele_bou, init_ref_elements
   PUBLIC   ::  init_ref_elements
   
   PRIVATE  ::  init_ref_segment, &
                init_ref_triangle, init_ref_quadrilater, &
                init_ref_tetrahedron, init_ref_brick
   !============================================================ 


!============================================================ 
CONTAINS
!============================================================ 


   !============================================================ 
   SUBROUTINE init_ref_elements(k_d, ref_ele_dom, ref_ele_bou) 
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER  :: k_d   ! Dimension of the space
      TYPE (reference_element), DIMENSION(:), POINTER  ::  ref_ele_dom, &
                                                           ref_ele_bou
      !------------------------------------------------------------ 


      SELECT CASE(k_d)


      CASE (1) ! 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D 1D
      ! ============================================================ 

         ! ------------------------------------------
         ! Initialization of volume reference element
         ! 1D: domain elements are only segments(2)
         ! ------------------------------------------
         ALLOCATE( ref_ele_dom(1:1) )
         CALL  init_ref_segment ( ref_ele_dom(1) ) 

         ! -----------------------------------------------
         ! NO initialization of boundary reference element
         ! 1D: boundary elements are nodes 
         ! -----------------------------------------------


      CASE (2) ! 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D
      ! ============================================================ 

         ! ------------------------------------------
         ! Initialization of volume reference element
         ! 2D: domain elements are either triangle(3) 
         ! or quadrilater(4)
         ! ------------------------------------------
         ALLOCATE( ref_ele_dom(2:3) )
         CALL  init_ref_triangle    ( ref_ele_dom(2) )
         CALL  init_ref_quadrilater ( ref_ele_dom(3) )

         ! -------------------------------------------
         ! Initialization of surface reference element
         ! 2D: surface elements are only segments(2)
         ! -------------------------------------------
         ALLOCATE( ref_ele_bou(1:1) )
         CALL  init_ref_segment ( ref_ele_bou(1) ) 


      CASE (3) ! 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D 3D
      ! ============================================================ 

         ! -----------------------------------------------------
         ! Initialization of volume reference elements 
         ! 3D: domain elements are: tetrahedron (4), pyramid (5)
         ! tri-prism (6) or bric (7)
         ! -----------------------------------------------------
         ALLOCATE( ref_ele_dom(4:7) )
         CALL  init_ref_tetrahedron ( ref_ele_dom(4) )
         ! CALL  init_ref_pyramid    ( ref_ele_dom(5) )
          CALL  init_ref_prism    ( ref_ele_dom(5) )
         CALL  init_ref_prism      ( ref_ele_dom(6) )
         CALL  init_ref_brick       ( ref_ele_dom(7) )


         ! --------------------------------------------
         ! Initialization of surface reference elements
         ! 3D: boundary elements are either triangle(3) 
         ! or quadrilater(4)
         ! --------------------------------------------
         ALLOCATE( ref_ele_bou(2:3) )
         CALL  init_ref_triangle    ( ref_ele_bou(2) )
         CALL  init_ref_quadrilater ( ref_ele_bou(3) )


      CASE DEFAULT
      ! ============================================================ 
         WRITE (*,*) 'I dont know how to handle'
         WRITE (*,*) 'a grid with spatial dimension', k_d
         STOP 


      END SELECT
      
      
   END SUBROUTINE init_ref_elements 
   !============================================================ 


!============================================================ 
!
!       11         DDDD    ONE-DIMENSIONAL ELEMENT
!     1111         D   D   -----------------------
!       11    ===  D   D  
!       11         D   D
!     111111       DDDD
!
!============================================================ 


   !============================================================ 
   SUBROUTINE  init_ref_segment (r_e)   
   !============================================================ 


   !------------------------------------------------------------
   !  one-dimensional element with linear interpolation
   !  and two Gauss integration points
   !------------------------------------------------------------


      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (reference_element), INTENT(OUT)  ::  r_e
      !------------------------------------------------------------ 
!      INTEGER, PARAMETER  ::  n_w  = 2,  l_G  = 2
      INTEGER, PARAMETER  ::  n_w  = 2,  l_G  = 3
      REAL(KIND=8), DIMENSION(   n_w, l_G)  :: w
      REAL(KIND=8), DIMENSION(1, n_w, l_G)  :: d
      REAL(KIND=8), DIMENSION(l_G)          :: p
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(l_G) :: xx
      REAL(KIND=8), DIMENSION(1, n_w)  :: d_C
      REAL(KIND=8)  ::  xx_C
      !------------------------------------------------------------ 
      REAL(KIND=8) ::  zero= 0.d0, one = 1.d0,  two = 2.d0,  &
                       three = 3.d0, five = 5.d0, eight = 8.d0, nine = 9.d0
      !------------------------------------------------------------ 


      xx(1) = - SQRT(three/five)
      xx(2) = zero
      xx(3) = + SQRT(three/five)

         w(1, :) = (one - xx)/two
      d(1, 1, :) = - one/two

         w(2, :) = (xx + one)/two
      d(1, 2, :) = + one/two

            p(:) = five/nine
            p(2) = eight/nine

      xx_C = 0
      d_C(1,1) = - one/two
      d_C(1,2) = + one/two
      
      ! Store the information in the reference element structure
      r_e % n_w  =  n_w
      r_e % l_G  =  l_G 

      ALLOCATE( r_e % rr_G(1,1:l_G) )
      r_e % rr_G(1,:)  =  xx 

      ALLOCATE( r_e % ww(1:n_w,1:l_G) )
      r_e % ww = w 
      ALLOCATE( r_e % dd(1:1,1:n_w,1:l_G) )
      r_e % dd = d 
      ALLOCATE( r_e % pp(1:l_G) )
      r_e % pp =  p

      ALLOCATE( r_e % rr_C(1) )
      r_e % rr_C(1)  =  xx_C 
      ALLOCATE( r_e % dd_C(1:1,1:n_w) )
      r_e % dd_C = d_C 

   END SUBROUTINE  init_ref_segment
   !============================================================ 


!============================================================ 
!
!     222222       DDDD    BI-DIMENSIONAL ELEMENTS
!         22       D   D   -----------------------
!     2222    ===  D   D  
!     2            D   D
!     222222       DDDD
!
!============================================================ 

   !============================================================ 
   SUBROUTINE  init_ref_triangle (r_e)   
   !============================================================ 


   !------------------------------------------------------------
   !  triangular element with linear interpolation
   !  and three Gauss integration points
   !------------------------------------------------------------

      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (reference_element), INTENT(OUT)  ::  r_e
      !------------------------------------------------------------ 
      INTEGER, PARAMETER  ::  n_w  = 3,  l_G  = 3
      REAL(KIND=8), DIMENSION(   n_w, l_G)  :: w
      REAL(KIND=8), DIMENSION(2, n_w, l_G)  :: d
      REAL(KIND=8), DIMENSION(l_G)          :: p
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(l_G) :: xx, yy
      REAL(KIND=8), DIMENSION(2, n_w)  :: d_C
      REAL(KIND=8)  ::  xx_C, yy_C
      !------------------------------------------------------------ 
      REAL(KIND=8) :: zero = 0.d0,  one = 1.d0,  two = 2.d0,  &
                     three = 3.d0,  six = 6.d0 
      !------------------------------------------------------------ 


      xx(1) = one/six;  xx(2) = two/three;  xx(3) = one/six
      yy(1) = one/six;  yy(2) = one/six;    yy(3) = two/three

         w(1, :) = one - xx - yy
      d(1, 1, :) = - one
      d(2, 1, :) = - one

         w(2, :) = xx
      d(1, 2, :) = one
      d(2, 2, :) = zero

         w(3, :) = yy
      d(1, 3, :) = zero
      d(2, 3, :) = one

            p(:) = one/six


      xx_C = one/three; yy_C = one/three

      d_C(1,1) = - one 
      d_C(2,1) = - one
      d_C(1,2) =   one
      d_C(2,2) =   zero
      d_C(1,3) =   zero
      d_C(2,3) =   one
    
      
      ! Store the information in the reference element structure
      r_e % n_w  =  n_w
      r_e % l_G  =  l_G 

      ALLOCATE( r_e % rr_G(1:2,l_G) )
      r_e % rr_G(1,:)  =  xx 
      r_e % rr_G(2,:)  =  yy 

      ALLOCATE( r_e % ww(1:n_w,1:l_G) )
      r_e % ww  =  w 
      ALLOCATE( r_e % dd(1:2,1:n_w,1:l_G) )
      r_e % dd  =  d 
      ALLOCATE( r_e % pp(1:l_G) )
      r_e % pp  =  p

      ALLOCATE( r_e % rr_C(2) )
      r_e % rr_C(1)  =  xx_C
      r_e % rr_C(2)  =  yy_C
      ALLOCATE( r_e % dd_C(1:2,1:n_w) )
      r_e % dd_C = d_C 


   END SUBROUTINE init_ref_triangle
   !============================================================ 



   !============================================================ 
   SUBROUTINE  init_ref_quadrilater (r_e)  
   !============================================================ 


   !------------------------------------------------------------
   !  quadrilateral element with bilinear interpolation
   !  and four Gauss integration points
   !------------------------------------------------------------


      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (reference_element), INTENT(OUT)  ::  r_e
      !------------------------------------------------------------ 
      INTEGER, PARAMETER  ::  n_w = 4,  l_G = 9
!      INTEGER, PARAMETER  ::  n_w = 4,  l_G = 4
      REAL(KIND=8), DIMENSION(   n_w, l_G)  :: w
      REAL(KIND=8), DIMENSION(2, n_w, l_G)  :: d
      REAL(KIND=8), DIMENSION(l_G)          :: p
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(l_G) :: xx, yy
      REAL(KIND=8), DIMENSION(2, n_w)  :: d_C
      REAL(KIND=8)  ::  xx_C, yy_C
      !------------------------------------------------------------ 
      REAL(KIND=8) ::  zero= 0.d0, one = 1.d0,  &
                       three = 3.d0, four=4.d0, five = 5.d0, eight = 8.d0, &
                       nine = 9.d0, p1, p2, p3
      !------------------------------------------------------------ 

      xx(1) = - SQRT(three/five);  yy(1) = - SQRT(three/five)
      xx(2) = zero;                yy(2) = - SQRT(three/five)
      xx(3) = + SQRT(three/five);  yy(3) = - SQRT(three/five)
      xx(4) = - SQRT(three/five);  yy(4) = zero
      xx(5) = zero;                yy(5) = zero
      xx(6) = + SQRT(three/five);  yy(6) = zero
      xx(7) = - SQRT(three/five);  yy(7) = + SQRT(three/five)
      xx(8) = zero;                yy(8) = + SQRT(three/five)
      xx(9) = + SQRT(three/five);  yy(9) = + SQRT(three/five)

      p1 = five/nine; p2 = eight/nine; p3 = p1

      p(1) = p1*p1
      p(2) = p2*p1
      p(3) = p3*p1
      p(4) = p1*p2
      p(5) = p2*p2
      p(6) = p3*p2
      p(7) = p1*p3
      p(8) = p2*p3
      p(9) = p3*p3
      


         w(1, :) = (one - xx) * (one - yy) / four 
      d(1, 1, :) = - (one - yy) / four 
      d(2, 1, :) = - (one - xx) / four 

         w(2, :) = (one + xx) * (one - yy) / four 
      d(1, 2, :) =   (one - yy) / four 
      d(2, 2, :) = - (one + xx) / four

         w(3, :) = (one + xx) * (one + yy) / four 
      d(1, 3, :) =   (one + yy) / four
      d(2, 3, :) =   (one + xx) / four

         w(4, :) = (one - xx) * (one + yy) / four 
      d(1, 4, :) = - (one + yy) / four
      d(2, 4, :) =   (one - xx) / four


      xx_C = zero; yy_C = zero

      d_C(1,1) = - (one - yy_C) / four 
      d_C(2,1) = - (one - xx_C) / four 
      
      d_C(1,2) =   (one - yy_C) / four 
      d_C(2,2) = - (one - xx_C) / four
       
      d_C(1,3) =   (one + yy_C) / four
      d_C(2,3) =   (one + xx_C) / four 

      d_C(1,4) = - (one - yy_C) / four
      d_C(2,4) =   (one + xx_C) / four 

      r_e % n_w  =  n_w
      r_e % l_G  =  l_G

      ALLOCATE( r_e % rr_G(1:2,1:l_G) )
      r_e % rr_G(1,:)  =  xx 
      r_e % rr_G(2,:)  =  yy 

      ALLOCATE( r_e % ww(1:n_w,1:l_G) )
      r_e % ww  =  w 
      ALLOCATE( r_e % dd(1:2,1:n_w,1:l_G) )
      r_e % dd  =  d 
      ALLOCATE( r_e % pp(1:l_G) )
      r_e % pp  =  p

      ALLOCATE( r_e % rr_C(2) )
      r_e % rr_C(1)  =  xx_C
      r_e % rr_C(2)  =  yy_C
      ALLOCATE( r_e % dd_C(1:2,1:n_w) )
      r_e % dd_C = d_C 


   END SUBROUTINE  init_ref_quadrilater
   !============================================================ 


!============================================================ 
!
!     333333       DDDD    TRI-DIMENSIONAL ELEMENTS
!         33       D   D   ------------------------
!       3333  ===  D   D  
!         33       D   D
!     333333       DDDD
!
!============================================================ 


   !============================================================ 
   SUBROUTINE  init_ref_tetrahedron (r_e) 
   !============================================================ 
   

   !------------------------------------------------------------
   ! tetraedral element with linear interpolation
   ! and four Gauss integration points
   !------------------------------------------------------------


      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (reference_element), INTENT(OUT)  ::  r_e
      !------------------------------------------------------------ 
      INTEGER, PARAMETER  ::  n_w = 4,  l_G = 4
      REAL(KIND=8), DIMENSION(   n_w, l_G)  :: w
      REAL(KIND=8), DIMENSION(3, n_w, l_G)  :: d
      REAL(KIND=8), DIMENSION(l_G)          :: p
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(l_G) :: xx, yy, zz
      REAL(KIND=8), DIMENSION(3, n_w)  :: d_C
      REAL(KIND=8)  ::  xx_C, yy_C, zz_C
      REAL(KIND=8) :: zero = 0.d0,  one = 1.d0, two = 2.d0, five = 5.d0
      REAL(KIND=8) :: a, b
      !------------------------------------------------------------ 


      a = (five -      SQRT(five))/20.d0
      b = (five + 3.d0*SQRT(five))/20.d0

      xx    = a;      yy = a;      zz = a
      xx(3) = b;   yy(2) = b;   zz(1) = b


         w(1, :) =  one - xx - yy - zz
      d(1, 1, :) = - one
      d(2, 1, :) = - one
      d(3, 1, :) = - one

         w(2, :) = xx
      d(1, 2, :) = one
      d(2, 2, :) = zero
      d(3, 2, :) = zero

         w(3, :) = yy
      d(1, 3, :) = zero
      d(2, 3, :) = one
      d(3, 3, :) = zero

         w(4, :) = zz
      d(1, 4, :) = zero
      d(2, 4, :) = zero
      d(3, 4, :) = one

            p(:) = one/24.d0

            xx_C = one/two
            yy_C = one/two
            zz_C = one/two

       d_C(1, 1) = - one
       d_C(2, 1) = - one
       d_C(3, 1) = - one
 
       d_C(1, 2) = one
       d_C(2, 2) = zero
       d_C(3, 2) = zero
 
       d_C(1, 3) = zero
       d_C(2, 3) = one
       d_C(3, 3) = zero
 
       d_C(1, 4) = zero
       d_C(2, 4) = zero
       d_C(3, 4) = one

      r_e % n_w  =  n_w
      r_e % l_G  =  l_G

      ALLOCATE( r_e % rr_G(1:3,l_G) )
      r_e % rr_G(1,:)  =  xx 
      r_e % rr_G(2,:)  =  yy 
      r_e % rr_G(3,:)  =  zz 

      ALLOCATE( r_e % ww(1:n_w,1:l_G) )
      r_e % ww  =  w 
      ALLOCATE( r_e % dd(1:3,1:n_w,1:l_G) )
      r_e % dd  =  d 
      ALLOCATE( r_e % pp(1:l_G) ) 
      r_e % pp  =  p
      ALLOCATE( r_e % rr_C(3) )
      r_e % rr_C(1)  =  xx_C
      r_e % rr_C(2)  =  yy_C
      r_e % rr_C(3)  =  zz_C
      ALLOCATE( r_e % dd_C(1:3,1:n_w) )
      r_e % dd_C = d_C 


   END SUBROUTINE init_ref_tetrahedron
   !============================================================ 



   !============================================================ 
   SUBROUTINE  init_ref_prism (r_e) 
   !============================================================ 
   

   !------------------------------------------------------------
   ! prism element (triangular base) 
   !------------------------------------------------------------


      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (reference_element), INTENT(OUT)  ::  r_e
      !------------------------------------------------------------ 
      INTEGER, PARAMETER  ::  n_w = 6,  l_G = 6
      REAL(KIND=8), DIMENSION(   n_w, l_G)  :: w
      REAL(KIND=8), DIMENSION(3, n_w, l_G)  :: d
      REAL(KIND=8), DIMENSION(l_G)          :: p
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(l_G) :: xx, yy, zz
      REAL(KIND=8), DIMENSION(3, n_w)  :: d_C
      REAL(KIND=8)  ::  xx_C, yy_C, zz_C
      REAL(KIND=8) :: zero = 0.d0,  one = 1.d0,   two = 2.d0, &
                     three = 3.d0,  six = 6.d0, sq3
      !------------------------------------------------------------ 


      sq3 = SQRT(three)

      xx(1) = one/six;   xx(2) = two/three;  xx(3) = one/six
      yy(1) = one/six;   yy(2) = one/six;    yy(3) = two/three
      zz(1) = -one/sq3;  zz(2) = -one/sq3;   zz(3) = -one/sq3

      xx(4) = one/six;   xx(5) = two/three;  xx(6) = one/six
      yy(4) = one/six;   yy(5) = one/six;    yy(6) = two/three
      zz(4) =  one/sq3;  zz(5) =  one/sq3;   zz(6) =  one/sq3


         w(1, :) = (one - xx - yy) * (one - zz) / two
      d(1, 1, :) =  - (one - zz) / two
      d(2, 1, :) =  - (one - zz) / two
      d(3, 1, :) =  - (one - xx - yy) / two

         w(2, :) =  xx * (one - zz) / two
      d(1, 2, :) =  (one - zz) / two
      d(2, 2, :) =  zero
      d(3, 2, :) =  - xx / two

         w(3, :) =  yy * (one - zz) / two
      d(1, 3, :) =  zero 
      d(2, 3, :) =  (one - zz) / two
      d(3, 3, :) =  - yy / two

         w(4, :) = (one - xx - yy) * (one + zz) / two
      d(1, 4, :) = - (one + zz) / two
      d(2, 4, :) = - (one + zz) / two
      d(3, 4, :) =   (one - xx - yy) / two

         w(5, :) =  xx * (one + zz) / two
      d(1, 5, :) =  (one + zz) / two
      d(2, 5, :) =  zero
      d(3, 5, :) =  xx / two

         w(6, :) =  yy * (one + zz) / two
      d(1, 6, :) =  zero 
      d(2, 6, :) =  (one + zz) / two
      d(3, 6, :) =  yy / two

            p(:) = one/six

            xx_C = one/three
            yy_C = one/three
            zz_C = zero

      d_C(1, 1) =  - (one - zz_C) / two
      d_C(2, 1) =  - (one - zz_C) / two
      d_C(3, 1) =  - (one - xx_C - yy_C) / two

      d_C(1, 2) =  (one - zz_C) / two
      d_C(2, 2) =  zero
      d_C(3, 2) =  - xx_C / two

      d_C(1, 3) =  zero 
      d_C(2, 3) =  (one - zz_C) / two
      d_C(3, 3) =  - yy_C / two

      d_C(1, 4) = - (one + zz_C) / two
      d_C(2, 4) = - (one + zz_C) / two
      d_C(3, 4) =   (one - xx_C - yy_C) / two

      d_C(1, 5) =  (one + zz_C) / two
      d_C(2, 5) =  zero
      d_C(3, 5) =  xx_C / two

      d_C(1, 6) =  zero 
      d_C(2, 6) =  (one + zz_C) / two
      d_C(3, 6) =  yy_C / two

      r_e % n_w  =  n_w
      r_e % l_G  =  l_G

      ALLOCATE( r_e % rr_G(1:3,l_G) )
      r_e % rr_G(1,:)  =  xx 
      r_e % rr_G(2,:)  =  yy 
      r_e % rr_G(3,:)  =  zz 

      ALLOCATE( r_e % ww(1:n_w,1:l_G) )
      r_e % ww  =  w 
      ALLOCATE( r_e % dd(1:3,1:n_w,1:l_G) )
      r_e % dd  =  d 
      ALLOCATE( r_e % pp(1:l_G) ) 
      r_e % pp  =  p
      ALLOCATE( r_e % rr_C(3) )
      r_e % rr_C(1)  =  xx_C
      r_e % rr_C(2)  =  yy_C
      r_e % rr_C(3)  =  zz_C
      ALLOCATE( r_e % dd_C(1:3,1:n_w) )
      r_e % dd_C = d_C 


   END SUBROUTINE init_ref_prism
   !============================================================ 



   !============================================================ 
   SUBROUTINE  init_ref_brick (r_e) 
   !============================================================ 
   

   !------------------------------------------------------------
   ! brick element 
   !------------------------------------------------------------


      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (reference_element), INTENT(OUT)  ::  r_e
      !------------------------------------------------------------ 
      INTEGER, PARAMETER  ::  n_w = 8,  l_G = 8
      REAL(KIND=8), DIMENSION(   n_w, l_G)  :: w
      REAL(KIND=8), DIMENSION(3, n_w, l_G)  :: d
      REAL(KIND=8), DIMENSION(l_G)          :: p
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(l_G) :: xx, yy, zz
      REAL(KIND=8), DIMENSION(3, n_w)  :: d_C
      REAL(KIND=8)  ::  xx_C, yy_C, zz_C
      REAL(KIND=8) ::  zero = 0.d0, one = 1.d0,  three = 3.d0, &
                      eight = 8.d0, sq3
      !------------------------------------------------------------ 


      sq3 = SQRT(three)

      xx(1) = -one/sq3;  xx(2) =  one/sq3;  xx(3) =  one/sq3;  xx(4) = -one/sq3
      yy(1) = -one/sq3;  yy(2) = -one/sq3;  yy(3) =  one/sq3;  yy(4) =  one/sq3
      zz(1) = -one/sq3;  zz(2) = -one/sq3;  zz(3) = -one/sq3;  zz(4) = -one/sq3

      xx(5) = -one/sq3;  xx(6) =  one/sq3;  xx(7) =  one/sq3;  xx(8) = -one/sq3
      yy(5) = -one/sq3;  yy(6) = -one/sq3;  yy(7) =  one/sq3;  yy(8) =  one/sq3
      zz(5) =  one/sq3;  zz(6) =  one/sq3;  zz(7) =  one/sq3;  zz(8) =  one/sq3


         w(1, :) =   (one - xx) * (one - yy) * (one - zz) / eight
      d(1, 1, :) = -              (one - yy) * (one - zz) / eight
      d(2, 1, :) = - (one - xx) *              (one - zz) / eight
      d(3, 1, :) = - (one - xx) * (one - yy)              / eight

         w(2, :) =   (one + xx) * (one - yy) * (one - zz) / eight
      d(1, 2, :) =                (one - yy) * (one - zz) / eight
      d(2, 2, :) = - (one + xx) *              (one - zz) / eight
      d(3, 2, :) = - (one + xx) * (one - yy)              / eight

         w(3, :) =   (one + xx) * (one + yy) * (one - zz) / eight
      d(1, 3, :) =                (one + yy) * (one - zz) / eight
      d(2, 3, :) =   (one + xx) *              (one - zz) / eight
      d(3, 3, :) = - (one + xx) * (one + yy)              / eight

         w(4, :) =   (one - xx) * (one + yy) * (one - zz) / eight
      d(1, 4, :) = -              (one + yy) * (one - zz) / eight
      d(2, 4, :) =   (one - xx) *              (one - zz) / eight
      d(3, 4, :) = - (one - xx) * (one + yy)              / eight

         w(5, :) =   (one - xx) * (one - yy) * (one + zz) / eight
      d(1, 5, :) = -              (one - yy) * (one + zz) / eight
      d(2, 5, :) = - (one - xx) *              (one + zz) / eight
      d(3, 5, :) =   (one - xx) * (one - yy)              / eight

         w(6, :) =   (one + xx) * (one - yy) * (one + zz) / eight
      d(1, 6, :) =                (one - yy) * (one + zz) / eight
      d(2, 6, :) = - (one + xx) *              (one + zz) / eight
      d(3, 6, :) =   (one + xx) * (one - yy)              / eight

         w(7, :) =   (one + xx) * (one + yy) * (one + zz) / eight
      d(1, 7, :) =                (one + yy) * (one + zz) / eight
      d(2, 7, :) =   (one + xx) *              (one + zz) / eight
      d(3, 7, :) =   (one + xx) * (one + yy)              / eight

         w(8, :) =   (one - xx) * (one + yy) * (one + zz) / eight
      d(1, 8, :) = -              (one + yy) * (one + zz) / eight
      d(2, 8, :) =   (one - xx) *              (one + zz) / eight
      d(3, 8, :) =   (one - xx) * (one + yy)              / eight

            p(:) = one

            xx_C = zero
            yy_C = zero
            zz_C = zero

       d_C(1, 1) = - one/eight
       d_C(2, 1) = - one/eight
       d_C(3, 1) = - one/eight

       d_C(1, 2) =   one/eight
       d_C(2, 2) = - one/eight
       d_C(3, 2) = - one/eight

       d_C(1, 3) =   one/eight
       d_C(2, 3) =   one/eight
       d_C(3, 3) = - one/eight

       d_C(1, 4) = - one/eight
       d_C(2, 4) =   one/eight
       d_C(3, 4) = - one/eight

       d_C(1, 5) = - one/eight
       d_C(2, 5) = - one/eight
       d_C(3, 5) =   one/eight

       d_C(1, 6) =   one/eight
       d_C(2, 6) = - one/eight
       d_C(3, 6) =   one/eight

       d_C(1, 7) =   one/eight
       d_C(2, 7) =   one/eight
       d_C(3, 7) =   one/eight

       d_C(1, 8) = - one/eight
       d_C(2, 8) =   one/eight
       d_C(3, 8) =   one/eight

      r_e % n_w  =  n_w
      r_e % l_G  =  l_G

      ALLOCATE( r_e % rr_G(1:3,l_G) )
      r_e % rr_G(1,:)  =  xx 
      r_e % rr_G(2,:)  =  yy 
      r_e % rr_G(3,:)  =  zz 

      ALLOCATE( r_e % ww(1:n_w,1:l_G) )
      r_e % ww  =  w 
      ALLOCATE( r_e % dd(1:3,1:n_w,1:l_G) )
      r_e % dd  =  d 
      ALLOCATE( r_e % pp(1:l_G) ) 
      r_e % pp  =  p
      ALLOCATE( r_e % rr_C(3) )
      r_e % rr_C(1)  =  xx_C
      r_e % rr_C(2)  =  yy_C
      r_e % rr_C(3)  =  zz_C
      ALLOCATE( r_e % dd_C(1:3,1:n_w) )
      r_e % dd_C = d_C 


   END SUBROUTINE init_ref_brick
   !============================================================ 


END MODULE  fem_ref_elements


