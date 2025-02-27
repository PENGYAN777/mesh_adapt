!============================================================ 
!
!      Module: fem_gauss_points
!
! Description: Procedures to perform the isoparametric
!              transformation on one-, two-, thre-dimensional 
!              P1 finite elements 
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


   ! Warning: the element structure can deal with element of
   ! order greater than one.  This is to be taken into account
   ! in the gauss points generation (reference element of order
   ! P2, P3 etc) and in the association of the node indeces
   ! with the degrees of freedom of each element.
   ! Here, only P1 elements are considered. 

   !============================================================ 
   ! NOTATION
   ! ========
   !
   !------------------------------------------------------------ 
   ! Coordinates and shape functions:
   ! --------------------------------
   !
   !   X   :  (vector) coordinate in the physical space
   !
   !   S   :  (vector) coordinate in the reference space
   !
   !   N   :  Shape function
   !
   ! Subscripts:
   ! -----------
   !
   !   k   : k-th component of a vector in the physical space
   !
   !   h   : h-th component of a vector in the reference space
   !
   !   n   : indices indicating the n-th degree of freedom 
   !
   !   l   : indices indicating the l-th Gauss point 
   !
   ! Examples:
   ! ---------
   !
   !   X_k_l : k-th coordinate of the l-th Gauss point 
   !           in the physical space
   !   S_h_l : k-th coordinate of the l-th Gauss point 
   !           in the reference space
   ! N_n(S_l): n-th shape function evaluated at the l-th 
   !           Gauss point
   !
   !------------------------------------------------------------ 
   ! Variables defined in the module:
   ! --------------------------------
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
   !
   !------------------------------------------------------------ 
   ! Remarks:
   ! --------
   !
   !     d         [  d (x,y)  ]      d    
   !  --------  =  [ --------- ]   --------
   !   d(r,s)      [  d (r,s)  ]    d(x,y) 
   !
   !
   !     d         [  d (x,y)  ]-1    d    
   !  --------  =  [ --------- ]   --------
   !   d(x,y)      [  d (r,s)  ]    d(x,y) 
   !
   !
   ! The isoparametric transormation is linear, so the
   ! following relations old:
   !
   !       X_l  =  SUM_n X_n N_n(S_l)             == rr_G
   !
   !      d X                 d N_n(S)                     [  d (x,y)  ]    
   !     -----  =  SUM_n X_n ----------|          == dr == [ --------- ]
   !      d S                   d S    |S = S_l            [  d (r,s)  ]
   !
   !
   !============================================================ 

MODULE fem_gauss_points


   !============================================================ 
   USE fem_ele_types
   USE fem_ref_elements
   USE dynamic_vector
   USE lin_algebra
   !============================================================ 


   !============================================================ 
   PUBLIC   ::  gen_gauss, gen_gauss_ele_dom, gen_gauss_ele_bou

!   PRIVATE  ::  gen_gauss_ele_dom, gen_gauss_ele_bou
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   SUBROUTINE gen_gauss(j_m_d, j_m_b, jd_jb, rr,      &
                        ele_d, ele_b, Nn_d, Nn_b, rr_n) 
   !============================================================ 



      !------------------------------------------------------------ 
      IMPLICIT NONE
      ! Input variables      
      TYPE( D_I_V ),            DIMENSION(:),   INTENT(IN)     ::  j_m_d
      TYPE( D_I_V ),            DIMENSION(:),   INTENT(IN)     ::  j_m_b
      INTEGER,                  DIMENSION(:),   INTENT(IN)     ::  jd_jb
      REAL(KIND=8),             DIMENSION(:,:), INTENT(IN)     ::  rr
      ! Output variables
      TYPE( element_domain ),   DIMENSION(:),   INTENT(INOUT)  ::  ele_d
      TYPE( element_boundary ), DIMENSION(:),   INTENT(INOUT)  ::  ele_b
      INTEGER,                                  INTENT(OUT)    ::  Nn_d, Nn_b
      REAL(KIND=8),             DIMENSION(:,:), POINTER        ::  rr_n                                                   
      !------------------------------------------------------------ 
      INTEGER  :: m, k_d
      TYPE (reference_element), DIMENSION(:), POINTER  ::  ref_ele_dom, &
                                                           ref_ele_bou
      !------------------------------------------------------------ 
      
      
      k_d = SIZE(rr(:,1))


      !------------------------------------------------------------ 
      ! Initialization of reference elements
      ! ------------------------------------

      CALL init_ref_elements(k_d, ref_ele_dom, ref_ele_bou)
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Compute the total number of degrees of freedom and 
      ! initialize the element/degree of freedom connectivity,  
      ! the element type and the vector rr_n 
      ! ------------------------------------------------------
      ! WARNING P1 element only :)))
      
      ALLOCATE( rr_n(SIZE(rr,1), SIZE(rr,2)) )
      rr_n = rr

      Nn_d = size_DIV(j_m_d,3)
      DO m = 1, SIZE(ele_d)
      
         ALLOCATE( ele_d(m)%nn(SIZE(j_m_d(m)%vec)) )
         ele_d(m)%nn = j_m_d(m)%vec

      ENDDO 
      
      Nn_b = size_DIV(j_m_b,3)
      DO m = 1, SIZE(ele_b)
      
         ALLOCATE( ele_b(m)%nn(SIZE(j_m_b(m)%vec)), &
                   ele_b(m)%nb(SIZE(j_m_b(m)%vec))  )

         ele_b(m)%nn = jd_jb(j_m_b(m)%vec)
         ele_b(m)%nb = j_m_b(m)%vec
         
      ENDDO 
      !------------------------------------------------------------ 
      
      

      !------------------------------------------------------------ 
      ! Domain element
      ! --------------
      DO m = 1, SIZE(ele_d)

         CALL gen_gauss_ele_dom( ele_d(m), rr_n(:,ele_d(m)%nn), &
                                 ref_ele_dom(ele_d(m)%ele_type)    )

      ENDDO 
      
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Boundary element
      ! ----------------
      IF ( SIZE(ele_b) > 0 ) THEN 

      DO m = 1, SIZE(ele_b)

         CALL gen_gauss_ele_bou( ele_b(m), rr_n(:,ele_b(m)%nn), &
                                 ref_ele_bou(ele_b(m)%ele_type) )

      ENDDO 
      
      ENDIF
      !------------------------------------------------------------ 


   END SUBROUTINE gen_gauss
   !============================================================ 



   !============================================================ 
   SUBROUTINE gen_gauss_ele_dom (ele_d, rr_nodes, r_e_d)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (element_domain),         INTENT(INOUT)  ::  ele_d
      REAL(KIND=8),  DIMENSION(:,:), INTENT(IN)     ::  rr_nodes
      TYPE (reference_element),      INTENT(IN)     ::  r_e_d
      !------------------------------------------------------------ 
      ! Parental element data      
      INTEGER ::  n_w,  l_G, k_d
      REAL(KIND=8), DIMENSION(:,:),   POINTER  :: w
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  :: d
      REAL(KIND=8), DIMENSION(:),     POINTER  :: p
      REAL(KIND=8), DIMENSION(:,:),   POINTER  :: d_C
      !------------------------------------------------------------ 
      ! Physical element data      
      INTEGER ::  l, n, k, h
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  :: rr_G
      REAL(KIND=8), DIMENSION(:),     ALLOCATABLE  :: r_C
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  :: dr, inv_dr
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE  :: dw
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  :: dw_C
      REAL(KIND=8), DIMENSION(:),     ALLOCATABLE  :: rjp
      !------------------------------------------------------------ 


      k_d = SIZE(rr_nodes, 1)

      !------------------------------------------------------------ 
      ! Retrieve the information on the parental element
      ! and allocates local structures
      ! ------------------------------------------------
      n_w = r_e_d % n_w;  l_G = r_e_d % l_G

        w => r_e_d % ww;  d => r_e_d % dd;  p => r_e_d % pp
      d_C => r_e_d % dd_C

      ALLOCATE ( rr_G(k_d,l_G), dr(k_d,k_d), inv_dr(k_d,k_d), & 
                 dw(k_d, n_w, l_G), rjp(l_G),                 &
                 r_C(k_d), dw_C(k_d, n_w)                     )        
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Isoparametric transformation in each Gauss point l
      ! --------------------------------------------------
      DO l = 1, l_G

         DO k = 1, k_d 
            DO h = 1, k_d
               dr(h, k) = SUM(rr_nodes(k,:) * d(h,:,l))
            ENDDO
         ENDDO

         inv_dr = inverse(dr)
         DO n = 1, n_w
            dw(:,n,l) = MATMUL( inv_dr, d(:,n,l) ) 
         ENDDO

         rjp(l) = determinant(dr) * p(l)

         DO k = 1, k_d 
            rr_G(k,l) = SUM( rr_nodes(k,:) * w(:,l) )
         ENDDO

      ENDDO
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Baricentrical quantities
      ! ------------------------
      DO k = 1, k_d 
         r_C(k) = SUM(rr_nodes(k,:))/SIZE(rr_nodes,2)
      ENDDO

      DO k = 1, k_d 
         DO h = 1, k_d
            dr(h, k) = SUM(rr_nodes(k,:) * d_C(h,:))
         ENDDO
      ENDDO

      inv_dr = inverse(dr)
      DO n = 1, n_w
         dw_C(:, n) = MATMUL( inv_dr, d_C(:,n) )
      ENDDO
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Store the information in the element structure
      ! ----------------------------------------------
      ele_d % n_w = n_w;  ele_d % l_G = l_G

      ALLOCATE( ele_d % dw(k_d, n_w, l_G), ele_d % rr_G(k_d,l_G), &
                ele_d % dw_C(k_d,n_w),     ele_d % rjp(l_G)       )

      ele_d % ww   => r_e_d % ww      
      ele_d % dw   = dw 
      ele_d % rr_G = rr_G
      ele_d % dw_C = dw_C;              
      ele_d % rjp  = rjp
      !------------------------------------------------------------ 

      DEALLOCATE ( rr_G, dr, inv_dr, dw, rjp, r_C, dw_C )        


   END SUBROUTINE gen_gauss_ele_dom
   !============================================================ 



   !============================================================ 
   SUBROUTINE gen_gauss_ele_bou(ele_b, rr_nodes, r_e_b)
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (element_boundary),                  INTENT(INOUT)  ::  ele_b
      REAL(KIND=8),             DIMENSION(:,:), INTENT(IN)     ::  rr_nodes
      TYPE (reference_element),                 INTENT(IN)     ::  r_e_b
      !------------------------------------------------------------ 
      ! Parental element data
      INTEGER ::  n_w, l_G, k_d
      REAL(KIND=8), DIMENSION(:,:),   POINTER  :: w
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  :: d
      REAL(KIND=8), DIMENSION(:),     POINTER  :: p
      !------------------------------------------------------------ 
      ! Physical element data
      INTEGER ::  l, k, h
      REAL(KIND=8) ::  rjac
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  :: rr_G
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  :: dr
      REAL(KIND=8), DIMENSION(:),     ALLOCATABLE  :: rjp
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  :: rnorms  
      !------------------------------------------------------------ 


      k_d = SIZE(rr_nodes, 1)
 
      !------------------------------------------------------------ 
      ! Retrieve the information on the parental element
      ! and allocates local structures
      ! ------------------------------------------------
      n_w = r_e_b % n_w;  l_G = r_e_b % l_G

      w => r_e_b % ww;  d => r_e_b % dd;  p => r_e_b % pp

      ALLOCATE ( rr_G(k_d,l_G), dr(k_d-1,k_d), &
                 rjp(l_G), rnorms(k_d,l_G)     )
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Isoparametric transformation in each Gauss point l
      ! --------------------------------------------------
      DO l = 1, l_G

         DO k = 1, k_d 
            DO h = 1, k_d - 1
               dr(h, k) = SUM( rr_nodes(k,:) * d(h,:,l) )
            ENDDO
         ENDDO
                      ! The minus sign is for Rebay's orientation
         rnorms(:,l) = - vec_prod(dr) 
         rjac = SQRT(SUM(rnorms(:,l)**2)) 
         rnorms(:,l) = rnorms(:,l)/rjac 
         rjp(l) = rjac * p(l)

         DO k = 1, k_d 
            rr_G(k,l) = SUM( rr_nodes(k,:) * w(:,l) )
         ENDDO


      ENDDO
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Store the information in the element strucure
      ! ---------------------------------------------
      ele_b % n_w = n_w;  ele_b % l_G = l_G
      
      ALLOCATE( ele_b % rr_G(k_d,l_G), ele_b % rjp(l_G), &
                ele_b % rnorms(k_d,l_G)                  )     

      ele_b % ww    => r_e_b % ww
      ele_b % rr_G   = rr_G
      ele_b % rjp    = rjp   
      ele_b % rnorms = rnorms 
      !------------------------------------------------------------ 
 
      DEALLOCATE ( rr_G, dr, rjp, rnorms )


   END SUBROUTINE gen_gauss_ele_bou
   !============================================================ 



END MODULE fem_gauss_points
