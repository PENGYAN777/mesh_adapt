!============================================================ 
!
!      Module: fem_quadrature
!
! Description: quadrature functions for the finite element
!              method
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

MODULE fem_quadrature



   !============================================================ 
   USE fem_ele_types
   USE csr
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



!============================================================ 
!************************************************************ 
!============================================================ 
!
!  =========================
!  INTEGRALS OVER THE DOMAIN 
!  =========================
!
!  BOTH THE UNKNOWN AND THE FUNCTIONS ARE REINTERPOLATED
!
!  
!  Scalar unknwon u
!  ----------------
!
!    [w_u]_i  ==  SUM_j ( N_i , N_j U_j )
!
!  [Gw_Gu]_i  ==  SUM_j ( Grad(N_i) , Grad(N_j) U_j )
!
!
!  Vector unknown v
!  ----------------
!
!   [w_Dv]_i  ==  SUM_j ( N_i , Div(N_j) V_j )
!
!   [Dw_v]_i  ==  SUM_j ( Grad(N_i) , N_j V_j )
!
!
!  Threefold products
!  ------------------
!
!  [w_vGu]_i  ==  SUM_j,h ( N_i , N_h V_h . Grad(N_j) U_j )
!
!  [Gw_vu]_i  ==  SUM_j,h ( Grad(N_i) , N_h V_h  N_j U_j )
!
!============================================================ 
!************************************************************ 
!============================================================ 


   !============================================================ 
   FUNCTION fem_quadr_w_u(ele_d, uu)  RESULT (pp)
   !============================================================ 


      ! < w , u >


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:), INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:), INTENT(IN)  ::  uu

      REAL(KIND=8),         DIMENSION(SIZE(uu))       ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),   POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:), POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:),   POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8)  :: u_l
      INTEGER  ::  m, l
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            u_l = SUM( uu(nn) * ww(:,l) ) 

            pp(nn) = pp(nn)  +  ww(:,l) * u_l * rjp(l)
         
         ENDDO

      ENDDO


   END FUNCTION fem_quadr_w_u
   !============================================================ 




   !============================================================ 
   FUNCTION fem_quadr_w_Dv(ele_d, vv)  RESULT (pp)
   !============================================================ 
 

      ! < w , Dv >


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:,:), INTENT(IN)  ::  vv

      REAL(KIND=8),         DIMENSION(SIZE(vv,2))       ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8)  ::  dv_l
      INTEGER  ::  m, l, k
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            dv_l = 0.d0
            DO k = 1, SIZE(dw,1)
               dv_l = dv_l + SUM( vv(k,nn) * dw(k,:,l) ) 
            ENDDO
            
            pp(nn) = pp(nn)  +  ww(:,l) * dv_l * rjp(l)
         
         ENDDO

      ENDDO


   END FUNCTION fem_quadr_w_Dv
   !============================================================ 



   !============================================================ 
   FUNCTION fem_quadr_Gw_v(ele_d, vv)  RESULT (pp)
   !============================================================ 


      ! << Gw , v >>


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:,:), INTENT(IN)  ::  vv

      REAL(KIND=8),         DIMENSION(SIZE(vv,2))       ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(SIZE(vv,1))  :: v_l
      INTEGER  ::  m, l, k, n
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            DO k = 1,SIZE(vv,1)
               v_l(k) = SUM( vv(k,nn) * ww(:,l) )
            ENDDO

            DO n = 1, ele_d(m) % n_w
               pp(nn(n)) = pp(nn(n))  +  SUM( dw(:,n,l) * v_l ) * rjp(l)
            ENDDO
         
         ENDDO
      ENDDO


   END FUNCTION fem_quadr_Gw_v
   !============================================================ 



   !============================================================ 
   FUNCTION fem_quadr_Gw_Gu(ele_d, uu)  RESULT (pp)
   !============================================================ 


      ! << Gw , Gu >>


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8),         DIMENSION(SIZE(uu))         ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(SIZE(ele_d(1)%dw,1))  :: du_l
      INTEGER  ::  m, l, k, n
      !------------------------------------------------------------ 

   
      pp = 0.d0
      
      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            DO k = 1, SIZE(dw,1)
               du_l(k) = SUM( uu(nn) * dw(k,:,l) ) 
            ENDDO

            DO n = 1, ele_d(m) % n_w
               pp(nn(n)) =  pp(nn(n))  +  SUM( dw(:,n,l) * du_l ) * rjp(l)
            ENDDO
         
         ENDDO
      ENDDO


   END FUNCTION fem_quadr_Gw_Gu
   !============================================================ 



   !============================================================ 
   FUNCTION fem_quadr_w_vGu(ele_d, vv, uu)  RESULT (pp)
   !============================================================ 


      ! < w , v.Gu >


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:,:), INTENT(IN)  ::  vv
      REAL(KIND=8),         DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8),         DIMENSION(SIZE(uu))         ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(SIZE(vv,1))  ::  Gu_l, v_l
      REAL(KIND=8)  :: vGu_l
      INTEGER  ::  m, l, k
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            DO k = 1,SIZE(vv,1)
               v_l(k)  =  SUM( vv(k,nn) * ww(  :,l) )
               Gu_l(k) =  SUM( uu(  nn) * dw(k,:,l) ) 
            ENDDO

            vGu_l =  SUM( v_l * Gu_l )

            pp(nn) = pp(nn)  +  ww(:,l) * vGu_l * rjp(l)
                              
         ENDDO

      ENDDO


   END FUNCTION fem_quadr_w_vGu
   !============================================================ 



   !============================================================ 
   FUNCTION fem_quadr_Gw_vu(ele_d, vv, uu)  RESULT (pp)
   !============================================================ 


      ! << Gw , vu >>


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),         DIMENSION(:,:), INTENT(IN)  ::  vv
      REAL(KIND=8),         DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8),         DIMENSION(SIZE(uu))         ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(SIZE(vv,1))  :: v_l
      REAL(KIND=8) :: u_l
      INTEGER  ::  m, l, k, n
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, SIZE(rjp)

            DO k = 1,SIZE(vv,1)
               v_l(k) = SUM( vv(k,nn) * ww(:,l) ) 
            ENDDO

            u_l = SUM( uu(nn) * ww(:,l) )

            DO n = 1, SIZE(nn)
               pp(nn(n)) = pp(nn(n))  &
                         +  SUM( dw(:,n,l) * v_l ) * u_l * rjp(l)
            ENDDO
         
         ENDDO

      ENDDO


   END FUNCTION fem_quadr_Gw_vu
   !============================================================ 




!============================================================ 
!************************************************************ 
!============================================================ 
!
!  =========================
!  INTEGRALS OVER THE DOMAIN 
!  =========================
!
!  FUNCTIONS ARE KNOWN AT GAUSS POINTS
!
!  
!  Scalar function fs
!  -----------------
!
!    [w_fs]_i  ==  ( N_i , fs )
!
!
!  Vector unknown v
!  ----------------
!
!   [Dw_fv]_i  ==  SUM_j ( Grad(N_i) , fv )
!
!
!  Threefold products
!  ------------------
!
!  [w_fvGu]_i  ==  SUM_j ( N_i , fv . Grad(N_j) U_j )
!
!  [Gw_fvu]_i  ==  SUM_j ( Grad(N_i) ,fv  N_j U_j )
!
!
!============================================================ 
!************************************************************ 
!============================================================ 


   !============================================================ 
   FUNCTION fem_quadr_w_fs(ele_d, fs, N_dof)  RESULT (pp)
   !============================================================ 


      ! < w , u >


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:), INTENT(IN)  ::  ele_d
      TYPE(scalar_func_gauss), DIMENSION(:), INTENT(IN)  ::  fs
      INTEGER  ::  N_dof

      REAL(KIND=8),            DIMENSION(N_dof)            ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),   POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:), POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:),   POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8)  :: f_l
      INTEGER  ::  m, l
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            f_l = fs(m)%ss(l) 

            pp(nn) = pp(nn)  +  ww(:,l) * f_l * rjp(l)
         
         ENDDO

      ENDDO


   END FUNCTION fem_quadr_w_fs
   !============================================================ 




   !============================================================ 
   FUNCTION fem_quadr_Gw_fv(ele_d, fv, N_dof)  RESULT (pp)
   !============================================================ 


      ! << Gw , v >>


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:), INTENT(IN)  ::  ele_d
      TYPE(vector_func_gauss), DIMENSION(:), INTENT(IN)  ::  fv
      INTEGER  ::  N_dof

      REAL(KIND=8),            DIMENSION(N_dof)          ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(SIZE(fv(1)%vv,1))  :: fv_l
      INTEGER  ::  m, l, n
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            fv_l = fv(m)%vv(:,l)

            DO n = 1, ele_d(m) % n_w
               pp(nn(n)) = pp(nn(n))  +  SUM( dw(:,n,l) * fv_l ) * rjp(l)
            ENDDO
         
         ENDDO
      ENDDO


   END FUNCTION fem_quadr_Gw_fv
   !============================================================ 



   !============================================================ 
   FUNCTION fem_quadr_w_fvGu(ele_d, fv, uu)  RESULT (pp)
   !============================================================ 


      ! < w , v.Gu >


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:), INTENT(IN)  ::  ele_d
      TYPE(vector_func_gauss), DIMENSION(:), INTENT(IN)  ::  fv
      REAL(KIND=8),            DIMENSION(:), INTENT(IN)  ::  uu

      REAL(KIND=8),            DIMENSION(SIZE(uu))       ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(SIZE(fv(1)%vv,1))  ::  Gu_l, fv_l
      REAL(KIND=8)  :: vGu_l
      INTEGER  ::  m, l, k
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G


            fv_l = fv(m)%vv(:,l)
            
            DO k = 1,SIZE(fv_l,1)
               Gu_l(k) =  SUM( uu(nn) * dw(k,:,l) ) 
            ENDDO

            vGu_l =  SUM( fv_l * Gu_l )

            pp(nn) = pp(nn)  +  ww(:,l) * vGu_l * rjp(l)
                              
         ENDDO

      ENDDO


   END FUNCTION fem_quadr_w_fvGu
   !============================================================ 



   !============================================================ 
   FUNCTION fem_quadr_Gw_fvu(ele_d, fv, uu)  RESULT (pp)
   !============================================================ 


      ! << Gw , vu >>


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:), INTENT(IN)  ::  ele_d
      TYPE(vector_func_gauss), DIMENSION(:), INTENT(IN)  ::  fv
      REAL(KIND=8),            DIMENSION(:), INTENT(IN)  ::  uu

      REAL(KIND=8),            DIMENSION(SIZE(uu))       ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8), DIMENSION(SIZE(fv(1)%vv,1))  :: fv_l
      REAL(KIND=8) :: u_l
      INTEGER  ::  m, l, n
      !------------------------------------------------------------ 


      pp = 0.d0

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, SIZE(rjp)

            fv_l = fv(m)%vv(:,l)

            u_l = SUM( uu(nn) * ww(:,l) )

            DO n = 1, SIZE(nn)
               pp(nn(n)) = pp(nn(n))  &
                         +  SUM( dw(:,n,l) * fv_l ) * u_l * rjp(l)
            ENDDO
         
         ENDDO

      ENDDO


   END FUNCTION fem_quadr_Gw_fvu
   !============================================================ 




!============================================================ 
!************************************************************ 
!============================================================ 
!
!  ===========================
!  INTEGRALS OVER THE BOUNDARY 
!  ===========================
!
!  BOTH THE UNKNOWN AND THE FUNCTIONS ARE REINTERPOLATED
!
!  [wv_nbu]_i  ==  SUM_j,h ( N_i v , norm u ) OUTLET 
!               +  SUM_j,h ( N_i v , norm b ) INLET
!
!
!============================================================ 
!************************************************************ 
!============================================================ 


   !============================================================ 
   FUNCTION fem_quadr_wv_nbu_inout(ele_b, vv_b, bb_b, uu_b)  RESULT (pp_b)
   !============================================================ 


      ! < w , vu.n > with in/out conditions


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_boundary), DIMENSION(:),   INTENT(IN)  ::  ele_b
      REAL(KIND=8),           DIMENSION(:,:), INTENT(IN)  ::  vv_b
      REAL(KIND=8),           DIMENSION(:),   INTENT(IN)  ::  bb_b, uu_b

      REAL(KIND=8),           DIMENSION(SIZE(uu_b))       ::  pp_b
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  rnorms
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8)  :: vn_l, u_l
      REAL(KIND=8), DIMENSION(SIZE(vv_b,1))  :: v_l
      INTEGER  ::  m, l, k
      !------------------------------------------------------------ 


      pp_b = 0.d0

      DO m = 1, SIZE(ele_b)

         nn     => ele_b(m)%nn
         ww     => ele_b(m)%ww
         rnorms => ele_b(m)%rnorms
         rjp    => ele_b(m)%rjp

         DO l = 1, ele_b(m) % l_G

            ! v(x_l) 
            DO k = 1,SIZE(vv_b,1)
               v_l(k) = SUM( vv_b(k,nn) * ww(:,l) ) 
            ENDDO

            ! (v(x_l) . n_l) 
            vn_l = SUM( v_l * rnorms(:,l) ) 

            ! In
            IF ( vn_l < 0.d0 ) THEN
               u_l =  SUM( bb_b(nn) * ww(:,l) )
            ! Out
            ELSE 
               u_l =  SUM( uu_b(nn) * ww(:,l) )
            ENDIF

            pp_b(nn) = pp_b(nn)  +  ww(:,l) * vn_l * u_l * rjp(l) 

         ENDDO
      ENDDO


   END FUNCTION fem_quadr_wv_nbu_inout
   !============================================================ 



!============================================================ 
!************************************************************ 
!============================================================ 
!
!  ======================
!  MATRICES IN CSR FORMAT 
!  ======================
!
!  Mass matrix
!  -----------
!   
!    [w_w]_ij  ==  ( N_i , N_j )          
!
!
!  Stiffness matrix
!  ----------------
!
!  [Gw_Gw]_ij  ==  ( Grad(N_i), Grad(N_j) )    
!
!  Lumping 
!  -------
!
!     lump_i   ==  SUM_j [ M_ij ]              
!
!      where M_ij is the (i,j) element of a generic matrix M 
!
!============================================================ 
!************************************************************ 
!============================================================ 



   !============================================================ 
   SUBROUTINE fem_quadr_w_w_CSR(ele_d, alpha,  AA)  
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:), INTENT(IN)  ::  ele_d
      REAL(KIND=8),                       INTENT(IN)  ::  alpha

      TYPE(CSR_matrix) :: AA
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8)  :: Mij_l
      INTEGER  ::  m, l, i_, i, j_, j, idx
      !------------------------------------------------------------ 


      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            DO i_ = 1, SIZE(nn);     i = nn(i_)
               DO j_ = 1, SIZE(nn);  j = nn(j_)
            
                  Mij_l = alpha * ww(i_,l) * ww(j_,l) * rjp(l)

                  idx = csridx_srt(i,j,AA%i, AA%j)
                  AA%e(idx) = AA%e(idx) + Mij_l
         
               ENDDO
            ENDDO

         ENDDO

      ENDDO


   END SUBROUTINE fem_quadr_w_w_CSR
   !============================================================ 



   !============================================================ 
   SUBROUTINE fem_quadr_Gw_Gw_CSR(ele_d, alpha, AA)  
   !============================================================ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:), INTENT(IN)  ::  ele_d
      REAL(KIND=8),                       INTENT(IN)  ::  alpha

      TYPE(CSR_matrix) :: AA
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8)  :: Kij_l
      INTEGER  ::  m, l, i_, i, j_, j, idx
      !------------------------------------------------------------ 


      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            DO i_ = 1, SIZE(nn);     i = nn(i_)
               DO j_ = 1, SIZE(nn);  j = nn(j_)
            
                  Kij_l = alpha * SUM( dw(:,i_,l) * dw(:,j_,l) ) * rjp(l)

                  idx = csridx_srt(i,j,AA%i,AA%j)
                  AA%e(idx) = AA%e(idx) + Kij_l   
         
               ENDDO
            ENDDO

         ENDDO

      ENDDO


   END SUBROUTINE fem_quadr_Gw_Gw_CSR
   !============================================================ 


   !============================================================ 
   FUNCTION fem_quadr_lump_CSR(AA)  RESULT (ll)
   !============================================================ 

      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(CSR_matrix) :: AA

      REAL(KIND=8), DIMENSION(SIZE(AA%i)-1)       ::  ll
      !------------------------------------------------------------ 
      INTEGER  ::  i, j, m
      !------------------------------------------------------------ 


      ll = 0.d0

      DO i = 1, SIZE(ll)
         DO m = AA%i(i), AA%i(i+1) - 1  
            j = AA%j(m)

            ll(i) = ll(i) + AA%e(csridx_srt(i,j,AA%i,AA%j))
         
         ENDDO
      ENDDO


   END FUNCTION fem_quadr_lump_CSR
   !============================================================ 

END MODULE fem_quadrature
