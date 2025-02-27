!============================================================ 
!
!      Module: fem_quadrature_tg
!
! Description: quadrature functions for Taylor-Galrking terms 
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

MODULE fem_quadrature_tg


   !============================================================ 
   USE fem_ele_types
   USE dynamic_vector
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
!  [vGw_vGu]_i  ==  SUM_j,h ( N_h V_h . Grad(N_i)     , 
!                             N_h V_h . Grad(N_j) U_j )
!
!
!============================================================ 
!************************************************************ 
!============================================================ 


   !============================================================ 
   FUNCTION fem_quadr_vGw_vGu(ele_d, vv, uu)  RESULT (pp)
   !============================================================ 



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
      INTEGER      ::  m, l, k, n
      REAL(KIND=8) :: vdu_l, vdw_l_n
      REAL(KIND=8), DIMENSION(SIZE(vv,1)) :: v_l, du_l
      !------------------------------------------------------------ 


      pp = 0   

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            DO k = 1, SIZE(vv,1)
               v_l(k) = SUM(vv(k,nn) * ww(  :,l))
              du_l(k) = SUM(uu(  nn) * dw(k,:,l))
            ENDDO
            vdu_l = SUM(v_l * du_l) 

            DO n = 1, ele_d(m) % n_w
               vdw_l_n   = SUM(v_l * dw(:,n,l))
               pp(nn(n)) = pp(nn(n))  +  vdw_l_n * vdu_l * rjp(l)
            ENDDO

         ENDDO
      ENDDO
      

   END FUNCTION fem_quadr_vGw_vGu
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
!  [fvGw_fvGu]_i  ==  SUM_j,h ( fv . Grad(N_i) , 
!                               fv . Grad(N_j) U_j )
!
!
!============================================================ 
!************************************************************ 
!============================================================ 



   !============================================================ 
   FUNCTION fem_quadr_fvGw_fvGu(ele_d, fv, uu)  RESULT (pp)
   !============================================================ 


   !------------------------------------------------------------    
   !
   !  < fv.Gw, fv.Gu >   ===>   pp
   !
   !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain),    DIMENSION(:), INTENT(IN)  ::  ele_d
      TYPE(vector_func_gauss), DIMENSION(:), INTENT(IN)  ::  fv
      REAL(KIND=8),            DIMENSION(:), INTENT(IN)  ::  uu

      REAL(KIND=8),            DIMENSION(SIZE(uu))       ::  pp
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      INTEGER      ::  m, l, k, n
      REAL(KIND=8) :: vdu_l, vdw_l_n
      REAL(KIND=8), DIMENSION(SIZE(fv(1)%vv,1))  :: fv_l, du_l
      !------------------------------------------------------------ 


      pp = 0   

      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            fv_l = fv(m)%vv(:,l)
            
            DO k = 1, SIZE(fv_l,1)
              du_l(k) = SUM(uu(  nn) * dw(k,:,l))
            ENDDO
            
            vdu_l = SUM(fv_l * du_l) 

            DO n = 1, ele_d(m) % n_w
               vdw_l_n   = SUM(fv_l * dw(:,n,l))
               pp(nn(n)) = pp(nn(n))  +  vdw_l_n * vdu_l * rjp(l)
            ENDDO

         ENDDO
         
      ENDDO
      

   END FUNCTION fem_quadr_fvGw_fvGu
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
!
!  [wnvv_gbgu]_i  ==  SUM_j,h ( N_i norm . v, v . Grad(u) ) OUTLET 
!                  +  SUM_j,h ( N_i norm . v, v . Grad(b) ) INLET
!
!                     ! NOTE: Grad(b) is in input, 
!                             Grad(u) is SUM_j Grad(N_j) U_j 
!
!
!============================================================ 
!************************************************************ 
!============================================================ 

   !============================================================ 
   FUNCTION fem_quadr_wnvv_gbgu_inout(ele_b, ele_d_shell, vv_b, &
                                      Gb_b, uu_b) RESULT (pp_b)
   !============================================================ 


      ! < w(n.v) , v.Gu > with in/out conditions


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_boundary), DIMENSION(:),   INTENT(IN)  ::  ele_b
      TYPE(element_domain),   DIMENSION(:),   INTENT(IN)  ::  ele_d_shell
      REAL(KIND=8),           DIMENSION(:,:), INTENT(IN)  ::  vv_b
      REAL(KIND=8),           DIMENSION(:,:), INTENT(IN)  ::  Gb_b
      REAL(KIND=8),           DIMENSION(:),   INTENT(IN)  ::  uu_b

      REAL(KIND=8),           DIMENSION(SIZE(vv_b,2))       ::  pp_b
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  rnorms
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      REAL(KIND=8)  :: vn_l, vGu_l
      REAL(KIND=8), DIMENSION(SIZE(vv_b,1))  :: v_l, Gu_l
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

            ! Gu_inout            
            ! In
            IF ( vn_l < 1.e-12 ) THEN
            
               Gu_l =  Gb_b(:,m) 
            
            ! Out
            ELSE 
           
               DO k = 1, SIZE(rnorms, 1)
                  Gu_l(k) = SUM( ele_d_shell(m)%dw_C(k,:) &
                                 * uu_b(ele_d_shell(m)%nn) )
               ENDDO               
            
            ENDIF

            ! v.Gu
            vGu_l = SUM(v_l * Gu_l)

            pp_b(nn) = pp_b(nn)  +  ww(:,l) * vn_l * vGu_l * rjp(l) 

         ENDDO
      ENDDO


   END FUNCTION fem_quadr_wnvv_gbgu_inout
   !============================================================ 

!============================================================ 
!************************************************************ 
!============================================================ 
!
!  ======================
!  MATRICES IN CSR FORMAT 
!  ======================
!
!  Modified mass matrix
!  --------------------
!   
!    [vGw_vGw]_ij  ==  SUM_h( N_h V_h . Grad(N_i) , 
!                             N_h V_h . Grad(N_j) )          
!
!
!============================================================ 
!************************************************************ 
!============================================================ 


   !============================================================ 
   SUBROUTINE fem_quadr_vGw_vGw_CSR(ele_d, alpha, vv, AA)  
   !============================================================ 


   !  alpha < v.Gw, v.Gw >   ===>   Matrix in CSR format


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),   INTENT(IN)  ::  ele_d
      REAL(KIND=8),                         INTENT(IN)  ::  alpha
      REAL(KIND=8),         DIMENSION(:,:), INTENT(IN)  ::  vv

      TYPE(CSR_matrix) :: AA
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:),   POINTER  ::  ww
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      INTEGER       ::  m, l, k, i, i_, j, j_, idx
      REAL(KIND=8)  ::  TGij_l
      REAL(KIND=8), DIMENSION(SIZE(vv,1)) :: v_l
      !------------------------------------------------------------ 


      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         ww  => ele_d(m)%ww
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            DO k = 1, SIZE(vv,1)
               v_l(k) = SUM( vv(k,nn) * ww(:,l) ) 
            ENDDO

            DO i_ = 1, SIZE(nn);    i = nn(i_)
               DO j_ = 1, SIZE(nn); j = nn(j_)
                  
                  TGij_l = SUM( v_l * dw(:,i_,l) ) &
                         * SUM( v_l * dw(:,j_,l) )
                  
                  TGij_l = alpha * TGij_l * rjp(l)
                  
                  idx = csridx_srt(i,j,AA%i,AA%j)
                  AA%e(idx) = AA%e(idx)  +  TGij_l       
         
               ENDDO
            ENDDO
            
         ENDDO
      ENDDO

      
   END SUBROUTINE fem_quadr_vGw_vGw_CSR
   !============================================================


!============================================================ 
!************************************************************ 
!============================================================ 
!
!  ======================
!  MATRICES IN CSR FORMAT 
!  ======================
!
!  Modified mass matrix
!  --------------------
!   
!    [fvGw_fvGw]_ij  ==  ( fv . Grad(N_i) , fv . Grad(N_j) )          
!
!
!============================================================ 
!************************************************************ 
!============================================================ 


   !============================================================ 
   SUBROUTINE fem_quadr_fvGw_fvGw_CSR(ele_d, alpha, fv, AA)  
   !============================================================ 


   !  alpha < v.Gw, v.Gw >   ===>   Matrix in CSR format


      !------------------------------------------------------------ 
      IMPLICIT NONE 

      TYPE(element_domain), DIMENSION(:),    INTENT(IN)  ::  ele_d
      REAL(KIND=8),                          INTENT(IN)  ::  alpha
      TYPE(vector_func_gauss), DIMENSION(:), INTENT(IN)  ::  fv

      TYPE(CSR_matrix) :: AA
      !------------------------------------------------------------ 
      ! Elemental quantities
      INTEGER,      DIMENSION(:),     POINTER  ::  nn
      REAL(KIND=8), DIMENSION(:,:,:), POINTER  ::  dw
      REAL(KIND=8), DIMENSION(:),     POINTER  ::  rjp
      !------------------------------------------------------------ 
      INTEGER       ::  m, l, i, i_, j, j_, idx
      REAL(KIND=8)  ::  TGij_l
      REAL(KIND=8), DIMENSION(SIZE(fv(1)%vv,1)) :: fv_l
      !------------------------------------------------------------ 


      DO m = 1, SIZE(ele_d)

         nn  => ele_d(m)%nn
         dw  => ele_d(m)%dw
         rjp => ele_d(m)%rjp

         DO l = 1, ele_d(m) % l_G

            fv_l = fv(m)%vv(:,l)
            
            DO i_ = 1, SIZE(nn);    i = nn(i_)
               DO j_ = 1, SIZE(nn); j = nn(j_)
                  
                  TGij_l = SUM( fv_l * dw(:,i_,l) ) &
                         * SUM( fv_l * dw(:,j_,l) )
                  
                  TGij_l = alpha * TGij_l * rjp(l)
                  
                  idx = csridx_srt(i,j,AA%i,AA%j)
                  AA%e(idx) = AA%e(idx)  +  TGij_l       
         
               ENDDO
            ENDDO
            
         ENDDO
      ENDDO

      
   END SUBROUTINE fem_quadr_fvGw_fvGw_CSR
   !============================================================ 
    

END MODULE fem_quadrature_tg
