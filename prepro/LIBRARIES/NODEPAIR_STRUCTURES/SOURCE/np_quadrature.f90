!============================================================ 
!
!      Module: np_quadrature
!
! Description: Finite Element quadrature in node-pair 
!              formulation 
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
!  Available functions:
!
!     np_quadr_w_u
!     np_quadr_w_Dv
!     np_quadr_w_Gu
!     np_quadr_w_Rv
!     np_quadr_Gw_Gu
!     np_quadr_Gw_nGu
!
!============================================================ 

  MODULE np_quadrature

  CONTAINS


   !============================================================ 
   FUNCTION np_quadr_w_u(j_c, mass_ii, mass_ij, uu)            &
                 RESULT (pp)
   !============================================================ 


      ! (w , u)


      !------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  mass_ii, mass_ij
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8), DIMENSION(SIZE(uu))         ::  pp
      !------------------------------------------------------------
      INTEGER  ::  i, j, c
      !------------------------------------------------------------


      pp = 0.d0
      
      DO c = 1, SIZE(mass_ij)

         i = j_c(1,c);  j = j_c(2,c)

         pp(i) = pp(i) + mass_ij(c) * uu(j)
         pp(j) = pp(j) + mass_ij(c) * uu(i)

      ENDDO

      pp  =  pp  +  mass_ii * uu


   END FUNCTION np_quadr_w_u
   !============================================================ 



   !============================================================ 
   FUNCTION np_quadr_w_Dv(j_c_d, j_c_b, jd_jb,                 &
                          eta,   chi_b, xi_bp,                 &
                          vv, vv_b)                            &
                    RESULT (pp)
   !============================================================ 
 
 
      !  (w , Dv)

      !------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c_d, j_c_b
      INTEGER,      DIMENSION(:),   INTENT(IN)  ::  jd_jb
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  eta, chi_b, xi_bp
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  vv
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN), OPTIONAL  ::  vv_b

      REAL(KIND=8), DIMENSION(SIZE(vv,2))       ::  pp
      !------------------------------------------------------------
      REAL(KIND=8), DIMENSION(SIZE(vv,1)) :: vv_ij
      REAL(KIND=8)                        :: pp_ij
      INTEGER  ::  i, j, c, i_, j_
      !------------------------------------------------------------


      pp = 0.d0
      
      !------------------------------------------------------------
      ! Domain contribution
      ! -------------------
      DO c = 1, SIZE(eta,2)

         i = j_c_d(1,c);  j = j_c_d(2,c)

         vv_ij = 0.5 * ( vv(:,i)  +  vv(:,j) )
 
         pp_ij = SUM( vv_ij * eta(:,c) )

         pp(i) = pp(i) + pp_ij
         pp(j) = pp(j) - pp_ij

      ENDDO
      !------------------------------------------------------------
 
      
      !------------------------------------------------------------
      ! Boundary contribution
      ! ---------------------
      IF ( .NOT. PRESENT(vv_b) ) THEN
         
         !------------------------------------------------------------
         ! No imposed boundary conditions
         ! ------------------------------
         
         ! Node-pair contribution
         DO c = 1, SIZE(chi_b,2)

            
            i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
            i  = jd_jb(i_);   j  = jd_jb(j_)

            vv_ij = 0.5 * ( vv(:,j)  -  vv(:,i) )
            pp_ij = SUM( vv_ij * chi_b(:,c) )

            pp(i) = pp(i) + pp_ij
            pp(j) = pp(j) - pp_ij

         ENDDO

         ! Nodal contribution
         DO i_ = 1, SIZE(xi_bp,2)

            i = jd_jb(i_)
            pp(i) = pp(i) + SUM( vv(:,i) * xi_bp(:,i_) )

         ENDDO
         !------------------------------------------------------------
   
      ELSE
      
         !------------------------------------------------------------         
         ! Imposed boundary conditions
         ! ---------------------------
         
         ! Node-pair contribution
         DO c = 1, SIZE(chi_b,2)

            i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
            i  = jd_jb(i_);   j  = jd_jb(j_)

            vv_ij = 0.5 * ( vv_b(:,j_)  -  vv_b(:,i_) )
            
            pp_ij = SUM( vv_ij * chi_b(:,c) )

            pp(i) = pp(i) + pp_ij
            pp(j) = pp(j) - pp_ij

         ENDDO

         ! Nodal contribution
         DO i_ = 1, SIZE(xi_bp,2)

            i = jd_jb(i_)

            pp(i) = pp(i) + SUM( vv_b(:,i_) * xi_bp(:,i_) )

         ENDDO
         !------------------------------------------------------------

      ENDIF
      !------------------------------------------------------------


   END FUNCTION np_quadr_w_Dv
   !============================================================ 



   !============================================================ 
   FUNCTION np_quadr_w_Gu(j_c_d, j_c_b, jd_jb,                 &
                          eta,   chi_b, xi_bp,                 &
                          uu, uu_b)                            &
                    RESULT (pp)
   !============================================================ 
 

      !  (w , Gu)

      !------------------------------------------------------------
      use mp_interface
      
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c_d, j_c_b
      INTEGER,      DIMENSION(:),   INTENT(IN)  ::  jd_jb
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  eta, chi_b, xi_bp
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  uu
      REAL(KIND=8), DIMENSION(:),   INTENT(IN), OPTIONAL  ::  uu_b

      REAL(KIND=8), DIMENSION(SIZE(eta,1), SIZE(uu))  ::  pp
      !------------------------------------------------------------
      REAL(KIND=8)                         :: uu_ij
      REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: pp_ij
      INTEGER  ::  i, j, c, i_, j_
      !------------------------------------------------------------

      pp = 0.d0
      
      !------------------------------------------------------------
      ! Domain contribution
      ! -------------------
      DO c = 1, SIZE(eta,2)

         i = j_c_d(1,c);  j = j_c_d(2,c)

         uu_ij = 0.5 * ( uu(i)  +  uu(j) )
 
         pp_ij = uu_ij * eta(:,c) 

         pp(:,i) = pp(:,i) + pp_ij
         pp(:,j) = pp(:,j) - pp_ij

      ENDDO
      !------------------------------------------------------------
 
      
      !------------------------------------------------------------
      ! Boundary contribution
      ! ---------------------
      IF (.NOT. PRESENT(uu_b)) THEN
         
         !------------------------------------------------------------
         ! No imposed boundary conditions
         ! ------------------------------
         
         ! Node-pair contribution
         DO c = 1, SIZE(chi_b,2)
            
            i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
            i  = jd_jb(i_);   j  = jd_jb(j_)

            uu_ij = 0.5 * ( uu(j)  -  uu(i) )
            pp_ij = uu_ij * chi_b(:,c) 

            pp(:,i) = pp(:,i) + pp_ij
            pp(:,j) = pp(:,j) - pp_ij

         ENDDO

         ! Nodal contribution
         DO i_ = 1, SIZE(xi_bp,2)

            i = jd_jb(i_)
            pp(:,i) = pp(:,i)  +  uu(i) * xi_bp(:,i_) 

         ENDDO
         !------------------------------------------------------------
   
      ELSE
      
         !------------------------------------------------------------         
         ! Imposed boundary conditions
         ! ---------------------------
         
         ! Node-pair contribution
         DO c = 1, SIZE(chi_b,2)

            i_ = j_c_b(1,c);  j_ = j_c_b(2,c)
            i  = jd_jb(i_);   j  = jd_jb(j_)

            uu_ij = 0.5 * ( uu_b(j_)  -  uu_b(i_) )
            pp_ij = uu_ij * chi_b(:,c) 

            pp(:,i) = pp(:,i) + pp_ij
            pp(:,j) = pp(:,j) - pp_ij

         ENDDO

         ! Nodal contribution
         DO i_ = 1, SIZE(xi_bp,2)

            i = jd_jb(i_)

            pp(:,i) = pp(:,i)  +  uu_b(i_) * xi_bp(:,i_) 

         ENDDO
         !------------------------------------------------------------

      ENDIF
      !------------------------------------------------------------


   END FUNCTION np_quadr_w_Gu
   !============================================================ 



   !============================================================ 
   FUNCTION np_quadr_Gw_Gu(j_c_d, stiff_ij, uu)  RESULT (pp)
   !============================================================ 

      ! (Gw, Gu)

      !------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c_d
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  stiff_ij
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  uu

      REAL(KIND=8), DIMENSION(SIZE(uu))         ::  pp
      !------------------------------------------------------------
      REAL(KIND=8) :: pp_ij
      INTEGER  ::  i, j, c
      !------------------------------------------------------------


      pp = 0.d0
      
      DO c = 1, SIZE(stiff_ij)

         i = j_c_d(1,c);  j = j_c_d(2,c)

         pp_ij =  (uu(j) - uu(i))  *  stiff_ij(c) 

         pp(i) = pp(i) + pp_ij
         pp(j) = pp(j) - pp_ij

      ENDDO
      

   END FUNCTION np_quadr_Gw_Gu
   !============================================================ 



   !============================================================ 
   FUNCTION np_quadr_Gw_nGu(j_c_d, stiff_ij, uu, nu)           &
                    RESULT (pp)
   !============================================================ 

      ! (Gw, Gu)

      !------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c_d
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  stiff_ij
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  ::  uu, nu

      REAL(KIND=8), DIMENSION(SIZE(uu))         ::  pp
      !------------------------------------------------------------
      REAL(KIND=8) :: pp_ij, nu_ij
      INTEGER  ::  i, j, c
      !------------------------------------------------------------


      pp = 0.d0
      
      DO c = 1, SIZE(stiff_ij)

         i = j_c_d(1,c);  j = j_c_d(2,c)


         nu_ij = 0.5 * ( nu(i) + nu(j) )
         pp_ij = nu_ij * (uu(j) - uu(i))  *  stiff_ij(c) 

         pp(i) = pp(i) + pp_ij
         pp(j) = pp(j) - pp_ij

      ENDDO
      

   END FUNCTION np_quadr_Gw_nGu
   !============================================================ 

   
END MODULE np_quadrature
