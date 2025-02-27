!============================================================ 
!
!      Module: np_quadrature_tg 
!
! Description: Finite Element quadrature in node-pair 
!              formulation for Taylor-Galrkin terms
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

MODULE np_quadrature_tg


!============================================================ 
CONTAINS
!============================================================ 


   !============================================================ 
   FUNCTION np_quadr_TG_aGw_aGu(j_c_d, Tstiff_ij, aa, uu)      &
                               RESULT (pp)
   !============================================================ 


      ! (aGw , aGu)


      !------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:),   INTENT(IN)  ::  j_c_d
      REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  ::  Tstiff_ij
      REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  ::  aa
      REAL(KIND=8), DIMENSION(:),     INTENT(IN)  ::  uu

      REAL(KIND=8), DIMENSION(SIZE(Tstiff_ij,3))  ::  TGT
      REAL(KIND=8), DIMENSION(SIZE(uu))           ::  pp
      !------------------------------------------------------------


      TGT = np_quadr_TGT(j_c_d, aa, Tstiff_ij)
      
      pp = np_quadr_TGT_u(j_c_d, TGT, uu)
      

   END FUNCTION np_quadr_TG_aGw_aGu
   !============================================================ 



   !============================================================ 
   FUNCTION np_quadr_TGT(j_c_d, aa, Tstiff_ij)         &
                         RESULT (TGT)
   !============================================================ 


      ! a : Tstiff


      !------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:),   INTENT(IN)  ::  j_c_d
      REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  ::  Tstiff_ij
      REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  ::  aa

      REAL(KIND=8), DIMENSION(SIZE(Tstiff_ij,3))  ::  TGT
      !------------------------------------------------------------
      REAL(KIND=8), DIMENSION(SIZE(Tstiff_ij,1))  ::  a_ij
      INTEGER  ::  i, j, c, k
      !------------------------------------------------------------


      TGT = 0.d0
      
      DO c = 1, SIZE(TGT)
      

         i = j_c_d(1,c);  j = j_c_d(2,c)  

         a_ij = (aa(:,i) + aa(:,j))/2

         DO k =1, SIZE(Tstiff_ij,1)
            
            TGT(c) =  TGT(c)  &
                   +  a_ij(k) * SUM( a_ij * Tstiff_ij(:,k,c) ) 
         
         ENDDO


      ENDDO
   
   
   END FUNCTION np_quadr_TGT
   !============================================================ 



   !============================================================ 
   FUNCTION np_quadr_TGT_u(j_c_d, TGT, uu)                  &
                      RESULT (pp)
   !============================================================ 



      !------------------------------------------------------------
      IMPLICIT NONE 

      INTEGER,      DIMENSION(:,:),  INTENT(IN)  ::  j_c_d
      REAL(KIND=8), DIMENSION(:),    INTENT(IN)  ::  TGT
      REAL(KIND=8), DIMENSION(:),    INTENT(IN)  ::  uu

      REAL(KIND=8), DIMENSION(SIZE(uu))         ::  pp
      !------------------------------------------------------------
      REAL(KIND=8)  ::  pp_ij
      INTEGER  ::  i, j, c
      !------------------------------------------------------------


      pp = 0.d0
      
      DO c = 1, SIZE(TGT)
      
         i = j_c_d(1,c);  j = j_c_d(2,c)  

         pp_ij = TGT(c) * (uu(j) - uu(i))

         pp(i) = pp(i) + pp_ij
         pp(j) = pp(j) - pp_ij

      ENDDO
   
   
   END FUNCTION np_quadr_TGT_u
   !============================================================ 


END MODULE np_quadrature_tg
