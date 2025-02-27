!============================================================ 
!
!      Module: euler_equations
!
! Description: Euler equations definition; eigenstructure,
!              fluxes, pressure etc.
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

   MODULE euler_equations

   USE structures
   USE thermodynamics

   CONTAINS


   !============================================================ 
   FUNCTION flux__ww(ww) RESULT(ff)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww)-2)  ::  ff
      
      REAL(KIND=8), DIMENSION(2:SIZE(ww)-1)  ::  II 
      REAL(KIND=8)  ::  Pi 
      INTEGER  ::  p, i 


      p = SIZE(ww) 

      Pi = Pi__ww(ww)

      ! Mass conservation law ( 1 equation )
      ! ------------------------------------
      ff(1,:) = ww(2:p-1) 

      ! Momentum conservation law ( p-2 equations )
      ! -------------------------------------------
      DO i = 2, p-1
         II = 0;  II(i) = 1
         ff(i,:) = (ww(i)*ww(2:p-1))/ww(1) + II*Pi   
      ENDDO 

      ! Energy conservation law ( 1 equation )
      ! --------------------------------------
      ff(p, :) = (ww(2:p-1)/ww(1))*(ww(p) + Pi) 


   END FUNCTION flux__ww
   !============================================================ 



   !============================================================ 
   FUNCTION flux__ww_eos(ww, eos) RESULT(ff)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      TYPE(eos_type),              INTENT(IN)       ::  eos

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww)-2)  ::  ff
      
      REAL(KIND=8), DIMENSION(2:SIZE(ww)-1)  ::  II 
      REAL(KIND=8)  ::  Pi 
      INTEGER  ::  p, i 


      p = SIZE(ww) 

      Pi = eos%P

      ! Mass conservation law ( 1 equation )
      ! ------------------------------------
      ff(1,:) = ww(2:p-1) 

      ! Momentum conservation law ( p-2 equations )
      ! -------------------------------------------
      DO i = 2, p-1
         II = 0;  II(i) = 1
         ff(i,:) = (ww(i)*ww(2:p-1))/ww(1) + II*Pi   
      ENDDO 

      ! Energy conservation law ( 1 equation )
      ! --------------------------------------
      ff(p, :) = (ww(2:p-1)/ww(1))*(ww(p) + Pi) 


   END FUNCTION flux__ww_eos
   !============================================================ 



   !============================================================ 
   FUNCTION jacobian__ww(ww) RESULT(AA)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)                ::  ww

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww),SIZE(ww)-2)  ::  AA
      
      REAL(KIND=8), DIMENSION(SIZE(ww))  ::  dPi
      REAL(KIND=8)  ::  ht
      REAL(KIND=8), DIMENSION(2:SIZE(ww)-1)  ::  II
      INTEGER  ::  p, i, j, k


      p = SIZE(ww) 

      dPi = dPi_dww(ww) 
      ht = (ww(p) + Pi__ww(ww)) / ww(1) 


      ! -------------------------------------------------------
      ! Build the Jacobian matrix for every spatial dimension k
      ! -------------------------------------------------------
      DO k = 1, p-2   
      
         ! Position in ww of the momentum along the k-th direction          
         i = k+1  

         II  = 0;   II(i) = 1   

         ! -----------------------------------------------------
         ! First column: partial derivatives with respect to rho
         ! -----------------------------------------------------
         AA(1,     1, k) = 0
         AA(2:p-1, 1, k) = -(ww(i)*ww(2:p-1))/(ww(1)**2) + II*dPi(1)
         AA(p,     1, k) = (ww(i)/ww(1)) * (dPi(1) - ht)

         ! -------------------------------------------------
         ! Second to p-2 column(s): partial derivatives with 
         ! respect to the momentum vector
         ! -------------------------------------------------

         ! Mass claw
         ! ---------
         AA(1, 2:p-1, k) = II

         ! Momentum claw
         ! -------------
         AA(2:p-1, 2:p-1, k) = 0
         DO j = 2, p-1
            AA(j,j,k) = ww(i)/ww(1)
         ENDDO
         AA(i,2:p-1, k) = AA(i,2:p-1,k) + dPi(2:p-1)
         AA(2:p-1,i,k) = AA(2:p-1,i,k) + ww(2:p-1)/ww(1)

         ! Energy claw
         ! -----------
         AA(p, 2:p-1, k) = (ww(i)/ww(1))*dPi(2:p-1) + II*ht

         ! ---------------------------------------------------
         ! Last column: partial derivatives with respect to Et
         ! ---------------------------------------------------
         AA(1,     p, k) = 0
         AA(2:p-1, p, k) = II*dPi(p)
         AA(p,     p, k) = (ww(i)/ww(1)) * (1 + dPi(p))

      ENDDO


   END FUNCTION jacobian__ww
   !============================================================ 



   !============================================================ 
   FUNCTION jacobian__ww_eos(ww, eos) RESULT(AA)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)                ::  ww
      TYPE(eos_ext_type),          INTENT(IN)                ::  eos

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww),SIZE(ww)-2)  ::  AA
      
      REAL(KIND=8), DIMENSION(SIZE(ww))  ::  dPi
      REAL(KIND=8)  ::  ht
      REAL(KIND=8), DIMENSION(2:SIZE(ww)-1)  ::  II
      INTEGER  ::  p, i, j, k


      p = SIZE(ww) 

      dPi = dPi_dww_eos(ww, eos) 
      ht = (ww(p) + eos%P) / ww(1) 


      ! -------------------------------------------------------
      ! Build the Jacobian matrix for every spatial dimension k
      ! -------------------------------------------------------
      DO k = 1, p-2   
      
         ! Position in ww of the momentum along the k-th direction          
         i = k+1  

         II  = 0;   II(i) = 1   

         ! -----------------------------------------------------
         ! First column: partial derivatives with respect to rho
         ! -----------------------------------------------------
         AA(1,     1, k) = 0
         AA(2:p-1, 1, k) = -(ww(i)*ww(2:p-1))/(ww(1)**2) + II*dPi(1)
         AA(p,     1, k) = (ww(i)/ww(1)) * (dPi(1) - ht)

         ! -------------------------------------------------
         ! Second to p-2 column(s): partial derivatives with 
         ! respect to the momentum vector
         ! -------------------------------------------------

         ! Mass claw
         ! ---------
         AA(1, 2:p-1, k) = II

         ! Momentum claw
         ! -------------
         AA(2:p-1, 2:p-1, k) = 0
         DO j = 2, p-1
            AA(j,j,k) = ww(i)/ww(1)
         ENDDO
         AA(i,2:p-1, k) = AA(i,2:p-1,k) + dPi(2:p-1)
         AA(2:p-1,i,k) = AA(2:p-1,i,k) + ww(2:p-1)/ww(1)

         ! Energy claw
         ! -----------
         AA(p, 2:p-1, k) = (ww(i)/ww(1))*dPi(2:p-1) + II*ht

         ! ---------------------------------------------------
         ! Last column: partial derivatives with respect to Et
         ! ---------------------------------------------------
         AA(1,     p, k) = 0
         AA(2:p-1, p, k) = II*dPi(p)
         AA(p,     p, k) = (ww(i)/ww(1)) * (1 + dPi(p))

      ENDDO


   END FUNCTION jacobian__ww_eos
   !============================================================ 



   !============================================================ 
   FUNCTION jacobian__ww_nn(ww, nn) RESULT(AA)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)     ::  ww
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)     ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww))  ::  AA
      
      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww),SIZE(nn))  ::  JAC
      INTEGER  ::  k

      JAC = jacobian__ww(ww)
      
      AA = 0.d0
      DO k = 1, SIZE(ww) - 2
         AA = AA + JAC(:,:,k)*nn(k)
      ENDDO


   END FUNCTION jacobian__ww_nn
   !============================================================ 



   !============================================================ 
   FUNCTION jacobian__ww_eos_nn(ww, eos, nn) RESULT(AA)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)     ::  ww
      TYPE(eos_ext_type),          INTENT(IN)                ::  eos
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)     ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww))  ::  AA
      
      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww),SIZE(nn))  ::  JAC
      INTEGER  ::  k

      JAC = jacobian__ww_eos(ww, eos)
      
      AA = 0.d0
      DO k = 1, SIZE(ww) - 2
         AA = AA + JAC(:,:,k)*nn(k)
      ENDDO


   END FUNCTION jacobian__ww_eos_nn
   !============================================================ 



   !============================================================ 
   FUNCTION eigenvalues__ww_nn(ww, nn) RESULT(lambda)
   !============================================================ 


      ! Ordering of the eigenvalues
      !
      ! l(1) = u.w - c,   l(2:p-1) = u.w,   l(p) = u.w + c

      ! Works with |nn| != 1 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww))             ::  lambda
      
      INTEGER  ::  p              ! Number of unknowns
      REAL(KIND=8)  :: c          ! Speed of sound
      REAL(KIND=8)  :: mod_nn  ! Length of vector nn
      REAL(KIND=8)  :: u_om       ! Shorthand for u.w


      c = c__ww(ww)

      p = SIZE(ww)
      mod_nn = SQRT( SUM(nn*nn) )
      u_om = SUM( (ww(2:p-1)/ww(1)) * nn )

       
      lambda(1)     = u_om - c*mod_nn  
      lambda(2:p-1) = u_om 
      lambda(p)     = u_om + c*mod_nn  
      

   END FUNCTION eigenvalues__ww_nn
   !============================================================ 


 
   !============================================================ 
   FUNCTION right_eigenvectors__ww_nn(ww, nn) RESULT(RR)
   !============================================================ 


      ! WARNING: WORKS IN 1 AND 2 DIMENSION ONLY

      ! The vector nn MUST have length equal to one

      ! Eigenvectors ordering corresponds to the following 
      ! choice for the eigenvalues l:
      !
      ! l(1) = u.w - c,   l(2:p-1) = u.w,   l(p) = u.w + c


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww))    ::  RR
      
      INTEGER  ::  p, i 
      REAL(KIND=8), DIMENSION(SIZE(nn))  ::  nnn !, bnn in 3D ?
      REAL(KIND=8), DIMENSION(SIZE(ww))  ::  dPi
      REAL(KIND=8)  ::  u_om, u_nom 
      REAL(KIND=8)  ::  u, v, w
      REAL(KIND=8)  ::  ht, c, cc


      ! Useful shorthands
      ! -----------------
      p = SIZE(ww)                               ! Dimension of the system
      u_om = SUM( (ww(2:p-1)/ww(1)) * nn )    ! Shorthand for u.w

      ! Thermodynamic-dependent quantities
      ! ----------------------------------
      dPi = dPi_dww(ww) 
      ht = (ww(p) + Pi__ww(ww)) / ww(1) 
      c   = c__ww(ww)
      cc  = c*c

      ! --------------------------------------
      ! Right eigenvectors (stored by columns)
      ! --------------------------------------
      RR(1,    1) = 1                          ! Pressure wave 1
      RR(2:p-1,1) = ww(2:p-1)/ww(1) - c*nn  ! Genuinely non-linear
      RR(p,    1) = ht - c*u_om

      RR(1,    2) = 1                          ! Entropy wave
      RR(2:p-1,2) = ww(2:p-1)/ww(1)            ! Linearly degenerate
      RR(p,    2) = ht - cc/dPi(p)

      RR(1,    p) = 1                          ! Pressure wave 2
      RR(2:p-1,p) = ww(2:p-1)/ww(1) + c*nn  ! Genuinely non-linear
      RR(p,    p) = ht + c*u_om

      DO i = 3, p-1  ! RRs corresponding to shear wave(s) 
                     ! (0 in 1D, 1 in 2D and 2 in 3D)
 
         nnn(1) = -nn(2);  nnn(2) = nn(1)  ! 2D only
         u_nom   = SUM( (ww(2:p-1)/ww(1)) * nnn )

         RR(1,    i) = 0                       ! Shear wave(s)
         RR(2:p-1,i) = nnn                  ! Linearly degenarate
         RR(p,    i) = u_nom

      ENDDO

      RR = RR / (2*c)     ! Normalization

      ! Three dimensional case.  WARNING: to be optimised
      IF ( p .NE. 5) RETURN
      
      u = ww(2)/ww(1);  v = ww(3)/ww(1);  w = ww(4)/ww(1)
      
      RR(1, 2) = nn(1)
      RR(2, 2) = u * nn(1)           
      RR(3, 2) = v * nn(1)  +  nn(3)
      RR(4, 2) = w * nn(1)  -  nn(2)
      RR(5, 2) = (ht - cc/dPi(p)) * nn(1)  + v * nn(3)  -  w * nn(2) 

      RR(1, 3) = nn(2)
      RR(2, 3) = u * nn(2)  -  nn(3)          
      RR(3, 3) = v * nn(2)  
      RR(4, 3) = w * nn(2)  +  nn(1)
      RR(5, 3) = (ht - cc/dPi(p)) * nn(2)  - u * nn(3)  +  w * nn(1) 
      
      RR(1, 4) = nn(3)
      RR(2, 4) = u * nn(3)  +  nn(2)          
      RR(3, 4) = v * nn(3)  -  nn(1)
      RR(4, 4) = w * nn(3)  
      RR(5, 4) = (ht - cc/dPi(p)) * nn(3)  + u * nn(2)  -  v * nn(1) 

      RR(:,2:4) = RR(:,2:4) / (2*c)     ! Normalization


   END FUNCTION right_eigenvectors__ww_nn
   !============================================================ 


 
   !============================================================ 
   FUNCTION left_eigenvectors__ww_nn(ww, nn) RESULT(LL)
   !============================================================ 


      ! WARNING: WORKS IN 1 AND 2 DIMENSION ONLY

      ! Actually, LL is the TRANSPOSED matrix of the left eigenvectors

      ! The vector nn MUST have length equal to one (for consistency)

      ! Eigenvectors ordering corresponds to the following 
      ! choice for the eigenvalues l:
      !
      ! l(1) = u.w - c,   l(2:p-1) = u.w,   l(p) = u.w + c


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww))    ::  LL
      
      INTEGER  ::  p, i 
      REAL(KIND=8), DIMENSION(SIZE(nn))  ::  nnn !, bnn in 3D ?
      REAL(KIND=8), DIMENSION(SIZE(ww))  ::  dPi
      REAL(KIND=8)  ::  u, v, w
      REAL(KIND=8)  ::  u_om, u_nom !, u_bom in 3D?
      REAL(KIND=8)  ::  ht, c, cc


      ! Useful shorthands
      ! -----------------
      p = SIZE(ww)                               ! Dimension of the system
      u_om = SUM( (ww(2:p-1)/ww(1)) * nn )    ! Shorthand for u.w

      ! Thermodynamic-dependent quantities
      ! ----------------------------------
      dPi = dPi_dww(ww) 
      ht = (ww(p) + Pi__ww(ww)) / ww(1) 
      c   = c__ww(ww)
      cc  = c*c

      ! ----------------------------------
      ! Left eigenvectors (stored by rows)
      ! ----------------------------------
      LL(1,1)     = dPi(1) + u_om*c           ! Pressure wave 1
      LL(1,2:p-1) = dPi(2:p-1) - c*nn      ! Genuinely non-linear
      LL(1,p)     = dPi(p)

      LL(2,1)     = -2*(dPi(1) - cc)          ! Entropy wave
      LL(2,2:p-1) = -2*dPi(2:p-1)             ! Linearly degenerate
      LL(2,p)     = -2*dPi(p)

      LL(p,1)     = dPi(1) - u_om*c           ! Pressure wave 2
      LL(p,2:p-1) = dPi(2:p-1) + c*nn      ! Genuinely non-linear
      LL(p,p)     = dPi(p)


      DO i = 3, p-1  ! LLs corresponding to shear wave(s) 
                     ! (0 in 1D, 1 in 2D and 2 in 3D)

         nnn(1) = -nn(2);  nnn(2) = nn(1) ! 2D only
         u_nom   = SUM( (ww(2:p-1)/ww(1)) * nnn )

         LL(i,1)     = -2*cc*u_nom            ! Shear wave(s) 
         LL(i,2:p-1) =  2*cc*nnn           ! Linearly degenerate
         LL(i,p)     =  0

      ENDDO

      LL = LL / c     ! Normalization

      ! Three dimensional case.  WARNING: to be optimised
      IF ( p .NE. 5) RETURN

      u = ww(2)/ww(1);  v = ww(3)/ww(1);  w = ww(4)/ww(1)
      
      LL(2,1) = - (dPi(1) - cc) * nn(1) - cc*( v*nn(3) - w*nn(2) )        
      LL(2,2) = -  dPi(2)       * nn(1)             
      LL(2,3) = -  dPi(3)       * nn(1) + cc*nn(3)            
      LL(2,4) = -  dPi(4)       * nn(1) - cc*nn(2)         
      LL(2,5) = -  dPi(5)       * nn(1)   
      
      LL(3,1) = - (dPi(1) - cc) * nn(2) - cc*( w*nn(1) - u*nn(3) )        
      LL(3,2) = -  dPi(2)       * nn(2) - cc*nn(3)           
      LL(3,3) = -  dPi(3)       * nn(2)             
      LL(3,4) = -  dPi(4)       * nn(2) + cc*nn(1)         
      LL(3,5) = -  dPi(5)       * nn(2)   
      
      LL(4,1) = - (dPi(1) - cc) * nn(3) - cc*( u*nn(2) - v*nn(1) )        
      LL(4,2) = -  dPi(2)       * nn(3) + cc*nn(2)           
      LL(4,3) = -  dPi(3)       * nn(3) - cc*nn(1)            
      LL(4,4) = -  dPi(4)       * nn(3)         
      LL(4,5) = -  dPi(5)       * nn(3)   
      
      LL(2:4,:) = 2 * LL(2:4,:) / c     ! Normalization
      

   END FUNCTION left_eigenvectors__ww_nn
   !============================================================ 


 
   !============================================================ 
   FUNCTION eigenvalues__ww_eos_nn(ww, eos, nn) RESULT(lambda)
   !============================================================ 


      ! Ordering of the eigenvalues
      !
      ! l(1) = u.w - c,   l(2:p-1) = u.w,   l(p) = u.w + c

      ! Works with |nn| != 1 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      TYPE(eos_ext_type),          INTENT(IN)       ::  eos
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww))             ::  lambda
      
      INTEGER  ::  p              ! Number of unknowns
      REAL(KIND=8)  :: c          ! Speed of sound
      REAL(KIND=8)  :: mod_nn  ! Length of vector nn
      REAL(KIND=8)  :: u_om       ! Shorthand for u.w


      c = eos%c

      p = SIZE(ww)
      mod_nn = SQRT( SUM(nn*nn) )
      u_om = SUM( (ww(2:p-1)/ww(1)) * nn )

       
      lambda(1)     = u_om - c*mod_nn  
      lambda(2:p-1) = u_om 
      lambda(p)     = u_om + c*mod_nn  
      

   END FUNCTION eigenvalues__ww_eos_nn
   !============================================================ 


 
   !============================================================ 
   FUNCTION right_eigenvectors__ww_eos_nn(ww, eos, nn) RESULT(RR)
   !============================================================ 


      ! WARNING: WORKS IN 1 AND 2 DIMENSION ONLY

      ! The vector nn MUST have length equal to one

      ! Eigenvectors ordering corresponds to the following 
      ! choice for the eigenvalues l:
      !
      ! l(1) = u.w - c,   l(2:p-1) = u.w,   l(p) = u.w + c


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      TYPE(eos_ext_type),          INTENT(IN)       ::  eos
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww))    ::  RR
      
      INTEGER  ::  p, i 
      REAL(KIND=8), DIMENSION(SIZE(nn))  ::  nnn !, bnn in 3D ?
      REAL(KIND=8), DIMENSION(SIZE(ww))  ::  dPi
      REAL(KIND=8)  ::  u_om, u_nom 
      REAL(KIND=8)  ::  u, v, w
      REAL(KIND=8)  ::  ht, c, cc


      ! Useful shorthands
      ! -----------------
      p = SIZE(ww)                               ! Dimension of the system
      u_om = SUM( (ww(2:p-1)/ww(1)) * nn )       ! Shorthand for u.w

      ! Thermodynamic-dependent quantities
      ! ----------------------------------
      dPi = dPi_dww_eos(ww, eos) 
      ht = (ww(p) + eos%P) / ww(1) 
      c = eos%c
      cc  = c*c

      ! --------------------------------------
      ! Right eigenvectors (stored by columns)
      ! --------------------------------------
      RR(1,    1) = 1                          ! Pressure wave 1
      RR(2:p-1,1) = ww(2:p-1)/ww(1) - c*nn     ! Genuinely non-linear
      RR(p,    1) = ht - c*u_om

      RR(1,    2) = 1                          ! Entropy wave
      RR(2:p-1,2) = ww(2:p-1)/ww(1)            ! Linearly degenerate
      RR(p,    2) = ht - cc/dPi(p)

      RR(1,    p) = 1                          ! Pressure wave 2
      RR(2:p-1,p) = ww(2:p-1)/ww(1) + c*nn  ! Genuinely non-linear
      RR(p,    p) = ht + c*u_om

      DO i = 3, p-1  ! RRs corresponding to shear wave(s) 
                     ! (0 in 1D, 1 in 2D and 2 in 3D)
 
         nnn(1) = -nn(2);  nnn(2) = nn(1)  ! 2D only
         u_nom   = SUM( (ww(2:p-1)/ww(1)) * nnn )

         RR(1,    i) = 0                       ! Shear wave(s)
         RR(2:p-1,i) = nnn                  ! Linearly degenarate
         RR(p,    i) = u_nom

      ENDDO

      RR = RR / (2*c)     ! Normalization

      ! Three dimensional case.  WARNING: to be optimised
      IF ( p .NE. 5) RETURN
      
      u = ww(2)/ww(1);  v = ww(3)/ww(1);  w = ww(4)/ww(1)
      
      RR(1, 2) = nn(1)
      RR(2, 2) = u * nn(1)           
      RR(3, 2) = v * nn(1)  +  nn(3)
      RR(4, 2) = w * nn(1)  -  nn(2)
      RR(5, 2) = (ht - cc/dPi(p)) * nn(1)  + v * nn(3)  -  w * nn(2) 

      RR(1, 3) = nn(2)
      RR(2, 3) = u * nn(2)  -  nn(3)          
      RR(3, 3) = v * nn(2)  
      RR(4, 3) = w * nn(2)  +  nn(1)
      RR(5, 3) = (ht - cc/dPi(p)) * nn(2)  - u * nn(3)  +  w * nn(1) 
      
      RR(1, 4) = nn(3)
      RR(2, 4) = u * nn(3)  +  nn(2)          
      RR(3, 4) = v * nn(3)  -  nn(1)
      RR(4, 4) = w * nn(3)  
      RR(5, 4) = (ht - cc/dPi(p)) * nn(3)  + u * nn(2)  -  v * nn(1) 

      RR(:,2:4) = RR(:,2:4) / (2*c)     ! Normalization




   END FUNCTION right_eigenvectors__ww_eos_nn
   !============================================================ 


 
   !============================================================ 
   FUNCTION left_eigenvectors__ww_eos_nn(ww, eos, nn) RESULT(LL)
   !============================================================ 


      ! WARNING: WORKS IN 1 AND 2 DIMENSION ONLY

      ! Actually, LL is the TRANSPOSED matrix of the left eigenvectors

      ! The vector nn MUST have length equal to one (for consistency)

      ! Eigenvectors ordering corresponds to the following 
      ! choice for the eigenvalues l:
      !
      ! l(1) = u.w - c,   l(2:p-1) = u.w,   l(p) = u.w + c


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      TYPE(eos_ext_type),          INTENT(IN)       ::  eos
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8), DIMENSION(SIZE(ww),SIZE(ww))    ::  LL
      
      INTEGER  ::  p, i 
      REAL(KIND=8), DIMENSION(SIZE(nn))  ::  nnn !, bnn in 3D ?
      REAL(KIND=8), DIMENSION(SIZE(ww))  ::  dPi
      REAL(KIND=8)  ::  u_om, u_nom !, u_bom in 3D?
      REAL(KIND=8)  ::  u, v, w
      REAL(KIND=8)  ::  ht, c, cc


      ! Useful shorthands
      ! -----------------
      p = SIZE(ww)                               ! Dimension of the system
      u_om = SUM( (ww(2:p-1)/ww(1)) * nn )    ! Shorthand for u.w

      ! Thermodynamic-dependent quantities
      ! ----------------------------------
      dPi = dPi_dww_eos(ww, eos) 
      ht = (ww(p) + eos%P) / ww(1) 
      c = eos%c 
      cc  = c*c

      ! ----------------------------------
      ! Left eigenvectors (stored by rows)
      ! ----------------------------------
      LL(1,1)     = dPi(1) + u_om*c           ! Pressure wave 1
      LL(1,2:p-1) = dPi(2:p-1) - c*nn      ! Genuinely non-linear
      LL(1,p)     = dPi(p)

      LL(2,1)     = -2*(dPi(1) - cc)          ! Entropy wave
      LL(2,2:p-1) = -2*dPi(2:p-1)             ! Linearly degenerate
      LL(2,p)     = -2*dPi(p)

      LL(p,1)     = dPi(1) - u_om*c           ! Pressure wave 2
      LL(p,2:p-1) = dPi(2:p-1) + c*nn      ! Genuinely non-linear
      LL(p,p)     = dPi(p)


      DO i = 3, p-1  ! LLs corresponding to shear wave(s) 
                     ! (0 in 1D, 1 in 2D and 2 in 3D)

         nnn(1) = -nn(2);  nnn(2) = nn(1) ! 2D only
         u_nom   = SUM( (ww(2:p-1)/ww(1)) * nnn )

         LL(i,1)     = -2*cc*u_nom            ! Shear wave(s) 
         LL(i,2:p-1) =  2*cc*nnn           ! Linearly degenerate
         LL(i,p)     =  0

      ENDDO

      LL = LL / c     ! Normalization

      ! Three dimensional case.  WARNING: to be optimised
      IF ( p .NE. 5) RETURN

      u = ww(2)/ww(1);  v = ww(3)/ww(1);  w = ww(4)/ww(1)
      
      LL(2,1) = - (dPi(1) - cc) * nn(1) - cc*( v*nn(3) - w*nn(2) )        
      LL(2,2) = -  dPi(2)       * nn(1)             
      LL(2,3) = -  dPi(3)       * nn(1) + cc*nn(3)            
      LL(2,4) = -  dPi(4)       * nn(1) - cc*nn(2)         
      LL(2,5) = -  dPi(5)       * nn(1)   
      
      LL(3,1) = - (dPi(1) - cc) * nn(2) - cc*( w*nn(1) - u*nn(3) )        
      LL(3,2) = -  dPi(2)       * nn(2) - cc*nn(3)           
      LL(3,3) = -  dPi(3)       * nn(2)             
      LL(3,4) = -  dPi(4)       * nn(2) + cc*nn(1)         
      LL(3,5) = -  dPi(5)       * nn(2)   
      
      LL(4,1) = - (dPi(1) - cc) * nn(3) - cc*( u*nn(2) - v*nn(1) )        
      LL(4,2) = -  dPi(2)       * nn(3) + cc*nn(2)           
      LL(4,3) = -  dPi(3)       * nn(3) - cc*nn(1)            
      LL(4,4) = -  dPi(4)       * nn(3)         
      LL(4,5) = -  dPi(5)       * nn(3)   
      
      LL(2:4,:) = 2 * LL(2:4,:) / c     ! Normalization
      

   END FUNCTION left_eigenvectors__ww_eos_nn
   !============================================================ 


 
   !============================================================ 
   FUNCTION c__ww(ww) RESULT(c)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8)                                  ::  c  ! Speed of sound
      
      REAL(KIND=8)  :: Pi, rho

      
      Pi  = Pi__ww(ww) 
      rho = ww(1)
      
      c = c__P_r(Pi, rho)
      

   END FUNCTION c__ww
   !============================================================ 



   !============================================================ 
   FUNCTION G__ww(ww) RESULT(G)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8)                                  ::  G  ! Fundamental
                                                           ! derivative
      REAL(KIND=8)  :: Pi, rho

      
      Pi  = Pi__ww(ww) 
      rho = ww(1)
      
      G = G__P_r(Pi, rho)
      

   END FUNCTION G__ww
   !============================================================ 



   !============================================================ 
   FUNCTION s__ww(ww) RESULT(s)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN) :: ww

      REAL(KIND=8) :: s  ! Entropy (up to a constant)
      
      REAL(KIND=8) :: Pi, rho

      
      Pi  = Pi__ww(ww) 
      rho = ww(1)
      
      s = s__P_r(Pi, rho)
      

   END FUNCTION s__ww
   !============================================================ 



   !============================================================ 
   FUNCTION T__ww(ww) RESULT(T)
   !============================================================ 


      ! T = T( e,                        1/rho ) 
      ! T = T( Et/rho - 0.5*|m|^2/rho^2, 1/rho )


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8)                                  ::  T  ! temperature
      
      INTEGER  ::  p
      REAL(KIND=8)  :: e, rho
      REAL(KIND=8)  :: mod_u_2

      p = SIZE(ww)

      mod_u_2 = SUM((ww(2:p-1)/ww(1))**2)
      
      e = ww(p) / ww(1)  -  0.5 * mod_u_2 
      rho = ww(1)

      T   = T__e_v(e, 1/rho)
      

   END FUNCTION T__ww
   !============================================================ 



   !============================================================ 
   FUNCTION Pi__ww(ww) RESULT(Pi)
   !============================================================ 


      ! P  = P( e,                      rho ) 
      ! Pi = P( Et/rho - 0.5*|m|^2/rho, rho )


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8)                                  ::  Pi  ! Pressure
      
      INTEGER  ::  p
      REAL(KIND=8)  :: e, rho
      REAL(KIND=8)  :: mod_u_2

      p = SIZE(ww)

      mod_u_2 = SUM((ww(2:p-1)/ww(1))**2)
      
      e = ww(p) / ww(1)  -  0.5 * mod_u_2 
      rho = ww(1)

      Pi   = P__e_r(e, rho)
      

   END FUNCTION Pi__ww
   !============================================================ 



   !============================================================ 
   FUNCTION eos__ww(ww) RESULT(eos)
   !============================================================ 


      ! P  = P( e,                      rho ) 
      ! Pi = P( Et/rho - 0.5*|m|^2/rho, rho )


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      TYPE(eos_type)                                ::  eos
      
      INTEGER  ::  p
      REAL(KIND=8)  :: e, rho
      REAL(KIND=8)  :: mod_u_2

      p = SIZE(ww)

      mod_u_2 = SUM( (ww(2:p-1)/ww(1))**2 )
      
      e = ww(p) / ww(1)  -  0.5 * mod_u_2 
      rho = ww(1)

      eos = eos__e_v(e, 1/rho)
      

   END FUNCTION eos__ww
   !============================================================ 



   !============================================================ 
   FUNCTION dPi_dww(ww) RESULT(dPi)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8), DIMENSION(SIZE(ww))             ::  dPi 
            
      INTEGER  :: p
      REAL(KIND=8), DIMENSION(2)  ::  dP 
      REAL(KIND=8)  :: e, rho
      REAL(KIND=8)  :: mod_u_2
      

      p = SIZE(ww)
      mod_u_2 = SUM( (ww(2:p-1)/ww(1))**2 )

      e = ww(p)/ww(1) - 0.5 * mod_u_2
      rho = ww(1)
      
      dP = dP__de_dr(e, rho)
      
      dPi(1)     =  dP(2)  -  (1/rho)*(ww(p)/rho - mod_u_2) * dP(1) 
      dPi(2:p-1) = - ( ww(2:p-1)/(rho**2) ) * dP(1)
      dPi(p)     = (1/rho) * dP(1)
      

   END FUNCTION dPi_dww
   !============================================================ 


  
   !============================================================ 
   FUNCTION dPi_dww_eos(ww, eos) RESULT(dPi)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      TYPE(eos_ext_type),          INTENT(IN)       ::  eos

      REAL(KIND=8), DIMENSION(SIZE(ww))             ::  dPi 
            
      REAL(KIND=8), DIMENSION(2)  ::  dP 
      REAL(KIND=8)  ::  rho
      REAL(KIND=8)  ::  mod_u_2
      INTEGER  ::  p
      

      p = SIZE(ww)
      mod_u_2 = SUM( (ww(2:p-1)/ww(1))**2 )

      rho = ww(1)
      dP = eos%dP
      
      dPi(1)     =  dP(2)  +  0.5* mod_u_2 * dP(1) 
      dPi(2:p-1) = - ( ww(2:p-1)/rho ) * dP(1)
      dPi(p)     = dP(1)
      

   END FUNCTION dPi_dww_eos
   !============================================================ 


  
   !============================================================ 
   FUNCTION Pi__vv(vv) RESULT(Pi)
   !============================================================ 


      ! P  = P( h,              rho )
      ! Pi = P( ht - 0.5*|v|^2, rho )

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  vv

      REAL(KIND=8)                                  ::  Pi  ! Pressure
      
      INTEGER  ::  p
      REAL(KIND=8)  :: h, rho
      REAL(KIND=8)  :: mod_u_2


      p = SIZE(vv)
      mod_u_2 = SUM( vv(2:p-1)**2 ) 

      h = vv(p) - 0.5 * mod_u_2
      rho = vv(1)
      
      Pi   = P__h_r(h, rho)
      

   END FUNCTION Pi__vv
   !============================================================ 



   !============================================================ 
   FUNCTION Etot__rho_m_P (rho, m, P) RESULT(Etot)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8),                INTENT(IN)  ::  rho
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)  ::  m
      REAL(KIND=8),                INTENT(IN)  ::  P

      REAL(KIND=8)                             ::  Etot
      
      
      Etot = rho * (e__P_r(P, rho) + 0.5*SUM((m/rho)**2 ))


   END FUNCTION Etot__rho_m_P 
   !============================================================ 



   !============================================================ 
   FUNCTION Etot__rho_m_T (rho, m, T) RESULT(Etot)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8),                INTENT(IN)  ::  rho
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)  ::  m
      REAL(KIND=8),                INTENT(IN)  ::  T

      REAL(KIND=8)                             ::  Etot
      
      Etot = rho  *  (e__T_v(T, 1/rho) + 0.5*SUM( (m/rho)**2 ) )
      

   END FUNCTION Etot__rho_m_T 
   !============================================================ 



   !============================================================ 
   FUNCTION Etot__rho_m_c (rho, m, c) RESULT(Etot)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8),                INTENT(IN)  ::  rho
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)  ::  m
      REAL(KIND=8),                INTENT(IN)  ::  c

      REAL(KIND=8)                             ::  Etot
      
      REAL(KIND=8)                             ::  P

      P = P__c_r(c, rho)
      
      Etot = rho  *  (e__P_r(P, rho) + 0.5*SUM( (m/rho)**2 ) )
      

   END FUNCTION Etot__rho_m_c 
   !============================================================ 



   !============================================================ 
   FUNCTION mod_u__ww(ww) RESULT(mod_u)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8)                                  ::  mod_u ! |u|
      
      INTEGER  ::  p


      p = SIZE(ww) 

      mod_u = SQRT( SUM(ww(2:p-1)**2) ) / ww(1)
      

   END FUNCTION mod_u__ww
   !============================================================ 



   !============================================================ 
   FUNCTION Mach__ww(ww) RESULT(Ma)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww

      REAL(KIND=8)                                  ::  Ma  ! Mach number
      
      REAL(KIND=8)  ::  mod_u, c
      INTEGER  ::  p


      p = SIZE(ww) 

      mod_u = SQRT( SUM(ww(2:p-1)**2) ) / ww(1)
          c = c__ww(ww)

      Ma = mod_u / c 
      

   END FUNCTION Mach__ww
   !============================================================ 



   !============================================================ 
   FUNCTION Mach_n__ww(ww, nn) RESULT(Ma_n)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8)                                  ::  Ma_n  ! Mach number
      
      REAL(KIND=8)  ::  u_n, c
      INTEGER  ::  p


      p = SIZE(ww) 

      u_n = SUM( (ww(2:p-1)*nn)/ww(1) ) 
        c = c__ww(ww)

      Ma_n = u_n / c 
      

   END FUNCTION Mach_n__ww
   !============================================================ 

   !============================================================ 
   FUNCTION Mach__ww_eos_nn(ww, eos, nn) RESULT(Ma_n)
   !============================================================ 


      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  ww
      TYPE(eos_ext_type),          INTENT(IN)       ::  eos
      REAL(KIND=8), DIMENSION(:),  INTENT(IN)       ::  nn

      REAL(KIND=8)                                  ::  Ma_n  ! Mach number
      
      REAL(KIND=8)  ::  u_n, c
      INTEGER  ::  p


      p = SIZE(ww) 

      u_n = SUM( (ww(2:p-1)*nn)/ww(1) ) 
        c = eos%c

      Ma_n = u_n / c 
      

   END FUNCTION Mach__ww_eos_nn
   !============================================================ 


END MODULE euler_equations
