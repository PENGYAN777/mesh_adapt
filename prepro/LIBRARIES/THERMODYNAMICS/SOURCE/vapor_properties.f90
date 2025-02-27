! !============================================================ 
! !
! !      Module: vapor_properties
! !
! ! Description: Procedure for the computation of the 
! !              saturation curve, G = 0 curve, isentropes
! !              and other loci pertaining the vapor-liquid
! !              phase transiction 
! !
! !      Author: Alberto Guardone
! !              Dipartimento di Ingegneria Aerospaziale
! !              Politecnico di Milano
! !              Via La Masa 34, 20156 Milano, ITALY
! !              e-mail: guardone@aero.polimi.it
! !
! !   Copyright: 1998-2003 Alberto Guardone
! !              See COPYING file for copyright notice
! !
! !============================================================ 
!  
 MODULE vapor_properties
! 
! 
!    !============================================================ 
!    USE thermodynamics
!    USE gas_properties, ONLY: gas
!    USE lin_algebra
!    USE derivatives
!    !============================================================ 
! 
! 
! 
!  !============================================================ 
!  CONTAINS
!  !============================================================ 
! 
!    
! 
!    !============================================================ 
!    SUBROUTINE saturation_curve( PP, vvg, vvl, TT)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::  PP
!       REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  vvg, vvl
!       REAL(KIND=8), DIMENSION(:), INTENT(OUT)    ::  TT
! 
!       REAL(KIND=8)  ::  P, vg, vl, Zc
! 
!       REAL(KIND=8), DIMENSION(2)  ::  FF, DD
!       REAL(KIND=8), DIMENSION(2, 2)  ::  JF
!       INTEGER   ::  i, it, N
!       REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-8
!       INTEGER ,     PARAMETER  ::  NEWTON_MAX_ITER = 100000
! 
!       
!       Zc = gas%Zc
!       N = SIZE(TT)
!       
!       DO i = 1, N
! 
!          P = PP(i) 
! 
!          ! Initial guess
!          IF (i == 1) THEN
!              vg = vvg(1);        vl = vvl(1)
!          ELSE
!              vg = 1.1*vvg(i-1);  vl = MAX(0.3d0, 0.9*vvl(i-1))
!          ENDIF
!          
!          DO it = 1, NEWTON_MAX_ITER + 1
!       
!             FF(1) = eq_chem_potential( (/ vg, vl /), (/ P /) )  
!             FF(2) = eq_temperature   ( (/ vg, vl /), (/ P /) )  
! 
!             IF ( MAXVAL(ABS(FF)) < NEWTON_EPS ) EXIT  
!          
!             JF(1,:) =  df_dxx_pp ( eq_chem_potential,                 &
!                                 (/ vg,   vl   /), (/ 1.d0, 1.d0 /), &
!                                 (/ P /)                             ) 
!             JF(2,:) =  df_dxx_pp ( eq_temperature,                       &
!                                 (/ vg,   vl   /), (/ 1.d0, 1.d0 /), &
!                                 (/ P /)                             )  
! 
!             DD = - MATMUL( inverse(JF), FF )
!                   
!             !   vg = vg +  DD(1);  vl = vl +  DD(2)
!             vg = vg + SIGN(1.d0,DD(1))* MIN(0.1d0,ABS(DD(1)))
!             vl = vl + SIGN(1.d0,DD(2))* MIN(0.1d0,ABS(DD(2)))
! 
!             IF(vg < 1) THEN   
!                IF (i == 1) THEN
!                   vg = vvg(1)
!                ELSE
!                   vg = 2*vvg(i-1)
!                ENDIF
!             ENDIF
!          
!             IF(vl > 1) THEN   
!                IF (i == 1) THEN
!                   vl = vvl(1)
!                ELSE
!                   vl = MAX(0.3d0, 0.5*vvl(i-1))
!                ENDIF
!             ENDIF
!          
! 
!         ENDDO
!       
!          IF ((it > NEWTON_MAX_ITER) .AND. &
!              (MAXVAL(ABS(FF)) > NEWTON_EPS) ) THEN
!             WRITE(*,*) ' Newton method failed to converge'
!             WRITE(*,*) ' in subroutine saturation_curve. STOP'
!             WRITE(*,*) FF
!             STOP
!          ENDIF
! 
!          vvg(i) = vg;  vvl(i) = vl
!          
!          TT(i) = T__P_v(P,vl)*gas%Zc  
!                
!       ENDDO
!       
!       
!    END SUBROUTINE saturation_curve
!    !============================================================ 
! 
!    
! 
!  
!    !============================================================ 
!    SUBROUTINE saturation_curve_T( TT, vvg, vvl, PP)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::  TT
!       REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  vvg, vvl
!       REAL(KIND=8), DIMENSION(:), INTENT(OUT)    ::  PP
! 
!       REAL(KIND=8)  ::  T, vg, vl, Zc
! 
!       REAL(KIND=8), DIMENSION(2)  ::  FF, DD
!       REAL(KIND=8), DIMENSION(2, 2)  ::  JF
!       INTEGER   ::  i, it, N
!       REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-6
!       INTEGER ,     PARAMETER  ::  NEWTON_MAX_ITER = 1000
! 
!       
!       Zc = gas%Zc
!       N = SIZE(TT)
!       
!       DO i = 1, N
! 
!          T = TT(i)/Zc 
! 
!          ! Initial guess
!          IF (i == 1) THEN
!              vg = vvg(1);        vl = vvl(1)
!          ELSE
!              vg = 1.1*vvg(i-1);  vl = MAX(0.3d0, 0.9*vvl(i-1))
!          ENDIF
!          
!          DO it = 1, NEWTON_MAX_ITER + 1
!       
!             FF(1) = eq_chem_potential_T( (/ vg, vl /), (/ T /) )  
!             FF(2) = eq_pressure        ( (/ vg, vl /), (/ T /) )  
! 
!             IF ( MAXVAL(ABS(FF)) < NEWTON_EPS ) EXIT  
!          
!             JF(1,:) =  df_dxx_pp ( eq_chem_potential_T,                 &
!                                 (/ vg,   vl   /), (/ 1.d0, 1.d0 /), &
!                                 (/ T /)                             ) 
!             JF(2,:) =  df_dxx_pp ( eq_pressure,                       &
!                                 (/ vg,   vl   /), (/ 1.d0, 1.d0 /), &
!                                 (/ T /)                             )  
! 
!             DD = - MATMUL( inverse(JF), FF )
!                   
!             !   vg = vg +  DD(1);  vl = vl +  DD(2)
!             vg = vg + SIGN(1.d0,DD(1))* MIN(0.01d0,ABS(DD(1)))
!             vl = vl + SIGN(1.d0,DD(2))* MIN(0.01d0,ABS(DD(2)))
! 
!             IF(vg < 1) THEN   
!                IF (i == 1) THEN
!                   vg = vvg(1)
!                ELSE
!                   vg = 2*vvg(i-1)
!                ENDIF
!             ENDIF
!          
!             IF(vl > 1) THEN   
!                IF (i == 1) THEN
!                   vl = vvl(1)
!                ELSE
!                   vl = MAX(0.3d0, 0.5*vvl(i-1))
!                ENDIF
!             ENDIF
!          
! 
!         ENDDO
!       
!          IF (it > NEWTON_MAX_ITER) THEN
!             WRITE(*,*) ' Newton method failed to converge'
!             WRITE(*,*) ' in subroutine saturation_curve. STOP'
!             STOP
!          ENDIF
! 
!          vvg(i) = vg;  vvl(i) = vl
!           PP(i) = P__T_v(T, vg) 
!                
!       ENDDO
!       
!       
!    END SUBROUTINE saturation_curve_T
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    SUBROUTINE iso_gamma(C, vv, PP)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8),               INTENT(IN)     ::  C
!       REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::  vv
!       REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  PP
! 
!       REAL(KIND=8)  ::  P, v
! 
!       REAL(KIND=8) ::  F, D
!       REAL(KIND=8) ::  JF
!       INTEGER   ::  i, it, N
!       REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-6
!       INTEGER ,     PARAMETER  ::  NEWTON_MAX_ITER = 1000
! 
!       
!       N = SIZE(vv)
!       
!       DO i = 1, N
!             
!          v = vv(i)   
!          IF( i == 1) THEN 
!             P = PP(1)
!          ELSE
!             P = PP(i-1)
!          ENDIF
!          
!          DO it = 1, NEWTON_MAX_ITER + 1
!      
!             F = eq_isogamma( P, (/ v, C /) )  
!             
!             IF ( ABS(F) < NEWTON_EPS ) EXIT        
!          
!          
!             JF =  df_dx_pp ( eq_isogamma, P , 1.d0, (/ v, C /) )  
! 
!             D  = - F/JF
!                   
!             P = P + SIGN(1.d0,D)* MIN(0.01d0,ABS(D))
! 
!             IF ( P < 0.d0 )  P = 0.1d0
! 
!          ENDDO
!       
!          IF (it > NEWTON_MAX_ITER) THEN
!             WRITE(*,*) ' Newton method failed to converge'
!             WRITE(*,*) ' in subroutine iso_gamma. STOP'
!             STOP
!          ENDIF
! 
!          PP(i) = P  
!                
!       ENDDO
!       
!       
!    END SUBROUTINE iso_gamma
!    !============================================================ 
!    
!    
!    
!    !============================================================ 
!    SUBROUTINE isentrope(C, vv, PP)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8),               INTENT(IN)     ::  C
!       REAL(KIND=8), DIMENSION(:), INTENT(IN)     ::  vv
!       REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  PP
! 
!       REAL(KIND=8)  ::  P, v
! 
!       REAL(KIND=8) ::  F, D
!       REAL(KIND=8) ::  JF
!       INTEGER   ::  i, it, N
!       REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-6
!       INTEGER ,     PARAMETER  ::  NEWTON_MAX_ITER = 1000
! 
!       
!       N = SIZE(vv)
!       
!       DO i = 1, N
!             
!          v = vv(i)   
!          IF( i == 1) THEN 
!             P = PP(1)
!          ELSE
!             P = PP(i-1)
!          ENDIF
!          
!          DO it = 1, NEWTON_MAX_ITER + 1
!      
!             F = eq_isentrope( P, (/ v, C /) )  
!             
!             IF ( ABS(F) < NEWTON_EPS ) EXIT        
!          
!          
!             JF =  df_dx_pp ( eq_isentrope, P , 1.d0, (/ v, C /) )  
! 
!             D  = - F/JF
!                   
!             ! P = P + SIGN(1.d0,D)* MIN(0.01d0,ABS(D))
!             P = P + D
! 
!             IF ( P < 0.d0 )  P = 0.1d0
! 
!          ENDDO
!       
!          IF (it > NEWTON_MAX_ITER) THEN
!             WRITE(*,*) ' Newton method failed to converge'
!             WRITE(*,*) ' in subroutine isentrope. STOP'
!             STOP
!          ENDIF
! 
!          PP(i) = P  
!                
!       ENDDO
!       
!       
!    END SUBROUTINE isentrope
!    !============================================================ 
! 
!  
! 
!    !============================================================ 
!    SUBROUTINE intersection_G_VD(C, vg, vl, P)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8),  INTENT(IN)      ::  C
!       REAL(KIND=8),  INTENT(INOUT)   ::  vg, vl, P
!       
!       REAL(KIND=8), DIMENSION(3)  ::  FF, DD
!       REAL(KIND=8), DIMENSION(3, 3)  ::  JF
!       REAL(KIND=8), DIMENSION(2)  ::  FF2, DD2
!       REAL(KIND=8), DIMENSION(2, 2)  ::  JF2
!       INTEGER   ::  it, it2
!       REAL(KIND=8), PARAMETER  ::  NEWTON_EPS      = 1.e-6
!       INTEGER ,     PARAMETER  ::  NEWTON_MAX_ITER = 1000
! 
!        
! 
!       DO it = 1, NEWTON_MAX_ITER + 1
! 
!          FF(1) = eq_chem_potential_vvP( (/ vg, vl, P /) )  
!          FF(2) = eq_temperature_vvP   ( (/ vg, vl, P /) )  
!          FF(3) = eq_isogamma_vvP      ( (/ vg, vl, P /), (/C/) )  
! 
!          IF ( MAXVAL(ABS(FF)) < NEWTON_EPS ) EXIT        
! 
! 
!          JF(1,:) =  df_dxx   ( eq_chem_potential_vvP,  &
!                                (/ vg, vl, P /), (/ 1.d0,1.d0,1.d0 /) ) 
!          JF(2,:) =  df_dxx   ( eq_temperature_vvP,     &
!                                (/ vg, vl, P /), (/ 1.d0,1.d0,1.d0 /) ) 
!          JF(3,:) =  df_dxx_pp ( eq_isogamma_vvP,        &   
!                                (/ vg, vl, P /), (/ 1.d0,1.d0,1.d0 /), (/C/) ) 
! 
!          DD = - MATMUL( inverse(JF), FF )
! 
!          vg = vg + SIGN(1.d0,DD(1))* MIN(0.01d0,ABS(DD(1)))
!          vl = vl + SIGN(1.d0,DD(2))* MIN(0.01d0,ABS(DD(2)))
!          P  = P  + SIGN(1.d0,DD(3))* MIN(0.01d0,ABS(DD(3)))
! 
!          DO it2 = 1, NEWTON_MAX_ITER + 1
!      
!             FF2(1) = eq_chem_potential( (/ vg, vl /), (/ P /) )  
!             FF2(2) = eq_temperature   ( (/ vg, vl /), (/ P /) )  
! 
!             IF ( MAXVAL(ABS(FF2)) < NEWTON_EPS ) EXIT  
!          
!             JF2(1,:) =  df_dxx_pp ( eq_chem_potential,                 &
!                                 (/ vg,   vl   /), (/ 1.d0, 1.d0 /), &
!                                 (/ P /)                             ) 
!             JF2(2,:) =  df_dxx_pp ( eq_temperature,                       &
!                                 (/ vg,   vl   /), (/ 1.d0, 1.d0 /), &
!                                 (/ P /)                             )  
! 
!             DD2 = - MATMUL( inverse(JF2), FF2 )
!                   
!             !   vg = vg +  DD(1);  vl = vl +  DD(2)
!             vg = vg + SIGN(1.d0,DD2(1))* MIN(0.01d0,ABS(DD2(1)))
!             vl = vl + SIGN(1.d0,DD2(2))* MIN(0.01d0,ABS(DD2(2)))
! 
!             IF(vg < 1)   vg = 2
!             IF(vl > 1)   vl = 0.4
! 
!          ENDDO
! 
!          IF(vg < 1)   vg = 2
!          IF(vl < 0.1) vl = 0.3
!          
!       ENDDO
! 
!       IF (it > NEWTON_MAX_ITER) THEN
!          WRITE(*,*) ' Newton method failed to converge'
!          WRITE(*,*) ' in subroutine intersection_G_VD. STOP'
!          STOP
!       ENDIF
! 
! 
!       
!    END SUBROUTINE intersection_G_VD
!    !============================================================ 
! 
!    
!    !============================================================ 
!    SUBROUTINE isotherm(T, vv, PP)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8),               INTENT(IN)     ::  T
!       REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  ::  vv, PP
! 
!       INTEGER   ::  i
! 
!       DO i = 1, SIZE(vv)
!       
!          PP(i) = P__T_v(T,vv(i))
!                
!       ENDDO
!       
!       
!    END SUBROUTINE isotherm
!    !============================================================ 
! 
! 
!    !============================================================ 
!    FUNCTION eq_isogamma( P, v_C ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8),               INTENT(IN)  ::  P  
!       REAL(KIND=8), DIMENSION(2), INTENT(IN)  ::  v_C  
!       REAL(KIND=8)                            ::  F     
! 
!       REAL(KIND=8)  ::  v, C     
! 
! 
!       v = v_C(1);  C = v_C(2)
! 
!       F =  G__P_r(P, 1/v) - C
! 
!  
!    END FUNCTION eq_isogamma
!    !============================================================ 
! 
!  
!   
!    !============================================================ 
!    FUNCTION eq_isentrope( P, v_C ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8),               INTENT(IN)  ::  P  
!       REAL(KIND=8), DIMENSION(2), INTENT(IN)  ::  v_C  
!       REAL(KIND=8)                            ::  F     
! 
!       REAL(KIND=8)  ::  v, C     
! 
! 
!       v = v_C(1);  C = v_C(2)
! 
!       F =  s__P_r(P, 1/v) - C
! 
!  
!    END FUNCTION eq_isentrope
!    !============================================================ 
! 
! 
!    !============================================================ 
!    FUNCTION eq_chem_potential( vv, PP ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(2), INTENT(IN)  ::  vv  
!       REAL(KIND=8), DIMENSION(1), INTENT(IN)  ::  PP  
!       REAL(KIND=8)              ::  F     
! 
!       REAL(KIND=8)              ::  vg, vl, P, Tg, Tl,eg, el     
! 
! 
!       vg = vv(1);  vl = vv(2);  P = PP(1)
!       
!       eg = e__P_r(P, 1/vg);  el = e__P_r(P, 1/vl)
!       Tg = T__P_v(P,   vg);  Tl = T__P_v(P,   vl)
! 
!       F =   eg + P*vg - Tg*s__T_v(Tg,vg) &
!         - ( el + P*vl - Tl*s__T_v(Tl,vl) )
! 
!  
!    END FUNCTION eq_chem_potential
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION eq_chem_potential_T( vv, TT ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(2), INTENT(IN)  ::  vv  
!       REAL(KIND=8), DIMENSION(1), INTENT(IN)  ::  TT  
!       REAL(KIND=8)              ::  F     
! 
!       REAL(KIND=8)              ::  vg, vl, T, Pg, Pl,eg, el     
! 
! 
!       vg = vv(1);  vl = vv(2);  T = TT(1)
!       
!       eg = e__T_v(T, vg);  el = e__T_v(T, vl)
!       Pg = P__T_v(T, vg);  Pl = P__T_v(T, vl)
! 
!       F =   eg + Pg*vg - T*s__T_v(T,vg) &
!         - ( el + Pl*vl - T*s__T_v(T,vl) )
! 
!  
!    END FUNCTION eq_chem_potential_T
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION eq_temperature( vv, PP ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(2), INTENT(IN)  ::  vv  
!       REAL(KIND=8), DIMENSION(1), INTENT(IN)  ::  PP  
!       REAL(KIND=8)              ::  F     
! 
!       REAL(KIND=8)              ::  vg, vl, P     
! 
! 
!       vg = vv(1);  vl = vv(2);  P = PP(1)
! 
!       F = T__P_v(P,vg) - T__P_v(P, vl)
! 
! 
!    END FUNCTION eq_temperature
!    !============================================================ 
! 
! 
!    
!    !============================================================ 
!    FUNCTION eq_pressure( vv, TT ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(2), INTENT(IN)  ::  vv  
!       REAL(KIND=8), DIMENSION(1), INTENT(IN)  ::  TT  
!       REAL(KIND=8)              ::  F     
! 
!       REAL(KIND=8)              ::  vg, vl, T    
! 
! 
!       vg = vv(1);  vl = vv(2);  T = TT(1)
! 
!       F = P__T_v(T,vg) - P__T_v(T, vl)
! 
! 
!    END FUNCTION eq_pressure
!    !============================================================ 
! 
! 
!    
!    !============================================================ 
!    FUNCTION eq_chem_potential_vvP( xx ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(3), INTENT(IN)  ::  xx  
!       REAL(KIND=8)              ::  F     
! 
!       REAL(KIND=8)              ::  vg, vl, P, Tg, Tl,eg, el     
! 
! 
!       vg = xx(1);  vl = xx(2);  P = xx(3)
!       
!       eg = e__P_r(P, 1/vg)
!       el = e__P_r(P, 1/vl)
!       Tg = T__e_v(eg,vg)
!       Tl = T__e_v(el,vl)
! 
!       F =   eg + P*vg - Tg*s__T_v(Tg,vg) &
!         - ( el + P*vl - Tl*s__T_v(Tl,vl) )
! 
!  
!    END FUNCTION eq_chem_potential_vvP
!    !============================================================ 
! 
! 
! 
!    !============================================================ 
!    FUNCTION eq_temperature_vvP( xx ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(3), INTENT(IN)  ::  xx  
!       REAL(KIND=8)              ::  F     
! 
!       REAL(KIND=8)              ::  vg, vl, P, Tg, Tl,eg, el     
! 
! 
!       
!       vg = xx(1);  vl = xx(2);  P = xx(3)
!       eg = e__P_r(P, 1/vg)
!       el = e__P_r(P, 1/vl)
!       Tg = T__e_v(eg,vg)
!       Tl = T__e_v(el,vl)
! 
!       F =   Tg - Tl
! 
! 
!    END FUNCTION eq_temperature_vvP
!    !============================================================ 
!    
!    
!    !============================================================ 
!    FUNCTION eq_isogamma_vvP( xx, C ) RESULT (F)
!    !============================================================ 
! 
! 
!       IMPLICIT NONE 
! 
!       REAL(KIND=8), DIMENSION(3), INTENT(IN)  ::  xx 
!       REAL(KIND=8), DIMENSION(1), INTENT(IN)  ::  C 
!       REAL(KIND=8)                            ::  F     
! 
!       REAL(KIND=8)  ::  P, v    
! 
! 
!       v = xx(1);  P = xx(3)
! 
!       F =  G__P_r(P, 1/v) - C(1)
! 
!  
!    END FUNCTION eq_isogamma_vvP
!    !============================================================ 
!    
! 
 END MODULE vapor_properties
