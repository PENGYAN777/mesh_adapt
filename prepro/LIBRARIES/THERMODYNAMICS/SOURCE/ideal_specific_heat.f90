!=========================================================================================
!
!      Module: ideal_specific_heat
!
! Description: Models for the specific heat cv_oo(T) at constant volume in the 
!              ideal gas limit.
!
!              The parameters read in procedure 'read_param_ideal_specific_heat'
!              are assumed to be either DIMENSIONAL OR DIMENSIONLESS with respect to
!              the CV_REF* according to the value of the logical
!              flag coeff_adim:
!
!                .TRUE.  = reading dimensionless data
!                .FALSE. = reading dimensional data
!
!              * NOTE: Before (Alberto) the values where nondimensional with respect
!                      to the Gas constant Rgas. Now (Marco) these are made dimensionless
!                      with respect to CV_REF. Note that CV_REF is choosen equal to 
!                      the reference gas constant RGAS_REF. If the latter is choosen equal
!                      to the Gas Constant Rgas the two approaches give the same result
!
!              The procedure 'init_ideal_specific_heat' makes these read parameters
!              dimensional if 'coeff_adim = .TRUE.'. This means that each function
!              in the following works with dimensional values of the read parameters.
!
!
!              Implemented curve fitting models are as follows:
!
!              1) Constant cv (polytropic)
!
!                  Parameter: cv_const (real)
!
!                    cv_oo(T) = cv_0
!
!                  Input field reads:
!                    cv_0
!
!
!              2) Power law for dense gases (power of the reduced temperature T/Tc)
!
!                  Parameters: cv_oo(Tc) (real), np (real)
!
!                    cv_oo(T) = cv_oo(Tc) (T/Tc)**np
!
!                  Input field reads:
!                    cv_oo(Tc), np
!
!
!              3) Polynomial fitting (function of the dimensional temperature)
!
!                  Parameters: N (integer), ai (vector, real), ni (vector, real)
!
!                    cv_oo(T) = cv_0 + SUM_i=1,N [ ai(i) * T**ni(i) ] 
!                    ni(i) .NE. 0 for all i = 1,N
!
!                  Input fields read:
!                    N
!                    cv_0
!                    ai(1), ai(2), ..... , ai(N)
!                    ni(1), ni(2), ..... , ni(N)
!
!
!              4) Piecewise polynomial fitting (function of the dimensional temperature)
!
!                  Parameters: M (integer), T_k (vector, integer), 
!                              Nk (vector, integer), cv_0 (vector, real), 
!                              aik (matrix, real), nik (matrix, real)  
!
!                    for T_k(k) < T < T_k(k+1), k = 1,M
!                    cv_oo(T) = cv_0(k) + SUM_i=1,Nk(k) [ aik(i, k) * T**nik(i,k) ]
!                    nik(i,k) .NE. 0 for all i = 1,Nk(k), k = 1,M
!
!                  with T_k(1) = 0 and the understanding that T_k(M+1) = oo.
!                  In the internal representation, T_k has dimension M+1, with
!                  T_k(M+1) = oo added for efficiency. 
!
!                  Input fields read:
!                    M
!                      T_k(1),   T_k(2), ...... ,       T_k(M)     
!                       Nk(1),    Nk(2), ...... ,        Nk(M)
!                     cv_0(1),  cv_0(2), ...... ,      cv_0(M)    
!                    aik(1,1), aik(2,1), ...... , aik(Nk(1),1)
!                    nik(1,1), nik(2,1), ...... , nik(Nk(1),1)
!                    aik(1,2), aik(2,2), ...... , aik(Nk(2),2)
!                    nik(1,2), nik(2,2), ...... , nik(Nk(2),2)
! 
!                     ...... ,  ...... , ...... , ...... 
! 
!                    aik(1,1), aik(2,1), ...... , aik(Nk(M),M)
!                    nik(1,1), nik(2,1), ...... , nik(Nk(M),M)
! 
!
!              5) Harmonic oscillator approximation
!
!                  Parameter: cv_const (real), N (integer), Tv(vector,real) 
!
!                    cv_oo(T) = cv_const 
!                             + Rgas * SUM( (Tv(:)/T)^2 * EXP(Tv(:)/T)/(EXP(Tv(:)/T-1)^2 )  
!
!                  Input field reads:
!                    cv_const
!                    N
!                    Tv(1), Tv(2), ...,  Tv(N)
!
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!              Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!   Copyright: 1998-2003 Alberto Guardone
!              See COPYING file for copyright notice
!
!=========================================================================================

   MODULE ideal_specific_heat

   !----------------------------------------------------------------------------
   IMPLICIT NONE
   PRIVATE

   ! Parameters to choose the curve fitting type
   LOGICAL :: coeff_adim
   INTEGER :: function_type

   INTEGER, PARAMETER :: POLYTROPIC             = 1, &
                         NP_cv_oo_Tr_np         = 2, &
                         NP_polynom_Tdim        = 3, &
                         NP_piecepolynom_Tdim   = 4, &
                         NP_harmonic_oscillator = 5

   ! Constant cv for polytropic (1) and polynomial fitting (3)                       
   REAL(KIND=8) :: cv_0

   ! Cv at the critical point and exponent for power law (2)
   REAL(KIND=8) :: cv_oo_Tc, np

   ! Coefficients a and exponents n for polynomial fitting (3) 
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: a_i, n_i
   ! Piecewise polynomial fitting (4):
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: T_k, cv_0_k
   INTEGER,      DIMENSION(:),   ALLOCATABLE :: N_k
   ! Coefficients a and exponents n  
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a_ik, n_ik

   ! Values of phi0 and psi0 at the temperature T_k  
   ! phi0_k = - INT_0^T_k cv_k(t) dt 
   !   + SUM i=1,k-1 [ int_0^Ti cv_i(t) dt - int_0^Ti+1 cv_i(t) dt ]
   ! where cv_k(T) == cv(T)   T_k < T < T_k+1 
   ! psi0_k = - INT_0^T_k cv_k(t)/t dt
   !   + SUM i=1,k-1 [ int_0^Ti cv_i(t)/t dt - int_0^Ti+1 cv_i(t)/t dt ]
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: phi0_k, psi0_k

   ! Vibrational temperature for the harmonic oscillator model
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Tv


   ! DIMENSIONAL Gas properties (stored locally)
   ! Reference temperature [K], pressure [Pa], 
   ! specific volume [m^3/Kg] and gas constant Rgas
   REAL(KIND=8) :: Tc, T_ref, P_ref, v_ref, Rgas,  &
                   cv_ref, cp_ref, Rgas_ref


   PUBLIC :: phi__T,                          &
             phip__T,                         & 
             phipp__T,                        & 
             psi__T,                          &
             cv_oo__Tdim,                     &
             init_ideal_specific_heat,        &
             read_param_ideal_specific_heat,  &
             write_param_ideal_specific_heat
   !----------------------------------------------------------------------------

   CONTAINS
   

   FUNCTION  cv_oo__Tdim(T)   RESULT(cv_oo)

   ! Computes the DIMENSIONAL ideal gas specific heat cv_oo(T)
   ! from the dimensional temperature
   !------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T     ! DIMENSIONAL temperature
   REAL(KIND=8)             :: cv_oo ! DIMENSIONAL cv_oo(T)  
   !------------------------------------------------------------
        
   !cv_oo = phip__T(T/T_ref) / (T_ref/(P_ref*v_ref))

   cv_oo = phip__T(T/T_ref) * cv_ref

   END FUNCTION  cv_oo__Tdim

   
   
   

   FUNCTION  phi__T(T)   RESULT(phi)

   ! Computes the ideal gas contribution to the energy, phi, 
   !
   !            _T
   !    phi = _/  cv_oo(t) dt
   !           T0
   !
   ! In this function, phi is computed from the dimensional 
   ! temperature  (NOTE  that  the  input  temperature  is 
   ! NONDIMENSIONAL it is made dimensional at the beginning
   ! of the function) and constants. Eventually the dimensionless
   ! value of phi is obtained by P_ref v_ref
   !
   ! The constant contribution due to T0 is not computed and
   ! can be retrieved as -phi__T(T0)
   !------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T    ! dimensionless temperature
   REAL(KIND=8)             :: phi  ! dimensionless phi(T)  

   REAL(KIND=8) :: Tr, Td
   INTEGER :: k
   !------------------------------------------------------------
        
   ! DIMENSIONAL temperature
   Td = T * T_ref

   SELECT CASE (function_type)
   
      CASE (POLYTROPIC)

        phi = cv_0*Td

      
      CASE (NP_cv_oo_Tr_np)

        ! Reduced temperature
        Tr = Td/Tc 
        phi = Tc * (cv_oo_Tc/(np+1)) * Tr**(np+1)

         
      CASE (NP_polynom_Tdim)

        phi = cv_0*Td + SUM( (a_i/(n_i + 1)) * Td**(n_i + 1) )

                                
      CASE (NP_piecepolynom_Tdim)

        DO k = 1, SIZE(T_k) - 1
           IF ( (Td > T_k(k)) .AND. (Td < T_k(k+1))) THEN  
              phi = phi0_k(k) + cv_0_k(k)*Td &
                  + SUM( (a_ik(:,k)/(n_ik(:,k) + 1))  &
                         * Td**(n_ik(:,k) + 1))
           ENDIF 
        ENDDO

      
      CASE (NP_harmonic_oscillator)

        phi = cv_0*Td + Rgas*SUM(Tv / (EXP(Tv/Td)-1))

   END SELECT 

   ! Dimensionless value
   phi = phi / (P_ref*v_ref)
   
   END FUNCTION  phi__T





   FUNCTION  phip__T(T)   RESULT(phip)

   !  Computes dphi/dT, namely, the ideal gas contribution to  
   !  the specific heat at constant volume, i.e.,
   !
   !    phi'(T) = cv_oo(T) 
   !
   ! In this function, phi' is computed from the dimensionless
   ! temperature using dimensional temperature and constants
   ! and its made dimensionless by (P_ref v_ref)/T_ref, 
   ! eventually
   !------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T     ! dimensionless temperature 
   REAL(KIND=8)             :: phip  ! dimensionless phi'(T)      

   REAL(KIND=8) :: Tr, Td
   INTEGER :: k
   !------------------------------------------------------------
        
   ! Dimensional temperature
   Td = T * T_ref

   SELECT CASE (function_type)
   
      CASE (POLYTROPIC)
      phip = cv_0
      
      CASE (NP_cv_oo_Tr_np)
      ! Reduced temperature
      Tr = Td/Tc 
      phip = cv_oo_Tc * Tr**np
     
      CASE (NP_polynom_Tdim)
      phip = cv_0 + SUM(a_i * Td**n_i)
                     
      CASE (NP_piecepolynom_Tdim)
      DO k = 1, SIZE(T_k) - 1
        IF ((Td > T_k(k))  .AND.  (Td < T_k(k+1))) THEN
          phip = cv_0_k(k) + SUM(a_ik(:,k) * Td**n_ik(:,k))
        ENDIF 
      ENDDO
      
      CASE (NP_harmonic_oscillator)
      phip = cv_0 &
           + Rgas*SUM( (Tv/Td)**2 * EXP(Tv/Td) / (EXP(Tv/Td)-1)**2 )
        
   END SELECT 

   ! Dimensionless value
   !phip = phip * T_ref/(P_ref*v_ref)
   phip = phip / cv_ref
   
   END FUNCTION  phip__T


   


   FUNCTION  phipp__T(T)   RESULT(phipp)

   ! Computes d2phi/dT2, namely, the derivative with respect to
   ! the temperature of the ideal gas contribution to  
   ! the specific heat at constant volume, i.e.,
   !
   !    phi''(T) = cv_oo'(T)
   !
   ! In this function, phi'' is computed from the dimensionless
   ! temperature using dimensional temperature and constants
   ! and its made dimensionless by (P_ref v_ref)/(T_ref^2), 
   ! eventually
   !------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T      ! dimensionless temperature 
   REAL(KIND=8)             :: phipp  ! dimensionless phi''(T)      

   REAL(KIND=8) :: Tr, Td
   INTEGER :: k
   !------------------------------------------------------------
        
   ! Dimensional temperature
   Td = T * T_ref

   SELECT CASE (function_type)
   
      CASE (POLYTROPIC)
      phipp = 0.d0
      
      CASE (NP_cv_oo_Tr_np)
      ! Reduced temperature
      Tr = Td/Tc 
      phipp = ((cv_oo_Tc * np)/Tc) * Tr**(np-1)
      
      CASE (NP_polynom_Tdim)
      phipp = SUM( a_i * n_i * Td**(n_i -1))
         
      CASE (NP_piecepolynom_Tdim)
      DO k = 1, SIZE(T_k)  - 1
         IF ((Td > T_k(k)) .AND. (Td < T_k(k+1))) THEN  
            phipp = SUM( a_ik(:,k) * n_ik(:,k) * Td**(n_ik(:,k) - 1))
         ENDIF 
      ENDDO
      
      CASE (NP_harmonic_oscillator)
      phipp = Rgas * SUM( &
            ((Tv**2)/(Td**3))*(EXP(Tv/Td)/(EXP(Tv/Td)-1)**2) &
            *(2*(Tv/Td)*EXP(Tv/Td)/(EXP(Tv/Td)-1) -2 - Tv/Td))

   END SELECT 

   ! Dimensionless value
   phipp = phipp * (T_ref**2)/(P_ref*v_ref)

   END FUNCTION  phipp__T





   FUNCTION  psi__T(T)   RESULT(psi)

   !  Computes the ideal gas contribution to the entropy, psi, 
   !
   !              _ T
   !             /   cv_oo(t)
   !    psi =   /   --------  dt
   !          _/        t
   !            T0
   !
   ! In this function, psi is computed from the dimensionless
   ! temperature using dimensional temperature and constants
   ! and its made dimensionless by P_ref v_ref, eventually
   !
   ! The constant contribution due to T0 is not computed and
   ! can be retrieved as -psi__T(T0)
   
   !------------------------------------------------------------
   IMPLICIT NONE 

   REAL(KIND=8), INTENT(IN) :: T     ! dimensionless temperature
   REAL(KIND=8)             :: psi   ! dimensionless psi(T)  

   REAL(KIND=8) :: Tr, Td
   INTEGER :: k
   !------------------------------------------------------------

   ! Dimensional temperature
   Td = T * T_ref

   SELECT CASE (function_type)
   
      CASE (POLYTROPIC)
      psi = cv_0*LOG(Td)
      
      CASE (NP_cv_oo_Tr_np)
      ! Reduced temperature
      Tr = Td/Tc 
      psi = (cv_oo_Tc/np) * Tr**np
      
      CASE (NP_polynom_Tdim)
      psi = cv_0*LOG(Td) + SUM((a_i/n_i) * Td**n_i)
         
      CASE (NP_piecepolynom_Tdim)
      DO k = 1, SIZE(T_k) - 1
         IF ( (Td > T_k(k)) .AND. (Td < T_k(k+1))) THEN  
            psi = psi0_k(k) + cv_0_k(k)*LOG(Td) &
                + SUM( (a_ik(:,k)/n_ik(:,k)) &
                       * Td**n_ik(:,k)       )
         ENDIF 
      ENDDO
      
      CASE (NP_harmonic_oscillator)
      psi = cv_0*LOG(Td) &
          + Rgas*SUM( (Tv/Td) / (EXP(Tv/Td)-1)  - LOG(1-EXP(-Tv/Td)) )
      
   END SELECT 

   ! Dimensionless value
   psi = psi / Rgas_ref

   END FUNCTION  psi__T

 
 


   SUBROUTINE  read_param_ideal_specific_heat(idf)
   !----------------------------------------------------------------------------
   IMPLICIT NONE 

   INTEGER, INTENT(IN) :: idf
   INTEGER :: M, N, k
   !----------------------------------------------------------------------------       
      
   READ(idf,*) coeff_adim
   READ(idf,*) function_type

   SELECT CASE (function_type)
   
      CASE (POLYTROPIC)

        READ(idf,*) cv_0

      
      CASE (NP_cv_oo_Tr_np)

        READ(idf,*) cv_oo_Tc, np


      CASE (NP_polynom_Tdim)

        READ(idf,*) N
        READ(idf,*) cv_0
	
        ALLOCATE (a_i(N), n_i(N))
        
	READ(idf,*) a_i
        READ(idf,*) n_i


      CASE (NP_piecepolynom_Tdim)

        READ(idf,*) M

        ALLOCATE (T_k(M+1), cv_0_k(M), N_k(M))
	
        READ(idf,*) T_k(1:M)
        T_k(M+1) = HUGE(T_k(1)) 

        READ(idf,*) cv_0_k

        READ(idf,*) N_k
        N = MAXVAL(N_k)

	ALLOCATE (a_ik(N,M), n_ik(N,M))
        a_ik = 0.d0
        n_ik = 0
     
	DO k = 1, M
          READ(idf,*) a_ik(1:N_k(k), k)
          READ(idf,*) n_ik(1:N_k(k), k)
        ENDDO
      
        ALLOCATE (phi0_k(M), psi0_k(M))


      CASE (NP_harmonic_oscillator)

        READ(idf,*) cv_0
        READ(idf,*) N

        ALLOCATE (Tv(N))

        READ(idf,*) Tv

   END SELECT 

   END SUBROUTINE  read_param_ideal_specific_heat




 
   SUBROUTINE  init_ideal_specific_heat
   !----------------------------------------------------------------------------
   USE structures,       ONLY: REF
   USE gas_properties,   ONLY: gas
   
   IMPLICIT NONE
   INTEGER :: k
   !----------------------------------------------------------------------------

   ! Reference values from nondimensionalization.
   ! Variables of the module, common to each procedure of the module
   T_ref = REF % T
   P_ref = REF % P
   v_ref = REF % v
   Rgas_ref = REF % Rgas

   ! Initialize reference (ideal) specific heats 
   ! in the gas properties
   REF % cv = REF % Rgas
   REF % cp = REF % Rgas

   cv_ref = REF % cv
   cp_ref = REF % cp

   ! Critical Temperature [K]
   Tc = gas % Tc
   ! Gas constant [J /(kg K)]
   Rgas = gas % Rgas


   SELECT CASE (function_type)
   
      CASE (POLYTROPIC)

        IF (coeff_adim)  cv_0 = cv_0 * cv_ref


      CASE (NP_cv_oo_Tr_np)

        IF(.NOT. gas % vapor) THEN
          WRITE(*,*) ''
          WRITE(*,*) 'ERROR. INIT_IDEAL_SPECIFIC_HEAT:'
          WRITE(*,*) 'The gas ', gas % name,'  is not a vapor,'
          WRITE(*,*) 'dense gas power law is not applicable.'
          WRITE(*,*) ''
          STOP
        ENDIF

        IF (coeff_adim)  cv_oo_Tc = cv_oo_Tc * cv_ref


      CASE (NP_polynom_Tdim)
      
        IF (coeff_adim) THEN
          cv_0 = cv_0 * cv_ref
          a_i  =  a_i * cv_ref
        ENDIF


      CASE (NP_piecepolynom_Tdim)

        IF (coeff_adim) THEN           
          cv_0_k = cv_0_k * cv_ref
          T_k = T_k * T_ref
          DO k = 1, SIZE(T_k)
            a_ik(:,k)  =  a_ik(:,k) * cv_ref 
          ENDDO
        ENDIF

        phi0_k = 0.d0;  psi0_k = 0.d0
     
     
        k = 2
        phi0_k(k) = phi0_k(k-1)                           &
                  - cv_0_k(k-1)*T_k(k)                    &
                  - SUM( (a_ik(:,k-1)/(n_ik(:,k-1) + 1))  &
                         * T_k(k)**(n_ik(:,k-1) + 1)) 

        psi0_k(k) = psi0_k(k-1)                     &
                  - cv_0_k(k-1)*LOG(T_k(k))         &
                  - SUM( (a_ik(:,k-1)/n_ik(:,k-1))  &
                         * T_k(k)**n_ik(:,k))

        DO k = 3, SIZE(T_k) - 1           

           phi0_k(k) = phi0_k(k-1)                           &
                     + cv_0_k(k-1)*T_k(k-1)                  &
                     + SUM( (a_ik(:,k-1)/(n_ik(:,k-1) + 1))  &
                            * T_k(k-1)**(n_ik(:,k-1) + 1))   &
                     - cv_0_k(k-1)*T_k(k)                    &
                     - SUM( (a_ik(:,k-1)/(n_ik(:,k-1) + 1))  &
                            * T_k(k)**(n_ik(:,k-1) + 1)) 

           psi0_k(k) = psi0_k(k-1)                     &
                     + cv_0_k(k-1)*LOG(T_k(k-1))       &
                     + SUM( (a_ik(:,k-1)/n_ik(:,k-1))  &
                            * T_k(k-1)**n_ik(:,k))     &
                     - cv_0_k(k-1)*LOG(T_k(k))         &
                     - SUM( (a_ik(:,k-1)/n_ik(:,k-1))  &
                            * T_k(k)**n_ik(:,k))
        ENDDO

        DO k = 3, SIZE(T_k) - 1

           phi0_k(k) = phi0_k(k)                         &
                     - cv_0_k(k)*T_k(k)                  &
                     - SUM( (a_ik(:,k)/(n_ik(:,k) + 1))  &
                            * T_k(k)**(n_ik(:,k) + 1))

           psi0_k(k) = psi0_k(k)                  &
                     - cv_0_k(k)*LOG(T_k(k))      &
                     - SUM( (a_ik(:,k)/n_ik(:,k)) &
                            * T_k(k)**n_ik(:,k))
        ENDDO

      
      CASE (NP_harmonic_oscillator)

        IF (coeff_adim) THEN
          cv_0 = cv_0 * cv_ref
          Tv   = Tv * T_ref
        ENDIF
      
   END SELECT 
   
   END SUBROUTINE  init_ideal_specific_heat





   SUBROUTINE  write_param_ideal_specific_heat(idf)
   !----------------------------------------------------------------------------
   IMPLICIT NONE 

   INTEGER, INTENT(IN) :: idf
   INTEGER :: M, k
   !----------------------------------------------------------------------------

   ! This writing operation occurs after the initialization 
   ! of the module, so the values of Cv0 and other parameters
   ! have already been transformed in dimensional values
   ! regardless what stated the input file. For this reason
   ! the parameter corersponding to the dimensional/nondimensional
   ! character of these values (coeff_adim) is forced to be false.
   WRITE(idf,*) '  .FALSE.'
   WRITE(idf,*) '  ', function_type

   SELECT CASE (function_type)
   
     CASE (POLYTROPIC)
       WRITE(idf,*) '', cv_0


     CASE (NP_cv_oo_Tr_np)
       WRITE(idf,*)'   ',  cv_oo_Tc, np


     CASE (NP_polynom_Tdim)
       WRITE(idf,*) '   ', SIZE(a_i)
       WRITE(idf,*) '   ', cv_0
       WRITE(idf,*) '   ', a_i
       WRITE(idf,*) '   ', n_i


     CASE (NP_piecepolynom_Tdim)
       M = SIZE(T_k - 1)
       WRITE(idf,*) '   ', M

       WRITE(idf,*) '   ', T_k(1:M)

       WRITE(idf,*) '   ', cv_0_k

       WRITE(idf,*) N_k
       DO k = 1, SIZE(T_k - 1)
          WRITE(idf,*) '   ', a_ik(1:N_k(k), k)
          WRITE(idf,*) '   ', n_ik(1:N_k(k), k)
       ENDDO


     CASE (NP_harmonic_oscillator)
       WRITE(idf,*) '   ', cv_0
       WRITE(idf,*) '   ', SIZE(Tv)
       WRITE(idf,*) '   ', Tv
        
   END SELECT 

   END SUBROUTINE  write_param_ideal_specific_heat


   END MODULE ideal_specific_heat
