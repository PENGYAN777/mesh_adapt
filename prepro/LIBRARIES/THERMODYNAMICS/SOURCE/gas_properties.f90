!===============================================================================
!
!      Module: gas_properties
!
! Description: Definition of the gas properties data 
!              structure and IO procedures
!
!      Author: Alberto Guardone, Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it,
!                      fossati@aero.polimi.it
!
!   Copyright: 1998-2007 Alberto Guardone, Marco Fossati
!              See COPYING file for copyright notice
!
!=============================================================================== 

   MODULE  gas_properties

   !----------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE gas_properties_type

     CHARACTER(LEN=32) :: name
     CHARACTER(LEN=32) :: chemic_formula
     
     ! Gas constant. Initialized while initializing
     ! gas model (i.e. for ideal gas it is initialized in
     ! init_ideal_gas)
     REAL(KIND=8) :: Rgas
     ! Molecular weight [kg/kmol]
     REAL(KIND=8) :: M_mole

     ! Energy and entropy at reference state. (NONDIMENSIONAL)
     REAL(KIND=8) :: e_0, s_0
     ! Reference (ideal) specific heat (not 
     ! scaling but reference state _0). Used for Newton
     ! initial guess in Thermodynamics. (NONDIMENSIONAL)
     REAL(KIND=8) :: cv_0
     ! Specific heats ratio
     REAL(KIND=8) :: gamma


     LOGICAL :: vapor
     ! Pressure [Pa], temperature [K], specific volume [m3/kg]
     ! and compressibility factor Zc = (Pc vc)/(Rgas Tc) at the critical point
     REAL(KIND=8) :: Pc, Tc, vc, Zc
     ! Boiling temperature at standard pressure, i.e., Pb = 1 atm.
     REAL(KIND=8) :: Tb
     ! Acentric factor, can be taken as -log (Ps/Pc) - 1, where Ps 
     ! is the vapor pressure evaluated at T/Tc = 0.7
     REAL(KIND=8) :: omega


     ! STRUCTURE OF THE MOLECULA:
     ! Number of thermal degrees of freedom 
     INTEGER :: thd
     ! Dimensionless dipole moment for polar moleculas
     REAL(KIND=8) :: mur
     ! Association factor for hydrogen bond moleculas
     REAL(KIND=8) :: kH

   END TYPE gas_properties_type


   TYPE (gas_properties_type) :: gas

   REAL(KIND=8), PARAMETER, PRIVATE :: Rgas_univ = 8314.47   ! [J/(kmole K)] 
   !----------------------------------------------------------------------------

   CONTAINS

 
   SUBROUTINE  read_param_gas_properties(idf)
   !----------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN)  ::  idf
   !----------------------------------------------------------------------------

   READ(idf,*) gas % name
   READ(idf,*) gas % chemic_formula
   READ(idf,*) gas % M_mole
   READ(idf,*) gas % vapor

   IF (gas % vapor) THEN
   
     READ(idf,*) gas % Tc
     READ(idf,*) gas % Pc

     ! [Pa]   =   [atm]  * 101325.d0
     gas % Pc = gas % Pc * 101325.d0

     ! gas % Zc can be changed depending on the EOS
     READ(idf,*) gas % Zc
     READ(idf,*) gas % Tb

     gas%omega = 3.d0*(LOG10(gas%Pc/101325)*(gas%Tb/gas%Tc) &
               / (1.d0 - gas%Tb/gas%Tc) /7.d0) - 1.d0

   ENDIF


   READ(idf,*) gas % thd
   READ(idf,*) gas % mur
   READ(idf,*) gas % kH

   gas % Rgas = Rgas_univ / gas % M_mole
   gas % gamma = ((gas % thd + 3.d0) + 2.d0)/(gas % thd + 3.d0)

   END SUBROUTINE  read_param_gas_properties 





   SUBROUTINE  write_param_gas_properties(idf)
   !----------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: idf            
   !----------------------------------------------------------------------------

   WRITE(idf,*) '   ', gas % name
   WRITE(idf,*) '   ', gas % chemic_formula
   WRITE(idf,*) ' ',   gas % M_mole
   WRITE(idf,*) '  ',  gas % vapor
   
   IF (gas % vapor) THEN
   
     WRITE(idf,*) gas % Tc
     WRITE(idf,*) gas % Pc / 101325.d0
     WRITE(idf,*) gas % Zc
     WRITE(idf,*) gas % Tb
     
   ENDIF

   WRITE(idf,*) '  ',  gas % thd
   WRITE(idf,*) ' ',   gas % mur
   WRITE(idf,*) ' ',   gas % kH 

   END SUBROUTINE  write_param_gas_properties 


END MODULE  gas_properties
