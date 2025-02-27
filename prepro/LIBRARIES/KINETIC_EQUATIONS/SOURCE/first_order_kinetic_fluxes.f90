!=============================================================================== 
!
!      Module: first_order_kinetic_fluxes
!
! Description: First order kinetic flux
!
!      Author: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!
!   Copyright: 2003-2006 Marco Fossati
!              See COPYING file for copyright notice
!
!===============================================================================

  MODULE  first_order_kinetic_fluxes

  CONTAINS


    FUNCTION  KFVS_FO_flux(ww, j_c_d, eta, dt)   RESULT (phi)

    ! In a local reference frame aligned with the node-pair*, the 
    ! flux results from the following Riemann problem for the BGK
    ! equation for the particle distribution function f:
    !
    !  ( f_t + f_x = 0
    !  ( 
    !  ( f(x,0) = f_0 = H(x)*g_l + (1-H(x))*g_r
    !
    ! with  H(x)  being the Heaviside function. g is a Maxwellian. 
    ! Subscripts  _r  and  _l indicate the states at the left and 
    ! the  right  side of the interface.  Subscript  _0 indicates 
    ! the intermediate (NO CONNECTION with Roe!!!!!) state at the 
    ! interface. g_l, g_r are assumed to be constant in space and
    ! time. The function g is assumed to have the form:
    !
    !  g = g_0
    !
    ! constant in space and time. g_0 is  Maxwellian distribution 
    ! function at the interface.  As a consequence of the choices
    ! ono f_0 and g the solution of Riemann problem at the inter-
    ! face is:
    !
    !  f(x_0,t,u,v,w,xi) = f_0(x - ut)
    !
    ! The macroscopic state at the interface is computed by means
    ! of the compatibility constraint for t = 0 and x = x_0. [See
    ! VKI_LS 1998-03 K.Xu]
    ! 
    ! Finally the flux at the interface is the result of the fol-
    ! lowing integral:
    !
    !           _dt  _
    !  1/dt * _/   _/  u*PSI*f(x_0,t,u,v,w,xi) dudvdwdxi dt
    !          0
    !
    ! where PSI is the vector of collisional invariants.
    !---------------------------------------------------------------------------
    USE kinetic_equations
    USE messages,              ONLY: terminate
    USE metric_coefficients,   ONLY: cosines_ETA_FV

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
    INTEGER,      DIMENSION(:,:), INTENT(IN) :: j_c_d
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: eta
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: dt

    REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,2)) :: phi
    REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,1)) :: Rot
    !REAL(KIND=8), DIMENSION(SIZE(eta,1),SIZE(eta,1)) :: versors

    REAL(KIND=8), DIMENSION(SIZE(ww,1))  :: ww_i, ww_j, phi_ij
    REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ij

    REAL(KIND=8) :: mod_eta_ij

    INTEGER  ::  c, i, j, p, space_dim
    !---------------------------------------------------------------------------

    space_dim = SIZE(eta,1)

    phi = 0.d0
    phi_ij = 0.d0

    ! Loop on the node-pairs
    DO c = 1, SIZE(j_c_d, 2)

       ! i is considered as the left state of Riemann problem
       ! while j is considered as the right state
       i = j_c_d(1,c)
       j = j_c_d(2,c)
       
       ww_i  = ww(:,i)
       ww_j  = ww(:,j)

       ! Integrated normal vector
       eta_ij = eta(:,c)
       mod_eta_ij = SQRT(SUM(eta_ij*eta_ij))

       !CALL cosine_matrix(eta_ij, versors)


       ! Rotation matrix to rotate from original frame
       ! into one aligned with normal vector eta_ij
       Rot = 0.d0
       Rot(1,1) = 1.d0
       Rot(SIZE(ww,1),SIZE(ww,1)) = 1.d0

       DO p = 1, space_dim
         Rot(p+1, 2:1+space_dim) = cosines_ETA_FV(p,:,c)
       ENDDO

       ww_i = MATMUL(Rot, ww_i)
       ww_j = MATMUL(Rot, ww_j)

       IF (L__ww(ww_i) <= 0.d0)  CALL terminate(6,'KFVS_FO','Negative temperature')
       IF (L__ww(ww_j) <= 0.d0)  CALL terminate(6,'KFVS_FO','Negative temperature')

       CALL compute_scalar_moments(ww_i,'L')
       CALL compute_scalar_moments(ww_j,'R')


       ! ij Node-pair numerical flux
       phi_ij = dt(1)*(ww_i(1)*IP__c_PSI_g((/1,0,0/), 'L')  +  ww_j(1)*IN__c_PSI_g((/1,0,0/), 'R'))

       ! Inverse of the rotation matrix
       !DO p = 1, space_dim
       !  Rot(p+1, 2:1+space_dim) = versors(:,p)
       !ENDDO

       phi_ij = MATMUL(Rot, phi_ij)

       ! Projection over edge ij
       phi_ij = mod_eta_ij * phi_ij/dt(1)


       ! Accumulate fluxes contribute in the two nodes.
       phi(:,i) = phi(:,i) + phi_ij
       phi(:,j) = phi(:,j) - phi_ij

    ENDDO

    END FUNCTION  KFVS_FO_flux





    FUNCTION  BGK_FO_flux(ww, j_c_d, eta, dt)   RESULT (phi)

    ! In a local reference frame aligned with the node-pair*, the 
    ! flux results from the following Riemann problem for the BGK
    ! equation for the particle distribution function f:
    !
    !  ( f_t + f_x = (g - f) / tau
    !  ( 
    !  ( f(x,0) = f_0 = H(x)*g_l + (1-H(x))*g_r
    !
    ! with  H(x)  being the Heaviside function. g is a Maxwellian. 
    ! Subscripts  _r  and  _l indicate the states at the left and 
    ! the  right  side of the interface.  Subscript  _0 indicates 
    ! the intermediate (NO CONNECTION with Roe!!!!!) state at the 
    ! interface. g_l, g_r are assumed to be constant in space and
    ! time. The function g is assumed to have the form:
    !
    !  g = g_0
    !
    ! constant in space and time. g_0 is  Maxwellian distribution 
    ! function at the interface.  As a consequence of the choices
    ! ono f_0 and g the solution of Riemann problem at the inter-
    ! face is:
    !
    !  f(x_0,t,u,v,w,xi) = (1-EXP(-t/tau))*g_0 + EXP(-t/tau)*f_0
    !
    ! The macroscopic state at the interface is computed by means
    ! of the compatibility constraint for t = 0 and x = x_0. [See
    ! VKI_LS 1998-03 K.Xu]
    ! 
    ! Finally the flux at the interface is the result of the fol-
    ! lowing integral:
    !
    !           _dt  _
    !  1/dt * _/   _/  u*PSI*f(x_0,t,u,v,w,xi) dudvdwdxi dt
    !          0
    !
    ! where PSI is the vector of collisional invariants. For the 
    ! case considered we have:
    !
    !                  _+oo                     _+oo             _0
    !   F_0 = gamma_1* _/  u PSI g_0  +  gamma_2*(_/  u PSI g_l  + _/  u PSI g_r  )
    !                 -oo                      0                -oo
    !
    ! where:
    !
    !   gamma_2 -> is the first term deriving from the initial 
    !            condition f_0.
    !   
    !   gamma_1 -> is the first term deriving from the integral 
    !            term in the form of f. It is associated to
    !            the form of function g.
    !
    !   These are the result of the time integration and averaging
    ! 
    !   *Since the flux is the result of integration over microsco-
    !    pic velocities the choice of the reference frame in which
    !    evaluate the Riemann problem is arbitrary. For convenience
    !    such a frame has been chosen with x-axis aligned with the
    !    node-pair
    !---------------------------------------------------------------------------
    USE kinetic_equations
    USE messages,              ONLY: terminate
    USE metric_coefficients,   ONLY: cosines_ETA_FV

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
    INTEGER,      DIMENSION(:,:), INTENT(IN) :: j_c_d
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: eta
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: dt

    REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,2)) :: phi
    REAL(KIND=8), DIMENSION(SIZE(ww,1), SIZE(ww,1)) :: Rot

    !REAL(KIND=8), DIMENSION(SIZE(eta,1),SIZE(eta,1)) :: versors

    REAL(KIND=8), DIMENSION(SIZE(ww,1))  :: ww_i, ww_j, phi_ij, ww_0
    REAL(KIND=8), DIMENSION(SIZE(eta,1)) :: eta_ij

    REAL(KIND=8) :: gamma_2, gamma_1, dt_

    REAL(KIND=8) :: mod_eta_ij, tau

    INTEGER, DIMENSION(3) :: order_1

    INTEGER  ::  c, i, j, p, space_dim

    LOGICAL :: local_time_step
    !---------------------------------------------------------------------------

    space_dim = SIZE(eta,1)

    phi = 0.d0
    phi_ij = 0.d0

    ! Loop on the node-pairs
    DO c = 1, SIZE(j_c_d, 2)

       ! i is considered as the left state of Riemann problem
       ! while j is considered as the right state
       i = j_c_d(1,c)
       j = j_c_d(2,c)
       
       ww_i  = ww(:,i)
       ww_j  = ww(:,j)

       ! Integrated normal vector
       eta_ij = eta(:,c)
       mod_eta_ij = SQRT(SUM(eta_ij*eta_ij))

       !CALL cosine_matrix(eta_ij, versors)

       ! Rotation matrix to rotate from original frame
       ! into one aligned with normal vector eta_ij
       Rot = 0.d0
       Rot(1,1) = 1.d0
       Rot(SIZE(ww,1),SIZE(ww,1)) = 1.d0

       DO p = 1, space_dim
         Rot(p+1, 2:1+space_dim) = cosines_ETA_FV(p,:,c)
       ENDDO

       ww_i = MATMUL(Rot, ww_i)
       ww_j = MATMUL(Rot, ww_j)

       IF (L__ww(ww_i) <= 0.d0)  CALL terminate(6,'BGK_FO','Negative temperature')
       IF (L__ww(ww_j) <= 0.d0)  CALL terminate(6,'BGK_FO','Negative temperature')
       
       CALL compute_scalar_moments(ww_i, 'L')
       CALL compute_scalar_moments(ww_j, 'R')

       ww_0 = ww_i(1)*IP__c_PSI_g((/0,0,0/), 'L')  +  ww_j(1)*IN__c_PSI_g((/0,0,0/), 'R')

       IF (L__ww(ww_0) <= 0.d0)  CALL terminate(6,'BGK_FO:','Negative temperature.')	
	
       CALL compute_scalar_moments(ww_0, 'I')
      
	
       ! Integrated-Averaged time dependent terms in the final
       ! form of distribution function.

       ! Recover information about local time stepping
       IF (ALL(dt == dt(1)))  local_time_step = .FALSE.

       IF (.NOT. local_time_step) THEN

         dt_ = dt(1)

         tau = mean_collision_time(dt_, ww_i, ww_j, ww_0)

         gamma_1 = (dt_ + tau*EXP(-dt_/tau) - tau)
         gamma_2 = tau*(1 - EXP(-dt_/tau))


         ! ij Node-pair numerical flux
	 order_1 = (/1, 0, 0/)
	 
         phi_ij = gamma_1 * ww_0(1) * I__c_PSI_g(order_1,'I')    +  &
                  gamma_2 * (ww_i(1) * IP__c_PSI_g(order_1,'L')  +  &
                             ww_j(1) * IN__c_PSI_g(order_1,'R'))

         ! Inverse of the rotation matrix
         !DO p = 1, space_dim
         !  Rot(p+1, 2:1+space_dim) = versors(:,p)
         !ENDDO

         phi_ij = MATMUL(Rot, phi_ij)

         ! Projection over edge ij
         phi_ij = mod_eta_ij * phi_ij/dt(1)


         ! Accumulate fluxes contribute in the two nodes
         phi(:,i) = phi(:,i) + phi_ij
         phi(:,j) = phi(:,j) - phi_ij

       ELSE

         PRINT*, ''
         PRINT*, 'ERROR. BGK_FO_FLUX:'
         PRINT*, 'Local time stepping not impemented yet'
         PRINT*, ''
         STOP

       ENDIF

    ENDDO

    END FUNCTION  BGK_FO_flux


  END MODULE  first_order_kinetic_fluxes
