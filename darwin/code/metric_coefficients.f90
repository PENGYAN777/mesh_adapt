!============================================================ 
!
!      Module: metric_coefficients
!
! Description: Definition of the metrics quantities 
!              pertaining to the domain and its boundary and 
!              of IO subroutines
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!   Copyright: 1998-2003 Alberto Guardone, Marco Fossati
!              See COPYING file for copyright notice
!
!============================================================ 

!=====================================================================
!
!   ---------------------------
!   Domain node-pair quantities
!   ---------------------------
!
!            mass(c)   Mass associated with the c-th node-pair
!
!           eta(k,c)   Metric vector (of the physical space)
!                      associated with the c-th node-pair
!
!           stiff(c)   Stiffness associated with the c-th node-pair
!
!       rot_stiff(c)   Rotated stiffness (antisymmetric) associated 
!                      with the c-th node-pair
!
!     stiff_T(x,y,c)   Stiffness tensor associate to the c-th node-pair
!
!           Drr(1,c)   Length of the node-pair  i--j
!           Drr(2,c)   Length of the extension i*--j
!           Drr(3,c)   Length of the extension  i--j*
!
!   ----------------------
!   Domain node quantities
!   ----------------------
!
!            cell(i)   Cell size == diagonally lumped mass matrix
!
!         mass_ii(i)   Diagonal element of the mass matrix
!
!        stiff_ii(i)   Diagonal element of the stiffness matrix
!
!
!   -----------------------------
!   Boundary node-pair quantities
!   -----------------------------
!
!          mass_b(c)   Mass (surface integral) associated with
!                      the c-th surface node-pair
!
!         chi_b(k,c)   Boundary metric vector (of the physical space)
!                      associated with the c-th surface node-pair
!
!         kappa_b(:,k,c)   Boundary metric vector (of the physical space)
!                          associated with the c-th surface node-pair
!                          first index for node i, second for node j
!        
!         Kb_ijs(:,:,c)  Symmetric term of the boundary metric tensor 
!                        associated to the c-th boundary node-pair
!
!         Kb_ija(:,:,c)  Antisymmetric term of the boundary metric tensor 
!                        associated to the c-th boundary node-pair

!   ------------------------
!   Boundary node quantities
!   ------------------------
!
!         xi_bp(k,i)   Boundary metric vector (of the physical space)
!                      associated with the i-th surface point
!
!
!=====================================================================

MODULE  metric_coefficients


   !=====================================================================
   USE nodes,               ONLY : Np_d, Np_b, k_d
   USE node_pair_structure, ONLY : Nc_d, Nc_b, Nc_fv
   !=====================================================================

   LOGICAL, PARAMETER   ::  fv_metric = .TRUE.

   !=====================================================================
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)     :: mass, stiff!, rot_stiff
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: eta
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: stiff_T
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: Drr

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)     :: cell   
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)     :: mass_ii, stiff_ii 

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)     :: mass_b
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: chi_b
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: kappa_b

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: xi_bp

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: eta_fv
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: Drr_fv
   
   TARGET  ::  eta, Drr, eta_fv, Drr_fv
   !=====================================================================


!=====================================================================
CONTAINS
!=====================================================================



   !=====================================================================
   SUBROUTINE save_metric(idf, name, name_length)
   !=====================================================================


      !---------------------------------------------------------------------
      IMPLICIT NONE
  
      INTEGER,           INTENT(IN)  ::  idf
      CHARACTER(LEN=64), INTENT(IN)  ::  name
      INTEGER,           INTENT(IN)  ::  name_length
      !---------------------------------------------------------------------
      INTEGER  ::  c, i, k
      !---------------------------------------------------------------------


      ! Writes header
      ! -------------
      WRITE(idf,1000)
      WRITE(idf,1010) name(1:name_length)

      WRITE(idf,1000)    
      WRITE(idf,1020)
      IF (fv_metric) THEN
         WRITE(idf,1021) k_d, Nc_d, Np_d, Nc_b, Np_b, Nc_fv 
      ELSE
         WRITE(idf,1021) k_d, Nc_d, Np_d, Nc_b, Np_b
      ENDIF

      WRITE(idf,1000)
      WRITE(idf,1025)
    
    
      ! ------------------------------------------------------------    
      ! Domain quantities
      ! -----------------

      ! Writes header
      ! -------------
      WRITE(idf,1000)
      WRITE(idf,1040)
      WRITE(idf,1044)
      WRITE(idf,1045)
      WRITE(idf,1036)
      WRITE(idf,1049)

      ! Volume node-pairs associated quantities
      ! ---------------------------------------
      DO c = 1, Nc_d
         WRITE(idf,1046) c, mass(c), stiff(c) 
         WRITE(idf,1048) eta(:,c)
         WRITE(idf,1048) Drr(:,c)
         DO k = 1, k_d
            WRITE(idf,1048) stiff_T(k,:,c)
         ENDDO      
      ENDDO   

      ! Writes header
      ! -------------
      WRITE(idf,1000)
      WRITE(idf,1080)
      WRITE(idf,1050)
      WRITE(idf,1051)

      ! Volume nodes associated quantities
      ! ----------------------------------
      DO i = 1, Np_d
         WRITE(idf,1046) i, mass_ii(i), cell(i)
         WRITE(idf,1048) stiff_ii(i)
      ENDDO   
      !------------------------------------------------------------    
    

      !------------------------------------------------------------    
      ! Boundary quantities
      ! -------------------

      ! Writes header
      WRITE(idf,1000)
      WRITE(idf,1055)
      
      WRITE(idf,1000)
      WRITE(idf,1040)
      WRITE(idf,1065)
      WRITE(idf,1066)
      WRITE(idf,1064)

      ! Boundary node-pairs associated quantities
      ! -----------------------------------------
      DO c = 1, Nc_b
         WRITE(idf,1067) c, mass_b(c) 
         WRITE(idf,1068) chi_b(:,c)
         WRITE(idf,1068) kappa_b(:,:,c)
      ENDDO   
    
      ! Writes header
      WRITE(idf,1000)
      WRITE(idf,1080)
      WRITE(idf,1090)

      ! Boundary nodes associated quantities
      ! ------------------------------------
      DO i = 1, Np_b
         WRITE(idf,1091) i, xi_bp(:,i) 
      ENDDO   
      !------------------------------------------------------------    

      IF (fv_metric) THEN

         ! FINITE VOLUME
         ! Writes header
         ! -------------
         WRITE(idf,1000)
         WRITE(idf,1043)
         WRITE(idf,1036)

         ! Volume node-pairs associated quantities
         ! ---------------------------------------
         DO c = 1, Nc_fv
            WRITE(idf,1047) c, eta_fv(:,c)
            WRITE(idf,1048) Drr_fv(:,c)
         ENDDO   

      ENDIF
      
      !DO c = 1, Nc_d
      !   WRITE(idf,*) c,  rot_stiff(c) 
      !ENDDO   

  
1000  FORMAT('######################################################################################')
1010  FORMAT('#    NAME:      ',a15,'                                                      #')
1020  FORMAT('#        K_D        NC_D        NP_D        NC_B        NP_B       NC_FV             #')
1021  FORMAT(6i12)
1025  FORMAT('#  **********  DOMAIN  **********                                                    #')
1040  FORMAT('#  ++++  NODE-PAIR ASSOCIATED METRICS QUANTITIES  ++++                               #')
1044  FORMAT('#        IDX                    MASS                   STIFF                         #')
1045  FORMAT('#                                ETA                                                 #')
1043  FORMAT('#        IDX                  ETA_FV                                                 #')
1036  FORMAT('#                           DRR i--j               DRR is--i               DRR j--js #')
1049  FORMAT('#                       STIFF TENSOR                                                 #')
1050  FORMAT('#        IDX                 MASS_II                    CELL                         #')
1051  FORMAT('#                           STIFF_II                                                 #')
1046  FORMAT(i12,2e24.16)
1048  FORMAT(12x,3e24.16)
1047  FORMAT(i12,3e24.16)
1055  FORMAT('#  **********  BOUNDARY  **********                                                  #')
1065  FORMAT('#        IDX                  MASS_B                                                 #')
1066  FORMAT('#                              CHI_B                                                 #')
1064  FORMAT('#                            KAPPA_B                                                 #')
1067  FORMAT(i12,1e24.16)
1068  FORMAT(12x,3e24.16)
1080  FORMAT('#  ++++  NODE ASSOCIATED METRICS QUANTITIES  ++++                                    #')
1090  FORMAT('#        IDX                   XI_BP                                                 #')
1091  FORMAT(i12,3e24.16)

  
   END SUBROUTINE save_metric
   !=====================================================================



   !=====================================================================
   SUBROUTINE read_metric(idf)
   !=====================================================================
  

      !---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER  ::  idf
      !---------------------------------------------------------------------
      INTEGER  ::  c, i, k, idx 
      CHARACTER(LEN=16)  ::  name
      !---------------------------------------------------------------------


      READ(idf,*) 
      READ(idf,'(16x,a15)') name
      !WRITE(*,*) '   Reading metric quantities: ',name

      READ(idf,*); READ(idf,*); 
      IF (fv_metric) THEN
         READ(idf,*) k_d, Nc_d, Np_d, Nc_b, Np_b, Nc_fv 
      ELSE
         READ(idf,*) k_d, Nc_d, Np_d, Nc_b, Np_b 
      ENDIF

      IF ( ALLOCATED(mass) )      DEALLOCATE( mass ) 
      IF ( ALLOCATED(stiff) )     DEALLOCATE( stiff ) 
      !IF ( ALLOCATED(rot_stiff) ) DEALLOCATE( rot_stiff ) 
      IF ( ALLOCATED(eta) )       DEALLOCATE( eta ) 
      IF ( ALLOCATED(Drr) )       DEALLOCATE( Drr ) 
      IF ( ALLOCATED(stiff_T) )   DEALLOCATE( stiff_T ) 
      IF ( ALLOCATED(mass_ii) )   DEALLOCATE( mass_ii ) 
      IF ( ALLOCATED(stiff_ii) )  DEALLOCATE( stiff_ii ) 
      IF ( ALLOCATED(cell) )      DEALLOCATE( cell ) 
      IF ( ALLOCATED(mass_b) )    DEALLOCATE( mass_b ) 
      IF ( ALLOCATED(chi_b) )     DEALLOCATE( chi_b ) 
      IF ( ALLOCATED(kappa_b) )   DEALLOCATE( kappa_b ) 
      IF ( ALLOCATED(xi_bp) )     DEALLOCATE( xi_bp )
      IF (fv_metric) THEN
         IF ( ALLOCATED(eta_fv) )   DEALLOCATE( eta_fv )
         IF ( ALLOCATED(Drr_fv) )   DEALLOCATE( Drr_fv )
      ENDIF

      ALLOCATE(mass(Nc_d),  stiff(Nc_d), & ! rot_stiff(Nc_d), 
               eta(k_d, Nc_d), stiff_T(k_d, k_d, Nc_d)) 
      ALLOCATE( Drr(3,Nc_d) )
      ALLOCATE(stiff_ii(Np_d),mass_ii(Np_d), cell(Np_d))

      ALLOCATE(mass_b(Nc_b), chi_b(k_d, Nc_b), kappa_b(2,k_d,Nc_b), &
               xi_bp(k_d, Np_b))
      IF (fv_metric) ALLOCATE( eta_fv(k_d, Nc_fv), Drr_fv(3,Nc_fv)  )

      !READ(idf,*); READ(idf,*); 
   
      ! Skips header
      READ(idf,*); READ(idf,*); READ(idf,*); READ(idf,*) 
      READ(idf,*); READ(idf,*); READ(idf,*); READ(idf,*) 

      ! Volume node-pairs associated quantities
      ! ---------------------------------------
      DO c = 1, Nc_d
         READ(idf,*) idx, mass(c), stiff(c) 
         READ(idf,*) eta(:,c)
         READ(idf,*) Drr(1,c),  Drr(2,c),  Drr(3,c)
         DO k = 1, k_d
            READ(idf,*) stiff_T(k,:,c)
         ENDDO      
      ENDDO   

      ! Skips header
      READ(idf,*); READ(idf,*); READ(idf,*); READ(idf,*)

      ! Volume nodes associated quantities
      ! ----------------------------------
      DO i = 1, Np_d
         READ(idf,*) idx, mass_ii(i), cell(i)
         READ(idf,*) stiff_ii(i)
      ENDDO   
   
      ! Skips header
      READ(idf,*); READ(idf,*); READ(idf,*); READ(idf,*)
      READ(idf,*); READ(idf,*); READ(idf,*)

      ! Boundary node-pairs associated quantities
      ! -----------------------------------------
      DO c = 1, Nc_b
         READ(idf,*) idx, mass_b(c) 
         READ(idf,*) chi_b(:,c)
         READ(idf,*) kappa_b(:,:,c)
      ENDDO   

      ! Skips header
      READ(idf,*); READ(idf,*); READ(idf,*) 

      ! Boundary nodes associated quantities
      ! ------------------------------------
      DO i = 1, Np_b
         READ(idf,*) idx, xi_bp(:,i) 
      ENDDO   


      ! Volume node-pairs associated quantities
      ! ---------------------------------------
      IF (fv_metric) THEN
      READ(idf,*);  READ(idf,*);  READ(idf,*)
      DO c = 1, Nc_fv
         READ(idf,*) idx, eta_fv(:,c)
         READ(idf,*) Drr_fv(:,c)
      ENDDO   
      ENDIF

      !DO c = 1, Nc_d
      !   READ(idf,*) idx,  rot_stiff(c) 
      !ENDDO   

  
   END SUBROUTINE read_metric
   !=====================================================================


END MODULE  metric_coefficients
