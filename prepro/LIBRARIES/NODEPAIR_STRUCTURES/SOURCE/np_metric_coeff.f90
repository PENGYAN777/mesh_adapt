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
!   Copyright: 1998-2003 Alberto Guardone
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
   
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: cosines_ETA_FV
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: cosines_XI_BP

   ! Axisymmetric metric quantities
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)     :: mass_y
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: eta_y
      
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)     :: cell_y   
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:)     :: mass_y_ii
   
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: xi_y_bp
   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: eta_y_fv   

   TARGET :: eta, Drr, eta_fv, eta_y_fv, Drr_fv   
   !REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)   :: chi_y_b
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
      INTEGER  ::  c, i, k, j
      !---------------------------------------------------------------------


      ! Writes header
      ! -------------
      WRITE(idf,1000)
      WRITE(idf,1010) name(1:name_length)

      WRITE(idf,1000)    
      WRITE(idf,1020)
      WRITE(idf,1021) k_d, Nc_d, Np_d, Nc_b, Np_b, Nc_fv 

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

      IF (k_d == 1) WRITE(idf,1036)      
      IF (k_d == 2) WRITE(idf,1037)      
      IF (k_d == 3) WRITE(idf,1038)

      WRITE(idf,1099)

      ! Volume node-pairs associated quantities
      ! ---------------------------------------
      DO c = 1, Nc_d
         WRITE(idf,1046) c, mass(c), stiff(c) 
         WRITE(idf,1048) eta(:,c)

         IF (k_d == 1) WRITE(idf,1049) Drr(:,c)
         IF (k_d == 2) WRITE(idf,1050) Drr(:,c)
         IF (k_d == 3) WRITE(idf,1051) Drr(:,c)

         DO k = 1, k_d
            WRITE(idf,1048) stiff_T(k,:,c)
         ENDDO      
      ENDDO   

      ! Writes header
      ! -------------
      WRITE(idf,1000)
      WRITE(idf,1080)
      WRITE(idf,1097)
      WRITE(idf,1098)

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

      ! FINITE VOLUME
      ! Writes header
      ! -------------
      WRITE(idf,1000)
      WRITE(idf,1043)

      IF (k_d == 1) WRITE(idf,1036)   
      IF (k_d == 2) WRITE(idf,1037)   
      IF (k_d == 3) WRITE(idf,1038)


      ! Volume node-pairs associated quantities
      ! ---------------------------------------
      DO c = 1, Nc_fv
         WRITE(idf,1047) c, eta_fv(:,c)
         IF (k_d == 1) WRITE(idf,1049) Drr_fv(:,c)
         IF (k_d == 2) WRITE(idf,1050) Drr_fv(:,c)
         IF (k_d == 3) WRITE(idf,1051) Drr_fv(:,c)
      ENDDO   

      
      ! ROTATION MATRICES IN FINITE VOLUMES
      ! Domain: consines_ETA_FV
      WRITE(idf,1000)
      WRITE(idf,1143)

      DO c = 1, Nc_fv
        IF (k_d == 1) THEN
          WRITE(idf,1144) c, cosines_ETA_FV(1,:,c)
	ELSEIF (k_d == 2) THEN
	  WRITE(idf,1145) c, cosines_ETA_FV(1,:,c)
	  WRITE(idf,1146)    cosines_ETA_FV(2,:,c)
	ELSE
          WRITE(idf,1147) c, cosines_ETA_FV(1,:,c)
	  WRITE(idf,1148)    cosines_ETA_FV(2,:,c)
          WRITE(idf,1149)    cosines_ETA_FV(3,:,c)
	ENDIF
      ENDDO

      ! Boundary: cosines_XI_BP
      WRITE(idf,1000)
      WRITE(idf,1190)

      DO j = 1, Np_b
        IF (k_d == 1) THEN
          WRITE(idf,1144) j, cosines_XI_BP(1,:,j)
	ELSEIF (k_d == 2) THEN
	  WRITE(idf,1145) j, cosines_XI_BP(1,:,j)
	  WRITE(idf,1146)    cosines_XI_BP(2,:,j)
	ELSE
          WRITE(idf,1147) j, cosines_XI_BP(1,:,j)
	  WRITE(idf,1148)    cosines_XI_BP(2,:,j)
          WRITE(idf,1149)    cosines_XI_BP(3,:,j)
	ENDIF
      ENDDO

  
1000  FORMAT('######################################################################################')
1010  FORMAT('#    NAME:      ',a15,'                                                      #')
1020  FORMAT('#        K_D        NC_D        NP_D        NC_B        NP_B       NC_FV             #')
1021  FORMAT(6i12)
1025  FORMAT('#  **********  DOMAIN  **********                                                    #')
1040  FORMAT('#  ++++  NODE-PAIR ASSOCIATED METRICS QUANTITIES  ++++                               #')
1044  FORMAT('#        IDX                    MASS                   STIFF                         #')
1045  FORMAT('#                                ETA                                                 #')
1043  FORMAT('#        IDX                  ETA_FV                                                 #')
1143  FORMAT('#        IDX                  COSINES_ETA_FV                                         #')
1036  FORMAT('#', 27x,'DRR i--j', 15x, 'DRR is--i', 15x, 'DRR j--js', 15x, '#')
1037  FORMAT('#', 27x,'DRR i--j', 15x, 'DRR is--i', 15x, 'DRR j--js', 15x, 'DRR extended' )
1038  FORMAT('#', 27x,'DRR i--j', 15x, 'DRR is--i', 15x, 'DRR j--js', 15x, 'DRR extended' )
1099  FORMAT('#                       STIFF TENSOR                                                 #')
1097  FORMAT('#        IDX                 MASS_II                    CELL                         #')
1098  FORMAT('#                           STIFF_II                                                 #')
1046  FORMAT(i12,2e24.16)
1048  FORMAT(12x,1e24.16)
1049  FORMAT(12x,3e24.16)
1050  FORMAT(12x,7e24.16)
1051  FORMAT(12x,11e24.16)
1047  FORMAT(i12,3e24.16)
1055  FORMAT('#  **********  BOUNDARY  **********                                                  #')
1065  FORMAT('#        IDX                  MASS_B                                                 #')
1066  FORMAT('#                              CHI_B                                                 #')
1064  FORMAT('#                            KAPPA_B                                                 #')
1067  FORMAT(i12,1e24.16)
1068  FORMAT(12x,3e24.16)
1080  FORMAT('#  ++++  NODE ASSOCIATED METRICS QUANTITIES  ++++                                    #')
1090  FORMAT('#        IDX                   XI_BP                                                 #')
1190  FORMAT('#        IDX                  COSINES_XI_BP                                          #')
1091  FORMAT(i12,3e24.16)
1144  FORMAT(i12,1e24.16)
1145  FORMAT(i12,2e24.16)
1146  FORMAT(12x,2e24.16)
1147  FORMAT(i12,3e24.16)
1148  FORMAT(12x,3e24.16)
1149  FORMAT(12x,3e24.16)

   END SUBROUTINE save_metric
   !=====================================================================



   !=====================================================================
   SUBROUTINE save_metric_axisymmetric(idf, name, name_length)
   !=====================================================================


      !---------------------------------------------------------------------
      IMPLICIT NONE
  
      INTEGER,           INTENT(IN)  ::  idf
      CHARACTER(LEN=64), INTENT(IN)  ::  name
      INTEGER,           INTENT(IN)  ::  name_length
      !---------------------------------------------------------------------
      INTEGER  ::  c, i
      !---------------------------------------------------------------------


      ! Writes header
      ! -------------
      WRITE(idf,2000)
      WRITE(idf,2010) name(1:name_length)

      WRITE(idf,2000)    
      WRITE(idf,2020)
      WRITE(idf,2021) k_d, Nc_d, Np_d, Nc_b, Np_b, Nc_fv 

      WRITE(idf,2000)
      WRITE(idf,2025)
    
    
      ! ------------------------------------------------------------    
      ! Domain quantities
      ! -----------------

      ! Writes header
      ! -------------
      WRITE(idf,2000)
      WRITE(idf,2040)
      WRITE(idf,2044)
      WRITE(idf,2045)

      IF (k_d == 1) WRITE(idf,2036)
      IF (k_d == 2) WRITE(idf,2037)
      IF (k_d == 3) WRITE(idf,2038)
      

      ! Volume node-pairs associated quantities
      ! ---------------------------------------
      DO c = 1, Nc_d
         WRITE(idf,2046) c, mass_y(c)
         WRITE(idf,2048) eta_y(:,c)

         IF (k_d == 1) WRITE(idf,2049) Drr(:,c)
         IF (k_d == 2) WRITE(idf,2050) Drr(:,c)
         IF (k_d == 3) WRITE(idf,2051) Drr(:,c)

      ENDDO   

      ! Writes header
      ! -------------
      WRITE(idf,2000)
      WRITE(idf,2080)
      WRITE(idf,2060)

      ! Volume nodes associated quantities
      ! ----------------------------------
      DO i = 1, Np_d
         WRITE(idf,2046) i, mass_y_ii(i), cell_y(i)
      ENDDO   
      !------------------------------------------------------------    
    

      !------------------------------------------------------------    
      ! Boundary quantities
      ! -------------------

      ! Writes header
      WRITE(idf,2000)
      WRITE(idf,2055)
      
      ! Writes header
      WRITE(idf,2000)
      WRITE(idf,2080)
      WRITE(idf,2090)

      ! Boundary nodes associated quantities
      ! ------------------------------------
      DO i = 1, Np_b
         WRITE(idf,2091) i, xi_y_bp(:,i) 
      ENDDO   
      !------------------------------------------------------------
      

      ! FINITE VOLUME
      ! Writes header
      ! -------------
      WRITE(idf,2000)
      WRITE(idf,2043)
      
      IF (k_d == 1) WRITE(idf,2036)   
      IF (k_d == 2) WRITE(idf,2037)   
      IF (k_d == 3) WRITE(idf,2038)


      ! Volume node-pairs associated quantities
      ! ---------------------------------------
      DO c = 1, Nc_fv
         WRITE(idf,2047) c, eta_y_fv(:,c)
         IF (k_d == 1) WRITE(idf,2049) Drr_fv(:,c)
         IF (k_d == 2) WRITE(idf,2050) Drr_fv(:,c)
         IF (k_d == 3) WRITE(idf,2051) Drr_fv(:,c)
      ENDDO   

  
2000  FORMAT('######################################################################################')
2010  FORMAT('#    NAME:      ',a15,'                                                      #')
2020  FORMAT('#        K_D        NC_D        NP_D        NC_B        NP_B       NC_FV             #')
2021  FORMAT(6i12)
2025  FORMAT('#  **********  DOMAIN  **********                                                    #')
2040  FORMAT('#  ++++  NODE-PAIR ASSOCIATED METRICS QUANTITIES  ++++                               #')
2044  FORMAT('#        IDX                  MASS_y                                          #')
2045  FORMAT('#                              ETA_y                                                 #')
2043  FORMAT('#        IDX                ETA_y_FV                                                 #')
2036  FORMAT('#', 27x,'DRR i--j', 15x, 'DRR is--i', 15x, 'DRR j--js', 15x, '#')
2037  FORMAT('#', 27x,'DRR i--j', 15x, 'DRR is--i', 15x, 'DRR j--js', 15x, 'DRR extended' )
2038  FORMAT('#', 27x,'DRR i--j', 15x, 'DRR is--i', 15x, 'DRR j--js', 15x, 'DRR extended' )
2060  FORMAT('#        IDX               MASS_y_II                  CELL_y                         #')
2046  FORMAT(i12,2e24.16)
2048  FORMAT(12x,3e24.16)
2047  FORMAT(i12,3e24.16)
2049  FORMAT(12x,3e24.16)
2050  FORMAT(12x,7e24.16)
2051  FORMAT(12x,11e24.16)
2055  FORMAT('#  **********  BOUNDARY  **********                                                  #')
2080  FORMAT('#  ++++  NODE ASSOCIATED METRICS QUANTITIES  ++++                                    #')
2090  FORMAT('#        IDX                 XI_y_BP                                                 #')
2091  FORMAT(i12,3e24.16)

  
   END SUBROUTINE save_metric_axisymmetric
   !=====================================================================




   SUBROUTINE read_metric(idf)
   !---------------------------------------------------------------------
   USE mp_interface
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf

   INTEGER :: c, i, k, idx 
   CHARACTER(LEN=16) :: name
   !---------------------------------------------------------------------

   READ(idf,*) 
   READ(idf,'(16x,a15)') name
   READ(idf,*); READ(idf,*);
    
   IF (.NOT. MP_job  .OR.  MP_master) THEN
     READ(idf,*) k_d, Nc_d, Np_d, Nc_b, Np_b, Nc_fv
     NcFE_G = Nc_d
     Nj_G = Np_d
     NcFE_b = Nc_b
     Njb_G = Np_b
     NcFV_G = Nc_fv
   ELSE
     READ(idf,*) k_d
     Nc_d = NcTFEP
     Np_d = NjTP
     Nc_b = NcFEP_b
     Np_b = 2*NjP_b
     Nc_fv = NcTFVP
   ENDIF  

   IF (ALLOCATED(mass))      DEALLOCATE(mass) 
   IF (ALLOCATED(stiff))     DEALLOCATE(stiff) 
   IF (ALLOCATED(eta))       DEALLOCATE(eta) 
   IF (ALLOCATED(Drr))       DEALLOCATE(Drr) 
   IF (ALLOCATED(stiff_T))   DEALLOCATE(stiff_T) 
   IF (ALLOCATED(mass_ii))   DEALLOCATE(mass_ii) 
   IF (ALLOCATED(stiff_ii))  DEALLOCATE(stiff_ii) 
   IF (ALLOCATED(cell))      DEALLOCATE(cell) 
   IF (ALLOCATED(mass_b))    DEALLOCATE(mass_b) 
   IF (ALLOCATED(chi_b))     DEALLOCATE(chi_b) 
   IF (ALLOCATED(kappa_b))   DEALLOCATE(kappa_b) 
   IF (ALLOCATED(xi_bp))     DEALLOCATE(xi_bp)
   IF (ALLOCATED(eta_fv))    DEALLOCATE(eta_fv)
   IF (ALLOCATED(Drr_fv))    DEALLOCATE(Drr_fv)

   IF (ALLOCATED(cosines_ETA_FV))   DEALLOCATE(cosines_ETA_FV)
   IF (ALLOCATED(cosines_XI_BP))    DEALLOCATE(cosines_XI_BP)
      
   !IF (ALLOCATED(rot_stiff)) DEALLOCATE(rot_stiff)
      

   ALLOCATE(mass(Nc_d), stiff(Nc_d), eta(k_d, Nc_d), stiff_T(k_d, k_d, Nc_d)) 
   ALLOCATE(Drr(4*k_d-1,Nc_d))
   ALLOCATE(stiff_ii(Np_d), mass_ii(Np_d), cell(Np_d))
   ALLOCATE(mass_b(Nc_b), chi_b(k_d, Nc_b), kappa_b(2,k_d,Nc_b), xi_bp(k_d, Np_b))
   ALLOCATE(eta_fv(k_d, Nc_fv), Drr_fv(4*k_d-1,Nc_fv))
   
   ALLOCATE(cosines_ETA_FV(k_d, k_d, Nc_d))
   ALLOCATE(cosines_XI_BP(k_d, k_d, Np_b))
   
   !ALLOCATE(rot_stiff(Nc_d))
   


   ! Volume node-pairs associated quantities   
   READ(idf,*); READ(idf,*); READ(idf,*); READ(idf,*) 
   READ(idf,*); READ(idf,*); READ(idf,*); READ(idf,*) 

   DO c = 1, NcFE_G
     
     IF (.NOT. MP_job  .OR.  MP_master) THEN
     
       READ(idf,*) idx, mass(c), stiff(c) 
       READ(idf,*) eta(:,c)
       READ(idf,*) Drr(1,c),  Drr(2,c),  Drr(3,c)
       DO k = 1, k_d
         READ(idf,*) stiff_T(k,:,c)
       ENDDO
     
     ELSE
     
       IF (ANY(cG_cP_FE == c)) THEN
       
         READ(idf,*) idx, mass(cP_cG_FE(c)), stiff(cP_cG_FE(c)) 
         READ(idf,*) eta(:,cP_cG_FE(c))
         READ(idf,*) Drr(1,cP_cG_FE(c)), Drr(2,cP_cG_FE(c)), Drr(3,cP_cG_FE(c))
         DO k = 1, k_d
           READ(idf,*) stiff_T(k,:,cP_cG_FE(c))
         ENDDO       
       
       ELSE
       
         READ(idf,*)
         READ(idf,*); READ(idf,*)
         READ(idf,*)
         DO k = 1, 2*k_d
           READ(idf,*)
         ENDDO       
       
       ENDIF
     
     ENDIF
     
   ENDDO   


   ! Volume nodes associated quantities
   READ(idf,*); READ(idf,*); 
   READ(idf,*); READ(idf,*)

   DO i = 1, Nj_G
   
     IF (.NOT. MP_job  .OR.  MP_master) THEN
     
       READ(idf,*) idx, mass_ii(i), cell(i)
       READ(idf,*) stiff_ii(i)
     
     ELSE
       
       IF (ANY(jG_jP == i)) THEN
         READ(idf,*) idx, mass_ii(jP_jG(i)), cell(jP_jG(i))
         READ(idf,*) stiff_ii(jP_jG(i))
       ELSE
         READ(idf,*)
	 READ(idf,*)
       ENDIF
       
     ENDIF  
   
   ENDDO   
   


   ! Boundary node-pairs associated quantities
   READ(idf,*); READ(idf,*); READ(idf,*);   READ(idf,*)
   READ(idf,*); READ(idf,*); READ(idf,*)

   DO c = 1, NcFE_b
      
      IF (.NOT. MP_job  .OR.  MP_master) THEN
      
        READ(idf,*) idx, mass_b(c) 
        READ(idf,*) chi_b(:,c)
        READ(idf,*) kappa_b(:,:,c)
      
      ELSE

        IF (cP_cG_FE_B(c) /= 0  .AND. NcFEP_b /= 0) THEN
          READ(idf,*) idx, mass_b(cP_cG_FE_B(c))
          READ(idf,*) chi_b(:,cP_cG_FE_B(c))
          READ(idf,*) kappa_b(:,:,cP_cG_FE_B(c))
	ELSE
	  READ(idf,*)
	  READ(idf,*)
	  READ(idf,*)
	  READ(idf,*)
	ENDIF
      
      ENDIF	
   
   ENDDO   



   ! Boundary nodes associated quantities   
   READ(idf,*); READ(idf,*); 
   READ(idf,*) 

   DO i = 1, Njb_G
   
     IF (.NOT. MP_job  .OR.  MP_master) THEN
       READ(idf,*) idx, xi_bp(:,i) 
     ELSE
       
       IF (jP_jG_B(i) /= 0 .AND. NjP_b /= 0) THEN                     
         READ(idf,*) idx, xi_bp(:,jP_jG_B(i))
       ELSE
         READ(idf,*)
       ENDIF
     
     ENDIF
   
   ENDDO   



   ! Volume node-pairs associated quantities
   READ(idf,*);  READ(idf,*);	 
   READ(idf,*)
   
   DO c = 1, NcFV_G
     
     IF (.NOT. MP_job  .OR.  MP_master) THEN
       
       READ(idf,*) idx, eta_fv(:,c)
       READ(idf,*) Drr_fv(:,c)
       
     ELSE
       
       IF (ANY(cG_cP_FV == c)) THEN
         READ(idf,*) idx, eta_fv(:,cP_cG_FV(c))
         READ(idf,*) Drr_fv(:,cP_cG_FV(c))
       ELSE
         READ(idf,*)
	 READ(idf,*)
       ENDIF
     
     ENDIF
       
   ENDDO   


   ! ROTATION MATRICES IN FINITE VOLUMES          
   ! Domain: consines_ETA_FV		          
   READ(idf,*); READ(idf,*)	          

   DO c = 1, NcFV_G	     
   
     IF (.NOT. MP_job  .OR.  MP_master) THEN		          
   
       IF (k_d == 1) THEN		          
         READ(idf,*) idx, cosines_ETA_FV(1,:,c)
       ELSEIF (k_d == 2) THEN		         
         READ(idf,*) idx, cosines_ETA_FV(1,:,c)
         READ(idf,*)	  cosines_ETA_FV(2,:,c)
       ELSE				         
         READ(idf,*) idx, cosines_ETA_FV(1,:,c)
         READ(idf,*)	  cosines_ETA_FV(2,:,c)
         READ(idf,*)	  cosines_ETA_FV(3,:,c)
       ENDIF
       
     ELSE
       
       IF (ANY(cG_cP_FV == c)) THEN
         
	 IF (k_d == 1) THEN			  
           READ(idf,*) idx, cosines_ETA_FV(1,:,cP_cG_FV(c))
         ELSEIF (k_d == 2) THEN 		 
           READ(idf,*) idx, cosines_ETA_FV(1,:,cP_cG_FV(c))
           READ(idf,*)      cosines_ETA_FV(2,:,cP_cG_FV(c))
         ELSE					 
           READ(idf,*) idx, cosines_ETA_FV(1,:,cP_cG_FV(c))
           READ(idf,*)      cosines_ETA_FV(2,:,cP_cG_FV(c))
           READ(idf,*)      cosines_ETA_FV(3,:,cP_cG_FV(c))
         ENDIF
       
       ELSE
       
    	 IF (k_d == 1) THEN			  
           READ(idf,*)
         ELSEIF (k_d == 2) THEN 		 
           READ(idf,*)
           READ(idf,*)
         ELSE					 
           READ(idf,*)
           READ(idf,*)
           READ(idf,*)
         ENDIF
   
       ENDIF
     
     ENDIF  
       
   ENDDO      


   ! Boundary: cosines_XI_BP
   READ(idf,*); READ(idf,*)
   DO i = 1, Njb_G	     
   
     IF (.NOT. MP_job  .OR.  MP_master) THEN		          
   
       IF (k_d == 1) THEN		          
         READ(idf,*) idx, cosines_XI_BP(1,:,i)
       ELSEIF (k_d == 2) THEN		         
         READ(idf,*) idx, cosines_XI_BP(1,:,i)
         READ(idf,*)	  cosines_XI_BP(2,:,i)
       ELSE				         
         READ(idf,*) idx, cosines_XI_BP(1,:,i)
         READ(idf,*)	  cosines_XI_BP(2,:,i)
         READ(idf,*)	  cosines_XI_BP(3,:,i)
       ENDIF
       
     ELSE
       
       IF (jP_jG_B(i) /= 0 .AND. Njp_b /= 0) THEN
         
	 IF (k_d == 1) THEN			  
           READ(idf,*) idx, cosines_XI_BP(1,:,jP_jG_B(i))
         ELSEIF (k_d == 2) THEN 		 
           READ(idf,*) idx, cosines_XI_BP(1,:,jP_jG_B(i))
           READ(idf,*)      cosines_XI_BP(2,:,jP_jG_B(i))
         ELSE					 
           READ(idf,*) idx, cosines_XI_BP(1,:,jP_jG_B(i))
           READ(idf,*)      cosines_XI_BP(2,:,jP_jG_B(i))
           READ(idf,*)      cosines_XI_BP(3,:,jP_jG_B(i))
         ENDIF
       
       ELSE
       
    	 IF (k_d == 1) THEN			  
           READ(idf,*)
         ELSEIF (k_d == 2) THEN 		 
           READ(idf,*)
           READ(idf,*)
         ELSE					 
           READ(idf,*)
           READ(idf,*)
           READ(idf,*)
         ENDIF
   
       ENDIF
     
     ENDIF  
       
   ENDDO      
  
   END SUBROUTINE read_metric
   


   

   SUBROUTINE read_metric_axisymmetric(idf)
   !---------------------------------------------------------------------
   USE mp_interface
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf

   INTEGER :: c, i, idx 
   CHARACTER(LEN=16) :: name
   !---------------------------------------------------------------------

   READ(idf,*) 
   READ(idf,'(16x,a15)') name
   
   !IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Reading metric quantities: ', TRIM(name)

   READ(idf,*); READ(idf,*)
   
   IF (.NOT. MP_job  .OR.  MP_master) THEN
     READ(idf,*) k_d, Nc_d, Np_d, Nc_b, Np_b, Nc_fv
     NcFE_G = Nc_d
     Nj_G = Np_d
     NcFE_b = Nc_b
     Njb_G = Np_b
     NcFV_G = Nc_fv
   ELSE
     READ(idf,*) k_d
     Nc_d = NcTFEP
     Np_d = NjTP
     Nc_b = NcFEP_b
     Np_b = 2*NjP_b
     Nc_fv = NcTFVP
   ENDIF  

   IF (ALLOCATED(mass_y))	DEALLOCATE(mass_y) 
   IF (ALLOCATED(eta_y))	DEALLOCATE(eta_y) 
   IF (ALLOCATED(Drr)) 	        DEALLOCATE(Drr) 
   IF (ALLOCATED(mass_y_ii))	DEALLOCATE(mass_y_ii) 
   IF (ALLOCATED(cell_y))	DEALLOCATE(cell_y) 
   IF (ALLOCATED(xi_y_bp))	DEALLOCATE(xi_y_bp)
   IF (ALLOCATED(eta_y_fv))	DEALLOCATE(eta_y_fv)
   IF (ALLOCATED(Drr_fv))	DEALLOCATE(Drr_fv)

   ALLOCATE(mass_y(Nc_d), eta_y(k_d, Nc_d))
   ALLOCATE(mass_y_ii(Np_d), cell_y(Np_d))
   ALLOCATE(Drr(4*k_d-1,Nc_d))
   ALLOCATE(xi_y_bp(k_d, Np_b))
   ALLOCATE(eta_y_fv(k_d, Nc_fv), Drr_fv(4*k_d-1,Nc_fv))

 

   ! Volume node-pairs associated quantities
   READ(idf,*); READ(idf,*); READ(idf,*); READ(idf,*) 
   READ(idf,*); READ(idf,*); READ(idf,*)
   
   DO c = 1, NcFE_G

     IF (.NOT. MP_job  .OR.  MP_master) THEN
   
       READ(idf,*) idx, mass_y(c)
       READ(idf,*) eta_y(:,c)
       READ(idf,*) Drr(1,c),  Drr(2,c),  Drr(3,c)
       
     ELSE
     
       IF (ANY(cG_cP_FE == c)) THEN
      
         READ(idf,*) idx, mass_y(cP_cG_FE(c))
         READ(idf,*) eta_y(:,cP_cG_FE(c))
         READ(idf,*) Drr(1,cP_cG_FE(c)),  Drr(2,cP_cG_FE(c)),  Drr(3,cP_cG_FE(c))
      
       ELSE
         
	 READ(idf,*)
	 READ(idf,*)
	 READ(idf,*)
	 	        
       ENDIF     
     
     ENDIF  
   
   ENDDO   



   ! Volume nodes associated quantities
   READ(idf,*); READ(idf,*); 
   READ(idf,*)

   DO i = 1, Nj_G

     IF (.NOT. MP_job  .OR.  MP_master) THEN   
       READ(idf,*) idx, mass_y_ii(i), cell_y(i)
     ELSE
     
       IF (ANY(jG_jP == i)) THEN
         READ(idf,*)  idx, mass_y_ii(jP_jG(i)), cell_y(jP_jG(i))
       ELSE
         READ(idf,*)
       ENDIF     
     
     ENDIF
   
   ENDDO   
   
   

   ! Boundary nodes associated quantities
   READ(idf,*); READ(idf,*)
   READ(idf,*); READ(idf,*); READ(idf,*) 

   DO i = 1, Njb_G
 
     IF (.NOT. MP_job  .OR.  MP_master) THEN   
       READ(idf,*) idx, xi_y_bp(:,i) 
     ELSE
       
       IF (jP_jG_B(i) /= 0) THEN
         READ(idf,*) idx, xi_y_bp(:,jP_jG_B(i))
       ELSE
         READ(idf,*)
       ENDIF        
       
     ENDIF
        
   ENDDO
   
   
   
   ! Volume node-pairs associated quantities
   READ(idf,*);  READ(idf,*);  
   READ(idf,*)

   DO c = 1, NcFV_G
   
     IF (.NOT. MP_job  .OR.  MP_master) THEN
      
       READ(idf,*) idx, eta_y_fv(:,c)
       READ(idf,*) Drr_fv(:,c)
     
     ELSE
     
       IF (ANY(cG_cP_FV == c)) THEN
         READ(idf,*) idx, eta_y_fv(:,cP_cG_FV(c))
         READ(idf,*) Drr_fv(:,cP_cG_FV(c))
       ELSE
         READ(idf,*)
	 READ(idf,*)
       ENDIF
     
     ENDIF
     
   ENDDO
  
   END SUBROUTINE read_metric_axisymmetric


END MODULE  metric_coefficients

