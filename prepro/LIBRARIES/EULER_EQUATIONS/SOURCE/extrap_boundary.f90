!============================================================ 
!
!      Module: extrap_boundary
!
! Description: extrapolation strategies at boundaries
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

MODULE extrap_boundary


   !============================================================ 
   USE euler_equations
   USE roe_average
   USE thermodynamics
   !============================================================ 


   !============================================================    
   INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  ext_type_bound    
   INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  ext_type_c    
   INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  c_ext    
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  ::  ww_ext    
   !------------------------------------------------------------ 
   INTEGER, PARAMETER  ::  ZERO_ORDER    = 0,   &
                           FIRST_ORDER   = 1,   &
                           SECOND_ORDER  = 2,   &
                           MUSCL_EXTRAP  = 3
   CHARACTER(*), DIMENSION(0:3), PARAMETER  ::  ext_names = &
                        (/ 'ZERO ORDER EXTRAPOLATION   ', &
                           'FIRST ORDER EXTRAPOLATION  ', &
                           'SECOND ORDER EXTRAPOLATION ', &
                           'LIMITED MUSCL EXTRAPOLATION' /)
   !============================================================ 



 !============================================================ 
 CONTAINS
 !============================================================ 

 
           
   !============================================================ 
   FUNCTION ww_ext_size() RESULT(ext_size)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  ::  ext_size
      !------------------------------------------------------------


      IF (ALLOCATED(ww_ext)) THEN
        ext_size = SIZE(ww_ext,2)
      ELSE
        ext_size = 0
      ENDIF

   
   END FUNCTION ww_ext_size
   !============================================================ 


   !============================================================ 
   SUBROUTINE read_param_extrap_boundary(idf)
   !============================================================ 


      !------------------------------------------------------------
      USE mp_interface
      
      IMPLICIT NONE

      INTEGER :: idf      
      INTEGER :: Nbound, b
      !------------------------------------------------------------

      IF (.NOT. MP_job  .OR.  MP_master) THEN
        READ(idf,*) Nbound
	Nb_G = Nbound
      ELSE
        READ(idf,*)
	Nbound = NbP
      ENDIF
      
      ALLOCATE (ext_type_bound(Nbound))
      
      DO b = 1, Nb_G
        
	IF (.NOT. MP_job  .OR.  MP_master) THEN
          READ(idf,*) ext_type_bound(b)
        ELSE
	
	  IF (bP_bG(b) /= 0) THEN
	    READ(idf,*) ext_type_bound(bP_bG(b))
	  ELSE
	    READ(idf,*)
	  ENDIF
	
	ENDIF
      
      ENDDO

   END SUBROUTINE read_param_extrap_boundary
   !============================================================ 



   !============================================================ 
   SUBROUTINE write_param_extrap_boundary(idf)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE

      INTEGER  ::  idf
      !------------------------------------------------------------
      INTEGER  ::  b


      WRITE(idf,*) '   PARAMETERS FOR HIGH ORDER EXTRAPOLATION AT BOUNDARIES'
      WRITE(idf,*) '   Number of boundaries: ', SIZE(ext_type_bound)
      DO b = 1, SIZE(ext_type_bound)
         WRITE(idf,*) '   Extrapolation of type: ', ext_names(ext_type_bound(b))
      ENDDO
      WRITE(idf,*)


   END SUBROUTINE write_param_extrap_boundary
   !============================================================ 



   !============================================================ 
   SUBROUTINE init_extrap_boundary(ww, j_c_d, cs_c_d, jd_jb, bound_p, Drr, fem_metric)
   !============================================================ 

      USE mp_interface
      
      !------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     :: ww
      INTEGER,      DIMENSION(:,:), INTENT(INOUT)  :: j_c_d
      INTEGER,      DIMENSION(:,:), INTENT(INOUT)  :: cs_c_d
      INTEGER,      DIMENSION(:),   INTENT(INOUT)  :: jd_jb
      INTEGER,      DIMENSION(:),   INTENT(INOUT)  :: bound_p
      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT)  :: Drr
      LOGICAL,                      INTENT(IN)     :: fem_metric
      !------------------------------------------------------------
      INTEGER, DIMENSION(SIZE(ww,2)) :: jb_jd
      INTEGER :: Nc_ext, c, c_e, i, j, Np_d, Nc_d


      Np_d = SIZE(ww,2);  Nc_d = SIZE(j_c_d,2)

      ! Computes the number of extensions
      Nc_ext = 0
      DO c = 1, SIZE(j_c_d,2)
        IF (j_c_d(4,c) == j_c_d(1,c))   Nc_ext = Nc_ext + 1
      ENDDO

      ALLOCATE (c_ext(Nc_ext), ww_ext(SIZE(ww,1), Nc_ext))

      ! Vector c_ext contains the indices
      ! of the node-pairs that require extension 
      ! at the bounary
      Nc_ext = 0
      DO c = 1, SIZE(j_c_d,2)

         IF (j_c_d(4,c) == j_c_d(1,c)) THEN
           Nc_ext = Nc_ext + 1
           c_ext(Nc_ext) = c
         ENDIF

      ENDDO

      ! Assign a progressive index
      ! to the extended nodes 
      DO c_e = 1, SIZE(c_ext)

         c = c_ext(c_e)

          j_c_d(4,c) = Np_d + c_e
         cs_c_d(2,c) = Nc_d + c_e
         Drr(3,c) = ABS(Drr(3,c))

      ENDDO
      
      
      ALLOCATE (ext_type_c(SIZE(c_ext)))

      jb_jd = 0      
      DO i = 1, SIZE(jd_jb)
        jb_jd(jd_jb(i)) = i
      ENDDO      
      
      DO c_e = 1, SIZE(c_ext)

        c = c_ext(c_e)
        j = j_c_d(2,c)

        IF (.NOT.MP_job  .OR.  (     fem_metric  .AND. c <= NcFEP)  &
	                 .OR.  (.NOT.fem_metric  .AND. c <= NcFVP)) THEN
	  ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
	  
        ELSE
	  ext_type_c(c_e) = 0
	ENDIF
	
      ENDDO
       
 
   END SUBROUTINE init_extrap_boundary
   !============================================================ 




   !============================================================ 
   FUNCTION ww_extrap( ww, j_c_d, Drr )  RESULT (ww_ext)
   !============================================================ 

      USE np_quadrature
      USE muscl_fluxes

      !------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  ww
      INTEGER,      DIMENSION(:,:), INTENT(IN)  ::  j_c_d
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  ::  Drr

      REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(c_ext))  ::  ww_ext
      !------------------------------------------------------------
      REAL(KIND=8) ::  Drr_ij, Drr_iis, Drr_jjs
      REAL(KIND=8), DIMENSION(SIZE(ww,1))    ::  Dww_ij, Dww_iis, s_i
      INTEGER  ::  p, i, j, is, c, c_e 
      !------------------------------------------------------------



      DO c_e = 1, SIZE(c_ext)
      
         c = c_ext(c_e)

         i  = j_c_d(1,c);  j  = j_c_d(2,c) 
         is = j_c_d(3,c) 

         Drr_ij = Drr(1,c);  Drr_iis = Drr(2,c);  Drr_jjs = Drr(3,c)
        
         SELECT CASE (ext_type_c(c_e))

            CASE (ZERO_ORDER)
            !-----------------
               ww_ext(:, c_e) = ww(:,j) 
              
            CASE (FIRST_ORDER)
            !-----------------
               ww_ext(:, c_e) = 2*ww(:,j) - ww(:,i) 
            
            CASE (SECOND_ORDER)
            !-----------------
               ww_ext(:, c_e) = ww(:,j)                                  &
                              - (Drr_jjs/Drr_ij) * (ww(:,j) - ww(:,i))   &
                              + ((Drr_ij + Drr_jjs)/(Drr_iis + Drr_ij))  &
                                * ((ww(:,j) - ww(:,i))/Drr_ij +          &
                                   (ww(:,i) - ww(:,is))/Drr_iis) 
        
            CASE (MUSCL_EXTRAP)
            !-----------------

               Dww_ij  =                      ww(:,j) - ww(:,i)  
               Dww_iis = (Drr_ij/Drr_iis) * ( ww(:,i) - ww(:,is) )
 
               DO p = 1, SIZE(ww(:,i))
                  s_i(p) = s_limiter(Dww_iis(p), Dww_ij(p), 1)
               ENDDO 

               ww_ext(:,c_e) = ww(:,j) + 0.5 * s_i &
                               * ( Dww_ij - onemk_2*(Dww_ij - Dww_iis) )    
            
            CASE DEFAULT
            !-----------------

            WRITE(*,*) ' Extrapolation of unknown type. STOP'
            STOP

         END SELECT
                      
      ENDDO


   END FUNCTION ww_extrap
   !============================================================ 


   !============================================================ 
   SUBROUTINE ww_extrap_roe(j_c_d, cs_c_d, wwe, eose,  wwe_, eose_)
   !============================================================ 


      !------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,            DIMENSION(:,:), INTENT(IN)     ::  j_c_d
      INTEGER,            DIMENSION(:,:), INTENT(IN)     ::  cs_c_d
      REAL(KIND=8),       DIMENSION(:,:), INTENT(IN)     ::  wwe
      TYPE(eos_type),     DIMENSION(:),   INTENT(IN)     ::  eose
      REAL(KIND=8),       DIMENSION(:,:), INTENT(INOUT)  ::  wwe_
      TYPE(eos_ext_type), DIMENSION(:),   INTENT(INOUT)  ::  eose_

      INTEGER  ::  j, js, c, c_e, cs 
      !------------------------------------------------------------



      DO c_e = 1, SIZE(c_ext)
      
         c = c_ext(c_e)

         j  = j_c_d(2,c);  js = j_c_d(4,c) 
           
         cs = cs_c_d(2,c)

         CALL intermediate_ww( wwe (:,j),  eose (j),  &
                               wwe (:,js), eose (js), &
                               wwe_(:,cs), eose_(cs))

      ENDDO


   END SUBROUTINE ww_extrap_roe
   !============================================================ 

END MODULE extrap_boundary
