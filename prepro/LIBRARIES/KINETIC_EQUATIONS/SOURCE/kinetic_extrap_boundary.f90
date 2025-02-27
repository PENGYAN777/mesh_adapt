!============================================================ 
!
!      Module: kinetic_extrap_boundary
!
! Description: Extrapolation strategies at boundaries.
!              This is an exact copy of the module for 
!              extrapolation at boundaries for the
!              Euler fluxes. Repeated for clarity.
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

MODULE  kinetic_extrap_boundary

 !----------------------------------------------------------------------------
 USE thermodynamics

 INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  ext_type_bound    
 INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  ext_bc
 INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  ext_type_c    
 INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  c_ext    
 INTEGER,      DIMENSION(:),     ALLOCATABLE  ::  j_c__ext
 
 REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  ::  normD
 REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE  ::  ww_ext    

 INTEGER, PARAMETER  ::  ZERO_ORDER    = 0,   &
                         FIRST_ORDER   = 1,   &
                         SECOND_ORDER  = 2,   &
                         MUSCL_EXTRAP  = 3

 CHARACTER(*), DIMENSION(0:3), PARAMETER  ::  ext_names = &
                      (/ 'ZERO ORDER EXTRAPOLATION   ', &
                         'FIRST ORDER EXTRAPOLATION  ', &
                         'SECOND ORDER EXTRAPOLATION ', &
                         'LIMITED MUSCL EXTRAPOLATION' /)
 !----------------------------------------------------------------------------

 CONTAINS
 

   SUBROUTINE  read_param_kin_extrap_boundary(idf)
   !----------------------------------------------------------------------------
   USE mp_interface
   
   IMPLICIT NONE
   
   INTEGER :: idf
   INTEGER :: Nbound, b
   !----------------------------------------------------------------------------

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
   
   END SUBROUTINE  read_param_kin_extrap_boundary





   SUBROUTINE  init_kin_extrap_boundary(ww, j_c_d, cs_c_d, bound_p, Drr)
   !----------------------------------------------------------------------------
   USE mp_interface
   USE nodes,   ONLY: k_d, jd_jb

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: ww
   INTEGER,      DIMENSION(:,:), INTENT(INOUT) :: j_c_d
   INTEGER,      DIMENSION(:,:), INTENT(INOUT) :: cs_c_d
   INTEGER,      DIMENSION(:),   INTENT(INOUT) :: bound_p
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: Drr

   INTEGER, DIMENSION(SIZE(ww,2)):: jb_jd

   INTEGER :: i, j, c, c_e, Np_d, Nc_d, Nc_ext  
   !----------------------------------------------------------------------------

   Np_d = SIZE(ww,2)
   Nc_d = SIZE(j_c_d,2)

   ! Computes the number of extensions. For each
   ! node-pair (c) counts the number of node-pairs
   ! for which boundary extrapolation is required
   Nc_ext = 0

   DO c = 1, SIZE(j_c_d,2)

      ! Normal
      IF (j_c_d(4,c) == j_c_d(1,c))  Nc_ext = Nc_ext + 1

      ! Tangent
      IF (k_d >= 2) THEN
        IF (j_c_d(5,c) == j_c_d(1,c))  Nc_ext = Nc_ext + 1
        IF (j_c_d(6,c) == j_c_d(2,c))  Nc_ext = Nc_ext + 1
        IF (j_c_d(7,c) == j_c_d(1,c))  Nc_ext = Nc_ext + 1
        IF (j_c_d(8,c) == j_c_d(2,c))  Nc_ext = Nc_ext + 1
      ENDIF 

      ! Binormal
      IF (k_d == 3) THEN
        IF (j_c_d(9,c)  == j_c_d(1,c))  Nc_ext = Nc_ext + 1
        IF (j_c_d(10,c) == j_c_d(2,c))  Nc_ext = Nc_ext + 1
        IF (j_c_d(11,c) == j_c_d(1,c))  Nc_ext = Nc_ext + 1
        IF (j_c_d(12,c) == j_c_d(2,c))  Nc_ext = Nc_ext + 1
      ENDIF

   ENDDO

   ALLOCATE (c_ext(Nc_ext), j_c__ext(Nc_ext), ww_ext(SIZE(ww,1), Nc_ext))


   ! Set the extensions' indices. To each node-pair extended outside the boundary
   ! associates the index of the node-pair (c) inside the domain
   Nc_ext = 0

   DO c = 1, SIZE(j_c_d,2)

      IF (j_c_d(4,c) == j_c_d(1,c)) THEN
        Nc_ext = Nc_ext + 1
        c_ext(Nc_ext) = c
      ENDIF


      IF (k_d >= 2) THEN

        IF (j_c_d(5,c) == j_c_d(1,c)  .OR.  &
            j_c_d(7,c) == j_c_d(1,c)) THEN

          Nc_ext = Nc_ext + 1
          c_ext(Nc_ext) = c

        ENDIF

        IF (j_c_d(6,c) == j_c_d(2,c)  .OR.  & 
            j_c_d(8,c) == j_c_d(2,c)) THEN 

          Nc_ext = Nc_ext + 1
          c_ext(Nc_ext) = c

        ENDIF

      ENDIF 

      IF (k_d == 3) THEN

        IF (j_c_d(9,c)  == j_c_d(1,c)  .OR.  &
            j_c_d(11,c) == j_c_d(1,c)) THEN 

          Nc_ext = Nc_ext + 1
          c_ext(Nc_ext) = c

        ENDIF

        IF (j_c_d(10,c) == j_c_d(2,c)  .OR.  &
            j_c_d(12,c) == j_c_d(2,c))  THEN

          Nc_ext = Nc_ext + 1
          c_ext(Nc_ext) = c

        ENDIF

      ENDIF 

   ENDDO
   
   
   ALLOCATE (ext_type_c(SIZE(c_ext)), &
             ext_bc(SIZE(c_ext)))

   ! Invert jd_jb
   jb_jd = 0
   DO i = 1, SIZE(jd_jb)
     jb_jd(jd_jb(i)) = i 
   ENDDO      

   ! For each "ouside-extended" node-pair sets the 
   ! order of extrapolation
   c_e = 1
   
   DO WHILE (c_e <= SIZE(c_ext))

     c = c_ext(c_e)

     ! NORMAL DIRECTION ------------------------------------
     IF (j_c_d(4,c) == j_c_d(1,c)) THEN
       j = j_c_d(2,c)  
       ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
       c_e = c_e + 1
       IF (c_e > SIZE(c_ext)) EXIT
       c = c_ext(c_e)
     ENDIF


     ! TANGENT DIRECTION ------------------------------------
     IF (k_d >= 2) THEN

       IF (j_c_d(5,c) == j_c_d(1,c)) THEN
         j = j_c_d(5,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

       IF (j_c_d(6,c) == j_c_d(2,c)) THEN
         j = j_c_d(6,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

       IF (j_c_d(7,c) == j_c_d(1,c)) THEN
         j = j_c_d(7,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

       IF (j_c_d(8,c) == j_c_d(2,c)) THEN
         j = j_c_d(8,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

     ENDIF 



     ! BINORMAL DIRECTION ----------------------------------
     IF (k_d == 3) THEN

       IF (j_c_d(9,c)  == j_c_d(1,c)) THEN
         j = j_c_d(9,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

       IF (j_c_d(10,c) == j_c_d(2,c)) THEN
         j = j_c_d(10,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

       IF (j_c_d(11,c) == j_c_d(1,c)) THEN
         j = j_c_d(11,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

       IF (j_c_d(12,c) == j_c_d(2,c)) THEN
         j = j_c_d(12,c)
         ext_type_c(c_e) = ext_type_bound(bound_p(jb_jd(j)))
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
         c = c_ext(c_e)
       ENDIF

     ENDIF 

   ENDDO


   ! Sets the indices of the nodes extended 
   ! outside the boundary
   c_e = 1
   DO WHILE (c_e <= SIZE(c_ext))

     c = c_ext(c_e)

     ! NORMAL DIRECTION ------------------------------------
     IF (j_c_d(4,c) == j_c_d(1,c)) THEN
       j_c_d(4,c) = Np_d + c_e
       j_c__ext(c_e) = Np_d + c_e
       cs_c_d(2,c) = Nc_d + c_e
       Drr(3,c) = ABS(Drr(3,c))
       c_e = c_e + 1
       IF (c_e > SIZE(c_ext)) EXIT
     ENDIF


     ! TANGENT DIRECTION ------------------------------------
     IF (k_d >= 2) THEN

       IF (j_c_d(5,c) == j_c_d(1,c)) THEN 
         j_c_d(5,c) = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(4,c) = Drr(6,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

       IF (j_c_d(6,c) == j_c_d(2,c)) THEN 
         j_c_d(6,c) = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(5,c) = Drr(7,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

       IF (j_c_d(7,c) == j_c_d(1,c)) THEN 
         j_c_d(7,c) = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(6,c) = Drr(4,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

       IF (j_c_d(8,c) == j_c_d(2,c)) THEN 
         j_c_d(8,c) = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(7,c) = Drr(5,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

     ENDIF 
     
     
     ! BINORMAL DIRECTION ----------------------------------
     IF (k_d == 3) THEN

       IF (j_c_d(9,c)  == j_c_d(1,c)) THEN 
         j_c_d(9,c)  = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(8,c) = Drr(10,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

       IF (j_c_d(10,c) == j_c_d(2,c)) THEN 
         j_c_d(10,c) = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(9,c) = Drr(11,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

       IF (j_c_d(11,c) == j_c_d(1,c)) THEN 
         j_c_d(11,c) = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(10,c) = Drr(8,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

       IF (j_c_d(12,c) == j_c_d(2,c)) THEN 
         j_c_d(12,c) = Np_d + c_e 
         j_c__ext(c_e) = Np_d + c_e
         Drr(11,c) = Drr(9,c)
         c_e = c_e + 1
         IF (c_e > SIZE(c_ext)) EXIT
       ENDIF

     ENDIF 

   ENDDO

   END SUBROUTINE  init_kin_extrap_boundary





   SUBROUTINE  write_param_kin_extrap_boundary(idf)
   !----------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER  ::  idf
   !----------------------------------------------------------------------------
   INTEGER  ::  b


   WRITE(idf,*) ' Parameters for high order extrapolation'
   WRITE(idf,*) ' at boundaries'
   WRITE(idf,*) ' Number of boundaries:          ', SIZE(ext_type_bound)
   DO b = 1, SIZE(ext_type_bound)
      WRITE(idf,*) '   Extrapolation of type: ', ext_names(ext_type_bound(b))
   ENDDO
   WRITE(idf,*)

   END SUBROUTINE  write_param_kin_extrap_boundary





   FUNCTION  ww_kin_ext_size()   RESULT(ext_size)
   !----------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER :: ext_size
   !----------------------------------------------------------------------------

   IF (ALLOCATED(ww_ext)) THEN
     ext_size = SIZE(ww_ext,2)
   ELSE
     ext_size = 0
   ENDIF

   END FUNCTION  ww_kin_ext_size



   FUNCTION  ww_kin_extrap(ww, j_c_d, Drr)   RESULT (ww_ext)
   !----------------------------------------------------------------------------
   USE np_quadrature
   USE nodes,  ONLY: k_d

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
   INTEGER,      DIMENSION(:,:), INTENT(IN) :: j_c_d
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Drr

   REAL(KIND=8), DIMENSION(SIZE(ww,1),SIZE(c_ext)) :: ww_ext

   REAL(KIND=8) :: Drr_ij, Drr_iis, Drr_jjs, &
                   Drr_iitp, Drr_iitm, Drr_jjtp, Drr_jjtm, &
		   Drr_iibp, Drr_iibm, Drr_jjbp, Drr_jjbm

   INTEGER :: i, j, is, js, itp, jtp, itm, jtm, &
              ibp, jbp, ibm, jbm, c, c_e, j_ext
   !----------------------------------------------------------------------------

   DO c_e = 1, SIZE(c_ext)
   
      c = c_ext(c_e)
      j_ext = j_c__ext(c_e)

      i  = j_c_d(1,c);  j  = j_c_d(2,c)
      is = j_c_d(3,c);  js = j_c_d(4,c)
      
      Drr_ij = Drr(1,c); Drr_iis = Drr(2,c); Drr_jjs = Drr(3,c)

      IF (k_d >= 2) THEN
        itp = j_c_d(5,c);    itm = j_c_d(7,c)
        Drr_iitp = Drr(4,c); Drr_iitm = Drr(6,c)
        jtp = j_c_d(6,c);    jtm = j_c_d(8,c)
        Drr_jjtp = Drr(5,c); Drr_jjtm = Drr(7,c)
      ENDIF

      IF (k_d == 3) THEN
        ibp = j_c_d(9,c);    ibm = j_c_d(11,c)
        Drr_iibp = Drr(8,c); Drr_iibm = Drr(10,c)
        jbp = j_c_d(10,c);   jbm = j_c_d(12,c)
        Drr_jjbp = Drr(9,c); Drr_jjbm = Drr(11,c)
      ENDIF


      IF (js == j_ext) THEN

        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:, c_e) = ww(:,j)
             
          CASE (FIRST_ORDER)
	  ww_ext(:,c_e) = 2*ww(:,j) - ww(:,i)

          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d >= 2  .AND.  itp == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:,c_e) = ww(:,i)
             
          CASE (FIRST_ORDER)
          ww_ext(:,c_e) = ww(:,i) + (ww(:,i) - ww(:,itm)) * Drr_iitp/Drr_iitm
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d >= 2  .AND.  jtp == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:,c_e) = ww(:,j)
             
          CASE (FIRST_ORDER)
          ww_ext(:,c_e) = ww(:,j) + (ww(:,j) - ww(:,jtm)) * Drr_jjtp/Drr_jjtm
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d >= 2  .AND.  itm == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:,c_e) = ww(:,i)
             
          CASE (FIRST_ORDER)	     
          ww_ext(:,c_e) = ww(:,i) + (ww(:,i) - ww(:,itp)) * Drr_iitm/Drr_iitp
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d >= 2  .AND.  jtm == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:,c_e) = ww(:,j)
             
          CASE (FIRST_ORDER)	     
          ww_ext(:,c_e) = ww(:,j) + (ww(:,j) - ww(:,jtp)) * Drr_jjtm/Drr_jjtp
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d == 3  .AND.  ibp == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:,c_e) = ww(:,i)
             
          CASE (FIRST_ORDER)
          ww_ext(:,c_e) = ww(:,i) + (ww(:,i) - ww(:,ibm)) * Drr_iibp/Drr_iibm
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d == 3  .AND.  jbp == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:,c_e) = ww(:,j)
             
          CASE (FIRST_ORDER)
          ww_ext(:,c_e) = ww(:,j) + (ww(:,j) - ww(:,jbm)) * Drr_jjbp/Drr_jjbm
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d == 3  .AND.  ibm == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:,c_e) = ww(:,i)
             
          CASE (FIRST_ORDER)
          ww_ext(:,c_e) = ww(:,i) + (ww(:,i) - ww(:,ibp)) * Drr_iibm/Drr_iibp
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

        CYCLE


      ELSEIF (k_d == 3  .AND.  jbm == j_ext) THEN


        SELECT CASE (ext_type_c(c_e))

          CASE (ZERO_ORDER)
          ww_ext(:, c_e) = ww(:,j)
             
          CASE (FIRST_ORDER)
          ww_ext(:, c_e) = ww(:,j) + (ww(:,j) - ww(:,jbp)) * Drr_jjbm/Drr_jjbp
          
          CASE DEFAULT
          WRITE(*,*) ' Extrapolation of unknown type. STOP'; STOP

        END SELECT 

      ENDIF

   ENDDO

   END FUNCTION ww_kin_extrap


   END MODULE  kinetic_extrap_boundary
