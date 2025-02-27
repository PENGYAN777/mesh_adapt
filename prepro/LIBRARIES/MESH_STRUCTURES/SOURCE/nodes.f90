!============================================================ 
!
!      Module: nodes
!
! Description: This module provide the data strucuture in 
!              which the geometry of the mesh is stored, and 
!              the procedures to read it and save it. 
!              Geometry is defined through the coordinates of 
!              the nodes and the connectivity between 
!              boundary nodes and domain nodes.
!              Both the data structures and subroutines 
!              defined in this modules can deal with one, 
!              two or three dimensional meshes 
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

!============================================================
!   This module provide the data strucuture in which the 
!   geometry of the mesh is stored, and the procedures to 
!   read it and save it. 
!   Geometry is defined through the coordinates of the nodes
!   and the connectivity between boundary nodes and domain
!   nodes.
!   Both the data structures and subroutines defined
!   in this modules can deal with one, two or three 
!   dimensional meshes
!============================================================

!============================================================
!
!     k_d     physical dimension of the space
!
!     Np_d    Total number of nodes of the domain
!
!     Np_b    Total number of nodes of the boundary
!    
!  rr(k,i)    k-th coordinate of the i-th node
!
!     jd_jb   Connectivity between the node jb in the 
!             boundary numeration and the corresponding node 
!             jd in the numeration of the domain 
!             [ jd = jd_jb(jb)]
!
!   bound_p   Connectivity between the node p and the 
!             boundary portion
!
!============================================================

MODULE  nodes


   !============================================================
   USE dynamic_vector
   !============================================================

   !============================================================
   INTEGER  ::  k_d, Np_d, Np_b
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  ::  rr
   INTEGER,      DIMENSION(:,:), ALLOCATABLE  ::  jdm_jb 
   INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  jd_jb
   INTEGER,      DIMENSION(:,:), ALLOCATABLE  ::  jb_jd
   INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  bound_p
   !============================================================

   CONTAINS


   SUBROUTINE  invert_jd_jb(jd_jb)
   !------------------------------------------------------------------------------------- 
   IMPLICIT NONE                                                                          
  
   INTEGER, DIMENSION(:), INTENT(IN)  :: jd_jb

   INTEGER :: dmn_idx, i                                                                  
   !------------------------------------------------------------------------------------- 

   ALLOCATE (jb_jd(2,Np_d))
   jb_jd = 0

   DO i = 1, Np_b

     dmn_idx = jd_jb(i)                                                            

     IF (jb_jd(1,dmn_idx) .EQ. 0) THEN                                           

       jb_jd(1,dmn_idx) = i                                                        

     ELSE                                                                                 

       jb_jd(2,dmn_idx) = i                                                        

     ENDIF                                                                                

   ENDDO                                                                                  

   END SUBROUTINE  invert_jd_jb                                                           




   !============================================================
   SUBROUTINE set_bound_p(j_m_b, bound_m)
   !============================================================


      IMPLICIT NONE

      TYPE(D_I_V), DIMENSION(:), INTENT(IN) :: j_m_b
      INTEGER,     DIMENSION(:), INTENT(IN) :: bound_m

      INTEGER   ::  m


      IF (ALLOCATED(bound_p)) DEALLOCATE(bound_p)
      ALLOCATE (bound_p(Np_b))

      DO m = 1, SIZE(j_m_b)
        bound_p(j_m_b(m)%vec) = bound_m(m)
      ENDDO


   END SUBROUTINE set_bound_p
   !============================================================




   
   SUBROUTINE  read_nodes(idf)
   !------------------------------------------------------------
   USE mp_interface
   
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf

   INTEGER :: j, k, Bidx, idx, bp
   CHARACTER(LEN=16) :: name
   !------------------------------------------------------------

   READ(idf,*)  
   READ(idf,'(16x,a15)') name 

   !IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Reading nodes: ', TRIM(name)

   READ(idf,*); READ (idf,*)

   IF (.NOT. MP_job  .OR.  MP_master) THEN
     READ(idf,*) k_d, Np_d, Np_b
     Nj_G  = Np_d
     Njb_G = Np_b
     NjP = Np_d
   ELSE
     READ(idf,*) k_d
     Np_d = NjTP
     Np_b = 2*NjP_b
   ENDIF  

   !IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Spatial dimension(s): ', k_d
   IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Domain nodes: ', Np_d
   IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Boundary nodes: ', Np_b
       
   IF (ALLOCATED(rr))  DEALLOCATE (rr)	    
   IF (ALLOCATED(jd_jb))  DEALLOCATE (jd_jb)
   IF (ALLOCATED(bound_p))  DEALLOCATE (bound_p)

   ALLOCATE (rr(k_d,Np_d))
   ALLOCATE (jd_jb(Np_b), bound_p(Np_b))


   ! Volume nodes' coordinates 
   READ(idf,*); READ(idf,*); READ(idf,*)
   READ(idf,*); READ(idf,*); READ(idf,*)

   DO j = 1, Nj_G

       READ(idf,*) idx

       IF (j .NE. idx) THEN
   	 WRITE(*,*) 'WARNING. Nodes numeration not ordered '
   	 WRITE(*,*) 'Using internal (sequential) numeration '
       ENDIF

       IF (.NOT. MP_job  .OR. MP_master) THEN
         READ(idf,*) rr(:,j)
       ELSE
         IF (ANY(jG_jP == idx)) THEN
           READ(idf,*) rr(:,jP_jG(idx))
	 ELSE
	   READ(idf,*)
	 ENDIF
       ENDIF 

   ENDDO


   ! Boundary nodes - domain nodes connectivity
   READ(idf,*); READ(idf,*); READ(idf,*)
   READ(idf,*); READ(idf,*)

   k = 1

   DO j = 1, Njb_G

     IF (.NOT. MP_job  .OR.  MP_master) THEN
  
       READ(idf,*) idx, jd_jb(j), bound_p(j)

       IF (j .NE. idx) THEN
   	 WRITE(*,*) 'WARNING. Nodes numeration not ordered '
   	 WRITE(*,*) 'Using internal (sequential) numeration '
       ENDIF
   
     ELSE

       READ(idf,*) Bidx, idx, bp

       IF (ANY(jG_jP == idx)) THEN
         jd_jb(k) = jP_jG(idx)
	 bound_p(k) = bP_bG(bp)
	 jP_jG_B(Bidx) = k
	 k = k + 1
       ENDIF
       
     ENDIF
   
       
   ENDDO

   CALL invert_jd_jb(jd_jb)

   END SUBROUTINE read_nodes





   !============================================================
   SUBROUTINE save_nodes (idf, name, name_length)
   !============================================================


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idf
      CHARACTER (LEN=64), INTENT(IN) :: name
      INTEGER ::  name_length 
 
      INTEGER :: i
    

      WRITE (idf,1000)
      WRITE (idf,1010) name(1:name_length)

      WRITE (idf,1000)
      WRITE (idf,1020)
      WRITE (idf,1021) k_d, Np_d, Np_b 

      WRITE (idf,1000)
      WRITE (idf,1025)

      WRITE (idf,1000)
      WRITE (idf,1048)
      WRITE (idf,1050)
      WRITE (idf,1051)

      ! Volume nodes coordinates
      ! ------------------------
      DO i=1,Np_d
         WRITE (idf,1052) i
         WRITE (idf,1053) rr(:,i)
      ENDDO

      WRITE (idf,1000)
      WRITE (idf,1055)

      WRITE (idf,1000)
      WRITE (idf,1048)     
      WRITE (idf,1075)

      ! Boundary nodes - domain nodes connectivity
      ! ------------------------------------------
      DO i = 1, Np_b
	 WRITE (idf,1076) i, jd_jb(i), bound_p(i)
      ENDDO
      

1000  FORMAT('###########################################################################')
1010  FORMAT('#    NAME:      ',a15,'                                           #')
1020  FORMAT('#         ND        NP_D        NP_B                                      #')
1021  FORMAT(5i12)
1025  FORMAT('#  **********  DOMAIN  **********                                         #')
1048  FORMAT('#  ++++   NODES   ++++                                                    #')
1050  FORMAT('#        IDX                                                              #')
1051  FORMAT('#         RR                                                              #')
1052  FORMAT(i12)
1053  FORMAT(3e18.9)
1055  FORMAT('#  **********  BOUNDARY  **********                                       #')
1075  FORMAT('#        IDX       JD_JB       BOUND                                      #')
1076  FORMAT(3i12)


   END SUBROUTINE save_nodes
   !============================================================



END MODULE  nodes
