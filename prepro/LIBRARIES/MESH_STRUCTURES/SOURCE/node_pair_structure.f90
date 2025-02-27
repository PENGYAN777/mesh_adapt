!============================================================ 
!
!      Module: node_pair_structure
!
! Description: In this module the node-pair data structure is 
!              defined and subroutines to initialize the 
!              node-pair based representation of the mesh, 
!              starting from the usual element/nodes one, are 
!              given.   
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
!
!   Quantities belonging to the domain
!   ----------------------------------
!   
!      Nc_d   Total number of node-pair in the domain
!
!     j_c_d   Connectivity between the node-pair c and its 
!             nodes j 
!
!     c_j_d   Connectivity between the node j and the 
!             node-pairs c it belongs to (DIV format) 
!
!     c_m_d   Connectivity between the element m and its 
!             node-pairs c (DIV format)
!
!
!   Quantities belonging to the boundary
!   ------------------------------------
!   
!      Nc_b   Total number of node-pair in the boundary
!
!     j_c_b   Connectivity between the node-pair c and its 
!             nodes j 
!
!     c_j_b   Connectivity between the node j and the 
!             node-pairs c it belongs to (DIV format) 
!
!     c_m_b   Connectivity between the element m and its 
!             node-pairs c (DIV format)
!   bound_c   Connectivity between the node-pair and the
!             boundaries
!
!============================================================ 

!============================================================ 
!
! Remarks:      
!
!   (1) Connectivity matrix  (other than j_c) are stored in 
!   DIV format to save memory, making therefore more 
!   difficult their manipulation.  Matrices stored in such a 
!   format has been explicitly indicated above, and has to be 
!   used as follows.  
!   For example, the vector c(:) of the indices of 
!   the node-pairs the node j belongs to is given by:  
!                     c(:) = c_j(j) % vec 
!   See the module dynamic_vector for details about the DIV 
!   data structure.
!
!============================================================ 

!============================================================ 
!
!         c1       c      c2
!       o------O======O------o
!      is      i      j      js
!
!============================================================
 
MODULE  node_pair_structure



   !============================================================ 
   USE nodes
   USE mesh_structure
   USE dynamic_vector

   !============================================================ 
   INTEGER  ::  Nc_d
   INTEGER,     DIMENSION(:,:), ALLOCATABLE  ::  j_c_d
   TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE  ::  c_m_d
   INTEGER,     DIMENSION(:,:), ALLOCATABLE  ::  cs_c_d

   INTEGER  ::  Nc_b
   INTEGER,     DIMENSION(:,:), ALLOCATABLE  ::  j_c_b
   TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE  ::  c_m_b
   INTEGER,     DIMENSION(:),   ALLOCATABLE  ::  cd_cb
   INTEGER,     DIMENSION(:),   ALLOCATABLE  ::  bound_c

   INTEGER  ::  Nc_fv
   INTEGER,     DIMENSION(:,:), ALLOCATABLE  ::  j_c_fv
   TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE  ::  c_m_fv
   INTEGER,     DIMENSION(:,:), ALLOCATABLE  ::  cs_c_fv
   
   TARGET  ::  j_c_d, cs_c_d, j_c_fv, cs_c_fv
   !============================================================ 

   CONTAINS


   !============================================================ 
   SUBROUTINE save_node_pair(idf, name, name_length)
   !============================================================ 


      !------------------------------------------------------------ 
      !
      ! Save the node-pair topology to disk
      !
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)  ::  idf
      CHARACTER(LEN=64), INTENT(IN) :: name
      INTEGER, INTENT(IN)  ::  name_length 
      !------------------------------------------------------------ 
      INTEGER  ::  c
      !------------------------------------------------------------ 
    
      WRITE(idf,1000)
      WRITE(idf,1010) name(1:name_length)

      WRITE(idf,1000)    
      WRITE(idf,1020)
      WRITE(idf,1021) k_d, Nc_d, Nc_b, Nc_fv 

      WRITE(idf,1000)
      WRITE(idf,1025)
    
      WRITE(idf,1000)
      WRITE(idf,1040)

      IF (k_d == 1)  WRITE(idf,1045) 
      IF (k_d == 2)  WRITE(idf,1046) 
      IF (k_d == 3)  WRITE(idf,1047) 


      ! Domain node-pairs
      ! -----------------
      DO c = 1, Nc_d
        IF (k_d == 1)  WRITE(idf,1048) c, j_c_d(:,c), cs_c_d(:,c) 
        IF (k_d == 2)  WRITE(idf,1049) c, j_c_d(:,c), cs_c_d(:,c) 
        IF (k_d == 3)  WRITE(idf,1050) c, j_c_d(:,c), cs_c_d(:,c)
      ENDDO   
        
      WRITE(idf,1000)
      WRITE(idf,1060)
    
      WRITE(idf,1000)
      WRITE(idf,1040)

      IF (k_d == 1)  WRITE(idf,1065) 
      IF (k_d == 2)  WRITE(idf,1066) 
      IF (k_d == 3)  WRITE(idf,1067) 

      ! Boundary node-pairs
      ! -------------------
      DO c = 1, Nc_b
        IF (k_d == 1)  WRITE(idf,1068) c, j_c_b(:,c), cd_cb(c), bound_c(c)
        IF (k_d == 2)  WRITE(idf,1069) c, j_c_b(:,c), cd_cb(c), bound_c(c)
        IF (k_d == 3)  WRITE(idf,1070) c, j_c_b(:,c), cd_cb(c), bound_c(c)
      ENDDO   


      ! Finite Volume domain node-pairs
      ! -------------------------------
      WRITE(idf,1000)
      WRITE(idf,1025)
    
      WRITE(idf,1000)
      WRITE(idf,1042)

      IF (k_d == 1)  WRITE(idf,1045) 
      IF (k_d == 2)  WRITE(idf,1046) 
      IF (k_d == 3)  WRITE(idf,1047) 

      DO c = 1, Nc_fv
        IF (k_d == 1)  WRITE(idf,1048) c, j_c_fv(:,c), cs_c_fv(:,c) 
        IF (k_d == 2)  WRITE(idf,1049) c, j_c_fv(:,c), cs_c_fv(:,c) 
        IF (k_d == 3)  WRITE(idf,1050) c, j_c_fv(:,c), cs_c_fv(:,c) 
      ENDDO   
        
  
1000  FORMAT('######################################################################################')
1010  FORMAT('#    NAME:      ',a15,'                                                      #')
1019  FORMAT('#        K_D        NC_D        NC_B                                                 #')
1020  FORMAT('#        K_D        NC_D        NC_B       NC_FV                                     #')
1021  FORMAT(4i12)
1025  FORMAT('#  **********  DOMAIN  **********                                                    #')
1040  FORMAT('#  ++++   NODE-PAIRS   ++++                                                          #')
1042  FORMAT('#  ++++  FV NODE_PAIRS ++++                                                          #')
1045  FORMAT('#        IDX           I           J          I*          J*       C_I*I       C_JJ* #')
1046  FORMAT('#',8x,'IDX', 11x, 'I', 11x, 'J', 11x, 'I*', 9x, 'J*', 9x, 'Itp', 9x, 'Jtp', 9x, 'Itm', 9x, 'Jtm', &
             7x, 'C_I*I       C_JJ* #')
1047  FORMAT('#',8x,'IDX', 11x, 'I', 11x, 'J', 11x, 'I*', 9x, 'J*', 9x, 'Itp', 9x, 'Jtp', 9x, 'Itm', 9x, 'Jtm', &
             9x, 'Ibp', 9x, 'Jbp', 9x, 'Ibm', 9x, 'Jbm', 7x, 'C_I*I       C_JJ* #')
1048  FORMAT(7i12)
1049  FORMAT(11i12)
1050  FORMAT(15i12)
1060  FORMAT('#  **********  BOUNDARY  **********                                                  #')
1065  FORMAT('#        IDX           I           J          I*          J*       C_DOM       BOUND #')
1066  FORMAT('#',8x,'IDX', 11x, 'I', 11x, 'J', 11x, 'I*', 9x, 'J*', 9x, 'Itp', 9x, 'Jtp', 9x, 'Itm', 9x, 'Jtm', &
             7x, 'C_DOM       BOUND #')
1067  FORMAT('#',8x,'IDX', 11x, 'I', 11x, 'J', 11x, 'I*', 9x, 'J*', 9x, 'Itp', 9x, 'Jtp', 9x, 'Itm', 9x, 'Jtm', &
             9x, 'Ibp', 9x, 'Jbp', 9x, 'Ibm', 9x, 'Jbm', 7x, 'C_DOM       BOUND #')
1068  FORMAT(7i12)
1069  FORMAT(11i12)
1070  FORMAT(15i12)

   END SUBROUTINE save_node_pair
   !============================================================ 




   SUBROUTINE read_node_pair(idf)
   !------------------------------------------------------------ 
   USE mp_interface

   !------------------------------------------------------------ 
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: idf

   INTEGER, DIMENSION(:), ALLOCATABLE :: nVec
   INTEGER, DIMENSION(2) :: cs_c
   INTEGER :: c, idx, bc, Bidx, k, p
   
   CHARACTER(LEN=16) :: name
   !------------------------------------------------------------ 
   
   READ(idf,*) 
   READ (idf,'(16x,a15)') name
   !IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Reading node-pair structure: ', TRIM(name)

   READ(idf,*); READ(idf,*)
   
   IF (.NOT. MP_job  .OR.  MP_master) THEN 
     READ(idf,*) k_d, Nc_d, Nc_b, Nc_fv
     NcFE_G = Nc_d
     NcFV_G = Nc_fv
     NcFE_b = Nc_b 
   ELSE
     READ(idf,*) k_d
     Nc_d = NcTFEP
     Nc_fv = NcTFVP
     Nc_b = NcFEP_b
   ENDIF     

   IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   FE Domain node-pairs: ', Nc_d
   IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Boundary node-pairs:  ', Nc_b
   IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   FV Domain node-pairs: ', Nc_fv

   IF (ALLOCATED(j_c_d))    DEALLOCATE (j_c_d)
   IF (ALLOCATED(cs_c_d))   DEALLOCATE (cs_c_d)
   IF (ALLOCATED(j_c_b))    DEALLOCATE (j_c_b)
   IF (ALLOCATED(cd_cb))    DEALLOCATE (cd_cb)
   IF (ALLOCATED(bound_c))  DEALLOCATE (bound_c)
   IF (ALLOCATED(j_c_fv))   DEALLOCATE (j_c_fv)
   IF (ALLOCATED(cs_c_fv))  DEALLOCATE (cs_c_fv)
   
   	 
   ALLOCATE (j_c_d(4*k_d,Nc_d),   cs_c_d(2,Nc_d))
   ALLOCATE (j_c_fv(4*k_d,Nc_fv), cs_c_fv(2,Nc_fv))

   ALLOCATE (j_c_b(4*k_d,Nc_b), cd_cb(Nc_b), bound_c(Nc_b))

   ALLOCATE (nVec(4*k_d))

   ! Domain node-pairs
   READ(idf,*); READ(idf,*); 
   READ(idf,*); READ(idf,*); 
   READ(idf,*) 

   DO c = 1, NcFE_G
     
     IF (.NOT. MP_job  .OR.  MP_master) THEN
       READ(idf,*) idx, j_c_d(:,c), cs_c_d(:,c) 
     ELSE

       READ(idf,*) idx, nVec, cs_c
       
       IF (ANY(cG_cP_FE == idx)) THEN
         
	 DO p = 1, SIZE(nVec)
           IF (nVec(p) /= 0) THEN
	     IF (jP_jG(nVec(p)) == 0) THEN
  	       j_c_d(p,cP_cG_FE(c)) = 1
	     ELSE
  	       j_c_d(p,cP_cG_FE(c)) = jP_jG(nVec(p))
	     ENDIF  
	   ELSE
             j_c_d(p,cP_cG_FE(c)) = 0
	   ENDIF
	 ENDDO  
	  
	 cs_c_d(:,cP_cG_FE(c)) = SIGN(cP_cG_FE(ABS(cs_c)), cs_c)

	 DO k = 1, SIZE(cs_c_d,1), 1
	   IF (cs_c_d(k, cP_cG_FE(c)) == 0) cs_c_d(k, cP_cG_FE(c)) = 1
	 ENDDO

       ENDIF
       
     ENDIF
       
   ENDDO   
   
   ! Boundary node-pairs    
   READ(idf,*); READ(idf,*); 
   READ(idf,*); READ(idf,*); 
   READ(idf,*) 
   
   k = 1
   DO c = 1, NcFE_b

     IF (.NOT. MP_job  .OR.  MP_master) THEN
       READ(idf,*) idx, j_c_b(:,c), cd_cb(c), bound_c(c)   
     ELSE

       READ(idf,*) Bidx, nVec, idx, bc

       IF (ANY(cG_cP_FE == idx)) THEN

         DO p = 1, SIZE(nVec)
           IF (nVec(p) /= 0) THEN
	     IF (jP_jG_B(nVec(p)) == 0) THEN
  	       j_c_b(p,k) = 1
	     ELSE
  	       j_c_b(p,k) = jP_jG_B(nVec(p))
	     ENDIF  	   
	   ELSE
	     j_c_b(p,k) = 0	   
	   ENDIF
	 ENDDO  

	 cd_cb(k) = cP_cG_FE(idx)
	 bound_c(k) = bc
	 cP_cG_FE_B(Bidx) = k
	 k = k + 1

       ENDIF
       
     ENDIF
   
   
   ENDDO   
 
 
   ! Finite Volume domain node-pairs
   READ(idf,*); READ(idf,*); 
   READ(idf,*); READ(idf,*); 
   READ(idf,*) 

   DO c = 1, NcFV_G

     IF (.NOT. MP_job  .OR.  MP_master) THEN	
       READ(idf,*) idx, j_c_fv(:,c), cs_c_fv(:,c)
     ELSE							      

       READ(idf,*) idx, nVec, cs_c

       IF (ANY(cG_cP_FV == idx)) THEN
         
	 DO p = 1, SIZE(nVec)
           IF (nVec(p) /= 0) THEN
	     IF (jP_jG(nVec(p)) == 0) THEN
  	       j_c_fv(p,cP_cG_FV(c)) = 1
	     ELSE
  	       j_c_fv(p,cP_cG_FV(c)) = jP_jG(nVec(p))
	     ENDIF  
	   ELSE
             j_c_fv(p,cP_cG_FV(c)) = 0
	   ENDIF
	 ENDDO  
	 
	 cs_c_fv(:,cP_cG_FV(c)) = SIGN(cP_cG_FV(ABS(cs_c)), cs_c)
	 
	 DO k = 1, SIZE(cs_c_fv,1), 1
	   IF (cs_c_fv(k, cP_cG_FV(c)) == 0) cs_c_fv(k, cP_cG_FV(c)) = 1
	 ENDDO
	 
       ENDIF

     ENDIF							      

   ENDDO

   END SUBROUTINE read_node_pair


END MODULE  node_pair_structure
