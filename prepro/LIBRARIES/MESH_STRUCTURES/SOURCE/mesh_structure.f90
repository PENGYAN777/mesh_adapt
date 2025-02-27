!============================================================ 
!
!      Module: mesh_structures
!
! Description: This module provide the definition of the 
!              variables in which the connectivity and the 
!              geometry of the mesh are stored. An HYBRID 
!              mesh is assumed, and hence the connectivity 
!              matrices are stored in DIV format to save
!              memory.  Both the data structures and 
!              subroutines defined in this modules can deal 
!              with one, two or three dimensional meshes. 
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
!   Domain elements
!   ---------------
!
!         Ne_d   Total number of elements
!
!        j_m_d   Connectivity between the element m and its 
!                nodes j (DIV format)
!
!        m_j_d   Connectivity between the node j and the 
!                elements m it belongs to (DIV format)
!
!       ma_m_d   Connectivity between the element m and the 
!                elements ma adjacent to the element faces 
!                (DIV format)
!
!   ele_type_d   Type of domain element
!
!
!   Boundary elements
!   -----------------
!
!         Ne_b   Total number of elements
!
!        j_m_b   Connectivity between the element m and its 
!                nodes j (DIV format)
!        m_j_b   Connectivity between the node j and the 
!                elements m it belongs to (DIV format)
!
!       ma_m_b   Connectivity between the element m and the 
!                elements ma adjacent to element m's faces 
!                (DIV format)
!
!      bound_m   Connectivity between nodes and the part of 
!                the boundary they belongs to
!    
!   ele_type_b   Type of boundary element
!    
!
!============================================================ 

!============================================================ 
!
! Remarks:   
!
!   (1) Connectivity matrix are stored in DIV format to 
!   save memory, making therefore more difficult their 
!   manipulation.  Matrices stored in such a format has been 
!   explicitly indicated above, and has to be used as 
!   follows.  For example, the vector (of variable size since 
!   we are dealing with ibrid meshes) of indices of the nodes 
!   belonging to the element m is given by the statement:  
!                     j(:) = j_m(m) % vec   
!   Similarly, the vector m(:) of the indices of the elements 
!   the node j belongs to is given by:  
!                     m(:) = m_j(j) % vec 
!   See the module dynamic_vector for details about the DIV 
!   data strucutre.
!
!   (2) The element type is defined in module 
!   element_topology
!
!============================================================ 

MODULE  mesh_structure


   !============================================================ 
   USE dynamic_vector
   USE element_topology
   USE nodes
   !============================================================ 


   !============================================================ 
   INTEGER  ::  Ne_d
   TYPE (D_I_V), DIMENSION(:), ALLOCATABLE  ::  j_m_d
   TYPE (D_I_V), DIMENSION(:), ALLOCATABLE  ::  m_j_d
   TYPE (D_I_V), DIMENSION(:), ALLOCATABLE  ::  ma_m_d 
   INTEGER,      DIMENSION(:), ALLOCATABLE  ::  ele_type_d

   INTEGER  ::  Ne_b
   TYPE (D_I_V), DIMENSION(:), ALLOCATABLE  ::  j_m_b
   TYPE (D_I_V), DIMENSION(:), ALLOCATABLE  ::  m_j_b
   TYPE (D_I_V), DIMENSION(:), ALLOCATABLE  ::  ma_m_b 
   INTEGER,      DIMENSION(:), ALLOCATABLE  ::  bound_m
   INTEGER,      DIMENSION(:), ALLOCATABLE  ::  ele_type_b
   !============================================================ 



 !============================================================ 
 CONTAINS
 !============================================================ 



   !============================================================ 
   FUNCTION md_mb_connectivity(j_m_d, j_m_b, jd_jb) RESULT (md_mb)
   !============================================================ 


      IMPLICIT NONE

      TYPE(D_I_V), DIMENSION(:), INTENT(IN)  ::  j_m_d, j_m_b 
      INTEGER,     DIMENSION(:), INTENT(IN)  ::  jd_jb
      INTEGER,     DIMENSION(:), POINTER     ::  md_mb
      
      TYPE(D_I_V), DIMENSION(:), ALLOCATABLE  ::  m_j_d
      INTEGER  ::  mb, md, md_, i, i_, j_, j, n
      LOGICAL  ::  found_ele, found_i


      ALLOCATE( m_j_d(size_DIV(j_m_d,3)) )
      ALLOCATE( md_mb(SIZE(j_m_b)) )

      md_mb = 0

      m_j_d = invert_DIV(j_m_d)

      DO mb = 1, SIZE(j_m_b)
                                       ! Bubble of one of the nodes of mb
         n =  jd_jb(j_m_b(mb)%vec(1))  ! (node 1 is always present) 

         DO md_ = 1, SIZE(m_j_d(n)%vec);  md = m_j_d(n)%vec(md_)
         found_ele = .true.           

            DO i_ = 1, SIZE(j_m_b(mb)%vec);  i = jd_jb(j_m_b(mb)%vec(i_))

               found_i = .false.

               DO j_ = 1, SIZE(j_m_d(md)%vec);  j = j_m_d(md)%vec(j_)
                  IF ( i == j ) THEN 
                     found_i = .true.
                  ENDIF
               ENDDO

               found_ele = found_ele .AND. found_i

            ENDDO              

            IF (found_ele) THEN 
               md_mb(mb) = md 
            ENDIF
   
         ENDDO

      ENDDO


   END FUNCTION md_mb_connectivity
   !============================================================ 



   !============================================================ 
   FUNCTION ma_m_connectivity(j_m, ele_type) RESULT (ma_m)
   !============================================================ 

      ! 1D and 2D only

      IMPLICIT NONE

      TYPE(D_I_V), DIMENSION(:), INTENT(IN)  ::  j_m
      INTEGER,     DIMENSION(:), INTENT(IN)  ::  ele_type
      TYPE(D_I_V), DIMENSION(:), POINTER     ::  ma_m
      
      TYPE(D_I_V), DIMENSION(:), ALLOCATABLE  ::  m_j
      INTEGER  ::  mb_, mb, m, i, i_, j_, j, f, fb, p
      LOGICAL  ::  found_ele, found_i
      TYPE(D_I_V), DIMENSION(:), POINTER ::  faces, faces_mb
      INTEGER, DIMENSION(:), ALLOCATABLE  ::  j_f, j_fb


      ALLOCATE( m_j(size_DIV(j_m,3)) )
      ALLOCATE( ma_m(SIZE(j_m)) )

      m_j = invert_DIV(j_m)

      DO m = 1, SIZE(j_m)
         ALLOCATE( ma_m(m)%vec(N_faces_ele(ele_type(m))) )
         ma_m(m)%vec = 0
      ENDDO

      DO m = 1, SIZE(j_m)
                             
         faces => ele_faces(ele_type(m))                     

         DO f = 1, SIZE(faces)
            
            IF ( ALLOCATED(j_f) ) DEALLOCATE (j_f)
            ALLOCATE( j_f(SIZE(faces(f)%vec)) ) 
            j_f = j_m(m)%vec(faces(f)%vec)
            
            p = j_f(1)  ! Node 1 always present
            
            
            DO mb_ = 1, SIZE(m_j(p)%vec);  mb = m_j(p)%vec(mb_)
               IF (mb == m) CYCLE
            
               faces_mb => ele_faces(ele_type(mb))
            
               DO fb = 1, SIZE(faces_mb)
            found_ele = .TRUE.           
            
                  IF ( ALLOCATED(j_fb) ) DEALLOCATE (j_fb)
                  ALLOCATE( j_fb(SIZE(faces_mb(fb)%vec)) ) 
            
                  j_fb = j_m(mb)%vec(faces_mb(fb)%vec)
                  
                  DO i_ = 1, SIZE(j_f); i = j_f(i_)
                   
                     found_i = .FALSE.      
                     DO j_ = 1, SIZE(j_fb);  j = j_fb(j_)
                        IF ( i == j) THEN
                           found_i = .TRUE.
                        ENDIF
                     ENDDO
                     found_ele = found_ele .AND. found_i
                  ENDDO
               
                
                  IF (found_ele) THEN 
                      ma_m(m) %vec(f)  = mb
                      ma_m(mb)%vec(fb) = m 
                  ENDIF
               
               ENDDO
               
            ENDDO
         
         ENDDO
            
      ENDDO



   END FUNCTION ma_m_connectivity
   !============================================================ 



   !============================================================ 
   SUBROUTINE del_double_nodes_boundary (jd_jb_double, bound_p, &
                                         j_m_b, bound_m,        &
                                         jd_jb_single,          &
                                         connect_boundaries     )
   !============================================================ 


      IMPLICIT NONE

      INTEGER,     DIMENSION(:), INTENT(IN)     ::  jd_jb_double
      INTEGER,     DIMENSION(:), INTENT(IN)     ::  bound_p, bound_m
      TYPE(D_I_V), DIMENSION(:), INTENT(INOUT)  ::  j_m_b
      INTEGER,     DIMENSION(:), POINTER        ::  jd_jb_single
      LOGICAL                                   ::  connect_boundaries

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: dn

      INTEGER :: Np_b, jd, jb, jb_, mb, b 


      ALLOCATE (dn(MAXVAL(jd_jb_double),MAXVAL(bound_p)))

      ! Computes the number of independent boundary
      ! nodes and creates the new numbering in vector dn
      dn = 0;  Np_b = 0
      DO jb = 1, SIZE(jd_jb_double); jd = jd_jb_double(jb)  
      
         b = bound_p(jb) 

         IF (dn(jd,b) == 0) THEN

            Np_b = Np_b + 1

            dn(jd, b) = Np_b

            IF(connect_boundaries) dn(jd,:) = Np_b

         ENDIF

      ENDDO

      ! Generate the boundary node --> domain node 
      ! connectivity jd_jb and update the 
      ! boundary element --> boundary node connecitvity j_m_b
      ALLOCATE (jd_jb_single(Np_b))
      DO mb = 1, SIZE(j_m_b)

         DO jb_ = 1, SIZE(j_m_b(mb)%vec);  jb = j_m_b(mb)%vec(jb_)   
            
            jd = jd_jb_double(jb) 
             b = bound_p(jb) 

            j_m_b(mb)%vec(jb_) = dn(jd,b)  

            jd_jb_single(dn(jd,b)) = jd     

         ENDDO

      ENDDO
      
      CALL set_bound_p(j_m_b, bound_m)
      
   END SUBROUTINE del_double_nodes_boundary
   !============================================================ 





   SUBROUTINE read_mesh (idf)
   !-----------------------------------------------------------
   USE mp_interface
   
   !-----------------------------------------------------------
   IMPLICIT NONE

   INTEGER, INTENT(IN)  ::  idf

   INTEGER, DIMENSION(:), ALLOCATABLE :: nVec, eVec, mbNodes

   INTEGER :: m, k   ! Generic element index
   INTEGER :: idx, typeD, typeB, idxP, bm
   INTEGER :: jb1, jb2
   
   LOGICAL, DIMENSION(Nm_G) :: mbFlag

   CHARACTER (LEN=16) :: name
   !-----------------------------------------------------------

   READ (idf,*)  
   READ (idf,'(16x,a15)') name 

   !IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Reading mesh: ', TRIM(name)

   READ (idf,*); READ (idf,*)

   IF (.NOT. MP_job  .OR.  MP_master) THEN
     READ (idf,*) Ne_d, Ne_b 
     Nm_G  = Ne_d
     Nmb_G = Ne_b
   ELSE
     READ(idf,*)
     Ne_d = NmP
     Ne_b = 0
   ENDIF  

   IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Domain elements: ', Ne_d 
   IF (.NOT. MP_job  .OR.  MP_master)  WRITE(*,*) '   Boundary elements: ', Ne_b 


   IF (ALLOCATED(j_m_d))      DEALLOCATE (j_m_d)
   IF (ALLOCATED(ma_m_d))     DEALLOCATE (ma_m_d)
   IF (ALLOCATED(ele_type_d)) DEALLOCATE (ele_type_d)

   IF (ALLOCATED(j_m_b))      DEALLOCATE (j_m_b)
   IF (ALLOCATED(ma_m_b))     DEALLOCATE (ma_m_b)    
   IF (ALLOCATED(bound_m))    DEALLOCATE (bound_m)	 
   IF (ALLOCATED(ele_type_b)) DEALLOCATE (ele_type_b)
 

   ALLOCATE (j_m_d(Ne_d), ma_m_d(Ne_d))
   ALLOCATE (ele_type_d(Ne_d))

   IF (.NOT. MP_job  .OR.  MP_master) THEN
     ALLOCATE (j_m_b(Ne_b), ma_m_b(Ne_b))
     ALLOCATE (bound_m(Ne_b))
     ALLOCATE (ele_type_b(Ne_b))
   ENDIF  


   ! Element of the domain
   READ (idf,*); READ (idf,*); READ (idf,*)
   READ (idf,*); READ (idf,*); READ (idf,*)
   READ (idf,*)

   DO m = 1, Nm_G   

     READ(idf,*) idx, typeD
     
     IF (m .NE. idx) THEN
       WRITE(*,*) 'WARNING. Grid numeration not ordered '
       WRITE(*,*) 'Using internal (sequential) numeration '
     ENDIF

   
     IF (.NOT. MP_job  .OR.  MP_master) THEN
     
       ele_type_d(m) = typeD
     
       ALLOCATE (j_m_d(m) % vec(N_nodes_ele(ele_type_d(m))))
       READ(idf,*) j_m_d(m) % vec
     
       ALLOCATE (ma_m_d(m) % vec(N_faces_ele(ele_type_d(m))))
       READ(idf,*)  ma_m_d(m) % vec


     ELSE
     
       IF (ANY(mG_mP == idx)) THEN

         idxP = mP_mG(idx)
         ele_type_d(idxP) = typeD

         ALLOCATE (j_m_d(idxP) % vec(N_nodes_ele(ele_type_d(idxP))))
         ALLOCATE (nVec(N_nodes_ele(ele_type_d(idxP))))

         READ(idf,*) nVec
         j_m_d(idxP) % vec = jP_jG(nVec) 
     
         ALLOCATE (ma_m_d(idxP) % vec(N_faces_ele(ele_type_d(idxP))))
         ALLOCATE (eVec(N_nodes_ele(ele_type_d(idxP))))
	 
         READ(idf,*) eVec
         ma_m_d(idxP) % vec = mP_mG(eVec)

         DEALLOCATE (nVec, eVec)

       ELSE
	 READ(idf,*)
	 READ(idf,*)
       ENDIF
	      
     ENDIF

   ENDDO


   ! Element of the boundary
   READ (idf,*); READ (idf,*); READ (idf,*) 
   READ (idf,*); READ (idf,*); READ (idf,*)
   READ (idf,*)

   NmbP = 0
   DO m = 1, Nmb_G

     READ(idf,*) idx, typeB, bm

     IF (m .NE. idx) THEN
       WRITE(*,*) ' Warning, grid numeration not ordered '
       WRITE(*,*) ' Using internal (sequential) numeration '
     ENDIF
     
     IF (.NOT. MP_job  .OR.  MP_master) THEN
     
       ele_type_B(m) = typeB
       bound_m(m) = bm
     
       ALLOCATE (j_m_b(m) % vec(N_nodes_ele(ele_type_b(m))))	     
       READ(idf,*)  j_m_b(m) % vec

       ALLOCATE (ma_m_b(m) % vec(N_faces_ele(ele_type_b(m))))
       READ(idf,*) ma_m_b(m) % vec
     
     ELSE
     
       IF (ALLOCATED(mbNodes))  DEALLOCATE(mbNodes)
       ALLOCATE (mbNodes(N_nodes_ele(typeB)))
       READ(idf,*) mbNodes
       READ(idf,*)
       
       ! Count the number boundary element for each partition
       IF (ALL(jP_jG_B(mbNodes) /= 0)) THEN 
         IF(ANY(jd_jb(jP_jG_B(mbNodes)) <= NjP)) THEN
       
           NmbP = NmbP + 1
	   mbFlag(m) = .TRUE.
	   mP_mG_B(idx) = NmbP

         ENDIF
       ENDIF	 
          
     ENDIF
      
   ENDDO


    ! Allocate and read  
   IF (MP_job  .AND.  MP_worker) THEN
   
     Ne_b = NmbP
   
     ALLOCATE (j_m_b(Ne_b), ma_m_b(Ne_b))
     ALLOCATE (bound_m(Ne_b))
     ALLOCATE (ele_type_b(Ne_b))

     ! Returns to the top of boundary elements' list
     DO m = 1, 3*Nmb_G
       BACKSPACE idf
     ENDDO
     
     k = 1
     DO m = 1, Nmb_G
       IF (mbFlag(m)) THEN
         
	 READ(idf,*) idx, ele_type_b(k), bound_m(k)
	 
         ALLOCATE (j_m_b(k) % vec(N_nodes_ele(ele_type_b(k))))	     
         
	 IF (ALLOCATED(mbNodes))  DEALLOCATE(mbNodes)
         ALLOCATE (mbNodes(N_nodes_ele(ele_type_b(k))))
         READ(idf,*) mbNodes
  	 
	 j_m_b(k) % vec = jP_jG_B(mbNodes)

         ALLOCATE (ma_m_b(k) % vec(N_faces_ele(ele_type_b(k))))
	 
	 IF (ALLOCATED(mbNodes))  DEALLOCATE(mbNodes)
         ALLOCATE (mbNodes(N_nodes_ele(ele_type_b(k))))
         READ(idf,*) mbNodes

         ma_m_b(k) % vec = mP_mG_B(mbNodes)
	 
	 k = k + 1
	 
       ELSE
         READ(idf,*)	 
         READ(idf,*)	 
         READ(idf,*)	 
       ENDIF
     ENDDO
   
   ENDIF
   
   IF (k_d == 2) THEN
  
     ALLOCATE (jdm_jb(2,Np_b))
  
     ! Two Dimensional only !!!!!!!
     ! Kinetic outflow boundary condition
     DO m = 1, Ne_b
     
       jb1 = j_m_b(m) % vec(1)
       jb2 = j_m_b(m) % vec(2)

       jdm_jb(1,jb1) = jd_jb(jb1)
       jdm_jb(2,jb1) = jd_jb(jb2)

       jdm_jb(1,jb2) = jd_jb(jb2)
       jdm_jb(2,jb2) = jd_jb(jb1)

     ENDDO

   ENDIF

   END SUBROUTINE read_mesh




   !============================================================ 
   SUBROUTINE save_mesh (idf, name, name_length)
   !============================================================ 


      IMPLICIT NONE

      INTEGER,            INTENT(IN) :: idf
      CHARACTER (LEN=64), INTENT(IN) :: name
      INTEGER,            INTENT(IN) :: name_length 
 
      INTEGER  ::  m   ! Generic element index      
      INTEGER  ::  f_  ! Generic element's face index in local coordinates  
      INTEGER  ::  j_  ! Generic node index in local coordinates  

    

      WRITE (idf,1000)
      WRITE (idf,1010) name(1:name_length)

      WRITE (idf,1000)
      WRITE (idf,1020)
      WRITE (idf,1021)  Ne_d, Ne_b

      WRITE (idf,1000)
      WRITE (idf,1025)

      WRITE (idf,1000)
      WRITE (idf,1040)
      WRITE (idf,1041) 
      WRITE (idf,1042) 
      WRITE (idf,1043) 

      ! Element of the domain
      ! ---------------------
      DO m=1,Ne_d
         WRITE (idf,1046) m, ele_type_d(m)
         DO j_ = 1, SIZE( j_m_d(m)%vec ) 
            WRITE (UNIT=idf,FMT=1047,ADVANCE='NO')  j_m_d(m)%vec(j_)
         ENDDO
         WRITE(idf,*)
         DO f_ = 1, SIZE( ma_m_d(m)%vec ) 
            WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') ma_m_d(m)%vec(f_)
         ENDDO
         WRITE(idf,*)
      ENDDO


      WRITE (idf,1000)
      WRITE (idf,1055)

      WRITE (idf,1000)
      WRITE (idf,1040)     
      WRITE (idf,1065)
      WRITE (idf,1042) 
      WRITE (idf,1043) 

      ! Element of the boundary
      ! -----------------------
      DO m=1,Ne_b
         WRITE (idf,1066) m, ele_type_b(m), bound_m(m)
         DO j_ = 1, SIZE( j_m_b(m)%vec ) 
            WRITE (UNIT=idf,FMT=1047,ADVANCE='NO')  j_m_b(m)%vec(j_)
         ENDDO
         WRITE(idf,*)
         DO f_ = 1, SIZE( ma_m_b(m)%vec ) 
            WRITE (UNIT=idf,FMT=1047,ADVANCE='NO') ma_m_b(m)%vec(f_)
         ENDDO
         WRITE(idf,*)
      ENDDO


1000  FORMAT('###########################################################################')
1010  FORMAT('#    NAME:      ',a15,'                                           #')
1020  FORMAT('#       NE_D        NE_B                                                  #')
1021  FORMAT(2i12)
1025  FORMAT('#  **********  DOMAIN  **********                                         #')
1040  FORMAT('#  ++++  ELEMENTS ++++                                                    #')
1041  FORMAT('#        IDX        TYPE                                                  #')
1042  FORMAT('#        J_M                                                              #')
1043  FORMAT('#       MA_M                                                              #')
1046  FORMAT(2i12)
1047  FORMAT(i12)
1055  FORMAT('#  **********  BOUNDARY  **********                                       #')
1065  FORMAT('#        IDX        TYPE       BOUND                                      #')
1066  FORMAT(3i12)


   END SUBROUTINE save_mesh
   !============================================================ 



END MODULE  mesh_structure
