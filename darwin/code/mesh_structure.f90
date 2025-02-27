!============================================================ 
!
!    Author:  Alberto Guardone
!             Dipartimento di Ingegneria Aerospaziale
!             Politecnico di Milano
!             Via La Msa 34, 20156 Milano, ITALY
!             e-mail: guardone@aero.polimi.it
!
! Copyright: 1998-2003 Alberto Guardone
!            See COPYING file for copyright notice
!
!============================================================ 

MODULE  mesh_structure

!============================================================ 
!   This module provide the definition of the variables in 
!   which the connectivity and the geometry of the mesh are 
!   stored. An IBRID mesh is assumed, and hence the 
!   connectivity matrices are stored in DIV format to save
!   memory.  Both the data structures and subroutines defined
!   in this modules can deal with one, two or three 
!   dimensional meshes.
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
!      bound_m   Connectivity between elements and the part of 
!                the boundary they belongs to
!    
!   ele_type_b   Type of boundary element
!    
!
!============================================================ 

!============================================================ 
!
! NOTES:   
! ------   
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

                  DEALLOCATE( faces_mb(fb)%vec )    !(S)
               
               ENDDO

               DEALLOCATE( faces_mb )    !(S)
               
            ENDDO

            DEALLOCATE( faces(f)%vec )    !(S)
         
         ENDDO

         DEALLOCATE( faces )    !(S)
            
      ENDDO



   END FUNCTION ma_m_connectivity
   !============================================================ 



   !============================================================ 
   SUBROUTINE del_double_nodes_boundary (jd_jb_DD,   &
                                         jb_jd_DD,   &
                                         j_m_b_DD,   &
                                         bound_p_DD, &
                                         bound_m_DD, &
                                         jd_jb_SS,   &
                                         j_m_b_SS,   &
                                         bound_p_SS)

   ! Heavily modifyed by marco fossati. 25 Sep 2006
   !============================================================ 

      IMPLICIT NONE

      INTEGER,     DIMENSION(:),   INTENT(IN)  :: jd_jb_DD
      INTEGER,     DIMENSION(:,:), INTENT(IN)  :: jb_jd_DD
      TYPE(D_I_V), DIMENSION(:),   INTENT(IN)  :: j_m_b_DD  
      INTEGER,     DIMENSION(:),   INTENT(IN)  :: bound_p_DD, &
                                                  bound_m_DD
      INTEGER,     DIMENSION(:), POINTER  :: jd_jb_SS
      TYPE(D_I_V), DIMENSION(:), POINTER  :: j_m_b_SS
      INTEGER,     DIMENSION(:), POINTER  :: bound_p_SS


      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: dn

      INTEGER :: Np_bound, jd, jb, jb_, mb, b 


      ALLOCATE (dn(MAXVAL(jd_jb_DD), MAXVAL(bound_p_DD)))

      ! Computes the number of independent boundary
      ! nodes and creates the new numbering in vector dn
      dn = 0;  Np_bound = 0

      DO jb = 1, SIZE(jd_jb_DD)

         jd = jd_jb_DD(jb)  
      
         b = bound_p_DD(jb)

         IF (dn(jd,b) == 0) THEN

            Np_bound = Np_bound + 1

            dn(jd, b) = Np_bound


            ! Verifies the fact that the node jd is at the
            ! conjunction of two different  boundary lines
            ! (marco fossati. Sept 22, 2006)
            ! old version   'IF (connect_boundaries) dn(jd,:) = Np_bound'
            IF (bound_p_DD(jb_jd_DD(1, jd))  /=   bound_p_DD(jb_jd_DD(2, jd)) ) THEN

              dn(jd, bound_p_DD(jb_jd_DD(1, jd))) = Np_bound
              dn(jd, bound_p_DD(jb_jd_DD(2, jd))) = Np_bound

            ENDIF

         ENDIF

      ENDDO

      ! Generate the boundary node --> domain node 
      ! connectivity jd_jb and update the 
      ! boundary element --> boundary node connecitvity j_m_b
      ALLOCATE (jd_jb_SS(Np_bound),       &
                j_m_b_SS(SIZE(j_m_b_DD)), &
                bound_p_SS(Np_bound))

      DO mb = 1, SIZE(j_m_b_SS)

         ALLOCATE (j_m_b_SS(mb) % vec( SIZE(j_m_b_DD(mb)%vec) ))

         DO jb_ = 1, SIZE(j_m_b_DD(mb)%vec)

            jb = j_m_b_DD(mb)%vec(jb_)   
            
            jd = jd_jb_DD(jb) 
             b = bound_p_DD(jb) 

            j_m_b_SS(mb)%vec(jb_) = dn(jd,b)  

            jd_jb_SS(dn(jd,b)) = jd

            bound_p_SS(j_m_b_SS(mb)%vec(jb_)) = bound_m_DD(mb)

         ENDDO

      ENDDO
      
   END SUBROUTINE del_double_nodes_boundary
   !============================================================ 



   !============================================================ 
   SUBROUTINE read_mesh (idf)
   !============================================================ 


      IMPLICIT NONE

      INTEGER, INTENT(IN)  ::  idf

      INTEGER  ::  m   ! Generic element index      

      INTEGER  ::  idx

      CHARACTER (LEN=16)  ::  name


      READ (idf,*)  
      READ (idf,'(16x,a15)') name 

      !WRITE(*,*) '   Reading mesh:  ',name 

      READ (idf,*); READ (idf,*)

      READ (idf,*) Ne_d, Ne_b 


      !WRITE(*,*) '  Domain elements:   ',Ne_d
      !WRITE(*,*) '  Boundary elements: ',Ne_b 
          
      IF ( ALLOCATED(j_m_d) )      DEALLOCATE (j_m_d)
      IF ( ALLOCATED(ma_m_d) )     DEALLOCATE (ma_m_d)
      IF ( ALLOCATED(ele_type_d) ) DEALLOCATE (ele_type_d)

      IF ( ALLOCATED(j_m_b) )      DEALLOCATE (j_m_b)
      IF ( ALLOCATED(ma_m_b) )     DEALLOCATE (ma_m_b)    
      IF ( ALLOCATED(bound_m) )    DEALLOCATE (bound_m)      
      IF ( ALLOCATED(ele_type_b) ) DEALLOCATE (ele_type_b)
 

      ALLOCATE ( j_m_d(Ne_d),  ma_m_d(Ne_d)  )
      ALLOCATE ( ele_type_d(Ne_d) )

      ALLOCATE ( j_m_b(Ne_b),  ma_m_b(Ne_b)  )
      ALLOCATE ( bound_m(Ne_b) )
      ALLOCATE ( ele_type_b(Ne_b) )


      READ (idf,*); READ (idf,*); READ (idf,*)
      READ (idf,*); READ (idf,*); READ (idf,*)
      READ (idf,*)


      ! Element of the domain
      ! ---------------------
      DO m = 1, Ne_d

          READ(idf,*) idx, ele_type_d(m)
          
          IF ( m .NE. idx ) THEN
              WRITE(*,*) ' Warning, grid numeration not ordered '
              WRITE(*,*) ' Using internal (sequential) numeration '
          ENDIF

          ALLOCATE ( j_m_d(m) % vec(N_nodes_ele(ele_type_d(m))) )         
          READ(idf,*) j_m_d(m) % vec
          
          ALLOCATE ( ma_m_d(m) % vec(N_faces_ele(ele_type_d(m))) )
          READ(idf,*)  ma_m_d(m) % vec

      ENDDO


      READ (idf,*); READ (idf,*); READ (idf,*) 
      READ (idf,*); READ (idf,*); READ (idf,*)
      READ (idf,*)

      ! Element of the boundary
      ! -----------------------
      DO m = 1, Ne_b

          READ(idf,*) idx, ele_type_b(m), bound_m(m)

          IF ( m .NE. idx ) THEN
              WRITE(*,*) ' Warning, grid numeration not ordered '
              WRITE(*,*) ' Using internal (sequential) numeration '
          ENDIF
 
          ALLOCATE ( j_m_b(m) % vec(N_nodes_ele(ele_type_b(m))) )         
          READ(idf,*)  j_m_b(m) % vec

          ALLOCATE ( ma_m_b(m) % vec(N_faces_ele(ele_type_b(m))) )
          READ(idf,*) ma_m_b(m) % vec

      ENDDO


   END SUBROUTINE read_mesh
   !============================================================ 



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
