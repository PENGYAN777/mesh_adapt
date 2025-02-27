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

MODULE  np_topology_gen

!============================================================ 
!   In this module the node-pair data structure is defined
!   and subroutines to initialize the node-pair based 
!   representation of the mesh, starting from the usual 
!   element/nodes one, are given.  
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
!   Quantities belonging to the domain
!   ----------------------------------
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
!
!============================================================ 

!============================================================ 
!
! REMARKS:   
! --------   
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
!   EXTENDED NODE-PAIR STRUCTURE:         
!   -----------------------------
!
!       o------O======O------o
!      i*      i      j      j*
!
!       3      1      2      4    <-- Local indices
!
!============================================================ 



   !============================================================ 
   USE node_pair_structure
   USE mesh_structure
   USE nodes
   USE element_topology
   USE dynamic_vector
   !============================================================ 
   


!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   SUBROUTINE node_pair_str_gen
   !============================================================ 


      IMPLICIT NONE

      INTEGER, DIMENSION(:), ALLOCATABLE  ::  fv_fe_np

      !============================================================ 
      ! Domain node-pairs
      ! ===================

      ! Connectivity matrix initialization
      ! ----------------------------------
      ALLOCATE ( c_m_d(Ne_d) )
      c_m_d = c_m_connectivity( j_m_d )

      Nc_d = size_DIV(c_m_d,3)
      !WRITE(*,*) '  I have found ',Nc_d,' node-pairs in the domain'

      ! Generation of the node-pairs structure
      ! --------------------------------------
      ALLOCATE ( j_c_d(4,Nc_d) )
      j_c_d = j_c_connectivity( j_m_d, c_m_d )   ! Node-pairs

      ALLOCATE( cs_c_d(2,Nc_d) )
      CALL extend_node_pair (j_c_d, cs_c_d, rr) ! Extended node-pairs
      IF ( Nc_D /= 0) THEN
      CALL set_extend_bou (j_c_d, cs_c_d, jd_jb, rr) ! Nullify extension on boundaries  
                          ! and reorder the node-pair outward
      ENDIF
      
      ! Selection of the node-pair to
      ! be kept in the FV approach
      ! -----------------------------
      ALLOCATE( fv_fe_np(Nc_d) )
      fv_fe_np = fv_fe_np_connectivity(j_c_d, j_m_d, c_m_d, ele_type_d, Nc_fv)

      ALLOCATE( j_c_fv(4,Nc_fv) ) 
      j_c_fv  = fv_node_pair( j_c_d, fv_fe_np, Nc_fv )
      
      ALLOCATE( cs_c_fv(2,Nc_d) )
      CALL extend_node_pair (j_c_fv, cs_c_fv, rr) ! Extended node-pairs
      IF ( Nc_d /= 0) THEN
      CALL set_extend_bou (j_c_fv, cs_c_fv, jd_jb, rr) ! Nullify extension on boundaries  
                          ! and reorder the node-pair outward
      ENDIF
      
      ALLOCATE ( c_m_fv(Ne_d) )
      c_m_fv  = fv_c_m_connectivity( c_m_d, fv_fe_np  )
      DEALLOCATE (fv_fe_np)
      !============================================================ 



      !============================================================ 
      ! Boundary node-pairs
      ! ===================

      ! Connectivity matrix initialization
      ! ----------------------------------
      IF (Ne_b > 0) THEN

      ALLOCATE ( c_m_b(Ne_b) )

      c_m_b = c_m_connectivity( j_m_b )

      Nc_b = size_DIV( c_m_b, 3 )
      !WRITE(*,*) '  I have found ',Nc_b,' node-pairs in the boundary'

      ! Generation of the node-pairs structure
      ! --------------------------------------
      ALLOCATE ( j_c_b(4,Nc_b) )

      j_c_b = j_c_connectivity( j_m_b, c_m_b )   ! Node-pairs
      
      ! No extension is needed on boundary node-pair
      ! CALL extend_node_pair (j_c_b, rr(:,jd_jb))     
      ! otherwise 
      j_c_b(3:4,:) = 0

      ! Boundary node-pair to domain node-pair connectivity.
      ! It may also change the order of the nodes in boundary 
      ! node-pairs to fit the orientation of the corresponding 
      ! domain node-pairs
      ! ------------------------------------------------------
      ALLOCATE ( cd_cb(Nc_b) )
      cd_cb = cd_cb_connectivity( j_c_d, jd_jb,  j_c_b )

      ! Node-pair to boundary patch connectivity
      ! ----------------------------------------
      ALLOCATE ( bound_c(Nc_b) )
      bound_c = bound_connectivity( c_m_b, bound_m, Nc_b )

      ELSE   ! One dimension

         Nc_b = 0

      ENDIF
      !============================================================ 


   END SUBROUTINE node_pair_str_gen
   !============================================================ 



   !============================================================ 
   FUNCTION c_m_connectivity (j_m) RESULT (c_m)
   !============================================================ 


      !------------------------------------------------------------ 
      !
      ! Starting from the element to node connectivity j_m, 
      ! retrieves the element to node-pair connectivity c_m         
      ! 
      ! NOTE: j_m and c_m are DIV vectors
      !
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      IMPLICIT NONE
      !------------------------------------------------------------ 
      TYPE(D_I_V), DIMENSION(:),         INTENT(IN)   ::  j_m

      TYPE(D_I_V), DIMENSION(SIZE(j_m))               ::  c_m
      !------------------------------------------------------------ 
      INTEGER  ::  Nc   ! Total number of node-pair
      TYPE(D_I_V), DIMENSION(:), ALLOCATABLE  ::  m_j ! node to element 
                                                      ! connectivity 
      INTEGER  ::  m,      &  ! General element 
                   i_, j_, &  ! Nodes in local element indexes
                   i,  j,  &  ! Nodes in local element indexes
                   c_,     &  ! Node-pair in local element coordinates
                   Nc_        ! Total number of node-pair in the element
                       
      ! Quantities pertaining the bubble of node i or j
      ! -----------------------------------------------
      INTEGER  ::  m_B         ! Element belonging to the bubble 
                               ! of node i or j
      INTEGER  ::  m_B_      
      INTEGER  ::  i_B_, j_B_  ! Nodes of element m_B in local indexes
      INTEGER  ::  i_B, j_B    ! Nodes of element m_B in local indexes
      INTEGER  ::  c_B_        ! Node-pair of element m_B
      !------------------------------------------------------------ 



      !------------------------------------------------------------ 
      ! Initialization of the connectivity matrix c_m
      ! ---------------------------------------------
      !
      ! For every element, the local connectivity is 
      ! set.  The total number of local node-pair is
      ! computed and the local connectivity vector
      ! is allocated 
           
      DO m = 1, SIZE(j_m)

         Nc_ = 0 

         ! Loop over all the couples ij to count them
         ! ------------------------------------------
         DO i_ = 1, SIZE( j_m(m)%vec )       
            DO j_ = i_+1, SIZE( j_m(m)%vec )  
               Nc_ = Nc_ + 1            
            ENDDO
         ENDDO

         ALLOCATE ( c_m(m)%vec(Nc_) )
         c_m(m)%vec = 0

      ENDDO
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Create the node to element connectivity 
      ! m_j inverting j_m
      ! ---------------------------------------
      ALLOCATE ( m_j(size_DIV(j_m,3)) )
      m_j = invert_DIV( j_m )
      !------------------------------------------------------------ 


      !============================================================ 
      ! Compute c_m 
      ! ===========
      
      ! The total number of node-pairs cn is set to zero
      ! and the loop on elements begins
      ! ------------------------------------------------

      Nc = 0     
      DO m = 1, SIZE(j_m)  
      
         !------------------------------------------------------------          
         ! >>> FOR EVERY ELEMENT
         ! ---------------------         
         !
         ! Local (of the element) index c_ of node-pair is set to zero 
         ! Loop on every possible different couple of nodes ij, i.e., 
         ! every possible node-pair
         
         c_ = 0
         DO i_ = 1, SIZE( j_m(m)%vec )       
            DO j_ = i_+1, SIZE( j_m(m)%vec )  

               i = j_m(m) % vec(i_);  j = j_m(m) % vec(j_)

               c_ = c_ + 1 
 
               !------------------------------------------------------------          
               ! >>> >>> FOR EVERY NEW NODE-PAIR
               ! -------------------------------         
               !
               ! If the connectivity ( element-nodepair ) is not 
               ! initialized we have a new node-pair
               
               IF ( c_m(m)%vec(c_) == 0) THEN  

                  ! The total number of node pair is increased
                  Nc = Nc + 1           
                  ! and the connectivity matrix for the element 
                  ! considered is initialized
                  c_m(m)%vec(c_) = Nc
 
                  ! To initialize c_m for the surrounding elements 
                  ! we have to loop on every element of the bubble
                  DO m_B_ = 1, SIZE( m_j(i)%vec ) 
                  
                  !------------------------------------------------------------          
                  ! >>> >>> >>> FOR EVERY ELEMENT OF THE BUBBLE
                  ! -------------------------------------------         
                  ! Same loop on all the node-pairs of each element m_B 
                  ! of the bubble, looking for the node-pair just inserted
                     m_B = m_j(i) % vec(m_B_)
                  
                     c_B_ = 0
                     DO i_B_ = 1, SIZE( j_m(m_B) % vec ) 
                        DO j_B_ = i_B_+1, SIZE( j_m(m_B) % vec )
                           
                           i_B = j_m(m_B)%vec(i_B_)
                           j_B = j_m(m_B)%vec(j_B_)
                           
                           c_B_ = c_B_ + 1
                           
                           ! If the connectivity matrix in m_B
                           ! is not initialized, finds the 
                           ! node-pair and sets it.  Other
                           ! node-pairs (not this one) may be
                           ! already initialized.
                           IF ( c_m(m_B)%vec(c_B_) == 0 ) THEN
                           
                              IF      ( (i == i_B) .AND. (j == j_B) ) THEN
                                 c_m(m_B)%vec(c_B_) = Nc
                              ELSE IF ( (i == j_B) .AND. (j == i_B) ) THEN
                                 c_m(m_B)%vec(c_B_) = Nc
                              ENDIF
                              
                           ENDIF
                           
                        ENDDO
                     ENDDO
                  
                  ENDDO
                  !------------------------------------------------------------          

               ENDIF
               !------------------------------------------------------------          

            ENDDO

         ENDDO
         ! -----------------------------------------------------------         

      ENDDO
      !============================================================ 


   END FUNCTION  c_m_connectivity
   !============================================================ 



   !============================================================ 
   FUNCTION  j_c_connectivity (j_m, c_m)  RESULT (j_c)
   !============================================================ 


      !------------------------------------------------------------ 
      !
      ! Starting from the element to node connectivity j_m and
      ! the element to node-pair connectivity c_m, retrieves the
      ! node-pair to node connectivity j_c         
      ! 
      ! REMARKS:
      ! -------- 
      ! 1) j_m and c_m are DIV vectors, j_c is a matrix of 
      ! dimension 4 (i,j,i*,j*) times the number of node-pairs
      ! 
      ! 2)The algorithm works on an element by element basis,
      ! so each node-pair is set more than once (in 2D, twice if
      ! at least one of the two nodes is not on the boundary, once
      ! otherwise) 
      !
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      IMPLICIT NONE
      !------------------------------------------------------------ 
      TYPE(D_I_V), DIMENSION(:),   INTENT(IN)  ::  j_m, c_m  

      INTEGER,     DIMENSION(:,:), POINTER     ::  j_c
      !------------------------------------------------------------ 
      INTEGER  ::  m       ! General element 
      INTEGER  ::  i, j    ! Points in global element coordinates
      INTEGER  ::  i_, j_  ! Points in global element coordinates
      INTEGER  ::  c       ! Node-pair in local element coordinates
      INTEGER  ::  c_      ! Node-pair in local element coordinates
      !------------------------------------------------------------ 


      ! The total number of node-pair is computed 
      ! by means of the DIV size function and
      ! the matrix j_c is allocated 

      ALLOCATE( j_c(4, size_DIV(c_m,3)) )


      !------------------------------------------------------------ 
      ! Loops on every element
      ! ----------------------
      DO m = 1, SIZE( j_m ) 

         c_ = 0
  
         !------------------------------------------------------------ 
         ! Find every possible couple ij
         ! in the element
         ! -----------------------------
         DO i_ = 1, SIZE( j_m(m)%vec )
            DO j_ = i_+1, SIZE( j_m(m)%vec ) 
               
               i = j_m(m) % vec(i_);  j = j_m(m) % vec(j_)
               
               c_ = c_ + 1  ! local index of the node-pair

               IF ( c_m(m)%vec(c_) .NE. 0) THEN

                  c = c_m(m)%vec(c_)

                  j_c(1,c) =  i;  j_c(2,c) =  j    ! <<<<<<

                  ! May 29 2006. MF. Added to avoid irrealistic
                  ! allocation in  cd_cb_connectivity  from div
                  ! grid % j_c_d_fem at line 268 of grid_utils.f90
                  j_c(3,c) =  0;  j_c(4,c) =  0

               ELSE

                  WRITE(*,*) 'Node pair',c_,' of element',m
                  WRITE(*,*) 'not initialized.'
                  WRITE(*,*) 'FUNCTION j_c_connectivity'
                  WRITE(*,*) 'in MODULE node_pair_generation. STOP'
                  STOP

               ENDIF

            ENDDO
         ENDDO
         !------------------------------------------------------------ 

      ENDDO
      !------------------------------------------------------------ 


   END FUNCTION j_c_connectivity
   !============================================================ 



   !============================================================ 
   SUBROUTINE  extend_node_pair (j_c, cs_c, rr)
   !============================================================ 


      !------------------------------------------------------------ 
      !
      ! Computes the extended node-pair structure         
      ! 
      !
      !                o---o---o
      !               / \ / \ / \
      !           i* o---o===o---o j*           
      !               \ /i\ /j\ /          
      !                o---o---o
      !
      ! The algorithm find i* [j*] by looking for the node n 
      ! to the bubble of i [j] for which the normalized scalar product
      ! 
      !     (X_i - X_j).(X_n - X_i)    [(X_j - X_i).(X_n - X_i)]
      !
      ! is maximum
      !
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER,      DIMENSION(:,:), INTENT(INOUT)  ::  j_c
      INTEGER,      DIMENSION(:,:), INTENT(OUT)    ::  cs_c
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  rr
      !------------------------------------------------------------ 
      TYPE(D_I_V), DIMENSION(:), ALLOCATABLE  :: j_c_DIV, c_j

      INTEGER  ::  c         ! General element 
      INTEGER  ::  c_B,c_B_  ! Node-pair belonging to the bubble 
      INTEGER  ::  i,j       ! Points in local element coordinates
      INTEGER  ::  is, js    ! Extended nodes 
      INTEGER  ::  orient

      REAL(KIND=8)  :: max_s_prod, s_prod
      REAL(KIND=8), DIMENSION(SIZE(rr,1)) :: DXij, DXji, DXiis, DXjjs
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Computes the node to node-pair connectivity 
      ! matrix c_j using DIV algorithms
      ! -------------------------------------------      

      j_c(3:4,:) = 0
      ALLOCATE ( j_c_DIV(SIZE(j_c,2)) )
      j_c_DIV = convert_matrix_to_DIV( j_c )
      
      ALLOCATE ( c_j(size_DIV(j_c_DIV,3)) )
      c_j = invert_DIV( j_c_DIV )
      DEALLOCATE( j_c_DIV )
      !------------------------------------------------------------ 
      


      !------------------------------------------------------------ 
      ! Loop oon every node-pair c
      ! to find i* and j*
      ! --------------------------
      
      DO c = 1, SIZE(j_c,2) 
      
         i = j_c(1,c);  j = j_c(2,c)

         !------------------------------------------------------------ 
         ! Find node i*         
         ! ------------
         !
         !      o   o
         !       \ / 
         ! i* X---o===o           c_B node-pair 
         !       /i\  j           bubble of i
         !      o   o
         !

         max_s_prod = - HUGE(max_s_prod)
         DO c_B_ = 1, SIZE( c_j(i)%vec )
         
            c_B = c_j(i)%vec(c_B_)             

            ! Selects the candidate node is
            IF      (( j_c(1,c_B) .NE. i ) .AND. &
                     ( j_c(2,c_B) .NE. j )) THEN
               
               is = j_c(1,c_B); orient = 1 

            ELSE IF (( j_c(2,c_B) .NE. i ) .AND. &
                     ( j_c(1,c_B) .NE. j )) THEN

               is = j_c(2,c_B); orient = -1 

            ELSE

               is = j   ! The node-pair is ij 

            ENDIF

            DXji  =           ( rr(:,i)  -  rr(:,j) )       &
                  / SQRT( SUM(( rr(:,i)  -  rr(:,j) )**2) )
            DXiis =           ( rr(:,is) -  rr(:,i) )       &  
                  / SQRT( SUM(( rr(:,is) -  rr(:,i) )**2) )

            s_prod = SUM( DXji * DXiis ) 

            IF (s_prod > max_s_prod) THEN
               j_c(3,c) = is
               max_s_prod = s_prod
               cs_c(1,c) = orient*c_B 
            ENDIF

         ENDDO
         !------------------------------------------------------------ 


         !------------------------------------------------------------ 
         ! Find node j*         
         ! ------------
         !
         !          o   o
         !           \ / 
         !        o===o---X j*    c_B node-pair 
         !        i  /j\          bubble of j
         !          o   o
         !

         max_s_prod = - HUGE(max_s_prod)
 
         DO c_B_ = 1, SIZE( c_j(j)%vec )
         
            c_B = c_j(j)%vec(c_B_)             
         
            ! Selects the candidate node js
            IF      (( j_c(1,c_B) .NE. j ) .AND. &
                     ( j_c(2,c_B) .NE. i )) THEN

               js = j_c(1,c_B); orient = -1 

            ELSE IF (( j_c(2,c_B) .NE. j ) .AND. &
                     ( j_c(1,c_B) .NE. i )) THEN

               js = j_c(2,c_B); orient = 1 

            ELSE

               js = j  ! The node-pair is ij 

            ENDIF

            DXij  =           ( rr(:,j)  -  rr(:,i) )       &
                  / SQRT( SUM(( rr(:,j)  -  rr(:,i) )**2) )
            DXjjs =           ( rr(:,js) -  rr(:,j) )       &
                  / SQRT( SUM(( rr(:,js) -  rr(:,j) )**2) )

            s_prod = SUM( DXij * DXjjs ) 

            IF (s_prod > max_s_prod) THEN
               j_c(4,c) = js
               max_s_prod = s_prod
               cs_c(2,c) = orient*c_B 
            ENDIF

         ENDDO
         !------------------------------------------------------------ 


      ENDDO
      !------------------------------------------------------------ 

      DEALLOCATE( c_j )
      

   END SUBROUTINE  extend_node_pair
   !============================================================ 



   !============================================================ 
   SUBROUTINE  set_extend_bou (j_c_d, cs_c_d, jd_jb, rr)
   !============================================================ 


      !------------------------------------------------------------ 
      ! 
      !  Nullifies the extended structure for node-pair
      !  coming into the boundary and reorient the node-pair 
      !  outward. Finally, sets j* == i for these node-pairs.
      !
      !           /Boundary                 /
      !  o---o---o/                o---o---o/+++o
      !  j*  j  i|/          ==>   i*  i   j/   i
      !          |/                         /  
      !        i*o/                         /     
      !           /                         /      
      !
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER, DIMENSION(:,:), INTENT(INOUT)  ::  j_c_d
      INTEGER, DIMENSION(:,:), INTENT(INOUT)  ::  cs_c_d
      INTEGER, DIMENSION(:),   INTENT(IN)     ::  jd_jb
      REAL(KIND=8), DIMENSION(:,:), INTENT(IN)     ::  rr
      !------------------------------------------------------------ 
      INTEGER  ::  c               ! General node-pair
      INTEGER  ::  i,j, is, js     ! Points of the node-pair
      LOGICAL,      DIMENSION(SIZE(rr,2))  ::  boundary 
      !REAL(KIND=8), DIMENSION(SIZE(rr,1))  ::  DXij, DXiis, &
      !                                         DXji, DXjjs
      !REAL(KIND=8)  ::  s_prod
      !------------------------------------------------------------ 

       

      boundary = .FALSE.;  boundary(jd_jb) = .TRUE.     

      DO c = 1, SIZE(j_c_d,2) 
      
         i  = j_c_d(1,c);  j  = j_c_d(2,c) 
         is = j_c_d(3,c);  js = j_c_d(4,c) 

         IF ( boundary(j) .AND. (.NOT. boundary(i)) ) THEN
            
            !------------------------------------------------------------ 
            ! I case:
            ! -------
            !
            !          /                      /
            !       j*o/                      /
            !         |/                      /
            !         |/                      /
            ! o---o---o/             o---o---o/+++o 
            ! i*  i   j/             i*  i   j/   i
            !          /                      /
            
! TEST: exstension only if angle is less that 90 degree
!            DXij  =           ( rr(:,j)  -  rr(:,i) )       &
!                  / SQRT( SUM(( rr(:,j)  -  rr(:,i) )**2) )
!            DXjjs =           ( rr(:,js) -  rr(:,j) )       &
!                  / SQRT( SUM(( rr(:,js) -  rr(:,j) )**2) )
!
!            s_prod = SUM( DXij * DXjjs ) 
!
!            IF (s_prod < 0.d0) 
                j_c_d(4,c) = i
               cs_c_d(2,c) = -c
!            ENDIF
            !------------------------------------------------------------ 

         ELSE IF ( boundary(i) .AND. (.NOT. boundary(j)) ) THEN
            
            !------------------------------------------------------------ 
            ! II case:
            ! --------
            !
            !          /                      /
            !       i*o/                      /
            !         |/                      /
            !         |/                      /
            ! o---o---o/             o---o---o/+++o 
            ! j*  j   i/             i*  i   j/   i
            !          /                      /
            

! TEST: exstension only if angle is less that 90 degree
!            DXji  =           ( rr(:,i)  -  rr(:,j) )       &
!                  / SQRT( SUM(( rr(:,i)  -  rr(:,j) )**2) )
!            DXiis =           ( rr(:,is) -  rr(:,i) )       &  
!                  / SQRT( SUM(( rr(:,is) -  rr(:,i) )**2) )
!
!            s_prod = SUM( DXji * DXiis ) 
!
!            IF (s_prod < 0.d0) THEN 
               j_c_d(1,c)  = j;   j_c_d(2,c) = i
               j_c_d(3,c)  = js;  j_c_d(4,c) = j
               cs_c_d(1,c) = - cs_c_d(2,c) 
               cs_c_d(2,c) = - c 
!            ENDIF
            !------------------------------------------------------------ 

         ENDIF

      ENDDO
     
      ! Node-pair orientation has been changed.
      ! Reorder the signs in the cs_c matrix
      DO c = 1, SIZE(j_c_d,2)

         i  = j_c_d(1,c);  j  = j_c_d(2,c) 
       
         IF ( j_c_d(2,ABS(cs_c_d(1,c))) .EQ. i) THEN
            cs_c_d(1,c) =   ABS(cs_c_d(1,c))
         ELSE 
            cs_c_d(1,c) = - ABS(cs_c_d(1,c))
         ENDIF 
      
         IF ( j_c_d(1,ABS(cs_c_d(2,c))) .EQ. j) THEN
            cs_c_d(2,c) =   ABS(cs_c_d(2,c))
         ELSE 
            cs_c_d(2,c) = - ABS(cs_c_d(2,c))
         ENDIF 
      
      
      ENDDO


   END SUBROUTINE set_extend_bou 
   !============================================================ 
   


   !============================================================ 
   FUNCTION cd_cb_connectivity ( j_c_d, jd_jb,  j_c_b ) RESULT (cd_cb)
   !============================================================ 




      !------------------------------------------------------------ 
      IMPLICIT NONE
      !------------------------------------------------------------ 
      INTEGER,     DIMENSION(:,:),       INTENT(IN)     ::  j_c_d
      INTEGER,     DIMENSION(:),         INTENT(IN)     ::  jd_jb
      INTEGER,     DIMENSION(:,:),       INTENT(INOUT)  ::  j_c_b

      INTEGER,     DIMENSION(SIZE(j_c_b,2))             ::  cd_cb
      !------------------------------------------------------------ 
      TYPE(D_I_V), DIMENSION(:), ALLOCATABLE  ::  j_c_d_DIV, c_j_d 
                                                       
      INTEGER  ::  ib, jb,      &  ! Nodes in local element indexes
                   id,  jd,     &  ! Nodes in local element indexes
                   cb, cd, cd_, &  ! Node-pair in local element coordinates
                   nd, tmp         ! Total number of node-pair in the element
                       
 
      !------------------------------------------------------------ 
      ! Create the node to node-pair connectivity 
      ! (domain) inverting j_c_d (in DIV format)
      ! -----------------------------------------
      ALLOCATE (j_c_d_DIV(SIZE(j_c_d,2)))
      j_c_d_DIV = convert_matrix_to_DIV(j_c_d)

      ALLOCATE (c_j_d(size_DIV(j_c_d_DIV,3)))
      c_j_d = invert_DIV(j_c_d_DIV)
      !------------------------------------------------------------ 

      cd_cb = 0

      DO cb = 1, SIZE(j_c_b,2)
                                   ! Bubble of one of the nodes of cb
         nd =  jd_jb(j_c_b(1,cb))  ! (node 1 is always present) 

         ib = jd_jb(j_c_b(1,cb));  jb = jd_jb(j_c_b(2,cb)) 

         DO cd_ = 1, SIZE(c_j_d(nd)%vec);  cd = c_j_d(nd)%vec(cd_)
                       
            id = j_c_d(1,cd);  jd = j_c_d(2,cd) 
            
            IF      ( (ib == id) .AND. (jb == jd) ) THEN 
            
               cd_cb(cb) = cd 
            
            ELSE IF ( (jb == id) .AND. (ib == jd) ) THEN
            
               cd_cb(cb) = cd 
               
               ! swap i and j
               tmp = j_c_b(2,cb); j_c_b(2,cb) = j_c_b(1,cb); j_c_b(1,cb) = tmp       
            
               ! swap is and js
               tmp = j_c_b(4,cb); j_c_b(4,cb) = j_c_b(3,cb); j_c_b(3,cb) = tmp       
           
            ENDIF
   
         ENDDO

      ENDDO

      IF ( ANY( cd_cb == 0) ) THEN
        WRITE(*,*) ' I was unable to find the vector cd_cb. '
        WRITE(*,*) ' domain node_pair ', cd,'   STOP'
        STOP
      ENDIF
      

   END FUNCTION  cd_cb_connectivity
   !============================================================ 

   
   !============================================================ 
   FUNCTION bound_connectivity(c_m_b, bound_m, Nc_b) RESULT(bound_c)
   !============================================================ 
   

      IMPLICIT NONE

      TYPE (D_I_V), DIMENSION(:),          INTENT(IN)  ::  c_m_b
      INTEGER,      DIMENSION(:),          INTENT(IN)  ::  bound_m
      INTEGER,                             INTENT(IN)  ::  Nc_b
      !INTEGER,      DIMENSION(SIZE(c_m_b))             ::  bound_c
      INTEGER,      DIMENSION(Nc_b)                    ::  bound_c
   
      INTEGER  ::  m

      DO m = 1, SIZE( c_m_b )
   
         bound_c(c_m_b(m)%vec) = bound_m(m)
   
      ENDDO

      
   END FUNCTION bound_connectivity
   !============================================================ 
   
   
   
!============================================================ 
!************************************************************  
!
!  FINITE VOLUME METHOD
!
!************************************************************  
!============================================================ 



   !============================================================ 
   FUNCTION fv_fe_np_connectivity( j_c, j_m, c_m, ele_type,  &
                                   Nc_fv ) RESULT( fv_fe )
   !============================================================ 
   

      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER,      DIMENSION(:,:), INTENT(IN)   ::  j_c
      TYPE (D_I_V), DIMENSION(:),   INTENT(IN)   ::  j_m, c_m
      INTEGER,      DIMENSION(:),   INTENT(IN)   ::  ele_type
      INTEGER,                      INTENT(OUT)  ::  Nc_fv

      INTEGER,      DIMENSION(SIZE(j_c,2))       ::  fv_fe
      !------------------------------------------------------------ 
      INTEGER, DIMENSION(:,:), POINTER  ::  edges
      LOGICAL, DIMENSION(SIZE(j_c,2))   ::  keep
      INTEGER  ::  m, i, j, i_e, j_e, c, c_, e
      !------------------------------------------------------------ 



      keep = .FALSE.

 
      !------------------------------------------------------------ 
      ! Loop on all the element to select the
      ! node-pair to be kept ( keep = true ) 
      ! ------------------------------------- 

      DO m = 1, SIZE(j_m)  ! <<< Loop on the elements

         ! List of the element's edges
         edges => ele_edges(ele_type(m)) 

         DO c_ = 1, SIZE(c_m(m)%vec) ! <<< Loop on the node-pairs
                                     !     of element m
            c = c_m(m)%vec(c_)
            
            i = j_c(1,c);  j = j_c(2,c)
            
            DO e = 1, SIZE(edges,1)  ! <<< Controls whether the 
                                     !     node-pair is an edge
            
               i_e = j_m(m)%vec(edges(e,1))
               j_e = j_m(m)%vec(edges(e,2))
            
               ! If the node-pair is also an edge, keep it
               ! -----------------------------------------   
               IF ( ((i == i_e ) .AND. (j == j_e)) .OR. &
                    ((j == i_e ) .AND. (i == j_e))      ) THEN
               
                  keep(c) = .TRUE.  ! <<<<<<
               
               ENDIF
            
            ENDDO
         
         ENDDO

         DEALLOCATE( edges )    !(S)

      ENDDO 
      !------------------------------------------------------------ 
      

      !------------------------------------------------------------
      ! Set the connectivity matrix fv_fe between
      ! FE node-pair indices and FV indices
      ! (0 if not kept) and compute total number
      ! of FV node-pairs (Nc_fv)
      ! -----------------------------------------

      fv_fe = 0;  Nc_fv = 0
      DO c = 1, SIZE(keep)

         IF ( keep(c) ) THEN

            Nc_fv = Nc_fv + 1
            fv_fe(c) = Nc_fv

         ENDIF

      ENDDO
      !------------------------------------------------------------


      
   END FUNCTION fv_fe_np_connectivity
   !============================================================ 



   !============================================================ 
   FUNCTION fv_node_pair( j_c, fv_fe, Nc_fv )  RESULT(j_c_fv)
                          
   !============================================================ 
   

      !------------------------------------------------------------ 
      IMPLICIT NONE

      INTEGER, DIMENSION(:,:), INTENT(IN)  ::  j_c
      INTEGER, DIMENSION(:),   INTENT(IN)  ::  fv_fe
      INTEGER,                 INTENT(IN)  ::  Nc_fv

      INTEGER, DIMENSION(4,Nc_fv)          ::  j_c_fv
      !------------------------------------------------------------ 
      INTEGER  ::  c, Nc
      !------------------------------------------------------------ 



      Nc = 0
      DO c = 1, SIZE(fv_fe)

         IF ( fv_fe(c) .NE. 0 ) THEN

            Nc = Nc + 1
            j_c_fv(:,Nc) = j_c(:,c)

         ENDIF

      ENDDO 
      !------------------------------------------------------------ 

      
   END FUNCTION fv_node_pair
   !============================================================ 



   !============================================================ 
   FUNCTION fv_c_m_connectivity( c_m, fv_fe ) RESULT(c_m_fv)
   !============================================================ 
   

      !------------------------------------------------------------ 
      IMPLICIT NONE

      TYPE (D_I_V), DIMENSION(:),   INTENT(IN)   ::  c_m
      INTEGER,      DIMENSION(:),   INTENT(IN)   ::  fv_fe

      TYPE (D_I_V), DIMENSION(SIZE(c_m))         ::  c_m_fv
      !------------------------------------------------------------ 
      INTEGER  ::  m, c, c_, Nc_ele
      !------------------------------------------------------------ 


      !------------------------------------------------------------ 
      ! Element to FV node-pair connectivity i
      ! matrix initialization
      DO m = 1, SIZE(c_m)

         ! Total number of node-pairs in this element
         ! ------------------------------------------
         Nc_ele = 0
         DO c_ = 1, SIZE(c_m(m)%vec)
         
            c = c_m(m)%vec(c_)
            IF ( fv_fe(c) .NE. 0 ) Nc_ele = Nc_ele + 1
         
         ENDDO

         ! Element/node-pair connectivity is stored
         ! ----------------------------------------
         ALLOCATE( c_m_fv(m)%vec(Nc_ele) )
         Nc_ele = 0
         DO c_ = 1, SIZE(c_m(m)%vec)
         
            c = c_m(m)%vec(c_)
            IF ( fv_fe(c) .NE. 0 ) THEN
               Nc_ele = Nc_ele + 1
               c_m_fv(m)%vec(Nc_ele) = fv_fe(c)
            ENDIF
         
         ENDDO

      ENDDO 
      !------------------------------------------------------------ 
      
   END FUNCTION fv_c_m_connectivity
   !============================================================ 


END MODULE  np_topology_gen
