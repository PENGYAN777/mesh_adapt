!============================================================ 
!
!      Module: np_topology_gen
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
!    Modified: Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
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
!   Extended node-pair structure (2D sample):         
!
!             itm(7)     jtp(6)
!              |          |
!              |          |
!     o--------O==========O---------o
!   i*(3)     i(1)       j(2)     j*(4)
!              |          | 
!              |          |
!             itp(5)     jtm(8)
!
!
!
!============================================================ 

  MODULE  np_topology_gen

   USE dynamic_vector
   USE element_topology
   USE mesh_structure
   USE node_pair_structure
   USE nodes

   TYPE(D_I_V), DIMENSION(:), POINTER, PUBLIC :: m_c

   CONTAINS


   SUBROUTINE node_pair_str_gen(grid_name)
   !---------------------------------------------------------
   IMPLICIT NONE

   CHARACTER(LEN=64), INTENT(IN) :: grid_name

   TYPE(D_I_V), DIMENSION(:), POINTER :: c_mfv!, m_c

   INTEGER, DIMENSION(:), ALLOCATABLE  ::  fv_fe_np
   INTEGER :: space_dim
   !---------------------------------------------------------

   space_dim = SIZE(rr,1)

   ! DOMAIN NODE-PAIRS

   ! Connectivity matrix initialization
   ALLOCATE (c_m_d(Ne_d))
   c_m_d = c_m_connectivity(j_m_d)

   Nc_d = size_DIV(c_m_d,3)
   WRITE(*,*) '   I have found ', Nc_d, ' node-pairs in the domain'


   ! GENERATION OF THE NODE-PAIRS STRUCTURE

   ! Node-pairs
   ALLOCATE (j_c_d(4*space_dim, Nc_d))
   j_c_d = j_c_connectivity(j_m_d, c_m_d)

   ! Extended node-pairs
   ALLOCATE (cs_c_d(2, Nc_d))
   CALL extend_node_pair(j_c_d, cs_c_d, rr, 'FE')


   IF (Nc_d /= 0) THEN
     ! Nullify extension on boundaries  
     ! and reorder the node-pair outward
     CALL set_extend_bou(j_c_d, cs_c_d, jd_jb, rr, 'FE', grid_name)
   ENDIF


   ! SELECTION OF THE NODE-PAIR TO BE KEPT IN THE FV APPROACH

   ALLOCATE (fv_fe_np(Nc_d))
   fv_fe_np = fv_fe_np_connectivity(j_c_d, j_m_d, c_m_d, ele_type_d, Nc_fv)

   ALLOCATE (j_c_fv(4*space_dim, Nc_fv))
   j_c_fv  = fv_node_pair(j_c_d, fv_fe_np, Nc_fv)


   ! C_M_fv
   ALLOCATE (c_mfv(Ne_d))
   c_mfv = fv_c_m_connectivity(c_m_d, fv_fe_np)
   
   ! M_C_fv 
   ALLOCATE (m_c(size_DIV(c_mfv, 3)))
   m_c = invert_DIV(c_mfv)
   
   DEALLOCATE (c_mfv)


   ! Extended node-pairs      
   ALLOCATE (cs_c_fv(2, Nc_d))
   CALL extend_node_pair (j_c_fv, cs_c_fv, rr, 'FV', m_c, j_m_d)

   IF (Nc_d /= 0) THEN
     ! Nullify extension on boundaries and reorder the node-pair outward
     CALL set_extend_bou (j_c_fv, cs_c_fv, jd_jb, rr, 'FV', grid_name)
   ENDIF
   
   ALLOCATE (c_m_fv(Ne_d))
   c_m_fv  = fv_c_m_connectivity(c_m_d, fv_fe_np)
   DEALLOCATE (fv_fe_np)

   !DEALLOCATE (m_c)




   ! BOUNDARY NODE-PAIRS

   ! Connectivity matrix initialization
   IF (Ne_b > 0) THEN

     ALLOCATE (c_m_b(Ne_b))

     c_m_b = c_m_connectivity(j_m_b)

     Nc_b = size_DIV(c_m_b, 3)
     WRITE(*,*) '   I have found ', Nc_b, ' node-pairs in the boundary'

     ! Generation of the node-pairs structure
     ALLOCATE (j_c_b(4*space_dim, Nc_b))

     ! Node-pairs
     j_c_b = j_c_connectivity(j_m_b, c_m_b)
   
     ! No extension is needed on boundary node-pair, otherwise 
     ! CALL extend_node_pair (j_c_b, rr(:,jd_jb))
     j_c_b(3:4*space_dim,:) = 0

     ! Boundary node-pair to domain node-pair connectivity.
     ! It may also change the order of the nodes in boundary 
     ! node-pairs to fit the orientation of the corresponding 
     ! domain node-pairs
     ALLOCATE (cd_cb(Nc_b))
     cd_cb = cd_cb_connectivity(j_c_d, jd_jb,  j_c_b)

     ! Node-pair to boundary patch connectivity
     ALLOCATE (bound_c(Nc_b))
     bound_c = bound_connectivity(c_m_b, bound_m, Nc_b)

   ELSE   

      Nc_b = 0

   ENDIF

   END SUBROUTINE node_pair_str_gen





   SUBROUTINE  ja_j_connectivity(j_c, ja_j)
   !--------------------------------------------------------------------
   USE dynamic_vector

   IMPLICIT NONE

   INTEGER,     DIMENSION(:,:), INTENT(IN)  :: j_c
   TYPE(D_I_V), DIMENSION(:),   INTENT(OUT) :: ja_j
   
   INTEGER, DIMENSION(SIZE(j_c,1),SIZE(j_c,2)) :: j_c_copy
         
   TYPE(D_I_V), DIMENSION(:), POINTER :: c_j      
   TYPE(D_I_V), DIMENSION(:), POINTER :: j_c_DIV

   INTEGER :: i, j, j_, c, c_
   !--------------------------------------------------------------------
   
   ! Computes the node to node-pair connectivity 
   ! matrix c_j using DIV algorithms
   j_c_copy = j_c
   j_c_copy(3:SIZE(j_c,1),:) = 0

   ALLOCATE (j_c_DIV(SIZE(j_c_copy,2)))
   j_c_DIV = convert_matrix_to_DIV(j_c_copy)
   
   ALLOCATE (c_j(size_DIV(j_c_DIV,3)))
   c_j = invert_DIV(j_c_DIV)
   DEALLOCATE (j_c_DIV)

   DO j_ = 1, SIZE(c_j)
   
     ALLOCATE (ja_j(j_)%vec(SIZE(c_j(j_)%vec)))
     
     DO c_ = 1, SIZE(c_j(j_)%vec)
     
       c = c_j(j_)%vec(c_)

       i = j_c(1,c)
       j = j_c(2,c)
       
       IF (i == j_)  ja_j(j_)%vec(c_) = j
       IF (j == j_)  ja_j(j_)%vec(c_) = i
       
     ENDDO     

   ENDDO

   END SUBROUTINE  ja_j_connectivity





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

   DEALLOCATE ( m_j )


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
      USE nodes, ONLY: k_d

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

      ALLOCATE( j_c(4*k_d, size_DIV(c_m,3)) )


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




   SUBROUTINE  extend_node_pair (j_c, cs_c, rr, feOfv, m_c, j_m)

   ! Computes the extended node-pair structure
   !
   !                 itm jtp
   !                o---o---o
   !               / \ / \ / \
   !           is o---o===o---o js           
   !               \ /i\ /j\ /          
   !                o---o---o
   !                 itp jtm
   !
   !  ibp, ibm binormal to plane associated to node i 
   !  jbp, jbm binormal to plane associated to node j
   !
   ! The algorithm find i* [j*] by looking for the node n 
   ! to the bubble of i [j] for which the normalized scalar product
   ! 
   !     (X_i - X_j).(X_n - X_i)    [(X_j - X_i).(X_n - X_i)]
   !
   ! is maximum
   !
   ! January 2007. Modified by Marco Fossati to include
   !               extended node-pairs and extended nodes
   !               in the directions normal and binormal 
   !               to the edge connecting two nodes i-j.
   !
   !       Legend:  - n = versor j-i [i-j]
   !                - t = perpendicular versor to n
   !                - b = binormal versor to n and t
   ! 
   !                - itp = index of node in direction t+
   !                - itm = index of node in direction t-
   !                - ibp = index of node in direction b+
   !                - ibm = index of node in direction b-
   !
   !                Similarly for indices jtp, jtm, jbp, jbm
   !-------------------------------------------------------------------------------
   USE metric_coefficients,   ONLY: Drr, Drr_fv
!   USE faces_gen,             ONLY: f_c_connectivity

   IMPLICIT NONE

   INTEGER,      DIMENSION(:,:), INTENT(INOUT) :: j_c
   INTEGER,      DIMENSION(:,:), INTENT(OUT)   :: cs_c
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: rr
   CHARACTER(LEN=2),             INTENT(IN)    :: feOfv
   TYPE(D_I_V),  DIMENSION(:),   INTENT(IN), OPTIONAL :: m_c, j_m

   TYPE(D_I_V), DIMENSION(:), ALLOCATABLE  :: j_c_DIV, c_j

   INTEGER, DIMENSION(:), ALLOCATABLE :: k, tgNodes

   INTEGER :: c                        ! General element 
   INTEGER :: c_B,c_B_                 ! Node-pair belonging to the bubble 
   INTEGER :: i,j                      ! Points in local element coordinates
   INTEGER :: is, itp, itm, ibp, ibm   ! extended nodes associated to node i
   INTEGER :: js, jtp, jtm, jbp, jbm   ! extended nodes associated to node j
   INTEGER :: orient
   INTEGER :: Nc, space_dim, Nm_c, m_, m, p, Nj_c, tgn

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: N_sProd, TP_sProd, TM_sProd,  &
                                              BP_sProd, BM_sProd

   REAL(KIND=8), DIMENSION(SIZE(rr,1)) :: DXij, DXji, DXiis, DXjjs
   REAL(KIND=8), DIMENSION(SIZE(rr,1), SIZE(rr,1)) :: cosines

   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Drr_   
   !-------------------------------------------------------------------------------

   space_dim = SIZE(rr,1)
   Nc = SIZE(j_c,2)

   ! Computes the node to node-pair connectivity 
   ! matrix c_j using DIV algorithms
   j_c(3:4*space_dim,:) = 0

   ALLOCATE (j_c_DIV(SIZE(j_c,2)))
   j_c_DIV = convert_matrix_to_DIV(j_c)
   
   ALLOCATE (c_j(size_DIV(j_c_DIV,3)))
   c_j = invert_DIV(j_c_DIV)
   DEALLOCATE (j_c_DIV)

   ALLOCATE (Drr_(4*space_dim-1, Nc))


   ! Loop on every node-pair c to find extended nodes
   DO c = 1, SIZE(j_c,2) 

      i = j_c(1,c)
      j = j_c(2,c)

      ! Length of edge i-j
      Drr_(1,c) = SQRT(SUM((rr(:,j) - rr(:,i))**2))

      ! Versor i->j
      DXij = (rr(:,j) - rr(:,i)) / Drr_(1,c)
      ! Versor j->i
      DXji = (rr(:,i) - rr(:,j)) / Drr_(1,c)


      IF (feOfv == 'FV') THEN

        ! Selecting candidates among which select for extended
        ! nodes in directions tangent to i--j normal
        IF (space_dim == 2) THEN

          Nm_c = SIZE(m_c(c)%vec)
          Nj_c = 0
	 
          ! Counting the number of nodes of the intersection
          ! between the bubble of i and j
          DO m_ = 1, Nm_c
            m = m_c(c)%vec(m_)
            Nj_c = Nj_c + SIZE(j_m(m)%vec) - 2
          ENDDO

          ALLOCATE (tgNodes(Nj_c))

          ! Selecting the nodes belonging to the intersection
          ! of the bubbles of the two nodes i and j excluded
          ! nodes i and j.
	  tgn = 1
          DO m_ = 1, Nm_c
	  
            m = m_c(c)%vec(m_)
	    
            DO p = 1, SIZE(j_m(m)%vec)
              IF (j_m(m)%vec(p) /= i  .AND.  &
        	  j_m(m)%vec(p) /= j) THEN
		
		tgNodes(tgn) = j_m(m)%vec(p)
		tgn = tgn + 1
		
              ENDIF
            ENDDO
	    
          ENDDO

        ELSEIF (space_dim == 3) THEN


!          CALL f_c_connectivity
          STOP

        ENDIF

      ENDIF





      !---------------------------------------------------------------------------- 
      ! Find extended nodes for I
      !      o   o
      !       \ / 
      ! i* X---o===o           c_B = bubble of node-pairs
      !       /i\  j                 issuing from i
      !      o   o

      ALLOCATE (k(SIZE(c_j(i) % vec)), N_sProd(SIZE(c_j(i) % vec)))

      IF (space_dim >= 2) THEN
        ALLOCATE (TP_sProd(SIZE(c_j(i) % vec)), TM_sProd(SIZE(c_j(i) % vec)))
        TP_sProd = 0.d0
        TM_sProd = 0.d0
      ENDIF

      IF (space_dim == 3)  ALLOCATE (BP_sProd(SIZE(c_j(i) % vec)),  &
                                     BM_sProd(SIZE(c_j(i) % vec)))


      DO c_B_ = 1, SIZE(c_j(i) % vec)

         c_B = c_j(i) % vec(c_B_)

         ! Selects the k node i*
         IF (( j_c(1,c_B) .NE. i ) .AND. &
             ( j_c(2,c_B) .NE. j )) THEN

            orient = 1
            k(c_B_) = j_c(1,c_B)

         ELSE IF (( j_c(2,c_B) .NE. i ) .AND. &
                  ( j_c(1,c_B) .NE. j )) THEN

            orient = -1
            k(c_B_) = j_c(2,c_B)

         ELSE

            ! The candidate is the node-pair ij 
            k(c_B_) = j

         ENDIF

         ! Computing versors in the two directions, normal 
         ! and binormal to versor DXji
         CALL np_versors(DXji, cosines)


         ! Versor i--k
         DXiis =           ( rr(:,k(c_B_)) -  rr(:,i) )       &  
               / SQRT( SUM(( rr(:,k(c_B_)) -  rr(:,i) )**2) )

         ! Scalar product with normal
         N_sProd(c_B_) = SUM(DXiis*cosines(1,:))

         ! Scalar product with tangent (+ & -)
	 IF (feOfv == 'FV') THEN
           IF (space_dim >= 2) THEN
	     IF(ANY(k(c_B_) == tgNodes)) THEN
               TP_sProd(c_B_) = SUM(DXiis*  cosines(2,:))
               TM_sProd(c_B_) = SUM(DXiis*(-cosines(2,:)))
             ENDIF
	   ENDIF

           ! Scalar product with binormal (+ & -)
           IF (space_dim == 3) THEN
             BP_sProd(c_B_) = SUM(DXiis*  cosines(3,:))
             BM_sProd(c_B_) = SUM(DXiis*(-cosines(3,:)))
           ENDIF
         ENDIF
	 
       ENDDO  
       

       ! Extended node in the NORMAL DIRECTION (In*) ---------------------------
       IF (MAXVAL(N_sProd) > 0.d0) THEN
         is = k(MAXLOC(N_sProd,1))
         j_c(3,c) = is
         cs_c(1,c) = orient * c_j(i) % vec(MAXLOC(N_sProd,1))
         Drr_(2,c) = ABS(SUM((rr(:,is) -  rr(:,i))*cosines(1,:)))
       ELSE
         j_c(3,c) = i
         cs_c(1,c) = orient * c
         Drr_(2,c) = Drr_(1,c)
       ENDIF


      
       ! Extended node in the TANGENT DIRECTION (It*) --------------------------
       IF (feOfv == 'FV'  .AND.  space_dim >= 2) THEN

         IF (MAXVAL(TP_sProd) > 0.d0) THEN
           itp = k(MAXLOC(TP_sProd,1))
           j_c(5,c) = itp
           Drr_(4,c) = ABS(SUM((rr(:,itp) -  rr(:,i))*cosines(2,:)))
         ELSE
           j_c(5,c) = i
           Drr_(4,c) = Drr_(1,c)
         ENDIF

         IF (MAXVAL(TM_sProd) > 0.d0) THEN
           itm = k(MAXLOC(TM_sProd,1))
           j_c(7,c) = itm
           Drr_(6,c) = ABS(SUM((rr(:,itm) -  rr(:,i))*(-cosines(2,:))))
         ELSE
           j_c(7,c) = i
           Drr_(6,c) = Drr_(1,c)
         ENDIF

       ENDIF



       ! Extended node in the BINORMAL DIRECTION (Ib*) -------------------------
       IF (feOfv == 'FV'  .AND.  space_dim == 3) THEN

         IF (MAXVAL(BP_sProd) > 0.d0) THEN
           ibp = k(MAXLOC(BP_sProd,1))
           j_c(9,c) = ibp
           Drr_(8,c) = ABS(SUM((rr(:,ibp) -  rr(:,i))*cosines(3,:)))
         ELSE
           j_c(9,c) = i
           Drr_(8,c) = Drr_(1,c)
         ENDIF

         IF (MAXVAL(BM_sProd) > 0.d0) THEN
           ibm = k(MAXLOC(BM_sProd,1))
           j_c(11,c) = ibm
           Drr_(10,c) = ABS(SUM((rr(:,ibm) -  rr(:,i))*(-cosines(3,:))))
         ELSE
           j_c(11,c) = i
           Drr_(10,c) = Drr_(1,c)
         ENDIF

       ENDIF

       DEALLOCATE (k, N_sProd)
       IF (space_dim >= 2)  DEALLOCATE (TP_sProd, TM_sProd)
       IF (space_dim == 3)  DEALLOCATE (BP_sProd, BM_sProd)





      !------------------------------------------------------------ 
      ! Find extended nodes for J         
      !
      !          o   o
      !           \ / 
      !        o===o---X j*    c_B = bubble of node-pairs
      !        i  /j\                issuing from node j
      !          o   o

      ALLOCATE (k(SIZE(c_j(j) % vec)), N_sProd(SIZE(c_j(j) % vec)))

      IF (space_dim >= 2) THEN
        ALLOCATE (TP_sProd(SIZE(c_j(j) % vec)), TM_sProd(SIZE(c_j(j) % vec)))
        TP_sProd = 0.d0
        TM_sProd = 0.d0
      ENDIF


      IF (space_dim == 3)  ALLOCATE (BP_sProd(SIZE(c_j(j) % vec)),  &
                                     BM_sProd(SIZE(c_j(j) % vec)))


      DO c_B_ = 1, SIZE(c_j(j) % vec)
      
         c_B = c_j(j) % vec(c_B_)             
      
         ! Selects the k node
         IF      (( j_c(1,c_B) .NE. j ) .AND. &
                  ( j_c(2,c_B) .NE. i )) THEN

            orient = -1 
            k(c_B_) = j_c(1,c_B)

         ELSE IF (( j_c(2,c_B) .NE. j ) .AND. &
                  ( j_c(1,c_B) .NE. i )) THEN

            orient = 1 
            k(c_B_) = j_c(2,c_B)

         ELSE

            ! The k is the node-pair ij 
            k(c_B_) = j  

         ENDIF


         ! Computing versors in the two directions, normal 
         ! and binormal to versor DXij
         CALL  np_versors(DXij, cosines)

         ! Versor j--k
         DXjjs =           ( rr(:,k(c_B_)) -  rr(:,j) )       &
               / SQRT( SUM(( rr(:,k(c_B_)) -  rr(:,j) )**2) )

         ! Scalar product with normal
         N_sProd(c_B_) = SUM(DXjjs*cosines(1,:))

         ! Scalar product with tangent (+ & -)
	 IF (feOfv == 'FV') THEN
           IF (space_dim >= 2) THEN
	     IF(ANY(k(c_B_) == tgNodes)) THEN
               TP_sProd(c_B_) = SUM(DXjjs*  cosines(2,:))
               TM_sProd(c_B_) = SUM(DXjjs*(-cosines(2,:)))
	     ENDIF  
           ENDIF

           ! Scalar product with binormal (+ & -)
           IF (feOfv == 'FV'  .AND.  space_dim == 3) THEN
             BP_sProd(c_B_) = SUM(DXjjs*  cosines(3,:))
             BM_sProd(c_B_) = SUM(DXjjs*(-cosines(3,:)))
           ENDIF
	 ENDIF

      ENDDO



       ! Extended node in the NORMAL DIRECTION (Jn*) ---------------------------
      IF (MAXVAL(N_sProd) > 0.d0) THEN
        js = k(MAXLOC(N_sProd,1))
        j_c(4,c) = js
        cs_c(2,c) = orient * c_j(j) % vec(MAXLOC(N_sProd,1))
        Drr_(3,c) = ABS(SUM((rr(:,js) -  rr(:,j))*cosines(1,:)))
      ELSE
        j_c(4,c) = j
        cs_c(2,c) = orient * c
        Drr_(3,c) = Drr_(1,c)
      ENDIF


       ! Extended node in the TANGENT DIRECTION (Jt*) --------------------------
      IF (feOfv == 'FV'  .AND.  space_dim >= 2) THEN

        IF (MAXVAL(TP_sProd) > 0.d0) THEN
          jtp = k(MAXLOC(TP_sProd,1))
          j_c(6,c) = jtp
          Drr_(5,c) = ABS(SUM((rr(:,jtp) -  rr(:,j))*cosines(2,:)))
        ELSE
          j_c(6,c) = j
          Drr_(5,c) = Drr_(1,c)
        ENDIF

        IF (MAXVAL(TM_sProd) > 0.d0) THEN
          jtm = k(MAXLOC(TM_sProd,1))
          j_c(8,c) = jtm
          Drr_(7,c) = ABS(SUM((rr(:,jtm) -  rr(:,j))*(-cosines(2,:))))
        ELSE
          j_c(8,c) = j
          Drr_(7,c) = Drr_(1,c)
        ENDIF

      ENDIF



       ! Extended node in the BINORMAL DIRECTION (Jb*) -------------------------
      IF (feOfv == 'FV'  .AND.  space_dim == 3) THEN

        IF (MAXVAL(BP_sProd) > 0.d0) THEN
          jbp = k(MAXLOC(BP_sProd,1))
          j_c(10,c) = jbp
          Drr_(9,c) = ABS(SUM((rr(:,jbp) -  rr(:,j))*cosines(3,:)))
        ELSE
          j_c(10,c) = j
          Drr_(9,c) = Drr_(1,c)
        ENDIF

        IF (MAXVAL(BM_sProd) > 0.d0) THEN
          jbm = k(MAXLOC(BM_sProd,1))
          j_c(12,c) = jbm
          Drr_(11,c) = ABS(SUM((rr(:,jbm) -  rr(:,j))*(-cosines(3,:))))
        ELSE
          j_c(12,c) = j
          Drr_(11,c) = Drr_(1,c)
        ENDIF

      ENDIF

      DEALLOCATE (k, N_sProd)
      IF (space_dim >= 2)  DEALLOCATE (TP_sProd, TM_sProd)
      IF (space_dim == 3)  DEALLOCATE (BP_sProd, BM_sProd)

      IF (ALLOCATED(tgNodes))  DEALLOCATE(tgNodes)

   ENDDO



   

   IF (feOfv == 'FE'  .AND.  ANY(Drr_(1:3,:) == 0)) THEN
     PRINT*, ''
     PRINT*, 'ERROR. EXTEND NODE_PAIR:'
     PRINT*, 'Zero projected length found.'
     PRINT*, ''
     STOP
   ENDIF


   IF (feOfv == 'FV'  .AND.  ANY(Drr_ == 0)) THEN
     PRINT*, ''
     PRINT*, 'ERROR. EXTEND NODE_PAIR:'
     PRINT*, 'Zero projected length found.'
     PRINT*, ''
     STOP
   ENDIF


   ! Allocate arrays for projected edges lenghts
   IF (feOfv == 'FE') THEN

     ALLOCATE (Drr(4*space_dim-1, Nc))
     Drr = Drr_

   ELSEIF (feOfv == 'FV') THEN

     ALLOCATE (Drr_fv(4*space_dim-1, Nc))
     Drr_fv = Drr_

   ENDIF

   DEALLOCATE (c_j, Drr_)
   
   END SUBROUTINE  extend_node_pair





   SUBROUTINE  set_extend_bou (j_c_d, cs_c_d, jd_jb, rr, feOfv, grid_name)

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
   !------------------------------------------------------------
   USE metric_coefficients,   ONLY: Drr, Drr_fv 
   USE nodes,                 ONLY: bound_p

   IMPLICIT NONE

   INTEGER,      DIMENSION(:,:), INTENT(INOUT) :: j_c_d
   INTEGER,      DIMENSION(:,:), INTENT(INOUT) :: cs_c_d
   INTEGER,      DIMENSION(:),   INTENT(IN)    :: jd_jb
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: rr
   CHARACTER(LEN=2),             INTENT(IN)    :: feOfv
   CHARACTER(LEN=64),            INTENT(IN)    :: grid_name

   LOGICAL, DIMENSION(SIZE(rr,2)) :: boundary, angle_node
   INTEGER, DIMENSION(SIZE(rr,2)) :: boundaryId
  
   INTEGER, DIMENSION(:), ALLOCATABLE :: angle_nodeID
   
   INTEGER :: c, i, j, is, js, sDim, itp, itm, jtp, jtm,  &
              ibp, ibm, jbp, jbm, p, Nangles
   
   REAL(KIND=8) :: Drr_4, Drr_5, Drr_6,  Drr_7,  &
                   Drr_8, Drr_9, Drr_10, Drr_11
                   
   LOGICAL :: FEX = .FALSE.                  
   !------------------------------------------------------------ 

   sDim = SIZE(rr,1)

   angle_node = .FALSE.
   boundary   = .FALSE.
   
   boundary(jd_jb) = .TRUE.     

   boundaryId = 0
   DO p = 1, SIZE(bound_p)
     boundaryId(jd_jb(p)) = bound_p(p)
   ENDDO


   INQUIRE (file='corners.'//TRIM(grid_name), exist=FEX)

   IF (FEX) THEN
   
     OPEN(UNIT=13, FILE='corners.'//TRIM(grid_name))
       READ(13,*) Nangles
       ALLOCATE (angle_nodeID(Nangles))
       DO p = 1, Nangles
         READ(13,*) angle_nodeID(p)
       ENDDO
     CLOSE(13)
   
     angle_node(angle_nodeID) = .TRUE.
   
     DEALLOCATE (angle_nodeID)
   
   ELSE

     DO p = 1, SIZE(jd_jb)
       
       ! Last (and implicitly first) node in the list of boundary nodes
       IF (p == SIZE(jd_jb)) THEN
         IF (bound_p(p) /= bound_p(1))  angle_node(jd_jb(p)) = .TRUE.
         EXIT
       ENDIF
        
       IF (bound_p(p) /= bound_p(p+1))  angle_node(jd_jb(p)) = .TRUE.
   
     ENDDO
     
     OPEN(UNIT=13, FILE='corners.'//TRIM(grid_name), STATUS='unknown')
       WRITE(13,*) COUNT(angle_node)
       DO p = 1, SIZE(angle_node) 
         IF (angle_node(p))  WRITE(13,*) p
       ENDDO
     CLOSE(13)     
     
   ENDIF


   DO c = 1, SIZE(j_c_d,2) 
   
      i  = j_c_d(1,c);  j  = j_c_d(2,c) 
      is = j_c_d(3,c);  js = j_c_d(4,c)

      ! CASE I ---------------------------------------------------------------------------
      IF (boundary(j)  .AND.  (.NOT. boundary(i))) THEN

         !          /                      /
         !       j*o/                      /
         !         |/                      /
         !         |/                      /
         ! o---o---o/             o---o---o/+++o 
         ! i*  i   j/             i*  i   j/   i
         !          /                      /
         
         j_c_d(4,c) = i
         cs_c_d(2,c) = -c

         IF (feOfv == 'FE')  Drr(3,c) = Drr(1,c)
         IF (feOfv == 'FV')  Drr_fv(3,c) = Drr_fv(1,c)


      ! CASE II --------------------------------------------------------------------------
      ELSE IF (boundary(i)  .AND.  (.NOT. boundary(j))) THEN

         !          /                      /
         !       i*o/                      /
         !         |/                      /
         !         |/                      /
         ! o---o---o/             o---o---o/+++o 
         ! j*  j   i/             i*  i   j/   i
         !          /                      /
         
         j_c_d(1,c)  = j;   j_c_d(2,c) = i
         j_c_d(3,c)  = js;  j_c_d(4,c) = j
         cs_c_d(1,c) = - cs_c_d(2,c) 
         cs_c_d(2,c) = - c 

         IF (feOfv == 'FE')  THEN
	   Drr(2,c) = Drr(3,c)
	   Drr(3,c) = Drr(1,c)
         ELSEIF (feOfv == 'FV') THEN
	   Drr_fv(2,c) = Drr_fv(3,c)
	   Drr_fv(3,c) = Drr_fv(1,c)
	 ENDIF  


         IF (sDim >= 2) THEN
	 
	   itp = j_c_d(5,c);   jtp = j_c_d(6,c)
 	   itm = j_c_d(7,c);   jtm = j_c_d(8,c)

           j_c_d(5,c) = jtp;   j_c_d(6,c) = itp
           j_c_d(7,c) = jtm;   j_c_d(8,c) = itm
	   
           IF (feOfv == 'FE')  THEN
             
	     Drr_4 = Drr(4,c);   Drr_5 = Drr(5,c)
             Drr_6 = Drr(6,c);   Drr_7 = Drr(7,c)
	     
	     Drr(4,c) = Drr_5;   Drr(5,c) = Drr_4
	     Drr(6,c) = Drr_7;   Drr(7,c) = Drr_6
	     	     
	   ELSEIF (feOfv == 'FV') THEN
             
	     Drr_4 = Drr_fv(4,c);   Drr_5 = Drr_fv(5,c)
             Drr_6 = Drr_fv(6,c);   Drr_7 = Drr_fv(7,c)

	     Drr_fv(4,c) = Drr_5;   Drr_fv(5,c) = Drr_4
	     Drr_fv(6,c) = Drr_7;   Drr_fv(7,c) = Drr_6

	   ENDIF
	   
	 ENDIF


         IF (sDim == 3) THEN
	 
	   ibp = j_c_d(9,c);    jbp = j_c_d(10,c)
 	   ibm = j_c_d(11,c);   jbm = j_c_d(12,c)

           j_c_d(9,c)  = jbp;   j_c_d(10,c) = ibp
           j_c_d(11,c) = jbm;   j_c_d(12,c) = ibm

           IF (feOfv == 'FE')  THEN
             
	     Drr_8 = Drr(8,c);     Drr_9 = Drr(9,c)
             Drr_10 = Drr(10,c);   Drr_11 = Drr(11,c)
	     
	     Drr(8,c) = Drr_9;     Drr(9,c) = Drr_8
	     Drr(10,c) = Drr_11;   Drr(11,c) = Drr_10
	     	    
	   ELSEIF (feOfv == 'FV') THEN
             
	     Drr_8 = Drr_fv(8,c);     Drr_9 = Drr_fv(9,c)
             Drr_10 = Drr_fv(10,c);   Drr_11 = Drr_fv(11,c)
	     
	     Drr_fv(8,c) = Drr_9;     Drr_fv(9,c) = Drr_8
	     Drr_fv(10,c) = Drr_11;   Drr_fv(11,c) = Drr_10

	   ENDIF

	 ENDIF

      ENDIF


      ! CASE III -------------------------------------------------------------------------
        IF (angle_node(i)  .AND.  boundaryId(j) == boundaryId(js)) THEN

	 j_c_d(1,c)  = j;   j_c_d(2,c) = i
         j_c_d(3,c)  = js;  j_c_d(4,c) = j
         cs_c_d(1,c) = - cs_c_d(2,c) 
         cs_c_d(2,c) = - c

         IF (feOfv == 'FE')  THEN
	   Drr(2,c) = Drr(3,c)
	   Drr(3,c) = Drr(1,c)
         ELSEIF (feOfv == 'FV') THEN
	   Drr_fv(2,c) = Drr_fv(3,c)
	   Drr_fv(3,c) = Drr_fv(1,c)
	 ENDIF  


         IF (sDim >= 2) THEN
	 
	   itp = j_c_d(5,c);   jtp = j_c_d(6,c)
 	   itm = j_c_d(7,c);   jtm = j_c_d(8,c)

           j_c_d(5,c) = jtp;   j_c_d(6,c) = itp
           j_c_d(7,c) = jtm;   j_c_d(8,c) = itm
	   
           IF (feOfv == 'FE')  THEN
             
	     Drr_4 = Drr(4,c);   Drr_5 = Drr(5,c)
             Drr_6 = Drr(6,c);   Drr_7 = Drr(7,c)
	     
	     Drr(4,c) = Drr_5;   Drr(5,c) = Drr_4
	     Drr(6,c) = Drr_7;   Drr(7,c) = Drr_6
	     	     
	   ELSEIF (feOfv == 'FV') THEN
             
	     Drr_4 = Drr_fv(4,c);   Drr_5 = Drr_fv(5,c)
             Drr_6 = Drr_fv(6,c);   Drr_7 = Drr_fv(7,c)

	     Drr_fv(4,c) = Drr_5;   Drr_fv(5,c) = Drr_4
	     Drr_fv(6,c) = Drr_7;   Drr_fv(7,c) = Drr_6

	   ENDIF
	   
	 ENDIF


         IF (sDim == 3) THEN
	 
	   ibp = j_c_d(9,c);    jbp = j_c_d(10,c)
 	   ibm = j_c_d(11,c);   jbm = j_c_d(12,c)

           j_c_d(9,c)  = jbp;   j_c_d(10,c) = ibp
           j_c_d(11,c) = jbm;   j_c_d(12,c) = ibm

           IF (feOfv == 'FE')  THEN
             
	     Drr_8 = Drr(8,c);     Drr_9 = Drr(9,c)
             Drr_10 = Drr(10,c);   Drr_11 = Drr(11,c)
	     
	     Drr(8,c) = Drr_9;     Drr(9,c) = Drr_8
	     Drr(10,c) = Drr_11;   Drr(11,c) = Drr_10
	     	    
	   ELSEIF (feOfv == 'FV') THEN
             
	     Drr_8 = Drr_fv(8,c);     Drr_9 = Drr_fv(9,c)
             Drr_10 = Drr_fv(10,c);   Drr_11 = Drr_fv(11,c)
	     
	     Drr_fv(8,c) = Drr_9;     Drr_fv(9,c) = Drr_8
	     Drr_fv(10,c) = Drr_11;   Drr_fv(11,c) = Drr_10

	   ENDIF

	 ENDIF

       ENDIF


     ! CASE IV --------------------------------------------------------------------------
     IF (angle_node(j)  .AND.  boundaryId(i) == boundaryId(is)) THEN

       j_c_d(4,c) = i
       cs_c_d(2,c) = -c
       
       IF (feOfv == 'FE')  Drr(3,c) = Drr(1,c)
       IF (feOfv == 'FV')  Drr_fv(3,c) = Drr_fv(1,c)

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

   END SUBROUTINE  set_extend_bou 

   



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
                                                       
      INTEGER  ::  ib, jb, &  ! Nodes in local element indexes
                   id,  jd,  &  ! Nodes in local element indexes
                   cb, cd, cd_,     &  ! Node-pair in local element coordinates
                   nd, tmp        ! Total number of node-pair in the element
                       


 
      !------------------------------------------------------------ 
      ! Create the node to node-pair connectivity 
      ! (domain) inverting j_c_d (in DIV format)
      ! -----------------------------------------
      ALLOCATE ( j_c_d_DIV(SIZE(j_c_d,2)) )
      j_c_d_DIV = convert_matrix_to_DIV ( j_c_d )
      ALLOCATE ( c_j_d(size_DIV(j_c_d_DIV,3)) )
      c_j_d = invert_DIV( j_c_d_DIV )
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
        WRITE(*,*) ' domain node_pair ', cd, id, jd, '   STOP'
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
      USE nodes, ONLY: k_d

      IMPLICIT NONE

      INTEGER, DIMENSION(:,:), INTENT(IN)  ::  j_c
      INTEGER, DIMENSION(:),   INTENT(IN)  ::  fv_fe
      INTEGER,                 INTENT(IN)  ::  Nc_fv

      INTEGER, DIMENSION(4*k_d,Nc_fv)          ::  j_c_fv
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





   SUBROUTINE  np_versors(n, versors)

   ! Versors is the cosine direction matrix for the
   ! direction aligned with normal vector (n) and
   ! two orthogonal directions tangent to the node pair
   ! interface (t, b)
   ! 1st row contains cosine directors of versors in
   ! the direction normal to the interface.
   ! The other two rows the versors in the 
   ! two tangential directions
   !---------------------------------------------------------------------------
   USE  lin_algebra,   ONLY: vec_prod

   IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), INTENT(IN) :: n

   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: versors

   REAL(KIND=8), DIMENSION(2, SIZE(n)) :: VV

   REAL(KIND=8), DIMENSION(SIZE(n)) :: t, b, temp

   REAL(KIND=8) :: mod_n, mod_t, mod_b

   INTEGER :: space_dim, z_pos
   !---------------------------------------------------------------------------

   space_dim = SIZE(n)

   mod_n = SQRT(SUM(n*n))

   versors(1,:) = n / mod_n

   ! Computing versors in directions 
   ! normal to n
   IF (space_dim == 2) THEN

     versors(2,:) = (/ -versors(1,2), versors(1,1) /)

   ELSEIF (space_dim == 3) THEN

     VV(1,:) = versors(1,:)

     IF (COUNT(versors(1,:) == 0.d0) == 2) THEN

       VV(2,:) = versors(1,:)

       temp = versors(1,:)*versors(1,:)
 
       z_pos = MINLOC(temp,1)

       VV(2,z_pos) = VV(2,z_pos) + 1.d0

     ELSE

       VV(2,:) = (/ versors(1,1),  versors(1,2), versors(1,3) + 1.d0 /)

     ENDIF  

     t = vec_prod(VV)

     VV(1,:) = t
     VV(2,:) = versors(1,:)

     b = vec_prod(VV)

     mod_t = SQRT(SUM(t*t))
     mod_b = SQRT(SUM(b*b))

     versors(2,:) = t / mod_t
     versors(3,:) = b / mod_b

   ENDIF

   END SUBROUTINE  np_versors



END MODULE  np_topology_gen
