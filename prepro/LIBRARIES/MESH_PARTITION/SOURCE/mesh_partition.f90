! ==============================================================================
!
!      Module: mesh_partition
!
! Description: Module devoted to mesh partitioning.
!              METIS library is used to compute partition
!              of the unstructured hybrid grid converted
!              to a graph.
!
!      Author: Marco Fossati
!              Department of Aerospace Engineering
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
! ==============================================================================

   MODULE  mesh_partition

   CONTAINS


   SUBROUTINE  partMesh ( grid_name, nparts, supp )
   !-------------------------------------------------------------------------
   USE dynamic_vector
   USE mesh_structure
   USE node_pair_structure 
   USE element_topology
   USE np_topology_gen 
   USE nodes

   !-------------------------------------------------------------------------
   IMPLICIT NONE

   CHARACTER(LEN=64), INTENT(IN) :: grid_name
   INTEGER,           INTENT(IN) :: nparts, supp

   TYPE(D_I_V), DIMENSION(:,:), ALLOCATABLE :: jFringes_FE, jFringes_FV
   TYPE(D_I_V), DIMENSION(:,:), ALLOCATABLE :: cFringes_FE, cFringes_FV
   TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE :: cE_FE, cE_FV
   TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE :: mFringes
   TYPE(D_I_V), DIMENSION(:),   ALLOCATABLE :: ja_j, boundPart_FE, boundPart_FV

   INTEGER, DIMENSION(SIZE(j_c_d,1), SIZE(j_c_d,2))  :: j_c_FE_copy
   INTEGER, DIMENSION(SIZE(j_c_fv,1),SIZE(j_c_fv,2)) :: j_c_FV_copy   
         
   TYPE(D_I_V), DIMENSION(:), POINTER :: c_j      
   TYPE(D_I_V), DIMENSION(:), POINTER :: j_c_DIV

 
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xm, ym
   REAL(KIND=8) :: x1, x2, y1, y2, xC, yC
 
   INTEGER, PARAMETER :: idf = 11
   INTEGER :: i, j, k, c, nF, p, q, NmP, c_
   INTEGER :: n, m, wgtflag, numflag, edgecut, nFringes_ij
   INTEGER :: nCe, m_, ck
   INTEGER :: boundIdx, loc_nBounds
 
   INTEGER, DIMENSION(:), ALLOCATABLE :: part, partM, partJ, partC_FE, partC_FV 
   INTEGER, DIMENSION(:), ALLOCATABLE :: N_part, fC
   INTEGER, DIMENSION(:), ALLOCATABLE :: xadj, adjncy 
   INTEGER, DIMENSION(:), ALLOCATABLE :: vwgt, adjwgt
   INTEGER, DIMENSION(:), ALLOCATABLE :: cSize, nJ, pJ
 
   INTEGER, DIMENSION(0:4) :: options
 
   LOGICAL, DIMENSION(:), ALLOCATABLE :: flagNodes, flagBound, cE_flag
   LOGICAL, PARAMETER :: debug = .FALSE.
   !-------------------------------------------------------------------------

   ALLOCATE (N_part(nparts))
   ALLOCATE (partM(Ne_d), partJ(Np_d))
   ALLOCATE (partC_FE(Nc_d), partC_FV(Nc_fv))

   IF (supp == 0) THEN
   ! Nodes of the graph are NODES of the mesh

       ! Allocate graph dimensions
       n = Np_d 
       ALLOCATE (xadj(n+1), vwgt(n))
       ALLOCATE (part(n))
     
       m = 2*Nc_fv
       ALLOCATE (adjncy(m), adjwgt(m))


       ! Build NODAL mesh graph connectivity in CSR format
       ALLOCATE (ja_j(Np_d))
       CALL ja_j_connectivity(j_c_fv, ja_j)

       xadj(1) = 1
       adjncy = 0
       j = 1
 
       DO i = 1, SIZE(xadj)-1
          
         adjncy(j:j+SIZE(ja_j(i)%vec)-1) = ja_j(i)%vec

         j = j + SIZE(ja_j(i)%vec)
         xadj(i+1) = j
         
       ENDDO

   ELSE
   ! Nodes of the graph are ELEMENTS of the mesh
 
       ! Allocate graph dimensions
       n = Ne_d
       ALLOCATE (xadj(n+1), vwgt(n))
       ALLOCATE (part(n))
       
       m = 0
       DO i = 1, Ne_d
         m = m + COUNT(ma_m_d(i)%vec /= 0)
       ENDDO
       ALLOCATE (adjncy(m), adjwgt(m))


       ! Build DUAL MESH graph connectivity in CSR format
       xadj(1) = 1
       adjncy = 0
       j = 0
 
       DO i = 1, SIZE(xadj)-1
       
         DO k = 1, SIZE(ma_m_d(i)%vec)
            IF (ma_m_d(i)%vec(k) /= 0) THEN
              j = j + 1
              adjncy(j) = ma_m_d(i)%vec(k)       
            ENDIF
         ENDDO
         
         xadj(i+1) = xadj(i) + COUNT(ma_m_d(i)%vec /= 0)
         
       ENDDO

   ENDIF
 
 
   ! Partition graph
   PRINT*, '   ...calling METIS library'
    
   numflag = 1
   wgtflag = 0
   vwgt = 0
   adjwgt = 0 
   options = 0
 
   IF (nparts <= 8) THEN

     CALL METIS_PartGraphRecursive(n, xadj, adjncy, vwgt, adjwgt,     &
                                   wgtflag, numflag, nparts, options, &
                                   edgecut, part)
   ELSE

     CALL METIS_PartGraphKway(n, xadj, adjncy, vwgt, adjwgt,     &
                              wgtflag, numflag, nparts, options, &
                              edgecut, part)
 
   ENDIF
 
 
 
   IF (supp == 0) THEN

     ! Partition ELEMENTS according to NODES partition
     partJ = part
     ALLOCATE (mFringes(nParts))
     
     ! Loop over elements to assign partition. An element
     ! is considered to belong to a partition if at least
     ! one of its nodes belongs to the partition
     DO j = 1, nparts
       
       NmP = 0
       
       DO i = 1, Ne_d
         IF (ANY(partJ(j_m_d(i)%vec) == j))  NmP = NmP + 1
       ENDDO

       ALLOCATE (mFringes(j)%vec(NmP))
       
       k = 1
       
       DO i = 1, Ne_d
         IF (ANY(partJ(j_m_d(i)%vec) == j)) THEN
	   mFringes(j)%vec(k) = i
	   k = k + 1
	 ENDIF  
       ENDDO
       
     ENDDO
 
   ELSE
 
     ! Partition NODES according to ELEMENTS partition
     partM = part
     
     ALLOCATE (m_j_d(size_DIV(j_m_d,3)))
     m_j_d = invert_DIV (j_m_d)

 
     ! Loop over nodes to assign partition
     DO i = 1, Np_d

       DO j = 1, nparts
         N_part(j) = COUNT(partM(m_j_d(i)%vec) == j)
       ENDDO
       
       partJ(i) = MAXLOC(N_part,1)
       
     ENDDO
 
   ENDIF
 
 
 
 
 
 
 
   ! Partition Node-Pairs according to nodes partition
 
   ! Finite element notation
   partC_FE = -1
   DO c = 1, Nc_d
     IF (partJ(j_c_d(1,c)) == partJ(j_c_d(2,c)))  partC_FE(c) = partJ(j_c_d(1,c))
   ENDDO
 
   ! Fringe node-pairs in FE notation
   ALLOCATE (cFringes_FE(nparts,nparts))

   DO i = 1, nparts
     ALLOCATE (cFringes_FE(i,i)%vec(0))
   ENDDO

   DO i = 1, nparts-1
     DO j = i+1, nparts
       
       nF = 1
       nFringes_ij = COUNT((partJ(j_c_d(1,:)) == i  .AND.  partJ(j_c_d(2,:)) == j) .OR. &
                           (partJ(j_c_d(1,:)) == j  .AND.  partJ(j_c_d(2,:)) == i))
       
       ALLOCATE (cFringes_FE(i,j)%vec(nFringes_ij))
       ALLOCATE (cFringes_FE(j,i)%vec(nFringes_ij))
       
       DO c = 1, Nc_d
         IF ((partJ(j_c_d(1,c)) == i  .AND.  partJ(j_c_d(2,c)) == j) .OR. &
             (partJ(j_c_d(1,c)) == j  .AND.  partJ(j_c_d(2,c)) == i))  THEN
           
           cFringes_FE(i,j)%vec(nF) = c
           nF = nF + 1
           
         ENDIF  
       ENDDO

       cFringes_FE(j,i)%vec = cFringes_FE(i,j)%vec

     ENDDO
   ENDDO
 
 
  ! Extended node-pairs in FE notation
  ALLOCATE (cE_FE(nparts))
  ALLOCATE (cE_flag(Nc_d))
  
  ! Computes the node to node-pair connectivity 
  ! matrix c_j using DIV algorithms
  j_c_FE_copy = j_c_d
  j_c_FE_copy(3:SIZE(j_c_d,1),:) = 0

  ALLOCATE (j_c_DIV(SIZE(j_c_FE_copy,2)))
  j_c_DIV = convert_matrix_to_DIV(j_c_FE_copy)
  
  ALLOCATE (c_j(size_DIV(j_c_DIV,3)))
  c_j = invert_DIV(j_c_DIV)
  DEALLOCATE (j_c_DIV)


  DO i = 1, nparts
  
    nCe = 0
    cE_flag = .FALSE.
  
    DO j = 1, nparts
    
      IF (i == j  .OR.  SIZE(cFringes_FE(i,j) % vec) == 0) CYCLE    
      
      DO c_ = 1, SIZE(cFringes_FE(i,j)%vec)
      
        c = cFringes_FE(i,j) % vec(c_)
      
        IF (partJ(j_c_d(1,c)) /= i) THEN	  
	  k = j_c_d(1,c)	  	  
	ELSEIF (partJ(j_c_d(2,c)) /= i) THEN
	  k = j_c_d(2,c)
	ENDIF
	
	DO m_ = 1, SIZE(c_j(k)%vec)

          m = c_j(k)%vec(m_)

          IF (.NOT.cE_flag(m)  .AND.  .NOT. ANY(cFringes_FE(i,j)%vec == m)) THEN
	    nCe = nCe + 1
	    cE_flag(m) = .TRUE.
	  ENDIF
	    
	ENDDO
      
      ENDDO    
    
    ENDDO
  
    
    ALLOCATE (cE_FE(i) % vec(nCe))
  
  
    ck = 1  
    cE_flag = .FALSE.
    
    DO j = 1, nparts
    
      IF (i == j  .OR.  SIZE(cFringes_FE(i,j) % vec) == 0) CYCLE    
      
      DO c_ = 1, SIZE(cFringes_FE(i,j)%vec)
      
        c = cFringes_FE(i,j) % vec(c_)
      
        IF (partJ(j_c_d(1,c)) /= i) THEN	  
	  k = j_c_d(1,c)	  	  
	ELSEIF (partJ(j_c_d(2,c)) /= i) THEN
	  k = j_c_d(2,c)
	ENDIF
	
	DO m_ = 1, SIZE(c_j(k)%vec)

          m = c_j(k)%vec(m_)

          IF (.NOT.cE_flag(m)  .AND.  .NOT. ANY(cFringes_FE(i,j)%vec == m)) THEN
	    cE_FE(i) % vec(ck) = m
	    ck = ck + 1
	    cE_flag(m) = .TRUE.
	  ENDIF
	    
	ENDDO
      
      ENDDO    
    
    ENDDO
        
  ENDDO  
  
 
  DEALLOCATE (cE_flag)
 
 
 
 
 
  ! Finite volumes notation
  partC_FV = -1
  DO c = 1, Nc_fv
    IF (partJ(j_c_fv(1,c)) == partJ(j_c_fv(2,c)))  partC_FV(c) = partJ(j_c_fv(1,c)) 
  ENDDO

  ! Fringe node-pairs in FV notation
  ALLOCATE (cFringes_FV(nparts,nparts))

  DO i = 1, nparts
    ALLOCATE (cFringes_FV(i,i)%vec(0))
  ENDDO


   DO i = 1, nparts-1
     DO j = i+1, nparts
       
       nF = 1
       nFringes_ij = COUNT((partJ(j_c_fv(1,:)) == i  .AND.  partJ(j_c_fv(2,:)) == j) .OR. &
                           (partJ(j_c_fv(1,:)) == j  .AND.  partJ(j_c_fv(2,:)) == i))
       
       ALLOCATE (cFringes_FV(i,j)%vec(nFringes_ij))
       ALLOCATE (cFringes_FV(j,i)%vec(nFringes_ij))
       
       DO c = 1, Nc_fv
         IF ((partJ(j_c_fv(1,c)) == i  .AND.  partJ(j_c_fv(2,c)) == j) .OR. &
             (partJ(j_c_fv(1,c)) == j  .AND.  partJ(j_c_fv(2,c)) == i))  THEN
           cFringes_FV(i,j)%vec(nF) = c
           nF = nF + 1
           
         ENDIF  
       ENDDO
       
       cFringes_FV(j,i)%vec = cFringes_FV(i,j)%vec
       
     ENDDO
   ENDDO

  DEALLOCATE (c_j)



  ALLOCATE (cE_FV(nparts))
  ALLOCATE (cE_flag(Nc_fv))
  
  ! Computes the node to node-pair connectivity 
  ! matrix c_j using DIV algorithms
  j_c_FV_copy = j_c_fv
  j_c_FV_copy(3:SIZE(j_c_fv,1),:) = 0

  ALLOCATE (j_c_DIV(SIZE(j_c_FV_copy,2)))
  j_c_DIV = convert_matrix_to_DIV(j_c_FV_copy)
  
  ALLOCATE (c_j(size_DIV(j_c_DIV,3)))
  c_j = invert_DIV(j_c_DIV)
  DEALLOCATE (j_c_DIV)


  DO i = 1, nparts
  
    nCe = 0
    cE_flag = .FALSE.
  
    DO j = 1, nparts
    
      IF (i == j  .OR.  SIZE(cFringes_FV(i,j) % vec) == 0) CYCLE    
      
      DO c_ = 1, SIZE(cFringes_FV(i,j)%vec)
      
        c = cFringes_FV(i,j) % vec(c_)
      
        IF (partJ(j_c_fv(1,c)) /= i) THEN	  
	  k = j_c_fv(1,c)	  	  
	ELSEIF (partJ(j_c_fv(2,c)) /= i) THEN
	  k = j_c_fv(2,c)
	ENDIF
	
	DO m_ = 1, SIZE(c_j(k)%vec)

          m = c_j(k)%vec(m_)

          IF (.NOT.cE_flag(m)  .AND.  .NOT. ANY(cFringes_FV(i,j)%vec == m)) THEN
	    nCe = nCe + 1
	    cE_flag(m) = .TRUE.
	  ENDIF
	    
	ENDDO
      
      ENDDO    
    
    ENDDO
  
    
    ALLOCATE (cE_FV(i) % vec(nCe))
  
  
    ck = 1  
    cE_flag = .FALSE.
    
    DO j = 1, nparts
    
      IF (i == j  .OR.  SIZE(cFringes_FV(i,j) % vec) == 0) CYCLE    
      
      DO c_ = 1, SIZE(cFringes_FV(i,j)%vec)
      
        c = cFringes_FV(i,j) % vec(c_)
      
        IF (partJ(j_c_fv(1,c)) /= i) THEN	  
	  k = j_c_fv(1,c)	  	  
	ELSEIF (partJ(j_c_fv(2,c)) /= i) THEN
	  k = j_c_fv(2,c)
	ENDIF
	
	DO m_ = 1, SIZE(c_j(k)%vec)

          m = c_j(k)%vec(m_)

          IF (.NOT.cE_flag(m)  .AND.  .NOT. ANY(cFringes_FV(i,j)%vec == m)) THEN
	    cE_FV(i) % vec(ck) = m
	    ck = ck + 1
	    cE_flag(m) = .TRUE.
	  ENDIF
	    
	ENDDO
      
      ENDDO    
    
    ENDDO
        
  ENDDO

  DEALLOCATE (cE_flag)




   ! Fringe Nodes Finite Elements Notation   
   ALLOCATE (jFringes_FE(nparts,nparts))

   ALLOCATE (cSize(nparts), nJ(nparts), pJ(nparts)) 
   ALLOCATE (flagNodes(SIZE(partJ)))
   
   DO i = 1, nparts
   
     flagNodes = .FALSE.
     nJ = 0
   
     ! Loop over fringe pairs
     DO j = 1, nparts
     
       IF (i == j  .OR.  SIZE(cFringes_FE(i,j)%vec) == 0)  CYCLE
     
       DO c_ = 1, SIZE(cFringes_FE(i,j)%vec)
       
         c = cFringes_FE(i,j)%vec(c_)
       
         IF (partJ(j_c_d(1,c)) /= i  .AND. .NOT.flagNodes(j_c_d(1,c))) THEN
           nJ(j) = nJ(j) + 1
	   flagNodes(j_c_d(1,c)) = .TRUE.
	 ELSEIF (partJ(j_c_d(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_d(2,c))) THEN
           nJ(j) = nJ(j) + 1
	   flagNodes(j_c_d(2,c)) = .TRUE.
	 ENDIF
       
       ENDDO
     
     ENDDO
   
   
     ! Loop over extended pairs
     DO c_ = 1, SIZE(cE_FE(i)%vec)
       
       c = cE_FE(i)%vec(c_)
       
       IF (partJ(j_c_d(1,c)) /= i  .AND.  .NOT.flagNodes(j_c_d(1,c))) THEN
         nJ(partJ(j_c_d(1,c))) = nJ(partJ(j_c_d(1,c))) + 1
	 flagNodes(j_c_d(1,c)) = .TRUE.
       ELSEIF (partJ(j_c_d(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_d(2,c))) THEN
         nJ(partJ(j_c_d(2,c))) = nJ(partJ(j_c_d(2,c))) + 1
	 flagNodes(j_c_d(2,c)) = .TRUE.
       ENDIF
       
     ENDDO
      
      
     ! Allocating
     DO j = 1, nparts
       ALLOCATE (jFringes_FE(i,j)%vec(nJ(j)))
     ENDDO
     
     flagNodes = .FALSE.
     pJ = 1

     DO j = 1, nparts
     
       IF (i == j  .OR.  SIZE(cFringes_FE(i,j)%vec) == 0)  CYCLE
     
       DO c_ = 1, SIZE(cFringes_FE(i,j)%vec)
       
         c = cFringes_FE(i,j)%vec(c_)
       
         IF (partJ(j_c_d(1,c)) /= i  .AND. .NOT.flagNodes(j_c_d(1,c))) THEN
           jFringes_FE(i,j) % vec(pJ(j)) = j_c_d(1,c)
           pJ(j) = pJ(j) + 1 
	   flagNodes(j_c_d(1,c)) = .TRUE.
	 ELSEIF (partJ(j_c_d(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_d(2,c))) THEN
           jFringes_FE(i,j) % vec(pJ(j)) = j_c_d(2,c)
           pJ(j) = pJ(j) + 1 
	   flagNodes(j_c_d(2,c)) = .TRUE.
	 ENDIF
       
       ENDDO
     
     ENDDO
   
   
     ! Loop over extended pairs
     DO c_ = 1, SIZE(cE_FE(i)%vec)
       
       c = cE_FE(i)%vec(c_)
       
       IF (partJ(j_c_d(1,c)) /= i  .AND.  .NOT.flagNodes(j_c_d(1,c))) THEN
         jFringes_FE(i,partJ(j_c_d(1,c))) % vec(pJ(partJ(j_c_d(1,c)))) = j_c_d(1,c)
         pJ(partJ(j_c_d(1,c))) = pJ(partJ(j_c_d(1,c))) + 1 
	 flagNodes(j_c_d(1,c)) = .TRUE.
       ELSEIF (partJ(j_c_d(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_d(2,c))) THEN
         jFringes_FE(i,partJ(j_c_d(2,c))) % vec(pJ(partJ(j_c_d(2,c)))) = j_c_d(2,c)
         pJ(partJ(j_c_d(2,c))) = pJ(partJ(j_c_d(2,c))) + 1 
	 flagNodes(j_c_d(2,c)) = .TRUE.
       ENDIF
       
     ENDDO
   
   ENDDO





   ! Fringe Nodes Finite Volumes Notation.
   ALLOCATE (jFringes_FV(nparts,nparts))
   
   DO i = 1, nparts
   
     flagNodes = .FALSE.
     nJ = 0
   
     ! Loop over fringe pairs
     DO j = 1, nparts
     
       IF (i == j  .OR.  SIZE(cFringes_FV(i,j)%vec) == 0)  CYCLE
     
       DO c_ = 1, SIZE(cFringes_FV(i,j)%vec)
       
         c = cFringes_FV(i,j)%vec(c_)
       
         IF (partJ(j_c_fv(1,c)) /= i  .AND. .NOT.flagNodes(j_c_fv(1,c))) THEN
           nJ(j) = nJ(j) + 1
	   flagNodes(j_c_fv(1,c)) = .TRUE.
	 ELSEIF (partJ(j_c_fv(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_fv(2,c))) THEN
           nJ(j) = nJ(j) + 1
	   flagNodes(j_c_fv(2,c)) = .TRUE.
	 ENDIF
       
       ENDDO
     
     ENDDO
   
   
     ! Loop over extended pairs
     DO c_ = 1, SIZE(cE_FV(i)%vec)
       
       c = cE_FV(i)%vec(c_)
       
       IF (partJ(j_c_fv(1,c)) /= i  .AND.  .NOT.flagNodes(j_c_fv(1,c))) THEN
         nJ(partJ(j_c_fv(1,c))) = nJ(partJ(j_c_fv(1,c))) + 1
	 flagNodes(j_c_fv(1,c)) = .TRUE.
       ELSEIF (partJ(j_c_fv(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_fv(2,c))) THEN
         nJ(partJ(j_c_fv(2,c))) = nJ(partJ(j_c_fv(2,c))) + 1
	 flagNodes(j_c_fv(2,c)) = .TRUE.
       ENDIF
       
     ENDDO
      
      
     ! Allocating
     DO j = 1, nparts
       ALLOCATE (jFringes_FV(i,j)%vec(nJ(j)))
     ENDDO
     
     flagNodes = .FALSE.
     pJ = 1

     DO j = 1, nparts
     
       IF (i == j  .OR.  SIZE(cFringes_FV(i,j)%vec) == 0)  CYCLE
     
       DO c_ = 1, SIZE(cFringes_FV(i,j)%vec)
       
         c = cFringes_FV(i,j)%vec(c_)
       
         IF (partJ(j_c_fv(1,c)) /= i  .AND. .NOT.flagNodes(j_c_fv(1,c))) THEN
           jFringes_FV(i,j) % vec(pJ(j)) = j_c_fv(1,c)
           pJ(j) = pJ(j) + 1 
	   flagNodes(j_c_fv(1,c)) = .TRUE.
	 ELSEIF (partJ(j_c_fv(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_fv(2,c))) THEN
           jFringes_FV(i,j) % vec(pJ(j)) = j_c_fv(2,c)
           pJ(j) = pJ(j) + 1 
	   flagNodes(j_c_fv(2,c)) = .TRUE.
	 ENDIF
       
       ENDDO
     
     ENDDO
   
   
     ! Loop over extended pairs
     DO c_ = 1, SIZE(cE_FV(i)%vec)
       
       c = cE_FV(i)%vec(c_)
       
       IF (partJ(j_c_fv(1,c)) /= i  .AND.  .NOT.flagNodes(j_c_fv(1,c))) THEN
         jFringes_FV(i,partJ(j_c_fv(1,c))) % vec(pJ(partJ(j_c_fv(1,c)))) = j_c_fv(1,c)
         pJ(partJ(j_c_fv(1,c))) = pJ(partJ(j_c_fv(1,c))) + 1 
	 flagNodes(j_c_fv(1,c)) = .TRUE.
       ELSEIF (partJ(j_c_fv(2,c)) /= i  .AND.  .NOT.flagNodes(j_c_fv(2,c))) THEN
         jFringes_FV(i,partJ(j_c_fv(2,c))) % vec(pJ(partJ(j_c_fv(2,c)))) = j_c_fv(2,c)
         pJ(partJ(j_c_fv(2,c))) = pJ(partJ(j_c_fv(2,c))) + 1 
	 flagNodes(j_c_fv(2,c)) = .TRUE.
       ENDIF
       
     ENDDO
   
   ENDDO   
   
   
   
   
   
   




   ! Identify the boundaries belonging to each partition
   ! FE notation
   ALLOCATE (boundPart_FE(nparts))
   ALLOCATE (flagBound(MAXVAL(bound_p)))

   DO p = 1, nparts
   
     loc_nBounds = 0
     flagBound = .FALSE.
     
     DO j = 1, Np_b

       DO m = 1, nparts

         IF (p == m  .OR.  SIZE(jFringes_FE(p,m)%vec) == 0)  CYCLE
     
         IF (partJ(jd_jb(j)) == p  .OR.  ANY(jFringes_FE(p,m)%vec == jd_jb(j))) THEN
       
           boundIdx = bound_p(j)
	   
	   IF (.NOT. flagBound(boundIdx)) THEN
	     loc_nBounds = loc_nBounds + 1
	     flagBound(boundIdx) = .TRUE.
	   ENDIF  
       
         ENDIF
       
       ENDDO
       
     ENDDO     
     
     ALLOCATE (boundPart_FE(p) % vec(loc_nBounds))
     
   ENDDO


   DO p = 1, nparts

     k = 1 
     flagBound = .FALSE.
     
     DO j = 1, Np_b

       DO m = 1, nparts

         IF (p == m  .OR.  SIZE(jFringes_FE(p,m)%vec) == 0)  CYCLE
      
         IF (partJ(jd_jb(j)) == p  .OR.  ANY(jFringes_FE(p,m)%vec == jd_jb(j))) THEN

           boundIdx = bound_p(j)
	   
	   IF (.NOT. flagBound(boundIdx)) THEN
             boundPart_FE(p) % vec(k) = boundIdx
	     flagBound(boundIdx) = .TRUE.
	     k = k + 1
	   ENDIF  
       
         ENDIF
	 
       ENDDO
       
     ENDDO
     
   ENDDO








   ! Identify the boundaries belonging to each partition
   ! FV notation
   ALLOCATE (boundPart_FV(nparts))

   DO p = 1, nparts
   
     loc_nBounds = 0
     flagBound = .FALSE.
     
     DO j = 1, Np_b

       DO m = 1, nparts

         IF (p == m  .OR.  SIZE(jFringes_FV(p,m)%vec) == 0)  CYCLE
     
         IF (partJ(jd_jb(j)) == p  .OR.  ANY(jFringes_FV(p,m)%vec == jd_jb(j))) THEN
       
           boundIdx = bound_p(j)
	   
	   IF (.NOT. flagBound(boundIdx)) THEN
	     loc_nBounds = loc_nBounds + 1
	     flagBound(boundIdx) = .TRUE.
	   ENDIF  
       
         ENDIF
       
       ENDDO
       
     ENDDO     
     
     ALLOCATE (boundPart_FV(p) % vec(loc_nBounds))
     
   ENDDO


   DO p = 1, nparts

     k = 1 
     flagBound = .FALSE.
     
     DO j = 1, Np_b

       DO m = 1, nparts

         IF (p == m  .OR.  SIZE(jFringes_FV(p,m)%vec) == 0)  CYCLE
      
         IF (partJ(jd_jb(j)) == p  .OR.  ANY(jFringes_FV(p,m)%vec == jd_jb(j))) THEN

           boundIdx = bound_p(j)
	   
	   IF (.NOT. flagBound(boundIdx)) THEN
             boundPart_FV(p) % vec(k) = boundIdx
	     flagBound(boundIdx) = .TRUE.
	     k = k + 1
	   ENDIF  
       
         ENDIF
	 
       ENDDO
       
     ENDDO
     
   ENDDO





 
   ! Saving Partitions
   PRINT*, '   ...saving partitions'

   OPEN (UNIT=idf, FILE='partitions.'//TRIM(grid_name))
     
     WRITE(idf,*) '-----------------------------------------------------------------------'
     WRITE(idf,*) ' GRID DIMENSIONS'
     WRITE(idf,*) '-----------------------------------------------------------------------'     
     WRITE(idf,*)  Np_d, Np_b, ' Nodes'
     WRITE(idf,*)  Ne_d, Ne_b, ' Elements'
     WRITE(idf,*)  Nc_d, Nc_b, ' FE node-pairs'
     WRITE(idf,*)  Nc_fv, ' FV node-pairs'
     WRITE(idf,*)  MAXVAL(bound_p), ' Number of boundaries'
     WRITE(idf,*)  nparts, ' Partitions'
     WRITE(idf,*) '-----------------------------------------------------------------------'
     
     WRITE(idf,*) ''
     
     DO i = 1, nparts
     
       WRITE(idf,*) '-----------------------------------------------------------------------'
       WRITE(idf,*) i, 'INDEX OF PARTITION'
       WRITE(idf,*) '', COUNT(partJ == i), ' nodes'
       DO j = 1, nparts         
         IF (i == j)  CYCLE
	 WRITE(idf,*) '', SIZE(jFringes_FE(i,j)%vec), ' nodes in FE notation from partition', j
       ENDDO
       DO j = 1, nparts
         IF (i == j)  CYCLE
	 WRITE(idf,*) '', SIZE(jFringes_FV(i,j)%vec), ' nodes in FV notation from partition', j
       ENDDO              
       
       WRITE(idf,*) '', SIZE(mFringes(i)%vec), ' elements'

       ! Finite elements              
       WRITE(idf,*) '', COUNT(partC_FE == i), ' FE node-pairs'
       DO j = 1, nparts       
         IF (i == j)  CYCLE
	 WRITE(idf,*) '', SIZE(cFringes_FE(i,j)%vec), ' FE node-pairs from partition', j
       ENDDO
       WRITE(idf,*) '', SIZE(cE_FE(i)%vec), '  FE extended node-pairs'

       ! Finite Volumes       
       WRITE(idf,*) '', COUNT(partC_FV == i), ' FV node-pairs'
       DO j = 1, nparts
         IF (i == j)  CYCLE   
         WRITE(idf,*) '', SIZE(cFringes_FV(i,j)%vec), ' FV node-pairs from partition', j
       ENDDO
       WRITE(idf,*) '', SIZE(cE_FV(i)%vec), '  FV extended node-pairs'
       WRITE(idf,*) '', SIZE(boundPart_FE(i) % vec), 'Number of boundaries in FE'
       WRITE(idf,*) '', SIZE(boundPart_FV(i) % vec), 'Number of boundaries in FV'
       WRITE(idf,*) '-----------------------------------------------------------------------'

       WRITE(idf,*) ' Nodes:'
       WRITE(idf,*) '-----------------------------------------------------------------------'
           
       DO j = 1, Np_d
         IF (partJ(j) == i  .AND.        ANY(jd_jb == j)) WRITE(idf,*) j, 'B'
         IF (partJ(j) == i  .AND.  .NOT. ANY(jd_jb == j)) WRITE(idf,*) j, 'D'
       ENDDO

       DO j = 1, nparts
         IF (i == j)  CYCLE
	 IF (SIZE(jFringes_FE(i,j)%vec)/= 0) THEN
           WRITE(idf,*) 'nodes from partition', j, 'in FE notation'
           DO k = 1, SIZE(jFringes_FE(i,j)%vec)
             IF (ANY(jd_jb == jFringes_FE(i,j)%vec(k))) THEN
	       WRITE(idf,*) jFringes_FE(i,j)%vec(k), 'B'
	     ELSE
	       WRITE(idf,*) jFringes_FE(i,j)%vec(k), 'D'
	     ENDIF  
           ENDDO
	 ENDIF  
       ENDDO

       DO j = 1, nparts
         IF (i == j)  CYCLE
	 IF (SIZE(jFringes_FV(i,j)%vec)/= 0) THEN
           WRITE(idf,*) 'nodes from partition', j, 'in FV notation'
           DO k = 1, SIZE(jFringes_FV(i,j)%vec)
             IF (ANY(jd_jb == jFringes_FV(i,j)%vec(k))) THEN
	       WRITE(idf,*) jFringes_FV(i,j)%vec(k), 'B'
	     ELSE
	       WRITE(idf,*) jFringes_FV(i,j)%vec(k), 'D'
	     ENDIF  
           ENDDO
         ENDIF
       ENDDO

       WRITE(idf,*) '-----------------------------------------------------------------------'
       WRITE(idf,*) ' Elements:'    
       WRITE(idf,*) '-----------------------------------------------------------------------'       
       DO j = 1, SIZE(mFringes(i)%vec)
         WRITE(idf,*) mFringes(i)% vec(j)
       ENDDO

       WRITE(idf,*) '-----------------------------------------------------------------------'       
       WRITE(idf,*) ' FINITE ELEMENTS node_pairs:'
       WRITE(idf,*) '-----------------------------------------------------------------------'
       DO j = 1, Nc_d
         IF (partC_FE(j) == i  .AND.        ANY(cd_cb == j))  WRITE(idf,*) j, 'B'
	 IF (partC_FE(j) == i  .AND.  .NOT. ANY(cd_cb == j))  WRITE(idf,*) j, 'D'
       ENDDO

       DO j = 1, nparts       
         IF (i == j)  CYCLE
	 IF (SIZE(cFringes_FE(i,j)%vec) /= 0) THEN
           WRITE(idf,*) 'node_pairs from partition', j    
           DO k = 1, SIZE(cFringes_FE(i,j)%vec)

             IF (ANY(cd_cb == cFringes_FE(i,j)%vec(k))) THEN
	       WRITE(idf,*) cFringes_FE(i,j)%vec(k), 'B'
	     ELSE
	       WRITE(idf,*) cFringes_FE(i,j)%vec(k), 'D'
	     ENDIF  
	     
           ENDDO
	 ENDIF  
       ENDDO
       
       WRITE(idf,*) ' extended node-pairs'       
       DO j = 1, SIZE(cE_FE(i)%vec)
       
         IF (ANY(cd_cb == cE_FE(i)%vec(j))) THEN
	   WRITE(idf,*) cE_FE(i)%vec(j), 'B'
	 ELSE
	   WRITE(idf,*) cE_FE(i)%vec(j), 'D'
	 ENDIF  
                    
       ENDDO



       WRITE(idf,*) '-----------------------------------------------------------------------'
       WRITE(idf,*) ' FINITE VOLUMES node_pairs:'
       WRITE(idf,*) '-----------------------------------------------------------------------'       
       DO j = 1, Nc_fv
         IF (partC_FV(j) == i  .AND.        ANY(cd_cb == j))  WRITE(idf,*) j, 'B'
	 IF (partC_FV(j) == i  .AND.  .NOT. ANY(cd_cb == j))  WRITE(idf,*) j, 'D'
       ENDDO

       DO j = 1, nparts       
         IF (i == j)  CYCLE       
	 IF (SIZE(cFringes_FV(i,j)%vec) /= 0) THEN
           WRITE(idf,*) 'node_pairs from partition', j    
           DO k = 1, SIZE(cFringes_FV(i,j)%vec)
	   
             IF (ANY(cd_cb == cFringes_FV(i,j)%vec(k))) THEN
	       WRITE(idf,*) cFringes_FV(i,j)%vec(k), 'B'
	     ELSE
	       WRITE(idf,*) cFringes_FV(i,j)%vec(k), 'D'
	     ENDIF  
           
	   ENDDO
         ENDIF
       ENDDO

       WRITE(idf,*) ' extended node-pairs'
       DO j = 1, SIZE(cE_FV(i)%vec)
       
         IF (ANY(cd_cb == cE_FV(i)%vec(j))) THEN
	   WRITE(idf,*) cE_FV(i)%vec(j), 'B'
	 ELSE
	   WRITE(idf,*) cE_FV(i)%vec(j), 'D'
	 ENDIF  
	               
       ENDDO

       WRITE(idf,*) '-----------------------------------------------------------------------'
       WRITE(idf,*) ' Boundaries:'
       WRITE(idf,*) '-----------------------------------------------------------------------'       
       WRITE(idf,*) 'boundaries for finite elements notation'
       DO j = 1, SIZE(boundPart_FE(i)%vec)
         WRITE(idf,*) boundPart_FE(i)%vec(j)
       ENDDO
       WRITE(idf,*) 'boundaries for finite volumes notation'
       DO j = 1, SIZE(boundPart_FV(i)%vec)
         WRITE(idf,*) boundPart_FV(i)%vec(j)
       ENDDO
     
     ENDDO
       
   CLOSE(idf)
 
   ! Plotting partitions
   PRINT*, '   ...saving plots'
 
    OPEN (UNIT=idf, FILE='partitions.'//TRIM(grid_name)//'.mtv')
 
     ! No Contour
     WRITE(idf,*) '$ DATA=CURVE2D'
     WRITE(idf,*) '% toplabel = "Partitions"'
     WRITE(idf,*) '% equalscale = true'
     WRITE(idf,*) ''
      
      ! Mesh
      DO m = 1, Nc_fv
        WRITE(idf,*) '% linecolor = 4'
        WRITE(idf,*) rr(1,j_c_fv(1,m)), rr(2,j_c_fv(1,m)), j_c_fv(1,m)
        WRITE(idf,*) rr(1,j_c_fv(2,m)), rr(2,j_c_fv(2,m)), j_c_fv(2,m)
        WRITE(idf,*) ''
      ENDDO

      ! Boundaries of the partitions
      DO m = 1, Ne_d
      
        IF (ALLOCATED(fC)) DEALLOCATE(fC) 
        ALLOCATE (fC(SIZE(c_m_fv(m)%vec)))
        
        fC = c_m_fv(m)%vec
        
        DO n = 1, nparts
          DO p = 1, nparts
            IF (n /= p) THEN
        
              nF = 0
              DO c = 1, SIZE(fC)
                IF (ANY(fC(c) == cFringes_FV(n,p)%vec)) nF = nF + 1
              ENDDO
        
              IF (nF == 2) THEN
                
                IF (ALLOCATED(xm)) DEALLOCATE (xm, ym)
                ALLOCATE (xm(2), ym(2))
                
                xC = 0
                yC = 0
                q = 1
                
                DO c = 1, SIZE(fC)
                  IF (ANY(fC(c) == cFringes_FV(n,p)%vec)) THEN
                  
                    x1 = rr(1, j_c_fv(1,fC(c)))
                    x2 = rr(1, j_c_fv(2,fC(c)))
                  
                    y1 = rr(2, j_c_fv(1,fC(c)))
                    y2 = rr(2, j_c_fv(2,fC(c)))

                    xm(q) = (x1 + x2)/2.0
                    ym(q) = (y1 + y2)/2.0
                    
                    q = q + 1
                    
                  ENDIF
                ENDDO
                
                IF (ele_type_d(m) == 2) THEN
                
                  WRITE(idf,*) '% linetype = 1'
                  WRITE(idf,*) '% markertype = 0'
                  WRITE(idf,*) '% linecolor = 0'
                  WRITE(idf,*) xm(1), ym(1)
                  WRITE(idf,*) xm(2), ym(2)
                  WRITE(idf,*)                   
                
                ELSE
                 
                  xC = SUM(rr(1,j_m_d(m)%vec)) / SIZE(j_m_d(m)%vec)
                  yC = SUM(rr(2,j_m_d(m)%vec)) / SIZE(j_m_d(m)%vec)

                  WRITE(idf,*) '% linetype = 1'
                  WRITE(idf,*) '% markertype = 0'
                  WRITE(idf,*) '% linecolor = 0'
                  WRITE(idf,*) xm(1), ym(1)
                  WRITE(idf,*) xC, yC
                  WRITE(idf,*)
 
                  WRITE(idf,*) '% linetype = 1'
                  WRITE(idf,*) '% markertype = 0'
                  WRITE(idf,*) '% linecolor = 0'
                  WRITE(idf,*) xm(2), ym(2)
                  WRITE(idf,*) xC, yC
                  WRITE(idf,*)
                
                ENDIF
              ENDIF
        
            ENDIF
          ENDDO
        ENDDO
        
      ENDDO !(m loop)

     DO m = 1, Ne_d

       IF (ele_type_d(m) == 2) THEN
    
          IF (partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(2)) .AND.  &
              partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(3)) .AND.  &
              partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(3)) ) THEN

            IF (ALLOCATED(xm)) DEALLOCATE (xm, ym)
            ALLOCATE (xm(3), ym(3))
 
            DO c = 1, 3
              xm(c) = rr(1, j_m_d(m)%vec(c))
              ym(c) = rr(2, j_m_d(m)%vec(c))
            ENDDO

            xC = SUM(xm)/SIZE(xm)
            yC = SUM(ym)/SIZE(ym)
 
            DO c = 1, 3
 
              x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
              y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

              x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
              y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
 
              xm(c) = (x1+x2)/2.0
              ym(c) = (y1+y2)/2.0
 
            ENDDO
 
            WRITE(idf,*)
            WRITE(idf,*) '% linetype = 1'
            WRITE(idf,*) '% markertype = 0'
            WRITE(idf,*) '% linecolor = 0'
            WRITE(idf,*) xm(1), ym(1)
            WRITE(idf,*) xC, yC
            WRITE(idf,*)
            WRITE(idf,*) '% linetype = 1'
            WRITE(idf,*) '% markertype = 0'
            WRITE(idf,*) '% linecolor = 0'
            WRITE(idf,*) xm(2), ym(2)
            WRITE(idf,*) xC, yC
            WRITE(idf,*)
            WRITE(idf,*) '% linetype = 1'
            WRITE(idf,*) '% markertype = 0'
            WRITE(idf,*) '% linecolor = 0'
            WRITE(idf,*) xm(3), ym(3)
            WRITE(idf,*) xC, yC
            WRITE(idf,*)

            IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)

          ENDIF
       ENDIF

       ! QUADRILATERS
       IF (ele_type_d(m) == 3) THEN
    
           ! Each node of the element belongs to a different partition
           IF (partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(2)) .AND.  & 
               partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(3)) .AND.  & 
               partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(4)) .AND.  & 
               partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(3)) .AND.  & 
               partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(4)) .AND.  & 
               partJ(j_m_d(m)%vec(3)) /= partJ(j_m_d(m)%vec(4))) THEN    

             ALLOCATE (xm(4), ym(4))                                     
                                                                         
             DO c = 1, 4                                                 
               xm(c) = rr(1, j_m_d(m)%vec(c))                            
               ym(c) = rr(2, j_m_d(m)%vec(c))                            
             ENDDO                                                       

             xC = SUM(xm)/SIZE(xm)                                       
             yC = SUM(ym)/SIZE(ym)                                       
                                                                         
             DO c = 1, 4                                                 
                                                                         
               x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))                     
               y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))                     

               x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))                     
               y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))                     
                                                                         
               xm(c) = (x1+x2)/2.0                                       
               ym(c) = (y1+y2)/2.0                                       
                                                                         
             ENDDO                                                       
                                                                         
             WRITE(idf,*)                                                
             WRITE(idf,*) '% linetype = 1'                               
             WRITE(idf,*) '% markertype = 0'                             
             WRITE(idf,*) '% linecolor = 0'                              
             WRITE(idf,*) xm(1), ym(1)                                   
             WRITE(idf,*) xC, yC                                         
             WRITE(idf,*)                                                
             WRITE(idf,*) '% linetype = 1'                               
             WRITE(idf,*) '% markertype = 0'                             
             WRITE(idf,*) '% linecolor = 0'                              
             WRITE(idf,*) xm(2), ym(2)                                   
             WRITE(idf,*) xC, yC                                         
             WRITE(idf,*)                                                
             WRITE(idf,*) '% linetype = 1'                               
             WRITE(idf,*) '% markertype = 0'                             
             WRITE(idf,*) '% linecolor = 0'                              
             WRITE(idf,*) xm(3), ym(3)                                   
             WRITE(idf,*) xC, yC                                         
             WRITE(idf,*)                                                
             WRITE(idf,*) '% linetype = 1'                               
             WRITE(idf,*) '% markertype = 0'                             
             WRITE(idf,*) '% linecolor = 0'                              
             WRITE(idf,*) xm(4), ym(4)                                   
             WRITE(idf,*) xC, yC                                         
             WRITE(idf,*)                                                

             IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                       

          ENDIF                                                         

          ! Count the number of nodes with different affiliations 
          nF = 0
          DO c = 1, 4
            IF (partJ(j_c_fv(1,c_m_fv(m)%vec(c))) /=  partJ(j_c_fv(2,c_m_fv(m)%vec(c)))) nF = nF + 1 
          ENDDO
       
          ! Three nodes belong to different partitions
          IF (nF == 3) THEN
            
            IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 
            ALLOCATE (xm(3), ym(3))
            
            n = 1
            DO c = 1, 4
              IF (partJ(j_c_fv(1,c_m_fv(m)%vec(c))) /=  partJ(j_c_fv(2,c_m_fv(m)%vec(c)))) THEN
              
                x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
                y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

                x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
                y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
            
                xm(n) = (x1+x2)/2.0
                ym(n) = (y1+y2)/2.0
               
                n = n + 1
               
              ENDIF
            ENDDO
            
            xC = SUM(rr(1,j_m_d(m)%vec))/SIZE(j_m_d(m)%vec)
            yC = SUM(rr(2,j_m_d(m)%vec))/SIZE(j_m_d(m)%vec) 

            
            WRITE(idf,*)                                          
            WRITE(idf,*) '% linetype = 1'                         
            WRITE(idf,*) '% markertype = 0'                       
            WRITE(idf,*) '% linecolor = 0'                        
            WRITE(idf,*) xm(1), ym(1)                             
            WRITE(idf,*) xC, yC                                   
            WRITE(idf,*)                                          
            WRITE(idf,*) '% linetype = 1'                         
            WRITE(idf,*) '% markertype = 0'                       
            WRITE(idf,*) '% linecolor = 0'                        
            WRITE(idf,*) xm(2), ym(2)                             
            WRITE(idf,*) xC, yC                                   
            WRITE(idf,*)                                          
            WRITE(idf,*) '% linetype = 1'                         
            WRITE(idf,*) '% markertype = 0'                       
            WRITE(idf,*) '% linecolor = 0'                        
            WRITE(idf,*) xm(3), ym(3)                             
            WRITE(idf,*) xC, yC                                   
            WRITE(idf,*)                                          

            IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)
       
          ENDIF
       ENDIF                                           
     ENDDO

     WRITE(idf,*) '$ DATA=CONTCURVE'
     WRITE(idf,*) '% toplabel = "Partitions Numbering"'   
     WRITE(idf,*) '% MESHPLOT = false'
     WRITE(idf,*) '% CONTSTYLE = 2'
     WRITE(idf,*) '% nsteps = 50'
     
     DO i = 1, Ne_d
       DO j = 1, SIZE(j_m_d(i)%vec)
         WRITE(idf,*) rr(1,j_m_d(i)%vec(j)), rr(2,j_m_d(i)%vec(j)), &
                      partJ(j_m_d(i)%vec(j)), j_m_d(i)%vec(j)
       ENDDO
       WRITE(idf,*) ''
     ENDDO
      
   CLOSE(idf)

!   OPEN (UNIT=idf, FILE='parts-Elems.'//TRIM(grid_name)//'.mtv')
!      
!     DO i = 1, nparts
!     
!       WRITE(idf,*) '$ DATA=CONTCURVE'
!       WRITE(idf,*) '% MESHPLOT = TRUE'
!       WRITE(idf,*) '% CONTSTYLE = 2'
!       WRITE(idf,*) '% nsteps = 50'
!       WRITE(idf,*) ''
!       
!       DO j = 1, SIZE(mFringes(i)%vec)
!         DO k = 1, SIZE(j_m_d(mFringes(i)%vec(j))%vec)
!           WRITE(idf,*) rr(1,j_m_d(mFringes(i)%vec(j))%vec(k)), &
!	                rr(2,j_m_d(mFringes(i)%vec(j))%vec(k)), &
!			i, j_m_d(mFringes(i)%vec(j))%vec(k)
!         ENDDO
!	 WRITE(idf,*)
!       ENDDO
!     ENDDO
!     
!     WRITE(idf,*)
!       
!   CLOSE(idf)

   IF (debug) THEN

     OPEN (UNIT=idf, FILE='interface-FE.'//TRIM(grid_name)//'.mtv')

       DO i = 1, nparts
         
	 WRITE(idf,*) ''
         WRITE(idf,*) '$ DATA=CURVE2D'
         WRITE(idf,*) '% toplabel = "Interface - Finite Element notation"'
         WRITE(idf,*) '% subtitle = "Partition', i,'"'
         WRITE(idf,*) '% equalscale = true'
         WRITE(idf,*) ''
         
         ! Mesh
         DO m = 1, Nc_fv
           IF (ALL(partJ(j_c_fv(1:2,m)) == i)) THEN
	     WRITE(idf,*) '% linetype = 1'
	   ELSE
	     WRITE(idf,*) '% linetype = 3'	   
	   ENDIF  
           WRITE(idf,*) '% linecolor = 1'
           WRITE(idf,*) rr(1,j_c_fv(1,m)), rr(2,j_c_fv(1,m)), j_c_fv(1,m)
           WRITE(idf,*) rr(1,j_c_fv(2,m)), rr(2,j_c_fv(2,m)), j_c_fv(2,m)
           WRITE(idf,*) ''
         ENDDO


         DO j = 1, nparts
         
	   IF (i == j)  CYCLE
	 
           DO k = 1, SIZE(cFringes_FE(i,j)%vec)
            
             c = cFringes_FE(i,j)%vec(k)

             WRITE(idf,*) '% linewidth = 2'			                    
             WRITE(idf,*) '% linecolor = 5'
             WRITE(idf,*) rr(1,j_c_d(1,c)), rr(2,j_c_d(1,c)), j_c_d(1,c)
             WRITE(idf,*) rr(1,j_c_d(2,c)), rr(2,j_c_d(2,c)), j_c_d(2,c)
             WRITE(idf,*)
            
           ENDDO

         ENDDO



         DO j = 1, SIZE(cE_FE(i)%vec)				        
         							        
           c = cE_FE(i)%vec(j)  				        
         							        
           WRITE(idf,*) '% linewidth = 2'			        
           WRITE(idf,*) '% linecolor = 4'			        
           WRITE(idf,*) rr(1,j_c_d(1,c)), rr(2,j_c_d(1,c)), j_c_d(1,c)  
           WRITE(idf,*) rr(1,j_c_d(2,c)), rr(2,j_c_d(2,c)), j_c_d(2,c)  
           WRITE(idf,*) 					        
         							        
         ENDDO  						        



         DO j = 1, nparts
           IF (i /= j) THEN

             ! Fringe nodes
             DO m = 1, SIZE(jFringes_FE(i,j)%vec)
               WRITE(idf,*) '% markercolor = 0' 
               WRITE(idf,*) '% markertype = 12'
               WRITE(idf,*) '% linetype = 0'
               WRITE(idf,*) rr(1,jFringes_FE(i,j)%vec(m)), rr(2,jFringes_FE(i,j)%vec(m)), jFringes_FE(i,j)%vec(m)
               WRITE(idf,*) ''
             ENDDO

             ! Boundaries of the partitions
             DO m = 1, Ne_d
             
               IF (ALLOCATED(fC)) DEALLOCATE(fC) 
               ALLOCATE (fC(SIZE(c_m_fv(m)%vec)))
               
               fC = c_m_fv(m)%vec
               
               DO n = 1, nparts
                 DO p = 1, nparts
                   IF (n /= p) THEN
               
                     nF = 0
                     DO c = 1, SIZE(fC)
                       IF (ANY(fC(c) == cFringes_FV(n,p)%vec)) nF = nF + 1
                     ENDDO
               
                     IF (nF == 2) THEN
                       
                       IF (ALLOCATED(xm)) DEALLOCATE (xm, ym)
                       ALLOCATE (xm(2), ym(2))
                       
                       xC = 0
                       yC = 0
                       q = 1
                       
                       DO c = 1, SIZE(fC)
                         IF (ANY(fC(c) == cFringes_FV(n,p)%vec)) THEN
                         
                           x1 = rr(1, j_c_fv(1,fC(c)))
                           x2 = rr(1, j_c_fv(2,fC(c)))
                         
                           y1 = rr(2, j_c_fv(1,fC(c)))
                           y2 = rr(2, j_c_fv(2,fC(c)))

                           xm(q) = (x1 + x2)/2.0
                           ym(q) = (y1 + y2)/2.0
                           
                           q = q + 1
                           
                         ENDIF
                       ENDDO
                       
                       IF (ele_type_d(m) == 2) THEN
                       
                         WRITE(idf,*) '% linetype = 1'
                         WRITE(idf,*) '% markertype = 0'
                         WRITE(idf,*) '% linecolor = 0'
                         WRITE(idf,*) xm(1), ym(1)
                         WRITE(idf,*) xm(2), ym(2)
                         WRITE(idf,*)                   
                       
                       ELSE
                        
                         xC = SUM(rr(1,j_m_d(m)%vec)) / SIZE(j_m_d(m)%vec)
                         yC = SUM(rr(2,j_m_d(m)%vec)) / SIZE(j_m_d(m)%vec)

                         WRITE(idf,*) '% linetype = 1'
                         WRITE(idf,*) '% markertype = 0'
                         WRITE(idf,*) '% linecolor = 0'
                         WRITE(idf,*) xm(1), ym(1)
                         WRITE(idf,*) xC, yC
                         WRITE(idf,*)
 
                         WRITE(idf,*) '% linetype = 1'
                         WRITE(idf,*) '% markertype = 0'
                         WRITE(idf,*) '% linecolor = 0'
                         WRITE(idf,*) xm(2), ym(2)
                         WRITE(idf,*) xC, yC
                         WRITE(idf,*)
                       
                       ENDIF
                     ENDIF
               
                   ENDIF
                 ENDDO
               ENDDO
               
             ENDDO !(m loop)  

       DO m = 1, Ne_d

         IF (ele_type_d(m) == 2) THEN
      
         IF (partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(2)) .AND.  &
             partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(3)) .AND.  &
             partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(3)) ) THEN

           IF (ALLOCATED(xm)) DEALLOCATE (xm, ym)
           ALLOCATE (xm(3), ym(3))                               
                                                                 
           DO c = 1, 3                                    
             xm(c) = rr(1, j_m_d(m)%vec(c))
             ym(c) = rr(2, j_m_d(m)%vec(c))
           ENDDO                                                 

           xC = SUM(xm)/SIZE(xm)                                 
           yC = SUM(ym)/SIZE(ym)                                 
           
           DO c = 1, 3
           
             x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
             y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

             x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
             y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
           
             xm(c) = (x1+x2)/2.0
             ym(c) = (y1+y2)/2.0
                    
           ENDDO
           
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(1), ym(1)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(2), ym(2)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(3), ym(3)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          

           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 

         ENDIF 
         ENDIF



         IF (ele_type_d(m) == 3) THEN
      
         IF (partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(2)) .AND.  &
             partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(3)) .AND.  &
             partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(4)) .AND.  &
             partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(3)) .AND.  &
             partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(4)) .AND.  &
             partJ(j_m_d(m)%vec(3)) /= partJ(j_m_d(m)%vec(4))) THEN

           ALLOCATE (xm(4), ym(4))                               
                                                                 
           DO c = 1, 4                                    
             xm(c) = rr(1, j_m_d(m)%vec(c))
             ym(c) = rr(2, j_m_d(m)%vec(c))
           ENDDO                                                 

           xC = SUM(xm)/SIZE(xm)                                 
           yC = SUM(ym)/SIZE(ym)                                 
           
           DO c = 1, 4
           
             x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
             y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

             x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
             y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
           
             xm(c) = (x1+x2)/2.0
             ym(c) = (y1+y2)/2.0
                    
           ENDDO
           
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(1), ym(1)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(2), ym(2)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(3), ym(3)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(4), ym(4)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          

           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 

         ENDIF 
         
         nF = 0
         DO c = 1, 4
           IF (partJ(j_c_fv(1,c_m_fv(m)%vec(c))) /=  partJ(j_c_fv(2,c_m_fv(m)%vec(c)))) nF = nF + 1
         ENDDO
         
         IF (nF == 3) THEN
           
           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 
           ALLOCATE (xm(3), ym(3))
           
           n = 1
           DO c = 1, 4
             IF (partJ(j_c_fv(1,c_m_fv(m)%vec(c))) /=  partJ(j_c_fv(2,c_m_fv(m)%vec(c)))) THEN
             
               x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
               y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

               x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
               y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
           
               xm(n) = (x1+x2)/2.0
               ym(n) = (y1+y2)/2.0
              
               n = n + 1
              
             ENDIF
           ENDDO
           
            xC = SUM(rr(1,j_m_d(m)%vec))/SIZE(j_m_d(m)%vec)
            yC = SUM(rr(2,j_m_d(m)%vec))/SIZE(j_m_d(m)%vec) 

           
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(1), ym(1)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(2), ym(2)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(3), ym(3)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          

           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 
           
           
         ENDIF
         
         
         ENDIF

                                                            
       ENDDO


           ENDIF !(if i/=j)
         ENDDO !(j loop)
       ENDDO  !(i loop)

     CLOSE(idf)



     OPEN (UNIT=idf, FILE='interface-FV.'//TRIM(grid_name)//'.mtv')

       DO i = 1, nparts
              
         WRITE(idf,*) '$ DATA=CURVE2D'
         WRITE(idf,*) '% toplabel = "Interface - Finite Volumes notation"'
         WRITE(idf,*) '% subtitle = "Partition', i, '"'
         WRITE(idf,*) '% equalscale = true'
         WRITE(idf,*) ''
        
         ! Mesh
         DO m = 1, Nc_fv
           IF (ALL(partJ(j_c_fv(1:2,m)) == i)) THEN
	     WRITE(idf,*) '% linetype = 1'
	   ELSE
	     WRITE(idf,*) '% linetype = 3'	   
	   ENDIF
           WRITE(idf,*) rr(1,j_c_fv(1,m)), rr(2,j_c_fv(1,m)), j_c_fv(1,m)
           WRITE(idf,*) rr(1,j_c_fv(2,m)), rr(2,j_c_fv(2,m)), j_c_fv(2,m)
           WRITE(idf,*) ''
         ENDDO
       
       
         DO j = 1, nparts
         
	   IF (i == j)  CYCLE
	 
           DO k = 1, SIZE(cFringes_FV(i,j)%vec)
            
             c = cFringes_FV(i,j)%vec(k)

             WRITE(idf,*) '% linewidth = 2'
             WRITE(idf,*) '% linecolor = 5'		    
             WRITE(idf,*) rr(1,j_c_fv(1,c)), rr(2,j_c_fv(1,c)), j_c_fv(1,c)
             WRITE(idf,*) rr(1,j_c_fv(2,c)), rr(2,j_c_fv(2,c)), j_c_fv(2,c)
             WRITE(idf,*)
            
           ENDDO

         ENDDO



         DO j = 1, SIZE(cE_FV(i)%vec)				        
         							        
           c = cE_FV(i)%vec(j)  				        
         							        
           WRITE(idf,*) '% linewidth = 2'			        
           WRITE(idf,*) '% linecolor = 4'			        
           WRITE(idf,*) rr(1,j_c_fv(1,c)), rr(2,j_c_fv(1,c)), j_c_fv(1,c)  
           WRITE(idf,*) rr(1,j_c_fv(2,c)), rr(2,j_c_fv(2,c)), j_c_fv(2,c)  
           WRITE(idf,*) 					        
         							        
         ENDDO
              
       
       
         DO j = 1, nparts	 
           IF (i /= j) THEN

             ! Fringe nodes
             DO m = 1, SIZE(jFringes_FV(i,j)%vec)
               WRITE(idf,*) '% markercolor = 0' 
               WRITE(idf,*) '% markertype = 12'
               WRITE(idf,*) '% linetype = 0'
               WRITE(idf,*) rr(1,jFringes_FV(i,j)%vec(m)), rr(2,jFringes_FV(i,j)%vec(m)), jFringes_FV(i,j)%vec(m)
               WRITE(idf,*) ''
             ENDDO

             ! Boundaries of the partitions
             DO m = 1, Ne_d
             
               IF (ALLOCATED(fC)) DEALLOCATE(fC) 
               ALLOCATE (fC(SIZE(c_m_fv(m)%vec)))
               
               fC = c_m_fv(m)%vec
               
               DO n = 1, nparts
                 DO p = 1, nparts
                   IF (n /= p) THEN
               
                     nF = 0
                     DO c = 1, SIZE(fC)
                       IF (ANY(fC(c) == cFringes_FV(n,p)%vec)) nF = nF + 1
                     ENDDO
               
                     IF (nF == 2) THEN
                       
                       IF (ALLOCATED(xm)) DEALLOCATE (xm, ym)
                       ALLOCATE (xm(2), ym(2))
                       
                       xC = 0
                       yC = 0
                       q = 1
                       
                       DO c = 1, SIZE(fC)
                         IF (ANY(fC(c) == cFringes_FV(n,p)%vec)) THEN
                         
                           x1 = rr(1, j_c_fv(1,fC(c)))
                           x2 = rr(1, j_c_fv(2,fC(c)))
                         
                           y1 = rr(2, j_c_fv(1,fC(c)))
                           y2 = rr(2, j_c_fv(2,fC(c)))

                           xm(q) = (x1 + x2)/2.0
                           ym(q) = (y1 + y2)/2.0
                           
                           q = q + 1
                           
                         ENDIF
                       ENDDO
                       
                       IF (ele_type_d(m) == 2) THEN
                       
                         WRITE(idf,*) '% linetype = 1'
                         WRITE(idf,*) '% markertype = 0'
                         WRITE(idf,*) '% linecolor = 0'
                         WRITE(idf,*) xm(1), ym(1)
                         WRITE(idf,*) xm(2), ym(2)
                         WRITE(idf,*)                   
                       
                       ELSE
                        
                         xC = SUM(rr(1,j_m_d(m)%vec)) / SIZE(j_m_d(m)%vec)
                         yC = SUM(rr(2,j_m_d(m)%vec)) / SIZE(j_m_d(m)%vec)

                         WRITE(idf,*) '% linetype = 1'
                         WRITE(idf,*) '% markertype = 0'
                         WRITE(idf,*) '% linecolor = 0'
                         WRITE(idf,*) xm(1), ym(1)
                         WRITE(idf,*) xC, yC
                         WRITE(idf,*)
 
                         WRITE(idf,*) '% linetype = 1'
                         WRITE(idf,*) '% markertype = 0'
                         WRITE(idf,*) '% linecolor = 0'
                         WRITE(idf,*) xm(2), ym(2)
                         WRITE(idf,*) xC, yC
                         WRITE(idf,*)
                       
                       ENDIF
                     ENDIF
               
                   ENDIF
                 ENDDO
               ENDDO
               
             ENDDO !(m loop)  

       DO m = 1, Ne_d

         IF (ele_type_d(m) == 2) THEN
      
         IF (partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(2)) .AND.  &
             partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(3)) .AND.  &
             partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(3)) ) THEN

           IF (ALLOCATED(xm)) DEALLOCATE (xm, ym)
           ALLOCATE (xm(3), ym(3))                               
                                                                 
           DO c = 1, 3                                    
             xm(c) = rr(1, j_m_d(m)%vec(c))
             ym(c) = rr(2, j_m_d(m)%vec(c))
           ENDDO                                                 

           xC = SUM(xm)/SIZE(xm)                                 
           yC = SUM(ym)/SIZE(ym)                                 
           
           DO c = 1, 3
           
             x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
             y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

             x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
             y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
           
             xm(c) = (x1+x2)/2.0
             ym(c) = (y1+y2)/2.0
                    
           ENDDO
           
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(1), ym(1)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(2), ym(2)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(3), ym(3)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          

           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 

         ENDIF 
         ENDIF



         IF (ele_type_d(m) == 3) THEN
      
         IF (partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(2)) .AND.  &
             partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(3)) .AND.  &
             partJ(j_m_d(m)%vec(1)) /= partJ(j_m_d(m)%vec(4)) .AND.  &
             partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(3)) .AND.  &
             partJ(j_m_d(m)%vec(2)) /= partJ(j_m_d(m)%vec(4)) .AND.  &
             partJ(j_m_d(m)%vec(3)) /= partJ(j_m_d(m)%vec(4))) THEN

           ALLOCATE (xm(4), ym(4))                               
                                                                 
           DO c = 1, 4                                    
             xm(c) = rr(1, j_m_d(m)%vec(c))
             ym(c) = rr(2, j_m_d(m)%vec(c))
           ENDDO                                                 

           xC = SUM(xm)/SIZE(xm)                                 
           yC = SUM(ym)/SIZE(ym)                                 
           
           DO c = 1, 4
           
             x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
             y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

             x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
             y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
           
             xm(c) = (x1+x2)/2.0
             ym(c) = (y1+y2)/2.0
                    
           ENDDO
           
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(1), ym(1)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(2), ym(2)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(3), ym(3)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(4), ym(4)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          

           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 

         ENDIF 
         
         nF = 0
         DO c = 1, 4
           IF (partJ(j_c_fv(1,c_m_fv(m)%vec(c))) /=  partJ(j_c_fv(2,c_m_fv(m)%vec(c)))) nF = nF + 1
         ENDDO
         
         IF (nF == 3) THEN
           
           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 
           ALLOCATE (xm(3), ym(3))
           
           n = 1
           DO c = 1, 4
             IF (partJ(j_c_fv(1,c_m_fv(m)%vec(c))) /=  partJ(j_c_fv(2,c_m_fv(m)%vec(c)))) THEN
             
               x1 = rr(1,j_c_fv(1,c_m_fv(m)%vec(c)))
               y1 = rr(2,j_c_fv(1,c_m_fv(m)%vec(c)))

               x2 = rr(1,j_c_fv(2,c_m_fv(m)%vec(c)))
               y2 = rr(2,j_c_fv(2,c_m_fv(m)%vec(c)))
           
               xm(n) = (x1+x2)/2.0
               ym(n) = (y1+y2)/2.0
              
               n = n + 1
              
             ENDIF
           ENDDO
           
           xC = SUM(xm)/SIZE(xm)                                 
           yC = SUM(ym)/SIZE(ym)                                 

           
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(1), ym(1)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(2), ym(2)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          
           WRITE(idf,*) '% linetype = 1'                         
           WRITE(idf,*) '% markertype = 0'                       
           WRITE(idf,*) '% linecolor = 0'                        
           WRITE(idf,*) xm(3), ym(3)                             
           WRITE(idf,*) xC, yC                                   
           WRITE(idf,*)                                          

           IF (ALLOCATED(xm)) DEALLOCATE (xm,ym)                 
           
           
         ENDIF
         
         
         ENDIF

                                                            
       ENDDO

           ENDIF !(if i/=j)
         ENDDO !(j loop)
       ENDDO  !(i loop)

     CLOSE(idf)
    
    
    
    ENDIF
    
  END SUBROUTINE  partMesh
  
  
  END MODULE  mesh_partition
