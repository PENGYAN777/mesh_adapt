! ==============================================================================
!
!      Module: mp_interface
!
! Description: Module devoted to multi-processor execution.
!              It contains definitions and subroutine
!              to handle multi-processor execution of
!              the code. 
!              The approach to parallelization
!              is to have ONE master process that drives everything
!              and N worker processes. So in terms of number of 
!              processes we have a total of N+1 machines involved.
!
!      Author: Marco Fossati
!              Department of Aerospace Engineering
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
! ==============================================================================
 
   MODULE mp_interface
   
   !------------------------------------------------------------------
   USE commons
   USE dynamic_vector
   
   IMPLICIT NONE
   
   INTEGER :: iErr
   INTEGER :: MP_nProc, MP_wProc, MP_pRank

   INTEGER :: MP_COMM_ALL, MP_COMM_MASTER, MP_COMM_WORKERS
   
   INTEGER :: Nj_G, Njb_G, Nm_G, Nmb_G, NcFE_G, NcFV_G, NcFE_b, Nb_G
   INTEGER :: iBou = 0

   REAL(KIND=8) :: loc_minDt = 0.0, &
                   loc_maxDt = 0.0


   TYPE partition
    ! Number of nodes of the partition
    INTEGER :: Nj_P
    ! Partial to global nodes connectivity
    INTEGER, DIMENSION(:), POINTER :: jG_jP
   END TYPE partition

   TYPE(partition), DIMENSION(:), POINTER :: part


   !----------------------------------------------------------------------
   ! GRID OF THE I-TH PARTITION (data visible to the single worker)
   ! Total number of nodes of the partition 
   INTEGER :: NjTP, NbP, NbP_FE, NbP_FV
   
   ! Internal nodes and element to the partition
   INTEGER :: NjP, NjbP, NmP, NmbP
   
   ! Total number of node-pairs of the partition    
   INTEGER :: NcTFEP, NcTFVP   

   ! Number of internal node-pairs
   INTEGER :: NcFEP, NcFVP   
   
   ! Interface nodes with other partitions
   INTEGER, DIMENSION(:), ALLOCATABLE :: NjI
   INTEGER, DIMENSION(:), ALLOCATABLE :: NcIfe, NcIfv
   INTEGER, DIMENSION(:), ALLOCATABLE :: NcEXT_FE, NcEXT_FV
      
   ! Part-Global and Global-Part nodes connectivity
   INTEGER, DIMENSION(:), POINTER :: jG_jP
   INTEGER, DIMENSION(:), POINTER :: jP_jG
   INTEGER, DIMENSION(:), POINTER :: jP_jG_B
      
   ! Part-Global and Global-Part elements connectivity
   INTEGER, DIMENSION(:), POINTER :: mG_mP
   INTEGER, DIMENSION(:), POINTER :: mP_mG
   INTEGER, DIMENSION(:), POINTER :: mP_mG_B
      
   ! Part-Global and Global-Part FE node_pairs connectivity
   INTEGER, DIMENSION(:), POINTER :: cG_cP_FE
   INTEGER, DIMENSION(:), POINTER :: cP_cG_FE
   INTEGER, DIMENSION(:), POINTER :: cP_cG_FE_B

   ! Part-Global and Global-Part FV node_pairs connectivity
   INTEGER, DIMENSION(:), POINTER :: cG_cP_FV
   INTEGER, DIMENSION(:), POINTER :: cP_cG_FV
   
   ! Number of boundary nodes and node-pair
   INTEGER :: NjP_b, NcFEP_b

   ! Part-Global and Global-Part boundaries connectivity
   INTEGER, DIMENSION(:), POINTER :: bG_bP
   INTEGER, DIMENSION(:), POINTER :: bP_bG
   
   TYPE(D_I_V), DIMENSION(:), POINTER :: jE
   !----------------------------------------------------------------------
   

   ! Include file for Message Parsing Interface (mpif77)
   INCLUDE 'mpif.h'
   !------------------------------------------------------------------
   
   CONTAINS
 
 
   SUBROUTINE  init_multiproc
   !------------------------------------------------------------------
   IMPLICIT NONE
   !------------------------------------------------------------------

   ! Initialize Message Parsing Interface
   CALL MP_initialize(iErr)

   ! Total number of machines involved in computation.
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MP_nProc, iErr)

   ! Each process gets its own identifier
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, MP_pRank, iErr)
 
   ! Creating communicators
   IF (MP_pRank == 0) THEN   
     MP_master = .TRUE.
     MP_worker = .FALSE.   
   ELSE   
     MP_master = .FALSE.
     MP_worker = .TRUE.   
   ENDIF


   IF (MP_nProc <= 2) THEN
   
     ! If MP_nProc parallel processing is not allowed
     ! since you have one master and at most one worker,
     ! no speed up is obtained.
     MP_job = .FALSE.
     CALL MP_finalize(iErr)
   
     IF (MP_worker)  STOP
   
   ELSE
 
     MP_job = .TRUE.
     
     ! Number of worker processes
     MP_wProc = MP_nProc - 1
     ALLOCATE (part(MP_wProc))    

     MP_COMM_ALL = MPI_COMM_WORLD
      
     IF (MP_master) THEN
       ! Communicator handling only the master process
       CALL MPI_COMM_SPLIT(MP_COMM_ALL, 1, MP_pRank, MP_COMM_MASTER, iErr)
     ELSE
       ! Communicator handling only worker processes
       CALL MPI_COMM_SPLIT(MP_COMM_ALL, 2, MP_pRank, MP_COMM_WORKERS, iErr)
     ENDIF   
   
   ENDIF 
 
   END SUBROUTINE  init_multiproc

  
  
  
  
   SUBROUTINE  MP_initialize(iErr)
   !---------------------------------------------------
   IMPLICIT NONE
   INTEGER :: iErr
   !---------------------------------------------------
   
   CALL MPI_INIT(iErr)
   
   END SUBROUTINE  MP_initialize



   SUBROUTINE  MP_bcast_path(path)
   !------------------------------------------------------------------- 
   IMPLICIT NONE
   
   CHARACTER(LEN=100), INTENT(IN) :: path
   !------------------------------------------------------------------- 

   CALL  MPI_BCAST(path, 100, MPI_CHARACTER, 0, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_bcast_path





   SUBROUTINE  read_partitions(idf)
   !------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: idf
   
   INTEGER, DIMENSION(:), ALLOCATABLE :: Nj_Ife, Nj_Ifv
   INTEGER, DIMENSION(MP_wProc) :: NjE
      
   INTEGER :: i, j, k
   INTEGER :: partIdx
   INTEGER :: Nidx, Ni, dumN, readN, formatN
   
   CHARACTER(LEN=1), DIMENSION(:), POINTER :: jSpec, cSpec
   !------------------------------------------------------------------
   
   READ(idf,*);  READ(idf,*); READ(idf,*)
   READ(idf,*) Nj_G, Njb_G
   READ(idf,*) Nm_G, Nmb_G
   READ(idf,*) NcFE_G, NcFE_b
   READ(idf,*) NcFV_G
   READ(idf,*) Nb_G         
   READ(idf,*);  READ(idf,*);  READ(idf,*)
   
   ALLOCATE (NjI(MP_wProc))
   NjI = 0
   
   ALLOCATE (NcIfe(MP_wProc), &
             NcIfv(MP_wProc))
  
   NcIfe = 0
   NcIfv = 0	     

   ALLOCATE (Nj_Ife(MP_wProc), &
             Nj_Ifv(MP_wProc))
   
   Nj_Ife = 0
   Nj_Ifv = 0
   
   
   ALLOCATE (jE(MP_wProc))
   
   ALLOCATE (NcEXT_FE(MP_wProc), &
             NcEXT_FV(MP_wProc))
   
   NcEXT_FE = 0
   NcEXT_FV = 0   
   
   DO i = 1, MP_wProc

     READ(idf,*)
     READ(idf,*) partIdx

     IF (partIdx == MP_pRank) THEN

       ! Number Nodes							   
       READ(idf,*) NjP  						   

       ! Number of extended nodes in FINITE ELEMENTS metrics		   
       DO j = 1, MP_wProc						   
         IF (j == partIdx) THEN 					   
           Nj_Ife(j) = 0						   
           CYCLE							   
         ENDIF  							   
         READ(idf,*) Nj_Ife(j)  					   
       ENDDO								   

       ! Number of extended nodes in FINITE VOLUMES metrics		   
       DO j = 1, MP_wProc						   
         IF (j == partIdx) THEN 					   
           Nj_Ifv(j) = 0						   
           CYCLE							   
         ENDIF  							   
         READ(idf,*) Nj_Ifv(j)  					   
       ENDDO								   
        								   
       ! Number of Elements						   
       READ(idf,*) NmP  						   
        								   
       ! Number of node-pairs in FINITE ELEMENTS metrics		   
       READ(idf,*) NcFEP						   
        								   
       DO j = 1, MP_wProc						   
         IF (j == partIdx) THEN 					   
           NcIfe(j) = 0 						   
           CYCLE							   
         ENDIF  							   
         READ(idf,*) NcIfe(j)						   
       ENDDO								  
        								   
       ! Extended node-pairs in FINITE ELEMENTS 			   
       READ(idf,*) NcEXT_FE(i)  					   

       ! Number of node-pairs in FINITE VOLUMES metrics 		   
       READ(idf,*) NcFVP						   
        								   
       DO j = 1, MP_wProc						   
         IF (j == partIdx) THEN 					   
           NcIfv(j) = 0 						   
           CYCLE							   
         ENDIF  							   
         READ(idf,*) NcIfv(j)						   
       ENDDO								   
        								   
       ! Extended node-pairs in FINITE VOLUMES  			   
       READ(idf,*) NcEXT_FV(i)  					   
        								   
       ! Number of boundaries
       READ(idf,*) NbP_FE
       READ(idf,*) NbP_FV

       IF (SUM(Nj_Ife) >= SUM(Nj_Ifv)) THEN
         NbP = NbP_FE
       ELSE
         NbP = NbP_FV
       ENDIF
       

       
       ! Nodes-related allocation
       !IF (fem_metric) THEN
       IF (SUM(Nj_Ife) >= SUM(Nj_Ifv)) THEN
         NjI = Nj_Ife
       ELSE
         NjI = Nj_Ifv
       ENDIF 

       NjTP = NjP + SUM(NjI)
       
       ALLOCATE (jSpec(NjTP))
       ALLOCATE (jG_jP(NjTP), jP_jG(Nj_G))
       jP_jG = 0
       ALLOCATE (jP_jG_B(Njb_G))
       jP_jG_B = 0
     
       ! Read Partition Nodes
       READ(idf,*); READ(idf,*); READ(idf,*)

       DO j = 1, NjP
         READ(idf,*) jG_jP(j), jSpec(j)
         jP_jG(jG_jP(j)) = j
       ENDDO
     
       Nidx = NjP
     
       ! Read extended nodes in FINITE ELEMENTS notation
       !IF (fem_metric) THEN
       IF (SUM(Nj_Ife) >= SUM(Nj_Ifv)) THEN
     
         DO j = 1, MP_wProc
           
	   IF (NjI(j) == 0)  CYCLE
	   
           READ(idf,*)
           DO k = Nidx + 1, Nidx + NjI(j), 1
             READ(idf,*) jG_jP(k), jSpec(k)
             jP_jG(jG_jP(k)) = k	     
           ENDDO
           
           Nidx = Nidx + NjI(j)

         ENDDO

         DO j = 1, MP_wProc
	   IF (Nj_Ifv(j) == 0)  CYCLE	 
           READ(idf,*) 
           DO k = 1, Nj_Ifv(j), 1
             READ(idf,*)    
           ENDDO
         ENDDO

       ! Read extended nodes in FINITE VOLUMES notation
       ELSE

         DO j = 1, MP_wProc
	   IF (Nj_Ife(j) == 0)  CYCLE
           READ(idf,*) 
           DO k = 1, Nj_Ife(j), 1
             READ(idf,*)    
           ENDDO
         ENDDO

         DO j = 1, MP_wProc

	   IF (NjI(j) == 0)  CYCLE
         
           READ(idf,*)
           DO k = Nidx + 1, Nidx + NjI(j), 1
             READ(idf,*) jG_jP(k), jSpec(k)
             jP_jG(jG_jP(k)) = k
           ENDDO
           
           Nidx = Nidx + NjI(j)
           
         ENDDO
     
       ENDIF	
     
       NjP_b = COUNT(jSpec == 'B')

       ! Elements-related allocation     
       ALLOCATE (mG_mP(NmP), mP_mG(0:Nm_G))
       mP_mG = 0
       ALLOCATE (mP_mG_B(0:Nmb_G))
       mP_mG_B = 0

       ! Read elements
       READ(idf,*); READ(idf,*); READ(idf,*)
            
       DO j = 1, NmP, 1
         READ(idf,*) mG_mP(j)
         mP_mG(mG_mP(j)) = j
       ENDDO
     
     
       ! FE node-pairs-related allocation     
       NcTFEP = NcFEP + SUM(NcIfe) + NcEXT_FE(i)
       
       ALLOCATE (cSpec(NcTFEP))
       ALLOCATE (cG_cP_FE(NcTFEP), cP_cG_FE(NcFE_G))
       cP_cG_FE = 0
       ALLOCATE (cP_cG_FE_B(NcFE_b))
       cP_cG_FE_B = 0
            
       ! Read Partition FE Node-Pairs
       READ(idf,*); READ(idf,*); READ(idf,*)

       DO j = 1, NcFEP
         READ(idf,*) cG_cP_FE(j), cSpec(j)     
         cP_cG_FE(cG_cP_FE(j)) = j
       ENDDO
     
       Nidx = NcFEP
     
       ! Read interface node-pair in FINITE ELEMENTS notation
       DO j = 1, MP_wProc

         IF (NcIfe(j) == 0)  CYCLE

         READ(idf,*)
         DO k = Nidx + 1, Nidx + NcIfe(j), 1
           READ(idf,*) cG_cP_FE(k), cSpec(k)
           cP_cG_FE(cG_cP_FE(k)) = k
         ENDDO
         
         Nidx = Nidx + NcIfe(j)
         
       ENDDO
       
       ! Extended node-pairs
       READ(idf,*)       
       DO k = Nidx + 1, Nidx + NcEXT_FE(i), 1
         READ(idf,*) cG_cP_FE(k), cSpec(k)
         cP_cG_FE(cG_cP_FE(k)) = k
       ENDDO
       
       ! Number of boundary node-pairs
       NcFEP_b = COUNT(cSpec == 'B')

     
       ! FV node-pairs-related allocation     
       NcTFVP = NcFVP + SUM(NcIfv) + NcEXT_FV(i)
     
       ALLOCATE (cG_cP_FV(NcTFVP), cP_cG_FV(NcFV_G))
       cP_cG_FV = 0

       ! Read Partition FV Node-Pairs
       READ(idf,*); READ(idf,*); READ(idf,*)

       DO j = 1, NcFVP
         READ(idf,*) cG_cP_FV(j)
         cP_cG_FV(cG_cP_FV(j)) = j
       ENDDO
     
       Nidx = NcFVP
     
       ! Read extended nodes in FINITE VOLUMES notation
       DO j = 1, MP_wProc

         IF (NcIfv(j) == 0)  CYCLE

         READ(idf,*)
         DO k = Nidx + 1, Nidx + NcIfv(j), 1
           READ(idf,*) cG_cP_FV(k)
           cP_cG_FV(cG_cP_FV(k)) = k
         ENDDO
         
         Nidx = Nidx + NcIFV(j)
         
       ENDDO

       ! Extended node-pairs
       READ(idf,*)       
       DO k = Nidx + 1, Nidx + NcEXT_FV(i), 1
         READ(idf,*) cG_cP_FV(k), cSpec(k)
         cP_cG_FV(cG_cP_FV(k)) = k
       ENDDO


       ! Reading boundaries index
       READ(idf,*);  READ(idf,*);  READ(idf,*)
     
       ALLOCATE (bG_bP(NbP), bP_bG(Nb_G))
     
       bP_bG = 0
       
       IF (SUM(Nj_Ife) >= SUM(Nj_Ifv)) THEN       
         
	 READ(idf,*)
         DO k = 1, NbP, 1
           READ(idf,*) bG_bP(k)
           bP_bG(bG_bP(k)) = k
         ENDDO
	 
         READ(idf,*)
         DO k = 1, NbP_FV, 1
           READ(idf,*)
         ENDDO
	 
       ELSE
       
	 READ(idf,*)
         DO k = 1, NbP_FE, 1
           READ(idf,*)
         ENDDO
	 
	 READ(idf,*)
         DO k = 1, NbP, 1
           READ(idf,*) bG_bP(k)
           bP_bG(bG_bP(k)) = k
         ENDDO

       ENDIF

            
       EXIT

     
     ELSE
         
	 formatN = 15 + 4*(MP_wProc-1) + 2
	 
	! Partition Nodes 
	READ(idf,*) dumN
	
	! FE nodes
	DO k = 1, MP_wProc - 1
	  READ(idf,*) readN
	  IF (readN == 0)  formatN = formatN - 1
	  dumN = dumN + readN
	ENDDO 
	     
	! FV nodes
	DO k = 1, MP_wProc - 1
	  READ(idf,*) readN
	  IF (readN == 0)  formatN = formatN - 1
	  dumN = dumN + readN
	ENDDO 

        ! Elements
	READ(idf,*) readN
	dumN = dumN + readN
	
        ! Partition FE node-pairs
	READ(idf,*) readN
	dumN = dumN + readN
	
	! FE node-pair
	DO k = 1, MP_wProc - 1
	  READ(idf,*) readN
	  IF (readN == 0)  formatN = formatN - 1
	  dumN = dumN + readN
	ENDDO 

        ! Extended FE node-pairs
	READ(idf,*) readN
	dumN = dumN + readN

        ! Partition FV node-pairs
	READ(idf,*) readN
	dumN = dumN + readN
	
	! FV node-pair
	DO k = 1, MP_wProc - 1
	  READ(idf,*) readN
	  IF (readN == 0)  formatN = formatN - 1	  
	  dumN = dumN + readN
	ENDDO 

        ! Extended FV node-pairs
	READ(idf,*) readN
	dumN = dumN + readN

        ! Partition boundaries
	READ(idf,*) readN
	dumN = dumN + readN + 1
	READ(idf,*) readN
	dumN = dumN + readN + 1

        ! Dummy reading
	DO k = 1, dumN + formatN
	  READ(idf,*)
	ENDDO

     ENDIF
     
   ENDDO


   ! Defining for each processor the nodes to give
   ! to the other processors (dual of the information
   ! present of file partitions.problem_name)
   Ni = NjP
   DO i = 1, MP_wProc
     IF (i == MP_pRank)  CYCLE
     CALL  MP_send_Isize(MP_pRank, i, NjI(i))
     CALL  MP_send_Inodes(MP_pRank, i, jG_jP(Ni+1:Ni+NjI(i)))
     Ni = Ni + NjI(i)
   ENDDO

   DO i = 1, MP_wProc
     IF (i == MP_pRank)  CYCLE
     CALL  MP_receive_Isize(i, NjE(i))
     ALLOCATE (jE(i) % vec(NjE(i)))
     CALL  MP_receive_Inodes(i, NjE(i), jE(i)%vec)
   ENDDO

      
   END SUBROUTINE  read_partitions  




   SUBROUTINE  MP_send_iSolution(proc, target, iSolution)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER,                      INTENT(IN) :: proc, target
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: iSolution
   !--------------------------------------------------------------------
   
   CALL  MPI_SEND(iSolution, SIZE(iSolution,1)*SIZE(iSolution,2), MPI_REAL8, &
                  target, proc*13, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_send_iSolution



   SUBROUTINE  MP_receive_iSolution(source, iEq, iN, iSolution)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER,                      INTENT(IN)    :: source, iEq, iN
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: iSolution
   
   INTEGER :: MP_STATUS(MPI_STATUS_SIZE)
   !--------------------------------------------------------------------

   CALL  MPI_RECV(iSolution, iEq*iN, MPI_REAL8, &
                  source, source*13, MP_COMM_ALL, MP_STATUS, iErr)
      
   END SUBROUTINE  MP_receive_iSolution





   SUBROUTINE  MP_send_Inodes(proc, target, iNodes)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER,               INTENT(IN) :: proc, target
   INTEGER, DIMENSION(:), INTENT(IN) :: iNodes
   !--------------------------------------------------------------------
   
   CALL  MPI_SEND(iNodes, SIZE(iNodes), MPI_INTEGER, &
                  target, proc*7, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_send_Inodes



   SUBROUTINE  MP_receive_iNodes(source, iN, iNodes)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER,               INTENT(IN)    :: source, iN
   INTEGER, DIMENSION(:), INTENT(INOUT) :: iNodes
   
   INTEGER, DIMENSION(iN) :: loc_iNodes
   
   INTEGER :: MP_STATUS(MPI_STATUS_SIZE), i
   !--------------------------------------------------------------------

   CALL  MPI_RECV(loc_iNodes, iN, MPI_INTEGER, &
                  source, source*7, MP_COMM_ALL, MP_STATUS, iErr)

   DO i = 1, iN
     iNodes(i) = jP_jG(loc_iNodes(i))
   ENDDO
      
   END SUBROUTINE  MP_receive_iNodes





   SUBROUTINE  MP_send_Isize(proc, target, iSize)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: proc, target
   INTEGER, INTENT(IN) :: iSize
   !--------------------------------------------------------------------
   
   CALL  MPI_SEND(iSize, 1, MPI_INTEGER, &
                  target, proc*3, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_send_Isize



   SUBROUTINE  MP_receive_iSize(source, Ne)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER, INTENT(IN)    :: source
   INTEGER, INTENT(INOUT) :: Ne
   
   INTEGER :: MP_STATUS(MPI_STATUS_SIZE)
   !--------------------------------------------------------------------

   CALL  MPI_RECV(Ne, 1, MPI_INTEGER, &
                  source, source*3, MP_COMM_ALL, MP_STATUS, iErr)
      
   END SUBROUTINE  MP_receive_iSize





   SUBROUTINE  MP_send_nodes_division(proc)
   !--------------------------------------------------------------------
   IMPLICIT NONE   
   INTEGER, INTENT(IN) :: proc   
   !--------------------------------------------------------------------

   ! Send number of total nodes 
   CALL MPI_SEND(NjTP, 1, MPI_INTEGER, &
                 0, proc*1, MP_COMM_ALL, iErr)

   ! Send nodes flag
   CALL MPI_SEND(jG_jP(1), NjTP, MPI_INTEGER, &
                 0, proc*2, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_send_nodes_division



   SUBROUTINE  MP_receive_nodes_division(proc)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: proc
   INTEGER :: MP_STATUS(MPI_STATUS_SIZE)
   !--------------------------------------------------------------------

   ! Receive number of nodes of partition 'proc'
   CALL MPI_RECV(part(proc)%Nj_P, 1, MPI_INTEGER, &
                 proc, proc*1, MP_COMM_ALL, MP_STATUS, iErr)

   ALLOCATE (part(proc)%jG_jP(part(proc)%Nj_P))

   ! Receive nodes connectivity
   CALL MPI_RECV(part(proc)%jG_jP(1), part(proc)%Nj_P, MPI_INTEGER, &
                 proc, proc*2, MP_COMM_ALL, MP_STATUS, iErr)
   
   PRINT*, '   Reading partition ', proc
   
   END SUBROUTINE  MP_receive_nodes_division





   SUBROUTINE  MP_send_bou_index(proc, loc_iBou)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: proc
   INTEGER, INTENT(IN) :: loc_iBou
   !--------------------------------------------------------------------
   
   CALL  MPI_SEND(loc_iBou, 1, MPI_INTEGER, &
                  proc, proc*3, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_send_bou_index



   SUBROUTINE  MP_receive_bou_index(proc)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: proc
   INTEGER :: MP_STATUS(MPI_STATUS_SIZE)
   !--------------------------------------------------------------------

   CALL  MPI_RECV(iBou, 1, MPI_INTEGER, &
                  0, proc*3, MP_COMM_ALL, MP_STATUS, iErr)
      
   END SUBROUTINE  MP_receive_bou_index





   SUBROUTINE  MP_send_solution(source, target, ww, ww_oo, b_types, dt)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER,                      INTENT(IN) :: source, target
   REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ww
   REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: ww_oo
   INTEGER,      DIMENSION(:),   INTENT(IN) :: b_types
   
   REAL(KIND=8), DIMENSION(:),   INTENT(IN), OPTIONAL :: dt
   !--------------------------------------------------------------------

   CALL MPI_SEND(SIZE(ww,1), 1, MPI_INTEGER, &
                 target, source*5, MP_COMM_ALL, iErr)

   IF (MP_pRank == 1) THEN
     CALL MPI_SEND(ww_oo, SIZE(ww_oo), MPI_REAL8, &
                   target, source*8, MP_COMM_ALL, iErr)
   ENDIF
   
   CALL MPI_SEND(NjP, 1, MPI_INTEGER, &
                 target, source*6, MP_COMM_ALL, iErr)   
   
   CALL MPI_SEND(ww(:,1:NjP), SIZE(ww,1)*NjP, MPI_REAL8, &
                 target, source*7, MP_COMM_ALL, iErr)


   CALL MPI_SEND(NbP, 1, MPI_INTEGER, &
                 target, source*95, MP_COMM_ALL, iErr)

   CALL MPI_SEND(bG_bP, NbP, MPI_INTEGER, &
                 target, source*97, MP_COMM_ALL, iErr)   


   CALL MPI_SEND(b_types, SIZE(b_types), MPI_INTEGER, &
                 target, source*93, MP_COMM_ALL, iErr)
   
   IF (PRESENT(dt)) THEN
     CALL MPI_SEND(dt, SIZE(dt), MPI_REAL8, &
                   target, source*15, MP_COMM_ALL, iErr)   
   ENDIF
   
   END SUBROUTINE  MP_send_solution



   SUBROUTINE  MP_receive_solution(source, ww, ww_oo, b_types, dt)
   !--------------------------------------------------------------------
   IMPLICIT NONE
   
   INTEGER,                      INTENT(IN)    :: source
   REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: ww
   REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: ww_oo
   INTEGER,      DIMENSION(:),   INTENT(INOUT) :: b_types

   REAL(KIND=8), DIMENSION(:), INTENT(INOUT), OPTIONAL :: dt
      
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: locWW
   REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: locDt
   
   INTEGER, DIMENSION(:), ALLOCATABLE :: loc_b_types, loc_bG_bP
   
   INTEGER :: MP_STATUS(MPI_STATUS_SIZE), locEq, locNj, locNb, p
   !--------------------------------------------------------------------

   CALL MPI_RECV(locEq, 1, MPI_INTEGER, &
                 source, source*5, MP_COMM_ALL, MP_STATUS, iErr)

   IF (source == 1) THEN
     CALL MPI_RECV(ww_oo, locEq, MPI_REAL8, &
                   source, source*8, MP_COMM_ALL, MP_STATUS, iErr)
   ENDIF

   CALL MPI_RECV(locNj, 1, MPI_INTEGER, &
                 source, source*6, MP_COMM_ALL, MP_STATUS, iErr)

   ALLOCATE (locWW(locEq,locNj), locDt(locNj))

   CALL MPI_RECV(locWW, locEq*locNj, MPI_REAL8, &
                 source, source*7, MP_COMM_ALL, MP_STATUS, iErr)

   DO p = 1, locNj, 1
     ww(:,part(source)%jG_jP(p)) = locWW(:,p)
   ENDDO

   CALL MPI_RECV(locNb, 1, MPI_INTEGER, &
                 source, source*95, MP_COMM_ALL, MP_STATUS, iErr)

   ALLOCATE (loc_b_types(locNb), loc_bG_bP(locNb))
   
   
   CALL MPI_RECV(loc_bG_bP, locNb, MPI_INTEGER, &
                 source, source*97, MP_COMM_ALL, MP_STATUS, iErr)   

   CALL MPI_RECV(loc_b_types, locNb, MPI_INTEGER, &
                 source, source*93, MP_COMM_ALL, MP_STATUS, iErr)

   DO p = 1, SIZE(loc_b_types), 1
     b_types(loc_bG_bP(p)) = loc_b_types(p)
   ENDDO
   
   IF (PRESENT(dt)) THEN
   
     CALL  MPI_RECV(locDt, locNj, MPI_REAL8, &
                    source, source*15, MP_COMM_ALL, MP_STATUS, iErr)
		    
     DO p = 1, locNj
       dt(part(source)%jG_jP(p)) = locDt(p)
     ENDDO
     
   ENDIF
   
   
   DEALLOCATE (locWW, locDt, loc_b_types, loc_bG_bP)
      
   END SUBROUTINE  MP_receive_solution





   SUBROUTINE  MP_reduce_resid(ww_resid_IN, ww_resid_OUT)
   !------------------------------------------------------------------- 
   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: ww_resid_IN
   REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: ww_resid_OUT
   !------------------------------------------------------------------- 

   CALL  MPI_REDUCE(ww_resid_IN, ww_resid_OUT, SIZE(ww_resid_IN,1), &
                    MPI_REAL8, MPI_SUM, 0, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_reduce_resid





   SUBROUTINE  MP_bcast_resid(resid)
   !------------------------------------------------------------------- 
   IMPLICIT NONE
   
   REAL(KIND=8), INTENT(IN) :: resid
   !------------------------------------------------------------------- 

   CALL  MPI_BCAST(resid, 1, MPI_REAL8, 0, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_bcast_resid




   SUBROUTINE  MP_bcast_dt_min(dtM)
   !------------------------------------------------------------------- 
   IMPLICIT NONE
   
   REAL(KIND=8), INTENT(IN) :: dtM
   !------------------------------------------------------------------- 

   CALL  MPI_BCAST(dtM, 1, MPI_REAL8, 0, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_bcast_dt_min





   SUBROUTINE  MP_send_MaxMin_dt(source, target, locDt)
   !------------------------------------------------------------------- 
   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: source, target
   REAL(KIND=8), DIMENSION(NjP) :: locDt
   !-------------------------------------------------------------------

   CALL MPI_SEND(MINVAL(locDt), 1, MPI_REAL8, &
                 target, source*91, MP_COMM_ALL, iErr)

   CALL MPI_SEND(MAXVAL(locDt), 1, MPI_REAL8, &
                 target, source*112, MP_COMM_ALL, iErr)

   END SUBROUTINE  MP_send_MaxMin_dt
   
   
   
   
   
   SUBROUTINE  MP_receive_MaxMin_dt(source, minDt, maxDt)
   !------------------------------------------------------------------- 
   IMPLICIT NONE
   
   INTEGER,      INTENT(IN) :: source
   REAL(KIND=8), INTENT(INOUT) :: minDt, maxDt
      
   INTEGER :: MP_STATUS(MPI_STATUS_SIZE)   
   !------------------------------------------------------------------- 
   
   CALL MPI_RECV(loc_minDt, 1, MPI_REAL8, &
                 source, source*91, MP_COMM_ALL, MP_STATUS, iErr)
   
   IF (minDt > loc_minDt)  minDt = loc_minDt   
   
   CALL MPI_RECV(loc_maxDt, 1, MPI_REAL8, &
                 source, source*112, MP_COMM_ALL, MP_STATUS, iErr)
   
   IF (maxDt < loc_maxDt)  maxDt = loc_maxDt
      
   END SUBROUTINE  MP_receive_MaxMin_dt





   SUBROUTINE  MP_gather_all_minDt(mDt, mDtV)
   !------------------------------------------------------------------- 
   IMPLICIT NONE
   
   REAL(KIND=8),                      INTENT(IN)  :: mDt
   REAL(KIND=8), DIMENSION(MP_wProc), INTENT(OUT) :: mDtV
   !------------------------------------------------------------------- 
   
   CALL MPI_ALLGATHER(mDt, 1, MPI_REAL8, mDtV, 1, MPI_REAL8, &
                      MP_COMM_WORKERS, iErr)
   
   END SUBROUTINE  MP_gather_all_minDt





   SUBROUTINE  MP_finalize(iErr)
   !---------------------------------------------------
   IMPLICIT NONE
   INTEGER :: iErr
   !---------------------------------------------------
   
   CALL MPI_FINALIZE(iErr)
   
   END SUBROUTINE  MP_finalize


   END MODULE mp_interface
