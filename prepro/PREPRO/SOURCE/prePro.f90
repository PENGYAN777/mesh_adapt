!============================================================ 
!
!     Program: PreProcessor
!
! Description: Computes the node-pair metric quantities
!              for both the finite element and the finite
!              volume methods starting from a grid in UMesh 
!              format.
!              Works in 1, 2 and 3 dimensions.  
!
!      Author: Alberto Guardone, Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!                      marco.fossati@polimi.it
!
!============================================================ 

  PROGRAM  geometric_kernel

! This program build the node-pair/node connectivity of the 
! mesh, computes the metric quantities and if required computes
! the mesh partition to allow for parallel processing of the solver
! 
! 
! INPUT FILES:
!   nodes.gd_name     Nodes coordinates and boundary connectivity
!   grid.gd_name      Element/Nodes connectivity
!
!
! OUTPUT FILES:
!   node_pair.gd_name         Node-pair / Nodes connectivity
!   metrics.gd_name           Metric quantities 
!   metrics_axial.gd_name     Axisymmetrical metric quantities
!
! La numerazione dei nodi di ciascun elemento e' sempre
! assegnata percorrendo il contorno dell'elemento in
! senso antiorario (2D)
!-----------------------------------------------------------------------

  use messages, only: write_prepro_logo
  use nodes 
  use mesh_structure
  use node_pair_structure
  use np_topology_gen
  use metric_coefficients
  use np_metric_gen
  use mesh_partition
  use faces_gen

  implicit none

  real*8 :: maxFs, t1, t2
  integer :: nCPU, support
  integer :: d_END, p_END, nBars

  character (len=64) :: directory = './'
  character (len=64) :: grid_name
  logical :: parallel, saveE, axisym
  logical :: filExist
  !-----------------------------------------------------------------------

  call cpu_time(t1)
  call write_prepro_logo(6)

  inquire (file='prepro.cfg', EXIST=filExist)
  if ( filExist ) then

    open ( unit=11, file='prepro.cfg' )
      read (11,*); read (11,*); read (11,*)
      read (11,*) grid_name
      read (11,*) axisym
      read (11,*) parallel, nCPU
      read (11,*) saveE
      read (11,*) nBars
      read (11,*) maxFs
    close (11)

  else

    write (*,*) ''
    write (*,*) 'ERROR. Missing config file prepro.cfg'
    write (*,*) ''
    stop

  end if

  d_END = last_c_leng (64, directory)
  p_END = last_c_leng (64, grid_name)

  ! Lettura dei nodi and of the elements of the grid
  WRITE(*,*)
  OPEN (1,file=directory(1:d_END)//'/nodes.'//grid_name(1:p_END))
  CALL read_nodes(1)
  CLOSE(1)
  OPEN (1,file=directory(1:d_END)//'/grid.'//grid_name(1:p_END))
  CALL read_mesh(1)
  CLOSE(1)

  ! Costruzione delle coppie di nodi
  WRITE(*,*)
  WRITE(*,*) '   Generating node-pairs...'
  CALL node_pair_str_gen(grid_name)

  ! Costruzione delle informazioni metriche
  WRITE(*,*) ''
  WRITE(*,*) '   Generating metric quantities...'
  CALL metrics_gen

  WRITE(*,*) ''
  WRITE(*,*) '   Checking metric consistency...'
  CALL metric_check
  WRITE(*,*) '   done'
    
  IF ( k_d .eq. 2 .and. axisym ) THEN 
    WRITE(*,*)
    WRITE(*,*) '   Checking AXISYMMETRIC metric consistency...'
    WRITE(*,*)
    CALL metric_check_axi
    WRITE(*,*) '   done';  WRITE(*,*)
  ENDIF
 
  ! Computes the cosines matrices at each interface
  ! relative to the integrated normals for the finite
  ! volumes notation (normals are eta_fv and xi_bp)
  WRITE(*,*)
  WRITE(*,*) '   Generating rotation matrices...'
  CALL rotation_matrices_gen
  WRITE(*,*) '   done'
  WRITE(*,*)
  WRITE(*,*) '   Saving...'

  ! Scrittura delle informazioni dei node-pairs su un file
  OPEN (1,file=directory(1:d_END)//'/node_pair.'//grid_name(1:p_END))
  CALL save_node_pair(1, grid_name,p_END)
  CLOSE(1)
  WRITE(*,*) '   node-pairs...'

  ! Scrittura delle informazioni metriche su un file
  OPEN (1,file=directory(1:d_end)//'/metrics.'//grid_name(1:p_end))
  CALL save_metric(1, grid_name,p_END)
  CLOSE(1)
  WRITE(*,*) '   metric quantities...'

  ! Scrittura delle informazioni metriche in simmetria assiale
  IF ( k_d .eq. 2 .and. axisym ) THEN
    OPEN (1,file=directory(1:d_end)//'/metrics_axial.'//grid_name(1:p_end))
    CALL save_metric_axisymmetric(1, grid_name,p_END)
    CLOSE(1)
    WRITE(*,*) '   axisymmetric metric quantities...'
  ENDIF  

  WRITE(*,*) '   done'
  WRITE(*,*)

  ! Partitioning mesh usign METIS library
  IF ( parallel ) THEN
    WRITE(*,'(a23,1x,i2,1x,a8)') '   Partitioning mesh in', nCPU, 'parts...'
    support = 0
    CALL partMesh ( grid_name, nCPU, support )
    PRINT*, '   done'
    PRINT*, ''  
  ENDIF

  ! Evaluating the quality of the grid elements
  IF ( k_d .gt. 1 .and. saveE) THEN
    WRITE(*,*) '   Evaluation of the grid quality...'      
    CALL grid_quality ( grid_name, nBars, maxFs )
    PRINT*, ''
  END IF  

  call cpu_time(t2)
  WRITE(*,'(4x,a13,1x,e11.4,1x,a1)') 'Elapsed time ', t2-t1, 's'
  WRITE(*,*)

  CONTAINS

  
  SUBROUTINE  grid_quality ( gName, nBars, maxFs )
  !---------------------------------------------------------------
  USE  dynamic_vector
  USE  node_pair_structure
  USE  mesh_structure
  USE  np_topology_gen
  USE  nodes
    
  IMPLICIT NONE
  
  CHARACTER(LEN=64), INTENT(IN) :: gName

  TYPE(D_I_V), DIMENSION(:), ALLOCATABLE  :: j_c_DIV, c_j
  
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: points, bubbleV, vectorCopy
  REAL(KIND=8), DIMENSION(:),	ALLOCATABLE :: area, aRatio, bubbleA, angleCopy
  REAL(KIND=8), DIMENSION(:),	ALLOCATABLE :: F_s
  REAL(KIND=8), DIMENSION(2,2,size(rr,2)) :: HM
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Inp, j_c
  INTEGER, DIMENSION(:),   ALLOCATABLE :: cIn, bubbleN, nodesCopy, ordered_idx
  
  REAL(KIND=8), DIMENSION(2) :: rDist 
  CHARACTER(LEN=1) :: info, saveE
  REAL(KIND=8) :: qEstimate, PIGR = 3.141592653589793
  REAL(KIND=8) :: sqLm, sqLn, bMN, gM, qFs, qFcng, qFo, maxFs
  INTEGER :: idf=11, i, j, j_, k, k_, c, m, m_, nP, c_B_, c_B, nBars
  LOGICAL, DIMENSION(:), ALLOCATABLE :: done
  LOGICAL :: cont = .TRUE., loop, fexist = .FALSE.
  LOGICAL :: filExist, riemann
  !---------------------------------------------------------------  
  
  INQUIRE (file='np-interface.'//TRIM(gName), EXIST=fexist)
  IF (fexist) THEN
  
    ! Structured/unstructured interface quality 
    OPEN(UNIT=idf, FILE='np-interface.'//TRIM(gName))
  
      DO WHILE (cont)
      
    	READ(idf,*) info
    	IF (info == 'N') THEN
    	  READ(idf,*) nP
    	  cont = .FALSE.
    	  EXIT    
    	ENDIF
      
      ENDDO

      ALLOCATE (Inp(2,nP))
      REWIND(idf)
      
      DO c = 1, nP
    	READ(idf,*) Inp(:,c)
      ENDDO

    CLOSE(idf)

    ALLOCATE (cIn(nP), aRatio(nP))

    DO k = 1, SIZE(Inp,2) 

      DO c = 1, SIZE(j_c_fv,2)
    	
    	i = j_c_fv(1,c)
    	j = j_c_fv(2,c)
      
    	IF ((i == Inp(1,k) .AND. j == Inp(2,k))  &
    			   .OR. 		 &
    	    (i == Inp(2,k) .AND. j == Inp(1,k))) THEN

    	  cIn(k) = c

    	ENDIF
    		  
      ENDDO
    ENDDO

    DO i = 1, nP
  
      c = cIn(i)
      ALLOCATE (area(SIZE(m_c(c)%vec)))
      area = 0.d0
  
      DO m_ = 1, SIZE(m_c(c)%vec)
      
    	m = m_c(c)%vec(m_)
      
    	ALLOCATE (points(2, SIZE(j_m_d(m)%vec) + 1))

    	! Loop over nodes of element m
    	DO k_ = 1, SIZE(j_m_d(m)%vec)
    	  k = j_m_d(m)%vec(k_)  	  
    	  points(:,k_) = rr(:,k)
    	ENDDO
    	
    	points(:,SIZE(points,2)) = points(:,1)
    	
    	DO k = 1, SIZE(points,2) - 1
    	  area(m_) = area(m_) + (points(2,k+1)+points(2,k))* &
	             (points(1,k+1)-points(1,k))/2
    	ENDDO
    	area(m_) = ABS(area(m_))
      
    	DEALLOCATE (points)
      
      ENDDO

      aRatio(i) = MAXVAL(area) / MINVAL(area)
  
      DEALLOCATE (area)
  
    ENDDO

    qEstimate = SQRT(SUM(aRatio**2)/nP)  
    WRITE(*,*) '   - Interface quality factor	 -->', qEstimate

  ENDIF

   
  ! Global surface mesh quality --------------------------------
  ! Computes the node to node-pair connectivity 
  ! matrix c_j using DIV algorithms
  ALLOCATE (j_c(SIZE(j_c_fv,1),Nc_fv))
  j_c = j_c_fv
  j_c(3:4*k_d,:) = 0
  
  ALLOCATE (j_c_DIV(SIZE(j_c,2)))
  j_c_DIV = convert_matrix_to_DIV(j_c)
  
  ALLOCATE (c_j(size_DIV(j_c_DIV,3)))
  c_j = invert_DIV(j_c_DIV)
  DEALLOCATE (j_c_DIV)  

  ALLOCATE( m_j_d(size_DIV(j_m_d,3)) )
  m_j_d = invert_DIV(j_m_d)

  
  ALLOCATE ( F_s(Np_d-Np_b/2) )
  F_s = 0.d0
  
  INQUIRE (file='hessian.'//trim(gName), EXIST=filExist)
  IF ( filExist ) THEN
    OPEN (UNIT=11, file='hessian.'//trim(gName))
    DO j = 1, size(rr,2)
      read (11,*) HM(1,:,j), HM(2,:,j)
    END DO
    CLOSE (11)
    riemann = .TRUE.
  ELSE
    riemann = .FALSE.
  END IF

  j = 0
  DO j_ = 1, Np_d

    IF (ANY(jd_jb == j_)) CYCLE
  
    j = j + 1
  
    ! Extract nodes of the bubble for node j and
    ! defines vectors P(:,m) - P(:,j) m = (1,...,size(bubble))
    ALLOCATE (bubbleN(     SIZE(c_j(j_)%vec)), &
              bubbleV(k_d, SIZE(c_j(j_)%vec)), &
              bubbleA(     SIZE(c_j(j_)%vec)))
  
    DO c_B_ = 1, SIZE(c_j(j_)%vec)
    
       c_B = c_j(j_)%vec(c_B_)            

       IF (j_c_fv(1,c_B) /= j_) THEN

         bubbleN(c_B_)   = j_c_fv(1,c_B)

         IF ( riemann ) then
           call riemann_distance ( rr, j_, j_c_fv(1,c_B), HM(:,:,j_),  &
             HM(:,:, j_c_fv(1,c_B)), rDist )
           bubbleV(:,c_B_) = rDist
         ELSE
           bubbleV(:,c_B_) = rr(:,j_c_fv(1,c_B)) - rr(:,j_)
         END IF

       ELSE 

         bubbleN(c_B_)   = j_c_fv(2,c_B)

         IF ( riemann ) THEN
           call riemann_distance ( rr, j_, j_c_fv(2,c_B),  HM(:,:,j_),  &
             HM(:,:, j_c_fv(2,c_B)), rDist )
           bubbleV(:,c_B_) = rDist
         ELSE
           bubbleV(:,c_B_) = rr(:,j_c_fv(2,c_B)) - rr(:,j_)
         END IF

       ENDIF  

       IF (bubbleV(1,c_B_) == 0.d0  .and.  bubbleV(2,C_B_) > 0.d0) THEN
         bubbleA(c_B_) = 90.d0
       ELSEIF (bubbleV(1,c_B_) == 0.d0  .and.  bubbleV(2,C_B_) < 0.d0) THEN
         bubbleA(c_B_) = 270.d0
       ELSE
         bubbleA(c_B_) = ATAN2(bubbleV(2,c_B_), bubbleV(1,c_B_)) * 180.0/PIGR
         IF (bubbleA(c_B_) < 0) bubbleA(c_B_) = 360.0 + bubbleA(c_B_)
       ENDIF

    ENDDO

    ! Order the nodes of the bubble counter-clockwise
    ALLOCATE (ordered_idx(SIZE(bubbleN)))
    ALLOCATE (angleCopy(SIZE(bubbleA)))
    ALLOCATE (nodesCopy(SIZE(bubbleN)))
    ALLOCATE (vectorCopy(SIZE(bubbleV,1), SIZE(bubbleV,2)))

    ordered_idx = 0
    angleCopy = - bubbleA
    
    i = 1
    k = 1
    loop = .TRUE.
    DO WHILE (loop)
    
      IF (angleCopy(k) == MAXVAL(angleCopy)) THEN
        ordered_idx(i) = k
        angleCopy(k) = -HUGE(angleCopy(k))
        i = i + 1
        k = 0
        IF (ALL(ordered_idx /= 0)) loop = .FALSE.
        
      ENDIF

      k = k + 1
      
    ENDDO


    DO k = 1, SIZE(bubbleN)     
      nodesCopy(k)    = bubbleN(ordered_idx(k))
      vectorCopy(:,k) = bubbleV(:,ordered_idx(k))
    ENDDO
    
    DEALLOCATE (bubbleN, bubbleV)
    ALLOCATE (bubbleN(  SIZE(nodesCopy)+1),  &
              bubbleV(2,SIZE(vectorCopy)+1))

    bubbleN(1:SIZE(nodesCopy)) = nodesCopy
    bubbleN(SIZE(bubbleN)) = bubbleN(1)
    
    bubbleV(:,1:SIZE(nodesCopy)) = vectorCopy
    bubbleV(:,SIZE(bubbleN)) = bubbleV(:,1)

    ! Compute quality functionals (IJNME vol 48, pp 401-420)
    F_s(j) = 0.d0
    
    IF (ALLOCATED(done)) DEALLOCATE (done)
    ALLOCATE (done(SIZE(m_j_d(j_)%vec)))
    done = .FALSE.

    k = 0
    k_ = 0 
    DO WHILE ( .NOT. ALL(done) )

      k = k + 1
      k_ = k_ + 1

      DO m = 1, SIZE(j_m_d)

        IF (ANY(j_m_d(m)%vec == bubbleN(k)) .AND. &
            ANY(j_m_d(m)%vec == bubbleN(k+1))) THEN

          sqLm = SUM(bubbleV(:,k)   * bubbleV(:,k))
          sqLn = SUM(bubbleV(:,k+1) * bubbleV(:,k+1))
          bMN  = SUM(bubbleV(:,k+1) * bubbleV(:,k))

          gM = sqLm * sqLn - bMN**2
          F_s(j) = F_s(j) + 0.5* (sqLm + sqLn) / SQRT(gM)

          done(k_) = .TRUE.
          EXIT
        ENDIF

      ENDDO

      k_ = COUNT(done)
      
    ENDDO
    
    F_s(j) = F_s(j) / (SIZE(bubbleN) - 1)
    
    DEALLOCATE ( bubbleN )
    DEALLOCATE ( bubbleA )
    DEALLOCATE ( bubbleV )
    DEALLOCATE ( ordered_idx )
    DEALLOCATE ( angleCopy )
    DEALLOCATE ( nodesCopy )
    DEALLOCATE ( vectorCopy )

  ENDDO

  qFs = SQRT(SUM(F_s**2) / SIZE(F_s))
  write (*,*) ''
  write (*,*) '   - Indicator range [', MINVAL(F_s), ' - ', MAXVAL(F_s), ']'
  write (*,*) '   - Indicator L2-norm = ', qFs
  write (*,*) ''

  OPEN (UNIT=15, FILE='quality-extrema.dat')
    WRITE(15,*) MINVAL(F_s), MAXVAL(F_s)
  CLOSE(15)
  CALL plot_quality_indicators( gName, F_s, rr, j_m_d, nBars, maxFs ) 

  return
  END SUBROUTINE  grid_quality




  SUBROUTINE riemann_distance ( rr, i, m, HM1, HM2, rDist )
  !--------------------------------------------------------------
  implicit none

  real*8, dimension(:,:), intent(in) :: rr
  integer,                intent(in) :: i, m
  real*8, dimension(2,2), intent(in) :: HM1, HM2
  real*8, dimension(2),   intent(inout) :: rDist

  real*8, dimension(2) :: gamma_p, vers
  real*8 :: l1, l2
  real*8 :: gamma_pT__M__gamma_p
  !--------------------------------------------------------------

  vers = (rr(:,m) - rr(:,i)) / DSQRT( SUM( (rr(:,m) - rr(:,i))**2 ) )

  gamma_p(1) = rr(1,m) - rr(1,i)
  gamma_p(2) = rr(2,m) - rr(2,i)

  gamma_pT__M__gamma_p = SUM( gamma_p*MATMUL(HM1(:,:), gamma_p) ) + 1.e-7
  l1 = DSQRT( gamma_pT__M__gamma_p )

  gamma_pT__M__gamma_p = SUM( gamma_p*MATMUL(HM2(:,:), gamma_p) ) + 1.e-7
  l2 = DSQRT( gamma_pT__M__gamma_p )

  rDist = (2.d0/3.d0) * DABS( (l1**2 + l1*l2 + l2**2) / (l1+l2) ) * vers

  RETURN
  END SUBROUTINE riemann_distance  



  SUBROUTINE plot_quality_indicators( gName, F_s, rr, j_m_d, ncol, MAXVALF_s ) 
  !--------------------------------------------------------------
  IMPLICIT NONE
  
  CHARACTER(LEN=64),            INTENT(IN) :: gName
  REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: F_s
  REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
  TYPE(D_I_V),  DIMENSION(:),   INTENT(IN) :: j_m_d
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: F_Ranges
  INTEGER,      DIMENSION(:), ALLOCATABLE :: nPoints
  REAL(KIND=8) :: dF, MINVALF_s, MAXVALF_s  
  
  INTEGER :: i, j, nCol
  !--------------------------------------------------------------
  
!   OPEN(UNIT=11, FILE='quality-field.'//TRIM(gName)//'.plt')
!   
!   WRITE(11,*) 'TITLE = "Quality indicators"'
!   WRITE(11,*) 'VARIABLES = "X" "Y" "Fsmooth" "Fcng" "Fodd"'
!   WRITE(11,*) 'ZONE T = "DOMAIN", N=', SIZE(rr,2), 'E=', SIZE(j_m_d), 
!     'DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
!   
!   DO i = 1, SIZE(F_s)
!     WRITE(11,*) rr(:,i), F_s(i), F_cng(i), F_o(i)
!   ENDDO
!   
!   DO i = 1, SIZE(j_m_d)
!     IF (SIZE(j_m_d(i)%vec) == 3) WRITE(11,*) j_m_d(i)%vec, j_m_d(i)%vec(1)
!     IF (SIZE(j_m_d(i)%vec) == 4) WRITE(11,*) j_m_d(i)%vec
!   ENDDO
!   
!   CLOSE(11)

!  OPEN(UNIT=11, FILE='quality-extrema.dat')
!    READ(11,*) MINVALF_s,     MAXVALF_s
!  CLOSE(11)

  !write (*,'(1x,a31)',advance='no') '   - Number of histogram bars? '
  !read (*,*) ncol
  !write (*,'(1x,a37)',advance='no') '   - Maximum value of the indicator? '
  !read(*,*) MAXVALF_s
  MINVALF_s = 0.0

  ALLOCATE ( F_Ranges(nCol+1) )
  ALLOCATE ( nPoints(nCol) )
  
  ! Hystogram data for F_smooth
  F_Ranges(1) = MINVALF_s
  F_Ranges(nCol+1) = MAXVALF_s

  if ( abs(MAXVALF_s - MINVALF_s) .lt. 1.d-5 ) then
    dF = 0.d0
  else
    dF = (MAXVALF_s - MINVALF_s) / nCol
  end if
  
  DO i = 2, ncol
    F_Ranges(i) = F_Ranges(i-1) + dF
  ENDDO
    
  nPoints = 0
  DO j = 1, nCol
    nPoints(j) = COUNT ( F_s >= F_ranges(j) .AND. F_s < F_ranges(j+1) )
  ENDDO

  OPEN(UNIT=11, FILE='histogram.'//TRIM(gName)//'.plt')
  
  WRITE(11,*) 'TITLE = "Quality indicators"'
  WRITE(11,*) 'VARIABLES = "Fsmooth" "Np_s"'
  WRITE(11,*) 'ZONE T = "DOMAIN", I=', nCol
  
  DO i = 1, nCol
    WRITE(11,*) (F_ranges(i+1)+F_ranges(i))/2.0, real(nPoints(i)*100/size(F_s))
  ENDDO
  
  CLOSE(11)
  

  
  END SUBROUTINE plot_quality_indicators
  

  

  
  FUNCTION  last_c_leng (len_str, string)   RESULT(leng)
  !--------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: len_str
  CHARACTER (LEN=len_str), INTENT(IN) :: string
  INTEGER :: leng

  INTEGER :: i
  !--------------------------------------------------------------

  leng = len_str

  DO i=1,len_str
    IF (string(i:i) == ' ') THEN
      leng = i-1
      EXIT
    ENDIF
  ENDDO

  END FUNCTION  last_c_leng


  END PROGRAM  geometric_kernel
