PROGRAM test_distances

   USE distances
   USE nodes
   USE mesh_structure
   USE element_topology
   USE dynamic_vector
   USE strings

   IMPLICIT NONE

   ! if is_solid_boundary(b)=TRUE, boundary b is a wall
   LOGICAL,      DIMENSION(:),  ALLOCATABLE  ::  is_solid_boundary
   ! wall_distance(i) is the distance of node i from the wall
   REAL(KIND=8), DIMENSION(:),  ALLOCATABLE  ::  wall_distance
   ! Interpolation nodes and interpolation weights for node i
   INTEGER,      DIMENSION(:,:), ALLOCATABLE  ::  interpolation_nodes
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  ::  interpolation_weights
   ! ib = nearest_node(i): ib is the boundary node nearest to node i                                             
   INTEGER,      DIMENSION(:),  ALLOCATABLE  ::  nearest_node
   ! local interpolation weights
   REAL(KIND=8), DIMENSION(:),  ALLOCATABLE  ::  weights
   ! Boundary node -> domain node connectivity for single boundary nodes
   INTEGER,      DIMENSION(:),  POINTER      ::  jd_jb_single
   REAL(KIND=8) :: d ! Distance
   INTEGER  ::  N_boundaries, & ! Number of boundaries
                idx             ! dummy
   INTEGER  ::  i, i1, i2, i3, &  ! Domain node indices
                ib,            &  ! Boundary node index
                m,             &  ! Boundary element index 
                m_                ! Bubble element local index
   ! Grid reading             
   CHARACTER(len=64)  ::  grid_name  
   INTEGER  ::  name_length          


   !======================================================
   ! Reads parameter
   OPEN(1, file='distances.param')
   READ(1,*) grid_name
   name_length = last_c_leng (64, grid_name)
   
   READ(1,*) N_boundaries
   ALLOCATE(is_solid_boundary(N_boundaries))
   DO i = 1, N_boundaries
      READ(1,*) idx, is_solid_boundary(i)
   ENDDO
   CLOSE(1)
   
   ! Reads nodes
   OPEN(2,file='nodes.'//grid_name(1:name_length))
   CALL  read_nodes(2)
   CLOSE(2)

   ! Reads grid
   OPEN(2,file='grid.'//grid_name(1:name_length))
   CALL  read_mesh(2)
   CLOSE(2)
   !======================================================


   !======================================================
   ! Delete double boundary nodes
   WRITE(*,*) '   Deleting double nodes'
   CALL del_double_nodes_boundary (jd_jb, bound_p,&
                                   j_m_b, bound_m,&
                                   jd_jb_single,  &
                                   .TRUE.         )
   
   DEALLOCATE(jd_jb)
   ALLOCATE(jd_jb(SIZE(jd_jb_single)))
   jd_jb = jd_jb_single
   Np_b = SIZE(jd_jb)
   !======================================================
  
     
   !======================================================
   ! Compute the nearest boundary node nearest_node(i)
   ! with respect to domain node i
   WRITE(*,*) '   Computing nearest boundary nodes on ', Np_b,' nodes'
   ALLOCATE(wall_distance(Np_d), nearest_node(Np_d))
   wall_distance = HUGE(d)

   DO i = 1, SIZE(rr,2)
      DO ib = 1, SIZE(jd_jb)
   
         IF (.NOT. is_solid_boundary(bound_p(ib))) CYCLE
   
         i1 = jd_jb(ib)

         d = distance_from_point( rr(:,i), rr(:,i1) )
         
         IF (d < wall_distance(i)) THEN
            wall_distance(i) = d
             nearest_node(i) = ib
         ENDIF
         
      ENDDO
   ENDDO
   !======================================================


   !======================================================
   ! Creates the boundary node --> boundary elements
   ! connectivity m_j_b and allocate interpolation nodes
   ! and weights
   WRITE(*,*) '   Computing distances'
   ALLOCATE(m_j_b(SIZE(jd_jb)))
   m_j_b = invert_DIV(j_m_b)
   
   ALLOCATE( weights(k_d) )
   ALLOCATE(   interpolation_nodes(k_d, Np_d), &
             interpolation_weights(k_d, Np_d)  )
   !======================================================
   
   !======================================================
   ! Computes the wall distance and the interpolation 
   ! nodes and weights for domain node i
   wall_distance = HUGE(d)
   DO i = 1, SIZE(rr,2)
      
      ! Retrive the nearest boundary node ib
      ib = nearest_node(i)
      
      ! For all boundary elements belonging to the bubble 
      ! of ib, computes the distance and the interpolation 
      ! weights
      DO m_ = 1, SIZE(m_j_b(ib)%vec); m = m_j_b(ib)%vec(m_)
   
   
         SELECT CASE (ele_type_b(m))
         
            !----------------------------------------------------------
            CASE(ELE_TYPE_SEGMENT)
            !---------------------   
               i1 = jd_jb(j_m_b(m)%vec(1))
               i2 = jd_jb(j_m_b(m)%vec(2))
               
               d = distance_from_segment( rr(:,i),            &
                                          rr(:,i1), rr(:,i2), &
                                          weights             )
               
               IF (d < wall_distance(i)) THEN
                    wall_distance(i) = d
                  interpolation_weights(:,i) = ABS(weights)
                  interpolation_nodes  (:,i) = jd_jb(j_m_b(m)%vec)
               ENDIF
            !----------------------------------------------------------
            
            !----------------------------------------------------------
            CASE(ELE_TYPE_TRIANGLE)
            !---------------------   
               i1 = jd_jb(j_m_b(m)%vec(1))
               i2 = jd_jb(j_m_b(m)%vec(2))
               i3 = jd_jb(j_m_b(m)%vec(3))
               
               d = distance_from_triangle( rr(:,i),                      & 
                                           rr(:,i1), rr(:,i2), rr(:,i3), &
                                           weights                       )
               
               IF (d < wall_distance(i)) THEN
                    wall_distance(i) = d
                  interpolation_weights(:,i) = ABS(weights) 
                  interpolation_nodes  (:,i) = jd_jb(j_m_b(m)%vec)
               ENDIF
            !----------------------------------------------------------
   
            !----------------------------------------------------------
            CASE(ELE_TYPE_QUADRILATER)
            !---------------------   
               i1 = jd_jb(j_m_b(m)%vec(1))
               i2 = jd_jb(j_m_b(m)%vec(2))
               i3 = jd_jb(j_m_b(m)%vec(3))
               
               d = distance_from_triangle( rr(:,i),                      & 
                                           rr(:,i1), rr(:,i2), rr(:,i3), &
                                           weights                       )
               
               IF (d < wall_distance(i)) THEN
                    wall_distance(i) = d
                  interpolation_weights(:,i) = ABS(weights) 
                  interpolation_nodes  (:,i) = jd_jb(j_m_b(m)%vec((/1,2,3/)))
               ENDIF
               
               i1 = jd_jb(j_m_b(m)%vec(3))
               i2 = jd_jb(j_m_b(m)%vec(4))
               i3 = jd_jb(j_m_b(m)%vec(1))
               
               d = distance_from_triangle( rr(:,i),                      & 
                                           rr(:,i1), rr(:,i2), rr(:,i3), &
                                           weights                       )
               
               IF (d < wall_distance(i)) THEN
                    wall_distance(i) = d
                  interpolation_weights(:,i) = ABS(weights) 
                  interpolation_nodes  (:,i) = jd_jb(j_m_b(m)%vec((/3,4,1/)))
               ENDIF
               
               
            !----------------------------------------------------------
   
        END SELECT
   
      ENDDO
   
   ENDDO
   !======================================================


    DO i = 1, SIZE(rr,2)
      IF (SUM(interpolation_weights(:,i)) < (1.d0 - 1d-9)) THEN
         WRITE(*,*) 'Error: the sum of the interpolation weights'
         WRITE(*,*) 'is less than 1.d0 - 1d-9'
         WRITE(*,*) ' at node  ', i
         WRITE(*,*) ' weights: ', interpolation_weights(:,i)
         WRITE(*,*) 'STOP'
       STOP
      ENDIF
   ENDDO


   OPEN(2,file='distances.'//grid_name(1:name_length))
   DO i = 1, SIZE(rr,2)
      WRITE(2,*) i, &
                 wall_distance(i), &
                 interpolation_nodes(:,i), &
                 interpolation_weights(:,i)
   ENDDO
   CLOSE(2)
   


END PROGRAM test_distances
