MODULE distances

   USE grid_utils
   USE lin_algebra
   USE nodes
   USE mesh_structure
   USE element_topology
   USE dynamic_vector
   USE structures


 CONTAINS


   FUNCTION Compute_Distance(mesh, is_solid_boundary)   RESULT(wall_distance)
   !-------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(grid_type),       INTENT(INOUT) :: mesh
   LOGICAL, DIMENSION(:), INTENT(IN) :: is_solid_boundary


   REAL(KIND=8), DIMENSION(:),  ALLOCATABLE :: wall_distance

   ! ib = nearest_node(i): ib is the boundary node nearest to node i                                             
   INTEGER,      DIMENSION(:),  ALLOCATABLE :: nearest_node

   ! Boundary node -> domain node connectivity for single boundary nodes
   INTEGER,      DIMENSION(:),  POINTER :: jd_jb_single,   &
                                           bound_p_single

   TYPE(D_I_V),  DIMENSION(:),  POINTER  :: j_m_b_single

   INTEGER,      DIMENSION(:,:),  POINTER :: jb_jd

   REAL(KIND=8) :: d ! Distance

   INTEGER :: Np_b_single,   &  ! Number of singled boundary points
              i, i1, i2, i3, &  ! Domain node indices
              ib,            &  ! Boundary node index
              m,             &  ! Boundary element index 
              m_                ! Bubble element local index

   TYPE(D_I_V), DIMENSION(:), ALLOCATABLE :: m_j_b_single
   !------------------------------------------------------------------------- 

     ALLOCATE (j_m_b_single(SIZE(mesh % j_m_b)))

     DO i = 1, SIZE(j_m_b_single)

       ALLOCATE (j_m_b_single(i) % vec(SIZE(mesh % j_m_b(i) % vec)))

       j_m_b_single(i) % vec = mesh % j_m_b(i) % vec

     ENDDO


     CALL Invert_jd_jb(mesh)


     ! Delete double boundary nodes
!      CALL del_double_nodes_boundary (mesh % jd_jb, mesh % bound_p, &
!                                      j_m_b_single, mesh % bound_m, &
!                                      mesh % jb_jd, jd_jb_single)
 
     CALL del_double_nodes_boundary (mesh % jd_jb,   &
                                     mesh % jb_jd,   &
                                     mesh % j_m_b,   &
                                     mesh % bound_p, &
                                     mesh % bound_m, &
                                     jd_jb_single,   &
                                     j_m_b_single,   &
                                     bound_p_single)


     Np_b_single = SIZE(jd_jb_single)

     ! Compute the nearest boundary node nearest_node(i)
     ! with respect to domain node i
     ALLOCATE (wall_distance(mesh % Nj_d), &
                nearest_node(mesh % Nj_d))

     wall_distance = HUGE(d)
     nearest_node = 0

     DO i = 1, SIZE(mesh % rr,2)

        DO ib = 1, SIZE(jd_jb_single)
 
           IF (.NOT. is_solid_boundary(bound_p_single(ib))) CYCLE
 
           i1 = jd_jb_single(ib)

           d = distance_from_point(mesh % rr(:,i), mesh % rr(:,i1))

           IF (d < wall_distance(i)) THEN

              wall_distance(i) = d
              nearest_node(i)  = ib

           ENDIF
 
        ENDDO
     ENDDO

     IF (ALL(nearest_node == 0)) THEN

       PRINT*, ''
       PRINT*, 'ERROR. COMPUTE_DISTANCE:'
       PRINT*, 'Computing nodes distance from wall failed.'
       PRINT*, ''
       STOP

     ENDIF


     ! Creates the boundary node --> boundary elements
     ! connectivity m_j_b and allocate interpolation nodes
     ! and weights
     ALLOCATE(m_j_b_single(SIZE(jd_jb_single)))

     m_j_b_single = invert_DIV(j_m_b_single)
 

     ! Computes the wall distance and the interpolation
     ! nodes and weights for domain node i
     wall_distance = HUGE(d)

     DO i = 1, SIZE(mesh % rr, 2)
 
        ! Retrive the nearest boundary node ib
        ib = nearest_node(i)
 
        ! For all boundary elements belonging to the bubble
        ! of ib, computes the distance and the interpolation
        ! weights
        DO m_ = 1, SIZE(m_j_b_single(ib) % vec)

           m = m_j_b_single(ib) % vec(m_)
 
           SELECT CASE (mesh % ele_type_b(m))
 

              CASE(ELE_TYPE_SEGMENT)

                 i1 = jd_jb_single(j_m_b_single(m) % vec(1))
                 i2 = jd_jb_single(j_m_b_single(m) % vec(2))
 
                 d = distance_from_segment(mesh % rr(:,i),  &
                                           mesh % rr(:,i1), &
                                           mesh % rr(:,i2))
 
                 IF (d < wall_distance(i))   wall_distance(i) = d
 

              CASE(ELE_TYPE_TRIANGLE)

                 i1 = jd_jb_single(j_m_b_single(m) % vec(1))
                 i2 = jd_jb_single(j_m_b_single(m) % vec(2))
                 i3 = jd_jb_single(j_m_b_single(m) % vec(3))
 
                 d = distance_from_triangle(mesh % rr(:,i),  &
                                            mesh % rr(:,i1), &
                                            mesh % rr(:,i2), &
                                            mesh % rr(:,i3))
 
                 IF (d < wall_distance(i))   wall_distance(i) = d
 

              CASE(ELE_TYPE_QUADRILATER)

                 i1 = jd_jb_single(j_m_b_single(m) % vec(1))
                 i2 = jd_jb_single(j_m_b_single(m) % vec(2))
                 i3 = jd_jb_single(j_m_b_single(m) % vec(3))
 
                 d = distance_from_triangle(mesh % rr(:,i),  &
                                            mesh % rr(:,i1), &
                                            mesh % rr(:,i2), &
                                            mesh % rr(:,i3))
 
                 IF (d < wall_distance(i))   wall_distance(i) = d
 
                 i1 = jd_jb_single(j_m_b_single(m) % vec(3))
                 i2 = jd_jb_single(j_m_b_single(m) % vec(4))
                 i3 = jd_jb_single(j_m_b_single(m) % vec(1))
 
                 d = distance_from_triangle(mesh % rr(:,i),  &
                                            mesh % rr(:,i1), &
                                            mesh % rr(:,i2), &
                                            mesh % rr(:,i3))
 
                 IF (d < wall_distance(i))   wall_distance(i) = d

 
          END SELECT
 
        ENDDO
 
   ENDDO


 CONTAINS


   FUNCTION distance_from_triangle(P, V0, V1, V2,  a_)   RESULT(d)

   ! Note: the coefficients of equation 
   !   ax + by + cz + d = 0
   ! defining the plane are given by
   ! a = n(1), b = n(2), c = n(3), d = -SUM(n*V0)

   !-------------------------------------------------------------------------
   USE lin_algebra

   IMPLICIT NONE
  
   REAL(KIND=8), DIMENSION(3), INTENT(IN)  ::  P
   REAL(KIND=8), DIMENSION(3), INTENT(IN)  ::  V0, V1, V2
   REAL(KIND=8), DIMENSION(3), OPTIONAL    ::  a_
  
   REAL(KIND=8)  ::  d
  
   REAL(KIND=8), DIMENSION(2,3)    ::  dn 
   REAL(KIND=8), DIMENSION(3)      ::  I  
   REAL(KIND=8), DIMENSION(3)      ::  n,  u, v, w
   REAL(KIND=8), DIMENSION(3)      ::  dd
   REAL(KIND=8), DIMENSION(0:2)    ::  a
   REAL(KIND=8), DIMENSION(0:2,3)  ::  as
   REAL(KIND=8), DIMENSION(0:2)    ::  as_
   REAL(KIND=8)  ::  uu, uv, vv, wu, wv, den
   REAL(KIND=8)  ::  s, t

   INTEGER, DIMENSION(1)  ::  is_minv
   INTEGER  ::  is_min
   !-------------------------------------------------------------------------
          
     u = V1 - V0;  v = V2 - V0
 
     dn(1,:) = u;  dn(2,:) = v
     n = vec_prod(dn)        ! Normal to plane Pi
     n = n / SQRT(SUM(n**2)) ! --> to unit normal
 

     ! Intersection point I(x,y,z)
     I = P - SUM( n*(P-V0) ) * n
 
     ! Parametric coordinates (s, t) of the intersection point I
     w = I - V0
 
     uu = SUM( u*u );  vv = SUM( v*v )
     uv = SUM( u*v );  wu = SUM( w*u );   wv = SUM( w*v )
 
     den = uv**2  -  uu*vv
 
     s = (uv*wv - vv*wu) / den;  t = (uv*wu - uu*wv) / den
 
     ! Areal baricentric coordinates
     a(0) = 1 - s - t;  a(1) = s;  a(2) = t
 
     is_min = 0; dd = 0
     IF (ANY(a < 0)) THEN ! Point I lies outside the triangle
 
        as = 0.d0
 
        dd(1) = distance_from_segment(P, V0, V1,  as_ )
        as(0,1) = as_(0);  as(1,1) = as_(1);  as(2,1) = 0.d0

        dd(2) = distance_from_segment(P, V0, V2,  as_ )
        as(0,2) = as_(0);  as(2,2) = as_(1);  as(1,2) = 0.d0

        dd(3) = distance_from_segment(P, V1, V2,  as_ )
        as(1,3) = as_(0);  as(2,3) = as_(1);  as(0,3) = 0.d0

 
        is_minv = MINLOC(dd); is_min = is_minv(1)
 
        d = dd(  is_min)
        a = as(:,is_min)
 
        ELSE ! Point I is inside the triangle
 
        d = SQRT( SUM(P - I)**2 )
 
     ENDIF
 
     IF ( PRESENT(a_) ) a_ = a
        
   END FUNCTION distance_from_triangle
  
  



   FUNCTION distance_from_segment(P, V0, V1,  a_) RESULT(d)
   !-------------------------------------------------------------------------
   USE lin_algebra

   IMPLICIT NONE
   
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  P
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  V0, V1
   REAL(KIND=8), DIMENSION(2), INTENT(OUT), OPTIONAL ::  a_
   
   REAL(KIND=8)  ::  d
   
   REAL(KIND=8), DIMENSION(SIZE(P))  ::  I
   REAL(KIND=8), DIMENSION(SIZE(P))  ::  v, w
   REAL(KIND=8)  ::  c1, c2, s
   !-------------------------------------------------------------------------
 
     v = V1 - V0;  w = P - V0
     
     c1 = SUM(w*v)
     IF (c1 <= 0) THEN
        d = distance_from_point(P, V0)
        IF (PRESENT(a_)) THEN
           a_(1) = 1;  a_(2) = 0  
        ENDIF      
        RETURN
     ENDIF
     c2 = SUM(v*v)
     IF (c2 <= c1) THEN
        d = distance_from_point(P, V1)
        IF (PRESENT(a_)) THEN
           a_(1) = 0;  a_(2) = 1  
        ENDIF      
        RETURN
     ENDIF
     
     ! Parametric coordinates s of the intersection point I
     s = c1 / c2
     I = V0  +  s*v
     
     d = distance_from_point(P, I)
      
     ! Areal baricentric coordinates
     IF (PRESENT(a_)) THEN
        a_(1) = 1-s;  a_(2) = s  
     ENDIF      
             
   END FUNCTION distance_from_segment
  



    
   FUNCTION distance_from_point(P, P0) RESULT(d)
   !-------------------------------------------------------------------------
   USE lin_algebra

   IMPLICIT NONE
      
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  P
   REAL(KIND=8), DIMENSION(:), INTENT(IN)  ::  P0
      
   REAL(KIND=8)  ::  d
   !-------------------------------------------------------------------------
    
     d = SQRT(SUM((P-P0)**2))
             
   END FUNCTION distance_from_point
  
  END FUNCTION Compute_Distance

END MODULE distances
