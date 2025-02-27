MODULE boundary_geometry

 USE dynamic_vector 
 USE structures        
       
 CONTAINS
 
 ! -------------------------------------------------------------------
 ! -  Add_b_Nodes( grid, s_new, jd_jsb, rj_b )
 ! -  SPL_interp ( i, curv, se, xe )
 ! -  evaluate_crds( n, s, x, xs, se, xe )
 ! -  patch (n, s, sx, ds, i1, i2)
 ! -------------------------------------------------------------------



  FUNCTION Bnd_Displ (grid, s_new, jd_jsb, rj_b)   RESULT(D)

  !----------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,       DIMENSION(:),   INTENT(IN) :: rj_b
  TYPE(crv_crd), DIMENSION(:),   INTENT(IN) :: s_new
  TYPE(cnv_idx), DIMENSION(:),   INTENT(IN) :: jd_jsb

  TYPE(D_I_V_R),  DIMENSION(:), POINTER :: D

  INTEGER                                    :: i, Nb, b_idx, j
  INTEGER                                    :: nnd  
  INTEGER                                    :: nd, ns, p
  REAL(KIND=8),  DIMENSION(:),   ALLOCATABLE :: s
  REAL(KIND=8),  DIMENSION(:,:), ALLOCATABLE :: crds   
  TYPE(grid_type)                            :: grid  
  TYPE(curve),   DIMENSION(:),   ALLOCATABLE :: kurve         
  !----------------------------------------------------------------------------------


   ALLOCATE (D(grid % Nj_d))
    
   DO i = 1, grid % Nj_d

     ALLOCATE (D(i) % vec(2))
     D(i) % vec = 0.d0

   ENDDO
  
   
   Nb = grid % Nb
   ALLOCATE (kurve(Nb))
 
   DO i = 1, Nb
 
    kurve(i) % nd = grid % line(i) % nd
    kurve(i) % ns = grid % line(i) % ns
 
    nd = grid % line(i) % nd
    ns = grid % line(i) % ns
 
    ALLOCATE (kurve(i) % s(ns),    kurve(i) % h(ns))
    ALLOCATE (kurve(i) % x(nd,ns), kurve(i) % xs(nd,ns))
 
    kurve(i) % s  = grid % line(i) % s
    kurve(i) % h  = grid % line(i) % h
    kurve(i) % x  = grid % line(i) % x
    kurve(i) % xs = grid % line(i) % xs
 
   ENDDO
   
   DO i = 1, Nb  !-Loop over the lines that make the boundary

    IF (rj_b(i) == 0) CYCLE  !-No point added to line 'i'

    b_idx = i                  !-boundary line index
    nnd = kurve(b_idx) % nd    !-number of spatial dimensions

    ALLOCATE (s(rj_b(b_idx)))      !-curve coordinates of points to add to 
                                    ! boundary line `b_idx`.
                                    ! rj_b(b_idx) says for each boundary line
                                    ! how many new points are added

    ALLOCATE (crds(nnd, SIZE(s)))    !-cartesian coordinates for curve coords 's'

    !-Assign 's' with s_new computed in new_nodes
    s = s_new(b_idx) % points(1:SIZE(s))

! print*, ''
! print*, 'Bound', i
! do p = 1, size(s)
! print*, p, s(p)
! enddo

    CALL SPL_interp(b_idx, kurve, s, crds)
  
    ! jd_jsb returns for each new point on the curve the domain node index
    ! it is a vector and each component is associated to a boundary line

    DO j = 1, SIZE(s)          

      D(jd_jsb(b_idx) % index(j)) % vec(1) = crds(1,j) - grid % rr( 1, jd_jsb(b_idx) % index(j) )
      D(jd_jsb(b_idx) % index(j)) % vec(2) = crds(2,j) - grid % rr( 2, jd_jsb(b_idx) % index(j) )

    ENDDO        
    
    DEALLOCATE (s)
    DEALLOCATE (crds) 
   
   ENDDO

  END FUNCTION Bnd_Displ







  SUBROUTINE Add_b_Nodes( grid, s_new, jd_jsb, rj_b )

  ! This subroutine evaluate cartesian coordinates from curve coordinates using
  ! SPLINE interpolation for each boundary curve.
  ! NOT USED ANYMORE ---
  !----------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER,       DIMENSION(:),   INTENT(IN) :: rj_b
   TYPE(crv_crd), DIMENSION(:),   INTENT(IN) :: s_new
   TYPE(cnv_idx), DIMENSION(:),   INTENT(IN) :: jd_jsb


   INTEGER                                    :: i, Nb, b_idx, j
   INTEGER                                    :: nnd  
   INTEGER                                    :: nd, ns
   REAL(KIND=8),  DIMENSION(:),   ALLOCATABLE :: s
   REAL(KIND=8),  DIMENSION(:,:), ALLOCATABLE :: crds   
   TYPE(grid_type)                            :: grid  
   TYPE(curve),   DIMENSION(:),   ALLOCATABLE :: kurve         
  !----------------------------------------------------------------------------------
   
   
   Nb = grid%Nb
   ALLOCATE( kurve(Nb) )
 
   DO i = 1, Nb
 
    kurve(i)%nd = grid%line(i)%nd
    kurve(i)%ns = grid%line(i)%ns
 
    nd = grid%line(i)%nd
    ns = grid%line(i)%ns
 
    ALLOCATE( kurve(i)%s(ns), kurve(i)%h(ns) )
    ALLOCATE( kurve(i)%x(nd,ns), kurve(i)%xs(nd,ns) )
 
    kurve(i)%s  = grid%line(i)%s
    kurve(i)%h  = grid%line(i)%h
    kurve(i)%x  = grid%line(i)%x
    kurve(i)%xs = grid%line(i)%xs
 
   ENDDO
   
   DO i = 1, Nb  !-Loop over the lines that make the boundary

    IF ( rj_b(i) .EQ. 0 ) CYCLE  !-No point added to line 'i'

    b_idx = i                !-boundary line index
    nnd = kurve(b_idx)%nd    !-number of spatial dimensions

    ALLOCATE( s(rj_b(b_idx)) )      !-curve coordinates of points to add to 
                                    ! boundary line `b_idx`.
                                    ! rj_b(b_idx) says for each boundary line
                                    ! how many new points are added

    ALLOCATE( crds( nnd,SIZE(s) ) )  !-cartesian coordinates for curve coords 's'


    s = s_new(b_idx)%points(1:SIZE(s))  !-Assign 's' with s_new computed in 
                                        ! new_nodes

    CALL SPL_interp( b_idx, kurve, s, crds )

  
    ! jd_jsb returns for each new point on the curve the domain node index
    ! it is a vector and each component is associated to a boundary line
  
    DO j = 1, SIZE(s)          
      grid % rr(1, jd_jsb(b_idx) % index(j)) = crds(1,j)
      grid % rr(2, jd_jsb(b_idx) % index(j)) = crds(2,j)            
    ENDDO        
    
    DEALLOCATE(s)
    DEALLOCATE(crds) 
   
   ENDDO

  END SUBROUTINE Add_b_Nodes      
      
  
  
    
 
  SUBROUTINE SPL_interp ( i, curv, se, xe )
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Compute cartesian coordinates 'xe' corresponding to curve coordinate 'se'    
  !=================================================================================!
   IMPLICIT NONE

   INTEGER,                     INTENT(IN) :: i
   REAL (KIND=8), DIMENSION(:), INTENT(IN) :: se
   TYPE(curve),   DIMENSION(:), INTENT(IN) :: curv
 
   REAL (KIND=8), DIMENSION(:,:), INTENT(OUT) :: xe
  !=================================================================================!
      
    CALL evaluate_crds( curv(i)%ns, curv(i)%s, curv(i)%x, curv(i)%xs, se, xe )

  END SUBROUTINE SPL_interp  
  
  
  


  SUBROUTINE evaluate_crds( n, s, x, xs, se, xe )
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Evaluate hermite univariate interpolation at point 'se'.'se' ranges from zero
  ! to one.        
  ! Warning: to be updatated in order to eliminate "n"  
  !=================================================================================!
   IMPLICIT NONE

   INTEGER,                       INTENT(IN)  :: n
   REAL (KIND=8), DIMENSION(:),   INTENT(IN)  :: s, se
   REAL (KIND=8), DIMENSION(:,:), INTENT(IN)  :: x, xs


   INTEGER       :: i1, i2, k, j
   REAL (KIND=8) :: uu, du, a1, b1, a2, b2, h01, h02, h03, h04


   REAL (KIND=8), DIMENSION(:,:), INTENT(OUT) :: xe
  !=================================================================================!
  
    DO j = 1, SIZE( se )
    
      uu  = se(j)

      CALL patch (n, s, uu, du, i1, i2)

      a1  = uu
      b1  = 1.-uu
      a2  = a1**2
      b2  = b1**2

      h01 = (b1+3.*a1)*b2
      h02 = +a1*b2
      h03 = -a2*b1
      h04 = (a1+3.*b1)*a2

      DO k = 1, SIZE(x,1)       
        xe(k,j) = x(k,i1)*h01 + du*xs(k,i1)*h02 + du*xs(k,i2)*h03 + x(k,i2)*h04
      ENDDO

    ENDDO

  END SUBROUTINE evaluate_crds
    
    
    
    
    
    
    
  SUBROUTINE patch (n, s, sx, ds, i1, i2)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Identify couple of nodes immediatly after and before the new interpolated point
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    IMPLICIT NONE
  !=================================================================================!
    INTEGER,                     INTENT(IN)    :: n
    REAL (KIND=8), DIMENSION(:), INTENT(IN)    :: s
    REAL (KIND=8),               INTENT(INOUT) :: sx
    REAL (KIND=8),               INTENT(OUT)   :: ds
    INTEGER,                     INTENT(OUT)   :: i1, i2

    INTEGER                                    :: l, h, r
    REAL (KIND=8), PARAMETER                   :: one=1.d0, zero=0.d0 
  !=================================================================================!

    IF ( ( sx .LT. zero ) .OR. ( sx .GT. one ) ) THEN

       WRITE (*,*) 'patch warning, parameter out of range'
       WRITE (*,*) 'parameter value = ',sx
       WRITE (*,*) 'parameter has been corrected'

       sx = MIN(one,MAX(zero,sx))

    ENDIF

    l = 1; r = n; sx = s(n)*sx

    DO; IF ( r-l == 1 ) EXIT
       h = (l+r)/2
       IF ( s(h) > sx ) THEN
          r = h
       ELSE
          l = h
       ENDIF
    ENDDO

    i1 = l; i2 = r; ds = s(r)-s(l); sx = (sx-s(l))/ds

  END SUBROUTINE patch      
       

END MODULE boundary_geometry
