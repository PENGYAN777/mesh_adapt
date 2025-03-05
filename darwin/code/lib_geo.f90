MODULE lib_geo

 CONTAINS

 SUBROUTINE locpt (x0, y0, x, y, n, l, m)

 ! Given a polygonal line connecting the vertices (x(i),y(i)) (i = 1,...,n)
 ! taken in this order.  It is assumed that the polygonal path is a loop,
 ! where (x(n),y(n)) = (x(1),y(1)) or there is an arc from (x(n),y(n)) to
 ! (x(1),y(1)).  N.B. The polygon may cross itself any number of times.

 ! (X0,Y0) Is an arbitrary point and l and m are variables.
 ! On output, L and M are assigned the following values ...

 !    L = -1   IF (X0,Y0) IS OUTSIDE THE POLYGONAL PATH
 !    L =  0   IF (X0,Y0) LIES ON THE POLYGONAL PATH
 !    L =  1   IF (X0,Y0) IS INSIDE THE POLYGONAL PATH

 ! M = 0 if (X0,Y0) is on or outside the path.  If (X0,Y0) is inside the
 ! path then M is the winding number of the path around the point (X0,Y0).

 ! Converted to ELF90 compatibility by Alan Miller, 15 February 1997
 ! Customized by Marco Fossati (PMI), 02 August 2006

 !----------------------------------------------------------------------------------------
 IMPLICIT NONE

 REAL(KIND=8),               INTENT(IN) :: x0, y0
 REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x, y
 INTEGER,                    INTENT(IN) :: n

 INTEGER, INTENT(OUT) :: l, m

 INTEGER :: i, n0

 REAL(KIND=8) :: angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v
 !----------------------------------------------------------------------------------------

 !     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
 !            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0
 eps = EPSILON(1.0)


 n0 = n

 IF (x(1) == x(n) .AND. y(1) == y(n)) n0 = n - 1

 pi = ATAN2(0.0, -1.0)
 pi2 = 2.0*pi
 tol = 4.0*eps*pi
 l = -1
 m = 0

 u = x(1) - x0
 v = y(1) - y0

 IF (u == 0.0 .AND. v == 0.0) GO TO 20

 IF (n0 < 2) RETURN

 theta1 = ATAN2(v, u)

 sum = 0.0
 theta = theta1

 DO i = 2, n0

   u = x(i) - x0
   v = y(i) - y0
   IF (u == 0.0 .AND. v == 0.0) GO TO 20
   thetai = ATAN2(v, u)
   
   angle = ABS(thetai - theta)

   IF (ABS(angle - pi) < tol) GO TO 20
   IF (angle > pi) angle = angle - pi2
   IF (theta > thetai) angle = -angle

   sum = sum + angle
   theta = thetai

 END DO

 angle = ABS(theta1 - theta)
 IF (ABS(angle - pi) < tol) GO TO 20
 IF (angle > pi) angle = angle - pi2
 IF (theta > theta1) angle = -angle
 sum = sum + angle

 !     SUM = 2*PI*M WHERE M IS THE WINDING NUMBER

 m = ABS(sum)/pi2 + 0.2
 IF (m == 0) RETURN
 l = 1
 IF (sum < 0.0) m = -m

 RETURN

 !     (X0, Y0) IS ON THE BOUNDARY OF THE PATH

 20 l = 0

 RETURN

 END SUBROUTINE locpt





 SUBROUTINE Check_Boxes_Overlap (boxes)
 
 USE dynamic_vector

 !----------------------------------------------------------------------------------------
 IMPLICIT NONE

 TYPE(D_I_V_R), DIMENSION(:), INTENT(IN) :: boxes

 REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x, y

 REAL(KIND=8) :: x0, y0

 INTEGER :: i, j, k, l, m, n
 !----------------------------------------------------------------------------------------

 DO i = 1, SIZE(boxes)

   DO j = 1, SIZE(boxes)

     IF (i == j) CYCLE

     ALLOCATE (x(SIZE(boxes(i) % vec)/2), &
               y(SIZE(boxes(i) % vec)/2))

     x = (/ boxes(i) % vec(1:SIZE(boxes(i) % vec):2) /)
     y = (/ boxes(i) % vec(2:SIZE(boxes(i) % vec):2) /)

     n = SIZE(boxes(i) % vec) / 2

     DO k = 1, SIZE(boxes(j) % vec), 2
     
       x0 = boxes(j) % vec(k)
       y0 = boxes(j) % vec(k+1)

       CALL locpt (x0, y0, x, y, n, l, m)

       IF (l /= -1) THEN       
         PRINT*, ''                                                             
         PRINT*, 'Error. CHECK_BOXES_OVERLAP:'
         PRINT*, 'Box n°', i, 'and box n°', j, 'overlap.'
         PRINT*, ''
         STOP
       ENDIF

     ENDDO

     DEALLOCATE (x, y)

   ENDDO
 ENDDO

 END SUBROUTINE Check_Boxes_Overlap



END MODULE lib_geo
