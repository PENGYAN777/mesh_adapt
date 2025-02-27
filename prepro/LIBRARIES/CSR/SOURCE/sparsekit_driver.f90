!============================================================ 
!
!      Module: solve_skit
!
! Description: driver for the SPARKIT iterative linear solver 
!
!     Authors: J.-L. Guermond
!              LIMSI
!              Orsay, FRANCE
!              Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!   Cosmetics: L. Quartapelle
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!
!============================================================ 

MODULE  solve_skit

   USE csr
   
   IMPLICIT NONE
 
   PRIVATE
 
!  From module csr:
!   TYPE  CSR_Matrix
!      INTEGER,      DIMENSION(:), POINTER :: i    
!      INTEGER,      DIMENSION(:), POINTER :: j
!      REAL(KIND=8), DIMENSION(:), POINTER :: e ! matrix elements
!   END TYPE  CSR_Matrix


   TYPE  MSR_Matrix
      INTEGER,      DIMENSION(:), POINTER :: i ! Different from CSR !   
      INTEGER,      DIMENSION(:), POINTER :: j
      REAL(KIND=8), DIMENSION(:), POINTER :: e ! matrix elements
   END TYPE  MSR_Matrix

   INTEGER,      DIMENSION(7) :: ipar
   REAL(KIND=8), DIMENSION(4) :: fpar
   INTEGER      ::  iout

   INTEGER,      PARAMETER :: max_save   = 10,  &
                              max_Krylov = 50,  &
                              std_Krylov = 15,  &
                              std_max_it = 500

   REAL(KIND=8), PARAMETER :: std_eps_rel = 1.d-6,  &
                              std_eps_abs = 1.d-9   
   
   !  Defined as static vector, to simplify acces and to improve speed

   TYPE(MSR_Matrix), DIMENSION(max_save), SAVE :: LU_stored
   
   LOGICAL  :: output_opened = .FALSE.
  
   PUBLIC :: init_spkit,  solve_spkit,  cancel_mem_precond, CSR_Matrix
   
   EXTERNAL cg


CONTAINS
!=======


!**************************************************************************

  
SUBROUTINE  init_spkit(idf)
  
   IMPLICIT NONE
  
   INTEGER, INTENT(IN) :: idf
      
   INTEGER :: i
      
   DO i = 1, SIZE(ipar);  READ (idf,*) ipar(i);  ENDDO
      
   DO i = 1, SIZE(fpar);  READ (idf,*) fpar(i);  ENDDO
     
END SUBROUTINE  init_spkit

  
!**************************************************************************


SUBROUTINE  solve_spkit(A, rhs, sol,  isave)

      !-----------------------------------------------------------------------
      !     Program for ilu preconditioned gmres.
      !     This program solves a linear system with a sparse matrix.
      !     Different methods available:
      !     ILU0, MILU0, ILUT, ILUK, ILUD  
      !     with different values of tol and lfil
      !     (from cheaper to more expensive preconditioners).
      !     The more accurate the preconditioner the fewer iterations
      !     are required in pgmres, in general.
      !-----------------------------------------------------------------------
      !
      !     ipar(1) = 1     ----> ILU0
      !             = 2     ----> MILU0
      !             = 3     ----> ILUT
      !             = 4     ----> ILUTP
      !             = 5     ----> ILUK
      !             = 6     ----> ILUD
      !             = 7     ----> GMRES or CG
      !
      !     ipar(2) = lfil  ----> If lfil < 0, default value is chosen.
      !                     ----> Needed for ILUT, ILUTP: lfil is the number
      !                           of non zero entries/row in the
      !                           LU decomposition (Default=15).
      !                     ----> Needed for ILUK: Level of fill used in
      !                           the LU decomposition (Default=6).
      !
      !     ipar(3) = n_work ---> Work space for the ILU decomposition
      !                           If n_work < 0, default value is chosen.
      !
      !     ipar(4) = iout  ----> Output unit number for printing intermediate
      !                           results. If iout <= 0 nothing is printed.
      !
      !     ipar(5) = im    ----> Size of Krylov space
      !
      !     ipar(6) = max_its --> Maximal number of iterations
      !
      !     ipar(7) = 1     ----> GMRES is used
      !             = 2     ----> CG is used
      !
      !     fpar(1) = eps_rel --> Stopping criterion:
      !                           ||current residual||/||initial residual|| <= eps
      !                           or ||current residual|| <= eps_absolute
      !                           Euclidian norm is used (If eps<0, default set).
      !
      !     fpar(2) = eps_abs --> Stopping criterion:
      !
      !     fpar(3) = tol   ----> Threshold for dropping term in ILU factorization
      !                           If tol < 0, default value chosen.
      !
      !     fpar(4) = alpha ----> Used only by ILUD: Diagonal conpensation parameter
      !                            0 =< alpha <= 1 (IF alpha < 0, default is chosen).
      !
      !     A%i             ----> Pointers to the beginning of each row of A.
      !     A%j             ----> Pointers to matrix elements of A stored in CSR.
      !     A%e             ----> Matrix elements of A stored in CSR.
      !
      !     rhs             ----> Right Hand Side
      !     sol             ----> Input: initial guess;  Output: solution.
      !
      !
      !     isave           ----> OPTIONAL: INTEGER
      !                           If isave < 0, the preconditioner is saved to
      !                                         be referred with index |isave|.
      !                           If isave > 0, the saved-preconditioner is
      !                                         reused.
      !-----------------------------------------------------------------------
            
   TYPE(CSR_Matrix),           INTENT(INOUT) :: A  !  A%j  only

   REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: rhs
   REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: sol
   INTEGER,      OPTIONAL,     INTENT(IN)    :: isave
                !========

   TYPE(MSR_Matrix) :: LU
  
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: vv
  
   INTEGER      :: k, index_save, method, ierr, lfil,  &
                   n_syst, n_work, im, max_its, n_bloc,  LUi_size
     
   REAL(KIND=8) :: tol, perm_tol, eps_rel, eps_abs, alpha

   INTEGER,      DIMENSION(16) :: ipar_loc
   REAL(KIND=8), DIMENSION(16) :: fpar_loc
  

   n_syst = SIZE(A%i) - 1  ! Dimension of the linear system
   
   !-----CHECK PARAMETERS IPAR - FPAR

   IF (ipar(1) <= 0  .OR.  ipar(1) > 6) THEN
         
      method = 3
      alpha = -1.
      tol =  -1.
      n_work = -1

      iout = ipar(4)
      WRITE (iout,*) 'Warning: method not existent -> default parameters used'
      WRITE (iout,*) 'Method:', method

   ELSE

      method = ipar(1)
      lfil   = ipar(2)
      n_work = ipar(3)

      iout = ipar(4)
      IF ( (.NOT. output_opened) .AND. (ipar(4) /= 6) ) THEN
         OPEN (iout, file = 'sparse.log')
         output_opened = .TRUE.
      ENDIF


      IF (ipar(5) <= 0  .OR.  ipar(5) >= max_Krylov) THEN
         WRITE (iout,*) 'Krylov space set to', std_Krylov
         im = std_Krylov
      ELSE
         im = ipar(5)
      ENDIF


      IF (ipar(6) <= 0) THEN
         WRITE (iout,*) 'MAXIMUM Iteration number set to', std_max_it
         max_its = std_max_it
      ELSE
         max_its = ipar(6)
      ENDIF


      IF (fpar(1) <= 0) THEN
          WRITE (iout,*) 'Relative epsilon set to', std_eps_rel
          eps_rel = std_eps_rel
      ELSE
         eps_rel = fpar(1)
      ENDIF


      IF (fpar(2) <= 0) THEN
         WRITE (iout,*) 'Absolute epsilon set to', std_eps_abs
         eps_abs = std_eps_abs
      ELSE
         eps_abs = fpar(2)
      ENDIF

      tol = fpar(3)
      alpha = fpar(4)

   ENDIF

   !  Test the RHS for NULL SOLUTION
   !IF (DOT_PRODUCT(rhs,rhs) <= eps_abs**2) THEN
   !   sol = 0
   !   IF (.NOT.PRESENT(isave)  .OR.  isave >= 0) RETURN
   !ENDIF

   ! Tolerance ratio used to determine wether or not to permute two columns

   perm_tol = 0.05

   ! Permuting is done within the diagonal blocs of size n_bloc. 
   ! Useful only if several degrees of freedom of the PDE are solved at once.

   n_bloc = n_syst

   !  Time to precondition


   !  Check the preconditioner status

   IF (PRESENT(isave)) THEN

      IF (ABS(isave) > max_save) THEN 
         WRITE (*,*) ''
         WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
         WRITE (*,*) 'isave out of range (maximum is ', max_save, ')' 
         WRITE (*,*) ''
         STOP
      ENDIF


      IF (isave >= 1) THEN

         ! Load saved preconditioner

         IF (ASSOCIATED(LU_stored(isave)%i)) THEN

            LUi_size = SIZE(LU_stored(isave)%i)

            IF (n_syst /= LUi_size) THEN
               WRITE (*,*) ''
               WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
               WRITE (*,*) 'Incoherence in the size of the linear system'
               WRITE (*,*) ' n_syst = ', n_syst, ' LU%i_size', LUi_size
               WRITE (*,*) ''
               STOP
            ENDIF

            n_work = SIZE(LU_stored(isave)%j)
            
            LU = LU_stored(isave)
            
         ELSE

            WRITE (*,*) ''
            WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
            WRITE (*,*) 'The problem referred to by isave = ', isave
            WRITE (*,*) 'has not yet been preconditioned.'
            WRITE (*,*) 'Use a negative integer first.'
            WRITE (*,*) ''
            STOP

         ENDIF

      ELSE

         ! Preconditioning

         index_save = ABS(isave)
            
         IF (ASSOCIATED(LU_stored(index_save)%i)) THEN
            WRITE (iout,*) 'Place already used! - Deleting old PRECONDITIONER'
            CALL cancel_mem_precond(index_save)
         ENDIF
            
         CALL precond_it(method)  !  it computes LU 

      ENDIF

   ELSE

      ! Preconditioning

      CALL precond_it(method)  !  it computes LU 

   ENDIF


   SELECT CASE (ipar(7))

      CASE (1)   !     == PGMRES ==

         ALLOCATE (vv(n_syst*(im+1)))
      
         CALL pgmres(n_syst, im, rhs, sol, vv, eps_rel, eps_abs, max_its, iout,  &  
                     A%e, A%j, A%i, LU%e, LU%j, LU%i, ierr)

      CASE (2)   !     == Conjugate Gradient ==
            
         ipar_loc(2) = 2
         ipar_loc(3) = 1
         ipar_loc(4) = 5*n_syst ! size of work space
         ipar_loc(5) = 2
         ipar_loc(6) = max_its
         fpar_loc(1) = eps_rel
         fpar_loc(2) = eps_abs
         
         ALLOCATE (vv(5*n_syst))
         
         CALL runrc(n_syst, rhs, sol, ipar_loc, fpar_loc, vv, sol, &
                    A%e, A%j, A%i, LU%e, LU%j, LU%i, cg)     !  guess
         
         IF (ipar_loc(1) == 0)  ierr = 0
         IF (ipar_loc(1) == -1) ierr = 1

      CASE DEFAULT

         WRITE (*,*) ''
         WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
         WRITE (*,*) 'Solver pointed to by ipar(7) not programmed yet'
         WRITE (*,*) ''
         STOP

   END SELECT


   SELECT CASE(ierr)

      CASE(1);       WRITE (*,*) ''
                     WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                     WRITE (*,*) 'Convergence not achieved in ', max_its,' iterations.'
                     WRITE (*,*) ''
                     STOP

      CASE(-1);      WRITE (iout,*) ' SOLVE_SPKIT: Initial guess was OK'

      CASE(0);       WRITE (iout,*) ' SOLVE_SPKIT: Convergence is OK'

      CASE DEFAULT;  WRITE (iout,*) ' SOLVE_SPKIT: Undetermined problem'

   END SELECT 


   WRITE (iout,*) 


   !  Save Preconditioner

   IF (PRESENT(isave)) THEN
         
      IF (isave < 0) THEN

         index_save = ABS(isave)
                                
         LU_stored(index_save) = LU
        
      ENDIF
     
   ELSE
   
      DEALLOCATE (LU%e, LU%j, LU%i)
              
   ENDIF

   DEALLOCATE (vv)

   !********************************************************************************
      
   CONTAINS
     
   !********************************************************************************
      
   SUBROUTINE  precond_it(method) ! internal subroutine

      INTEGER, INTENT(IN) :: method

      INTEGER,      DIMENSION(:), ALLOCATABLE :: iw, perm, levs
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: wk


      ALLOCATE (LU%i(n_syst))


      SELECT CASE (method)
        
!__________________________________________________________________________________
        
      CASE (1)

         WRITE (iout,*) ' +++++ ILU(0) Preconditioner ++++ '
         n_work = SIZE(A%e) + n_syst
         ALLOCATE (LU%e(n_work), LU%j(n_work), iw(n_syst))
          
         CALL ilu0 (n_syst, A%e, A%j, A%i, LU%e, LU%j, LU%i, iw, ierr)
          
         IF (ierr /= 0) THEN
           WRITE (*,*) ''
           WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
           WRITE (*,*) 'Zero pivot at step ', ierr
           WRITE (*,*) ''
           STOP
         ENDIF
           
!__________________________________________________________________________________
        
      CASE (2)

         WRITE (iout,*) ' +++++ MILU(0) Preconditioner ++++ '
         n_work = SIZE(A%e) + n_syst
         ALLOCATE (LU%e(n_work), LU%j(n_work), iw(n_syst))
           
         CALL milu0 (n_syst, A%e, A%j, A%i, LU%e, LU%j, LU%i, iw, ierr)
           
         IF (ierr /= 0) THEN
           WRITE (*,*) ''
           WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
           WRITE (*,*) 'Zero pivot at step ', ierr
           WRITE (*,*) ''
           STOP
         ENDIF
           
!__________________________________________________________________________________
        
      CASE (3)

         WRITE (iout,*) ' +++++ ILUT Preconditioner ++++ '
         IF (lfil < 0) lfil = 30
         IF (tol  < 0) tol = 0.0001d0
         WRITE (iout,*) ' +++++ tol = ', tol, ' lfil = ', lfil,'++++ '
         n_work = (2*lfil + 1) * n_syst
         ALLOCATE (wk(n_syst), iw(2*n_syst))

         DO
            ALLOCATE (LU%e(n_work), LU%j(n_work))
              
            CALL ilut (n_syst, A%e, A%j, A%i, lfil, tol,  &
                              LU%e,LU%j,LU%i, n_work, wk, iw, ierr)
              
            SELECT CASE (ierr)
              
               CASE(0);   EXIT

               CASE(1:);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Zero pivot at step ', ierr
                          WRITE (*,*) ''
                          STOP

               CASE(-1);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Input matrix is wrong'
                          WRITE (*,*) ''
                          STOP

                          
                          
                 
               CASE(-3:-2) 
                  WRITE (iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorization'
                  WRITE (iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
                  n_work = 2 * n_work
                  DEALLOCATE (LU%e, LU%j)
                  CYCLE
                 
               CASE(-4);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Illegal value of lfil', lfil
                          WRITE (*,*) ''
                          STOP

  
               CASE DEFAULT;  WRITE (*,*) ''
                              WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                              WRITE (*,*) 'Zero row encountered'
                              WRITE (*,*) ''
                              STOP

                 
            END SELECT
              
         ENDDO
           
!__________________________________________________________________________________
    
      CASE (4)

         WRITE (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         IF (lfil < 0) lfil = 30
         IF (tol  < 0) tol = 0.0001d0
         WRITE (iout,*) ' +++++ tol = ',tol,' lfil = ',lfil,'++++ '
         n_work = (2*lfil + 1) * n_syst
         ALLOCATE (wk(n_syst), iw(2*n_syst), perm(2*n_syst))

         DO
            ALLOCATE (LU%e(n_work), LU%j(n_work))
              
            CALL ilutp (n_syst, A%e, A%j, A%i, lfil, tol, perm_tol, n_bloc,  &
                               LU%e,LU%j,LU%i, n_work, wk, iw, perm, ierr)

            ! TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH 
            ! LU-SOLVE, THE MATRIX A IS PERMUTED ON RETURN. 
            ! [all column indices are changed]. 
            ! SIMILARLY FOR THE LU MATRIX.
            ! To permute the matrix back to its original state use the loop:

            DO k = A%i(1),  A%i(n_syst+1) - 1
               A%j(k) = perm(A%j(k))
            ENDDO

            SELECT CASE(ierr)

               CASE(0);   EXIT 
              
               CASE(1:);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Zero pivot at step ', ierr
                          WRITE (*,*) ''
                          STOP

                  
               CASE(-1);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Input matrix is wrong'
                          WRITE (*,*) ''
                          STOP

                 
               CASE(-3:-2)
                  WRITE (iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorization'
                  WRITE (iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
                  n_work = 2 * n_work
                  DEALLOCATE (LU%e, LU%j)
                  CYCLE
              
               CASE(-4);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Illegal value of lfil', lfil
                          WRITE (*,*) ''
                          STOP

                   
               CASE(:-5); WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Zero row encountered'
                          WRITE (*,*) ''
                          STOP

            
            END SELECT
              
         ENDDO
           
!__________________________________________________________________________________
        
      CASE (5)

         WRITE (iout,*) ' +++++ ILUK Preconditioner ++++ '
         IF (lfil < 1  .OR.  lfil > 20) lfil = 6
         WRITE (iout,*) ' +++++  lfil = ',lfil,'++++ '
         IF (n_work <= 0) n_work = 60 * lfil * n_syst
         ! ALLOCATE (LU%i(n_syst), wk(n_syst), iw(3*n_syst))  !  error by JLG (?)
         ALLOCATE (wk(n_syst), iw(3*n_syst))

         DO
            ALLOCATE (LU%e(n_work+n_syst), LU%j(n_work), levs(n_work))

            CALL iluk (n_syst, A%e, A%j, A%i, lfil,  &
                              LU%e,LU%j,LU%i, levs, n_work, wk, iw, ierr)

            SELECT CASE(ierr)

               CASE(0);   EXIT

               CASE(1:);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Zero pivot at step ', ierr
                          WRITE (*,*) ''
                          STOP

                  
               CASE(-1);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Input matrix is wrong'
                          WRITE (*,*) ''
                          STOP
               
               CASE(-3:-2)
                  WRITE (iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorization'
                  WRITE (iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
                  n_work = 2 * n_work
                  DEALLOCATE (LU%e, LU%j, levs)
                  CYCLE
             
               CASE(-4);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Illegal value of lfil', lfil
                          WRITE (*,*) ''
                          STOP

                   
               CASE(:-5); WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Zero row encountered'
                          WRITE (*,*) ''
                          STOP
              
            END SELECT 
          
         ENDDO
         
!__________________________________________________________________________________
        
      CASE (6)

         WRITE (iout,*) ' +++++ ILUD Preconditioner ++++ '
         IF (tol < 0) tol = 0.01
         IF (alpha < 0  .OR.  alpha > 1) alpha = 0
         WRITE (iout,*) ' +++++ tol = ', tol,' alpha = ', alpha,'++++ '
         IF (n_work <= 0) n_work = 100 * n_syst
         ! ALLOCATE (LU%i(n_syst), wk(n_syst), iw(2*n_syst))  !  error by JLG (?)
         ALLOCATE (wk(n_syst), iw(2*n_syst))

         DO
            ALLOCATE (LU%e(n_work+n_syst), LU%j(n_work))
          
            CALL ilud (n_syst, A%e, A%j, A%i, alpha, tol,  &
                              LU%e,LU%j,LU%i, n_work, wk, iw, ierr)

            SELECT CASE(ierr)

               CASE(0);   EXIT

               CASE(1:);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Zero pivot at step ', ierr
                          WRITE (*,*) ''
                          STOP

                  
               CASE(-1);  WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Input matrix is wrong'
                          WRITE (*,*) ''
                          STOP
              
               CASE(-2)  
                  WRITE (iout,*) ' SOLVE_SPKIT: Insufficient storage for LU factorization'
                  WRITE (iout,*) ' SOLVE_SPKIT: n_work is set to 2*n_work'
                  n_work = 2 * n_work
                  DEALLOCATE (LU%e, LU%j)
                  CYCLE
              
               CASE(:-3); WRITE (*,*) ''
                          WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
                          WRITE (*,*) 'Zero row encountered'
                          WRITE (*,*) ''
                          STOP
               
            END SELECT
           
         ENDDO
           
!__________________________________________________________________________________

      CASE DEFAULT

         WRITE (*,*) ''
         WRITE (*,*) 'ERROR. SOLVE_SPKIT:'
         WRITE (*,*) 'Preconditioning method', method, ' not existent!'
         WRITE (*,*) ''
         STOP

      END SELECT

      IF (ALLOCATED(iw))   DEALLOCATE (iw)
      IF (ALLOCATED(wk))   DEALLOCATE (wk)
      IF (ALLOCATED(perm)) DEALLOCATE (perm)
      IF (ALLOCATED(levs)) DEALLOCATE (levs)

   END SUBROUTINE  precond_it ! internal subroutine


END SUBROUTINE  solve_spkit
 
 
!********************************************************************************
 
 
SUBROUTINE  cancel_mem_precond(index)

   INTEGER, INTENT(IN) :: index

   INTEGER :: well_done


   IF (ASSOCIATED(LU_stored(index)%i)) THEN

      DEALLOCATE (LU_stored(index)%i,  STAT = well_done)
      IF (well_done /= 0) THEN
         WRITE (*,*) ''
         WRITE (*,*) 'ERROR. CANCEL_MEM_PRECOND:'
         WRITE (*,*) 'Deallocation of stored preconditioner failed (LU%i)'
         WRITE (*,*) ''
         STOP
      ENDIF
         
      DEALLOCATE (LU_stored(index)%j,  STAT = well_done)
      IF (well_done /= 0) THEN
         WRITE (*,*) ''
         WRITE (*,*) 'ERROR. CANCEL_MEM_PRECOND:'
         WRITE (*,*) 'Deallocation of stored preconditioner failed (LU%j)'
         WRITE (*,*) ''
         STOP
      ENDIF
         
      DEALLOCATE (LU_stored(index)%e,  STAT = well_done)
      IF (well_done /= 0) THEN
         WRITE (*,*) ''
         WRITE (*,*) 'ERROR. CANCEL_MEM_PRECOND:'
         WRITE (*,*) 'Deallocation of stored preconditioner failed (LU%e)'
         WRITE (*,*) ''
         STOP
      ENDIF

   ENDIF


END SUBROUTINE  cancel_mem_precond


!********************************************************************************


END MODULE  solve_skit

