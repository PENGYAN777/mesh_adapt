!============================================================ 
!
!      Module: umfpack_driver
!
! Description: driver for the UMFPACK direct linear solver 
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: guardone@aero.polimi.it
!
!   Copyright: 1998-2003 Alberto Guardone
!              See COPYING file for copyright notice
!
!============================================================ 

MODULE umfpack_driver



   !============================================================ 
   USE csr
   !============================================================ 

   IMPLICIT NONE

   INTEGER  ::  N_sys, N_nz,  N_index, N_value, N_work, &
                job_FA, job_RF, job_SO
   LOGICAL  ::  transp

   INTEGER, DIMENSION(20)       ::  keep, icntl
   INTEGER, DIMENSION(40)       ::  info
   REAL(KIND=8), DIMENSION(10)  ::  cntl
   REAL(KIND=8), DIMENSION(20)  ::  rinfo

   INTEGER,      DIMENSION(:), ALLOCATABLE  ::  index, index_
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE  ::  value
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE  ::  workspace
   
   LOGICAL  ::  factorize

   !============================================================ 
!   PUBLIC :: read_param_umfpack, init_umfpack, solve_umfpack
   PUBLIC ::  init_umfpack, solve_umfpack
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 



   !============================================================ 
   SUBROUTINE  init_umfpack(CSR_M) 
   !============================================================ 



      IMPLICIT NONE
      
      TYPE (CSR_matrix), INTENT(IN) :: CSR_M
      
      EXTERNAL   UMD21I

    
      INTEGER  ::  wlen, alen, acopy, add, dn, dne
      INTEGER  ::  i, j, e
      
      
      job_FA = 1   ! Preserve the matrix for iterative refinement
      job_RF = 1   ! Preserve the matrix for iterative refinement
      job_SO = 0   ! Allow for iterative refinement 
      transp = .FALSE.
      
      ! Allocation of the data vectors
      N_sys   = SIZE(CSR_M%i) - 1
      N_nz    = SIZE(CSR_M%j)
      N_value = 16*N_nz 
      N_work  =  4*N_sys
      
      dn      =    N_sys
      dne     =              N_nz
      wlen    = 11*N_sys          +  3*dn + 8
      alen    = 11*N_sys + 2*N_nz + 11*dn
      acopy   =    N_sys +   N_nz         + 1
      add     =  7*N_sys      
      
      N_index = wlen + alen + acopy + add

      
      OPEN (11, file = 'umfpack.log')
        WRITE (11,*) '   DIRECT SOLVER UMFPACK'
        WRITE (11,*) '   N_sys   = ', N_sys
        WRITE (11,*) '   N_nz    = ', N_nz
        WRITE (11,*) '   N_value = ', N_value
        WRITE (11,*) '   N_work  = ', N_work
        WRITE (11,*) '   N_index = ', N_index
      CLOSE (11)

      
      
      ALLOCATE ( index(N_index), value(N_value), workspace(N_work) )
      ALLOCATE ( index_(2*N_nz))
      
      ! Initialization of the index vector
      DO i = 1, SIZE(CSR_M%i) - 1   
         DO e = CSR_M%i(i), CSR_M%i(i+1) - 1;  j = CSR_M%j(e) 

             index_(e       ) = i  
             index_(e + N_nz) = j  

         ENDDO
      ENDDO
            
      ! Initialiazion of the UMFPACK library      
      CALL UMD21I( keep, cntl, icntl )
      icntl(4) = 0  ! Matrix is not reducible to block triangular
      icntl(6) = 1  ! Matrix has nearly synmmetric nonzero pattern
      icntl(8) = 10 ! Number of steps of the iterative refinement
      
      cntl(1)  = 0.d0  ! Pivoting: this way entries of A are 
                       ! ALL acceptable. Default is 0.1d0
      
      
      ! Factorization is to be computed      
      factorize = .TRUE. 
      
   END SUBROUTINE init_umfpack
   !============================================================ 



   !============================================================ 
   SUBROUTINE solve_umfpack(MM, rhs, sol) 
   !============================================================ 


      IMPLICIT NONE
      
      TYPE( CSR_matrix), INTENT(IN)  ::  MM

      REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: rhs
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: sol
      
      EXTERNAL   UMD2FA, UMD2SO

     
      value(1:N_nz) = MM%e

      IF (factorize) THEN 
     
          index(1:2*N_nz) = index_
          ! Factorization
          CALL UMD2FA( N_sys, N_nz, job_FA, transp, N_value, N_index, value, index, &
                       keep, cntl, icntl, info, rinfo                            )
      
          IF (info(1) < 0) THEN
             WRITE (*,*) ''
             WRITE (*,*) 'ERROR. SOLVE_UMFPACK:'
             WRITE (*,*) 'Error in factoring the matrix.'
             WRITE (*,*) ''
             STOP
          ENDIF
     
          factorize = .FALSE.
       
      ELSE 

          ! Factorization using the old profile
          index(1:2*N_nz) = index_
          CALL UMD2RF( N_sys, N_nz, job_RF, transp, N_value, N_index, value, index, &
                       keep, cntl, icntl, info, rinfo                            )
      
          IF ( info(1) < 0 ) THEN
          
!          WRITE(*,*) ' Warning: refactorizaion failed. Factorizing matrix...'
             value(1:N_nz) = MM%e
             index(1:2*N_nz) = index_
             ! Factorization
             CALL UMD2FA( N_sys, N_nz, job_FA, transp, N_value, N_index, value, index, &
                          keep, cntl, icntl, info, rinfo                            )

             IF ( info(1) < 0 ) THEN
               WRITE (*,*) ''
               WRITE (*,*) 'ERROR. SOLVE_UMFPACK:'
               WRITE (*,*) 'Error in factoring the matrix.'
               WRITE (*,*) ''
               STOP
             ENDIF

             factorize = .FALSE.
          
          ENDIF
            
      ENDIF
      
      ! Solution               
      CALL UMD2SO( N_sys, job_SO, transp, N_value, N_index, value, index, &
                    keep, rhs, sol, workspace, cntl, icntl, info, rinfo )
                    
      IF ( info(1) < 0 ) THEN
        WRITE (*,*) ''
        WRITE (*,*) 'ERROR. SOLVE_UMFPACK:'
        WRITE (*,*) 'Error in solving the matrix.'
        WRITE (*,*) ''
        STOP
      ENDIF
      
      
   END SUBROUTINE solve_umfpack
   !============================================================ 



END MODULE umfpack_driver
