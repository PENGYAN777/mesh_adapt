!=============================================================================== 
!
!      Module: plot_procedures
!
! Description: Procedures for posprocessing results. Two different
!              plotting programs can be selected:
!
!               - Tecplot
!               - Plotmtv
!
!      Author: Alberto Guardone
!              Dipartimento di Ingegneria Aerospaziale
!	       Politecnico di Milano
! 	       Via La Masa 34, 20156 Milano, ITALY
! 	       e-mail: guardone@aero.polimi.it
!
!              Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!   Copyright: 2003-2006 Marco Fossati
!              See COPYING file for copyright notice
!
!===============================================================================

   MODULE  plot_procedures

   !-----------------------------------------------------------------------------
   USE dynamic_vector

   IMPLICIT NONE
   
   TYPE plot_data_type
   
      CHARACTER(LEN=64) :: plot_name
      CHARACTER(LEN=64) :: file_name
      INTEGER :: pb_l, var_l, fl_l

      REAL(KIND=8),     DIMENSION(:,:), POINTER :: variables
      CHARACTER(LEN=6), DIMENSION(:),   POINTER :: var_names
      
      TYPE(D_I_V),  DIMENSION(:),   POINTER :: j_m_d	  
      TYPE(D_I_V),  DIMENSION(:),   POINTER :: j_m_b	  
      INTEGER,      DIMENSION(:),   POINTER :: jd_jb
      TYPE(D_I_V),  DIMENSION(:),   POINTER :: p_bound      
      TYPE(D_I_V),  DIMENSION(:),   POINTER :: m_bound      
      TYPE(D_I_V),  DIMENSION(:),   POINTER :: p_type	   
      INTEGER,      DIMENSION(:),   POINTER :: boundary_types
      
      REAL(KIND=8), DIMENSION(:,:), POINTER :: rr
      REAL(KIND=8), DIMENSION(:),   POINTER :: time	       
   
   END TYPE plot_data_type
   

   ! Tecplot Options
   INTEGER, PARAMETER  ::  TEC_STANDARD = 0, &
                           TEC_OPEN     = 1, &
                           TEC_WRITE    = 2, &
                           TEC_CLOSE    = 3   

   INTEGER,      PARAMETER  ::  Debug = 0,      &
                                VIsDouble = 1,  &
                                DIsDouble = 1

   CHARACTER(*), PARAMETER  ::  scrtch_dir = '.'
   CHARACTER(1), PARAMETER  ::  NULLCHAR = CHAR(0)   
   !-----------------------------------------------------------------------------

   CONTAINS


   !-----------------------------------------------------------------------------
   ! PLOTMTV
   !-----------------------------------------------------------------------------
     
   FUNCTION  write_plot_file_plotmtv(plot_data)   RESULT(err)
   !-----------------------------------------------------------------------------	
   IMPLICIT NONE 
   
   TYPE(plot_data_type), INTENT(IN) :: plot_data

   INTEGER ::  err
   
   INTEGER :: idf = 11, space_dim, m, i, i_, nFields
   !-----------------------------------------------------------------------------	

   ! Number of space dimensions
   space_dim = SIZE(plot_data % rr, 1)


   OPEN (UNIT=idf, file=plot_data%file_name)
   
     SELECT CASE (space_dim)
     
       CASE (1)
       
         DO nFields = 1, SIZE(plot_data%var_names)-1
	 
	   WRITE(idf,*) '$ DATA = CURVE2D'
           WRITE(idf,*) '% toplabel = "'//plot_data%var_names(nFields+1)//'"'
	   WRITE(idf,*) '% markertype = 12'
	   
           ! Domain
           DO i = 1, SIZE(plot_data%rr,2)
             WRITE (idf,101) plot_data%rr(1,i), plot_data%variables(nFields,i), i
	   ENDDO

           101 FORMAT(2(e15.9,3x),i5)
           WRITE(idf,*) ''
       
         ENDDO



       CASE (2)

         DO nFields = 1, SIZE(plot_data%var_names)-2

           WRITE(idf,*) '$ DATA = CONTCURVE'
           WRITE(idf,*) '% toplabel = "'//plot_data%var_names(nFields+2)//'"'
           WRITE(idf,*) '% contstyle = 2'
           WRITE(idf,*) '% nsteps = 50'
           WRITE(idf,*) '% meshplot = false'

           ! Domain
           DO m = 1, SIZE(plot_data%j_m_d)
             DO i_ = 1, SIZE(plot_data%j_m_d(m)%vec)

	       i = plot_data%j_m_d(m)%vec(i_)
               WRITE (idf,100) plot_data%rr(1,i), plot_data%rr(2,i), plot_data%variables(nFields,i), i
          
	     ENDDO
	     WRITE(idf,*) ''
	   ENDDO

           100 FORMAT(3(e15.9,3x),i5)
           WRITE(idf,*) ''

         ENDDO



       CASE (3)
              
     END SELECT
   
   CLOSE (idf)

   err = 0

   END FUNCTION  write_plot_file_plotmtv







   !-----------------------------------------------------------------------------
   ! TECPLOT
   !----------------------------------------------------------------------------- 

   FUNCTION  write_plot_file_tecplot(plot_data, tec_mode)   RESULT(err)
   !-----------------------------------------------------------------------------	
   IMPLICIT NONE 
   
   TYPE(plot_data_type), INTENT(IN) :: plot_data
   INTEGER, OPTIONAL :: tec_mode
   
   INTEGER ::  err

   INTEGER, EXTERNAL :: TecIni, TecZne, TecDat, TecNod, TecEnd
   INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  jm

   INTEGER :: i, i_, m, m_, k, b, var
   INTEGER :: space_dim, element_size_d, element_size_b
   INTEGER :: Ne_d, Ne_b, Np_d, Np_b, Np_bound, Ne_bound
   
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE  ::  bdata
   CHARACTER(LEN=18) :: bname
   CHARACTER(LEN=180) :: vcorr, var_names
   INTEGER  ::  mode
   !-------------------------------------------------------------------------------

   ! Opens TecPlot file
   IF (PRESENT(tec_mode)) THEN
     mode = tec_mode
   ELSE
     mode = TEC_STANDARD
   ENDIF

   WRITE(var_names,*) plot_data % var_names

   IF ((mode == TEC_STANDARD) .OR. (mode == TEC_OPEN)) THEN
   
     err = TecIni( plot_data%plot_name(1:plot_data%pb_l)//NULLCHAR,  &
     		   TRIM(var_names)//NULLCHAR,                        &
     		   plot_data%file_name(1:plot_data%fl_l)//NULLCHAR,  &
     		   scrtch_dir//NULLCHAR,			     &
     		   Debug,					     &
     		   VIsDouble					     )

     IF (err == -1) RETURN
     IF (mode == TEC_OPEN) RETURN
   
   ENDIF


   ! Initialization
   IF ((mode==TEC_STANDARD) .OR. (mode==TEC_WRITE)) THEN

     ! Number of space dimensions
     space_dim = SIZE(plot_data%rr,1)
   
     IF(space_dim > 1) THEN

     	! Total number of domain and boundary nodes
     	Np_d = SIZE(plot_data%rr,2)
     	Np_b = SIZE(plot_data%jd_jb)
   
     	! Total number of domain and boundary elements
     	Ne_d = SIZE(plot_data%j_m_d)
     	Ne_b = SIZE(plot_data%j_m_b)
   
     	! Maximum number of nodes per element
     	element_size_d = 0  ! Domain
     	DO m = 1, Ne_d
     	   element_size_d = MAX( element_size_d, & 
     				 SIZE(plot_data%j_m_d(m)%vec) )
     	ENDDO
     	element_size_b = 0  ! Boundary
     	DO m = 1, Ne_b
     	   element_size_b = MAX( element_size_b, & 
     				 SIZE(plot_data%j_m_b(m)%vec) )
     	ENDDO
	
     ELSE

     	Np_d = SIZE(plot_data%rr,2)
     	Np_b = 0
     	
     	Ne_d = 0
     	Ne_b = 0
     	
     	element_size_d = 0
     	element_size_b = 0
     	
     ENDIF



     ! Writes the plot file
     SELECT CASE (space_dim)

     	CASE(1)  ! One-dimensional plot
     	
     	   ! Domain
     	   ! Initialize the TecPlot zone
     	   err = TecZne( 'Domain'//NULLCHAR, &
     			 Np_d, 1, 1,	     &
     			 'BLOCK'//NULLCHAR,  &
     			 NULLCHAR	     )
     	   
     	   IF (err == -1) RETURN

     	   ! Writes the coordinates
     	   err = TecDat( Np_d,  	    &
     			 plot_data%rr(1,:), &
     			 DIsDouble	    )
     	      
     	   IF (err == -1) RETURN

     	   
     	   ! Writes the variables
     	   DO var = 1, SIZE(plot_data%variables,1)
     	      
     	      err = TecDat( Np_d,			&
     			    plot_data%variables(var,:), &
     			    DIsDouble			)
     	      
     	      IF (err == -1) RETURN
     	       
     	   ENDDO
     	   
     	   ! No boundary data is written


     	CASE(2)  ! Two-dimensional plot
     	
     	   ! Domain
     	   ! Initialize the TecPlot zone
     	   err = TecZne( 'Domain'//NULLCHAR,		 &
     			 Np_d, Ne_d, element_size_d - 3, &
     			 'FEBLOCK'//NULLCHAR,		 &
     			 NULLCHAR			 )
     	   
     	   IF (err == -1) RETURN

     	   ! Writes the coordinates
     	   DO k = 1, space_dim
     	   
     	      err = TecDat( Np_d,	       &
     			    plot_data%rr(k,:), &
     			    DIsDouble	       )
     	      
     	      IF (err == -1) RETURN
     	   
     	   ENDDO

     	   ! Writes the variables
     	   DO var = 1, SIZE(plot_data%variables,1)

     	      err = TecDat( Np_d,			&
     			    plot_data%variables(var,:), &
     			    DIsDouble			)
     	      
     	      IF (err == -1) RETURN
     	       
     	   ENDDO
     	
     	   ! Computes and writes the connectivity list
     	   ALLOCATE (jm(element_size_d, Ne_d))
     	   
     	   DO m = 1, Ne_d
     	      
     	      DO i_ = 1, SIZE(plot_data%j_m_d(m)%vec);  
     		 i  = plot_data%j_m_d(m)%vec(i_) 
     		 jm(i_,m) = i
     	      ENDDO
     	      
     	      ! If the elements has fewer nodes than the 
     	      ! maximum allowed, fills the blanks with node 1
     	      DO i_ = SIZE(plot_data%j_m_d(m)%vec) + 1, element_size_d
     		 i  = plot_data%j_m_d(m)%vec(1)
     		 jm(i_,m) = i
     	      ENDDO
     	   
     	   ENDDO
     	   
     	   err = TecNod(jm)

     	   IF (err == -1) RETURN
     	   
     	   DEALLOCATE (jm)
     	   
     	   
     	   ! Boundary
     	   ! For each boundary
     	   DO b = 1, SIZE(plot_data%p_bound)
     	      
     	      Np_bound = SIZE(plot_data%p_bound(b)%vec)
     	      ALLOCATE( bdata(Np_bound) )
     	   
     	      ! Initialize the TecPlot zone
     	      SELECT CASE (plot_data%boundary_types(b))

     		 CASE(0)   ! Numerical boundary (no bc)
     		 WRITE(bname,20) b
     	      
     		 CASE(1)   ! Solid wall
     		 WRITE(bname,10) b
     	      
     		 CASE(2)   ! Numerical boundary (Riemann)
     		 WRITE(bname,20) b

     		 CASE DEFAULT
     		 WRITE(bname,30) b
     	      
     	      END SELECT
     	      
     	      10 FORMAT('Solid boundary ',i3)
     	      20 FORMAT('Numer boundary ',i3)
     	      30 FORMAT('      boundary ',i3)
     			 
     	      err = TecZne( bname//NULLCHAR, &
     			    Np_bound, 1, 1,	    &
     			    'BLOCK'//NULLCHAR,  &
     			    NULLCHAR		)
     	   
     	      IF (err == -1) RETURN

     	      ! Writes the coordinates
     	      DO k = 1, space_dim
     	      
     		 bdata = plot_data%rr(k, &
     			    plot_data%jd_jb(plot_data%p_bound(b)%vec))
     			    
     		 err = TecDat( Np_bound, &
     			       bdata,	 &
     			       DIsDouble )
     	      
     		 IF (err == -1) RETURN
     	   
     	      ENDDO
     			  
     	      ! Writes the variables
     	      DO var = 1, SIZE(plot_data%variables,1)

     		 bdata  = plot_data%variables(var, &
     			    plot_data%jd_jb(plot_data%p_bound(b)%vec))

     		 err = TecDat( Np_bound, &
     			       bdata,	 &
     			       DIsDouble )

     		 IF (err == -1) RETURN

     	      ENDDO

     	      DEALLOCATE (bdata)
     	   
     	   ENDDO
     	

     	CASE(3)  ! Three-dimensional plot
     	
     	   ! Domain
     	   WRITE(*,*) ' Domain '

     	   ! Initialize the TecPlot zone
     	   err = TecZne( 'Domain'//NULLCHAR,		 &
     			 Np_d, Ne_d, 3, 		 &
     			 'FEBLOCK'//NULLCHAR,		 &
     			 NULLCHAR			 )
     	   
     	   IF (err == -1) RETURN

     	   ! Writes the coordinates
     	   DO k = 1, space_dim
     	   
     	      err = TecDat( Np_d,	       &
     			    plot_data%rr(k,:), &
     			    DIsDouble	       )
     	      
     	      IF (err == -1) RETURN
     	   
     	   ENDDO

     	   ! Writes the variables
     	   DO var = 1, SIZE(plot_data%variables,1)
     	      err = TecDat( Np_d,			&
     			    plot_data%variables(var,:), &
     			    DIsDouble			)
     	      
     	      IF (err == -1)  RETURN
     	   ENDDO
     	
     	   ! Computes and writes the connectivity list  
     	   ! ALLOCATE (jm(element_size_d, Ne_d))
     	   ALLOCATE (jm(8, Ne_d))
     	   
     	   DO m = 1, Ne_d
     	      
     	      DO i_ = 1, SIZE(plot_data%j_m_d(m)%vec);  
     		 i  = plot_data%j_m_d(m)%vec(i_)
     		 jm(i_,m) = i
     	      ENDDO
     	      
     	      ! If the elements has fewer nodes than the 
     	      ! maximum allowed, fills the blanks with node 1
     	      DO i_ = SIZE(plot_data%j_m_d(m)%vec) + 1, 8
     		 i  = plot_data%j_m_d(m)%vec(1)
     		 
     		 jm(i_,m) = i
     	      ENDDO
     	   
     	   ENDDO
     	   
     	   err = TecNod(jm)

     	   IF (err == -1) RETURN
     	   
     	   DEALLOCATE (jm)


     	   ! Boundary
     	   vcorr = ""
     	   CALL variable_correspondance(3+SIZE(plot_data%variables,1),vcorr)
     	   
     	   ! For each boundary
     	   DO b = 1, SIZE(plot_data%p_bound)

     	      WRITE(*,*) ' Boundary ', b

     	      Np_bound = SIZE(plot_data%p_bound(b)%vec)
     	      Ne_bound = SIZE(plot_data%m_bound(b)%vec)
	      
     	      ! ALLOCATE (bdata(Np_bound))
     	   
     	      ! Initialize the TecPlot zone
     	      SELECT CASE(plot_data%boundary_types(b))

     		 CASE(0)   ! Numerical boundary (no bc)
                 WRITE(bname,20) b
     	      
     		 CASE(1)   ! Solid wall
                 WRITE(bname,13) b
     	      
     		 CASE(2)   ! Numerical boundary (Riemann)
                 WRITE(bname,23) b

     		 CASE DEFAULT
                 WRITE(bname,33) b
     	      
     	      END SELECT
     	      
     	      13 FORMAT('Solid boundary ',i3)
     	      23 FORMAT('Numer boundary ',i3)
     	      33 FORMAT('      boundary ',i3)
     	      
     	      err = TecZne( bname//NULLCHAR,	  &
     			    Np_d, Ne_bound, 1,    &
     			    'FEBLOCK'//NULLCHAR,  &
     			    TRIM(vcorr)//NULLCHAR )
     	   
     	      IF (err == -1) RETURN

     	      ! Computes and writes the connectivity list   
     	      ALLOCATE ( jm(4, Ne_bound) )
     	    
     	      DO m_ = 1, Ne_bound
     	      
     		 m = plot_data%m_bound(b)%vec(m_)
     		 DO i_ = 1, SIZE(plot_data%j_m_b(m)%vec)
     		    i  = plot_data%j_m_b(m)%vec(i_) 
     		    jm(i_, m_) = plot_data%jd_jb(i)
     		 ENDDO

     		 ! If the elements has fewer nodes than the 
     		 ! maximum allowed, fills the blanks with node 1
     		 DO i_ = SIZE(plot_data%j_m_b(m)%vec) + 1, 4
   		    i  = plot_data%j_m_b(m)%vec(1)
     		    jm(i_,m_) = plot_data%jd_jb(i)
     		 ENDDO

     	      ENDDO

     	      err = TecNod(jm)

     	      IF (err == -1) RETURN

     	      DEALLOCATE (jm)

     	      ENDDO


     	CASE DEFAULT

     	   WRITE(*,*) ''      
     	   WRITE(*,*) ' ERROR. WRITE_PLOT_FILE_TECPLOT'
     	   WRITE(*,*) ' space dimension ', space_dim, ' not implemented'
     	   WRITE(*,*) ''
     	   STOP

     END SELECT 

   ENDIF

   IF (mode == TEC_WRITE) RETURN
  
   ! Closes TecPlot file
   IF ((mode==TEC_STANDARD) .OR. (mode==TEC_CLOSE))  err = TecEnd()

   END FUNCTION write_plot_file_tecplot





   SUBROUTINE variable_correspondance(Nvar, var) 
   !---------------------------------------------------------------------------------------   
   IMPLICIT NONE

   INTEGER  :: Nvar
   CHARACTER(LEN=64)  :: var
   CHARACTER(LEN=1)  :: dummy_1
   CHARACTER(LEN=2)  :: dummy_2
   CHARACTER(LEN=3)  :: dummy_3
   
   INTEGER  :: i, pos
   !---------------------------------------------------------------------------------------
         
   pos = 1
   DO i = 1, Nvar

      SELECT CASE (i)
   	 CASE(1:9)
   	    WRITE(dummy_1,10) i
   	    var = TRIM(var)//dummy_1
   	 CASE(10:99)
   	    WRITE(dummy_2,100) i
   	    var = TRIM(var)//dummy_2
   	 CASE(100:999)
   	    WRITE(dummy_3,1000) i
   	    var = TRIM(var)//dummy_3
      END SELECT
   
      IF (i < Nvar) var = TRIM(var)//','
    
   ENDDO
   
   10	FORMAT(i1.1)
   100  FORMAT(i2.2)
   1000 FORMAT(i2.2)
   
   END SUBROUTINE  variable_correspondance


   END MODULE  plot_procedures
