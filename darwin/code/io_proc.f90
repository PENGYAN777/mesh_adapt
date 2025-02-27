MODULE io_proc

 USE structures
 USE dynamic_vector
 USE lib_geo
 USE mesh_structure
 USE metric_coefficients
 USE nodes
 USE np_topology_gen 
 USE node_pair_structure
 USE grid_utils


 CONTAINS

 
   FUNCTION Load_Adaption_Param()   RESULT(params)
   !--------------------------------------------------------------------------------------
   USE thermodynamic,  ONLY: mu_ref, T_ref, T_char, T_scaling, L_scaling, P_scaling
   
   IMPLICIT NONE
 
   INTEGER              :: i, j, nv
   TYPE(adaption_param) :: params
   !--------------------------------------------------------------------------------------

    OPEN (UNIT=11, FILE='darwin.cfg', STATUS='old')
 
      READ(11,*)
      READ(11,*) params % prob_name
      READ(11,*) params % grid_name
      READ(11,*) params % sol_fmt
      READ(11,*)
      READ(11,*)
      READ(11,*) params % number_adaptation_levels
      READ(11,*) params % entry_level      
      READ(11,*)
      READ(11,*)
      READ(11,*) params % ref_type
      READ(11,*) params % nbbox


      IF ( params % nbbox .ge. 1 ) THEN
        ALLOCATE (params % box(params % nbbox))

        DO i = 1, params % nbbox

          READ(11,'(i1)', ADVANCE='NO') nv
          ALLOCATE (params % box(i) % vec(2*nv))
          READ(11,*) params % box(i) % vec

        ENDDO
      END IF

      READ(11,*)
      READ(11,*)
 
      READ(11,*) params % N_estimators

      ALLOCATE (params % estimator(params % N_estimators))

      DO i = 1, params % N_estimators
        READ(11,*) params % estimator(i) % error_variable

        
        IF (params % estimator(i) % error_variable == 6) THEN
        
          OPEN(UNIT=12, FILE='solution.'//params % prob_name, STATUS='old')
            
            DO j = 1, 38
              READ(12,*)
            ENDDO
            
            READ(12,*) mu_ref, T_ref, T_char
            
            DO j = 1, 6
              READ(12,*)
            ENDDO
            
            READ(12,*) l_scaling
            READ(12,*) p_scaling            
            READ(12,*) T_scaling
            
          CLOSE(12)
        
        ENDIF

        READ(11,*) params % estimator(i) % estimator_function
        READ(11,*) params % estimator(i) % passages
      ENDDO


      READ(11,*)
      READ(11,*)
      READ(11,*) params % Nb_lines

      ALLOCATE (params % solid_bou(params % Nb_lines))

      DO i = 1, params % Nb_lines
        READ(11,*) params % solid_bou(i)
      ENDDO

      READ(11,*) params % bou_tre_typ
      READ(11,*) params % jou_mod_fml 

      IF (params % jou_mod_fml == 0) READ(11,*) params % jmf_const
      IF (params % jou_mod_fml == 1) READ(11,*) params % jmf_alpha

      READ(11,*) params % Poi_Coe

      READ(11,*)
      READ(11,*)
      READ(11,*) params % min_angle
      READ(11,*) params % max_number_QL
      READ(11,*) params % min_length
      READ(11,*) params % swap
      READ(11,*) params % smooth
      READ(11,*) params % l_its, params % relax_1, params % relax_2
      
    CLOSE (11)
   


   END FUNCTION Load_Adaption_Param
 




   FUNCTION Load_Grid( level, g_name, ext ) RESULT(loaded_grid)     
   !===================================================================================! 
   IMPLICIT NONE
   
   INTEGER,           INTENT(IN)  :: level
   CHARACTER(LEN=7) , INTENT(IN)  :: ext
   CHARACTER(LEN=30), INTENT(IN)  :: g_name
   
   
   TYPE(grid_type)                        :: loaded_grid                               
   INTEGER                                :: idf, i
   TYPE(D_I_V), DIMENSION(:), ALLOCATABLE :: m_j_d
   TYPE(D_I_V), DIMENSION(:), ALLOCATABLE :: m_j_b
   !===================================================================================! 

    idf = 10 
         
    OPEN (idf, FILE = 'nodes.'//TRIM(g_name)//'_'//TRIM(ext), STATUS='old') 
      CALL read_nodes(idf)
    CLOSE (idf) 
    
    OPEN (idf, FILE = 'grid.'//TRIM(g_name)//'_'//TRIM(ext), STATUS = 'old') 
      CALL read_mesh(idf) 
    CLOSE (idf)

    OPEN (idf, FILE = 'metrics.'//TRIM(g_name)//'_'//TRIM(ext), STATUS = 'old') 
      CALL read_metric(idf) 
    CLOSE (idf)

    OPEN (idf, FILE = 'node_pair.'//TRIM(g_name)//'_'//TRIM(ext), STATUS = 'old') 
      CALL read_node_pair(idf) 
    CLOSE (idf)        
    
    loaded_grid%adaptation_level = level
    loaded_grid%k_d  = k_d
    loaded_grid%Nb   = MAXVAL(bound_p)                   
    
    ! Nodes
    loaded_grid%Nj_d = Np_d;  loaded_grid%Nj_b = Np_b
    ALLOCATE( loaded_grid%rr(k_d,Np_d)  )
    ALLOCATE( loaded_grid%jd_jb(Np_b)   )
    ALLOCATE( loaded_grid%bound_p(Np_b) )
                              
    loaded_grid%rr      = rr
    loaded_grid%jd_jb   = jd_jb
    loaded_grid%bound_p = bound_p           
    
    CALL Invert_jd_jb(loaded_grid)
    
    ! Domain elements  
    loaded_grid % Nm_d = Ne_d
    loaded_grid % Nm_b = Ne_b                
    
    ALLOCATE( loaded_grid%m_j_d(Np_d)  )
    ALLOCATE( loaded_grid%j_m_d(Ne_d)  )
    ALLOCATE( loaded_grid%ma_m_d(Ne_d) )
   
    ! Inverting j_m_d
    ALLOCATE( m_j_d(size_DIV(j_m_d,3)) )
    m_j_d = invert_DIV( j_m_d )
   
    DO i = 1, Np_d        
      ALLOCATE( loaded_grid%m_j_d(i)%vec(SIZE(m_j_d(i)%vec)) )      
      loaded_grid%m_j_d(i)%vec = m_j_d(i)%vec   
    ENDDO    
    
    DO i = 1, Ne_d    
      
      ALLOCATE( loaded_grid%j_m_d(i)%vec(SIZE(j_m_d(i)%vec)) )
      loaded_grid%j_m_d(i)%vec = j_m_d(i)%vec
      
      ALLOCATE( loaded_grid%ma_m_d(i)%vec(SIZE(ma_m_d(i)%vec)) )
      loaded_grid%ma_m_d(i)%vec = ma_m_d(i)%vec
      
    ENDDO

    ALLOCATE( loaded_grid%ele_type_d(Ne_d) )
    loaded_grid%ele_type_d = ele_type_d
    
    ! Boundary elements    
    loaded_grid%Nm_b = Ne_b;  loaded_grid%Nm_b = Ne_b                
    ALLOCATE( loaded_grid%m_j_b(Np_b)  )
    ALLOCATE( loaded_grid%j_m_b(Ne_b)  )
    ALLOCATE( loaded_grid%ma_m_b(Ne_b) )
    
    ! Inverting m_j_d
    ALLOCATE( m_j_b(size_DIV(j_m_b,3)) )
    m_j_b = invert_DIV( j_m_b )

    DO i = 1, Np_b        
      ALLOCATE( loaded_grid%m_j_b(i)%vec(SIZE(m_j_b(i)%vec)) )      
      loaded_grid%m_j_b(i)%vec = m_j_b(i)%vec 
    ENDDO    
    
    DO i = 1, Ne_b    
      
      ALLOCATE( loaded_grid%j_m_b(i)%vec(SIZE(j_m_b(i)%vec)) )
      loaded_grid%j_m_b(i)%vec = j_m_b(i)%vec     
      
      ALLOCATE( loaded_grid%ma_m_b(i)%vec(SIZE(ma_m_b(i)%vec)) )
      loaded_grid%ma_m_b(i)%vec = ma_m_b(i)%vec      
      
    ENDDO

    ALLOCATE( loaded_grid%bound_m(Ne_b)    )
    loaded_grid%bound_m = bound_m
    
    ALLOCATE( loaded_grid%ele_type_b(Ne_b) )
    loaded_grid%ele_type_b = ele_type_b

    ! Node-Pairs          
    loaded_grid%Nc_d = Nc_fv
    loaded_grid%Nc_b = Nc_b
    
    ALLOCATE (loaded_grid%jcd_fem(4*k_d,Nc_d))
    loaded_grid%jcd_fem = j_c_d

    ALLOCATE (loaded_grid%jcb_fem(4*k_d,Nc_b))
    loaded_grid%jcb_fem = j_c_b
    
    ALLOCATE (loaded_grid%j_c(2,Nc_fv))
    loaded_grid%j_c = j_c_fv(1:2,:)
    
    ALLOCATE (loaded_grid%j_c_b(2,Nc_b))
    loaded_grid%j_c_b = j_c_b(1:2,:)
    
    ALLOCATE (loaded_grid%cd_cb(Nc_b))
    loaded_grid%cd_cb = cd_cb
    
    ALLOCATE (loaded_grid%bound_c(Nc_b))
    loaded_grid%bound_c = bound_c
        
    CALL Edge_Structure_Gen(loaded_grid, 'P')                        
    

    ! Metrics       
    ALLOCATE (loaded_grid%cell(SIZE(cell)))
    loaded_grid%cell = cell
    
    ALLOCATE (loaded_grid%eta(SIZE(eta,1),SIZE(eta,2))) 
    loaded_grid%eta = eta    

    ALLOCATE (loaded_grid%chi_b(SIZE(chi_b,1),SIZE(chi_b,2)))
    loaded_grid%chi_b = chi_b    
 
    ALLOCATE (loaded_grid%xi_bp(SIZE(xi_bp,1),SIZE(xi_bp,2)))
    loaded_grid%xi_bp = xi_bp

    ALLOCATE (loaded_grid%stiff_T(SIZE(stiff_T,1), SIZE(stiff_T,2), SIZE(stiff_T,3)))
    loaded_grid%stiff_T = stiff_T

!    ALLOCATE (loaded_grid%Kb_ijs(SIZE(Kb_ijs,1), SIZE(Kb_ijs,2), SIZE(Kb_ijs,3)))
!    loaded_grid%Kb_ijs = Kb_ijs

!    ALLOCATE (loaded_grid%Kb_ija(SIZE(Kb_ija,1), SIZE(Kb_ija,2), SIZE(Kb_ija,3)))
!    loaded_grid%Kb_ija = Kb_ija


    ! Loading Boundary Curve Coordinates and Geometry
    IF (ext == '0  ') THEN

      ALLOCATE (loaded_grid % curve_coordinates(loaded_grid % Nj_b))
      loaded_grid % curve_coordinates = Load_Curve_Coords(g_name)

      ALLOCATE (loaded_grid%line(loaded_grid%Nb))
      loaded_grid%line = Load_Curves_Geometry(g_name)

    ENDIF
                            
 END FUNCTION Load_grid
 
 
 
 
 
 FUNCTION Load_Curve_Coords(g_name)   RESULT (curve_coordinates)
 !=================================================================================
  IMPLICIT NONE
 
  INTEGER           :: j, dum, ns
  CHARACTER(LEN=30) :: g_name
 
  REAL(KIND=8), DIMENSION(:), POINTER :: curve_coordinates
 !=================================================================================
    
  OPEN (UNIT=10, FILE='abscissa.'//TRIM(g_name)//'_0', STATUS='old')

    READ(10,*)
    READ(10,*)
    READ(10,*) ns
    READ(10,*)
    READ(10,*)
    READ(10,*)

    ALLOCATE (curve_coordinates(ns))

    DO j = 1, ns
      READ(10,*) dum, curve_coordinates(j)
    ENDDO

  CLOSE (10)

 END FUNCTION Load_Curve_Coords
 
 
 
 
 
 
 FUNCTION Load_Curves_Geometry(g_name) RESULT(line)
 !=================================================================================!
  IMPLICIT NONE
 
  INTEGER           :: j, k, b_nmr
  INTEGER           :: cc_dim, sp_dim
  CHARACTER(LEN=30) :: g_name
 
  TYPE(curve), DIMENSION(:), POINTER :: line
 !=================================================================================!
 
  k = 0

  OPEN(UNIT=11, FILE='geometry.'//TRIM(g_name)//'_0', STATUS='old')
 
    READ(11,*)
    READ(11,*) b_nmr
 
    ALLOCATE( line(b_nmr) )
 
1   k = k + 1

    READ(11,*)
    READ(11,*) line(k)%nd, line(k)%ns
    READ(11,*)
 
    sp_dim = line(k)%nd
    cc_dim = line(k)%ns
 
    ALLOCATE( line(k)%s(cc_dim), line(k)%h(cc_dim) )
 
    ALLOCATE( line(k)%x(sp_dim,cc_dim) )
    ALLOCATE( line(k)%xs(sp_dim,cc_dim) )
 
    DO j = 1, cc_dim
      READ(11,*) line(k)%s(j), line(k)%x(1,j),  line(k)%x(2,j),  &
                               line(k)%xs(1,j), line(k)%xs(2,j)
    ENDDO
 
    IF (k < b_nmr)  GOTO 1
 
  CLOSE (11)
 
 END FUNCTION Load_Curves_Geometry
 
 
 
 
 
 FUNCTION Load_Solution(sol_fmt, Sd, Np, Nb, p_name)   RESULT(solution)  
 !===================================================================================!   
  IMPLICIT NONE
 
  INTEGER,           INTENT(IN) :: sol_fmt, Sd, Np, Nb
  CHARACTER(LEN=30), INTENT(IN) :: p_name
  CHARACTER(LEN=80) :: aaa
  
  REAL(KIND=8)  ::  x, y
  INTEGER :: i, idx
  
  TYPE(solution_type) :: solution
 !===================================================================================! 

  solution % k_d  = Sd
  solution % Nj_d = Np
  solution % sol_fmt = sol_fmt
  
  OPEN (UNIT=11, FILE='solution.'//TRIM(p_name), STATUS='old')

    SELECT CASE (sol_fmt)

      CASE (PIG_EU)

           PRINT*, ' Reading solution'
           PRINT*, ' - equations:  Euler'
           PRINT*, ' - thermodynamics:  Perfect Ideal Gas'

           DO i = 1, 65+Nb   !npcode 65   !solver 51
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution%ww(Sd+2,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)  
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
           ENDDO


      CASE (PIG_NSL)

           PRINT*, ' Reading solution'
           PRINT*, ' - equations:  Navier-Stokes Laminar'
           PRINT*, ' - thermodynamics:  Perfect Ideal Gas'

           DO i = 1, 68+Nb
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution%ww(Sd+2,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)  
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
           ENDDO


      CASE (PIG_NST)

           PRINT*, ' Reading solution'
           PRINT*, ' - equations:  RANS'
           PRINT*, ' - thermodynamics:  Perfect Ideal Gas'

           DO i = 1, 55 + Nb !59+Nb
             READ(11,*) aaa
           ENDDO

    
           ALLOCATE (solution%ww(Sd+3,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)  
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
             READ(11,*) solution % ww(Sd+3,i)  ! Turbulent viscosity
           ENDDO


      CASE (VdW_EU)

           PRINT*, ' Reading solution'
           PRINT*, ' - equations:  Euler'
           PRINT*, ' - thermodynamics:  Van der Waals'

           DO i = 1, 52+Nb   !npCode 56   solver 52 
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution % ww(Sd+2,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
           ENDDO


      CASE (VdW_NSL)

           PRINT*, ' Reading solution'
           PRINT*, ' - equations:  Navier-Stokes Laminar'
           PRINT*, ' - thermodynamics:  Van der Waals'

           DO i = 1, 52+Nb   !npCode 56   solver 52 
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution % ww(Sd+2,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
           ENDDO


      CASE (VdW_NST)

           PRINT*, ' Reading solution'
           PRINT*, ' - equations:  RANS'
           PRINT*, ' - thermodynamics:  Van der Waals'

           DO i = 1, 61+Nb
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution % ww(Sd+3,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)  
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
             READ(11,*) solution % ww(Sd+3,i)   ! Turbulent viscosity
           ENDDO

      CASE (FENSAP)
      
           PRINT*, '   Reading FENSAP solution'
           PRINT*, '   - equations:  RANS'
           PRINT*, '   - thermodynamics:  Perfect Ideal Gas'

           DO i = 1, 68 + Nb
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution % ww(5,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)  
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
             READ(11,*) solution % ww(Sd+3,i)   ! Turbulent viscosity
           ENDDO
      
      CASE (FENSAPDROP)
      
           PRINT*, '   Reading FENSAP+DROP solution'
           PRINT*, '   - equations:  RANS'
           PRINT*, '   - thermodynamics:  Perfect Ideal Gas'

           DO i = 1, 68 + Nb
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution % ww(7,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)  
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
             READ(11,*) solution % ww(Sd+3,i)   ! Turbulent viscosity
             READ(11,*) solution % ww(Sd+4,i)   ! Liquid Water Content
             READ(11,*) solution % ww(Sd+5,i)   ! Collection efficiency
           ENDDO
      
      CASE (FENSAPDROPICE)

           PRINT*, '   Reading FENSAP+DROP+ICE solution'
           PRINT*, '   - equations:  RANS'
           PRINT*, '   - thermodynamics:  Perfect Ideal Gas'

           DO i = 1, 68 + Nb
             READ(11,*) aaa
           ENDDO
    
           ALLOCATE (solution % ww(9,Np))
              
           DO i = 1, Np    
             READ(11,*) idx, solution % ww(1,i)  
             READ(11,*) solution % ww(2:Sd+1,i)
             READ(11,*) solution % ww(Sd+2,i)
             READ(11,*) solution % ww(Sd+3,i)   ! Turbulent viscosity
             READ(11,*) solution % ww(Sd+4,i)   ! Liquid Water Content
             READ(11,*) solution % ww(Sd+5,i)   ! Collection efficiency
             READ(11,*) solution % ww(Sd+6,i)   ! Liquid Water Content
             READ(11,*) solution % ww(Sd+7,i)   ! Collection efficiency
           ENDDO

      CASE (SU2_EU, SU2_NSL, SU2_NSSA, SU2_NSSST)

           PRINT*, '   Reading SU2 solution'
           SELECT CASE (sol_fmt)
              CASE (SU2_EU) 
                 PRINT*, '   - equations:  Euler'
              CASE (SU2_NSL) 
                 PRINT*, '   - equations:  NS laminar'
              CASE (SU2_NSSA)
                 PRINT*, '   - equations:  RANS Spalart-Allmaras model'
              CASE (SU2_NSSST)
                 PRINT*, '   - equations:  RANS SST model'
           END SELECT
           PRINT*, '   - thermodynamics:  not used'

           SELECT CASE (sol_fmt)
              CASE (SU2_EU) 
                 ALLOCATE (solution % ww(6+Sd,Np))
                 !ALLOCATE (solution % ww(10,Np))
              CASE (SU2_NSL) 
                 ALLOCATE (solution % ww(11+Sd,Np))
                 !ALLOCATE (solution % ww(14,Np))
              CASE (SU2_NSSA)
                 ALLOCATE (solution % ww(13+Sd,Np))
                 !ALLOCATE (solution % ww(15,Np))
              CASE (SU2_NSSST)
                 ALLOCATE (solution % ww(14+Sd,Np))
                 !ALLOCATE (solution % ww(15,Np))
           END SELECT

           READ(11,*) aaa
              
           DO i = 1, Np    
               READ(11,*) idx, x, y, solution % ww(1:,i)

               ! SELECT CASE (sol_fmt)
               !    CASE (SU2_EU) 
               !       ! READ(11,*) idx, x, y,    &      ! ID and x,y coordinates 
               !       !    solution % ww(1:Sd+2,i), &   ! Conservative variables
               !       !    solution % ww(Sd+5:10,i)     ! P, T, CP, Mach
               !       ! solution % ww(Sd+3:Sd+4,i) = 0.d0  ! Dummy
               !       READ(11,*) idx, x, y, solution % ww(1:,i)
               !    CASE (SU2_NSL) 
               !       ! READ(11,*) idx, x, y,    &      ! ID and x,y coordinates 
               !       !    solution % ww(1:Sd+2,i), &   ! Conservative variables
               !       !    solution % ww(Sd+5:14,i)     ! P, T, CP, Mach, 
               !       !                                 ! mu, Cf, qw, yp
               !       ! solution % ww(Sd+3:Sd+4,i) = 0.d0  ! Dummy
               !       READ(11,*) idx, x, y, solution % ww(1:,i)
               !    CASE (SU2_NSSA)
               !       ! READ(11,*) idx, x, y,    &   ! ID and x,y coordinates 
               !       !    solution % ww(1:Sd+2,i), &   ! Conservative variables
               !       !    solution % ww(Sd+3,i),   &   ! SA eddy viscosity 
               !       !    solution % ww(Sd+5:15,i)     ! P, T, CP, Mach, 
               !       !                                 ! mu, Cf, qw, yp, mu_turb
               !       ! solution % ww(Sd+4,i) = 0.d0       ! Dummy
               !       READ(11,*) idx, x, y, solution % ww(1:,i)
               !    CASE (SU2_NSSST)
               !       ! READ(11,*) idx, x, y,    &   ! ID and x,y coordinates 
               !       !    solution % ww(1:Sd+2,i), &   ! Conservative variables
               !       !    solution % ww(Sd+3,i),   &   ! SA eddy viscosity 
               !       !    solution % ww(Sd+5:15,i)     ! P, T, CP, Mach, 
               !       !                                 ! mu, Cf, qw, yp, mu_turb
               !       ! solution % ww(Sd+4,i) = 0.d0       ! Dummy
               !       READ(11,*) idx, x, y, solution % ww(1:,i)
               ! END SELECT

           ENDDO

    END SELECT
    
  CLOSE(11)
    
 END FUNCTION Load_Solution
 
 
 
 
 
 FUNCTION Load_Adaption_History(level)   RESULT(history)
 !===================================================================================! 
 IMPLICIT NONE
 
 INTEGER, INTENT(IN) :: level
 
 INTEGER                            :: i, j
 INTEGER                            :: idx
 INTEGER, DIMENSION(:), ALLOCATABLE :: dims  
 CHARACTER(LEN=15)                  :: string
 
 TYPE(E_R_P), DIMENSION(0:level) :: history 
 !===================================================================================! 
  
  OPEN( UNIT=11, FILE='ADAPTION_HISTORY.dat', STATUS='old' )
  
   DO i = 0, level-1
 
    READ(11,*)
    READ(11,*) string, history(i) % N_refined_edges
 
    ALLOCATE( history(i) % inserted_nodes(history(i) % N_refined_edges) )
    ALLOCATE( history(i) % edge_nodes( 2, history(i) % N_refined_edges) )
    ALLOCATE( history(i) % adjacent_nodes(history(i) % N_refined_edges) )
    ALLOCATE( dims(history(i) % N_refined_edges) )
 
    READ(11,*) string, history(i) % inserted_nodes
 
    DO j = 1, history(i) % N_refined_edges
      READ(11,*) idx, string, history(i) % edge_nodes(:,j)
    ENDDO
 
    READ(11,*) dims
 
    DO j = 1, history(i) % N_refined_edges

      ALLOCATE (history(i) % adjacent_nodes(j) % vec(dims(j)))

      READ(11,*) idx, string, history(i) % adjacent_nodes(j)%vec

    ENDDO
 
    DEALLOCATE(dims)
 
   ENDDO
  
  CLOSE(11)

 END FUNCTION Load_Adaption_History





 SUBROUTINE Write_Adaption_History(history, level, entry_level)
 
 !===================================================================================! 
  IMPLICIT NONE
 
  INTEGER                                      :: i, j, level, entry_level
  INTEGER,     DIMENSION(:),       ALLOCATABLE :: dimensions
  
  TYPE(E_R_P), DIMENSION(0:entry_level), INTENT(IN)  :: history
 !===================================================================================! 
  
  
  OPEN(UNIT=11, FILE='ADAPTION_HISTORY.dat', STATUS='unknown')
    
    DO i = 0, level-1
      
      WRITE(11,*) 'Level', i
      WRITE(11,*) 'N_refined_edges:', history(i)%N_refined_edges    
      WRITE(11,*) 'Inserted_nodes:', history(i)%inserted_nodes 
      
      DO j = 1, history(i)%N_refined_edges
        WRITE(11,*) history(i)%inserted_nodes(j), 'Edge_nodes:', history(i)%edge_nodes(:,j)      
      ENDDO      
      
      IF ( ALLOCATED(dimensions) ) DEALLOCATE(dimensions)
      ALLOCATE( dimensions(history(i)%N_refined_edges) )
      
      DO j = 1, history(i)%N_refined_edges
        dimensions(j) = SIZE(history(i)%adjacent_nodes(j)%vec) 
      ENDDO      
      
      WRITE(11,*) dimensions
      
      DO j = 1, history(i)%N_refined_edges
        WRITE(11,*) history(i)%inserted_nodes(j), &
                    'Adjacent_nodes:', history(i)%adjacent_nodes(j)%vec      
      ENDDO      
      
    ENDDO  
    
  CLOSE(11)

 END SUBROUTINE Write_Adaption_History




 SUBROUTINE Write_Grid(grid, g_name, ext)
 !==================================================================================!
  IMPLICIT NONE
  
  TYPE(grid_type),   INTENT(IN) :: grid
  CHARACTER(LEN=7),  INTENT(IN) :: ext
  CHARACTER(LEN=30), INTENT(IN) :: g_name
  
  INTEGER, PARAMETER :: idf=11
  INTEGER            :: i
  INTEGER            :: m   ! Generic element index      
  INTEGER            :: f_  ! Generic element's face index in local coordinates  
  INTEGER            :: j_  ! Generic node index in local coordinates  
  
 !==================================================================================!
  
  OPEN( UNIT=idf, FILE='nodes.'//TRIM(g_name)//'_'//TRIM(ext), FORM='formatted', STATUS='unknown' )
  
    WRITE (idf,1000)
    WRITE (idf,1010) g_name

    WRITE (idf,1000)
    WRITE (idf,1020)
    WRITE (idf,1021) grid%k_d, grid%Nj_d, grid%Nj_b 

    WRITE (idf,1000)
    WRITE (idf,1025)

    WRITE (idf,1000)
    WRITE (idf,1048)
    WRITE (idf,1050)
    WRITE (idf,1051)

    ! Volume nodes coordinates
    ! ------------------------
    DO i=1,grid%Nj_d
     WRITE (idf,1052) i
     WRITE (idf,1053) grid%rr(:,i)
    ENDDO

    WRITE (idf,1000)
    WRITE (idf,1055)

    WRITE (idf,1000)
    WRITE (idf,1048)     
    WRITE (idf,1075)

    ! Boundary nodes - domain nodes connectivity
    ! ------------------------------------------
    DO i=1,grid%Nj_b
     WRITE (idf,1076) i, grid%jd_jb(i), grid%bound_p(i)
    ENDDO
  
  CLOSE(idf)    

1000  FORMAT('###########################################################################')
1010  FORMAT('#    NAME:      ',a15,'                                           #')
1020  FORMAT('#         ND        NP_D        NP_B                                      #')
1021  FORMAT(5i12)
1025  FORMAT('#  **********  DOMAIN  **********                                         #')
1048  FORMAT('#  ++++   NODES   ++++                                                    #')
1050  FORMAT('#        IDX                                                              #')
1051  FORMAT('#         RR                                                              #')
1052  FORMAT(i12)
1053  FORMAT(3e18.9)
1055  FORMAT('#  **********  BOUNDARY  **********                                       #')
1075  FORMAT('#        IDX       JD_JB       BOUND                                      #')
1076  FORMAT(3i12) 


 OPEN (UNIT=idf, FILE='grid.'//TRIM(g_name)//'_'//TRIM(ext), FORM='formatted', STATUS='unknown')  

   WRITE (idf,2000)
   WRITE (idf,2010) g_name(1:LEN_TRIM(g_name))

   WRITE (idf,2000)
   WRITE (idf,2020)
   WRITE (idf,2021)  grid % Nm_d, grid % Nm_b

   WRITE (idf,2000)
   WRITE (idf,2025)

   WRITE (idf,2000)
   WRITE (idf,2040)
   WRITE (idf,2041) 
   WRITE (idf,2042) 
   WRITE (idf,2043) 

   ! Element of the domain
   ! ---------------------
   DO m=1,grid%Nm_d
      WRITE (idf,2046) m, grid%ele_type_d(m)
      DO j_ = 1, SIZE( grid%j_m_d(m)%vec ) 
         WRITE (UNIT=idf,FMT=2047,ADVANCE='NO') grid%j_m_d(m)%vec(j_)
      ENDDO
      WRITE(idf,*)
      DO f_ = 1, SIZE( grid%ma_m_d(m)%vec ) 
         WRITE (UNIT=idf,FMT=2047,ADVANCE='NO') grid%ma_m_d(m)%vec(f_)
      ENDDO
      WRITE(idf,*)
   ENDDO


   WRITE (idf,2000)
   WRITE (idf,2055)

   WRITE (idf,2000)
   WRITE (idf,2040)     
   WRITE (idf,2065)
   WRITE (idf,2042) 
   WRITE (idf,2043) 

   ! Element of the boundary
   ! -----------------------
   DO m=1,grid%Nm_b
      WRITE (idf,2066) m, grid%ele_type_b(m), grid%bound_m(m)
      DO j_ = 1, SIZE( grid%j_m_b(m)%vec ) 
         WRITE (UNIT=idf,FMT=2047,ADVANCE='NO')  grid%j_m_b(m)%vec(j_)
      ENDDO
      WRITE(idf,*)
      DO f_ = 1, SIZE( grid%ma_m_b(m)%vec ) 
         WRITE (UNIT=idf,FMT=2047,ADVANCE='NO') grid%ma_m_b(m)%vec(f_)
      ENDDO
      WRITE(idf,*)
   ENDDO
   
 CLOSE (idf)  

2000  FORMAT('###########################################################################')
2010  FORMAT('#    NAME:      ',a15,'                                           #')
2020  FORMAT('#       NE_D        NE_B                                                  #')
2021  FORMAT(2i12)
2025  FORMAT('#  **********  DOMAIN  **********                                         #')
2040  FORMAT('#  ++++  ELEMENTS ++++                                                    #')
2041  FORMAT('#        IDX        TYPE                                                  #')
2042  FORMAT('#        J_M                                                              #')
2043  FORMAT('#       MA_M                                                              #')
2046  FORMAT(2i12)
2047  FORMAT(i12)
2055  FORMAT('#  **********  BOUNDARY  **********                                       #')
2065  FORMAT('#        IDX        TYPE       BOUND                                      #')
2066  FORMAT(3i12)

  CALL Write_Curve_Coordinates(g_name, ext, grid % curve_coordinates,  &
                               grid % bound_p)

 
 END SUBROUTINE Write_Grid





 SUBROUTINE Write_Curve_Coordinates(g_name, ext, curve_coordinates, b_p)
 !=================================================================================!
 IMPLICIT NONE

 INTEGER           :: i
 CHARACTER(LEN=3)  :: ext
 CHARACTER(LEN=30) :: g_name
 REAL(KIND=8), DIMENSION(:) :: curve_coordinates
 INTEGER,      DIMENSION(:), INTENT(IN) :: b_p
 !=================================================================================!
    
  OPEN ( UNIT=12, FILE='abscissa.'//TRIM(g_name)//'_'//TRIM(ext), &
         FORM='formatted', STATUS='unknown' )    
   
   WRITE(12,1108) TRIM(g_name)//'_'//TRIM(ext)
   WRITE(12,1109)
   WRITE(12,1110) SIZE(curve_coordinates)       
   WRITE(12,1111)
   WRITE(12,1112)
   WRITE(12,1113)

   DO i = 1, SIZE(curve_coordinates)
     WRITE(12,1114) i, curve_coordinates(i), b_p(i)
   ENDDO

1108 FORMAT(' Grid Name: ',a64)
1109 FORMAT(' Number of coordinates')
1110 FORMAT(8x,i6)    
1111 FORMAT('')
1112 FORMAT('      Node Index          Curve Coordinate          Boundary')
1113 FORMAT(' (boundary notation)')
1114 FORMAT(4x, i7, 15x, e15.7, 15x, i3)

    
  CLOSE(12) 
  
 END SUBROUTINE Write_Curve_Coordinates







 SUBROUTINE Write_Edge_Error(estimator)
 !============================================================================
 IMPLICIT NONE

 TYPE(estimator_type), DIMENSION(:), INTENT(IN) :: estimator 

 INTEGER :: i, j
 !============================================================================
  
  OPEN(UNIT=11, FILE='EDGE_ERRORS.dat', FORM='formatted', STATUS='unknown')

   WRITE(11,1000) 
   WRITE(11,1001)
   
   DO i = 1, SIZE(estimator)

    WRITE(11,1002) estimator(i)%error_variable, estimator(i)%estimator_function
    WRITE(11,1009) SIZE(estimator(i)%errors,2)
    WRITE(11,1010) SIZE(estimator(i)%errors,1)
    WRITE(11,1003)
    WRITE(11,1004)
    WRITE(11,1005)
 
    DO j = 1, SIZE(estimator(i)%errors,2)
      IF (estimator(i)%estimator_function == 1) WRITE(11,1006) j, estimator(i)%errors(1,j)
      IF (estimator(i)%estimator_function == 2) WRITE(11,1007) j, estimator(i)%errors(1,j), estimator(i)%errors(2,j)
      IF (estimator(i)%estimator_function == 3) WRITE(11,1007) j, estimator(i)%errors(1,j), estimator(i)%errors(2,j)
      IF (estimator(i)%estimator_function == 4) WRITE(11,1007) j, estimator(i)%errors(1,j)
      IF (estimator(i)%estimator_function == 5) WRITE(11,1007) j, estimator(i)%errors(1,j)
    ENDDO
 
    WRITE(11,1008)

   ENDDO
     
  CLOSE(11) 
 
  1000 FORMAT(' ------------------------------------------------------------------')
  1001 FORMAT(' Edge Error ')
  1002 FORMAT(' Estimator (variable,function): ', i3,',', i3)
  1009 FORMAT(' Edges: ', i7)
  1010 FORMAT(' Error Values: ', i3)
  1003 FORMAT(' ------------------------------------------------------------------')
  1004 FORMAT('    Edge Index             Error(1)                Error(2)        ')
  1005 FORMAT(' (domain notation)                                                 ')
  1006 FORMAT( 3x, i7, 10x, e15.7 ) 
  1007 FORMAT( 3x, i7, 10x, e15.7, 10x, e15.7 ) 
  1008 FORMAT(' ------------------------------------------------------------------')
  
  CLOSE(11)

 END SUBROUTINE Write_Edge_Error     





 SUBROUTINE Write_Interpolated_Solution(grid, sol_fmt, Nb, solution, problem_name)
 !===================================================================================!   
  IMPLICIT NONE

  TYPE(grid_type),     INTENT(IN) :: grid
 
  INTEGER, INTENT(IN) :: sol_fmt, Nb
  
  TYPE(solution_type), INTENT(IN) :: solution  
  
  CHARACTER(LEN=30),   INTENT(IN) :: problem_name

  CHARACTER(LEN=400) :: head_sol
  
  INTEGER :: i, p  
 !===================================================================================! 
 
  p = solution % k_d + 2
  !OPEN (UNIT=11, FILE='solution.'//TRIM(problem_name),                  STATUS='old')
  OPEN (UNIT=12, FILE='solution.'//TRIM(problem_name)//'.interpolated', FORM='formatted', STATUS='unknown')

    
    SELECT CASE (sol_fmt)


      CASE (PIG_EU)

           CALL file_header(Nb, solution%k_d, sol_fmt)
    
           DO i = 1, SIZE(solution%ww,2)
             WRITE(12,1085) i, solution % ww(1,i)              
             WRITE(12,1086)    solution % ww(2:p-1, i)        
             WRITE(12,1087)    solution % ww(p, i)
           ENDDO                                               
    

      CASE (PIG_NSL)
           
           CALL file_header(Nb, solution%k_d, sol_fmt)
           
           DO i = 1, SIZE(solution%ww,2)
             WRITE(12,1085) i, solution % ww(1,i)              
             WRITE(12,1086)    solution % ww(2:p-1, i)        
             WRITE(12,1087)    solution % ww(p, i)
           ENDDO                                               
    

      CASE (PIG_NST)
           
           !CALL file_header(Nb, solution%k_d, sol_fmt)
           
           DO i = 1, SIZE(solution%ww,2)
             WRITE(12,1085) i, solution % ww(1,i)              
             WRITE(12,1086)    solution % ww(2:p-1, i)        
             WRITE(12,1087)    solution % ww(p, i)
             WRITE(12,1087)    solution % ww(p+1,i)
           ENDDO

     CASE (SU2_EU, SU2_NSL, SU2_NSSA, SU2_NSSST)

           OPEN (UNIT=101, FILE='solution.'//TRIM(problem_name), STATUS='old')
           READ(101,'(a)') head_sol
           WRITE(12,*) TRIM(head_sol)
           CLOSE(101)
          
           DO i = 1, SIZE(solution%ww,2)
             WRITE(12,*) i-1, grid % rr(:,i),  solution % ww(1:,i) 
           ENDDO  

     ! CASE (SU2_EU)

           
     !       !CALL file_header(Nb, solution%k_d, sol_fmt)
           
     !       OPEN (UNIT=101, FILE='solution.'//TRIM(problem_name),                  STATUS='old')
     !       READ(101,'(a)') head_sol
     !       WRITE(12,*) TRIM(head_sol)
     !       CLOSE(101)

     !       !WRITE(12,*) 'Dummy'
     !       DO i = 1, SIZE(solution%ww,2)
     !         WRITE(12,1088) i-1, grid % rr(:,i),    &      ! ID and x,y coordinates 
     !                    solution % ww(1:p,i), &   ! Conservative variables
     !                    solution % ww(p+3:10,i)     ! P, T, CP, Mach
     !       ENDDO                                                       

     !  CASE (SU2_NSL)
           
     !       !CALL file_header(Nb, solution%k_d, sol_fmt)
 
     !       OPEN (UNIT=101, FILE='solution.'//TRIM(problem_name),                  STATUS='old')
     !       READ(101,'(a)') head_sol
     !       WRITE(12,*) TRIM(head_sol)
     !       CLOSE(101)
           
     !       !WRITE(12,*) 'Dummy'
     !       DO i = 1, SIZE(solution%ww,2)
     !         WRITE(12,1089) i-1, grid % rr(:,i),    &      ! ID and x,y coordinates 
     !                    solution % ww(1:p,i), &   ! Conservative variables
     !                    solution % ww(p+3:14,i)     ! P, T, CP, Mach, 
     !                                                 ! mu, Cf, qw, yp
     !       ENDDO     

     !  CASE (SU2_NSSA)
           
     !       !CALL file_header(Nb, solution%k_d, sol_fmt)

     !       OPEN (UNIT=101, FILE='solution.'//TRIM(problem_name),                  STATUS='old')
     !       READ(101,'(a)') head_sol
     !       WRITE(12,*) TRIM(head_sol)
     !       CLOSE(101)

     !       !WRITE(12,*) 'Dummy'
     !       DO i = 1, SIZE(solution%ww,2)
     !         WRITE(12,1090) i-1, grid % rr(:,i),    &   ! ID and x,y coordinates 
     !                    solution % ww(1:p,i), &   ! Conservative variables
     !                    solution % ww(p+1,i),   &   ! SA eddy viscosity 
     !                    solution % ww(p+3:15,i)     ! P, T, CP, Mach, 
     !                                                 ! mu, Cf, qw, yp, mu_turb
     !       ENDDO           
      
    
      CASE DEFAULT 

           WRITE(*,*) ''
           WRITE(*,*) 'ERROR. WRITE_INTERPOLATED_SOLUTION'
           WRITE(*,*) 'Invalid solution format.'
           WRITE(*,*) ''
           STOP

    END SELECT


  !CLOSE (11)
  CLOSE (12)

  
  1085  FORMAT(i7,1e24.16)
  1086  FORMAT(7x,3e24.16)
  1087  FORMAT(7x,1e24.16)
  1088  FORMAT(i7,2x,10(1e24.16,1x))
  1089  FORMAT(i7,2x,14(1e24.16,1x))
  1090  FORMAT(i7,2x,16(1e24.16,1x)) 

  
  CONTAINS
  
  SUBROUTINE  file_header(Nb, sD, sol_fmt)
  !--------------------------------------------------------------------------------
  IMPLICIT NONE
    
  INTEGER, INTENT(IN) :: Nb, sD, sol_fmt
  
  CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: string
  
  REAL(KIND=8), DIMENSION(sD) :: mom
  REAL(KIND=8), DIMENSION(2) :: time_resid
  REAL(KIND=8) :: R1, R2, R3
  LOGICAL :: L1, L2
  
  INTEGER :: i, it, bi, bt
  !--------------------------------------------------------------------------------

  IF (sol_fmt == PIG_EU)  ALLOCATE (string(65+Nb))
  IF (sol_fmt == PIG_NSL) ALLOCATE (string(68+Nb))
  
  DO i = 1, 10
    READ(11,*)  string(i)
    WRITE(12,*) string(i)
  ENDDO
  
  READ(11,*)  it, time_resid
  WRITE(12,*) it, time_resid
  
  READ(11,*)  string(1)
  READ(11,*)  string(2)
  WRITE(12,*) string(1)
  WRITE(12,*) string(2)

  READ(11,*)  it, L1, L2, R1
  WRITE(12,*) it, L1, L2, R1
  
  IF (sol_fmt == PIG_EU) THEN
    DO i = 1, 36
      READ(11,*)  string(i)
      WRITE(12,*) string(i)
    ENDDO
  ELSEIF (SOL_FMT == PIG_NSL) THEN
    DO i = 1, 23
      READ(11,*)  string(i)
      WRITE(12,*) string(i)
    ENDDO
    READ(11,*) it
    WRITE(12,*) it
    IF (it == 1) THEN
      READ(11,*) R1
      WRITE(12,*) R1
    ELSE 
      READ(11,*) R1, R2, R3
      WRITE(12,*) R1, R2, R3  
    ENDIF
    DO i = 1, 14
      READ(11,*)  string(i)
      WRITE(12,*) string(i)
    ENDDO    
  ENDIF
  

  DO i = 1, Nb
    READ(11,*)  bi, bt
    WRITE(12,*) bi, bt
  ENDDO
  
  DO i = 1, 7
    READ(11,*)  string(i)
    WRITE(12,*) string(i)
  ENDDO
  
  READ(11,*)  mom
  WRITE(12,*) mom
  
  DO i = 1, 7
    READ(11,*)  string(i)
    WRITE(12,*) string(i)
  ENDDO 
  
  DEALLOCATE (string)
  
  END SUBROUTINE file_header
  
  
 END SUBROUTINE Write_Interpolated_Solution
 
 
 
 
 
  SUBROUTINE  Post_Interpolated_Solution(grid, sol_fmt, solution, problem_name)
  !----------------------------------------------------------------------------------------------------------------
  IMPLICIT NONE
 
  TYPE(grid_type),     INTENT(IN) :: grid
  INTEGER,             INTENT(IN) :: sol_fmt
  TYPE(solution_type), INTENT(IN) :: solution 
  CHARACTER(LEN=30),   INTENT(IN) :: problem_name
  
  INTEGER :: i, p
  !----------------------------------------------------------------------------------------------------------------
 
  p = SIZE(solution%ww,1)
 
  OPEN (UNIT = 11, FILE = 'intsol.'//TRIM(problem_name)//'.plt', STATUS = 'unknown')
  
    WRITE(11,*) 'TITLE = "Interpolated Solution"'

    IF (sol_fmt == 2  .OR.  sol_fmt == 5) THEN
      WRITE(11,*) 'VARIABLES = "X", "Y", "rho", "Momentum X", "Momentum Y", "Total Energy", "mu_tilde"'
    ELSE
      WRITE(11,*) 'VARIABLES = "X", "Y", "rho", "Momentum X", "Momentum Y", "Total Energy"'
    ENDIF

    WRITE(11,*) 'ZONE T = "Domain", N=', grid % Nj_d, 'E=', grid % Nm_d, 'DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
    
    DO i = 1, grid % Nj_d
      WRITE(11,*) grid % rr(:,i), solution % ww(:, i) 
    ENDDO
  
    DO i = 1, grid % Nm_d
      IF (grid % ele_type_d(i) == 2)  WRITE(11,*) grid % j_m_d(i) % vec, grid % j_m_d(i) % vec(3) 
      IF (grid % ele_type_d(i) == 3)  WRITE(11,*) grid % j_m_d(i) % vec
    ENDDO   
         
  CLOSE(11)
  
  WRITE(*,*) ' - intsol.'//TRIM(problem_name)//'.plt'
     
 END SUBROUTINE  Post_Interpolated_Solution





 SUBROUTINE Post_Error(grid, estimator, name)
 !------------------------------------------------------------------------------------
  IMPLICIT NONE
 
  CHARACTER(LEN=30),    INTENT(IN) :: name
  TYPE(grid_type),      INTENT(IN) :: grid
  TYPE(estimator_type), DIMENSION(:), INTENT(IN) :: estimator
 

  INTEGER :: i, j, idf
  INTEGER :: n1, n2

  CHARACTER(LEN=12), DIMENSION(:), ALLOCATABLE :: variable_name
 
  LOGICAL,      DIMENSION(:),   ALLOCATABLE :: mask
  INTEGER,      DIMENSION(:),   ALLOCATABLE :: Nc_onB
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: middle_point
 !------------------------------------------------------------------------------------
  
  idf = 11 
    
  ALLOCATE( middle_point(grid%k_d, grid%Nc_d) )

  DO i = 1, grid%Nc_d
 
   n1 = grid%j_c(1,i)
   n2 = grid%j_c(2,i)
 
   middle_point(1,i) = ( grid%rr(1,n1) + grid%rr(1,n2) ) / 2.d0
   middle_point(2,i) = ( grid%rr(2,n1) + grid%rr(2,n2) ) / 2.d0
 
  ENDDO
 
 
 
  OPEN( UNIT = idf, FILE = 'errors.'//TRIM(name)//'.plt', FORM = 'formatted', STATUS = 'unknown' )

 
     WRITE(idf,*) 'TITLE = "Edge Error"'
 
     ALLOCATE ( variable_name( SIZE(estimator(1)%errors,1) ) )
 
     DO i = 1, SIZE( variable_name )
      WRITE(variable_name(i),'(a9,i1.1,a2)') '"Estimate ',i,'" '
     ENDDO
 
     WRITE(idf,*)  'VARIABLES = "X", "Y",', variable_name
 

 
     !--Domain
     WRITE(idf,*) 'ZONE T = "Domain", I=', grid%Nc_d, ' DATAPACKING=POINT, ZONETYPE=ORDERED'
 
     DO i = 1, grid%Nc_d
      WRITE(idf,*) middle_point(:,i), estimator(1)%errors(:,i)
     ENDDO



     !--Boundaries
     ALLOCATE( mask(grid%Nc_b) )
     ALLOCATE( Nc_onB( MAXVAL(grid%bound_c) ) )

     DO i = 1, MAXVAL(grid%bound_c)
 
      mask = .FALSE.
      WHERE ( grid%bound_c == i )
       mask = .TRUE.
      END WHERE
      Nc_onB(i) = COUNT(mask)
 
      WRITE(idf,*) 'ZONE T = "Boundary ', i, '", I=', Nc_onB(i), ' DATAPACKING=POINT, ZONETYPE=ORDERED'
 
      DO j = 1, grid%Nc_b
       IF ( grid%bound_c(j) == i )  &
         WRITE(idf,*) middle_point( :, grid%cd_cb(j) ), estimator(1)%errors( :, grid%cd_cb(j) )
      ENDDO
 
     ENDDO
     
     WRITE(idf,*) 'TEXT X=15, Y=95, H=2, T="Refinement Thresholds"'
     DO i = 1, estimator(1)%passages
      WRITE(idf,*) 'TEXT X=15, Y=', 95-3*i, ' H=2 T="p= ', i, ' -> ', estimator(1)%R_Thres(i,:), '"'
     ENDDO
     
     !IF ( ASSOCIATED(estimator(1)%C_Thres) ) THEN
     ! WRITE(idf,*) 'TEXT X=15, Y=75, H=2, T="Derefinement Thresholds"'
     ! DO i = 1, estimator(1)%passages
     !  WRITE(idf,*) 'TEXT X=15, Y=', 75-3*i, ' H=2 T="p= ', i, ' -> ', estimator(1)%C_Thres(i,:), '"'
     ! ENDDO
     !ENDIF 


  CLOSE(idf)

  WRITE(*,*) '   - errors.'//TRIM(name)//'.plt'

 END SUBROUTINE Post_Error





 SUBROUTINE Post_Grid (grid, grid_name, ref_type, boxes)
 !----------------------------------------------------------------------------------------
  IMPLICIT NONE
  
  TYPE(grid_type),               INTENT(IN) :: grid
  CHARACTER(LEN=30),             INTENT(IN) :: grid_name
  INTEGER,                       INTENT(IN) :: ref_type  
  TYPE(D_I_V_R), DIMENSION(:),   INTENT(IN) :: boxes
  
  INTEGER           :: i, j, idf, counter    
  
  LOGICAL, DIMENSION(:), ALLOCATABLE :: mask_j, mask_m
  INTEGER, DIMENSION(:), ALLOCATABLE :: Nj_onB, Nm_onB
  INTEGER, DIMENSION(:), ALLOCATABLE :: reorder_Np
 !----------------------------------------------------------------------------------------
  
   idf = 11 
     
   OPEN( UNIT = idf, FILE = 'admesh.'//TRIM(grid_name)//'.plt' )     
 
     
     WRITE(idf,*) 'TITLE = "Adapted Mesh"'
     WRITE(idf,*) 'VARIABLES = "X", "Y"'
     
     
     !--Domain     
     WRITE(idf,*) 'ZONE T = "Domain", N=', grid%Nj_d, 'E=', grid%Nm_d, 'DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
     
     DO i = 1, grid%Nj_d
      WRITE(idf,*) grid%rr(:,i)
     ENDDO
        
     DO i = 1, grid%Nm_d
      IF ( grid%ele_type_d(i) == 2 ) WRITE(idf,*) grid%j_m_d(i)%vec, grid%j_m_d(i)%vec(3) 
      IF ( grid%ele_type_d(i) == 3 ) WRITE(idf,*) grid%j_m_d(i)%vec
     ENDDO          

     WRITE(idf,*)

     !--Boundaries
     ALLOCATE( mask_j(grid%Nj_b) )
     ALLOCATE( Nj_onB( MAXVAL(grid%bound_p) ) )
     ALLOCATE( mask_m(grid%Nm_b) )
     ALLOCATE( Nm_onB( MAXVAL(grid%bound_m) ) )
     
     ALLOCATE(reorder_Np(grid%Nj_b))
     
     DO i = 1, MAXVAL(grid%bound_p)  ! Loop over boundaries
 
      mask_j = .FALSE.
      WHERE ( grid%bound_p == i )
       mask_j = .TRUE.
      END WHERE
      Nj_onB(i) = COUNT(mask_j)
      
      mask_m = .FALSE.
      WHERE ( grid%bound_m == i )
       mask_m = .TRUE.
      END WHERE
      Nm_onB(i) = COUNT(mask_m)
      
      counter = 0      
 
      WRITE(idf,*) 'ZONE T = "Boundary ', i, '", N=', Nj_onB(i), 'E=', Nm_onB(i), 'DATAPACKING=POINT, ZONETYPE=FELINESEG' 

      DO j = 1, grid%Nj_b
      
       IF ( grid%bound_p(j) == i ) THEN
        WRITE(idf,*) grid%rr( :, grid%jd_jb(j) )
        counter = counter + 1
        reorder_Np(j) = counter
       ENDIF       
              
      ENDDO
      
         
      DO j = 1, grid%Nm_b
       IF ( grid%bound_m(j) == i ) WRITE(idf,*) reorder_Np( grid%j_m_b(j)%vec(1) ), reorder_Np( grid%j_m_b(j)%vec(2) )
      ENDDO      
          
     ENDDO

     IF (ref_type == 1) THEN

       DO i = 1, SIZE(boxes)

         WRITE(idf,*) 'ZONE T = "Box ', i, '", N=', SIZE(boxes(i) % vec)/2, &
                      'E=', SIZE(boxes(i) % vec)/2, 'DATAPACKING=POINT, ZONETYPE=FELINESEG'

         DO j = 1, SIZE(boxes(i) % vec), 2
           WRITE(idf,*) boxes(i) % vec(j), boxes(i) % vec(j+1) 
         ENDDO

         DO j = 1, SIZE(boxes(i) % vec)/2 - 1
           WRITE(idf,*) j, j+1
         ENDDO
         WRITE(idf,*) j, 1  

         WRITE(idf,*) ''

       ENDDO

     ENDIF

   CLOSE(idf)
 
   WRITE(*,*) '   - admesh.'//TRIM(grid_name)//'.plt'

   OPEN( UNIT = idf, FILE = 'admesh.'//TRIM(grid_name)//'.su2' )    !(S)

     WRITE(idf,'(a1)') '%'
     WRITE(idf,'(a21)') '% Name = Adapted Mesh'
     WRITE(idf,'(a1)') '%'
     WRITE(idf,'(a19)') '% Problem dimension'
     WRITE(idf,'(a1)') '%'
     WRITE(idf,'(a7,i7)') 'NDIME= ', grid%k_d
     WRITE(idf,'(a1)') '%'
     WRITE(idf,'(a28)') '% Inner element connectivity'
     WRITE(idf,'(a1)') '%'
     WRITE(idf,'(a7,i7)') 'NELEM= ', grid%Nm_d

     DO i=1, grid%Nm_d
       IF(grid%ele_type_d(i) .eq. 2) THEN
         WRITE(idf,'(5(i7))') grid%ele_type_d(i)+3, grid%j_m_d(i)%vec-1, i-1
       ELSE IF(grid%ele_type_d(i) .eq. 3) THEN
         WRITE(idf,'(6(i7))') grid%ele_type_d(i)+6, grid%j_m_d(i)%vec-1, i-1
       ENDIF
     ENDDO

     WRITE(idf,'(a1)') '%'
     WRITE(idf,'(a18)') '% Node coordinates'
     WRITE(idf,'(a1)') '%'
     WRITE(idf,'(a7,i7)') 'NPOIN= ', grid%Nj_d

     DO i=1, grid%Nj_d
       WRITE (idf,'(2(1x,e22.16),i7)') grid%rr(:,i), i-1
     ENDDO

     WRITE(idf,*) '%'
     WRITE(idf,*) '% Boundary elements'
     WRITE(idf,*) '%'
     WRITE(idf,'(a7,i7)') 'NMARK= ', grid%Nb

     counter=0
     DO i = 1, grid%Nm_b
       IF ( grid%bound_m(i) .ne. counter ) THEN
         counter=grid%bound_m(i)
         WRITE(idf,'(a12,i7)') 'MARKER_TAG= ', counter
         WRITE(idf,'(a14,i7)') 'MARKER_ELEMS= ', COUNT(grid%bound_m==counter)
       END IF
       WRITE(idf,'(3(i7))') grid%ele_type_b(i)+2, grid%jd_jb(grid%j_m_b(i)%vec)-1
     ENDDO

   WRITE(*,*) '   - admesh.'//TRIM(grid_name)//'.su2'

   CLOSE(idf)                                                       !(S)

 END SUBROUTINE Post_Grid


END MODULE io_proc
