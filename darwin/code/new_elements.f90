MODULE new_elements

 USE structures


 CONTAINS

 SUBROUTINE Build_Elements(     e_type,       elem_to_nds, elem_adj, patt,    &
                            svd_e_type,   svd_elem_to_nds,   svd_elem_adj,    &
                            np_to_rnd, np_flag, b_m, svd_b_m, nds, crop, rr,  &
                            c_m, zone, sol_ww, add_ww, Nnodes)

 ! IMPORTANT !!!
 ! Differently from triangles, edges of quadrilateral elements are not ordered 
 ! in clockwise order.  To recover clockwise order the array c_m_cw (c_m_[c]lock[w]ise)
 ! is used to swap edge 3 with edge 4.
 !------------------------------------------------------------------------------------------- 
 IMPLICIT NONE

 INTEGER,      DIMENSION(:),        INTENT(INOUT) :: e_type
 TYPE(D_I_V),  DIMENSION(:),        INTENT(INOUT) :: elem_to_nds, &
                                                     elem_adj                                                  
 TYPE(R_P),    DIMENSION(:),        INTENT(IN)    :: patt
 INTEGER,      DIMENSION(:),        INTENT(INOUT) :: svd_e_type 
 TYPE(D_I_V),  DIMENSION(:),        INTENT(INOUT) :: svd_elem_to_nds, &
                                                     svd_elem_adj
 INTEGER,      DIMENSION(:,:),      INTENT(IN)    :: np_to_rnd
 LOGICAL,      DIMENSION(:),        INTENT(IN)    :: np_flag
 INTEGER,      DIMENSION(:),        INTENT(INOUT) :: b_m
 INTEGER,      DIMENSION(:),        INTENT(IN)    :: svd_b_m
 INTEGER                                          :: nds
 INTEGER,                           INTENT(INOUT) :: crop
 REAL(KIND=8), DIMENSION(:,:),      INTENT(INOUT) :: rr
 TYPE(D_I_V),  DIMENSION(:),        INTENT(IN)    :: c_m
 CHARACTER(LEN=1),                  INTENT(IN)    :: zone
 REAL(KIND=8),      DIMENSION(:,:), INTENT(IN)    :: sol_ww
 REAL(KIND=8),      DIMENSION(:,:), INTENT(INOUT) :: add_ww
 INTEGER,                           INTENT(IN)    :: Nnodes

 TYPE(D_I_V), DIMENSION(:), ALLOCATABLE :: c_m_cw

 REAL(KIND=8) :: Xnd, Ynd

 INTEGER, DIMENSION(:), ALLOCATABLE :: node_pairs

 INTEGER :: j, l, elem_idx,          &
            rnd1, rnd2, rnd3, rnd4,  &
            np_idx, aux, nd1, nd2,   &
            nd3, nd4, nd5, nd6,      &
            nd7, nd8, nd9

 CHARACTER(LEN=1) :: config
 !-------------------------------------------------------------------------------------------

 elem_idx = 1

 ALLOCATE (c_m_cw(SIZE(svd_e_type)))
 
 DO j = 1, SIZE(svd_e_type)
  
   IF (patt(j) % ref_flag) THEN

      SELECT CASE (patt(j) % strategy)

        CASE ('S1 ') !-----------------------------------------------------------

          ! Writing nodes idx
          nd1 = svd_elem_to_nds(j) % vec(1)
          nd2 = np_to_rnd(j,1)
          nd3 = np_to_rnd(j,2)
          nd4 = svd_elem_to_nds(j) % vec(2)


          ! Writing new element  (1) 
          e_type(elem_idx) = 1
          b_m(elem_idx) = svd_b_m(j)

          ALLOCATE (elem_to_nds(elem_idx) % vec(2))
          ALLOCATE (elem_adj(elem_idx) % vec(2))
          elem_to_nds(elem_idx) % vec(1) = nd1
          elem_to_nds(elem_idx) % vec(2) = nd2
          
          elem_idx = elem_idx + 1         

          ! Writing new element  (2) 
          e_type(elem_idx) = 1
          b_m(elem_idx) = svd_b_m(j)

          ALLOCATE (elem_to_nds(elem_idx) % vec(2))
          ALLOCATE (elem_adj(elem_idx) % vec(2))
          elem_to_nds(elem_idx) % vec(1) = nd3
          elem_to_nds(elem_idx) % vec(2) = nd4

          elem_idx = elem_idx + 1                 

          
          
        CASE ('T1 ') !-----------------------------------------------------------

          ! Rearranging node pairs and nodes in order 
          ! to build new element. Edges are rearranged
          ! so that edge 1 is flagged for refinement
          l = 1

          ALLOCATE (node_pairs(SIZE(c_m(j)%vec)))
          node_pairs = c_m(j) % vec
1         np_idx = node_pairs(l)

          IF (.NOT. np_flag(np_idx)) THEN       
          
            ! Node pairs 
            aux             = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = aux
            
            ! Nodes 
            aux                         = svd_elem_to_nds(j)%vec(l+2)
            svd_elem_to_nds(j)%vec(l+2) = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = aux

            GOTO 1
          ENDIF
          
          ! New node's idx associated to 'old' j element
          rnd1 = np_to_rnd(node_pairs(1),1)

          ! Writing nodes idx 
          nd1 = svd_elem_to_nds(j)%vec(1)
          nd2 = rnd1
          nd3 = svd_elem_to_nds(j)%vec(2)
          nd4 = svd_elem_to_nds(j)%vec(3)


          ! Interpolating solution
          add_ww(:,nd2 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd3)) / 2


          ! Writing new element  (1) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )
          elem_to_nds(elem_idx)%vec(1) = nd1
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd4
          
          elem_idx = elem_idx + 1

          ! Writing new element  (2)
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )
          elem_to_nds(elem_idx)%vec(1) = nd2
          elem_to_nds(elem_idx)%vec(2) = nd3
          elem_to_nds(elem_idx)%vec(3) = nd4
          
          elem_idx = elem_idx + 1
          
          DEALLOCATE (node_pairs)                



        CASE ('T2 ') !-----------------------------------------------------------
        
          ! Rearranging node pairs  and nodes in order 
          ! to build new element. Edges are rearranged so that
          ! edge 1 and 3 are flagged for refinement
          l = 1  
          ALLOCATE (node_pairs(SIZE(c_m(j) % vec)))
          node_pairs = c_m(j) % vec

2         np_idx = node_pairs(l)
          IF (.NOT. np_flag(np_idx)) THEN       

            ! Node pairs 
            aux             = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = aux
            
            ! Nodes             
            aux                         = svd_elem_to_nds(j)%vec(l+2)
            svd_elem_to_nds(j)%vec(l+2) = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = aux

            GOTO 2

          ENDIF
        
          l = l + 1
          np_idx = node_pairs(l)
          IF (np_flag(np_idx)) THEN     
                    
            ! Node pairs 
            aux             = node_pairs(l-1)
            node_pairs(l-1) = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = aux
            
            ! Nodes             
            aux                         = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = svd_elem_to_nds(j)%vec(l-1)
            svd_elem_to_nds(j)%vec(l-1) = aux
            
            l = 1
            GOTO 2      
           
           ENDIF            

          ! New node's idx associated to 'old' j element 
          rnd1 = np_to_rnd( node_pairs(1),1 )
          rnd2 = np_to_rnd( node_pairs(3),1 )
                                                                          
          ! Writing nodes idx 
          nd1 = svd_elem_to_nds(j)%vec(1)
          nd2 = rnd1
          nd3 = svd_elem_to_nds(j)%vec(2)
          nd4 = rnd2
          nd5 = svd_elem_to_nds(j)%vec(3)       


          ! Interpolating solution
          add_ww(:,nd2 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd3)) / 2
          add_ww(:,nd4 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd5)) / 2


          ! Writing new element  (1) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )
          elem_to_nds(elem_idx)%vec(1) = nd1
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd5
          
          elem_idx = elem_idx + 1       

          ! Writing new element  (2) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )
          elem_to_nds(elem_idx)%vec(1) = nd2
          elem_to_nds(elem_idx)%vec(2) = nd3
          elem_to_nds(elem_idx)%vec(3) = nd4
          
          elem_idx = elem_idx + 1               
        
          ! Writing new element  (3) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )
          elem_to_nds(elem_idx)%vec(1) = nd2
          elem_to_nds(elem_idx)%vec(2) = nd4
          elem_to_nds(elem_idx)%vec(3) = nd5
          
          elem_idx = elem_idx + 1
                                
          DEALLOCATE (node_pairs)       
        
        
        
        CASE ('T3 ')  !-----------------------------------------------------------

          ! New nodes' idx associated to 'old' j element 
          rnd1 = np_to_rnd( c_m(j)%vec(1),1 )
          rnd2 = np_to_rnd( c_m(j)%vec(3),1 )
          rnd3 = np_to_rnd( c_m(j)%vec(2),1 )     

          ! Writing nodes idx 
          nd1 = svd_elem_to_nds(j)%vec(1)
          nd2 = rnd1
          nd3 = svd_elem_to_nds(j)%vec(2)
          nd4 = rnd2
          nd5 = svd_elem_to_nds(j)%vec(3)
          nd6 = rnd3


          ! Interpolating solution
          add_ww(:,nd2 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd3)) / 2
          add_ww(:,nd4 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd5)) / 2
          add_ww(:,nd6 - Nnodes) = (sol_ww(:,nd5) + sol_ww(:,nd1)) / 2


          ! Writing new element  (1) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )          
          elem_to_nds(elem_idx)%vec(1) = nd1
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd6        
          
          elem_idx = elem_idx + 1
          
          ! Writing new element  (2) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )          
          elem_to_nds(elem_idx)%vec(1) = nd2
          elem_to_nds(elem_idx)%vec(2) = nd3
          elem_to_nds(elem_idx)%vec(3) = nd4
          
          elem_idx = elem_idx + 1

          ! Writing new element  (3) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )          
          elem_to_nds(elem_idx)%vec(1) = nd4
          elem_to_nds(elem_idx)%vec(2) = nd5
          elem_to_nds(elem_idx)%vec(3) = nd6
          
          elem_idx = elem_idx + 1
          
          ! Writing new element  (4) 
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )          
          elem_to_nds(elem_idx)%vec(1) = nd6
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd4
          
          elem_idx = elem_idx + 1         



        CASE('Q1 ')  !-----------------------------------------------------------

          l = 1  
          ALLOCATE (c_m_cw(j) % vec(SIZE(c_m(j) % vec)))
          ALLOCATE (     node_pairs(SIZE(c_m(j) % vec)))
          
          c_m_cw(j) % vec(1) = c_m(j) % vec(1)
          c_m_cw(j) % vec(2) = c_m(j) % vec(2)
          c_m_cw(j) % vec(3) = c_m(j) % vec(4)
          c_m_cw(j) % vec(4) = c_m(j) % vec(3)
          
          node_pairs = c_m_cw(j) % vec

4         np_idx = node_pairs(l)
          IF (.NOT. np_flag(np_idx)) THEN       

            ! Node pairs
            aux             = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = node_pairs(l+3)
            node_pairs(l+3) = aux
            
            ! Nodes
            aux                           = svd_elem_to_nds(j) % vec(l+3)
            svd_elem_to_nds(j) % vec(l+3) = svd_elem_to_nds(j) % vec(l+2)
            svd_elem_to_nds(j) % vec(l+2) = svd_elem_to_nds(j) % vec(l+1)
            svd_elem_to_nds(j) % vec(l+1) = svd_elem_to_nds(j) % vec(l)
            svd_elem_to_nds(j) % vec(l)   = aux

            GOTO 4

          ENDIF 

          ! New node's idx associated to 'old' j element 
          rnd1 = np_to_rnd(node_pairs(1),1)

          ! Writing nodes idx 
          nd1 = svd_elem_to_nds(j) % vec(1)
          nd2 = svd_elem_to_nds(j) % vec(2)
          nd3 = svd_elem_to_nds(j) % vec(3)
          nd4 = svd_elem_to_nds(j) % vec(4)
          nd5 = rnd1      
          nd6 = nds

          ! Computing element subdivision configuration (T or Q).
          ! T if the node along the median line follows outside
          !   the quadrilater
          ! Q in the other case          
          rr(:,nds) = node_coords(rr, nd1, nd2, nd3, nd4, nd5, crop, config)

          IF (config == 'T') THEN

            ! Interpolating solution
            add_ww(:,nd5 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2)) / 2


            ! Writing new element  (1) 
            e_type(elem_idx) = 2

            ALLOCATE (elem_to_nds(elem_idx) % vec(3))
            ALLOCATE (elem_adj(elem_idx) % vec(3))        
            elem_to_nds(elem_idx) % vec(1) = nd1
            elem_to_nds(elem_idx) % vec(2) = nd5
            elem_to_nds(elem_idx) % vec(3) = nd4
            
            elem_idx = elem_idx + 1

            ! Writing new element  (2) 
            e_type(elem_idx) = 2

            ALLOCATE (elem_to_nds(elem_idx) % vec(3))
            ALLOCATE (elem_adj(elem_idx) % vec(3))
            elem_to_nds(elem_idx) % vec(1) = nd5
            elem_to_nds(elem_idx) % vec(2) = nd2
            elem_to_nds(elem_idx) % vec(3) = nd3
            
            elem_idx = elem_idx + 1             

            ! Writing new element  (3) 
            e_type(elem_idx) = 2

            ALLOCATE (elem_to_nds(elem_idx) % vec(3))
            ALLOCATE (elem_adj(elem_idx) % vec(3))
            elem_to_nds(elem_idx) % vec(1) = nd5
            elem_to_nds(elem_idx) % vec(2) = nd3
            elem_to_nds(elem_idx) % vec(3) = nd4

            elem_idx = elem_idx + 1                               
          
          ELSEIF (config == 'Q') THEN               

            ! Interpolating solution
            add_ww(:,nd5 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2)) / 2
            add_ww(:,nd6 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2) +  &
                                              sol_ww(:,nd3) + sol_ww(:,nd4)) / 4
                 
            nds = nds + 1 
          
            ! Writing new element  (1) 
            e_type(elem_idx) = 3

            ALLOCATE (elem_to_nds(elem_idx) % vec(4))
            ALLOCATE (elem_adj(elem_idx) % vec(4))  
            elem_to_nds(elem_idx) % vec(1) = nd1
            elem_to_nds(elem_idx) % vec(2) = nd5
            elem_to_nds(elem_idx) % vec(3) = nd6
            elem_to_nds(elem_idx) % vec(4) = nd4
            
            elem_idx = elem_idx + 1

            ! Writing new element  (2) 
            e_type(elem_idx) = 3

            ALLOCATE (elem_to_nds(elem_idx) % vec(4))
            ALLOCATE (elem_adj(elem_idx) % vec(4))  
            elem_to_nds(elem_idx) % vec(1) = nd5
            elem_to_nds(elem_idx) % vec(2) = nd2
            elem_to_nds(elem_idx) % vec(3) = nd3
            elem_to_nds(elem_idx) % vec(4) = nd6
            
            elem_idx = elem_idx + 1             
        
            ! Writing new element  (3) 
            e_type(elem_idx) = 2

            ALLOCATE (elem_to_nds(elem_idx) % vec(3))
            ALLOCATE (elem_adj(elem_idx) % vec(3))
            elem_to_nds(elem_idx) % vec(1) = nd6
            elem_to_nds(elem_idx) % vec(2) = nd3
            elem_to_nds(elem_idx) % vec(3) = nd4

            elem_idx = elem_idx + 1                     

          ENDIF

          DEALLOCATE (node_pairs, c_m_cw(j)%vec)


        
        CASE ('Q2 ')   !-----------------------------------------------------------

          ! Rearranging node pairs  and nodes in order 
          ! to build new element 
          l = 1  
          ALLOCATE (c_m_cw(j)%vec(SIZE(c_m(j)%vec)))
          ALLOCATE (   node_pairs(SIZE(c_m(j)%vec)))
          
          c_m_cw(j) % vec(1) = c_m(j) % vec(1)
          c_m_cw(j) % vec(2) = c_m(j) % vec(2)
          c_m_cw(j) % vec(3) = c_m(j) % vec(4)
          c_m_cw(j) % vec(4) = c_m(j) % vec(3)

          node_pairs = c_m_cw(j)%vec

5         np_idx = node_pairs(l)
          IF ( .NOT. np_flag(np_idx) ) THEN     

            ! Node pairs 
            aux             = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = node_pairs(l+3)
            node_pairs(l+3) = aux
            
            ! Nodes             
            aux                         = svd_elem_to_nds(j)%vec(l+3)
            svd_elem_to_nds(j)%vec(l+3) = svd_elem_to_nds(j)%vec(l+2)
            svd_elem_to_nds(j)%vec(l+2) = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = aux

            GOTO 5

          ENDIF 
          l = l + 1
          np_idx = node_pairs(l)
          IF ( np_flag(np_idx) ) THEN   
                    
            !--Node pairs---------------------------------------------!
            aux             = node_pairs(l-1)
            node_pairs(l-1) = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = aux
            
            !--Nodes--------------------------------------------------!            
            aux                         = svd_elem_to_nds(j)%vec(l+2)
            svd_elem_to_nds(j)%vec(l+2) = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = svd_elem_to_nds(j)%vec(l-1)
            svd_elem_to_nds(j)%vec(l-1) = aux
            
            l = 1
            GOTO 5      
           
           ENDIF            

          !-------New node's idx associated to 'old' j element--------------------------!
          rnd1 = np_to_rnd( node_pairs(1),1 )
          rnd2 = np_to_rnd( node_pairs(4),1 )

          !------Writing nodes idx------------------------------------------------------!
          nd1 = svd_elem_to_nds(j)%vec(1)
          nd2 = svd_elem_to_nds(j)%vec(2)
          nd3 = svd_elem_to_nds(j)%vec(3)
          nd4 = svd_elem_to_nds(j)%vec(4)
          nd5 = rnd1      
          nd6 = rnd2
          nd7 = nds          

          !-------Adding new node-------------------------------------------------------!
          Xnd = ( rr(1,nd1) + rr(1,nd2) + rr(1,nd3) + rr(1,nd4) ) / 4
          Ynd = ( rr(2,nd1) + rr(2,nd2) + rr(2,nd3) + rr(2,nd4) ) / 4
           
          rr(1,nds) = Xnd          
          rr(2,nds) = Ynd

          nds = nds + 1  


          ! Interpolating solution
          add_ww(:,nd5 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2)) / 2
          add_ww(:,nd6 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd2)) / 2
          add_ww(:,nd7 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2) +  &
                                            sol_ww(:,nd3) + sol_ww(:,nd4)) / 4


          !-------Writing new element  (1)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd1
          elem_to_nds(elem_idx)%vec(2) = nd5
          elem_to_nds(elem_idx)%vec(3) = nd7
          elem_to_nds(elem_idx)%vec(4) = nd4
          
          elem_idx = elem_idx + 1

          !-------Writing new element  (2)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd5
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd6
          elem_to_nds(elem_idx)%vec(4) = nd7
          
          elem_idx = elem_idx + 1               
        
          !-------Writing new element  (3)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd6
          elem_to_nds(elem_idx)%vec(2) = nd3
          elem_to_nds(elem_idx)%vec(3) = nd4
          elem_to_nds(elem_idx)%vec(4) = nd7
          
          elem_idx = elem_idx + 1                                 
          
          DEALLOCATE ( node_pairs, c_m_cw(j)%vec )                       


        CASE('Q2p')

          !-------Rearranging node pairs  and nodes in order to build new element-------!
          l = 1  
          ALLOCATE ( c_m_cw(j)%vec( SIZE(c_m(j)%vec) ) )
          ALLOCATE (    node_pairs( SIZE(c_m(j)%vec) ) )
          
          c_m_cw(j)%vec(1) = c_m(j)%vec(1)
          c_m_cw(j)%vec(2) = c_m(j)%vec(2)
          c_m_cw(j)%vec(3) = c_m(j)%vec(4)
          c_m_cw(j)%vec(4) = c_m(j)%vec(3)
          
          node_pairs = c_m_cw(j)%vec
          
6         np_idx = node_pairs(l)
          IF ( .NOT. np_flag(np_idx) ) THEN     

          !--Node pairs---------------------------------------------!
            aux             = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = node_pairs(l+3)
            node_pairs(l+3) = aux
            
          !--Nodes--------------------------------------------------!            
            aux                         = svd_elem_to_nds(j)%vec(l+3)
            svd_elem_to_nds(j)%vec(l+3) = svd_elem_to_nds(j)%vec(l+2)
            svd_elem_to_nds(j)%vec(l+2) = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = aux

            GOTO 6

          ENDIF  

!-------New node's idx associated to 'old' j element--------------------------!
          rnd1 = np_to_rnd( node_pairs(1),1 )
          rnd2 = np_to_rnd( node_pairs(3),1 )

!------Writing nodes idx------------------------------------------------------!
          nd1 = svd_elem_to_nds(j)%vec(1)
          nd2 = svd_elem_to_nds(j)%vec(2)
          nd3 = svd_elem_to_nds(j)%vec(3)
          nd4 = svd_elem_to_nds(j)%vec(4)
          nd5 = rnd1      
          nd6 = rnd2


          ! Interpolating solution
          add_ww(:,nd5 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2)) / 2
          add_ww(:,nd6 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd4)) / 2


!-------Writing new element  (1)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd1
          elem_to_nds(elem_idx)%vec(2) = nd5
          elem_to_nds(elem_idx)%vec(3) = nd6
          elem_to_nds(elem_idx)%vec(4) = nd4
          
          elem_idx = elem_idx + 1

!-------Writing new element  (2)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd5
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd3
          elem_to_nds(elem_idx)%vec(4) = nd6
          
          elem_idx = elem_idx + 1               
        
          DEALLOCATE ( node_pairs, c_m_cw(j)%vec )                      

        
        CASE('Q3 ')
        
!-------Rearranging node pairs  and nodes in order to build new element-------!
          l = 1  
          ALLOCATE ( c_m_cw(j)%vec( SIZE(c_m(j)%vec) ) )
          ALLOCATE (    node_pairs( SIZE(c_m(j)%vec) ) )
          
          c_m_cw(j)%vec(1) = c_m(j)%vec(1)
          c_m_cw(j)%vec(2) = c_m(j)%vec(2)
          c_m_cw(j)%vec(3) = c_m(j)%vec(4)
          c_m_cw(j)%vec(4) = c_m(j)%vec(3)
          
          node_pairs = c_m_cw(j)%vec

7         np_idx = node_pairs(l)
          IF ( .NOT. np_flag(np_idx) ) THEN     

          !--Node pairs---------------------------------------------!
            aux             = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = node_pairs(l+3)
            node_pairs(l+3) = aux
            
          !--Nodes--------------------------------------------------!            
            aux                         = svd_elem_to_nds(j)%vec(l+3)
            svd_elem_to_nds(j)%vec(l+3) = svd_elem_to_nds(j)%vec(l+2)
            svd_elem_to_nds(j)%vec(l+2) = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = aux

            GOTO 7

          ENDIF 
          l = l + 1
          np_idx = node_pairs(l)
          IF ( np_flag(np_idx) ) THEN   
                    
          !--Node pairs---------------------------------------------!
            aux             = node_pairs(l-1)
            node_pairs(l-1) = node_pairs(l)
            node_pairs(l)   = node_pairs(l+1)
            node_pairs(l+1) = node_pairs(l+2)
            node_pairs(l+2) = aux
            
          !--Nodes--------------------------------------------------!            
            aux                         = svd_elem_to_nds(j)%vec(l+2)
            svd_elem_to_nds(j)%vec(l+2) = svd_elem_to_nds(j)%vec(l+1)
            svd_elem_to_nds(j)%vec(l+1) = svd_elem_to_nds(j)%vec(l)
            svd_elem_to_nds(j)%vec(l)   = svd_elem_to_nds(j)%vec(l-1)
            svd_elem_to_nds(j)%vec(l-1) = aux
            
            l = 1
            GOTO 7      
           
           ENDIF            

!-------New node's idx associated to 'old' j element--------------------------!
          rnd1 = np_to_rnd( node_pairs(1),1 )
          rnd2 = np_to_rnd( node_pairs(4),1 )
          rnd3 = np_to_rnd( node_pairs(3),1 )
        
!------Writing nodes idx------------------------------------------------------!
          nd1 = svd_elem_to_nds(j)%vec(1)
          nd2 = svd_elem_to_nds(j)%vec(2)
          nd3 = svd_elem_to_nds(j)%vec(3)
          nd4 = svd_elem_to_nds(j)%vec(4)
          nd5 = rnd1      
          nd6 = rnd2
          nd7 = rnd3
          nd8 = nds          
                  
!-------Adding new node-------------------------------------------------------!
          Xnd = ( rr(1,nd1) + rr(1,nd2) + rr(1,nd3) + rr(1,nd4) ) / 4
          Ynd = ( rr(2,nd1) + rr(2,nd2) + rr(2,nd3) + rr(2,nd4) ) / 4
           
          rr(1,nds) = Xnd          
          rr(2,nds) = Ynd

          nds = nds + 1  


          ! Interpolating solution
          add_ww(:,nd5 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2)) / 2
          add_ww(:,nd6 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd2)) / 2
          add_ww(:,nd7 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd4)) / 2
          add_ww(:,nd8 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2) +  &
                                            sol_ww(:,nd3) + sol_ww(:,nd4)) / 4


!-------Writing new element  (1)----------------------------------------------!

          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )          
          elem_to_nds(elem_idx)%vec(1) = nd1
          elem_to_nds(elem_idx)%vec(2) = nd5
          elem_to_nds(elem_idx)%vec(3) = nd8
          
          elem_idx = elem_idx + 1

!-------Writing new element  (2)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd5
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd6
          elem_to_nds(elem_idx)%vec(4) = nd8
          
          elem_idx = elem_idx + 1               
        
!-------Writing new element  (3)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd8
          elem_to_nds(elem_idx)%vec(2) = nd6
          elem_to_nds(elem_idx)%vec(3) = nd3
          elem_to_nds(elem_idx)%vec(4) = nd7
          
          elem_idx = elem_idx + 1                                 
                
!-------Writing new element  (4)----------------------------------------------!
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )          
          elem_to_nds(elem_idx)%vec(1) = nd8
          elem_to_nds(elem_idx)%vec(2) = nd7
          elem_to_nds(elem_idx)%vec(3) = nd4
          
          elem_idx = elem_idx + 1                                 
          
!-------Writing new element  (5)----------------------------------------------!
          e_type(elem_idx) = 2
          ALLOCATE ( elem_to_nds(elem_idx)%vec(3) )
          ALLOCATE ( elem_adj(elem_idx)%vec(3) )          
          elem_to_nds(elem_idx)%vec(1) = nd8
          elem_to_nds(elem_idx)%vec(2) = nd4
          elem_to_nds(elem_idx)%vec(3) = nd1
          
          elem_idx = elem_idx + 1                                               
        
          DEALLOCATE ( node_pairs, c_m_cw(j)%vec )                      
        
        
        CASE('Q4 ')
!-------New nodes' idx associated to 'old' j element--------------------------!
          rnd1 = np_to_rnd( c_m(j)%vec(1),1 )
          rnd2 = np_to_rnd( c_m(j)%vec(3),1 )
          rnd3 = np_to_rnd( c_m(j)%vec(4),1 )
          rnd4 = np_to_rnd( c_m(j)%vec(2),1 )
          
!-------Writing nodes idx-----------------------------------------------------!
          nd1 = svd_elem_to_nds(j)%vec(1)
          nd2 = svd_elem_to_nds(j)%vec(2)
          nd3 = svd_elem_to_nds(j)%vec(3)
          nd4 = svd_elem_to_nds(j)%vec(4)
          nd5 = rnd1
          nd6 = rnd2
          nd7 = rnd3
          nd8 = rnd4
          nd9 = nds
 
!-------Adding new node-------------------------------------------------------!
          Xnd = ( rr(1,nd1) + rr(1,nd2) + rr(1,nd3) + rr(1,nd4) ) / 4
          Ynd = ( rr(2,nd1) + rr(2,nd2) + rr(2,nd3) + rr(2,nd4) ) / 4
           
          rr(1,nds) = Xnd          
          rr(2,nds) = Ynd

          nds = nds + 1                   


          ! Interpolating solution
          add_ww(:,nd5 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2)) / 2
          add_ww(:,nd6 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd2)) / 2
          add_ww(:,nd7 - Nnodes) = (sol_ww(:,nd3) + sol_ww(:,nd4)) / 2
          add_ww(:,nd8 - Nnodes) = (sol_ww(:,nd4) + sol_ww(:,nd1)) / 2
          add_ww(:,nd9 - Nnodes) = (sol_ww(:,nd1) + sol_ww(:,nd2) +  &
                                            sol_ww(:,nd3) + sol_ww(:,nd4)) / 4



!-------Writing new element  (1)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd1
          elem_to_nds(elem_idx)%vec(2) = nd5
          elem_to_nds(elem_idx)%vec(3) = nd9
          elem_to_nds(elem_idx)%vec(4) = nd8
          
          elem_idx = elem_idx + 1
          
!-------Writing new element  (2)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd5
          elem_to_nds(elem_idx)%vec(2) = nd2
          elem_to_nds(elem_idx)%vec(3) = nd6
          elem_to_nds(elem_idx)%vec(4) = nd9
          
          elem_idx = elem_idx + 1

!-------Writing new element  (3)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd9
          elem_to_nds(elem_idx)%vec(2) = nd6
          elem_to_nds(elem_idx)%vec(3) = nd3
          elem_to_nds(elem_idx)%vec(4) = nd7
          
          elem_idx = elem_idx + 1
          
!-------Writing new element  (4)----------------------------------------------!
          e_type(elem_idx) = 3
          ALLOCATE ( elem_to_nds(elem_idx)%vec(4) )
          ALLOCATE ( elem_adj(elem_idx)%vec(4) )          
          elem_to_nds(elem_idx)%vec(1) = nd8
          elem_to_nds(elem_idx)%vec(2) = nd9
          elem_to_nds(elem_idx)%vec(3) = nd7
          elem_to_nds(elem_idx)%vec(4) = nd4
          
          elem_idx = elem_idx + 1         
             
      END SELECT
           
   ELSE      
          
      e_type(elem_idx) = svd_e_type(j)
      
      IF ( zone == 'B' ) THEN
        b_m(elem_idx) = svd_b_m(j) 
      ENDIF
      
      elem_to_nds(elem_idx) = svd_elem_to_nds(j)
      elem_adj(elem_idx)    = svd_elem_adj(j)
      elem_idx              = elem_idx + 1
      
   ENDIF

 ENDDO

 DEALLOCATE (c_m_cw)

 CONTAINS


 
 FUNCTION  node_coords(rr, nd1, nd2, nd3, nd4, nd5, crop, config)
 
 IMPLICIT NONE
 !==Variables Definition======================================================!
 INTEGER,                      INTENT(IN)    :: nd1, nd2, nd3, nd4, nd5
 INTEGER,                      INTENT(INOUT) :: crop
 
 REAL(KIND=8)                                :: alpha
 REAL(KIND=8)                                :: beta_3, beta_4, beta, beta_half
 REAL(KIND=8)                                :: gamma, delta, t
 REAL(KIND=8)                                :: num, den, pi = 3.141592654
 REAL(KIND=8)                                :: x_m43, y_m43, lenght
 REAL(KIND=8)                                :: xx1, xx2, xx3, xx4, xx5, xxm
 REAL(KIND=8)                                :: yy1, yy2, yy3, yy4, yy5, yym
 REAL(KIND=8)                                :: xx_1, xx_2, xx_3, &
                                                xx_4, xx_5, xx_m
 REAL(KIND=8)                                :: yy_1, yy_2, yy_3, &
                                                yy_4, yy_5, yy_m
 REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: rr
 REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: points
 REAL(KIND=8), DIMENSION(2)                  :: node_coords, coords
 
 CHARACTER(LEN=1),            INTENT(INOUT)  :: config
 !============================================================================!

  x_m43 = (rr(1,nd4) + rr(1,nd3)) / 2.d0
  y_m43 = (rr(2,nd4) + rr(2,nd3)) / 2.d0

  num = rr(2,nd5) - y_m43
  den = rr(1,nd5) - x_m43

  IF (ABS(den) < 1.E-8)  den = 1.E-8
    
  alpha = ATAN2(num,den) 

  lenght = SQRT((rr(2,nd5) - y_m43)**2 + (rr(1,nd5) - x_m43)**2)


  ! Transforming coordinates in local reference system
  ! Origin of local system is point '5'. 'X' axis passes for points '5' 
  ! and 'm43'
  xx1 = rr(1,nd1) - rr(1,nd5);  yy1 = rr(2,nd1) - rr(2,nd5)
  xx2 = rr(1,nd2) - rr(1,nd5);  yy2 = rr(2,nd2) - rr(2,nd5)
  xx3 = rr(1,nd3) - rr(1,nd5);  yy3 = rr(2,nd3) - rr(2,nd5)
  xx4 = rr(1,nd4) - rr(1,nd5);  yy4 = rr(2,nd4) - rr(2,nd5)
  xxm =     x_m43 - rr(1,nd5);  yym =     y_m43 - rr(2,nd5)
  xx5 = 0.0;                    yy5 = 0.0
          
  xx_1 =  xx1*COS( pi-alpha ) - yy1*SIN( pi-alpha )
  yy_1 =  xx1*SIN( pi-alpha ) + yy1*COS( pi-alpha )

  xx_2 =  xx1*COS( pi-alpha ) - yy2*SIN( pi-alpha )
  yy_2 =  xx2*SIN( pi-alpha ) + yy2*COS( pi-alpha )

  xx_3 =  xx3*COS( pi-alpha ) - yy3*SIN( pi-alpha )
  yy_3 =  xx3*SIN( pi-alpha ) + yy3*COS( pi-alpha )

  xx_4 =  xx4*COS( pi-alpha ) - yy4*SIN( pi-alpha )
  yy_4 =  xx4*SIN( pi-alpha ) + yy4*COS( pi-alpha )

  xx_m =  xxm*COS( pi-alpha ) - yym*SIN( pi-alpha )       
  yy_m =  xxm*SIN( pi-alpha ) + yym*COS( pi-alpha )

  xx_5 = xx5
  yy_5 = yy5


  ! Computing quadrilater's angles in order to fix configuration:
  ! Configuration is determined by the intersection  of  two  straight  lines. 
  ! Line one passes for points '5' and 'm', line two passes  for  the  vertex
  ! associated to the minumun angle between '4' and '3', with slope  equal to
  ! half that angle
  ALLOCATE (points(2,4))
  points(1,1) = xx_1;         points(2,1) = yy_1
  points(1,2) = xx_2;         points(2,2) = yy_2
  points(1,3) = xx_3;         points(2,3) = yy_3
  points(1,4) = xx_4;         points(2,4) = yy_4
                  
  beta_3 = angle(points, 3)
  beta_4 = angle(points, 4)

  beta      = MIN(beta_3, beta_4)
  beta_half = beta/2.d0 

  IF (beta == beta_3) THEN
          
    points(1,1) = xx_5;         points(2,1) = yy_5
    points(1,2) = xx_2;         points(2,2) = yy_2
    points(1,3) = xx_3;         points(2,3) = yy_3
    points(1,4) = xx_m;         points(2,4) = yy_m        
     
    gamma = angle(points, 4)
    delta = pi - (beta_half + gamma)
       
    t = ( TAN(-delta)*(xx_5 - xx_3) + yy_3 - yy_5 ) / &
        ( yy_m - yy_5 - TAN(-delta)*(xx_m - xx_5) )                                                       
  ELSE

    points(1,1) = xx_1;         points(2,1) = yy_1
    points(1,2) = xx_5;         points(2,2) = yy_5
    points(1,3) = xx_m;         points(2,3) = yy_m
    points(1,4) = xx_4;         points(2,4) = yy_4
    

    ! delta is the angle of the line that bisects minimum angle with respect
    ! the local reference frame
    gamma = angle( points, 3)
    delta = pi - ( beta_half + gamma )
    
    t = ( TAN(delta)*(xx_5 - xx_4) + yy_4 - yy_5 ) / &
        ( yy_m - yy_5 - TAN(delta)*(xx_m - xx_5) )          

  ENDIF


  IF (t > 1.d0) THEN
  
    node_coords(1) = xx_5  
    node_coords(2) = yy_5                
  
  ELSEIF (t < 0.d0) THEN

    node_coords(1) = xx_5
    node_coords(2) = yy_5

  ELSEIF (t > 0.d0  .AND.  t < 1.d0) THEN

    node_coords(1) = xx_5 + t*(xx_m - xx_5)
    node_coords(2) = yy_5 + t*(yy_m - yy_5)

  ENDIF                           

  IF (SQRT((node_coords(2) - yy_5)**2 + &
           (node_coords(1) - xx_5)**2 ) <  1.E-1*lenght) THEN
      
    config = 'T' 
    crop   = crop + 1
    
  ELSE

    config = 'Q'

  ENDIF  


  ! Transforming coordinates in global reference system
  coords(1) =  node_coords(1)*COS( pi-alpha ) + &
               node_coords(2)*SIN( pi-alpha )
           
  coords(2) = -node_coords(1)*SIN( pi-alpha ) + &
               node_coords(2)*COS( pi-alpha )

  node_coords(1) = coords(1) + rr(1,nd5)
  node_coords(2) = coords(2) + rr(2,nd5)

 END FUNCTION node_coords


 END SUBROUTINE Build_Elements





 FUNCTION angle(points, which)

 ! Computes angles of polygon specyfied by its vertices coordinates and
 ! returns the angle (in radiants) specified by 'which' index
 !-----------------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER                                   :: dim, l, mx, mn
 INTEGER                                   :: which
 
 REAL(KIND=8)                              :: angle
 REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: points 
 REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: vrs_x, vrs_y, alpha
 !-----------------------------------------------------------------------------

 !--WHICH identifies which angle the function should return:------------------!
 !  which = 0   return minimum angle
 !  which = 1   return alpha(1)
 !  which = 2   return alpha(2)
 !  which = 3   return alpha(3)
 !  which = 4   return alpha(4)
 
 !--POINTS is the vector containing the coordinates of the vertices ordered---!
 !  according to user needs---------------------------------------------------!

 !--vrs_(x,y) = vector containing edges' versors:-----------------------------!
 !  vrs_(x,y)(1) = edge 2-1
 !  vrs_(x,y)(2) = edge 3-2
 !
 !                 _ edge 3-1 (if triangle) 
 !  vrs_(x,y)(3) |_ 
 !                  edge 4-3 (if quadrilater)
 !
 !  vrs_(x,y)(4) = edge 4-1 (if quadrilater)----------------------------------!
 
 !--alpha = vector containing angles:-----------------------------------------!
 !  alpha(1) = angle between edges crossing in vertex 1
 !  alpha(2) = angle between edges crossing in vertex 2
 !  alpha(3) = angle between edges crossing in vertex 3
 !  alpha(4) = angle between edges crossing in vertex 4 (if quadrilater)------!
 
 
 dim = SIZE(points,2)
 ALLOCATE( vrs_x(dim), vrs_y(dim) )
 ALLOCATE( alpha(dim) )
 
 DO l = 1, dim
 
   IF ( l+1 .GT. dim )THEN
      mx = 1
      mn = dim
   ELSE  
      mx = l+1
      mn = l
   ENDIF  

   vrs_x(l) = ( points(1,mx) - points(1,mn) ) / SQRT( ( points(2,mx)-points(2,mn) )**2 &
                                                    + ( points(1,mx)-points(1,mn) )**2 )

   vrs_y(l) = ( points(2,mx) - points(2,mn) ) / SQRT( ( points(2,mx)-points(2,mn) )**2 &
                                                    + ( points(1,mx)-points(1,mn) )**2 )
 ENDDO
 
 DO l = 1, dim
 
   IF ( l+1 .GT. dim )THEN
      mx = 1
      mn = dim
   ELSE  
      mx = l+1
      mn = l
   ENDIF    
   
   alpha(mx) = -vrs_x(mx)*vrs_x(mn) - vrs_y(mx)*vrs_y(mn)
   alpha(mx) = ACOS(alpha(mx))
     
 ENDDO  
 
 IF (which == 0) THEN
   angle = MINVAL(alpha)                                       
 ELSE
   angle = alpha(which)
 ENDIF    
 
 END FUNCTION angle


END MODULE new_elements
