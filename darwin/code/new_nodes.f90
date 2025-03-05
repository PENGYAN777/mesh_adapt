MODULE new_nodes

 USE structures
 USE grid_utils

 CONTAINS
 
 SUBROUTINE Add_Nodes(grid, old_grid, d_np_to_rnd, b_np_to_rnd, nds, &
                      s_new, jd_jsb, updt_sb, rj_b, ref_history)

!==================================================================================!
 IMPLICIT NONE
 
 TYPE(grid_type),               INTENT(IN)    :: old_grid
 TYPE(E_R_P),                   INTENT(IN)    :: ref_history 
 
 TYPE(grid_type),               INTENT(INOUT) :: grid
 INTEGER,                       INTENT(INOUT) :: nds
 INTEGER,       DIMENSION(:),   INTENT(INOUT) :: rj_b
 TYPE(crv_crd), DIMENSION(:),   INTENT(INOUT) :: s_new 
 TYPE(cnv_idx), DIMENSION(:),   INTENT(INOUT) :: jd_jsb
 INTEGER,       DIMENSION(:,:), INTENT(INOUT) :: d_np_to_rnd
 INTEGER,       DIMENSION(:,:), INTENT(INOUT) :: b_np_to_rnd 
 REAL(KIND=8),  DIMENSION(:,:), INTENT(INOUT) :: updt_sb  

 INTEGER                                      :: edge_index, new_node                                            
 INTEGER                                      :: nd1, nd2
 INTEGER                                      :: nd1a, nd1b, nd2a, nd2b
 INTEGER                                      :: j, d_idx, b_idx
 INTEGER                                      :: b_c, indice
 INTEGER,       DIMENSION(2)                  :: b_n1, b_n2  
 INTEGER,       DIMENSION(:),   ALLOCATABLE   :: jsb            
 REAL(KIND=8),  DIMENSION(:),   ALLOCATABLE   :: sb
 
 REAL(KIND=8)                                 :: s_sx, s_dx 
!==================================================================================!

 ALLOCATE (sb(grid%Nj_b))
 sb = grid % curve_coordinates

 ALLOCATE (jsb(old_grid%Nb))
 jsb = 0
 indice = 1

 d_idx = old_grid % Nj_d + 1  
 b_idx = old_grid % Nj_b + 1
 
 
 DO j = 1, ref_history % N_refined_edges
 
   new_node = ref_history % inserted_nodes(j)
    
   nd1 = ref_history % edge_nodes(1,j)
   nd2 = ref_history % edge_nodes(2,j)

   edge_index = Extract_np_index(old_grid % c_j, nd1, nd2)
  
   
   IF (old_grid % cb_cd(edge_index) .EQ. 0) THEN 
        
    !--DOMAIN: 
      grid % rr(:,new_node) = (old_grid%rr(:,nd1) + old_grid%rr(:,nd2)) / 2.d0

    !--For each refined edge assing new nd idx    
      d_np_to_rnd(edge_index,1) = new_node           
         
   ELSEIF (old_grid%cb_cd(edge_index) .NE. 0) THEN 




    !_______________________________________________________________________________


     grid % rr(:,new_node) = ( old_grid%rr(:,nd1) + old_grid%rr(:,nd2) ) / 2.d0

    !_______________________________________________________________________________





    !--BOUNDARY:
    !--From domain notation to boundary notation
      nd1a = old_grid%jb_jd( 1, nd1 )  
      nd1b = old_grid%jb_jd( 2, nd1 )  

      nd2a = old_grid%jb_jd( 1, nd2 )  
      nd2b = old_grid%jb_jd( 2, nd2 ) 
         
    !--From nodes boundary index to boundary line index
      b_n1(1) = old_grid%bound_p( nd1a )
      b_n1(2) = old_grid%bound_p( nd1b )

      b_n2(1) = old_grid%bound_p( nd2a )
      b_n2(2) = old_grid%bound_p( nd2b )        
    
    !--Extract boundary line index to which edge considered belongs
      b_c = old_grid%bound_c( old_grid%cb_cd(edge_index) )
      
    !--For each boundary line counts the number of nodes to add to
    !--that line. It is an index for the new nodes added over a
    !--specified boundary line
      jsb(b_c) = jsb(b_c) + 1 


    !--Selecting left and right curve coordinate in order to compute 
    !--new node curve coordinate
       IF ( ( b_n1(1) .EQ. b_n1(2) ) .AND. sb(nd1a) .EQ. 0.d0  ) THEN

         IF ( ABS(1.d0-sb(nd2a)) .GT. ABS(0.d0-sb(nd2a)) ) THEN      
           s_sx = 0.d0
           s_dx = sb(nd2a)           
         ELSE      
           s_sx = 1.d0
           s_dx = sb(nd2a)         
         ENDIF
           
       ELSEIF ( ( b_n2(1) .EQ. b_n2(2) ) .AND. sb(nd2a) .EQ. 0.d0  ) THEN
        
         IF ( ABS(1.d0-sb(nd1a)) .GT. ABS(0.d0-sb(nd1a)) ) THEN      
           s_sx = 0.d0
           s_dx = sb(nd1a)           
         ELSE      
           s_sx = 1.d0
           s_dx = sb(nd1a)         
         ENDIF  
        
       ELSE       

         IF ( b_n1(1) .EQ. b_c ) THEN
           s_sx = sb(nd1a)
         ELSE
           s_sx = sb(nd1b)  
         ENDIF

         IF ( b_n2(1) .EQ. b_c ) THEN
           s_dx = sb(nd2a)
         ELSE
           s_dx = sb(nd2b)  
         ENDIF
                          
       ENDIF     
 
      !--Informations necessary for spline interpolation:
       s_new(b_c)%points(jsb(b_c)) = ( s_sx + s_dx ) / 2.d0  !-curve coordinate                          
       jd_jsb(b_c)%index(jsb(b_c)) = new_node                !-domain index
       rj_b(b_c) = jsb(b_c)                                  !-number of nodes added

    
      !PRINT*,  'node couple', nd1, nd2
      !PRINT*,  'left right and middle curve coord', s_sx, s_dx, ( s_sx + s_dx ) / 2.d0      
      !PRINT*,  'new node index', new_node
    

      !--For each refined np assing new nd idx
       d_np_to_rnd(edge_index,1) = new_node

       grid%jd_jb( b_idx ) = new_node
       grid%bound_p( b_idx ) = old_grid%bound_c(old_grid%cb_cd(edge_index))     
       b_np_to_rnd( grid%cb_cd(edge_index) ,1 ) = b_idx

       !-Copying informations for updating file curve_coords:
       ! curve coordinate
       updt_sb(b_idx-old_grid%Nj_b,1) = s_new(b_c)%points(jsb(b_c))     
       ! boundary line index
       updt_sb(b_idx-old_grid%Nj_b,2) = b_c*1.0

       b_idx = b_idx + 1



       grid%jd_jb( b_idx ) = new_node
       grid%bound_p( b_idx ) = old_grid%bound_c(old_grid%cb_cd(edge_index))
       b_np_to_rnd( grid%cb_cd(edge_index), 2 ) = b_idx

       !-Copying informations for updating file curve_coords    
       ! curve coordinate
       updt_sb(b_idx-old_grid%Nj_b,1) = s_new(b_c)%points(jsb(b_c))
       ! boundary line index    
       updt_sb(b_idx-old_grid%Nj_b,2) = b_c*1.0

       b_idx = b_idx + 1


    ENDIF      
                        
    d_idx = d_idx + 1    
 
 ENDDO   
   
 nds = d_idx


 END SUBROUTINE Add_Nodes
        


END MODULE new_nodes
