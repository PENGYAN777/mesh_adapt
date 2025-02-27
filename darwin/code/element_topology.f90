!============================================================ 
!
!    Author:  Alberto Guardone
!             Dipartimento di Ingegneria Aerospaziale
!             Politecnico di Milano
!             Via La Msa 34, 20156 Milano, ITALY
!             e-mail: guardone@aero.polimi.it
!
! Copyright: 1998-2003 Alberto Guardone
!            See COPYING file for copyright notice
!
!============================================================ 

MODULE element_topology


!============================================================ 
!
!        TYPE   CORRESPONDS TO   
!
!           1   segment 
!
!           2   triangle
!           3   quadrilater 
!
!           4   tetrahedron
!           5   pyramid 
!           6   tri-prism 
!           7   quad-prism 
!
!============================================================ 

   !============================================================ 
   USE dynamic_vector
   !============================================================ 

   !============================================================ 
   INTEGER, DIMENSION(1:7), PARAMETER  :: N_nodes_ele &
                                          = (/ 2,3,4,4,5,6,8 /)
   INTEGER, DIMENSION(1:7), PARAMETER  :: N_faces_ele &
                                          = (/ 2,3,4,4,5,5,6 /)
   INTEGER, DIMENSION(1:7), PARAMETER  :: N_edges_ele &
                                          = (/ 1,3,4,6,8,9,12 /)
                                          
   INTEGER, PARAMETER  ::  ELE_TYPE_SEGMENT     = 1, &
                           ELE_TYPE_TRIANGLE    = 2, &
                           ELE_TYPE_QUADRILATER = 3, &
                           ELE_TYPE_TETRAHEDRON = 4, &
                           ELE_TYPE_PYRAMID     = 5, &
                           ELE_TYPE_TRIPRISM    = 6, &
                           ELE_TYPE_QUADPRISM   = 7
   !============================================================ 



!============================================================ 
CONTAINS
!============================================================ 


   !============================================================ 
   FUNCTION ele_faces(ele_type) RESULT (faces)
   !============================================================ 

   
      !------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)  ::  ele_type
      
      TYPE(D_I_V), DIMENSION(:), POINTER  :: faces 
      !------------------------------------------------------------
      
      

      SELECT CASE(ele_type)
      
      
         !---------------------
         CASE(ELE_TYPE_SEGMENT)
         !---------------------
         
         ALLOCATE( faces(2) )   

         ALLOCATE( faces(1)%vec(1) )
         ALLOCATE( faces(2)%vec(1) )

         faces(1)%vec(1) = 2
         faces(2)%vec(1) = 1
            
         !----------------------
         CASE(ELE_TYPE_TRIANGLE)
         !----------------------

            ALLOCATE( faces(3) )

            ALLOCATE( faces(1)%vec(2) )
            ALLOCATE( faces(2)%vec(2) )
            ALLOCATE( faces(3)%vec(2) )
            
            faces(1)%vec  =  (/ 2, 3 /)
            faces(2)%vec  =  (/ 3, 1 /)
            faces(3)%vec  =  (/ 1, 2 /)
      
         !-------------------------
         CASE(ELE_TYPE_QUADRILATER)
         !-------------------------
          
            ALLOCATE( faces(4) )

            ALLOCATE( faces(1)%vec(2) )
            ALLOCATE( faces(2)%vec(2) )
            ALLOCATE( faces(3)%vec(2) )
            ALLOCATE( faces(4)%vec(2) )
            
            faces(1)%vec  =  (/ 1, 2 /)
            faces(2)%vec  =  (/ 2, 3 /)
            faces(3)%vec  =  (/ 3, 4 /)
            faces(4)%vec  =  (/ 4, 1 /)
            
         !-------------------------
         CASE(ELE_TYPE_TETRAHEDRON) 
         !-------------------------
          
            ALLOCATE( faces(4) )

            ALLOCATE( faces(1)%vec(3) )
            ALLOCATE( faces(2)%vec(3) )
            ALLOCATE( faces(3)%vec(3) )
            ALLOCATE( faces(4)%vec(3) )
            
            faces(1)%vec  =  (/ 1, 2, 3/)
            faces(2)%vec  =  (/ 1, 4, 2/)
            faces(3)%vec  =  (/ 2, 4, 3/)
            faces(4)%vec  =  (/ 1, 3, 4/)
            
         !---------------------
         CASE(ELE_TYPE_PYRAMID) 
         !---------------------
          
            ALLOCATE( faces(5) )

            ALLOCATE( faces(1)%vec(4) )
            ALLOCATE( faces(2)%vec(3) )
            ALLOCATE( faces(3)%vec(3) )
            ALLOCATE( faces(4)%vec(3) )
            ALLOCATE( faces(5)%vec(3) )
            
            faces(1)%vec  =  (/ 1, 2, 3, 4 /)
            faces(2)%vec  =  (/ 1, 5, 2 /)
            faces(3)%vec  =  (/ 2, 5, 3 /)
            faces(4)%vec  =  (/ 3, 5, 4 /)
            faces(5)%vec  =  (/ 4, 5, 1 /)
            
             
         !----------------------
         CASE(ELE_TYPE_TRIPRISM)         
         !----------------------
          
            ALLOCATE( faces(5) )

            ALLOCATE( faces(1)%vec(3) )
            ALLOCATE( faces(2)%vec(3) )
            ALLOCATE( faces(3)%vec(4) )
            ALLOCATE( faces(4)%vec(4) )
            ALLOCATE( faces(5)%vec(4) )
            
            faces(1)%vec  =  (/ 1, 2, 3 /)
            faces(2)%vec  =  (/ 4, 6, 5 /)
            faces(3)%vec  =  (/ 1, 4, 5, 2 /)
            faces(4)%vec  =  (/ 2, 5, 6, 3 /)
            faces(5)%vec  =  (/ 3, 6, 4, 1 /)
            
             
         !-----------------------
         CASE(ELE_TYPE_QUADPRISM)  
         !-----------------------
        
            ALLOCATE( faces(6) )

            ALLOCATE( faces(1)%vec(4) )
            ALLOCATE( faces(2)%vec(4) )
            ALLOCATE( faces(3)%vec(4) )
            ALLOCATE( faces(4)%vec(4) )
            ALLOCATE( faces(5)%vec(4) )
            ALLOCATE( faces(6)%vec(4) )
            
            faces(1)%vec  =  (/ 1, 2, 3, 4 /)
            faces(2)%vec  =  (/ 5, 8, 7, 6 /)
            faces(3)%vec  =  (/ 1, 5, 6, 2 /)
            faces(4)%vec  =  (/ 2, 6, 7, 3 /)
            faces(5)%vec  =  (/ 3, 7, 8, 4 /)
            faces(6)%vec  =  (/ 1, 4, 8, 5 /)

      END SELECT
   
      
   END FUNCTION ele_faces
   !============================================================ 



   !============================================================ 
   FUNCTION ele_edges(ele_type) RESULT (edges)
   !============================================================ 

   
      !------------------------------------------------------------
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)  ::  ele_type
      
      INTEGER, DIMENSION(:,:), POINTER  :: edges 
      !------------------------------------------------------------
      
      

      SELECT CASE(ele_type)
      
      
         !---------------------
         CASE(ELE_TYPE_SEGMENT)
         !---------------------
         
            ALLOCATE( edges(1,2) )
            
            edges(1,:)  =  (/ 1, 2 /)
            
            
         !----------------------
         CASE(ELE_TYPE_TRIANGLE)
         !----------------------
            ALLOCATE( edges(3,2) )
            
            edges(1,:)  =  (/ 1, 2 /)
            edges(2,:)  =  (/ 2, 3 /)
            edges(3,:)  =  (/ 3, 1 /)
      
      
         !-------------------------
         CASE(ELE_TYPE_QUADRILATER)
         !-------------------------
          
            ALLOCATE( edges(4,2) )
            
            edges(1,:)  =  (/ 1, 2 /)
            edges(2,:)  =  (/ 2, 3 /)
            edges(3,:)  =  (/ 3, 4 /)
            edges(4,:)  =  (/ 4, 1 /)
            
            
         !-------------------------
         CASE(ELE_TYPE_TETRAHEDRON) 
         !-------------------------
          
            ALLOCATE( edges(6,2) )
            
            edges(1,:)  =  (/ 1, 2 /)
            edges(2,:)  =  (/ 2, 3 /)
            edges(3,:)  =  (/ 3, 1 /)
            edges(4,:)  =  (/ 1, 4 /)
            edges(5,:)  =  (/ 2, 4 /)
            edges(6,:)  =  (/ 3, 4 /)
            
            
         !---------------------
         CASE(ELE_TYPE_PYRAMID) 
         !---------------------
          
            ALLOCATE( edges(8,2) )
            
            edges(1,:)  =  (/ 1, 2 /)
            edges(2,:)  =  (/ 2, 3 /)
            edges(3,:)  =  (/ 3, 1 /)
            edges(4,:)  =  (/ 4, 1 /)
            edges(5,:)  =  (/ 1, 5 /)
            edges(6,:)  =  (/ 2, 5 /)
            edges(7,:)  =  (/ 3, 5 /)
            edges(8,:)  =  (/ 4, 5 /)
            
             
         !----------------------
         CASE(ELE_TYPE_TRIPRISM)         
         !----------------------
         
            ALLOCATE( edges(9,2) )
            
            edges(1,:)  =  (/ 1, 2 /)
            edges(2,:)  =  (/ 2, 3 /)
            edges(3,:)  =  (/ 3, 1 /)
            edges(4,:)  =  (/ 1, 4 /)
            edges(5,:)  =  (/ 2, 5 /)
            edges(6,:)  =  (/ 3, 6 /)
            edges(7,:)  =  (/ 4, 5 /)
            edges(8,:)  =  (/ 5, 6 /)
            edges(9,:)  =  (/ 6, 4 /)
            
             
         !-----------------------
         CASE(ELE_TYPE_QUADPRISM)  
         !-----------------------
        
            ALLOCATE( edges(12,2) )
            
            edges(1,:)  = (/ 1, 2 /)
            edges(2,:)  = (/ 2, 3 /)
            edges(3,:)  = (/ 3, 4 /)
            edges(4,:)  = (/ 4, 1 /)
            edges(5,:)  = (/ 1, 5 /)
            edges(6,:)  = (/ 2, 6 /)
            edges(7,:)  = (/ 3, 7 /)
            edges(8,:)  = (/ 4, 8 /)
            edges(9,:)  = (/ 5, 6 /)
            edges(10,:) = (/ 6, 7 /)
            edges(11,:) = (/ 7, 8 /)
            edges(12,:) = (/ 8, 5 /)
            

      END SELECT
   
   
   END FUNCTION ele_edges
   !============================================================ 





END MODULE element_topology
