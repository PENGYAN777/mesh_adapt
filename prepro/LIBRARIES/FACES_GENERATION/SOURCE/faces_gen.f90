!============================================================ 
!
!      Module: faces_gen
!
! Description: In this module the faces data structure is 
!              defined.  
!
!    Author:   Marco Fossati
!              Dipartimento di Ingegneria Aerospaziale
!              Politecnico di Milano
!              Via La Masa 34, 20156 Milano, ITALY
!              e-mail: fossati@aero.polimi.it
!
!============================================================ 

!============================================================ 
!
!   Quantities belonging to the domain
!   ----------------------------------
!   
!      Nf_d   Total number of faces in the domain
!
!     j_f_d   Connectivity between the face f and its 
!             nodes j 
!
!     f_c_d   Connectivity between the node-pair c and its 
!             faces f 
!
!============================================================ 

   MODULE  faces_gen

   USE dynamic_vector
   USE element_topology
   USE mesh_structure
   USE node_pair_structure
   USE nodes

   IMPLICIT NONE
   
   TYPE (D_I_V), DIMENSION(:), POINTER :: j_f_d
   TYPE (D_I_V), DIMENSION(:), POINTER :: f_c_d
   
   INTEGER :: Nf_d
   

   CONTAINS


   SUBROUTINE  faces_str_gen
   !---------------------------------------------------------
   IMPLICIT NONE
   
   TYPE(D_I_V), DIMENSION(:), POINTER :: m_faces

   INTEGER :: m, jC, f_, mf, f, jf
   
   LOGICAL, DIMENSION(:), ALLOCATABLE :: present_face_A
   LOGICAL :: write_face
   !---------------------------------------------------------
   
   ! Count the number of internal faces twice + the number of
   ! boundary faces that are counted once. Then divide the number
   ! faces counted twice by two and sum the number of boundary faces
   Nf_d = ((COUNT(ele_type_d == 4)*4 + COUNT(ele_type_d == 5)*5 + &
            COUNT(ele_type_d == 6)*5 + COUNT(ele_type_d == 7)*6) - Ne_b) / 2 &
	   + Ne_b
   
   ALLOCATE (j_f_d(Nf_d))
   
   f_ = 0
   DO m = 1, Ne_d

     SELECT CASE(ele_type_d(m))

       CASE(ELE_TYPE_TETRAHEDRON)	    
       ALLOCATE (m_faces(4))

       CASE(ELE_TYPE_PYRAMID) 
       ALLOCATE (m_faces(5))
          
       CASE(ELE_TYPE_TRIPRISM)         
       ALLOCATE (m_faces(5))

       CASE(ELE_TYPE_QUADPRISM)        
       ALLOCATE (m_faces(6))

     END SELECT
   
     ! Returns the nodes of each face in a local numbering     
     m_faces = ele_faces(ele_type_d(m))

     DO mf = 1, SIZE (m_faces)
     
       ! Convert the local numbering of the face's nodes
       ! into the numbering of the nodes of the grid
       m_faces(mf) % vec = j_m_d(m)% vec((m_faces(mf) % vec))

       IF (ALLOCATED(present_face_A))  DEALLOCATE (present_face_A)
       ALLOCATE (present_face_A(SIZE(m_faces(mf) % vec)))

       ! Loop over the already defined faces to check
       ! if face mf is already written in j_f_d
       write_face = .TRUE.
       current_face: DO f = 1, f_, 1

         DO jf = 1, SIZE(m_faces(mf) % vec)
	 
	   jC = m_faces(mf) % vec(jf)
	   
	   IF (ANY(jC == j_f_d(f)%vec)) THEN
	     present_face_A(jf) = .TRUE.
	   ELSE
	     present_face_A(jf) = .FALSE.
	     EXIT
	   ENDIF           

	 ENDDO
	 
	 IF (ALL(present_face_A)  .AND.  &
	     SIZE(m_faces(mf) % vec) ==  SIZE(j_f_d(f)%vec)) THEN
	   write_face = .FALSE.
	   EXIT
	 ENDIF  
     
       ENDDO  current_face


       ! Write face to j_f_d
       IF (write_face) THEN

         f_ = f_ + 1
         ALLOCATE (j_f_d(f_) % vec(SIZE(m_faces(mf) % vec)))
       
         j_f_d(f_) % vec = m_faces(mf) % vec	   

       ENDIF

     ENDDO

   ENDDO

   END SUBROUTINE  faces_str_gen





   SUBROUTINE  f_c_connectivity
   !------------------------------------------------------------------------------   
   IMPLICIT NONE
   
   INTEGER :: c, f, i, j, Nfc, f_
   LOGICAL :: connected = .FALSE.
   !------------------------------------------------------------------------------

   ALLOCATE (f_c_d(Nc_fv))

   DO c = 1, Nc_fv
   
     i = j_c_fv(1,c)
     j = j_c_fv(2,c)
 
     ! Count the number of faces associated
     ! to each single edge
     Nfc = 0
     DO f = 1, Nf_d
     
       IF (ANY(i == j_f_d(f)%vec)  .AND. &
           ANY(j == j_f_d(f)%vec))  connected = .TRUE.

       IF (connected) THEN
         Nfc = Nfc + 1
         connected = .FALSE.
       ENDIF	 

     ENDDO

     ALLOCATE (f_c_d(c)%vec(Nfc))

     ! Fill with the index of the face
     f_ = 1
     DO f = 1, Nf_d
     
       IF (ANY(i == j_f_d(f)%vec)  .AND. &
           ANY(j == j_f_d(f)%vec))  connected = .TRUE.

       IF (connected) THEN
         f_c_d(c) % vec(f_) = f
	 f_ = f_ + 1
         connected = .FALSE.	 
       ENDIF

     ENDDO

   ENDDO


! do c = 1, Nc_fv
! print*, c, f_c_d(c)%vec
! enddo



   END SUBROUTINE  f_c_connectivity





   SUBROUTINE  plot_faces(name)
   !------------------------------------------------------------------------------
   IMPLICIT NONE
   
   CHARACTER(LEN=64) :: name
   INTEGER :: i, j, Np_f, m
   !------------------------------------------------------------------------------
      
   OPEN (UNIT=100, FILE='faces.'//TRIM(name)//'.plt')                                           

   ! Writes header and data
   WRITE(100,*) 'VARIABLES = "x", "y", "z"'

   WRITE(100,*) 'ZONE T= "Outline", N=', Np_d, ', E=', Ne_d, ' F=FEPOINT, ET=BRICK'
   DO i = 1, Np_d									  
     WRITE(100,*) rr(:,i)								  
   ENDDO

   DO m = 1, Ne_d										    
   												    
     SELECT CASE (ele_type_d(m))								    

       CASE (4)  ! Tetrahedron									    
       WRITE(100,*) j_m_d(m) % vec(1:3), j_m_d(m) % vec(3), &					    
   		    j_m_d(m) % vec(4),   j_m_d(m) % vec(4), j_m_d(m) % vec(4), j_m_d(m) % vec(4)    
   												    
       CASE (5)  ! Pyramid									    
       WRITE(100,*) j_m_d(m) % vec(1:4), j_m_d(m) % vec(5), &					    
   		    j_m_d(m) % vec(5),   j_m_d(m) % vec(5), j_m_d(m) % vec(5)			    

       CASE (6)  ! Prism (triangular base)							    
       WRITE(100,*) j_m_d(m) % vec(1:3), j_m_d(m) % vec(3), j_m_d(m) % vec(4:6), j_m_d(m) % vec(6)  

       CASE (7)  ! Prism (quadrilateral base)							    
       WRITE(100,*) j_m_d(m) % vec								    

       CASE DEFAULT										    
       WRITE(*,*)										    
       WRITE(*,*) 'ERROR. Unknown element type' 						    
       WRITE(*,*) 'Element', m, 'unknown flag', ele_type_d(m)					    
       WRITE(*,*)										    
       STOP											    

     END SELECT 										    

   ENDDO											    

 
   DO i = 1, Nf_d									  

     Np_f = SIZE(j_f_d(i) % vec)

     WRITE(100,*) 'ZONE T= "Face', i, '", N=', Np_f, ', E=1 F=FEPOINT, ET=BRICK'

     DO j = 1, Np_f
       WRITE(100,*) rr(:,j_f_d(i)%vec(j))
     ENDDO
     
     IF (SIZE(j_f_d(i)%vec) == 3)  WRITE(100,*) 1, 2, 3, 3, 1, 2, 3, 3
     IF (SIZE(j_f_d(i)%vec) == 4)  WRITE(100,*) 1, 2, 3, 4, 1, 2, 3, 4
     
   ENDDO

   CLOSE (100)

   END SUBROUTINE  plot_faces


   END MODULE  faces_gen
