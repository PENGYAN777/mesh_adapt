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

MODULE  nodes

!============================================================
!   This module provide the data strucuture in which the 
!   geometry of the mesh is stored, and the procedures to 
!   read it and save it. 
!   Geometry is defined through the coordinates of the nodes
!   and the connectivity between boundary nodes and domain
!   nodes.
!   Both the data structures and subroutines defined
!   in this modules can deal with one, two or three 
!   dimensional meshes
!============================================================

!============================================================
!
!     k_d     physical dimension of the space
!
!     Np_d    Total number of nodes of the domain
!
!     Np_b    Total number of nodes of the boundary
!    
!  rr(k,i)    k-th coordinate of the i-th node
!
!     jd_jb   Connectivity between the node jb in the 
!             boundary numeration and the corresponding node 
!             jd in the numeration of the domain 
!             [ jd = jd_jb(jb)]
!
!   bound_p   Connectivity between the node p and the 
!             boundary portion
!
!============================================================


   !============================================================
   USE dynamic_vector
   !============================================================

   !============================================================
   INTEGER  ::  k_d, Np_d, Np_b
   REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  ::  rr
   INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  jd_jb 
   INTEGER,      DIMENSION(:),   ALLOCATABLE  ::  bound_p
   !============================================================


!============================================================
CONTAINS
!============================================================



   !============================================================
   SUBROUTINE set_bound_p(j_m_b, bound_m)
   !============================================================


      IMPLICIT NONE


      TYPE(D_I_V), DIMENSION(:), INTENT(IN)  ::  j_m_b
      INTEGER,     DIMENSION(:), INTENT(IN)  ::  bound_m

      INTEGER   ::  m


      IF(ALLOCATED(bound_P)) DEALLOCATE(bound_p)
      ALLOCATE( bound_p(Np_b))

      DO m = 1, SIZE(j_m_b)

         bound_p(j_m_b(m)%vec) = bound_m(m)

      ENDDO

      do m = 1, size(bound_p) 
       print*, m, bound_p(m)
      enddo

   END SUBROUTINE set_bound_p
   !============================================================



   !============================================================
   SUBROUTINE read_nodes (idf)
   !============================================================


      IMPLICIT NONE

      INTEGER, INTENT(IN)  ::  idf

      INTEGER  ::  j   ! Generic node index   
      INTEGER  ::  idx

      CHARACTER (LEN=16)  ::  name


      READ (idf,*)  
      READ (idf,'(16x,a15)') name 

      !WRITE(*,*) '   Reading nodes:  ',name 

      READ (idf,*); READ (idf,*)

      READ (idf,*) k_d, Np_d, Np_b 

      !WRITE(*,*) '  Spatial dimension(s): ',k_d
      !WRITE(*,*) '  Domain nodes:   ',Np_d
      !WRITE(*,*) '  Boundary nodes: ',Np_b
          
      IF ( ALLOCATED(rr) )       DEALLOCATE (rr)       
      IF ( ALLOCATED(jd_jb) )    DEALLOCATE (jd_jb)
      IF ( ALLOCATED(bound_p) )  DEALLOCATE (bound_p)

      ALLOCATE ( rr(k_d,Np_d) )
      ALLOCATE ( jd_jb(Np_b), bound_p(Np_b) )


      ! Volume nodes' coordinates 
      ! -------------------------
      READ (idf,*); READ (idf,*); READ (idf,*)
      READ (idf,*); READ (idf,*); READ (idf,*)

      DO j = 1, Np_d

          READ(idf,*) idx

          IF ( j .NE. idx ) THEN
              WRITE(*,*) ' Warning, nodes numeration not ordered '
              WRITE(*,*) ' Using internal (sequential) numeration '
          ENDIF

          READ(idf,*) rr(:,j)

      ENDDO


      ! Boundary nodes - domain nodes connectivity
      ! ------------------------------------------
      READ (idf,*); READ (idf,*); READ (idf,*)
      READ (idf,*); READ (idf,*)

      DO j = 1, Np_b

          READ(idf,*) idx, jd_jb(j), bound_p(j)

          IF ( j .NE. idx ) THEN
              WRITE(*,*) ' Warning, nodes numeration not ordered '
              WRITE(*,*) ' Using internal (sequential) numeration '
          ENDIF
          
      ENDDO


   END SUBROUTINE read_nodes
   !============================================================



   !============================================================
   SUBROUTINE save_nodes (idf, name, name_length)
   !============================================================


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idf
      CHARACTER (LEN=64), INTENT(IN) :: name
      INTEGER ::  name_length 
 
      INTEGER :: i
    
      WRITE (idf,1000)
      WRITE (idf,1010) name(1:name_length)

      WRITE (idf,1000)
      WRITE (idf,1020)
      WRITE (idf,1021) k_d, Np_d, Np_b 

      WRITE (idf,1000)
      WRITE (idf,1025)

      WRITE (idf,1000)
      WRITE (idf,1048)
      WRITE (idf,1050)
      WRITE (idf,1051)

      ! Volume nodes coordinates
      ! ------------------------
      DO i=1,Np_d
         WRITE (idf,1052) i
         WRITE (idf,1053) rr(:,i)
      ENDDO

      WRITE (idf,1000)
      WRITE (idf,1055)

      WRITE (idf,1000)
      WRITE (idf,1048)     
      WRITE (idf,1075)

      ! Boundary nodes - domain nodes connectivity
      ! ------------------------------------------
      DO i=1,Np_b
         WRITE (idf,1076) i, jd_jb(i), bound_p(i)
      ENDDO
      

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


   END SUBROUTINE save_nodes
   !============================================================



END MODULE  nodes
