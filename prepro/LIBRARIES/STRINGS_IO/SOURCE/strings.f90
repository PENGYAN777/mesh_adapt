!============================================================ 
!
!      Module: strings
!
! Description: Procedures for manipulating strings 
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

MODULE strings


   !============================================================ 
   ! Procedures: private/public policy
   ! ---------------------------------
   PUBLIC   ::  last_c_leng
   !============================================================ 



!============================================================ 
 CONTAINS
!============================================================ 



   !============================================================ 
   FUNCTION last_c_leng (len_str, string) RESULT (leng)
   !============================================================ 

      ! Returns the length of the string

      IMPLICIT NONE                                                               

      INTEGER, INTENT(IN) :: len_str 
      CHARACTER (LEN=len_str), INTENT(IN) :: string 
      INTEGER :: leng     

      INTEGER :: i  

      leng = len_str
      DO i=1,len_str 

         IF ( string(i:i) .EQ. ' ' ) THEN 
            leng = i-1; EXIT
         ENDIF  

      ENDDO


   END FUNCTION last_c_leng 
   !============================================================ 



END MODULE strings
