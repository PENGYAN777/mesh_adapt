!============================================================ 
!
!      Module: messages
!
! Description: Procedures for printing messages on screen 
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

MODULE messages

   INTEGER, PARAMETER  ::  MSG_FILE = 6, &
                           ERR_FILE = 6

   CONTAINS

   SUBROUTINE  write_prepro_logo(idf)

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: idf

   WRITE(idf,*) ''
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   NODE-PAIR CODE   -   PREPROCESSING PHASE'
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   Department of Aerospace Engineering - Politecnico di Milano'   
   WRITE(idf,*) '   authors: Alberto Guardone, Marco Fossati'
   WRITE(idf,*) ''
   WRITE(idf,*) '   Warning: Check the status of file CORNERS.grid_name'
   WRITE(idf,*) ''
   
   END SUBROUTINE  write_prepro_logo





   SUBROUTINE  write_prescale_logo(idf)

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: idf

   WRITE(idf,*) ''
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   NODE-PAIR CODE   -   PRE-SCALE UTILITY'
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   Department of Aerospace Engineering - Politecnico di Milano'   
   WRITE(idf,*) '   author: Marco Fossati'
   WRITE(idf,*) ''
   WRITE(idf,*) ''

   END SUBROUTINE  write_prescale_logo





   SUBROUTINE  write_npcode_logo(idf)

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: idf

   WRITE(idf,*) ''
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   NODE-PAIR CODE'
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   Department of Aerospace Engineering - Politecnico di Milano'   
   WRITE(idf,*) '   authors: Alberto Guardone, Marco Fossati'
   WRITE(idf,*) ''

   END SUBROUTINE  write_npcode_logo





   SUBROUTINE  write_pospro_logo(idf)

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: idf

   WRITE(idf,*) ''
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   NODE-PAIR CODE   -   POSTPROCESSING PHASE'
   WRITE(idf,*) ' -------------------------------------------------------------------'
   WRITE(idf,*) '   Department of Aerospace Engineering - Politecnico di Milano'   
   WRITE(idf,*) '   authors: Alberto Guardone, Marco Fossati'
   WRITE(idf,*) ''
   WRITE(idf,*) ''

   END SUBROUTINE  write_pospro_logo





   SUBROUTINE  terminate(idf, name, msg)

   IMPLICIT NONE
   CHARACTER(*), INTENT(IN) :: name, msg
   INTEGER,      INTENT(IN) :: idf

   WRITE(idf,*) ''
   WRITE(idf,*) 'ERROR. ', TRIM(name)
   WRITE(idf,*) msg
   WRITE(idf,*) ''

   STOP

   END SUBROUTINE  terminate





   SUBROUTINE  write_line(ph)

   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: ph

   WRITE(MSG_FILE, * ) 
   WRITE(MSG_FILE, * )  ph
   WRITE(MSG_FILE, * ) 

100 FORMAT('------------------------------------------------------------')

   END SUBROUTINE write_line





   SUBROUTINE  write_message(msg)

   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: msg

   WRITE(MSG_FILE, * ) 
   WRITE(MSG_FILE,100)
   WRITE(MSG_FILE, * )  msg
   WRITE(MSG_FILE,100)
   WRITE(MSG_FILE, * ) 

100 FORMAT('------------------------------------------------------------')

   END SUBROUTINE write_message





   SUBROUTINE write_title(tle)

   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: tle

   WRITE(MSG_FILE, * ) 
   WRITE(MSG_FILE,100)
   WRITE(MSG_FILE,100)
   WRITE(MSG_FILE, * ) 
   WRITE(MSG_FILE, * )  tle
   WRITE(MSG_FILE, * ) 
   WRITE(MSG_FILE,100)
   WRITE(MSG_FILE,100)
   WRITE(MSG_FILE, * ) 

100 FORMAT('============================================================')

   END SUBROUTINE write_title





   SUBROUTINE write_error(err)

   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: err

   WRITE(ERR_FILE, * ) 
   WRITE(ERR_FILE,100)
   WRITE(ERR_FILE,110)
   WRITE(ERR_FILE, * ) 
   WRITE(ERR_FILE, * )  err
   WRITE(ERR_FILE, * ) 
   WRITE(ERR_FILE,110)
   WRITE(ERR_FILE,100)
   WRITE(ERR_FILE, * ) 

100 FORMAT('============================================================')
110 FORMAT('--ERROR--X--ERROR--X--ERROR--X--ERROR--X--ERROR--X--ERROR---')

   END SUBROUTINE write_error





   SUBROUTINE write_warning(war)

   IMPLICIT NONE

   CHARACTER(*), INTENT(IN) :: war

   WRITE(ERR_FILE, * ) 
   WRITE(ERR_FILE,100)
   WRITE(ERR_FILE,110)
   WRITE(ERR_FILE, * ) 
   WRITE(ERR_FILE, * )  war
   WRITE(ERR_FILE, * ) 
   WRITE(ERR_FILE,110)
   WRITE(ERR_FILE,100)
   WRITE(ERR_FILE, * ) 

100 FORMAT('------------------------------------------------------------')
110 FORMAT('--WARNING---WARNING---WARNING---WARNING---WARNING---WARNING--')

   END SUBROUTINE write_warning


END MODULE messages
