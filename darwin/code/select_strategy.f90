MODULE select_strategy

 CONTAINS
 
  SUBROUTINE select_refinement(k, flag, e_type, nmr_ref_np, e_to_add, &
                               nds_to_add, bnds_to_add, strategy_id, np_idxs)
  !-------------------------------------------------------------------------------
  IMPLICIT NONE
  
  LOGICAL, DIMENSION(:), INTENT(IN) :: flag
  
  INTEGER, INTENT(IN) :: k
  INTEGER, INTENT(INOUT) :: e_to_add, nds_to_add, bnds_to_add
  INTEGER, DIMENSION(:), INTENT(IN) :: e_type, nmr_ref_np
  INTEGER, DIMENSION(:), INTENT(IN) :: np_idxs 
  
  CHARACTER(LEN=3), DIMENSION(:), INTENT(INOUT) :: strategy_id  
  !-------------------------------------------------------------------------------
 
   IF (flag(k)) THEN
   
      IF (e_type(k) == 1) THEN 
      
        strategy_id(k) = 'S1 '
        e_to_add = e_to_add + 1         
   
      ELSEIF (e_type(k) == 2 .AND. nmr_ref_np(k) == 1) THEN 
      
        strategy_id(k) = 'T1 '
        e_to_add = e_to_add + 1
              
      ELSEIF (e_type(k) == 2 .AND. nmr_ref_np(k) == 2) THEN
      
        strategy_id(k) = 'T2 '
        e_to_add = e_to_add + 2
        
      ELSEIF (e_type(k) == 2 .AND. nmr_ref_np(k) == 3) THEN
      
        strategy_id(k) = 'T3 '
        e_to_add = e_to_add + 3

      ELSEIF (e_type(k) == 3 .AND. nmr_ref_np(k) == 1) THEN     

        strategy_id(k) = 'Q1 '  
        e_to_add = e_to_add + 2
        nds_to_add = nds_to_add + 1
        bnds_to_add = bnds_to_add + 1

      ELSEIF (e_type(k) == 3 .AND. nmr_ref_np(k) == 2) THEN     

        IF ((np_idxs(2) /= 0) .AND. (np_idxs(3) /= 0) .OR.  &
            (np_idxs(1) /= 0) .AND. (np_idxs(4) /= 0)) THEN          

          strategy_id(k) = 'Q2p'
          e_to_add = e_to_add + 1       

        ELSE    

          strategy_id(k) = 'Q2 '
          e_to_add = e_to_add + 2
          nds_to_add = nds_to_add + 1
          bnds_to_add = bnds_to_add + 1

        ENDIF    
        
      ELSEIF (e_type(k) == 3 .AND. nmr_ref_np(k) == 3) THEN     

        strategy_id(k) = 'Q3 '
        e_to_add = e_to_add + 4
        nds_to_add = nds_to_add + 1
        bnds_to_add = bnds_to_add + 1
        
      ELSEIF (e_type(k) == 3 .AND. nmr_ref_np(k) == 4) THEN     

        strategy_id(k) = 'Q4 '
        e_to_add = e_to_add + 3
        nds_to_add = nds_to_add + 1
        bnds_to_add = bnds_to_add + 1
              
      ENDIF
                   
   ENDIF
   
  END SUBROUTINE select_refinement

END MODULE select_strategy
