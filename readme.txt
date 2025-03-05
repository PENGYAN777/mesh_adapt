Modification:

1. in code/grid_util.f90, inside EDGE_STRUCTURE_GEN function, for ASSOCAITED, we should ues NULLIY not DEALLOCATED, since the variable is a pointer to static ,not dynamics memory.
2. in code/error_estimaor.f90, insisde Anisotropic, copy several line of codes from Second_Derivetive to include nodes on the boundary.
