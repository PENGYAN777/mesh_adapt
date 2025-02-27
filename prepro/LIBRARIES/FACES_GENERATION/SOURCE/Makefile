################################################################################
#
#  NODE-PAIR CODE   -   Faces structures Library
#
#  Department of Aerospace Engineering
#  Politecnico di Milano. 1995-2008
#
################################################################################

SHELL	 = /bin/sh

.SUFFIXES:
.SUFFIXES: .f .f90 .o

#-------------------------------------------------------------------------------
# SETTING COMPILERS VARIABLES
#-------------------------------------------------------------------------------
VOID=
CMPL_CFG_DIR=../../../PREPRO/SOURCE

ifneq ($(CMPL),$(VOID))
    
  ifeq ($(DEBUG),yes)
    include $(CMPL_CFG_DIR)/Makefile.$(CMPL).DBG
  else
    include $(CMPL_CFG_DIR)/Makefile.$(CMPL)    
  endif
    
endif  

#-------------------------------------------------------------------------------
# FOLDERS
#-------------------------------------------------------------------------------
DV_DIR	   = ../../DYNAMIC_VECTORS/SOURCE
MESH_DIR   = ../../MESH_STRUCTURES/SOURCE
FEM_DIR    = ../../FINITE_ELEMENT_METHOD/SOURCE
NP_DIR	   = ../../NODEPAIR_STRUCTURES/SOURCE
ALG_DIR    = ../../ALGEBRA/SOURCE
MPROC_DIR  = ../../MULTI_PROCESSOR/SOURCE
COMM_DIR   = ../../COMMONS/SOURCE

DIRS = ./          \
       $(DV_DIR)   \
       $(MESH_DIR) \
       $(FEM_DIR)  \
       $(NP_DIR)   \
       $(ALG_DIR)  \
       $(PLOT_DIR) \
       $(COMM_DIR) \
       $(MPROC_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
DV_LIB     = $(DV_DIR)/libdyn_vec.a
MESH_LIB   = $(MESH_DIR)/libmesh.a
FEM_LIB    = $(FEM_DIR)/libfem.a
NP_LIB     = $(NP_DIR)/libnp.a
ALG_LIB    = $(ALG_DIR)/libalgebra.a
MPROC_LIB  = $(MPROC_DIR)/libmproc.a
COMM_LIB   = $(COMM_DIR)/libcomm.a

LIBS = $(DV_LIB)   \
       $(MESH_LIB) \
       $(FEM_LIB)  \
       $(NP_LIB)   \
       $(ALG_LIB)  \
       $(COMM_LIB) \
       $(MPROC_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = faces_gen.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)   \
            $(INCL_CMD)$(MESH_DIR) \
            $(INCL_CMD)$(FEM_DIR)  \
            $(INCL_CMD)$(NP_DIR)   \
            $(INCL_CMD)$(ALG_DIR)  \
            $(INCL_CMD)$(COMM_DIR) \
            $(INCL_CMD)$(MPROC_DIR) 

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)

#----------------------------------------
# MAIN RULES
#----------------------------------------
libface_gen.a: $(LIBS) $(OBJS)
	($(AR) libface_gen.a $(OBJS); ranlib libface_gen.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
faces_gen.o: $(LIBS) \
             faces_gen.f90 
	($(COMPILE) faces_gen.f90)


#-------------------------------------------------------------------------------
# CLEANING
#-------------------------------------------------------------------------------
.PHONY: cleanLoc cleanAll

cleanLoc:
	rm -f *.o *.mod *.exe *.a

cleanAll:
	@echo
	@echo ' --------------------------------------------------------------------------'
	@echo '   Cleaning objects and libs in PREREQUISITE folders'
	@echo ' --------------------------------------------------------------------------'
	@for dir in $(DIRS); do    \
	echo Cleaning $$dir;       \
	(cd $$dir; make cleanLoc); \
	done
	@echo
