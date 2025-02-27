################################################################################
#
#  NODE-PAIR CODE   -   Axial symmetry Library
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
DV_DIR    = ../../DYNAMIC_VECTORS/SOURCE
NP_DIR    = ../../NODEPAIR_STRUCTURES/SOURCE
CSR_DIR   = ../../CSR/SOURCE
EUL_DIR   = ../../EULER_EQUATIONS/SOURCE
MESH_DIR  = ../../MESH_STRUCTURES/SOURCE
MPROC_DIR = ../../MULTI_PROCESSOR/SOURCE
THER_DIR  = ../../THERMODYNAMICS/SOURCE
SV_DIR    = ../../STRUCTURED_VARIABLES/SOURCE
COMM_DIR  = ../../COMMONS/SOURCE

DIRS = ./           \
       $(NP_DIR)    \
       $(CSR_DIR)   \
       $(EUL_DIR)   \
       $(MESH_DIR)  \
       $(DV_DIR)    \
       $(MPROC_DIR) \
       $(THER_DIR)  \
       $(SV_DIR)    \
       $(COMM_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
NP_LIB    = $(NP_DIR)/libnp.a
CSR_LIB   = $(CSR_DIR)/libcsr.a
EUL_LIB   = $(EUL_DIR)/libeuler.a
MESH_LIB  = $(MESH_DIR)/libmesh.a
DV_LIB    = $(DV_DIR)/libdyn_vec.a
MPROC_LIB = $(MPROC_DIR)/libmproc.a
THER_LIB  = $(THER_DIR)/libthermodyn.a
SV_LIB    = $(SV_DIR)/libsv.a
COMM_LIB  = $(COMM_DIR)/libcomm.a

LIBS = $(NP_LIB)    \
       $(CSR_LIB)   \
       $(EUL_LIB)   \
       $(MESH_LIB)  \
       $(DV_LIB)    \
       $(MPROC_LIB) \
       $(THER_LIB)  \
       $(SV_LIB)    \
       $(COMM_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = axial_symmetry.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(NP_DIR)    \
            $(INCL_CMD)$(CSR_DIR)   \
            $(INCL_CMD)$(EUL_DIR)   \
	    $(INCL_CMD)$(MESH_DIR)  \
	    $(INCL_CMD)$(DV_DIR)    \
	    $(INCL_CMD)$(MPROC_DIR) \
	    $(INCL_CMD)$(THER_DIR)  \
	    $(INCL_CMD)$(SV_DIR)    \
	    $(INCL_CMD)$(COMM_DIR)
	    
COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS) 



#----------------------------------------
# MAIN RULES
#----------------------------------------
libaxis.a: $(LIBS) $(OBJS)
	($(AR) libaxis.a $(OBJS); ranlib libaxis.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
axial_symmetry.o: $(LIBS) axial_symmetry.f90
	($(COMPILE) axial_symmetry.f90)



#-------------------------------------------------------------------------------
# CLEANING
#-------------------------------------------------------------------------------
.PHONY: subs all clean cleanskit

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
