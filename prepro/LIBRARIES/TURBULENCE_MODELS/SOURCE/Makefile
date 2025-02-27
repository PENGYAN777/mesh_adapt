################################################################################
#
#  NODE-PAIR CODE   -   Turbulence models Library
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

#----------------------------------------
# FOLDERS
#----------------------------------------
EUL_DIR   = ../../EULER_EQUATIONS/SOURCE
DV_DIR    = ../../DYNAMIC_VECTORS/SOURCE
MESH_DIR  = ../../MESH_STRUCTURES/SOURCE
NP_DIR    = ../../NODEPAIR_STRUCTURES/SOURCE
THER_DIR  = ../../THERMODYNAMICS/SOURCE
MPROC_DIR = ../../MULTI_PROCESSOR/SOURCE
SV_DIR    = ../../STRUCTURED_VARIABLES/SOURCE
COMM_DIR = ../../COMMONS/SOURCE

DIRS = ./	    \
       $(DV_DIR)    \
       $(NP_DIR)    \
       $(EUL_DIR)   \
       $(THER_DIR)  \
       $(MPROC_DIR) \
       $(SV_DIR)    \
       $(COMM_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
DV_LIB    = $(DV_DIR)/libdyn_vec.a
EUL_LIB   = $(EUL_DIR)/libeuler.a
THER_LIB  = $(THER_DIR)/libthermodyn.a
MPROC_LIB = $(MPROC_DIR)/libmproc.a
SV_LIB    = $(SV_DIR)/libsv.a
COMM_LIB = $(COMM_DIR)/libcomm.a

LIBS = $(DV_DIR)    \
       $(NS_LIB)    \
       $(EUL_LIB)   \
       $(THER_LIB)  \
       $(MPROC_LIB) \
       $(SV_LIB)    \
       $(COMM_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = turbulence_models.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)    \
            $(INCL_CMD)$(NP_DIR)    \
	    $(INCL_CMD)$(EUL_DIR)   \
            $(INCL_CMD)$(MESH_DIR)  \
	    $(INCL_CMD)$(THER_DIR)  \
	    $(INCL_CMD)$(MPROC_DIR) \
	    $(INCL_CMD)$(SV_DIR)    \
            $(INCL_CMD)$(COMM_DIR)

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)



#----------------------------------------
# MAIN RULES
#----------------------------------------
libturb.a: $(LIBS) $(OBJS)
	($(AR) libturb.a $(OBJS); ranlib libturb.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
turbulence_models.o: $(LIBS) turbulence_models.f90
	($(COMPILE) turbulence_models.f90)



#----------------------------------------
# CLEANING
#----------------------------------------
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
