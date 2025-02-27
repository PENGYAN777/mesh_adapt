################################################################################
#
#  NODE-PAIR CODE   -   Node-Pair Structures Library
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
MESH_DIR  = ../../MESH_STRUCTURES/SOURCE
MPROC_DIR = ../../MULTI_PROCESSOR/SOURCE
COMM_DIR  = ../../COMMONS/SOURCE

DIRS = ./           \
       $(DV_DIR)    \
       $(MESH_DIR)  \
       $(MPROC_DIR) \
       $(COMM_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
DV_LIB    = $(DV_DIR)/libdyn_vec.a
MESH_LIB  = $(MESH_DIR)/libmesh.a
MPROC_LIB = $(MPROC_DIR)/libmproc.a
COMM_LIB  = $(COMM_DIR)/libcomm.a

LIBS = $(DV_LIB)    \
       $(MESH_LIB)  \
       $(MPROC_LIB) \
       $(COMM_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS =  np_metric_coeff.o \
        np_quadrature.o   \
        np_quadrature_tg.o 

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)    \
            $(INCL_CMD)$(MESH_DIR)  \
            $(INCL_CMD)$(MPROC_DIR) \
	    $(INCL_CMD)$(COMM_DIR)

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)


#----------------------------------------
# EXECUTABLES
#----------------------------------------
libnp.a: $(LIBS) $(OBJS)
	($(AR) libnp.a $(OBJS); ranlib libnp.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
np_metric_coeff.o: $(MESH_LIB) np_metric_coeff.f90
	($(COMPILE) np_metric_coeff.f90)

np_quadrature.o: np_quadrature.f90
	($(COMPILE) np_quadrature.f90)
	
np_quadrature_tg.o: np_quadrature_tg.f90
	($(COMPILE) np_quadrature_tg.f90)



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
