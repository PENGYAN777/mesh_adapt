################################################################################
#
#  NODE-PAIR CODE   -   Mesh Strucutre Library
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
MPROC_DIR = ../../MULTI_PROCESSOR/SOURCE
COMM_DIR  = ../../COMMONS/SOURCE

DIRS = $(DV_DIR)    \
       $(MPROC_DIR) \
       $(COMM_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
DV_LIB    = $(DV_DIR)/libdyn_vec.a
MPROC_LIB = $(MPROC_DIR)/libmproc.a
COMM_LIB  = $(COMM_DIR)/libcomm.a

LIBS = $(DV_LIB)    \
       $(MPROC_LIB) \
       $(COMM_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = nodes.o               \
       mesh_structure.o      \
       node_pair_structure.o \
       element_topology.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)    \
            $(INCL_CMD)$(MPROC_DIR) \
	    $(INCL_CMD)$(COMM_DIR)

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS) 



#----------------------------------------
# MAIN RULES
#----------------------------------------
libmesh.a: $(LIBS) $(OBJS)
	($(AR) libmesh.a $(OBJS); ranlib libmesh.a)

#----------------------------------------
# OBJECTS AND LIBRARIES
#----------------------------------------
mesh_structure.o: $(LIBS) mesh_structure.f90 \
                  element_topology.o
	($(COMPILE) mesh_structure.f90)

nodes.o: $(LIBS) nodes.f90
	($(COMPILE) nodes.f90)

node_pair_structure.o: $(LIBS) node_pair_structure.f90 \
                       mesh_structure.o
	($(COMPILE) node_pair_structure.f90)

element_topology.o: $(LIBS) element_topology.f90 
	($(COMPILE) element_topology.f90)



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
