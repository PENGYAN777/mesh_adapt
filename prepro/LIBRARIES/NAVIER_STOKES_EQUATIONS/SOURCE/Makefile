################################################################################
#
#  NODE-PAIR CODE   -   Navier-Stokes Equations Library
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
EUL_DIR   = ../../EULER_EQUATIONS/SOURCE
THER_DIR  = ../../THERMODYNAMICS/SOURCE
DV_DIR    = ../../DYNAMIC_VECTORS/SOURCE
MESH_DIR  = ../../MESH_STRUCTURES/SOURCE
NP_DIR    = ../../NODEPAIR_STRUCTURES/SOURCE
SV_DIR    = ../../STRUCTURED_VARIABLES/SOURCE
MPROC_DIR = ../../MULTI_PROCESSOR/SOURCE
COMM_DIR  = ../../COMMONS/SOURCE

DIRS = ./           \
       $(DV_DIR)    \
       $(NP_DIR)    \
       $(EUL_DIR)   \
       $(THER_DIR)  \
       $(SV_DIR)    \
       $(MPROC_DIR) \
       $(MESH_DIR)  \
       $(COMM_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
THER_LIB  = $(THER_DIR)/libthermodyn.a
DV_LIB    = $(DV_DIR)/libdyn_vec.a
EUL_LIB   = $(EUL_DIR)/libeuler.a
SV_LIB    = $(SV_DIR)/libsv.a
MPROC_LIB = $(MPROC_DIR)/libmproc.a
NP_LIB    = $(NP_DIR)/libnp.a
MESH_LIB  = $(MESH_DIR)/libmesh.a
COMM_LIB  = $(COMM_DIR)/libcomm.a

LIBS = $(THER_LIB)  \
       $(SV_LIB)    \
       $(DV_LIB)    \
       $(EUL_LIB)   \
       $(MPROC_LIB) \
       $(NP_LIB)    \
       $(MESH_LIB)  \
       $(COMM_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS  = ns_viscous_bc.o  \
        ns_thermal_bc.o  \
        ns_num_fluxes.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)    \
            $(INCL_CMD)$(NP_DIR)    \
	    $(INCL_CMD)$(EUL_DIR)   \
	    $(INCL_CMD)$(THER_DIR)  \
            $(INCL_CMD)$(MESH_DIR)  \
	    $(INCL_CMD)$(SV_DIR)    \
	    $(INCL_CMD)$(MPROC_DIR) \
            $(INCL_CMD)$(COMM_DIR)
	    
COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS) 


#----------------------------------------
# LIBRARIES AND EXECTUABLES
#----------------------------------------
libns.a: $(LIBS) $(OBJS)
	($(AR) libns.a $(OBJS); ranlib libns.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
ns_viscous_bc.o: ns_viscous_bc.f90 \
                 $(DV_LIB) 
	($(COMPILE) ns_viscous_bc.f90 )

ns_thermal_bc.o: ns_thermal_bc.f90 \
                 $(DV_LIB) \
                 $(EUL_LIB) 
	($(COMPILE) ns_thermal_bc.f90 )

ns_num_fluxes.o: ns_num_fluxes.f90 \
                 $(NP_DIR)/np_quadrature.o \
                 $(EUL_LIB)
	($(COMPILE) ns_num_fluxes.f90 )



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
