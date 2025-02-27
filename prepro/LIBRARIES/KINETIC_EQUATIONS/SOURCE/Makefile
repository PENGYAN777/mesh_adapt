################################################################################
#
#  NODE-PAIR CODE   -   Kinetic Equations Library
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
THER_DIR  = ../../THERMODYNAMICS/SOURCE
DV_DIR    = ../../DYNAMIC_VECTORS/SOURCE
MESH_DIR  = ../../MESH_STRUCTURES/SOURCE
NP_DIR    = ../../NODEPAIR_STRUCTURES/SOURCE
NPGEN_DIR = ../../NODEPAIR_GENERATION/SOURCE
ALG_DIR   = ../../ALGEBRA/SOURCE
STRIO_DIR = ../../STRINGS_IO/SOURCE
SV_DIR    = ../../STRUCTURED_VARIABLES/SOURCE
MPROC_DIR = ../../MULTI_PROCESSOR/SOURCE
COMM_DIR  = ../../COMMONS/SOURCE
AXI_DIR   = ../../AXIAL_SYMMETRY/SOURCE
CSR_DIR   = ../../CSR/SOURCE



# Temporary using euler equations for bc
EUL_DIR = ../../EULER_EQUATIONS/SOURCE

DIRS = ./	    \
       $(THER_DIR)  \
       $(CSR_DIR)   \
       $(DV_DIR)    \
       $(NP_DIR)    \
       $(MESH_DIR)  \
       $(ALG_DIR)   \
       $(NPGEN_DIR) \
       $(EUL_DIR)   \
       $(STRIO_DIR) \
       $(SV_DIR)    \
       $(MPROC_DIR) \
       $(COMM_DIR)  \
       $(AXI_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
THER_LIB  = $(THER_DIR)/libthermodyn.a
CSR_LIB   = $(CSR_DIR)/libcsr.a
DV_LIB    = $(DV_DIR)/libdyn_vec.a
MESH_LIB  = $(MESH_DIR)/libmesh.a
NP_LIB    = $(NP_DIR)/libnp.a
NPGEN_LIB = $(NPGEN_DIR)/libnp_gen.a
ALG_LIB   = $(ALG_DIR)/libalgebra.a
STRIO_LIB = $(STRIO_DIR)/libio.a
SV_LIB    = $(SV_DIR)/libsv.a
MPROC_LIB = $(MPROC_DIR)/libmproc.a
COMM_LIB  = $(COMM_DIR)/libcomm.a
AXI_LIB   = $(AXI_DIR)/libaxis.a
# Temporary euler equations
EUL_LIB   = $(EUL_DIR)/libeuler.a
# System libraries
LAP_LIB   = /usr/lib/liblapack.a

LIBS = $(THER_LIB)  \
       $(CSR_LIB)   \
       $(DV_LIB)    \
       $(MESH_LIB)  \
       $(NP_LIB)    \
       $(LAP_LIB)   \
       $(ALG_LIB)   \
       $(NPGEN_LIB) \
       $(EUL_LIB)   \
       $(STRIO_LIB) \
       $(SV_LIB)    \
       $(MPROC_LIB) \
       $(COMM_LIB)  \
       $(AXI_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS  = derfc.o                        \
        kinetic_equations.o            \
        kinetic_extrap_boundary.o      \
        kinetic_boundary_cond.o        \
        kinetic_num_fluxes.o           \
        first_order_kinetic_fluxes.o   \
        high_resolution_kinetic_fluxes.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(THER_DIR)  \
            $(INCL_CMD)$(CSR_DIR)   \
            $(INCL_CMD)$(DV_DIR)    \
            $(INCL_CMD)$(NP_DIR)    \
            $(INCL_CMD)$(MESH_DIR)  \
            $(INCL_CMD)$(ALG_DIR)   \
            $(INCL_CMD)$(NPGEN_DIR) \
            $(INCL_CMD)$(EUL_DIR)   \
            $(INCL_CMD)$(STRIO_DIR) \
            $(INCL_CMD)$(SV_DIR)    \
	    $(INCL_CMD)$(MPROC_DIR) \
	    $(INCL_CMD)$(COMM_DIR)  \
	    $(INCL_CMD)$(AXI_DIR)  
	    
COMPILEF = $(CMP90) $(OPT90) $(INCL_DIRS)
COMPILEC = $(CC) $(COPT)



#----------------------------------------
# MAIN RULES
#----------------------------------------
libkin.a: $(LIBS) $(OBJS) 
	($(AR) libkin.a $(OBJS); ranlib libkin.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
derfc.o: derfc.c
	($(COMPILEC) derfc.c)

kinetic_equations.o: $(LIBS)\
                     kinetic_equations.f90\
                     derfc.o
	($(COMPILEF) kinetic_equations.f90 )

kinetic_extrap_boundary.o: $(LIBS)\
                           kinetic_extrap_boundary.f90
	($(COMPILEF) kinetic_extrap_boundary.f90 )

kinetic_boundary_cond.o: $(LIBS)\
                         kinetic_boundary_cond.f90\
                         kinetic_equations.o
	($(COMPILEF) kinetic_boundary_cond.f90 )

first_order_kinetic_fluxes.o: $(LIBS)\
                            first_order_kinetic_fluxes.f90\
                            kinetic_equations.o
	($(COMPILEF) first_order_kinetic_fluxes.f90 )

high_resolution_kinetic_fluxes.o: $(LIBS)\
                                  high_resolution_kinetic_fluxes.f90\
                                  kinetic_equations.o
	($(COMPILEF) high_resolution_kinetic_fluxes.f90 )

kinetic_num_fluxes.o: $(LIBS)\
                      kinetic_num_fluxes.f90\
                      kinetic_equations.o\
                      first_order_kinetic_fluxes.o\
                      high_resolution_kinetic_fluxes.o\
                      kinetic_extrap_boundary.o
	($(COMPILEF) kinetic_num_fluxes.f90 )



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
