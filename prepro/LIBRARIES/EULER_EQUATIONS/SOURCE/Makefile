################################################################################
#
#  NODE-PAIR CODE   -   Euler Equations Library
#
#  department of Aerospace Engineering
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
CSR_DIR   = ../../CSR/SOURCE
MESH_DIR  = ../../MESH_STRUCTURES/SOURCE
NP_DIR    = ../../NODEPAIR_STRUCTURES/SOURCE
ALG_DIR   = ../../ALGEBRA/SOURCE
SV_DIR    = ../../STRUCTURED_VARIABLES/SOURCE
MPROC_DIR = ../../MULTI_PROCESSOR/SOURCE
COMM_DIR  = ../../COMMONS/SOURCE

DIRS = ./	    \
       $(THER_DIR)  \
       $(ALG_DIR)   \
       $(DV_DIR)    \
       $(MESH_DIR)  \
       $(CSR_DIR)   \
       $(NP_DIR)    \
       $(SV_DIR)    \
       $(MPROC_DIR) \
       $(COMM_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
THER_LIB  = $(THER_DIR)/libthermodyn.a
ALG_LIB   = $(ALG_DIR)/libalgebra.a
DV_LIB    = $(DV_DIR)/libdyn_vec.a
MESH_LIB  = $(MESH_DIR)/libmesh.a
NP_LIB    = $(NP_DIR)/libnp.a
CSR_LIB   = $(CSR_DIR)/libcsr.a
SV_LIB    = $(SV_DIR)/libsv.a
MPROC_LIB = $(MPROC_DIR)/libmproc.a
COMM_LIB  = $(COMM_DIR)/libcomm.a 

LIBS = $(THER_LIB)  \
       $(ALG_LIB)   \
       $(DV_LIB)    \
       $(MESH_LIB)  \
       $(NP_LIB)    \
       $(CSR_LIB)   \
       $(SV_LIB)    \
       $(MPROC_LIB) \
       $(COMM_LIB)

#----------------------------------------
#  LOCAL OBJECTS
#----------------------------------------
OBJS  = euler_equations.o        \
        euler_num_fluxes.o	 \
        euler_boundary_cond.o	 \
        first_order_fluxes.o	 \
        high_resolution_fluxes.o \
        muscl_fluxes.o  	 \
        extrap_boundary.o	 \
        roe_average.o

#-------------------------------------------------------------------------------
# COMMANDS
#-------------------------------------------------------------------------------
INCL_DIRS = $(INCL_CMD)$(THER_DIR)  \
            $(INCL_CMD)$(ALG_DIR)   \
            $(INCL_CMD)$(DV_DIR)    \
            $(INCL_CMD)$(NP_DIR)    \
            $(INCL_CMD)$(MESH_DIR)  \
            $(INCL_CMD)$(CSR_DIR)   \
            $(INCL_CMD)$(SV_DIR)    \
            $(INCL_CMD)$(MPROC_DIR) \
	    $(INCL_CMD)$(COMM_DIR)
	   
COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)



#----------------------------------------
# MAIN RULES
#----------------------------------------
libeuler.a: $(LIBS) $(OBJS)
	($(AR) libeuler.a $(OBJS); ranlib libeuler.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
euler_num_fluxes.o: $(LIBS) \
                    euler_num_fluxes.f90\
                    euler_equations.o\
                    first_order_fluxes.o\
                    high_resolution_fluxes.o\
                    muscl_fluxes.o\
                    extrap_boundary.o
	($(COMPILE) euler_num_fluxes.f90 )

roe_average.o: $(LIBS) \
               roe_average.f90\
               euler_equations.o
	($(COMPILE) roe_average.f90 )

first_order_fluxes.o: $(LIBS) \
                      first_order_fluxes.f90\
                      roe_average.o\
                      euler_equations.o
	($(COMPILE) first_order_fluxes.f90 )

high_resolution_fluxes.o: $(LIBS) \
                          high_resolution_fluxes.f90\
                          roe_average.o\
                          euler_equations.o
	($(COMPILE) high_resolution_fluxes.f90 )

euler_boundary_cond.o: $(LIBS) \
                       euler_boundary_cond.f90\
                       euler_equations.o
	($(COMPILE) euler_boundary_cond.f90 )

euler_equations.o: $(LIBS) \
                   euler_equations.f90
	($(COMPILE) euler_equations.f90 )

muscl_fluxes.o: $(LIBS) \
                muscl_fluxes.f90\
                euler_equations.o\
                roe_average.o\
                high_resolution_fluxes.o
	($(COMPILE) muscl_fluxes.f90 )

extrap_boundary.o: $(LIBS) \
                   extrap_boundary.f90\
                   euler_equations.o\
                   roe_average.o
	($(COMPILE) extrap_boundary.f90 )



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
