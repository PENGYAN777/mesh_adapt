################################################################################
#
#  NODE-PAIR CODE   -   Finite Element Method Library
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
CSR_DIR   = ../../CSR/SOURCE
ALG_DIR   = ../../ALGEBRA/SOURCE

DIRS = ./         \
       $(DV_DIR)  \
       $(ALG_DIR) \
       $(CSR_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
DV_LIB   = $(DV_DIR)/libdyn_vec.a
CSR_LIB  = $(CSR_DIR)/libcsr.a
ALG_LIB  = $(ALG_DIR)/libalgebra.a

LIBS = $(DV_LIB)  \
       $(CSR_LIB) \
       $(ALG_LIB)
	    
#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = fem_ele_types.o     \
       csr_elem.o          \
       fem_ref_elements.o  \
       fem_gauss_points.o  \
       fem_quadrature.o    \
       fem_quadrature_tg.o \
       fem_norms.o         \
       fem_inlet.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)  \
            $(INCL_CMD)$(CSR_DIR) \
            $(INCL_CMD)$(ALG_DIR)

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)



#----------------------------------------
# MAIN RULES
#----------------------------------------
libfem.a: $(LIBS) $(OBJS)
	($(AR) libfem.a $(OBJS); ranlib libfem.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
fem_ele_types.o: fem_ele_types.f90
	( $(COMPILE) fem_ele_types.f90)

csr_elem.o: csr_elem.f90 fem_ele_types.o 
	($(COMPILE) csr_elem.f90 )
  
fem_ref_elements.o: fem_ref_elements.f90
	( $(COMPILE) fem_ref_elements.f90)

fem_gauss_points.o: fem_gauss_points.f90 \
                    fem_ele_types.o \
                    fem_ref_elements.o \
                    $(LIBS)
	( $(COMPILE) fem_gauss_points.f90)

fem_quadrature.o: fem_quadrature.f90 \
                  fem_ele_types.o \
                  fem_gauss_points.o \
                  $(LIBS)
	( $(COMPILE) fem_quadrature.f90)

fem_quadrature_tg.o: fem_quadrature_tg.f90 \
                     fem_ele_types.o \
                     $(LIBS)
	( $(COMPILE) fem_quadrature_tg.f90)

fem_norms.o: fem_norms.f90 \
             fem_ele_types.o \
             $(LIBS)
	( $(COMPILE) fem_norms.f90)

fem_inlet.o: fem_inlet.f90
	( $(COMPILE) fem_inlet.f90)



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
