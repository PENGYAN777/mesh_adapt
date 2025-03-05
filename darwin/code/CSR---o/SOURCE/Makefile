################################################################################
#
#  NODE-PAIR CODE   -   CSR Library
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
CMPL_CFG_DIR=../../

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
DV_DIR = ../../DYNAMIC_VECTORS/SOURCE

DIRS = ./        \
       $(DV_DIR)

#-------------------------------------------------------------------------------
# EXTERNAL LIBRARIES
#-------------------------------------------------------------------------------
DV_LIB = $(DV_DIR)/libdyn_vec.a

LIBS = $(DV_LIB) 

LIBS_DIRS = $(DV_DIR)

#-------------------------------------------------------------------------------
# LOCAL OBJECTS
#-------------------------------------------------------------------------------
OBJS  = csr.o              \
        coo.o              \
        csr_pair.o         \
        csr_pair_sys.o     \
        csr_dirichlet.o    \
        sparsekit_driver.o \
        umfpack_driver.o

#-------------------------------------------------------------------------------
# COMMANDS
#-------------------------------------------------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)
COMPILE   = $(CMP90) $(OPT90) $(INCL_DIRS)



#-------------------------------------------------------------------------------
# MAIN RULES
#-------------------------------------------------------------------------------
libcsr.a: $(LIBS) $(OBJS)
	($(AR) libcsr.a $(OBJS); ranlib libcsr.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
csr.o: csr.f90
	($(COMPILE) csr.f90)

coo.o: coo.f90 csr.o
	($(COMPILE) coo.f90)

csr_pair.o: $(LIBS) \
            csr_pair.f90 \
            csr.o
	($(COMPILE) csr_pair.f90)

csr_pair_sys.o: $(LIBS) \
                csr_pair_sys.f90 \
                csr.o
	($(COMPILE) csr_pair_sys.f90)

csr_dirichlet.o: $(LIBS) \
                 csr_dirichlet.f90 \
                 csr.o
	($(COMPILE) csr_dirichlet.f90)

sparsekit_driver.o: sparsekit_driver.f90 \
                    csr.o
	($(COMPILE) sparsekit_driver.f90)

umfpack_driver.o: umfpack_driver.f90 \
                  csr.o
	($(COMPILE) umfpack_driver.f90)



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
