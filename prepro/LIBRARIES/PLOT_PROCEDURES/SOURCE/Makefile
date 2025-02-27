################################################################################
#
#  NODE-PAIR CODE   -   Plot Procedures Library
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
DV_DIR = ../../DYNAMIC_VECTORS/SOURCE

DIRS = ./       \
       $(DV_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
DV_LIB = $(DV_DIR)/libdyn_vec.a

LIBS = $(DV_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = plot_procedures.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(DV_DIR)

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)


#----------------------------------------
# LIBRARIES AND EXECTUABLES
#----------------------------------------
libplot.a: $(LIBS) $(OBJS)
	($(AR) libplot.a $(OBJS); ranlib libplot.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
plot_procedures.o: $(LIBS) \
                   plot_procedures.f90
	($(COMPILE) plot_procedures.f90)



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
