################################################################################
#
#  NODE-PAIR CODE   -   Structured Variables Library
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
DIRS = ./

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
LIBS =
LIBS_DIR = 

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS  = structures.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = 
COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)



#----------------------------------------
# MAIN RULES
#----------------------------------------
libsv.a: $(LIBS) $(OBJS)
	($(AR) libsv.a $(OBJS); ranlib libsv.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
structures.o: structures.f90
	($(COMPILE) structures.f90)



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
