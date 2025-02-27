################################################################################
#
#  NODE-PAIR CODE   -   Strings Input/Output Library
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
DIRS = ./

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = messages.o \
       strings.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS =
COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)



#----------------------------------------
# MAIN RULES
#----------------------------------------
libio.a: $(OBJS)
	($(AR) libio.a $(OBJS); ranlib libio.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
messages.o:  messages.f90  
	($(COMPILE) messages.f90)

strings.o:  strings.f90
	($(COMPILE) strings.f90)



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
