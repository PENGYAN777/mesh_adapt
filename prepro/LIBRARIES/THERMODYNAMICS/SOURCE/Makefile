################################################################################
#
#  NODE-PAIR CODE   -   Thermodynamics Library
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
SV_DIR   = ../../STRUCTURED_VARIABLES/SOURCE
ALG_DIR  = ../../ALGEBRA/SOURCE

DIRS = ./         \
       $(ALG_DIR) \
       $(SV_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
ALG_LIB = $(ALG_DIR)/libalgebra.a
SV_LIB  = $(SV_DIR)/libsv.a

LIBS = $(ALG_LIB) \
       $(SV_LIB)

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS  = vapor_properties.o        \
        thermodynamics.o	  \
        ideal_gas.o		  \
        van_der_waals_gas.o	  \
        martin_hou_gas.o	  \
        soave_redlich_kwong_gas.o \
        peng_robinson_gas.o	  \
        redlich_kwong_gas.o	  \
        clausius_II_gas.o	  \
        ideal_specific_heat.o	  \
        transport_properties.o    \
        gas_properties.o

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS = $(INCL_CMD)$(ALG_DIR) \
            $(INCL_CMD)$(SV_DIR)

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)



#----------------------------------------
# MAIN RULES
#----------------------------------------
libthermodyn.a: $(LIBS) $(OBJS)
	($(AR) libthermodyn.a $(OBJS); ranlib libthermodyn.a)

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
vapor_properties.o: $(LIBS) \
                    vapor_properties.f90 \
                    thermodynamics.o \
                    transport_properties.o 
	($(COMPILE) vapor_properties.f90 )

thermodynamics.o: thermodynamics.f90 \
                  ideal_gas.o \
                  van_der_waals_gas.o \
                  martin_hou_gas.o \
                  soave_redlich_kwong_gas.o \
                  peng_robinson_gas.o \
                  redlich_kwong_gas.o \
                  clausius_II_gas.o   \
                  gas_properties.o \
                  ideal_specific_heat.o 
	($(COMPILE) thermodynamics.f90 )

ideal_gas.o: ideal_gas.f90 \
             gas_properties.o \
             ideal_specific_heat.o 
	($(COMPILE) ideal_gas.f90 )

van_der_waals_gas.o: van_der_waals_gas.f90 \
                     gas_properties.o \
                     ideal_specific_heat.o 
	($(COMPILE) van_der_waals_gas.f90 )

martin_hou_gas.o: martin_hou_gas.f90 \
                  gas_properties.o \
                  ideal_specific_heat.o 
	($(COMPILE) martin_hou_gas.f90 )

soave_redlich_kwong_gas.o: soave_redlich_kwong_gas.f90 \
                  gas_properties.o \
                  ideal_specific_heat.o 
	($(COMPILE) soave_redlich_kwong_gas.f90 )

peng_robinson_gas.o: peng_robinson_gas.f90 \
                  gas_properties.o \
                  ideal_specific_heat.o 
	($(COMPILE) peng_robinson_gas.f90 )

redlich_kwong_gas.o: redlich_kwong_gas.f90 \
                  gas_properties.o \
                  ideal_specific_heat.o 
	($(COMPILE) redlich_kwong_gas.f90 )

clausius_II_gas.o: clausius_II_gas.f90 \
                  gas_properties.o \
                  ideal_specific_heat.o 
	($(COMPILE) clausius_II_gas.f90 )

ideal_specific_heat.o: ideal_specific_heat.f90 \
                       gas_properties.o
	($(COMPILE) ideal_specific_heat.f90 )

gas_properties.o: gas_properties.f90 
	($(COMPILE) gas_properties.f90 )

transport_properties.o: transport_properties.f90 \
                        thermodynamics.o \
                        gas_properties.o
	($(COMPILE) transport_properties.f90 )



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
