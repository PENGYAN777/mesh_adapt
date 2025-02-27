#
#------------------------------------------------------------------------------
#     Definitions
#------------------------------------------------------------------------------
#
SYS=$(shell uname)_$(shell uname -m)
#
gfortran:
#
	@echo
	@echo "---------------------------------------------------------"
	@echo "Configuring for ${SYS}..."
	@echo "CMP_F90 = gfortran -c" > Makefile.config
	@echo "LNK_F90 = gfortran " >> Makefile.config
#	@echo "FFLAGS = -O3" >> Makefile.config
	@echo "FFLAGS = " >> Makefile.config
	@echo "LIBS = ../libs/liblapack.a ../libs/libskit.a ../libs/libblas.a" >> Makefile.config
	@echo "INCL = -lstdc++" >> Makefile.config
#
	@echo "Building DARWIN for ${SYS}..."
	@echo "---------------------------------------------------------"
	@echo
	@cd code/ ; $(MAKE)
	@echo
#
gfortran_dbg:
#
	@echo
	@echo "---------------------------------------------------------"
	@echo "Configuring for ${SYS}..."
	@echo "CMP_F90 = gfortran -c" > Makefile.config
	@echo "LNK_F90 = gfortran " >> Makefile.config
	@echo "FFLAGS = -g -Wall -fbacktrace" >> Makefile.config
	@echo "LIBS = ../libs/liblapack.a ../libs/libskit.a ../libs/libblas.a" >> Makefile.config
	@echo "INCL = -lstdc++" >> Makefile.config
#
	@echo "Building DARWIN for ${SYS}..."
	@echo "---------------------------------------------------------"
	@echo
	@cd code/ ; $(MAKE)
	@echo
#
ifort:
	@echo
	@echo "---------------------------------------------------------"
	@echo "Configuring for ${SYS}..."
	@echo "CMP_F90 = ifort -c -heap-arrays 100" > Makefile.config
	@echo "LNK_F90 = ifort" >> Makefile.config
	@echo "FFLAGS = -Bstatic -DCONVERT_BIG -fpp -DMAIN_UCASE -O3 -traceback -check bounds" >> Makefile.config
	@echo "LIBS = ../libs/liblapack.a ../libs/libskit.a ../libs/libblas.a ../libs/libkdtree2.a" >> Makefile.config
	@echo "INCL = -lg2c -lstdc++ -lgfortran -Ikdtree2_mod/" >> Makefile.config
#
	@echo "Building DARWIN for ${SYS}..."
	@echo "---------------------------------------------------------"
	@echo
	@cd code/ ; $(MAKE)
	@echo
#
#------------------------------------------------------------------------------
#     Cleaning definitions
#------------------------------------------------------------------------------
#
clean:
	@echo
	@echo "---------------------------------------------------------"
	@echo "Cleaning DARWIN ..."
	@echo
	@cd code/ ; $(MAKE) clean
	rm -f Makefile.config
	@echo "---------------------------------------------------------"
	@echo
