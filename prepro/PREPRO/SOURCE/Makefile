################################################################################
#
#  NODE-PAIR CODE   -   PREPOCESSOR
#
#  Department of Aerospace Engineering
#  Politecnico di Milano. 1995-2010
#
################################################################################

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .f .f90 .o

.PHONY: libs

#-------------------------------------------------------------------------------
# COMPILER SETTINGS AND COMMANDS
#-------------------------------------------------------------------------------
#
# CMPL = VOID --->  No compile nor link commands are specified. The system 
#                   default commands are used. These are specified in the 
#                   system variable ENV.
# CMPL /= VOID ---> Commands to compile and link are explicitly specified
#                   through the included file Makefile.CMPL

VOID=
ifneq ($(CMPL),$(VOID))
  ifeq ($(DEBUG),yes)
    include Makefile.$(CMPL).DBG
  else
    include Makefile.$(CMPL)    
  endif
else
  include Makefile.ERR
endif

#-------------------------------------------------------------------------------
#  FOLDERS
#-------------------------------------------------------------------------------
DV_DIR     = ../../LIBRARIES/DYNAMIC_VECTORS/SOURCE
CSR_DIR    = ../../LIBRARIES/CSR/SOURCE
MESH_DIR   = ../../LIBRARIES/MESH_STRUCTURES/SOURCE
NP_DIR     = ../../LIBRARIES/NODEPAIR_STRUCTURES/SOURCE
NP_GEN_DIR = ../../LIBRARIES/NODEPAIR_GENERATION/SOURCE
FEM_DIR    = ../../LIBRARIES/FINITE_ELEMENT_METHOD/SOURCE
ALG_DIR    = ../../LIBRARIES/ALGEBRA/SOURCE
IO_DIR     = ../../LIBRARIES/STRINGS_IO/SOURCE
MESHP_DIR  = ../../LIBRARIES/MESH_PARTITION/SOURCE
MPROC_DIR  = ../../LIBRARIES/MULTI_PROCESSOR/SOURCE
COMM_DIR   = ../../LIBRARIES/COMMONS/SOURCE
FACES_DIR  = ../../LIBRARIES/FACES_GENERATION/SOURCE
PRECOMP_DIR  = ../../LIBRARIES/PRECOMPILED

DIRS = ./ $(DV_DIR)     $(CSR_DIR) $(MESH_DIR) $(NP_DIR) \
          $(NP_GEN_DIR) $(FEM_DIR) $(ALG_DIR)  $(IO_DIR) \
          $(MESHP_DIR)  $(MPROC_DIR) $(COMM_DIR) $(FACES_DIR)

#----------------------------------------
# EXTERNAL LIBRARIES
#----------------------------------------
DV_LIB     = $(DV_DIR)/libdyn_vec.a
CSR_LIB    = $(CSR_DIR)/libcsr.a
MESH_LIB   = $(MESH_DIR)/libmesh.a
NP_LIB     = $(NP_DIR)/libnp.a
NP_GEN_LIB = $(NP_GEN_DIR)/libnp_gen.a
FEM_LIB    = $(FEM_DIR)/libfem.a
ALG_LIB    = $(ALG_DIR)/libalgebra.a
IO_LIB     = $(IO_DIR)/libio.a
MESHP_LIB  = $(MESHP_DIR)/libmeshP.a
MPROC_LIB  = $(MPROC_DIR)/libmproc.a
COMM_LIB   = $(COMM_DIR)/libcomm.a
FACES_LIB  = $(FACES_DIR)/libface_gen.a
METIS_LIB  = $(PRECOMP_DIR)/libmetis.a

LIBS = $(DV_LIB)  $(CSR_LIB) $(ALG_LIB)      \
       $(FEM_LIB) $(IO_LIB)  $(METIS_LIB)    \
       $(MESH_LIB) $(NP_LIB) $(NP_GEN_LIB)   \
       $(MESHP_LIB) $(MPROC_LIB) $(COMM_LIB) \
       $(FACES_LIB)

LIBS_DIRS = $(DV_DIR)     $(CSR_DIR) $(COMM_DIR) $(ALG_DIR)  $(IO_DIR) \
            $(MPROC_DIR)  $(FEM_DIR) $(MESH_DIR) $(NP_DIR) \
            $(FACES_DIR) $(NP_GEN_DIR) $(MESHP_DIR) 

#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJS = 

#----------------------------------------
# COMMANDS
#----------------------------------------
INCL_DIRS =  $(INCL_CMD)$(DV_DIR)     \
             $(INCL_CMD)$(CSR_DIR)    \
             $(INCL_CMD)$(MESH_DIR)   \
             $(INCL_CMD)$(NP_DIR)     \
             $(INCL_CMD)$(NP_GEN_DIR) \
             $(INCL_CMD)$(FEM_DIR)    \
             $(INCL_CMD)$(ALG_DIR)    \
             $(INCL_CMD)$(IO_DIR)     \
             $(INCL_CMD)$(MPROC_DIR)  \
             $(INCL_CMD)$(COMM_DIR)   \
             $(INCL_CMD)$(MESHP_DIR)  \
             $(INCL_CMD)$(FACES_DIR)	     

LNK_LIBS   = $(LIBDIR_CMD)$(MESH_DIR)   $(LNKLIB_CMD)mesh     \
             $(LIBDIR_CMD)$(NP_DIR)     $(LNKLIB_CMD)np       \
             $(LIBDIR_CMD)$(NP_GEN_DIR) $(LNKLIB_CMD)np_gen   \
             $(LIBDIR_CMD)$(FEM_DIR)    $(LNKLIB_CMD)fem      \
             $(LIBDIR_CMD)$(ALG_DIR)    $(LNKLIB_CMD)algebra  \
             $(LIBDIR_CMD)$(CSR_DIR)    $(LNKLIB_CMD)csr      \
             $(LIBDIR_CMD)$(DV_DIR)     $(LNKLIB_CMD)dyn_vec  \
             $(LIBDIR_CMD)$(IO_DIR)     $(LNKLIB_CMD)io       \
             $(LIBDIR_CMD)$(MESHP_DIR)  $(LNKLIB_CMD)meshP    \
             $(LIBDIR_CMD)$(MPROC_DIR)  $(LNKLIB_CMD)mproc    \
             $(LIBDIR_CMD)$(COMM_DIR)   $(LNKLIB_CMD)comm     \
             $(LIBDIR_CMD)$(FACES_DIR)  $(LNKLIB_CMD)face_gen \
             $(LIBDIR_CMD)$(PRECOMP_DIR)  $(LNKLIB_CMD)metis

COMPILE = $(CMP90) $(OPT90) $(INCL_DIRS)
#LINK    = $(LNK) $(LNK_OPT)
LINK    = $(LNK) $(LNK_OPT) -lpthread


#----------------------------------------
# MAIN RULES
#----------------------------------------
all: msg libs prP   prePro.exe 

msg:
	@echo
	@echo ' --------------------------------------------------------------------------'
	@echo '   BUILD COMMANDS FOR THE NODE-PAIR CODE: PRE-PROCESSOR'
	@echo '   copyright 1995-2010 - Alberto Guardone, Marco Fossati'
	@echo ' --------------------------------------------------------------------------'
	@echo
	@echo '   - Libraries:'

libs:
	@for dir in $(LIBS_DIRS) ;\
	do \
	(echo; echo 'MAKING ' $$dir; cd $$dir; $(MAKE) CMPL=$(CMPL) DEBUG=$(DEBUG))  ;\
	done
	@echo

prP:
	@echo
	@echo '   - PreProcessor:'

 
prePro.exe: $(LIBS) $(OBJS) prePro.o 
	$(LINK) -o prePro.exe prePro.o $(OBJS) $(LNK_LIBS)

 
#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
prePro.o:  $(LIBS) $(OBJS) prePro.f90
	($(COMPILE) prePro.f90 )

#-------------------------------------------------------------------------------
#  CLEANING COMMANDS
#-------------------------------------------------------------------------------
.PHONY: cleanLoc cleanAll

cleanLoc:
	rm -f *.o *.mod

cleanAll:
	@echo
	@echo ' --------------------------------------------------------------------------'
	@echo '   Cleaning objects, libs and exes in PRE-PROCESSOR dependency folders'
	@echo ' --------------------------------------------------------------------------'
	@for dir in $(DIRS) ;          \
	do       \
	echo ''; echo 'Cleaning' $$dir ;          \
	(cd $$dir ; $(MAKE) cleanLoc) ; \
	done
	@echo
