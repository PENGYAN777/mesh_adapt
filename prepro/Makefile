################################################################################
#
#  MANAGE PROJECT FILES
#
#  department of Aerospace Engineering
#  Politecnico di Milano. 1995-2010
#
################################################################################

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .f .f90 .o

#-------------------------------------------------------------------------------
#  LIBRARIES FOLDERS
#-------------------------------------------------------------------------------
ALG_DIR    = LIBRARIES/ALGEBRA
CSR_DIR    = LIBRARIES/CSR
DV_DIR     = LIBRARIES/DYNAMIC_VECTORS
EUL_DIR    = LIBRARIES/EULER_EQUATIONS
FEM_DIR    = LIBRARIES/FINITE_ELEMENT_METHOD
KIN_DIR    = LIBRARIES/KINETIC_EQUATIONS
MPART_DIR  = LIBRARIES/MESH_PARTITION
MESH_DIR   = LIBRARIES/MESH_STRUCTURES
MPROC_DIR  = LIBRARIES/MULTI_PROCESSOR
NS_DIR     = LIBRARIES/NAVIER_STOKES_EQUATIONS
NPGEN_DIR  = LIBRARIES/NODEPAIR_GENERATION
NP_DIR     = LIBRARIES/NODEPAIR_STRUCTURES
PLT_DIR    = LIBRARIES/PLOT_PROCEDURES
SCALE_DIR  = LIBRARIES/SCALING
SKIT_DIR   = LIBRARIES/SPARSKIT/
IO_DIR     = LIBRARIES/STRINGS_IO
SV_DIR     = LIBRARIES/STRUCTURED_VARIABLES
THER_DIR   = LIBRARIES/THERMODYNAMICS
UMFP_DIR   = LIBRARIES/UMFPACK/

#-------------------------------------------------------------------------------
#  SOURCE FOLDERS
#-------------------------------------------------------------------------------
ALG_SRC_DIR    = LIBRARIES/ALGEBRA/SOURCE
CSR_SRC_DIR    = LIBRARIES/CSR/SOURCE
DV_SRC_DIR     = LIBRARIES/DYNAMIC_VECTORS/SOURCE
EUL_SRC_DIR    = LIBRARIES/EULER_EQUATIONS/SOURCE
FEM_SRC_DIR    = LIBRARIES/FINITE_ELEMENT_METHOD/SOURCE
KIN_SRC_DIR    = LIBRARIES/KINETIC_EQUATIONS/SOURCE
MPART_SRC_DIR  = LIBRARIES/MESH_PARTITION/SOURCE
MESH_SRC_DIR   = LIBRARIES/MESH_STRUCTURES/SOURCE
MPROC_SRC_DIR  = LIBRARIES/MULTI_PROCESSOR/SOURCE
NS_SRC_DIR     = LIBRARIES/NAVIER_STOKES_EQUATIONS/SOURCE
NPGEN_SRC_DIR  = LIBRARIES/NODEPAIR_GENERATION/SOURCE
NP_SRC_DIR     = LIBRARIES/NODEPAIR_STRUCTURES/SOURCE
PLT_SRC_DIR    = LIBRARIES/PLOT_PROCEDURES/SOURCE
SCALE_SRC_DIR  = LIBRARIES/SCALING/SOURCE
SKIT_SRC_DIR   = LIBRARIES/SPARSKIT/
IO_SRC_DIR     = LIBRARIES/STRINGS_IO/SOURCE
SV_SRC_DIR     = LIBRARIES/STRUCTURED_VARIABLES/SOURCE
THER_SRC_DIR   = LIBRARIES/THERMODYNAMICS/SOURCE
UMFP_SRC_DIR   = LIBRARIES/UMFPACK/

#-------------------------------------------------------------------------------
#  MAIN RULES
#-------------------------------------------------------------------------------
.PHONY: msg openMakes

all: msg

msg:
	@echo
	@echo ' -------------------------------------------------------------------------------------'
	@echo '   MANAGE PROJECT`S FILES'
	@echo ' -------------------------------------------------------------------------------------'
	@echo
	@echo ' type:   make openMakes                to open ALL project`s Makefiles'
	@echo ' type:   make openSources              to open ALL project`s Fortran90 source files'
	@echo ' type:   make LIB=library  openLib     to open only LIB Makefile and Fortran90 source files'
	@echo ' type:   make rmDirCont CONT=file_name to recursively remove file_name from project`s folders'
	@echo ' type:   make rmSrcCont CONT=file_name to recursively remove file_name from project`s source folders'
	@echo ''

openMakes:
	(nedit $(ALG_SRC_DIR)/Makefile \
	$(CSR_SRC_DIR)/Makefile \
	$(DV_SRC_DIR)/Makefile \
	$(EUL_SRC_DIR)/Makefile \
	$(FEM_SRC_DIR)/Makefile \
	$(KIN_SRC_DIR)/Makefile \
	$(MPART_SRC_DIR)/Makefile \
	$(MESH_SRC_DIR)/Makefile \
	$(MPROC_SRC_DIR)/Makefile \
	$(NS_SRC_DIR)/Makefile \
	$(NPGEN_SRC_DIR)/Makefile \
	$(NP_SRC_DIR)/Makefile \
	$(PLT_SRC_DIR)/Makefile \
	$(SCALE_SRC_DIR)/Makefile \
	$(SKIT_SRC_DIR)/Makefile \
	$(IO_SRC_DIR)/Makefile \
	$(SV_SRC_DIR)/Makefile \
	$(THER_SRC_DIR)/Makefile \
	$(UMFP_SRC_DIR)/Makefile \
	SOLVER/SOURCE/Makefile \
	PREPRO/SOURCE/Makefile &)
	@echo


openSources:
	(nedit $(ALG_SRC_DIR)/*.f90 \
	$(CSR_SRC_DIR)/*.f90 \
	$(DV_SRC_DIR)/*.f90 \
	$(EUL_SRC_DIR)/*.f90 \
	$(FEM_SRC_DIR)/*.f90 \
	$(KIN_SRC_DIR)/*.f90 \
	$(MPART_SRC_DIR)/*.f90 \
	$(MESH_SRC_DIR)/*.f90 \
	$(MPROC_SRC_DIR)/*.f90 \
	$(NS_SRC_DIR)/*.f90 \
	$(NPGEN_SRC_DIR)/*.f90 \
	$(NP_SRC_DIR)/*.f90 \
	$(PLT_SRC_DIR)/*.f90 \
	$(SCALE_SRC_DIR)/*.f90 \
	$(IO_SRC_DIR)/*.f90 \
	$(SV_SRC_DIR)/*.f90 \
	$(THER_SRC_DIR)/*.f90 \
	SOLVER/SOURCE/*.f90 \
	PREPRO/SOURCE/*.f90 &)
	@echo

openLib:
	(nedit LIBRARIES/$(LIB)/SOURCE/*.f90 \
	LIBRARIES/$(LIB)/SOURCE/Makefile)

rmDirCont:
	rm -rf $(ALG_DIR)/$(CONT) \
	$(CSR_DIR)/$(CONT) \
	$(DV_DIR)/$(CONT) \
	$(EUL_DIR)/$(CONT) \
	$(FEM_DIR)/$(CONT) \
	$(KIN_DIR)/$(CONT) \
	$(MPART_DIR)/$(CONT) \
	$(MESH_DIR)/$(CONT) \
	$(MPROC_DIR)/$(CONT) \
	$(NS_DIR)/$(CONT) \
	$(NPGEN_DIR)/$(CONT) \
	$(NP_DIR)/$(CONT) \
	$(PLT_DIR)/$(CONT) \
	$(SCALE_DIR)/$(CONT) \
	$(IO_DIR)/$(CONT) \
	$(SV_DIR)/$(CONT) \
	$(THER_DIR)/$(CONT) \
	SOLVER/SOURCE/$(CONT) \
	PREPRO/SOURCE/$(CONT)
	@echo

rmSrcCont:
	rm -rf $(ALG_SRC_DIR)/$(CONT) \
	$(CSR_SRC_DIR)/$(CONT) \
	$(DV_SRC_DIR)/$(CONT) \
	$(EUL_SRC_DIR)/$(CONT) \
	$(FEM_SRC_DIR)/$(CONT) \
	$(KIN_SRC_DIR)/$(CONT) \
	$(MPART_SRC_DIR)/$(CONT) \
	$(MESH_SRC_DIR)/$(CONT) \
	$(MPROC_SRC_DIR)/$(CONT) \
	$(NS_SRC_DIR)/$(CONT) \
	$(NPGEN_SRC_DIR)/$(CONT) \
	$(NP_SRC_DIR)/$(CONT) \
	$(PLT_SRC_DIR)/$(CONT) \
	$(SCALE_SRC_DIR)/$(CONT) \
	$(IO_SRC_DIR)/$(CONT) \
	$(SV_SRC_DIR)/$(CONT) \
	$(THER_SRC_DIR)/$(CONT) \
	SOLVER/SOURCE/$(CONT) \
	PREPRO/SOURCE/$(CONT)
	@echo
