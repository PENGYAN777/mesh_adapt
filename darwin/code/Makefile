#
SHELL = /bin/sh
UNAME = /bin/uname
OPERSYS = $(shell $(UNAME))
ARCH = $(shell $(UNAME) -m)
#
#----------------------------------------
# COMPILER SETTINGS AND COMMANDS
#----------------------------------------
#
include ../Makefile.config
COMPILE = $(CMP_F90) $(FFLAGS) $(INCL)
LINK = $(LNK_F90) $(FFLAGS)
EXETAG = $(OPERSYS)_$(ARCH)
#
#-------------------------------------
#  SET DIRECTORIES
#-------------------------------------
BINDIR = .
LIBDIR = .
OBJDIR = .
SRCDIR = .

#-------------------------------------
#  LIBRARIES
#-------------------------------------
#LIBBLAS = /usr/lib64/libblas.a
#LIBS = $(LIBBLAS)

#----------------------------------------
# OBJECTS
#----------------------------------------
OBJ = dynamic_vector.o\
      lin_algebra.o\
      csr.o \
      solve_skit.o \
      csr_pair_sys.o \
      np_quadrature.o\
      select_strategy.o\
      nodes.o\
      element_topology.o\
      structures.o\
      boundary_geometry.o\
      mesh_structure.o\
      node_pair_str.o\
      metric_coefficients.o\
      np_topology_gen.o\
      lib_geo.o\
      grid_utils.o\
      derivatives.o\
      thermodynamic.o\
      io_proc.o\
      solution_interpolation.o\
      error_estimator.o\
      marking_proc.o\
      rh_setting.o\
      new_nodes.o\
      new_elements.o\
      triangulate_proc.o\
      reconstruction_proc.o \
      grid_adaption.o\
      darwin.o

all: darwin.exe move 

darwin.exe: darwin.o $(OBJ)
	($(LINK) -o darwin.$(EXETAG) $(OBJ) $(LIBS))
move:
	mv darwin.$(EXETAG) ../bins        

#----------------------------------------
#  OBJECTS
#----------------------------------------
dynamic_vector.o: dynamic_vector.f90
	($(COMPILE) dynamic_vector.f90)

lin_algebra.o: lin_algebra.f90
	($(COMPILE)  lin_algebra.f90)

csr.o: csr.f90
	($(COMPILE)  csr.f90)

csr_pair_sys.o: csr_pair_sys.f90\
                          dynamic_vector.o\
                          csr.o
	($(COMPILE)  csr_pair_sys.f90)

solve_skit.o: solve_skit.f90\
                        csr.o 
	($(COMPILE)  solve_skit.f90)

np_quadrature.o: np_quadrature.f90\
                           lin_algebra.o
	($(COMPILE) np_quadrature.f90)

select_strategy.o: select_strategy.f90
	($(COMPILE) select_strategy.f90)

nodes.o: nodes.f90\
	           dynamic_vector.o
	($(COMPILE) nodes.f90)

element_topology.o: element_topology.f90\
        	    dynamic_vector.o
	($(COMPILE)   element_topology.f90)

structures.o: structures.f90\
                        dynamic_vector.o
	($(COMPILE) structures.f90)

boundary_geometry.o: boundary_geometry.f90\
                               dynamic_vector.o\
                               structures.o
	($(COMPILE) boundary_geometry.f90)	

mesh_structure.o: mesh_structure.f90\
	                    dynamic_vector.o\
	                    element_topology.o\
	                    nodes.o
	($(COMPILE) mesh_structure.f90)

node_pair_str.o: node_pair_str.f90\
		           nodes.o\
		           mesh_structure.o\
		           dynamic_vector.o
	($(COMPILE)  node_pair_str.f90)

metric_coefficients.o: metric_coefficients.f90\
                       node_pair_str.o
	($(COMPILE) metric_coefficients.f90)

np_topology_gen.o: np_topology_gen.f90\
		             node_pair_str.o\
		             mesh_structure.o\
		             nodes.o\
		             element_topology.o\
		             dynamic_vector.o
	($(COMPILE)  np_topology_gen.f90)

lib_geo.o: lib_geo.f90\
                     dynamic_vector.o
	($(COMPILE)   lib_geo.f90)

derivatives.o: derivatives.f90\
                         metric_coefficients.o\
	                 node_pair_str.o\
	                 nodes.o\
	                 np_quadrature.o\
	                 structures.o
	($(COMPILE) derivatives.f90)

thermodynamic.o: thermodynamic.f90\
                           derivatives.o\
                           metric_coefficients.o\
		           node_pair_str.o\
                           nodes.o\
	                   np_quadrature.o
	($(COMPILE)  thermodynamic.f90)

io_proc.o: io_proc.f90\
                     structures.o\
	             dynamic_vector.o\
	             mesh_structure.o\
                     metric_coefficients.o\
                     nodes.o\
                     np_topology_gen.o\
                     node_pair_str.o\
	             grid_utils.o\
	             lib_geo.o\
		     thermodynamic.o
	($(COMPILE)  io_proc.f90)

solution_interpolation.o: solution_interpolation.f90\
                          structures.o
	($(COMPILE) solution_interpolation.f90)

error_estimator.o: error_estimator.f90\
                             grid_utils.o\
                             structures.o\
		             derivatives.o\
		             thermodynamic.o\
                             np_quadrature.o
	($(COMPILE) error_estimator.f90)

marking_proc.o: marking_proc.f90\
                          structures.o\
                          dynamic_vector.o\
		          select_strategy.o\
		          error_estimator.o\
		          lib_geo.o\
		          grid_utils.o
	($(COMPILE)  marking_proc.f90)

rh_setting.o: rh_setting.f90\
                        structures.o\
                        error_estimator.o\
	                marking_proc.o
	($(COMPILE)   rh_setting.f90)

new_nodes.o: new_nodes.f90\
                       structures.o\
	               grid_utils.o
	($(COMPILE) new_nodes.f90)

new_elements.o: new_elements.f90\
                          structures.o\
		          solution_interpolation.o
	($(COMPILE)  new_elements.f90)

grid_utils.o: grid_utils.f90\
                        dynamic_vector.o\
		        mesh_structure.o\
		        lib_geo.o\
		        metric_coefficients.o\
		        nodes.o\
		        np_topology_gen.o\
		        node_pair_str.o\
		        new_elements.o
	($(COMPILE)  grid_utils.f90)

triangulate_proc.o: triangulate_proc.f90\
                              structures.o\
		              grid_utils.o\
		              marking_proc.o\
		              boundary_geometry.o\
		              new_nodes.o\
		              new_elements.o\
		              nodes.o
	($(COMPILE) triangulate_proc.f90)

reconstruction_proc.o: reconstruction_proc.f90\
		                 structures.o\
		                 dynamic_vector.o\
		                 grid_utils.o\
		                 triangulate_proc.o\
		                 mesh_structure.o\
		                 rh_setting.o\
		                 io_proc.o\
                                 csr.o\
                                 csr_pair_sys.o\
                                 solve_skit.o
	($(COMPILE)  reconstruction_proc.f90)

grid_adaption.o: grid_adaption.f90\
                           io_proc.o\
                           structures.o\
                           rh_setting.o\
		           reconstruction_proc.o
	($(COMPILE)  grid_adaption.f90)

darwin.o: darwin.f90\
                    structures.o\
	            io_proc.o\
	            grid_utils.o\
	            solution_interpolation.o\
	            grid_adaption.o
	($(COMPILE) darwin.f90)
#
#-----------------------------------------
# CLEANING
#-----------------------------------------
.PHONY: clean
clean:
	rm -f ../bins/darwin.* *.o *.mod

