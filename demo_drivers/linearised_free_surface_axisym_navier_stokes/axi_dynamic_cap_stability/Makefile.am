#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

CC=ccache gcc
CXX=ccache g++


#Temporary!
all: run_tests
	./run_tests --run_test=free_surface_fixed_c

perturbed.pvd: run_tests
	./run_tests --run_test=mode_0_fix_c
	../../../bin/oomph-convert.py RESLT/perturbed0.dat
	../../../bin/oomph-convert.py RESLT/perturbed1.dat
	$(MAKE) add-labels
	../../../bin/makePvd RESLT/perturbed perturbed.pvd
	paraview perturbed.pvd

iwyu.out: ${run_test_SOURCES}
	include-what-you-use -Xiwyu $(AM_CPPFLAGS) run_tests.cc 2> iwyu.out

add-labels: labels.csv add-labels-to-vtu.sh 
	./add-labels-to-vtu.sh labels.csv 'RESLT/perturbed*.vtu'
	

fix_include: iwyu.out
	fix_include < iwyu.out

# Name of executables
check_PROGRAMS= \
	linear_stability_run \
	neutral_stability_run \
	run_tests

# Sources for second executable
linear_sources = linearised_axisymmetric_fluid_interface_elements.cc \
				 full_linearised_axisym_navier_stokes_elements.cc \
				 Tfull_linearised_axisym_navier_stokes_elements.cc \
				 decomposed_pvd_elements.cc \
				 decomposed_linear_elasticity_elements.cc \
				 base_element.cc \
				 base_problem.cc \
				 overlaying_my_linear_element.cc \
				 singular_overlaying_my_linear_elements.cc \
				 perturbed_linear_stability_cap_problem_base.cc \
				 singular_perturbed_linear_stability_cap_problem.cc


linear_stability_run_SOURCES = $(linear_sources) linear_stability_run.cc 
neutral_stability_run_SOURCES = $(linear_sources) neutral_stability_run.cc 
run_tests_SOURCES = $(linear_sources) run_tests.cc 

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
LDADD = -L@libdir@ \
		../../axisym_navier_stokes/axi_dynamic_cap/parameters.o \
		../../axisym_navier_stokes/axi_dynamic_cap/parameter_functions.o \
	-llinear_elasticity \
	-llinearised_axisym_navier_stokes \
	-laxisym_navier_stokes \
	-lfluid_interface \
	-lsolid \
	-lconstitutive \
	-lnavier_stokes \
	-lgeneric \
	$(EXTERNAL_LIBS) $(FLIBS)# -lboost_unit_test_framework -lgcov

# Add the paranoid and range checking flags
AM_CPPFLAGS += -g0 -O0 -w -Wno-misleading-indentation #-DPARANOID  -DRANGE_CHECKING #-Wfatal-errors # --coverage -O0
