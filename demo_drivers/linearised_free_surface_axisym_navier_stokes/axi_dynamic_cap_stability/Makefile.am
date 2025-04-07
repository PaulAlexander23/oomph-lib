#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

%.h.pch: %.h
	$(CXX) $(AM_CPPFLAGS) -x c++-header -o $@ $<

iwyu.out: ${run_test_SOURCES}
	include-what-you-use -Xiwyu $(AM_CPPFLAGS) run_tests.cc 2> iwyu.out

add-labels: labels.csv add-labels-to-vtu.sh 
	./add-labels-to-vtu.sh labels.csv 'RESLT/perturbed*.vtu'

fix_include: iwyu.out
	fix_include < iwyu.out

# Name of executables
check_PROGRAMS= \
	linear_stability_run \
	neutral_stability_run

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

linear_precompiled_headers = $(patsubst %.cc,%.h.pch,$(linear_sources))

echo-headers:
	echo $(linear_precompiled_headers)

compile-headers: $(linear_precompiled_headers)

linear_stability_run_SOURCES = $(linear_sources) linear_stability_run.cc 
neutral_stability_run_SOURCES = $(linear_sources) neutral_stability_run.cc 

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
