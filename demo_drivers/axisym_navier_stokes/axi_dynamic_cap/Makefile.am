#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Variables
EXEC:= \
	axi_dynamic_cap \
	run_eigenproblem \
	create_parameter_files \
	continuation_run \
	run_tests

# Makefile rules
run: ${EXEC}

%.iwyu: %.cc
	include-what-you-use $(AM_CPPFLAGS) $< 2> $@

# Run parallel axi_dynamic_cap with the parameters in the validata directory
# test: axi_dynamic_cap
# 	./generate_tests.zsh
# 	./run_tests.zsh
# 	cat tests.log
iwyu: $(IWYU)

%.fixed: %.iwyu
	fix_include < $< > $@

fix_include: $(IWYU_FIXED)

test.log: run_tests
	./run_tests -t extra > test.log

update-regression-tests-validata: 
	gzip Validation/*.dat
	cp Validation/*.dat.gz validata/


# Automake
#continuation_run_SOURCES = continuation_run.cc

# Name of executables
check_PROGRAMS= ${EXEC}
# Name of library used for linking to within the linearised free surface
# axisymmetric Navier-Stokes demos which performs a linear stability about the
# solutions in this folder.
noinst_LIBRARIES = libaxi_dynamic_cap.a

sources = parameters.cc parameter_functions.cc base_element.cc \
		  singular_axisym_dynamic_cap_problem.cc

# Pattern replace sources cc with iwyu
IWYU=$(patsubst %.cc,%.iwyu,$(sources))
IWYU_FIXED=$(patsubst %.iwyu,%.fixed,$(IWYU))

axi_dynamic_cap_SOURCES = $(sources) axi_dynamic_cap.cc 
run_eigenproblem_SOURCES = $(sources) run_eigenproblem.cc 
create_parameter_files_SOURCES = $(sources) create_parameter_files.cc
continuation_run_SOURCES = $(sources) continuation_run.cc
run_tests_SOURCES = $(sources) run_tests.cc
libaxi_dynamic_cap_a_SOURCES = $(sources)

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
LDADD = -L@libdir@ \
	-laxisym_navier_stokes \
	-lnavier_stokes \
	-lfluid_interface \
	-lsolid \
	-lconstitutive \
	-lmeshes \
	-lgeneric \
	$(EXTERNAL_LIBS) $(FLIBS)
#-lboost_unit_test_framework \
#-lgcov

# Add the paranoid and range checking flags
AM_CPPFLAGS += -Wno-misleading-indentation -Wfatal-errors #-I../../../build/include/ #-DPARANOID -DRANGE_CHECKING --coverage -O0 # 
