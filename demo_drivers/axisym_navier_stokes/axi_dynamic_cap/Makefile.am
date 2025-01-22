#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Variables
EXEC:= \
	axi_dynamic_cap \
	create_parameter_files \
	continuation_run \
	run_tests

# Makefile rules
run: ${EXEC}

# Run parallel axi_dynamic_cap with the parameters in the validata directory
test: axi_dynamic_cap
	./generate_tests.zsh
	./run_tests.zsh
	cat tests.log

# Automake
#continuation_run_SOURCES = continuation_run.cc

# Name of executables
check_PROGRAMS= ${EXEC}

axi_dynamic_cap_SOURCES = axi_dynamic_cap.cc 
create_parameter_files_SOURCES = create_parameter_files.cc
continuation_run_SOURCES = continuation_run.cc
run_tests_SOURCES = run_tests.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
LDADD = -L@libdir@ \
	-laxisym_navier_stokes \
	-lfluid_interface \
	-lsolid \
	-lconstitutive \
	-lgeneric \
	$(EXTERNAL_LIBS) $(FLIBS)
#-lboost_unit_test_framework \
#-lgcov

# Add the paranoid and range checking flags
AM_CPPFLAGS += -Wno-misleading-indentation -Wfatal-errors #-I../../../build/include/ #-DPARANOID -DRANGE_CHECKING --coverage -O0 # 
