#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
	test_case_1 \
	run_tests

# Sources for second executable
test_case_1_SOURCES = test_case_1.cc
run_tests_SOURCES = run_tests.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
test_case_1_LDADD = -L@libdir@ -lnavier_stokes -laxisym_navier_stokes -lfluid_interface \
									-lsolid -lconstitutive \
									-lgeneric -lgcov \
									$(EXTERNAL_LIBS) $(FLIBS)
run_tests_LDADD = -L@libdir@ -lnavier_stokes -laxisym_navier_stokes -lfluid_interface \
									-lsolid -lconstitutive \
									-lgeneric -lgcov -lboost_unit_test_framework \
									$(EXTERNAL_LIBS) $(FLIBS)

# Add the paranoid and range checking flags
#AM_CPPFLAGS += -DPARANOID  -DRANGE_CHECKING -ggdb -Og -Wno-misleading-indentation --coverage
