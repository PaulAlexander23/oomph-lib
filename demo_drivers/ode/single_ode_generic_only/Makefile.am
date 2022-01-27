#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=single_ode_generic_only

# Sources for executable
single_ode_generic_only_SOURCES = single_ode.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
single_ode_generic_only_LDADD = -L@libdir@ -lode -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST+=zeros
