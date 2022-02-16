#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
	integral

# Sources for executable
integral_SOURCES = integral.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
integral_LDADD = -L@libdir@ -lintegral -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
