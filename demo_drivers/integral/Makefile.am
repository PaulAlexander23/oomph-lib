#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
	integral \
	integral_tri

# Sources for executable
integral_SOURCES = integral.cc
integral_tri_SOURCES = integral_tri.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
integral_LDADD = -L@libdir@ -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
integral_tri_LDADD = -L@libdir@ -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
