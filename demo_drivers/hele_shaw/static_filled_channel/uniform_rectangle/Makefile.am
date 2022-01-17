#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=\
uniform_rectangle_quad \
uniform_rectangle_tri

# Sources for executable
uniform_rectangle_quad_SOURCES = uniform_rectangle_quad.cc validate.sh
uniform_rectangle_tri_SOURCES =  uniform_rectangle_tri.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
uniform_rectangle_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
uniform_rectangle_tri_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
