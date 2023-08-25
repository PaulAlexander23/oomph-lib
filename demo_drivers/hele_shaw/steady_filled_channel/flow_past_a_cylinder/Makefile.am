#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=\
flow_past_a_cylinder_tri

# Sources for executable
flow_past_a_cylinder_tri_SOURCES = flow_past_a_cylinder_tri.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
flow_past_a_cylinder_tri_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
