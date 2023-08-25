#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=\
depth_perturbed_rectangle_tri

# Sources for executable
depth_perturbed_rectangle_tri_SOURCES = depth_perturbed_rectangle_tri.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
depth_perturbed_rectangle_tri_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
