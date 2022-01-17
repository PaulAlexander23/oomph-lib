#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
variable_channel_depth_quad

# Sources for executable
variable_channel_depth_quad_SOURCES = variable_channel_depth_quad.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
variable_channel_depth_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
