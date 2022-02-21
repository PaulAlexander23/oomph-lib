#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
variable_flux_channel_flow_quad

# Sources for executable
variable_flux_channel_flow_quad_SOURCES = variable_flux_channel_flow_quad.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
variable_flux_channel_flow_quad_LDADD = -L@libdir@ -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
