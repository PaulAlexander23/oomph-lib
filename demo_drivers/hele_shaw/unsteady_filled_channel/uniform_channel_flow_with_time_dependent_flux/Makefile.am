#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
uniform_channel_flow_with_time_dependent_flux

# Sources for executable
uniform_channel_flow_with_time_dependent_flux_SOURCES = uniform_channel_flow_with_time_dependent_flux.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
uniform_channel_flow_with_time_dependent_flux_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
