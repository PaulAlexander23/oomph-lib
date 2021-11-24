#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
time_dependent_channel_height \
uniform_channel_flow_with_time_dependent_flux \
introduction_of_tape \
variable_flux_channel_flow_quad \
variable_channel_depth_quad

# Sources for executable
uniform_channel_flow_with_time_dependent_flux_SOURCES = uniform_channel_flow_with_time_dependent_flux.cc validate.sh
time_dependent_channel_height_SOURCES = time_dependent_channel_height.cc validate.sh
introduction_of_tape_SOURCES = introduction_of_tape.cc validate.sh
variable_flux_channel_flow_quad_SOURCES = variable_flux_channel_flow_quad.cc validate.sh
variable_channel_depth_quad_SOURCES = variable_channel_depth_quad.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
uniform_channel_flow_with_time_dependent_flux_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
time_dependent_channel_height_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
introduction_of_tape_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
variable_flux_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
variable_channel_depth_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
