#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
ode_controlled_channel_height \
hele_shaw_perturbation_moving_reference_frame \
hele_shaw_perturbation_moving_reference_frame_flux

# Sources for executable
ode_controlled_channel_height_SOURCES = ode_controlled_channel_height.cc validate.sh
hele_shaw_perturbation_moving_reference_frame_SOURCES = hele_shaw_perturbation_moving_reference_frame.cc validate.sh
hele_shaw_perturbation_moving_reference_frame_flux_SOURCES = hele_shaw_perturbation_moving_reference_frame_flux.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
ode_controlled_channel_height_LDADD = -L@libdir@ -lode -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_perturbation_moving_reference_frame_LDADD = -L@libdir@ -lode -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_perturbation_moving_reference_frame_flux_LDADD = -L@libdir@ -lode -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
