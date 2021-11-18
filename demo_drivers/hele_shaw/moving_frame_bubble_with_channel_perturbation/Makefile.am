#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=moving_frame_bubble_with_channel_perturbation

# Sources for executable
moving_frame_bubble_with_channel_perturbation_SOURCES = moving_frame_bubble_with_channel_perturbation.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
moving_frame_bubble_with_channel_perturbation_LDADD = -L@libdir@ -lsolid -lconstitutive -lfluid_interface -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
