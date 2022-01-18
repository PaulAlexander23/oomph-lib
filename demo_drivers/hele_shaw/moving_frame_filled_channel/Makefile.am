#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
moving_frame_filled_channel

# Sources for executable
moving_frame_filled_channel_SOURCES = moving_frame_filled_channel.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
moving_frame_filled_channel_LDADD = -L@libdir@ -lode -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
