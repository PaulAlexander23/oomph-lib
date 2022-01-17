#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
time_dependent_channel_height

# Sources for executable
time_dependent_channel_height_SOURCES = time_dependent_channel_height.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
time_dependent_channel_height_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
