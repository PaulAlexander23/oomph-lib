#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
uniform_channel

# Sources for executable
uniform_channel_SOURCES = uniform_channel.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
uniform_channel_LDADD = -L@libdir@ -lode -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
