#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=bubble_in_a_uniform_channel

# Sources for executable
bubble_in_a_uniform_channel_SOURCES = bubble_in_a_uniform_channel.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
bubble_in_a_uniform_channel_LDADD = -L@libdir@ -lsolid -lconstitutive -lfluid_interface -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
