#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = \
relaxing_bubble

# Sources for executable
relaxing_bubble_SOURCES = relaxing_bubble.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
relaxing_bubble_LDADD = -L@libdir@ -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
