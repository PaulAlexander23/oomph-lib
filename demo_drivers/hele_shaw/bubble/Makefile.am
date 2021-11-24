#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=bubble_unsteady

# Sources for executable
bubble_unsteady_SOURCES = bubble_unsteady.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
bubble_unsteady_LDADD = -L@libdir@ -lsolid -lconstitutive -lfluid_interface -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
