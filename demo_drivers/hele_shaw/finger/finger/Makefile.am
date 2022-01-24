#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=finger

# Sources for executable
finger_SOURCES=finger.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
finger_LDADD = -L@libdir@ -lhele_shaw -lintegral -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
