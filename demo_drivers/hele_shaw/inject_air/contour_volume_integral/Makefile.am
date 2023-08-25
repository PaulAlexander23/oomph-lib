#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = \
inject_air

# Sources for executable
inject_air_SOURCES = inject_air.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
inject_air_LDADD = -L@libdir@ -lhele_shaw -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
