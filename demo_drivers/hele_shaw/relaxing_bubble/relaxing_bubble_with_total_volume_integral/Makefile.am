#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = \
relaxing_bubble_fixed_volume \
relaxing_bubble_traded

# Sources for executable
relaxing_bubble_fixed_volume_SOURCES = relaxing_bubble_fixed_volume.cc
relaxing_bubble_traded_SOURCES = relaxing_bubble_traded.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
relaxing_bubble_fixed_volume_LDADD = -L@libdir@ -lhele_shaw -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
relaxing_bubble_traded_LDADD = -L@libdir@ -lhele_shaw -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
