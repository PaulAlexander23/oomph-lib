#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = \
relaxing_bubble \
relaxing_bubble_dirichlet_side \
relaxing_bubble_fixed_volume

# Sources for executable
relaxing_bubble_SOURCES = relaxing_bubble.cc
relaxing_bubble_dirichlet_side_SOURCES = relaxing_bubble_dirichlet_side.cc
relaxing_bubble_fixed_volume_SOURCES = relaxing_bubble_fixed_volume.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
relaxing_bubble_LDADD = -L@libdir@ -lhele_shaw -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
relaxing_bubble_dirichlet_side_LDADD = -L@libdir@ -lhele_shaw -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
relaxing_bubble_fixed_volume_LDADD = -L@libdir@ -lhele_shaw -lfluid_interface -lconstitutive -lsolid -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
