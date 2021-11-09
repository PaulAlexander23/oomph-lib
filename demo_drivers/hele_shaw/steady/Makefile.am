#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=hele_shaw_uniform_channel_flow_quad hele_shaw_uniform_channel_flow_tri hele_shaw_taped_channel_flow_quad hele_shaw_taped_channel_flow_quad_with_integral hele_shaw_perturbed_channel_flow_quad 

# Sources for executable
hele_shaw_uniform_channel_flow_quad_SOURCES = hele_shaw_uniform_channel_flow_quad.cc validate.sh
hele_shaw_uniform_channel_flow_tri_SOURCES = hele_shaw_uniform_channel_flow_tri.cc validate.sh
hele_shaw_taped_channel_flow_quad_SOURCES = hele_shaw_taped_channel_flow_quad.cc validate.sh
hele_shaw_taped_channel_flow_quad_with_integral_SOURCES = hele_shaw_taped_channel_flow_quad_with_integral.cc validate.sh
hele_shaw_perturbed_channel_flow_quad_SOURCES = hele_shaw_perturbed_channel_flow_quad.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hele_shaw_uniform_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_uniform_channel_flow_tri_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_taped_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_taped_channel_flow_quad_with_integral_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_perturbed_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
