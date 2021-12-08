#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=\
hele_shaw_uniform_channel_flow_quad \
hele_shaw_uniform_channel_flow_tri \
hele_shaw_taped_channel_flow_quad \
hele_shaw_taped_channel_flow_quad_with_integral \
hele_shaw_perturbed_channel_flow_quad \
filled_channel_flow_quad \
filled_channel_flow_with_integral_quad \
filled_channel_flow_tri \
filled_channel_flow_import_problem_quad \
filled_channel_flow_point_bc_quad \
flow_past_a_cylinder_tri

# Sources for executable
hele_shaw_uniform_channel_flow_quad_SOURCES = hele_shaw_uniform_channel_flow_quad.cc validate.sh
hele_shaw_uniform_channel_flow_tri_SOURCES = hele_shaw_uniform_channel_flow_tri.cc validate.sh
hele_shaw_taped_channel_flow_quad_SOURCES = hele_shaw_taped_channel_flow_quad.cc validate.sh
hele_shaw_taped_channel_flow_quad_with_integral_SOURCES = hele_shaw_taped_channel_flow_quad_with_integral.cc validate.sh
hele_shaw_perturbed_channel_flow_quad_SOURCES = hele_shaw_perturbed_channel_flow_quad.cc validate.sh
filled_channel_flow_quad_SOURCES = filled_channel_flow_quad.cc validate.sh
filled_channel_flow_with_integral_quad_SOURCES = filled_channel_flow_with_integral_quad.cc validate.sh
filled_channel_flow_tri_SOURCES = filled_channel_flow_tri.cc validate.sh
filled_channel_flow_import_problem_quad_SOURCES = filled_channel_flow_import_problem_quad.cc validate.sh
filled_channel_flow_point_bc_quad_SOURCES = filled_channel_flow_point_bc_quad.cc validate.sh
flow_past_a_cylinder_tri_SOURCES = flow_past_a_cylinder_tri.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hele_shaw_uniform_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_uniform_channel_flow_tri_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_taped_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_taped_channel_flow_quad_with_integral_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
hele_shaw_perturbed_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
filled_channel_flow_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
filled_channel_flow_with_integral_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
filled_channel_flow_tri_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
filled_channel_flow_import_problem_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
filled_channel_flow_point_bc_quad_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
flow_past_a_cylinder_tri_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
