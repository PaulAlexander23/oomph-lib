#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
				qtaylor_hood_impose_impenetratability_elements\
				multi_face_elements_test\
				elastic_qtaylor_hood_impose_impenetratability_elements\
				slip

#----------------------------------------------------------------------

# Sources for executable
qtaylor_hood_impose_impenetratability_elements_SOURCES = qtaylor_hood_impose_impenetratability_elements.cc
multi_face_elements_test_SOURCES = multi_face_elements_test.cc
elastic_qtaylor_hood_impose_impenetratability_elements_SOURCES = elastic_qtaylor_hood_impose_impenetratability_elements.cc
slip_SOURCES = slip.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
qtaylor_hood_impose_impenetratability_elements_LDADD = -L@libdir@ -lgeneric -lnavier_stokes $(EXTERNAL_LIBS) \
								  $(FLIBS)
multi_face_elements_test_LDADD = -L@libdir@ -lgeneric -lnavier_stokes $(EXTERNAL_LIBS) \
								  $(FLIBS)
elastic_qtaylor_hood_impose_impenetratability_elements_LDADD = -L@libdir@ -lgeneric -lnavier_stokes -lsolid -lconstitutive $(EXTERNAL_LIBS) \
								  $(FLIBS)
slip_LDADD = -L@libdir@ -lgeneric -lnavier_stokes -lsolid -lconstitutive $(EXTERNAL_LIBS) \
								  $(FLIBS)

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS += -DPARANOID -DRANGE_CHECKING
