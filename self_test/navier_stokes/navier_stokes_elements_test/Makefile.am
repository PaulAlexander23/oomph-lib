#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= qtaylor_hood_elements

#----------------------------------------------------------------------

# Sources for executable
qtaylor_hood_elements_SOURCES = qtaylor_hood_elements.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
qtaylor_hood_elements_LDADD = -L@libdir@ -lgeneric -lnavier_stokes $(EXTERNAL_LIBS) \
								  $(FLIBS)

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS += -DPARANOID -DRANGE_CHECKING
