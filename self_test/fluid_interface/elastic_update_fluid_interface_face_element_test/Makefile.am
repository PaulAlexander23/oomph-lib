#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= qtaylor_hood_fluid_interface_element ttaylor_hood_axi_contact_line

#----------------------------------------------------------------------

# Sources for executable
qtaylor_hood_fluid_interface_element_SOURCES = qtaylor_hood_fluid_interface_element.cc
ttaylor_hood_axi_contact_line_SOURCES = ttaylor_hood_axi_contact_line.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
qtaylor_hood_fluid_interface_element_LDADD = -L@libdir@ -lgeneric -lnavier_stokes -lfluid_interface -lsolid -lconstitutive $(EXTERNAL_LIBS) \
								  $(FLIBS)
ttaylor_hood_axi_contact_line_LDADD = -L@libdir@ -lgeneric -lnavier_stokes -laxisym_navier_stokes -lfluid_interface -lsolid -lconstitutive $(EXTERNAL_LIBS) \
								  $(FLIBS)

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS += -DPARANOID -DRANGE_CHECKING
