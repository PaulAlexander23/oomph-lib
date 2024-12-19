#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Variables
EXEC:= \
	structured_no_correction \
	structured_no_correction_region \
	structured_with_correction \
	structured_with_correction_region \
	unstructured_no_correction \
	unstructured_with_correction \
	hijack_boundary

FOLDERS:= \
	RESLT_axi_fix \
	RESLT_axi_fix_unstr \
	RESLT_axi_fix_region \
	RESLT_axi_no_fix \
	RESLT_axi_no_fix_unstr \
	RESLT_axi_no_fix_region \
	RESLT_hijack

# Makefile rules
run: ${EXEC} ${FOLDERS} parameters.dat
	echo ${EXEC} | xargs -d " " -I {} bash -c "./{}"

$(FOLDERS):
	mkdir -p $@

parameters.dat: 
	cp default_parameters.dat parameters.dat

# Automake

# Name of executables
check_PROGRAMS= ${EXEC}

# Sources for second executable
structured_no_correction_SOURCES = structured_no_correction.cc
structured_no_correction_region_SOURCES = structured_no_correction_region.cc
structured_with_correction_SOURCES = structured_with_correction.cc
structured_with_correction_region_SOURCES = structured_with_correction_region.cc
unstructured_no_correction_SOURCES = unstructured_no_correction.cc
unstructured_with_correction_SOURCES = unstructured_with_correction.cc
hijack_boundary_SOURCES = hijack_boundary.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
LDADD = -L@libdir@ \
		-lnavier_stokes \
		-laxisym_navier_stokes \
		-lfluid_interface \
		-lsolid \
		-lconstitutive \
		-lgeneric \
		$(EXTERNAL_LIBS) $(FLIBS) \
		-lgcov 

# Add the paranoid and range checking flags
AM_CPPFLAGS +=-Wno-misleading-indentation -ggdb -O0 --coverage #-DPARANOID  -DRANGE_CHECKING -ggdb -O0 -Wno-misleading-indentation 
