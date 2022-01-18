#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
introduction_of_tape

# Sources for executable
introduction_of_tape_SOURCES = introduction_of_tape.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
introduction_of_tape_LDADD = -L@libdir@ -lhele_shaw -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
