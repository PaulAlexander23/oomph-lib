#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
<<<<<<< HEAD
check_PROGRAMS= dense_complex_matrix_test cr_complex_matrix_test cc_complex_matrix_test eigensolver_test
=======
check_PROGRAMS= dense_complex_matrix_test cr_complex_matrix_test cc_complex_matrix_test
>>>>>>> 90e9d245128a6ebfaab03149f249eeafe51bd657

#----------------------------------------------------------------------

# Sources for executable
<<<<<<< HEAD
dense_complex_matrix_test_SOURCES = dense_complex_matrix_test.cc complex_matrices_test_utility_functions.cc
cr_complex_matrix_test_SOURCES = cr_complex_matrix_test.cc complex_matrices_test_utility_functions.cc
cc_complex_matrix_test_SOURCES = cc_complex_matrix_test.cc complex_matrices_test_utility_functions.cc
eigensolver_test_SOURCES = eigensolver_test.cc
=======
dense_complex_matrix_test_SOURCES = dense_complex_matrix_test.cc
cr_complex_matrix_test_SOURCES = cr_complex_matrix_test.cc
cc_complex_matrix_test_SOURCES = cc_complex_matrix_test.cc
>>>>>>> 90e9d245128a6ebfaab03149f249eeafe51bd657

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
dense_complex_matrix_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
cr_complex_matrix_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
cc_complex_matrix_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
<<<<<<< HEAD
eigensolver_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
=======
>>>>>>> 90e9d245128a6ebfaab03149f249eeafe51bd657

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

