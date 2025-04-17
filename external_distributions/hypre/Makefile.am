EXTRA_DIST = \
build_scripts

#makefile variable to hold the configure command (part 1)
run_command1 = " ./configure "
if WANT_MPI
run_command1 += " --with-MPI  \\" 
else
run_command1 += " --without-MPI  \\"  
endif

#makefile variable to hold the configure command (part 2)
run_command2 = " "
if ALREADY_HAVE_BLAS
run_command2 += " --with-blas-libs=\"$(BLAS_LIB)\" --with-blas-libs-dir=\"/home/paul/apps/open-blas/lib/\" "
else
run_command2 += "  --with-blas=dummy  "
endif
if ALREADY_HAVE_LAPACK
run_command2 += " --with-lapack-libs=\"$(LAPACK_LIB)\" --with-lapack-libs-dir=\"/home/paul/apps/open-blas/lib/\" "
else
run_command2 += " --with-lapack=dummy" 
endif
run_command2 += " F77=\"$(F77)\"  CC=\"$(CC)\" CXX=\"$(CXX)\" FFLAGS=\"$(FFLAGS)\"  CFLAGS=\"$(CFLAGS)\" CXXFLAGS=\"$(CXXFLAGS)\" LDFLAGS=\"$(LDFLAGS)\" FLIBS=\"-L/opt/AMD/aocc-compiler-5.0.0/lib -lflang -lflangrti -lompstub -lpgmath -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi\" " 

all:
	echo "#! /bin/bash" > build_script.bash
	echo " " >> build_script.bash
	echo "# WARNING: THIS IS AN AUTOMATICALLY GENERATED FILE" >> build_script.bash
	echo "#          DO NOT EDIT!" >> build_script.bash
	echo " " >> build_script.bash
if ALREADY_HAVE_BUILT_HYPRE                                                                 
	echo "ALREADY_HAVE_BUILT_HYPRE=1" >> build_script.bash
else
	echo "ALREADY_HAVE_BUILT_HYPRE=0" >> build_script.bash
endif	
	cat $(srcdir)/build_scripts/build_script.part1 >> build_script.bash
	echo " " >> build_script.bash
	echo $(run_command1) >> build_script.bash
	cat $(srcdir)/build_scripts/build_script.part2 >> build_script.bash
	echo $(run_command2) >> build_script.bash
	echo " " >> build_script.bash 
	echo "# Build and install" >> build_script.bash 
	echo "make">> build_script.bash 
	echo "make install" >> build_script.bash 
	echo " ">> build_script.bash 
	echo "#Cleanup">> build_script.bash 
	echo "cd ../..">> build_script.bash 
	echo "# rm -rf hypre-2.0.0">> build_script.bash 

	chmod a+x build_script.bash

	./build_script.bash 

clean-local:
	rm -f build_script.bash
	rm -rf hypre-2.0.0

distclean-local:
	rm -f build_script.bash
	rm -rf hypre-2.0.0
	rm -rf hypre_default_installation

