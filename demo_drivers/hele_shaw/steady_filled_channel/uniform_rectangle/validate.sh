#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo hele shaw
#----------------------------
cd Validation

echo "Running Hele-Shaw steady uniform rectangle quad validation "
mkdir RESLT
../uniform_rectangle_quad > OUTPUT_uniform_rectangle_quad
echo "done"
echo " " >> validation.log
echo "Hele-Shaw steady uniform rectangle quad validation " >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > uniform_rectangle_quad.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/uniform_rectangle_quad.dat.gz   \
    uniform_rectangle_quad.dat 0.1 1e-11  >> validation.log
fi

echo "Running Hele-Shaw steady uniform rectangle tri validation "
mkdir RESLT
../uniform_rectangle_tri > OUTPUT_uniform_rectangle_tri
echo "done"
echo " " >> validation.log
echo "Hele-Shaw steady uniform rectangle tri validation " >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > uniform_rectangle_tri.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../../bin/fpdiff.py ../validata/uniform_rectangle_tri.dat.gz   \
    uniform_rectangle_tri.dat 0.1 1e-11  >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
