#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Validation for complex eigensolver
#-----------------------------------------
echo "Running complex eigensolver validation "
../axisym_ttaylor_hood_pvd_elements > axisym_ttaylor_hood_pvd_elements.dat
echo "done"
echo " " >> validation.log
echo "Complex eigensolver validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat test.dat >> axisym_ttaylor_hood_pvd_elements.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/axisym_ttaylor_hood_pvd_elements.dat.gz  \
         axisym_ttaylor_hood_pvd_elements.dat >> validation.log 0.1 1e-9
fi
rm -rf RESLT

#-----------------------------------------

# Append log to main validation log
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