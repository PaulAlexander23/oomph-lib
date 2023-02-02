#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=3

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Validation for complex eigensolver
#-----------------------------------------
echo "Running impose impenetrability elements validation "
../qtaylor_hood_impose_impenetratability_elements > qtaylor_hood_impose_impenetratability_elements.dat
cat test.dat test_face.dat >> qtaylor_hood_impose_impenetratability_elements.dat
echo "done"
echo " " >> validation.log
echo "Impose impenetrability validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/qtaylor_hood_impose_impenetratability_elements.dat.gz  \
         qtaylor_hood_impose_impenetratability_elements.dat >> validation.log
fi

# Validation for complex eigensolver
#-----------------------------------------
echo "Running impose impenetrability elements validation "
../elastic_qtaylor_hood_impose_impenetratability_elements > elastic_qtaylor_hood_impose_impenetratability_elements.dat
cat test.dat test_face.dat >> elastic_qtaylor_hood_impose_impenetratability_elements.dat
echo "done"
echo " " >> validation.log
echo "Impose impenetrability validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/elastic_qtaylor_hood_impose_impenetratability_elements.dat.gz  \
         elastic_qtaylor_hood_impose_impenetratability_elements.dat >> validation.log
fi
#-----------------------------------------
# Validation for complex eigensolver
#-----------------------------------------
echo "Running impose impenetrability elements validation "
../slip > slip.dat
cat test.dat test_face.dat >> slip.dat
echo "done"
echo " " >> validation.log
echo "Impose impenetrability validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/slip.dat.gz  \
         slip.dat >> validation.log
fi

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
