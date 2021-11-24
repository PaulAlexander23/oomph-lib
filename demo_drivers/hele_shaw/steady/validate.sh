#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=5


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo hele shaw
#----------------------------
cd Validation

echo "Running Hele-Shaw uniform channel flow quad validation "
mkdir RESLT
../hele_shaw_uniform_channel_flow_quad > OUTPUT_hele_shaw_uniform_channel_flow_quad
echo "done"
echo " " >> validation.log
echo "Hele-Shaw uniform channel flow quad validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > hele_shaw_uniform_channel_flow_quad_results.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/hele_shaw_uniform_channel_flow_quad_results.dat.gz   \
    hele_shaw_uniform_channel_flow_quad_results.dat 0.1 1e-13  >> validation.log
fi

echo "Running Hele-Shaw uniform channel flow tri validation "
mkdir RESLT
../hele_shaw_uniform_channel_flow_tri > OUTPUT_hele_shaw_uniform_channel_flow_tri
echo "done"
echo " " >> validation.log
echo "Hele-Shaw uniform channel flow tri validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > hele_shaw_uniform_channel_flow_tri_results.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/hele_shaw_uniform_channel_flow_tri_results.dat.gz   \
    hele_shaw_uniform_channel_flow_tri_results.dat 0.1 1e-13  >> validation.log
fi

echo "Running Hele-Shaw perturbed channel flow quad validation "
mkdir RESLT
../hele_shaw_perturbed_channel_flow_quad > OUTPUT_hele_shaw_perturbed_channel_flow_quad
echo "done"
echo " " >> validation.log
echo "Hele-Shaw perturbed channel flow quad validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > hele_shaw_perturbed_channel_flow_quad_results.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/hele_shaw_perturbed_channel_flow_quad_results.dat.gz   \
    hele_shaw_perturbed_channel_flow_quad_results.dat 0.1 1e-13  >> validation.log
fi

echo "Running Hele-Shaw taped channel flow quad validation "
mkdir RESLT
../hele_shaw_taped_channel_flow_quad > OUTPUT_hele_shaw_taped_channel_flow_quad
echo "done"
echo " " >> validation.log
echo "Hele-Shaw taped channel flow quad validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > hele_shaw_taped_channel_flow_quad_results.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/hele_shaw_taped_channel_flow_quad_results.dat.gz   \
    hele_shaw_taped_channel_flow_quad_results.dat 0.1 1e-13  >> validation.log
fi

echo "Running Hele-Shaw taped channel flow quad with integral validation "
mkdir RESLT
../hele_shaw_taped_channel_flow_quad_with_integral > OUTPUT_hele_shaw_taped_channel_flow_quad_with_integral
echo "done"
echo " " >> validation.log
echo "Hele-Shaw taped channel flow quad with integral validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > hele_shaw_taped_channel_flow_quad_with_integral_results.dat
rm RESLT -r -f

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/hele_shaw_taped_channel_flow_quad_with_integral_results.dat.gz   \
    hele_shaw_taped_channel_flow_quad_with_integral_results.dat 0.1 1e-11  >> validation.log
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
