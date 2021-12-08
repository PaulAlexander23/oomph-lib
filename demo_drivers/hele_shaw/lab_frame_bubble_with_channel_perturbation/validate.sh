#!/usr/bin/sh

# Get the OOPMH-LIB root directory from a makefile
# OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
EXPECTED_NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Validation for lab frame bubble with channel perturbation
#-----------------------------------------
mkdir RESLT
echo "Running lab frame validation "
../lab_frame_bubble_with_channel_perturbation -q 0.01 -n 100 > \
    OUTPUT_lab_frame_bubble_with_channel_perturbation
echo "done"
echo " " >> validation.log
echo "Lab frame validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/* > lab_frame_bubble_with_channel_perturbation.dat

diff=$(zcmp ../validata/lab_frame_bubble_with_channel_perturbation.dat.gz lab_frame_bubble_with_channel_perturbation.dat )
if [ $? != 0 ]; then
    echo "[ERROR] Compare failed to run." >> validation.log
elif [ $diff ]; then
    echo "[FAILED]" >> validation.log
else
    echo "[OK]" >> validation.log
fi
#-----------------------------------------

cd ..

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
ACTUAL_NUM_TESTS_PASSED=$(grep OK Validation/validation.log | wc -l)
ACTUAL_NUM_TESTS_FAILED=$(grep FAILED Validation/validation.log | wc -l)
ACTUAL_NUM_TESTS_ERROR=$(grep ERROR Validation/validation.log | wc -l)
if [ $ACTUAL_NUM_TESTS_ERROR -gt 0 ]; then
    echo "Error in $ACTUAL_NUM_TESTS_ERROR testing scripts. Check validation.log."
    return 1
fi
if [ $ACTUAL_NUM_TESTS_FAILED -gt 0 ]; then
    echo "Failed $ACTUAL_NUM_TESTS_FAILED tests. Check validation.log."
    return 1
fi
if [ $ACTUAL_NUM_TESTS_PASSED -gt $EXPECTED_NUM_TESTS ]; then
    echo "Passed more tests than expected! Check EXPECTED_NUM_TESTS and validation.log."
    return 2
fi
if [ $ACTUAL_NUM_TESTS_PASSED -eq $EXPECTED_NUM_TESTS ]; then
    echo "Passed all $ACTUAL_NUM_TESTS_PASSED tests."
    return 0
fi
if [ $ACTUAL_NUM_TESTS_PASSED -lt $EXPECTED_NUM_TESTS ]; then
    echo "Passed less tests than expected! Check EXPECTED_NUM_TESTS and validation.log."
    return 1
fi

# Never get here
exit 10
