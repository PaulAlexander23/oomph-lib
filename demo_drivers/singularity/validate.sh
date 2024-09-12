#!/bin/sh

PASSED=0
FAILED=1
MORETEST=2

./test_case_1

ACTUAL=$?

if [ $ACTUAL -ne 0 ]; then
    exit $FAILED
fi

exit $PASSED
