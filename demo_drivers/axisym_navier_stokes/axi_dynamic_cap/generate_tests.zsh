#!/usr/bin/zsh
# Description: Generate test parameters for the test cases
#
# Usage: ./generate_tests.zsh 

CONTACT_ANGLE=(45 90 135)
BO=(0 1)
CA=(0 0.1)

for angle in ${CONTACT_ANGLE[@]}; do
    for bo in ${BO[@]}; do
        for ca in ${CA[@]}; do
            echo ./generate_single_test.zsh $angle $bo $ca
            ./generate_single_test.zsh $angle $bo $ca
        done
    done
done
