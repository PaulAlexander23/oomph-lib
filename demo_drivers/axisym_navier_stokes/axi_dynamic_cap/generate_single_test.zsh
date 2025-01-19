#!/usr/bin/zsh
# Description: Generate test parameters for the test cases
#
# Usage: ./generate_test_parameters.sh contact_angle bo ca

CONTACT_ANGLE=$1
BO=$2
CA=$3

mkdir -p TEST
# Create a unique filename based on inputs
PARAMETERS_FOLDER=TEST/$CONTACT_ANGLE
mkdir -p $PARAMETERS_FOLDER
PARAMETERS_FOLDER=$PARAMETERS_FOLDER/$BO
mkdir -p $PARAMETERS_FOLDER
PARAMETERS_FOLDER=$PARAMETERS_FOLDER/$CA
mkdir -p $PARAMETERS_FOLDER

cp default_parameters.dat $PARAMETERS_FOLDER/parameters.dat
SPECIAL=$(echo ${PARAMETERS_FOLDER} | perl -pe 's/\//\\\//g')
sed -i "s/90/$CONTACT_ANGLE/g" $PARAMETERS_FOLDER/parameters.dat
sed -i "s/0 # Reynolds Inverse/$BO # Reynolds Inverse/g" $PARAMETERS_FOLDER/parameters.dat
sed -i "s/0 # Wall velocity/$CA # Wall velocity/g" $PARAMETERS_FOLDER/parameters.dat
sed -i "s/RESLT/$SPECIAL/g" $PARAMETERS_FOLDER/parameters.dat

