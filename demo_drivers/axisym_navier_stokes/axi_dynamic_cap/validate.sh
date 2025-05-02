#! /bin/sh

# Default: serial
MAX_JOBS=1

# Get the number of cpu cores
NUM_CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -j) 
            shift
            if [[ "$1" =~ ^[0-9]+$ ]]; then
                MAX_JOBS="$1"
                shift
            else
                MAX_JOBS="$NUM_CORES"
            fi
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=23

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for linear stability run with restart
#--------------------------------------------------
validate(){
    # Expects to be called as
    # validate executable input actual_output_filename expected_output_filename

  TEMP_DIR=$(mktemp)
  cd $TEMP_DIR

  mkdir -p Validation/RESLT

  var="../$1 $2"
  echo $var
  eval $var 
  if test "$?" = "0"; then
      echo "done"
      LOG="Validation/"$(mktemp)
      echo " " >> $LOG 
      echo "Validation run" >> $LOG
      echo "---------------------------------------------" >> $LOG
      echo " " >> $LOG
      echo "Validation directory: " >> $LOG
      echo " " >> $LOG
      echo "  " `pwd` >> $LOG
      echo " " >> $LOG

      ## Compute expected file name 
      #TEMP=$(basename -- $2)
      #FILE=${1%% *}"_"${TEMP%.*}"_results.dat"
      FILE=${4%.*}
      echo Validation/$FILE
      cat  Validation/RESLT/$3 > Validation/$FILE

      if test "$1" = "no_fpdiff"; then
        echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> $LOG
      else
        ../../../../bin/fpdiff.py validata/$4  \
               Validation/$FILE 0.1 2e-7 >> $LOG
      fi
  else
      echo "done - returned error"
      LOG="Validation/validation.log"
      echo " " >> $LOG 
      echo "Validation run" >> $LOG
      echo "---------------------------------------------" >> $LOG
      echo " " >> $LOG
      echo "Validation directory: " >> $LOG
      echo " " >> $LOG
      echo "  " `pwd` >> $LOG
      echo " " >> $LOG
      FILE=${4%.*}
      echo Validation/$FILE
      echo "[FAILED] -- Returned and exit code" >> $LOG
  fi


  # Copy temp log into the main log
  cat $LOG > "Validation/validation.log"

  # Exit temporary directory and quit
  cd ..
  rm $TEMP_DIR -r
}

run_tests(){
    touch RESLT
    rm -r -f RESLT
    mkdir RESLT

    var="./run_tests -t normal > Validation/OUTPUT"
    echo $var
    eval $var
    echo "done"
    LOG="Validation/validation.log"
    echo " " >> $LOG 
    echo "Validation run" >> $LOG
    echo "---------------------------------------------" >> $LOG
    echo " " >> $LOG
    echo "Validation directory: " >> $LOG
    echo " " >> $LOG
    echo "  " `pwd` >> $LOG
    echo " " >> $LOG
    if grep "fatal error" Validation/OUTPUT; then
        cat Validation/OUTPUT >> $LOG
        echo "[FAILED] -- Unit tests failed see validation log." >> $LOG
    else
        echo "[OK] -- Unit tests passed" >> $LOG
    fi

}

# Throttle to N parrallel jobs
run_with_limit() {
    while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
        sleep 0.2
    done
    "$@" &
}



# Utility scripts
run_with_limit validate "create_parameter_files --folder Validation/RESLT --overwrite --parameters" validata/unsteady-parameters-with-restart.dat parameters.dat create_parameter_files_unsteady-parameters-with-restart_results.dat.gz

# Base state scripts
run_with_limit validate "axi_dynamic_cap --parameters " validata/parameters.dat trace.dat axi_dynamic_cap_parameters_results.dat.gz
run_with_limit validate "axi_dynamic_cap --parameters " validata/parameters-with-restart.dat trace.dat axi_dynamic_cap_parameters-with-restart_results.dat.gz
run_with_limit validate "axi_dynamic_cap --parameters " validata/unsteady-parameters.dat trace.dat axi_dynamic_cap_unsteady-parameters_results.dat.gz
run_with_limit validate "axi_dynamic_cap --parameters " validata/unsteady-parameters-with-restart.dat trace.dat axi_dynamic_cap_unsteady-parameters-with-restart_results.dat.gz

# Obtuse runs
run_with_limit validate "axi_dynamic_cap --parameters " validata/obtuse-parameters.dat trace.dat axi_dynamic_cap_obtuse-parameters_results.dat.gz
run_with_limit validate "axi_dynamic_cap --parameters " validata/obtuse-parameters-with-restart.dat trace.dat axi_dynamic_cap_obtuse-parameters-with-restart_results.dat.gz
run_with_limit validate "axi_dynamic_cap --parameters " validata/obtuse-unsteady-parameters.dat trace.dat unsteady_run_parameters_results.dat.gz
run_with_limit validate "axi_dynamic_cap --parameters " validata/obtuse-unsteady-parameters-with-restart.dat trace.dat unsteady_run_parameters-with-restart_results.dat.gz

# Continuation runs
run_with_limit validate "continuation_run --Bo 0.1 --parameters " validata/unsteady-parameters-with-restart.dat trace.dat cont-bo-results.dat.gz
run_with_limit validate "continuation_run --wall_velocity 0.01 --parameters" validata/unsteady-parameters-with-restart.dat trace.dat cont-ca-results.dat.gz
run_with_limit validate "continuation_run --arc --Bo 0.01 --parameters" validata/unsteady-parameters-with-restart.dat trace.dat arc-cont-bo-results.dat.gz
run_with_limit validate "continuation_run --arc --wall_velocity 0.01 --parameters" validata/unsteady-parameters-with-restart.dat trace.dat continuation_run_unsteady-parameters-with-restart_results.dat.gz
run_with_limit validate "continuation_run --height_control --Bo 0.01 --parameters" validata/unsteady-parameters-with-restart.dat trace.dat height-cont-bo-results.dat.gz
run_with_limit validate "continuation_run --height_control --wall_velocity 0.1 --parameters" validata/unsteady-parameters-with-restart.dat trace.dat continuation_run_height-continuation-parameters-with-restart_results.dat.gz

# Test exception handling.
run_with_limit validate "continuation_run --wall_velocity 0.1 --parameters" validata/unsteady-parameters-with-restart-reduced-tolerance.dat exit_code.txt exit_code.txt.gz
 
# Obtuse continuation runs
run_with_limit validate "continuation_run --Bo 0.1 --parameters " validata/obtuse-unsteady-parameters-with-restart.dat trace.dat obtuse-cont-bo-results.dat.gz
run_with_limit validate "continuation_run --wall_velocity 0.01 --parameters" validata/obtuse-unsteady-parameters-with-restart.dat trace.dat obtuse-cont-ca-results.dat.gz
run_with_limit validate "continuation_run --arc --Bo 0.1 --parameters" validata/obtuse-unsteady-parameters-with-restart.dat trace.dat obtuse-arc-cont-bo-results.dat.gz
run_with_limit validate "continuation_run --arc --wall_velocity 0.01 --parameters" validata/obtuse-unsteady-parameters-with-restart.dat trace.dat continuation_run_obtuse-unsteady-parameters-with-restart_results.dat.gz
run_with_limit validate "continuation_run --height_control --Bo 0.01 --parameters" validata/obtuse-unsteady-parameters-with-restart.dat trace.dat obtuse-height-cont-bo-results.dat.gz
run_with_limit validate "continuation_run --height_control --wall_velocity 0.01 --parameters" validata/obtuse-unsteady-parameters-with-restart.dat trace.dat continuation_run_obtuse-height-continuation-parameters-with-restart_results.dat.gz

run_with_limit run_tests

wait
#######################################################################

# Append log to main validation log
LOG="Validation/validation.log"
cat $LOG >> ../../../validation.log

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
