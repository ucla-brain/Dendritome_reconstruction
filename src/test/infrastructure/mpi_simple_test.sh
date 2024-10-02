#!/bin/bash
#$ -V
#$ -o $HOME/log/$JOB_NAME.o$JOB_ID
#$ -e $HOME/log/$JOB_NAME.e$JOB_ID

project_dir="$(bash "../../util/get_project_dir.sh")"
build_dir="$project_dir/grid_build"
util_dir="$project_dir/util"
verify_caller_script="$util_dir/verify_caller.sh"
if [ ! -e "$verify_caller_script" ]
then
    echo "$verify_caller_script does not exist"
    exit 1
fi
source "$verify_caller_script"

mpi_simple_test_bin="$build_dir/test/infrastructure/mpi_simple_test"
if [ ! -e "$mpi_simple_test_bin" ]
then
    echo "executable $mpi_simple_test_bin not found"
    exit 1
fi

/usr/local/mpich-3.1.3/bin/mpiexec  -n $NSLOTS "$mpi_simple_test_bin"
exit 0
