#!/bin/bash

#$ -V
#$ -o $HOME/log/$JOB_NAME.o$JOB_ID
#$ -e $HOME/log/$JOB_NAME.e$JOB_ID

if [[ "$#" < 2 ]]
then
    echo "usage: remove_path_mpi.sh caller root_dir"
    exit 1
fi

project_dir="$(bash "../../util/get_project_dir.sh")"
build_dir="$project_dir/grid_build_master"
util_dir="$project_dir/util"
verify_caller_script="$util_dir/verify_caller.sh"
source "$verify_caller_script"

async_remove_path_mpi_bin="$build_dir/parallel/async_remove_path_mpi"

if [ ! -e "$async_remove_path_mpi_bin" ]
then
    echo "executable $async_remove_path_mpi_bin not found"
    exit 1
fi

root_dir="$2"
if [ ! -e "$root_dir" ]; then
    echo "$root_dir does not exist"
    exit 1
elif [ -f "$root_dir" ]; then
    rm "$root_dir"
    echo "removed file $root_dir"
    exit 0
else
    /usr/local/mpich-3.1.3/bin/mpiexec -n "$NSLOTS" "$async_remove_path_mpi_bin" "$root_dir"
fi

