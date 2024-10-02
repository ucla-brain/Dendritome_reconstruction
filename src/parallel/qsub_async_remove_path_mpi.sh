#!/bin/bash

if [[ "$#" < 2 ]]
then
    echo "usage: qsub_remove_path_mpi.sh n_procs root_dir"
    exit 1
fi

n_procs="$1"
root_dir="$2"

script_dir="$(dirname "${BASH_SOURCE[0]}")"
script_name="$(basename "${BASH_SOURCE[0]}")"
cd "$script_dir"
project_dir="$(bash "../../util/get_project_dir.sh")"
cd "$project_dir/grid_build_master/parallel"

async_remove_path_mpi="async_remove_path_mpi.sh"

if [ ! -e "$async_remove_path_mpi" ]
then
    echo "executable $async_remove_path_mpi not found"
    exit 1
fi

job_name="async_remove_path_mpi"
qsub -cwd -c n -N "$job_name" -pe mpich "$n_procs" "$async_remove_path_mpi" "$script_name" "$root_dir"