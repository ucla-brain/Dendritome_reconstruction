#!/bin/bash

if [[ "$#" -lt 1 ]]
then
    echo "usage: bash ${BASH_SOURCE[0]} n_procs"
    exit 1
fi
n_procs="${1}"

script_dir="$(dirname "${BASH_SOURCE[0]}")"
script_name="$(basename "${BASH_SOURCE[0]}")"
cd "$script_dir"
project_dir="$(bash "../../util/get_project_dir.sh")"
cd "$project_dir/test/infrastructure"
job_name="mpi_simple_test"
qsub -cwd -c n -N "$job_name" -pe mpich "$n_procs" "$job_name.sh" "$script_name"
