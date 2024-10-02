#!/usr/bin/env bash

if [[ "$#" -lt 1 ]]
then
    echo "usage: bash qsub_mpi_io_test.sh n_writers"
    exit 1
fi
n_writers="${1}"
script_dir="$(dirname "${BASH_SOURCE[0]}")"
script_name="$(basename "${BASH_SOURCE[0]}")"
cd "$script_dir"
project_dir="$(bash "../../util/get_project_dir.sh")"
cd "$project_dir/test/infrastructure"
export QSUB_BINDING_STRIDING=1
job_name="mpi_io_test"
qsub -cwd -c n -N "$job_name" -pe mpich "$n_writers" "$job_name.sh" "$script_name"
#qsub -cwd -c n -binding striding:"$1":"$QSUB_BINDING_STRIDING" -N "$job_name" -pe mpich "$n_writers" "$job_name.sh" "$script_name"
