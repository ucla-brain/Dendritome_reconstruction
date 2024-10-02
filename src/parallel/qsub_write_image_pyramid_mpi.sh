#!/bin/bash

if [[ "$#" < 4 ]]
then
    echo "usage: qsub_write_image_pyramid_mpi.sh n_procs img_root_dir parent_level output_format [abort_all_on_failure]"
    exit 1
fi

n_procs="$1"
img_root_dir="$2"
parent_level="$3"
output_format="$4"
if [[ "$#" < 5 ]]
then
    abort_all_on_failure="true"
else
    abort_all_on_failure="$5"
fi

script_dir="$(dirname "${BASH_SOURCE[0]}")"
script_name="$(basename "${BASH_SOURCE[0]}")"
cd "$script_dir"
project_dir="$(bash "../../util/get_project_dir.sh")"
cd "$project_dir/grid_build_master/parallel"
write_image_pyramid_mpi="write_image_pyramid_mpi.sh"
if [ ! -e "$write_image_pyramid_mpi" ]
then
    echo "$write_image_pyramid_mpi not found"
    exit 1
fi

job_name="write_image_pyramid_mpi"
qsub -cwd -c n -N "$job_name" -pe mpich "$n_procs" "$write_image_pyramid_mpi" "$script_name" "$img_root_dir" "$parent_level" "$output_format" "$abort_all_on_failure"