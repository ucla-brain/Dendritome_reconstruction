#!/bin/bash

#$ -V
#$ -o $HOME/log/$JOB_NAME.o$JOB_ID
#$ -e $HOME/log/$JOB_NAME.e$JOB_ID

if [[ "$#" < 4 ]]
then
    echo "usage: write_image_pyramid_mpi.sh caller img_root_dir parent_level output_format [abort_all_on_failure]"
    exit 1
fi

project_dir="$(bash "../../util/get_project_dir.sh")"
build_dir="$project_dir/grid_build_master"
util_dir="$project_dir/util"
verify_caller_script="$util_dir/verify_caller.sh"
source "$verify_caller_script"

write_pyramid_bin="$build_dir/parallel/write_image_pyramid_mpi"

if [ ! -e "$write_pyramid_bin" ]
then
    echo "executable $write_pyramid_bin not found"
    exit 1
fi

img_root_dir="$2"
parent_level="$3"
output_format="$4"
if [[ "$#" < 5 ]]
then
    abort_all_on_failure="true"
else
    abort_all_on_failure="$5"
fi

/usr/local/mpich-3.1.3/bin/mpiexec -n "$NSLOTS" "$write_pyramid_bin" "$img_root_dir" "$parent_level" "$output_format" "$abort_all_on_failure"