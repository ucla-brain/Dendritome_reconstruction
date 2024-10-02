#!/bin/bash
#$ -V
#$ -o ~/log/conv_thread_number.o$JOB_ID
#$ -e ~/log/conv_thread_number.e$JOB_ID
#$ -N conv_thread_number
#$ -pe mpich 768

script_name="$(basename "${BASH_SOURCE[0]}")"
mother_script_name="qsub_$script_name"
# command used to call this current script
caller_command="$(ps -o args= $PPID)"
if [[ ! "$caller_command" = *"$mother_script_name"* ]]
then
    echo "$script_name must be called through $mother_script_name"
    exit 1
fi

export MPD_CON_EXT="sge_$JOB_ID.$SGE_TASK_ID"
echo "Got $NSLOTS slots."
echo "current working directory: $(pwd)"

project_dir="$(bash "../../util/get_project_dir.sh")"
build_dir="$project_dir/grid_build"
conv_thread_number_bin="$build_dir/benchmark/cluster_runtime/conv_thread_number"

/usr/local/mpich-3.1.3/bin/mpiexec  -n 192 -bind-to core:4 -map-by socket "$conv_thread_number_bin" parallel
exit 0