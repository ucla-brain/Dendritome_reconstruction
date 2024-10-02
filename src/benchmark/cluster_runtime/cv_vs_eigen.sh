#!/bin/bash

# Export all environment variables
#$ -V
#$ -o ~/log/cv_vs_eigen.o$JOB_ID
#$ -e ~/log/cv_vs_eigen.e$JOB_ID
#$ -N cv_vs_eigen
#$ -pe mpich 400

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
echo "working directory: $(pwd)"

project_dir="$(bash "../../util/get_project_dir.sh")"
build_dir="$project_dir/grid_build"
cv_vs_eigen_bin="$build_dir/benchmark/cluster_run_time/cv_vs_eigen"

echo "getting qhost command output"
qhost > "$project_dir/qhost.out"
echo "generating hostfile for mpiexec"
# generate a hostfile suitable for 100 processes, with 4 threads per process
python "$project_dir/util/gen_hostfile.py" 100 4
hostfile="$project_dir/hostfile"
# mpiexec with 100 processes (-n 100), each process spawning 4 threads (-bind-to core:4)
/usr/local/mpich-3.1.3/bin/mpiexec  -f "$hostfile" -n 100 -bind-to socket -map-by socket "$cv_vs_eigen_bin" 4 parallel
exit 0