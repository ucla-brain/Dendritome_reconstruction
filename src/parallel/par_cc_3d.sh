#!/bin/bash
#$ -V
#$ -o ~/log/par_cc3d.o$JOB_ID
#$ -e ~/log/par_cc3d.e$JOB_ID
#$ -N par_cc3d
#$ -pe mpich 52

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
par_cc_3d_bin="$build_dir/parallel/par_cc_3d"

echo "getting qhost command output"
qhost > "$project_dir/qhost.out"
echo "generating hostfile for mpiexec"
# generate a hostfile suitable for 100 processes, with 4 threads per process
python "$project_dir/util/gen_hostfile.py" 52 1
# The order of arguments is important. First global, then local options.
tissue_dir="/$project_dir/test/cc_seg"
/usr/local/mpich-3.1.3/bin/mpiexec  -n 52 -bind-to socket -map-by socket "$par_cc_3d_bin" "$tissue_dir" 22528 512 2
exit 0