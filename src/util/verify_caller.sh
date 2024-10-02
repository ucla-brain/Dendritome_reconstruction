#!/usr/bin/env bash

# this script should only be sourced
# when a script is submitted through SGE,
# script_name="$(basename "${BASH_SOURCE[0]}")" does not recover script name.
# job id is used as name by SGE during execution
# caller_command="$(ps -o args= $PPID)" does not recover calling command of
# parent process as intended aka bash qsub_hdf5_parallel_test.sh
# will get a SGE command such as sge_shepherd-3146640 -bg
# every qsub script should pass its script name as first argument to
# corresponding execution script containing mpiexec command

if [[ "${BASH_SOURCE[0]}" = "${0}" ]]
then
    echo "${BASH_SOURCE[0]} should only be sourced"
    exit 1
fi

export MPD_CON_EXT="sge_$JOB_ID.$SGE_TASK_ID"
echo "Got $NSLOTS slots."
echo "working directory: $(pwd)"
if [ ! -d "$project_dir" ]
then
    echo "project root directory does not exist: $project_dir"
    exit 1
fi
echo "project root directory: $project_dir"

script_name="$JOB_NAME"
if [[ -z "$script_name" ]]
then
    echo "$(basename "${BASH_SOURCE[0]}") is intended to be used through SGE"
    exit 1
fi

mother_script_name="qsub_$script_name"
if [[ $# < 1 ]]
then
    echo "too few arguments. required arguments: calling script name + [script specific arguments]"
    exit 1
else
    caller_sciprt="$1"
    if [[ ! "$caller_sciprt" = "$mother_script_name.sh" ]]
    then
        echo "$script_name.sh must be called through $mother_script_name.sh"
        exit 1
    fi
fi
