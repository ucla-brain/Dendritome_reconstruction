#!/usr/bin/env bash

#$ -V
#$ -o $HOME/log/$JOB_NAME.o$JOB_ID
#$ -e $HOME/log/$JOB_NAME.e$JOB_ID

project_dir="$(bash "../../util/get_project_dir.sh")"
build_dir="$project_dir/grid_build"
util_dir="$project_dir/util"
verify_caller_script="$util_dir/verify_calling_script.sh"
source "$verify_caller_script"


mpi_io_test_bin="$build_dir/test/infrastructure/mpi_io_test"
if [ ! -e "$mpi_io_test_bin" ]
then
    echo "executable $mpi_io_test_bin not found"
    exit 1
fi

export MPICH_MPIIO_HINTS_DISPLAY=1
python "$util_dir/gen_hostfile.py" "$NSLOTS" 1
romio_json_path="$project_dir/configs/romio_hints.json"
if [ ! -e "$romio_json_path" ]
then
    echo "romio hints json not found: $romio_json_path"
    exit 1
fi
echo "MPICH io hints: $MPICH_MPIIO_HINTS"
/usr/local/mpich-3.1.3/bin/mpiexec -prepend-rank -print-all-exitcodes -verbose -n "$NSLOTS" "$mpi_io_test_bin" "$romio_json_path"
mpi_io_dir="$project_dir/test/infrastructure/mpi_io"
if [ ! -d "$mpi_io_dir" ]
then
    echo "$mpi_io_dir does not exist"
    exit 1
fi

nslot_log="$mpi_io_dir"'/'"$NSLOTS"'writers.master_log'
execution_log="$mpi_io_dir"'/'"$NSLOTS"'writers/output/logs/mpi_io_ints.master_log'
if [ ! -f "$execution_log" ]
then
    echo "$execution_log does not exist"
    exit 1
fi
blank=$(cat <<EOF



EOF
)
if [ -e "$nslot_log" ]
then
    echo "$blank" >> "$nslot_log"
fi
echo "qsub binding stride: $QSUB_BINDING_STRIDING" >> "$nslot_log"
cat "$execution_log" >> "$nslot_log"
exit 0
